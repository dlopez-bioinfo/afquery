import json
import sqlite3
import warnings
from pathlib import Path

import duckdb
from pyroaring import BitMap

from .bitmaps import deserialize
from .capture import CaptureIndex, load_capture_indices
from .constants import normalize_chrom
from .models import AfqueryWarning, QueryParams, QueryResult, VariantKey, Sample, Technology, SampleFilter
from .ploidy import compute_AN, split_ploidy

BATCH_IN_THRESHOLD = 10_000


class QueryEngine:
    def __init__(self, db_path: str):
        self._db = Path(db_path)
        self._manifest = json.loads((self._db / "manifest.json").read_text())
        self._genome_build = self._manifest["genome_build"]

        con = sqlite3.connect(self._db / "metadata.sqlite")

        # Load precomputed bitmaps
        self._bitmaps: dict[str, dict[str, BitMap]] = {}
        for btype, bkey, bdata in con.execute(
            "SELECT bitmap_type, bitmap_key, bitmap_data FROM precomputed_bitmaps"
        ).fetchall():
            self._bitmaps.setdefault(btype, {})[bkey] = deserialize(bytes(bdata))

        # Load samples and technologies
        samples = [
            Sample(r[0], r[1], r[2], r[3])
            for r in con.execute(
                "SELECT sample_id, sample_name, sex, tech_id FROM samples"
            ).fetchall()
        ]
        techs = [
            Technology(r[0], r[1], r[2])
            for r in con.execute(
                "SELECT tech_id, tech_name, bed_path FROM technologies"
            ).fetchall()
        ]
        con.close()

        self._samples = samples
        self._tech_map = {t.tech_id: t for t in techs}
        sex_bms = self._bitmaps.get("sex", {})
        self._male_bm = sex_bms.get("male", BitMap())
        self._female_bm = sex_bms.get("female", BitMap())
        self._capture = load_capture_indices(techs, str(self._db / "capture"))
        self._tech_bitmaps = self._bitmaps.get("tech", {})
        self._all_samples_bm = BitMap(s.sample_id for s in self._samples)
        self._tech_name_to_id: dict[str, str] = {
            t.tech_name: str(t.tech_id) for t in self._tech_map.values()
        }

        # Cache which chroms use partitioned storage (variants/{chrom}/ directory)
        self._partitioned_chroms: set[str] = set()
        variants_dir = self._db / "variants"
        if variants_dir.exists():
            for p in variants_dir.iterdir():
                if p.is_dir():
                    self._partitioned_chroms.add(p.name)

        self._flat_chroms: set[str] = set()
        if variants_dir.exists():
            for p in variants_dir.iterdir():
                if p.is_file() and p.suffix == ".parquet":
                    self._flat_chroms.add(p.stem)
        self._all_known_chroms: set[str] = self._partitioned_chroms | self._flat_chroms

        # Precompute covered bitmap for positions where ALL technologies cover (common case)
        # For each tech, store its bitmap; also precompute the union of all WGS tech bitmaps
        self._covered_all_bm: BitMap = BitMap()
        for tech_id, capture_idx in self._capture.items():
            if capture_idx._always_covered:
                bm = self._tech_bitmaps.get(str(tech_id), BitMap())
                self._covered_all_bm |= bm

        # Cache parquet glob patterns per chrom (avoids filesystem scan on every query)
        self._glob_cache: dict[str, str | None] = {}

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _validate_sample_filter(self, sf: SampleFilter) -> None:
        ph_bms = self._bitmaps.get("phenotype", {})
        for code in sf.phenotype_include:
            if code not in ph_bms:
                warnings.warn(
                    f"Phenotype {code!r} not in database — include will match 0 samples.",
                    AfqueryWarning, stacklevel=4,
                )
        for code in sf.phenotype_exclude:
            if code not in ph_bms:
                warnings.warn(
                    f"Phenotype {code!r} not in database — exclude has no effect.",
                    AfqueryWarning, stacklevel=4,
                )
        for name in sf.tech_include:
            if name not in self._tech_name_to_id:
                warnings.warn(
                    f"Technology {name!r} not in database — include will match 0 samples.",
                    AfqueryWarning, stacklevel=4,
                )
        for name in sf.tech_exclude:
            if name not in self._tech_name_to_id:
                warnings.warn(
                    f"Technology {name!r} not in database — exclude has no effect.",
                    AfqueryWarning, stacklevel=4,
                )
        if sf.sex not in ("male", "female", "both"):
            raise ValueError(f"Invalid sex {sf.sex!r}. Must be 'male', 'female', or 'both'.")

    def _build_sample_bitmap(self, sf: SampleFilter) -> BitMap:
        """Construye bitmap de samples elegibles según filtros de fenotipo, sexo y tecnología."""
        self._validate_sample_filter(sf)
        # --- Fenotipo ---
        ph_bms = self._bitmaps.get("phenotype", {})
        if sf.phenotype_include:
            ph_bm = BitMap()
            for code in sf.phenotype_include:
                if code in ph_bms:
                    ph_bm |= ph_bms[code]
        else:
            ph_bm = BitMap(self._all_samples_bm)  # copia, todas las muestras
        for code in sf.phenotype_exclude:
            if code in ph_bms:
                ph_bm -= ph_bms[code]

        # --- Sexo ---
        sex_bms = self._bitmaps.get("sex", {})
        if sf.sex == "both":
            sex_bm = sex_bms.get("male", BitMap()) | sex_bms.get("female", BitMap())
        else:
            sex_bm = sex_bms.get(sf.sex, BitMap())

        # --- Tecnología (nivel de sample; cobertura por posición sigue en _compute_eligible) ---
        tech_bms = self._bitmaps.get("tech", {})
        if sf.tech_include or sf.tech_exclude:
            if sf.tech_include:
                tech_bm = BitMap()
                for name in sf.tech_include:
                    tid = self._tech_name_to_id.get(name)
                    if tid and tid in tech_bms:
                        tech_bm |= tech_bms[tid]
            else:
                tech_bm = BitMap(self._all_samples_bm)
            for name in sf.tech_exclude:
                tid = self._tech_name_to_id.get(name)
                if tid and tid in tech_bms:
                    tech_bm -= tech_bms[tid]
        else:
            tech_bm = self._all_samples_bm  # sin copia — solo lectura en la intersección

        result_bm = ph_bm & sex_bm & tech_bm
        if len(result_bm) == 0:
            warnings.warn(
                "Sample filter produces an empty eligible set — all queries will return AN=0. "
                "Check that include/exclude filters are not contradictory.",
                AfqueryWarning, stacklevel=3,
            )
        return result_bm

    def _compute_eligible(
        self,
        chrom: str,
        pos: int,
        sample_bitmap: BitMap,
    ) -> tuple[BitMap, int]:
        """Return (eligible_bitmap, AN) for a single position."""
        covered = BitMap(self._covered_all_bm)  # start with WGS samples (always covered)
        for tech_id, capture_idx in self._capture.items():
            if capture_idx._always_covered:
                continue  # already included in _covered_all_bm
            if capture_idx.covers(chrom, pos):
                covered |= self._tech_bitmaps.get(str(tech_id), BitMap())
        eligible = sample_bitmap & covered
        AN = compute_AN(
            eligible,
            self._male_bm,
            self._female_bm,
            chrom, pos, self._genome_build,
        )
        return eligible, AN

    def _parquet_path(self, chrom: str, pos: int) -> Path | None:
        """Resolve Parquet path for a point query. Partitioned takes priority over flat."""
        if chrom in self._partitioned_chroms:
            bucket = pos // 1_000_000
            p = self._db / "variants" / chrom / f"bucket_{bucket}.parquet"
            return p if p.exists() else None
        flat = self._db / "variants" / f"{chrom}.parquet"
        return flat if flat.exists() else None

    def _parquet_glob(self, chrom: str) -> str | None:
        """Return path/glob pattern for batch/region queries. None if no data for chrom."""
        if chrom in self._glob_cache:
            return self._glob_cache[chrom]
        if chrom in self._partitioned_chroms:
            chrom_dir = self._db / "variants" / chrom
            result = str(chrom_dir / "bucket_*.parquet") if any(chrom_dir.glob("bucket_*.parquet")) else None
        else:
            flat = self._db / "variants" / f"{chrom}.parquet"
            result = str(flat) if flat.exists() else None
        self._glob_cache[chrom] = result
        return result

    # ------------------------------------------------------------------
    # Public query methods
    # ------------------------------------------------------------------

    def query(self, params: QueryParams) -> list[QueryResult]:
        chrom = normalize_chrom(params.chrom)
        pos = params.pos

        if chrom not in self._all_known_chroms:
            warnings.warn(
                f"Chromosome {chrom!r} has no data in this database. "
                f"Available: {sorted(self._all_known_chroms)!r}",
                AfqueryWarning, stacklevel=3,
            )
            return []

        sample_bm = self._build_sample_bitmap(params.filter)
        eligible, AN = self._compute_eligible(chrom, pos, sample_bm)

        if AN == 0:
            return []

        parquet_path = self._parquet_path(chrom, pos)
        if parquet_path is None:
            return []

        con = duckdb.connect(config={"temp_directory": "/tmp"})
        rows = con.execute(
            "SELECT ref, alt, het_bitmap, hom_bitmap, fail_bitmap FROM read_parquet(?) WHERE pos = ?",
            [str(parquet_path), pos],
        ).fetchall()
        con.close()

        if not rows:
            return []

        results = []
        for row in rows:
            ref, alt, het_bytes, hom_bytes, fail_bytes = row
            het_bm = deserialize(bytes(het_bytes))
            hom_bm = deserialize(bytes(hom_bytes))
            haploid_elig, diploid_elig = split_ploidy(
                eligible, self._male_bm, self._female_bm, chrom, pos, self._genome_build
            )
            het_elig = het_bm & eligible
            hom_elig = hom_bm & eligible
            AC = (len((het_elig | hom_elig) & haploid_elig)
                  + len(het_elig & diploid_elig)
                  + 2 * len(hom_elig & diploid_elig))
            N_HET = len(het_elig & diploid_elig)
            N_HOM_ALT = len(hom_elig & diploid_elig) + len((het_elig | hom_elig) & haploid_elig)
            AF = AC / AN if AN > 0 else None
            fail_bm = deserialize(bytes(fail_bytes))
            N_FAIL = len(fail_bm & eligible)
            N_HOM_REF = len(eligible) - N_HET - N_HOM_ALT - N_FAIL
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC, AN=AN, AF=AF, n_samples_eligible=len(eligible),
                N_HET=N_HET, N_HOM_ALT=N_HOM_ALT, N_HOM_REF=N_HOM_REF, N_FAIL=N_FAIL,
            ))
        if params.ref is not None:
            results = [r for r in results if r.variant.ref == params.ref]
        if params.alt is not None:
            results = [r for r in results if r.variant.alt == params.alt]
        return results

    def query_batch(
        self,
        chrom: str,
        variants: list[tuple[int, str, str]],
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query multiple variants (pos, ref, alt) in a single Parquet read."""
        chrom = normalize_chrom(chrom)

        if chrom not in self._all_known_chroms:
            warnings.warn(
                f"Chromosome {chrom!r} has no data in this database. "
                f"Available: {sorted(self._all_known_chroms)!r}",
                AfqueryWarning, stacklevel=3,
            )
            return []

        # Deduplicate by full variant, preserve order
        unique_variants = list(dict.fromkeys(variants))

        # AN is per-position; collect unique positions
        unique_positions = list(dict.fromkeys(pos for pos, _ref, _alt in unique_variants))

        sample_bm = self._build_sample_bitmap(sf)

        pos_data: dict[int, tuple[BitMap, int]] = {}
        for pos in unique_positions:
            eligible, AN = self._compute_eligible(chrom, pos, sample_bm)
            pos_data[pos] = (eligible, AN)

        valid_positions = [p for p in unique_positions if pos_data[p][1] > 0]
        if not valid_positions:
            return []

        # Only keep variants whose position has AN > 0
        requested_variants: set[tuple[int, str, str]] = {
            (pos, ref, alt) for pos, ref, alt in unique_variants if pos in set(valid_positions)
        }

        parquet_glob = self._parquet_glob(chrom)
        if parquet_glob is None:
            return []

        escaped_glob = parquet_glob.replace("'", "''")

        con = duckdb.connect(config={"temp_directory": "/tmp"})
        if len(valid_positions) < BATCH_IN_THRESHOLD:
            placeholders = ", ".join("?" * len(valid_positions))
            rows = con.execute(
                f"SELECT pos, ref, alt, het_bitmap, hom_bitmap, fail_bitmap"
                f" FROM read_parquet('{escaped_glob}') WHERE pos IN ({placeholders})",
                valid_positions,
            ).fetchall()
        else:
            con.execute("CREATE TEMP TABLE pos_filter (pos UINTEGER)")
            con.executemany("INSERT INTO pos_filter VALUES (?)", [(p,) for p in valid_positions])
            rows = con.execute(
                f"SELECT v.pos, v.ref, v.alt, v.het_bitmap, v.hom_bitmap, v.fail_bitmap"
                f" FROM read_parquet('{escaped_glob}') v JOIN pos_filter pf ON v.pos = pf.pos"
            ).fetchall()
        con.close()

        results = []
        for row in rows:
            pos, ref, alt, het_bytes, hom_bytes, fail_bytes = row
            if (pos, ref, alt) not in requested_variants:
                continue
            eligible, AN = pos_data[pos]
            het_bm = deserialize(bytes(het_bytes))
            hom_bm = deserialize(bytes(hom_bytes))
            haploid_elig, diploid_elig = split_ploidy(
                eligible, self._male_bm, self._female_bm, chrom, pos, self._genome_build
            )
            het_elig = het_bm & eligible
            hom_elig = hom_bm & eligible
            AC = (len((het_elig | hom_elig) & haploid_elig)
                  + len(het_elig & diploid_elig)
                  + 2 * len(hom_elig & diploid_elig))
            N_HET = len(het_elig & diploid_elig)
            N_HOM_ALT = len(hom_elig & diploid_elig) + len((het_elig | hom_elig) & haploid_elig)
            AF = AC / AN if AN > 0 else None
            fail_bm = deserialize(bytes(fail_bytes))
            N_FAIL = len(fail_bm & eligible)
            N_HOM_REF = len(eligible) - N_HET - N_HOM_ALT - N_FAIL
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC, AN=AN, AF=AF, n_samples_eligible=len(eligible),
                N_HET=N_HET, N_HOM_ALT=N_HOM_ALT, N_HOM_REF=N_HOM_REF, N_FAIL=N_FAIL,
            ))
        results.sort(key=lambda r: (r.variant.pos, r.variant.alt))
        return results

    def query_region_multi(
        self,
        regions: list[tuple[str, int, int]],
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query variants across multiple genomic regions (may span chromosomes).

        Overlapping regions are deduplicated by variant key. Results are sorted
        by (chrom, pos, alt).

        Args:
            regions: List of ``(chrom, start, end)`` tuples, 1-based inclusive.
            sf: Sample filter.

        Returns:
            List of :class:`~afquery.models.QueryResult` objects.
        """
        seen: set[tuple[str, int, str, str]] = set()
        results: list[QueryResult] = []
        for chrom, start, end in regions:
            for r in self.query_region(chrom, start, end, sf):
                key = (r.variant.chrom, r.variant.pos, r.variant.ref, r.variant.alt)
                if key not in seen:
                    seen.add(key)
                    results.append(r)
        results.sort(key=lambda r: (r.variant.chrom, r.variant.pos, r.variant.alt))
        return results

    def query_batch_multi(
        self,
        variants: list[tuple[str, int, str, str]],
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query variants across multiple chromosomes, preserving input order."""
        from collections import defaultdict
        indexed = list(enumerate(variants))
        by_chrom: dict[str, list[tuple[int, tuple[int, str, str]]]] = defaultdict(list)
        for idx, (chrom, pos, ref, alt) in indexed:
            by_chrom[chrom].append((idx, (pos, ref, alt)))
        tagged: list[tuple[int, QueryResult]] = []
        for chrom, idx_variants in by_chrom.items():
            idxs = [i for i, _ in idx_variants]
            per_chrom = [v for _, v in idx_variants]
            results = self.query_batch(chrom, per_chrom, sf)
            result_map = {(r.variant.pos, r.variant.ref, r.variant.alt): r for r in results}
            for i, (pos, ref, alt) in zip(idxs, per_chrom):
                if (pos, ref, alt) in result_map:
                    tagged.append((i, result_map[(pos, ref, alt)]))
        tagged.sort(key=lambda x: x[0])
        return [r for _, r in tagged]

    def query_region(
        self,
        chrom: str,
        start: int,
        end: int,
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query all variants in a 1-based inclusive range [start, end]."""
        chrom = normalize_chrom(chrom)

        if chrom not in self._all_known_chroms:
            warnings.warn(
                f"Chromosome {chrom!r} has no data in this database. "
                f"Available: {sorted(self._all_known_chroms)!r}",
                AfqueryWarning, stacklevel=3,
            )
            return []

        parquet_glob = self._parquet_glob(chrom)
        if parquet_glob is None:
            return []

        escaped_glob = parquet_glob.replace("'", "''")

        sample_bm = self._build_sample_bitmap(sf)

        con = duckdb.connect(config={"temp_directory": "/tmp"})
        rows = con.execute(
            f"SELECT pos, ref, alt, het_bitmap, hom_bitmap, fail_bitmap"
            f" FROM read_parquet('{escaped_glob}') WHERE pos BETWEEN ? AND ?",
            [start, end],
        ).fetchall()
        con.close()

        if not rows:
            return []

        pos_data: dict[int, tuple[BitMap, int]] = {}
        for row in rows:
            pos = row[0]
            if pos not in pos_data:
                eligible, AN = self._compute_eligible(chrom, pos, sample_bm)
                pos_data[pos] = (eligible, AN)

        results = []
        for row in rows:
            pos, ref, alt, het_bytes, hom_bytes, fail_bytes = row
            eligible, AN = pos_data[pos]
            if AN == 0:
                continue
            het_bm = deserialize(bytes(het_bytes))
            hom_bm = deserialize(bytes(hom_bytes))
            haploid_elig, diploid_elig = split_ploidy(
                eligible, self._male_bm, self._female_bm, chrom, pos, self._genome_build
            )
            het_elig = het_bm & eligible
            hom_elig = hom_bm & eligible
            AC = (len((het_elig | hom_elig) & haploid_elig)
                  + len(het_elig & diploid_elig)
                  + 2 * len(hom_elig & diploid_elig))
            N_HET = len(het_elig & diploid_elig)
            N_HOM_ALT = len(hom_elig & diploid_elig) + len((het_elig | hom_elig) & haploid_elig)
            AF = AC / AN if AN > 0 else None
            fail_bm = deserialize(bytes(fail_bytes))
            N_FAIL = len(fail_bm & eligible)
            N_HOM_REF = len(eligible) - N_HET - N_HOM_ALT - N_FAIL
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC, AN=AN, AF=AF, n_samples_eligible=len(eligible),
                N_HET=N_HET, N_HOM_ALT=N_HOM_ALT, N_HOM_REF=N_HOM_REF, N_FAIL=N_FAIL,
            ))
        results.sort(key=lambda r: (r.variant.pos, r.variant.alt))
        return results
