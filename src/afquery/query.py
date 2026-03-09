import json
import sqlite3
from pathlib import Path

import duckdb
from pyroaring import BitMap

from .bitmaps import deserialize
from .capture import CaptureIndex, load_capture_indices
from .constants import normalize_chrom
from .models import QueryParams, QueryResult, VariantKey, Sample, Technology, SampleFilter
from .ploidy import compute_AN

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

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _build_sample_bitmap(self, sf: SampleFilter) -> BitMap:
        """Construye bitmap de samples elegibles según filtros de fenotipo, sexo y tecnología."""
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

        return ph_bm & sex_bm & tech_bm

    def _compute_eligible(
        self,
        chrom: str,
        pos: int,
        sample_bitmap: BitMap,
    ) -> tuple[BitMap, int]:
        """Return (eligible_bitmap, AN) for a single position."""
        covered = BitMap()
        for tech_id, capture_idx in self._capture.items():
            if capture_idx.covers(chrom, pos):
                covered |= self._tech_bitmaps.get(str(tech_id), BitMap())
        eligible = sample_bitmap & covered
        sex_bms = self._bitmaps.get("sex", {})
        AN = compute_AN(
            eligible,
            sex_bms.get("male", BitMap()),
            sex_bms.get("female", BitMap()),
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
        if chrom in self._partitioned_chroms:
            chrom_dir = self._db / "variants" / chrom
            if any(chrom_dir.glob("bucket_*.parquet")):
                return str(chrom_dir / "bucket_*.parquet")
            return None
        flat = self._db / "variants" / f"{chrom}.parquet"
        return str(flat) if flat.exists() else None

    # ------------------------------------------------------------------
    # Public query methods
    # ------------------------------------------------------------------

    def query(self, params: QueryParams) -> list[QueryResult]:
        chrom = normalize_chrom(params.chrom)
        pos = params.pos

        sample_bm = self._build_sample_bitmap(params.filter)
        eligible, AN = self._compute_eligible(chrom, pos, sample_bm)

        if AN == 0:
            return []

        parquet_path = self._parquet_path(chrom, pos)
        if parquet_path is None:
            return []

        con = duckdb.connect(config={"temp_directory": "/tmp"})
        rows = con.execute(
            "SELECT ref, alt, het_bitmap, hom_bitmap FROM read_parquet(?) WHERE pos = ?",
            [str(parquet_path), pos],
        ).fetchall()
        con.close()

        if not rows:
            return []

        results = []
        for ref, alt, het_bytes, hom_bytes in rows:
            het_bm = deserialize(bytes(het_bytes))
            hom_bm = deserialize(bytes(hom_bytes))
            AC = len(het_bm & eligible) + 2 * len(hom_bm & eligible)
            AF = AC / AN
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC,
                AN=AN,
                AF=AF,
                n_samples_eligible=len(eligible),
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
                f"SELECT pos, ref, alt, het_bitmap, hom_bitmap"
                f" FROM read_parquet('{escaped_glob}') WHERE pos IN ({placeholders})",
                valid_positions,
            ).fetchall()
        else:
            con.execute("CREATE TEMP TABLE pos_filter (pos UINTEGER)")
            con.executemany("INSERT INTO pos_filter VALUES (?)", [(p,) for p in valid_positions])
            rows = con.execute(
                f"SELECT v.pos, v.ref, v.alt, v.het_bitmap, v.hom_bitmap"
                f" FROM read_parquet('{escaped_glob}') v JOIN pos_filter pf ON v.pos = pf.pos"
            ).fetchall()
        con.close()

        results = []
        for pos, ref, alt, het_bytes, hom_bytes in rows:
            if (pos, ref, alt) not in requested_variants:
                continue
            eligible, AN = pos_data[pos]
            het_bm = deserialize(bytes(het_bytes))
            hom_bm = deserialize(bytes(hom_bytes))
            AC = len(het_bm & eligible) + 2 * len(hom_bm & eligible)
            AF = AC / AN
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC,
                AN=AN,
                AF=AF,
                n_samples_eligible=len(eligible),
            ))
        results.sort(key=lambda r: (r.variant.pos, r.variant.alt))
        return results

    def query_region(
        self,
        chrom: str,
        start: int,
        end: int,
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query all variants in a 1-based inclusive range [start, end]."""
        chrom = normalize_chrom(chrom)

        parquet_glob = self._parquet_glob(chrom)
        if parquet_glob is None:
            return []

        escaped_glob = parquet_glob.replace("'", "''")

        sample_bm = self._build_sample_bitmap(sf)

        con = duckdb.connect(config={"temp_directory": "/tmp"})
        rows = con.execute(
            f"SELECT pos, ref, alt, het_bitmap, hom_bitmap"
            f" FROM read_parquet('{escaped_glob}') WHERE pos BETWEEN ? AND ?",
            [start, end],
        ).fetchall()
        con.close()

        if not rows:
            return []

        pos_data: dict[int, tuple[BitMap, int]] = {}
        for pos, *_ in rows:
            if pos not in pos_data:
                eligible, AN = self._compute_eligible(chrom, pos, sample_bm)
                pos_data[pos] = (eligible, AN)

        results = []
        for pos, ref, alt, het_bytes, hom_bytes in rows:
            eligible, AN = pos_data[pos]
            if AN == 0:
                continue
            het_bm = deserialize(bytes(het_bytes))
            hom_bm = deserialize(bytes(hom_bytes))
            AC = len(het_bm & eligible) + 2 * len(hom_bm & eligible)
            AF = AC / AN
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC,
                AN=AN,
                AF=AF,
                n_samples_eligible=len(eligible),
            ))
        results.sort(key=lambda r: (r.variant.pos, r.variant.alt))
        return results
