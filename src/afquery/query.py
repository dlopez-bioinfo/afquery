import json
import sqlite3
from pathlib import Path

import duckdb
from pyroaring import BitMap

from .bitmaps import deserialize
from .capture import CaptureIndex, load_capture_indices
from .constants import normalize_chrom
from .models import QueryParams, QueryResult, VariantKey, Sample, Technology
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

    def _build_icd10_bitmap(self, icd10_codes: list[str]) -> BitMap:
        bm = BitMap()
        for code in icd10_codes:
            if code in self._bitmaps.get("icd10", {}):
                bm |= self._bitmaps["icd10"][code]
        return bm

    def _build_sex_bitmap(self, sex_filter: str) -> BitMap:
        sex_bms = self._bitmaps.get("sex", {})
        if sex_filter == "both":
            return sex_bms.get("male", BitMap()) | sex_bms.get("female", BitMap())
        return sex_bms.get(sex_filter, BitMap())

    def _compute_eligible(
        self,
        chrom: str,
        pos: int,
        icd10_bitmap: BitMap,
        sex_bitmap: BitMap,
    ) -> tuple[BitMap, int]:
        """Return (eligible_bitmap, AN) for a single position."""
        covered = BitMap()
        for tech_id, capture_idx in self._capture.items():
            if capture_idx.covers(chrom, pos):
                covered |= self._tech_bitmaps.get(str(tech_id), BitMap())
        eligible = icd10_bitmap & sex_bitmap & covered
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

        icd10_bitmap = self._build_icd10_bitmap(params.icd10_codes)
        sex_bitmap = self._build_sex_bitmap(params.sex_filter)
        eligible, AN = self._compute_eligible(chrom, pos, icd10_bitmap, sex_bitmap)

        if AN == 0:
            return []

        parquet_path = self._parquet_path(chrom, pos)
        if parquet_path is None:
            return []

        con = duckdb.connect()
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
        return results

    def query_batch(
        self,
        chrom: str,
        positions: list[int],
        icd10_codes: list[str],
        sex_filter: str = "both",
    ) -> list[QueryResult]:
        """Query multiple positions in a single Parquet read."""
        chrom = normalize_chrom(chrom)
        unique_positions = list(dict.fromkeys(positions))  # dedup, preserve order

        icd10_bitmap = self._build_icd10_bitmap(icd10_codes)
        sex_bitmap = self._build_sex_bitmap(sex_filter)

        pos_data: dict[int, tuple[BitMap, int]] = {}
        for pos in unique_positions:
            eligible, AN = self._compute_eligible(chrom, pos, icd10_bitmap, sex_bitmap)
            pos_data[pos] = (eligible, AN)

        valid_positions = [p for p in unique_positions if pos_data[p][1] > 0]
        if not valid_positions:
            return []

        parquet_glob = self._parquet_glob(chrom)
        if parquet_glob is None:
            return []

        escaped_glob = parquet_glob.replace("'", "''")

        con = duckdb.connect()
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
        icd10_codes: list[str],
        sex_filter: str = "both",
    ) -> list[QueryResult]:
        """Query all variants in a 1-based inclusive range [start, end]."""
        chrom = normalize_chrom(chrom)

        parquet_glob = self._parquet_glob(chrom)
        if parquet_glob is None:
            return []

        escaped_glob = parquet_glob.replace("'", "''")

        icd10_bitmap = self._build_icd10_bitmap(icd10_codes)
        sex_bitmap = self._build_sex_bitmap(sex_filter)

        con = duckdb.connect()
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
                eligible, AN = self._compute_eligible(chrom, pos, icd10_bitmap, sex_bitmap)
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
