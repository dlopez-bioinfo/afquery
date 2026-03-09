import glob as glob_module
import json
import logging
import os
import sqlite3
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from pyroaring import BitMap

from ..bitmaps import deserialize, serialize, build_icd10_bitmaps, build_sex_bitmaps, build_tech_bitmaps
from ..constants import VALID_GENOME_BUILDS
from ..models import Sample, Technology
from .build import PARQUET_SCHEMA, get_chroms_in_temp_files
from .ingest import ingest_all
from .manifest import parse_manifest
from .regions import build_capture_indices

logger = logging.getLogger(__name__)


class UpdateError(RuntimeError):
    pass


@dataclass
class CheckResult:
    severity: str   # 'error' | 'warning' | 'info'
    message: str


# ---- Internal helpers ----

def _read_manifest(db_dir: str) -> dict:
    path = os.path.join(db_dir, "manifest.json")
    with open(path) as f:
        return json.load(f)


def _update_manifest(
    db_dir: str,
    sample_count: int,
    next_sample_id: int | None = None,
) -> None:
    path = os.path.join(db_dir, "manifest.json")
    try:
        with open(path) as f:
            manifest = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        manifest = {}
    manifest["sample_count"] = sample_count
    manifest["updated_at"] = datetime.now(timezone.utc).isoformat()
    if next_sample_id is not None:
        manifest["next_sample_id"] = next_sample_id
    tmp_path = path + ".tmp"
    with open(tmp_path, "w") as f:
        json.dump(manifest, f, indent=2)
    os.replace(tmp_path, path)


def _get_next_sample_id(con: sqlite3.Connection) -> int:
    row = con.execute("SELECT MAX(sample_id) FROM samples").fetchone()
    return 0 if row[0] is None else row[0] + 1


def _get_next_tech_id(con: sqlite3.Connection) -> int:
    row = con.execute("SELECT MAX(tech_id) FROM technologies").fetchone()
    return 0 if row[0] is None else row[0] + 1


def _regenerate_precomputed_bitmaps(con: sqlite3.Connection) -> None:
    samples_rows = con.execute(
        "SELECT sample_id, sample_name, sex, tech_id FROM samples"
    ).fetchall()
    samples = [Sample(r[0], r[1], r[2], r[3]) for r in samples_rows]

    icd10_pairs = con.execute(
        "SELECT sample_id, icd10_code FROM sample_icd10"
    ).fetchall()

    con.execute("DELETE FROM precomputed_bitmaps")

    for sex, bm in build_sex_bitmaps(samples).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("sex", sex, serialize(bm)),
        )

    for code, bm in build_icd10_bitmaps(icd10_pairs).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("icd10", code, serialize(bm)),
        )

    for tech_id, bm in build_tech_bitmaps(samples).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("tech", str(tech_id), serialize(bm)),
        )


def _merge_chromosome_parquet(
    chrom: str,
    db_dir: str,
    update_tmp_dir: str,
    row_group_size: int = 100_000,
) -> tuple[int, int]:
    """Merge new temp files into existing chrom Parquet. Returns (new_variants, updated_variants)."""
    variants_dir = os.path.join(db_dir, "variants")
    out_path = os.path.join(variants_dir, f"{chrom}.parquet")

    # Read existing Parquet via pyarrow (NOT DuckDB — need Python bitmap deserialization)
    existing: dict[tuple, tuple[BitMap, BitMap]] = {}
    if os.path.exists(out_path):
        table = pq.read_table(out_path)
        for i in range(len(table)):
            pos = table["pos"][i].as_py()
            ref = table["ref"][i].as_py()
            alt = table["alt"][i].as_py()
            het_bm = deserialize(table["het_bitmap"][i].as_py())
            hom_bm = deserialize(table["hom_bitmap"][i].as_py())
            existing[(pos, ref, alt)] = (het_bm, hom_bm)

    # Check if there are any new temp files
    parquet_files = glob_module.glob(os.path.join(update_tmp_dir, "sample_*.parquet"))
    if not parquet_files:
        return (0, 0)

    # Aggregate new rows via DuckDB
    glob_pattern = os.path.join(update_tmp_dir, "sample_*.parquet").replace("'", "''")
    con = duckdb.connect()
    try:
        rows = con.execute(
            f"""
            SELECT pos, ref, alt,
                list(sample_id ORDER BY sample_id),
                list(gt_ac     ORDER BY sample_id)
            FROM read_parquet('{glob_pattern}')
            WHERE chrom = ?
            GROUP BY pos, ref, alt
            ORDER BY pos, alt
            """,
            [chrom],
        ).fetchall()
    finally:
        con.close()

    # No new data for this chrom — skip (don't touch the file)
    if not rows:
        return (0, 0)

    new_variants = 0
    updated_variants = 0

    for pos, ref, alt, sample_ids, gt_acs in rows:
        key = (pos, ref, alt)
        het_ids = [sid for sid, ac in zip(sample_ids, gt_acs) if ac == 1]
        hom_ids = [sid for sid, ac in zip(sample_ids, gt_acs) if ac == 2]
        new_het = BitMap(het_ids)
        new_hom = BitMap(hom_ids)

        if key in existing:
            old_het, old_hom = existing[key]
            existing[key] = (old_het | new_het, old_hom | new_hom)
            updated_variants += 1
        else:
            existing[key] = (new_het, new_hom)
            new_variants += 1

    # Sort by (pos, alt) and write atomically
    sorted_keys = sorted(existing.keys(), key=lambda k: (k[0], k[2]))

    positions = [k[0] for k in sorted_keys]
    refs = [k[1] for k in sorted_keys]
    alts = [k[2] for k in sorted_keys]
    het_bitmaps = [serialize(existing[k][0]) for k in sorted_keys]
    hom_bitmaps = [serialize(existing[k][1]) for k in sorted_keys]

    table = pa.table(
        {
            "pos":        pa.array(positions,   type=pa.uint32()),
            "ref":        pa.array(refs,        type=pa.large_utf8()),
            "alt":        pa.array(alts,        type=pa.large_utf8()),
            "het_bitmap": pa.array(het_bitmaps, type=pa.large_binary()),
            "hom_bitmap": pa.array(hom_bitmaps, type=pa.large_binary()),
        },
        schema=PARQUET_SCHEMA,
    )

    os.makedirs(variants_dir, exist_ok=True)
    tmp_path = out_path + ".tmp"
    pq.write_table(table, tmp_path, row_group_size=row_group_size)
    os.replace(tmp_path, out_path)

    return (new_variants, updated_variants)


def _clear_bits_from_parquet(parquet_file: str, removal_ids: BitMap) -> None:
    """Clear removal_ids bits from het/hom bitmaps in a Parquet file. Rewrites atomically if dirty."""
    table = pq.read_table(parquet_file)
    dirty = False

    new_het_list = []
    new_hom_list = []

    for i in range(len(table)):
        het_bytes = table["het_bitmap"][i].as_py()
        hom_bytes = table["hom_bitmap"][i].as_py()
        het_bm = deserialize(het_bytes)
        hom_bm = deserialize(hom_bytes)

        if removal_ids & (het_bm | hom_bm):
            het_bm = het_bm - removal_ids
            hom_bm = hom_bm - removal_ids
            dirty = True

        new_het_list.append(serialize(het_bm))
        new_hom_list.append(serialize(hom_bm))

    if not dirty:
        return

    new_table = pa.table(
        {
            "pos":        table["pos"],
            "ref":        table["ref"],
            "alt":        table["alt"],
            "het_bitmap": pa.array(new_het_list, type=pa.large_binary()),
            "hom_bitmap": pa.array(new_hom_list, type=pa.large_binary()),
        },
        schema=PARQUET_SCHEMA,
    )

    tmp_path = parquet_file + ".tmp"
    pq.write_table(new_table, tmp_path)
    os.replace(tmp_path, parquet_file)


# ---- Public functions ----

def add_samples(
    db_dir: str,
    manifest_path: str,
    threads: int = 8,
    tmp_dir: str | None = None,
    bed_dir: str | None = None,
    genome_build: str | None = None,
) -> dict:
    """Add new samples to the database without full rebuild."""
    # 1. Read DB manifest
    manifest = _read_manifest(db_dir)
    db_genome_build = manifest["genome_build"]

    # 2. Validate genome build if provided
    if genome_build is not None and genome_build != db_genome_build:
        raise UpdateError(
            f"genome_build mismatch: DB has '{db_genome_build}', got '{genome_build}'"
        )

    # 3. Parse new manifest
    samples_raw, techs_raw = parse_manifest(manifest_path, bed_dir)

    # 4. Open SQLite connection
    db_path = os.path.join(db_dir, "metadata.sqlite")
    con = sqlite3.connect(db_path)

    total_new = 0
    total_updated = 0

    try:
        # 5. Check for duplicate sample names
        existing_names = {
            r[0] for r in con.execute("SELECT sample_name FROM samples").fetchall()
        }
        new_names = [s.sample_name for s in samples_raw]
        duplicates = [n for n in new_names if n in existing_names]
        if duplicates:
            raise UpdateError(
                f"Sample(s) already in database: {', '.join(duplicates)}"
            )

        # 6. Assign IDs sequentially from next available.
        # Use manifest's next_sample_id when present (survives removals) so IDs
        # are never reused even after a sample has been deleted.
        manifest_next = manifest.get("next_sample_id")
        db_next = _get_next_sample_id(con)
        starting_id = max(manifest_next if manifest_next is not None else 0, db_next)
        new_samples = [
            Sample(
                sample_id=starting_id + i,
                sample_name=ps.sample_name,
                sex=ps.sex,
                tech_id=-1,  # filled below
            )
            for i, ps in enumerate(samples_raw)
        ]

        # 7. Handle technologies
        capture_dir = os.path.join(db_dir, "capture")
        os.makedirs(capture_dir, exist_ok=True)

        existing_techs = {
            r[1]: (r[0], r[2])
            for r in con.execute(
                "SELECT tech_id, tech_name, bed_path FROM technologies"
            ).fetchall()
        }

        tech_name_to_id: dict[str, int] = {}
        next_tech_id = _get_next_tech_id(con)

        for pt in techs_raw:
            if pt.tech_name in existing_techs:
                tech_name_to_id[pt.tech_name] = existing_techs[pt.tech_name][0]
            else:
                tech_id = next_tech_id
                next_tech_id += 1
                con.execute(
                    "INSERT INTO technologies VALUES (?, ?, ?)",
                    (tech_id, pt.tech_name, pt.bed_path),
                )
                tech = Technology(
                    tech_id=tech_id, tech_name=pt.tech_name, bed_path=pt.bed_path
                )
                build_capture_indices([tech], capture_dir)
                tech_name_to_id[pt.tech_name] = tech_id

        # Assign tech_ids to new samples
        for i, ps in enumerate(samples_raw):
            new_samples[i].tech_id = tech_name_to_id[ps.tech_name]

        vcf_paths = [ps.vcf_path for ps in samples_raw]

        # 8. Ingest VCFs into fresh tmp_dir
        auto_tmp = tmp_dir is None
        if auto_tmp:
            tmp_dir = tempfile.mkdtemp(prefix="afquery_update_")

        try:
            ingest_all(new_samples, vcf_paths, tmp_dir, n_workers=threads)

            # 9. Collect chroms from new temp files
            chroms = get_chroms_in_temp_files(tmp_dir)

            # 10. Merge Parquet files
            for chrom in chroms:
                n, u = _merge_chromosome_parquet(chrom, db_dir, tmp_dir)
                total_new += n
                total_updated += u
        finally:
            if auto_tmp:
                import shutil
                shutil.rmtree(tmp_dir, ignore_errors=True)

        # 11. Insert new samples and ICD10 pairs into SQLite
        con.executemany(
            "INSERT INTO samples VALUES (?, ?, ?, ?)",
            [(s.sample_id, s.sample_name, s.sex, s.tech_id) for s in new_samples],
        )
        sample_icd10_pairs = [
            (starting_id + i, code)
            for i, ps in enumerate(samples_raw)
            for code in ps.icd10_codes
        ]
        con.executemany("INSERT INTO sample_icd10 VALUES (?, ?)", sample_icd10_pairs)

        # 12. Regenerate precomputed bitmaps
        _regenerate_precomputed_bitmaps(con)
        con.commit()

    finally:
        con.close()

    # 13. Update manifest (persist next_sample_id so future adds don't reuse IDs)
    next_id = starting_id + len(new_samples)
    _update_manifest(db_dir, next_id, next_sample_id=next_id)

    return {
        "new_samples": len(new_samples),
        "new_variants": total_new,
        "updated_variants": total_updated,
    }


def remove_samples(db_dir: str, sample_names: list[str]) -> dict:
    """Remove samples from the database by clearing their bitmap bits."""
    if not sample_names:
        return {"removed": [], "not_found": []}

    manifest = _read_manifest(db_dir)

    db_path = os.path.join(db_dir, "metadata.sqlite")
    con = sqlite3.connect(db_path)

    try:
        # 1. Look up sample IDs for given names
        placeholders = ",".join("?" * len(sample_names))
        rows = con.execute(
            f"SELECT sample_id, sample_name FROM samples WHERE sample_name IN ({placeholders})",
            sample_names,
        ).fetchall()

        found_names = {r[1] for r in rows}
        missing = [n for n in sample_names if n not in found_names]
        if missing:
            raise UpdateError(f"Sample(s) not found: {', '.join(missing)}")

        removal_ids = BitMap([r[0] for r in rows])
        id_list = list(removal_ids)

        # 2-3. Clear bits from all Parquet files (flat and partitioned)
        variants_dir = os.path.join(db_dir, "variants")
        if os.path.exists(variants_dir):
            for pq_file in sorted(
                glob_module.glob(os.path.join(variants_dir, "*.parquet"))
            ):
                _clear_bits_from_parquet(pq_file, removal_ids)
            # Also handle partitioned format (variants/{chrom}/bucket_*.parquet)
            for chrom_dir in sorted(Path(variants_dir).iterdir()):
                if chrom_dir.is_dir():
                    for pq_file in sorted(
                        glob_module.glob(str(chrom_dir / "bucket_*.parquet"))
                    ):
                        _clear_bits_from_parquet(pq_file, removal_ids)

        # 4-5. Delete from SQLite
        ph2 = ",".join("?" * len(id_list))
        con.execute(f"DELETE FROM sample_icd10 WHERE sample_id IN ({ph2})", id_list)
        con.execute(f"DELETE FROM samples WHERE sample_id IN ({ph2})", id_list)

        # 6. Regenerate precomputed bitmaps
        _regenerate_precomputed_bitmaps(con)
        con.commit()

        # 7. Get new count and preserve next_sample_id for future adds
        new_count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        max_id_row = con.execute("SELECT MAX(sample_id) FROM samples").fetchone()
        max_id = max_id_row[0] if max_id_row[0] is not None else -1
        # next_sample_id should always be at least max_id + 1 (from manifest or DB)
        next_id = max(manifest.get("next_sample_id", max_id + 1), max_id + 1)

    finally:
        con.close()

    # 8. Update manifest (preserve next_sample_id to prevent reuse after removal)
    _update_manifest(db_dir, new_count, next_sample_id=next_id)

    return {"removed": list(sample_names), "not_found": []}


def check_database(db_dir: str) -> list[CheckResult]:
    """Validate database integrity. Returns list of CheckResult."""
    results: list[CheckResult] = []

    def err(msg: str) -> None:
        results.append(CheckResult("error", msg))

    def warn(msg: str) -> None:
        results.append(CheckResult("warning", msg))

    def info(msg: str) -> None:
        results.append(CheckResult("info", msg))

    # Check 1: manifest.json exists
    manifest_path = os.path.join(db_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        err("manifest.json not found")
        return results

    # Check 2: manifest.json parses as JSON
    try:
        with open(manifest_path) as f:
            manifest = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        err(f"manifest.json parse error: {e}")
        return results

    # Check 3: required keys present
    required_keys = {"genome_build", "sample_count"}
    missing_keys = required_keys - set(manifest.keys())
    if missing_keys:
        err(f"manifest.json missing keys: {', '.join(sorted(missing_keys))}")
        return results

    # Check 4: genome_build is valid
    if manifest["genome_build"] not in VALID_GENOME_BUILDS:
        err(f"Invalid genome_build: '{manifest['genome_build']}'")

    # Check 5: metadata.sqlite exists
    sqlite_path = os.path.join(db_dir, "metadata.sqlite")
    if not os.path.exists(sqlite_path):
        err("metadata.sqlite not found")
        return results

    # Check 6: all 4 tables exist in SQLite
    required_tables = {"samples", "technologies", "sample_icd10", "precomputed_bitmaps"}
    try:
        con = sqlite3.connect(sqlite_path)
        existing_tables = {
            r[0]
            for r in con.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
        }
        missing_tables = required_tables - existing_tables
        if missing_tables:
            err(f"SQLite missing tables: {', '.join(sorted(missing_tables))}")
            con.close()
            return results
    except Exception as e:
        err(f"SQLite open error: {e}")
        return results

    n_techs = 0
    max_sample_id = -1

    try:
        # Check 7: sample_count matches COUNT(*) FROM samples
        db_count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        if manifest["sample_count"] != db_count:
            warn(
                f"sample_count mismatch: manifest says {manifest['sample_count']}, "
                f"database has {db_count}"
            )

        # Check 8: capture pickle files exist for each tech_id
        tech_ids = [
            r[0] for r in con.execute("SELECT tech_id FROM technologies").fetchall()
        ]
        n_techs = len(tech_ids)
        capture_dir = os.path.join(db_dir, "capture")
        for tech_id in tech_ids:
            pickle_path = os.path.join(capture_dir, f"tech_{tech_id}.pickle")
            if not os.path.exists(pickle_path):
                err(f"Missing capture index: tech_{tech_id}.pickle")

        max_sid_row = con.execute("SELECT MAX(sample_id) FROM samples").fetchone()
        max_sample_id = max_sid_row[0] if max_sid_row[0] is not None else -1

    finally:
        con.close()

    # Check 9: variants/ directory exists
    variants_dir = os.path.join(db_dir, "variants")
    if not os.path.exists(variants_dir):
        err("variants/ directory not found")
        return results

    # Collect flat parquet files
    flat_parquets = sorted(glob_module.glob(os.path.join(variants_dir, "*.parquet")))
    # Collect partitioned parquets (variants/{chrom}/bucket_*.parquet)
    chrom_dirs = []
    bucket_parquets = []
    for entry in sorted(os.scandir(variants_dir), key=lambda e: e.name):
        if entry.is_dir():
            chrom_dirs.append(entry.path)
            bucket_parquets.extend(
                sorted(glob_module.glob(os.path.join(entry.path, "bucket_*.parquet")))
            )

    parquet_files = flat_parquets + bucket_parquets
    n_chroms = len(flat_parquets) + len(chrom_dirs)

    expected_fields = {
        "pos":        pa.uint32(),
        "ref":        pa.large_utf8(),
        "alt":        pa.large_utf8(),
        "het_bitmap": pa.large_binary(),
        "hom_bitmap": pa.large_binary(),
    }

    for pq_path in parquet_files:
        basename = os.path.basename(pq_path)
        if basename.startswith("bucket_"):
            chrom = os.path.basename(os.path.dirname(pq_path))
            label = f"{chrom}/{basename}"
        else:
            chrom = basename.replace(".parquet", "")
            label = f"{chrom}.parquet"

        # Check 10: each Parquet file readable
        try:
            table = pq.read_table(pq_path)
        except Exception as e:
            err(f"{label}: read error: {e}")
            continue

        # Check 11: correct schema
        for fname, ftype in expected_fields.items():
            if fname not in table.schema.names:
                err(f"{label}: missing field '{fname}'")
            elif table.schema.field(fname).type != ftype:
                err(
                    f"{label}: field '{fname}' has wrong type "
                    f"'{table.schema.field(fname).type}' (expected '{ftype}')"
                )

        # Check 12: rows sorted by (pos, alt)
        if len(table) > 1:
            pos_col = table["pos"].to_pylist()
            alt_col = table["alt"].to_pylist()
            pairs = list(zip(pos_col, alt_col))
            if pairs != sorted(pairs):
                warn(f"{label}: rows not sorted by (pos, alt)")

        # Check 13: bitmap max element ≤ max(sample_id) (sample first 1000 rows)
        sample_limit = min(1000, len(table))
        for i in range(sample_limit):
            for field in ("het_bitmap", "hom_bitmap"):
                bm_bytes = table[field][i].as_py()
                if bm_bytes:
                    bm = deserialize(bm_bytes)
                    if bm and max(bm) > max_sample_id:
                        err(
                            f"{label} row {i}: {field} contains sample_id "
                            f"{max(bm)} > max sample_id {max_sample_id}"
                        )

    # Check 14: schema_version key absent → warning
    if "schema_version" not in manifest:
        warn("manifest.json missing 'schema_version' key")

    # Check 15: info summary
    info(
        f"genome_build={manifest['genome_build']}, "
        f"sample_count={manifest['sample_count']}, "
        f"chroms={n_chroms}, "
        f"techs={n_techs}"
    )

    return results
