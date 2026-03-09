import json
import logging
import os
import sqlite3
import tempfile
from datetime import datetime, timezone

from ..bitmaps import serialize, build_sex_bitmaps, build_phenotype_bitmaps, build_tech_bitmaps
from ..constants import VALID_GENOME_BUILDS
from ..models import Sample, Technology
from .build import build_all_parquets
from .ingest import ingest_all
from .manifest import parse_manifest
from .regions import build_capture_indices
from .update import add_samples, remove_samples, check_database, UpdateError, CheckResult

logger = logging.getLogger(__name__)


def run_preprocess(
    manifest_path: str,
    output_dir: str,
    genome_build: str,
    bed_dir: str | None = None,
    n_threads: int = 8,
    tmp_dir: str | None = None,
    n_workers: int | None = None,
) -> None:
    if genome_build not in VALID_GENOME_BUILDS:
        raise ValueError(
            f"Invalid genome_build '{genome_build}', must be one of {sorted(VALID_GENOME_BUILDS)}"
        )

    samples_raw, techs_raw = parse_manifest(manifest_path, bed_dir)

    os.makedirs(os.path.join(output_dir, "variants"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "capture"), exist_ok=True)

    # Assign tech_ids (0-indexed, first-seen order) and sample_ids (0-indexed, manifest order)
    tech_name_to_id: dict[str, int] = {}
    for idx, pt in enumerate(techs_raw):
        tech_name_to_id[pt.tech_name] = idx

    technologies = [
        Technology(
            tech_id=tech_name_to_id[pt.tech_name],
            tech_name=pt.tech_name,
            bed_path=pt.bed_path,
        )
        for pt in techs_raw
    ]

    samples = [
        Sample(
            sample_id=idx,
            sample_name=ps.sample_name,
            sex=ps.sex,
            tech_id=tech_name_to_id[ps.tech_name],
        )
        for idx, ps in enumerate(samples_raw)
    ]

    vcf_paths = [ps.vcf_path for ps in samples_raw]

    sample_phenotype_pairs: list[tuple[int, str]] = [
        (idx, code)
        for idx, ps in enumerate(samples_raw)
        for code in ps.phenotype_codes
    ]

    _write_sqlite(output_dir, samples, technologies, sample_phenotype_pairs)

    capture_dir = os.path.join(output_dir, "capture")
    build_capture_indices(technologies, capture_dir)

    auto_tmp = tmp_dir is None
    if auto_tmp:
        tmp_dir = tempfile.mkdtemp(prefix="afquery_preprocess_")

    try:
        ingest_all(samples, vcf_paths, tmp_dir, n_workers=n_threads)
        variants_dir = os.path.join(output_dir, "variants")
        build_all_parquets(tmp_dir, variants_dir, n_workers=n_workers)
    finally:
        if auto_tmp:
            import shutil
            shutil.rmtree(tmp_dir, ignore_errors=True)

    _write_manifest(output_dir, genome_build, len(samples))


def _write_sqlite(
    output_dir: str,
    samples: list[Sample],
    technologies: list[Technology],
    sample_phenotype_pairs: list[tuple[int, str]],
) -> None:
    db_path = os.path.join(output_dir, "metadata.sqlite")
    con = sqlite3.connect(db_path)

    con.executescript("""
        CREATE TABLE samples (
            sample_id   INTEGER PRIMARY KEY,
            sample_name TEXT NOT NULL,
            sex         TEXT NOT NULL,
            tech_id     INTEGER NOT NULL
        );
        CREATE TABLE technologies (
            tech_id   INTEGER PRIMARY KEY,
            tech_name TEXT NOT NULL,
            bed_path  TEXT
        );
        CREATE TABLE sample_phenotype (
            sample_id    INTEGER NOT NULL,
            phenotype_code   TEXT NOT NULL,
            PRIMARY KEY (sample_id, phenotype_code)
        );
        CREATE TABLE precomputed_bitmaps (
            bitmap_type TEXT NOT NULL,
            bitmap_key  TEXT NOT NULL,
            bitmap_data BLOB NOT NULL,
            PRIMARY KEY (bitmap_type, bitmap_key)
        );
    """)

    con.executemany(
        "INSERT INTO samples VALUES (?, ?, ?, ?)",
        [(s.sample_id, s.sample_name, s.sex, s.tech_id) for s in samples],
    )
    con.executemany(
        "INSERT INTO technologies VALUES (?, ?, ?)",
        [(t.tech_id, t.tech_name, t.bed_path) for t in technologies],
    )
    con.executemany(
        "INSERT INTO sample_phenotype VALUES (?, ?)",
        sample_phenotype_pairs,
    )

    for sex, bm in build_sex_bitmaps(samples).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("sex", sex, serialize(bm)),
        )

    for code, bm in build_phenotype_bitmaps(sample_phenotype_pairs).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("phenotype", code, serialize(bm)),
        )

    for tech_id, bm in build_tech_bitmaps(samples).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("tech", str(tech_id), serialize(bm)),
        )

    con.commit()
    con.close()


def _write_manifest(output_dir: str, genome_build: str, sample_count: int) -> None:
    manifest = {
        "genome_build": genome_build,
        "version": "0.1.0",
        "sample_count": sample_count,
        "schema_version": "1.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    with open(os.path.join(output_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
