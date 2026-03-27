#!/usr/bin/env python3
"""Build AFQuery databases: 1 WGS ground truth + 3 WES scenarios.

Inputs:
  - manifests/manifest_wgs.tsv
  - manifests/manifest_{balanced,skewed,extreme}.tsv
  - beds/afquery/*.bed  (for WES databases)

Outputs:
  - dbs/db_wgs/         ground truth (all samples as WGS)
  - dbs/db_balanced/    technology-aware, balanced kit assignment
  - dbs/db_skewed/      technology-aware, realistic skew
  - dbs/db_extreme/     technology-aware, extreme skew
"""

import logging
import shutil
import sys
import time
from pathlib import Path

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from config import (
    AFQUERY_BED_DIR,
    DB_DIR,
    GENOME_BUILD,
    MANIFEST_DIR,
    SCENARIOS,
    THREADS,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def build_db(manifest_path: Path, db_dir: Path, bed_dir: Path | None = None):
    """Build one AFQuery database using run_preprocess()."""
    from afquery.preprocess import run_preprocess

    if (db_dir / "manifest.json").exists():
        logger.info("DB already exists: %s", db_dir)
        return

    if db_dir.exists():
        shutil.rmtree(db_dir)

    t0 = time.perf_counter()
    run_preprocess(
        manifest_path=str(manifest_path),
        output_dir=str(db_dir),
        genome_build=GENOME_BUILD,
        bed_dir=str(bed_dir) if bed_dir else None,
        threads=THREADS,
    )
    elapsed = time.perf_counter() - t0

    # Validate that the build produced variant data
    variants_dir = db_dir / "variants"
    parquet_files = list(variants_dir.rglob("*.parquet"))
    if not parquet_files:
        raise RuntimeError(
            f"Build of {db_dir.name} produced no variant parquet files. "
            "Check that VCF chromosomes match the genome build and BED file naming."
        )

    logger.info("Built %s in %.1fs (%d parquet file(s))", db_dir.name, elapsed, len(parquet_files))


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--scenario", choices=["wgs"] + list(SCENARIOS.keys()))
    args = parser.parse_args()

    ensure_dirs()

    if args.scenario:
        # Single-scenario CLI mode
        scenario = args.scenario
        if scenario == "wgs":
            logger.info("=== Building WGS ground truth database ===")
            build_db(MANIFEST_DIR / "manifest_wgs.tsv", DB_DIR / "db_wgs")
        else:
            logger.info("=== Building WES database: %s ===", scenario)
            build_db(
                MANIFEST_DIR / f"manifest_{scenario}.tsv",
                DB_DIR / f"db_{scenario}",
                bed_dir=AFQUERY_BED_DIR,
            )
        return

    # Full-sweep mode (backward compatible)
    logger.info("=== Building WGS ground truth database ===")
    build_db(MANIFEST_DIR / "manifest_wgs.tsv", DB_DIR / "db_wgs")

    for scenario_name in SCENARIOS:
        logger.info("=== Building WES database: %s ===", scenario_name)
        build_db(
            MANIFEST_DIR / f"manifest_{scenario_name}.tsv",
            DB_DIR / f"db_{scenario_name}",
            bed_dir=AFQUERY_BED_DIR,
        )

    logger.info("All databases built.")


if __name__ == "__main__":
    main()
