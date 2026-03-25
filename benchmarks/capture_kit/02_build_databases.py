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

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from capture_kit_bench.config import (
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
    logger.info("Built %s in %.1fs", db_dir.name, elapsed)


def main():
    ensure_dirs()

    # DB1: WGS ground truth (no BED files)
    logger.info("=== Building WGS ground truth database ===")
    build_db(MANIFEST_DIR / "manifest_wgs.tsv", DB_DIR / "db_wgs")

    # DB2: WES databases for each scenario (with BED files)
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
