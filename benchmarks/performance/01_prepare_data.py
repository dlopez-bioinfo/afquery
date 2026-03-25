#!/usr/bin/env python3
"""Build AFQuery databases for benchmarking.

Creates:
  - 1KG databases at subset sizes (500, 1000, 2504 samples)
  - Synthetic databases at scales (1K, 5K, 10K, 25K, 50K samples)

Reads the 1KG manifest from 00_download_1kg.sh output and subsamples it.
"""

import json
import logging
import random
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from config import (
    DATA_DIR,
    GENOME_BUILD,
    ONEKG_DIR,
    ONEKG_DB_DIR,
    ONEKG_SUBSETS,
    RESULTS_DIR,
    SEED,
    SYNTH_CHROMS,
    SYNTH_DB_DIR,
    SYNTH_DIR,
    SYNTH_SCALES,
    SYNTH_VARIANTS_PER_CHROM,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def _subsample_manifest(
    manifest_path: Path, n_samples: int, seed: int, output_path: Path
) -> Path:
    """Create a manifest with a random subset of samples."""
    lines = manifest_path.read_text().strip().split("\n")
    header, rows = lines[0], lines[1:]

    rng = random.Random(seed)
    if n_samples >= len(rows):
        subset = rows
    else:
        subset = rng.sample(rows, n_samples)

    output_path.write_text(header + "\n" + "\n".join(subset) + "\n")
    return output_path


def build_1kg_databases() -> dict:
    """Build AFQuery databases from 1KG manifest at different subset sizes."""
    from afquery.preprocess import run_preprocess

    manifest_path = ONEKG_DIR / "manifest.tsv"
    if not manifest_path.exists():
        logger.error("1KG manifest not found. Run 00_download_1kg.sh first.")
        return {}

    results = {}
    for n in ONEKG_SUBSETS:
        db_dir = ONEKG_DB_DIR / f"1kg_{n}"
        if (db_dir / "manifest.json").exists():
            logger.info("1KG DB already exists for %d samples: %s", n, db_dir)
            results[n] = str(db_dir)
            continue

        logger.info("Building 1KG DB with %d samples...", n)
        if db_dir.exists():
            shutil.rmtree(db_dir)

        sub_manifest = ONEKG_DIR / f"manifest_{n}.tsv"
        _subsample_manifest(manifest_path, n, SEED, sub_manifest)

        t0 = time.perf_counter()
        run_preprocess(
            manifest_path=str(sub_manifest),
            output_dir=str(db_dir),
            genome_build=GENOME_BUILD,
        )
        elapsed = time.perf_counter() - t0
        logger.info("  Built %d-sample DB in %.1fs", n, elapsed)
        results[n] = str(db_dir)

    return results


def build_synthetic_databases() -> dict:
    """Build synthetic AFQuery databases at multiple scales."""
    from afquery.preprocess import run_preprocess
    from afquery.preprocess.synth import generate_synthetic_manifest

    results = {}
    for n in SYNTH_SCALES:
        db_dir = SYNTH_DB_DIR / f"synth_{n}"
        if (db_dir / "manifest.json").exists():
            logger.info("Synth DB already exists for %d samples: %s", n, db_dir)
            results[n] = str(db_dir)
            continue

        logger.info("Generating synthetic data for %d samples...", n)
        synth_out = SYNTH_DIR / f"synth_{n}"
        if synth_out.exists():
            shutil.rmtree(synth_out)

        manifest_path = generate_synthetic_manifest(
            output_dir=synth_out,
            n_samples=n,
            n_variants_per_chrom=SYNTH_VARIANTS_PER_CHROM,
            chroms=SYNTH_CHROMS,
            seed=SEED,
        )

        logger.info("Building synthetic DB for %d samples...", n)
        if db_dir.exists():
            shutil.rmtree(db_dir)

        t0 = time.perf_counter()
        run_preprocess(
            manifest_path=str(manifest_path),
            output_dir=str(db_dir),
            genome_build=GENOME_BUILD,
        )
        elapsed = time.perf_counter() - t0
        logger.info("  Built %d-sample synth DB in %.1fs", n, elapsed)
        results[n] = str(db_dir)

    return results


def main():
    ensure_dirs()

    inventory = {}

    logger.info("=== Building 1KG databases ===")
    inventory["1kg"] = build_1kg_databases()

    logger.info("=== Building synthetic databases ===")
    inventory["synth"] = build_synthetic_databases()

    # Save inventory for other scripts
    out = RESULTS_DIR / "data_inventory.json"
    out.write_text(json.dumps(inventory, indent=2))
    logger.info("Data inventory saved to %s", out)


if __name__ == "__main__":
    main()
