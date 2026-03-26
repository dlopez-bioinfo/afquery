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

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

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


def _build_one_1kg(n: int, output_path: Path):
    """Build one 1KG database and write a per-run JSON."""
    from afquery.preprocess import run_preprocess

    manifest_path = ONEKG_DIR / "manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError("1KG manifest not found. Run download_1kg first.")

    db_dir = ONEKG_DB_DIR / f"1kg_{n}"
    if (db_dir / "manifest.json").exists():
        logger.info("1KG DB already exists for %d samples: %s", n, db_dir)
    else:
        import shutil
        if db_dir.exists():
            shutil.rmtree(db_dir)
        sub_manifest = ONEKG_DIR / f"manifest_{n}.tsv"
        _subsample_manifest(manifest_path, n, SEED, sub_manifest)
        logger.info("Building 1KG DB with %d samples...", n)
        t0 = time.perf_counter()
        run_preprocess(
            manifest_path=str(sub_manifest),
            output_dir=str(db_dir),
            genome_build=GENOME_BUILD,
        )
        elapsed = time.perf_counter() - t0
        logger.info("  Built %d-sample DB in %.1fs", n, elapsed)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps({"mode": "1kg", "n_samples": n, "db_path": str(db_dir)}, indent=2)
    )


def _build_one_synth(n: int, output_path: Path):
    """Build one synthetic database and write a per-run JSON."""
    from afquery.preprocess import run_preprocess
    from afquery.preprocess.synth import generate_synthetic_manifest

    db_dir = SYNTH_DB_DIR / f"synth_{n}"
    if (db_dir / "manifest.json").exists():
        logger.info("Synth DB already exists for %d samples: %s", n, db_dir)
    else:
        import shutil
        synth_out = SYNTH_DIR / f"synth_{n}"
        if synth_out.exists():
            shutil.rmtree(synth_out)
        logger.info("Generating synthetic data for %d samples...", n)
        manifest_path = generate_synthetic_manifest(
            output_dir=synth_out,
            n_samples=n,
            n_variants_per_chrom=SYNTH_VARIANTS_PER_CHROM,
            chroms=SYNTH_CHROMS,
            seed=SEED,
        )
        if db_dir.exists():
            shutil.rmtree(db_dir)
        logger.info("Building synthetic DB for %d samples...", n)
        t0 = time.perf_counter()
        run_preprocess(
            manifest_path=str(manifest_path),
            output_dir=str(db_dir),
            genome_build=GENOME_BUILD,
        )
        elapsed = time.perf_counter() - t0
        logger.info("  Built %d-sample synth DB in %.1fs", n, elapsed)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps({"mode": "synth", "n_samples": n, "db_path": str(db_dir)}, indent=2)
    )


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["1kg", "synth"])
    parser.add_argument("--n-samples", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    ensure_dirs()

    if args.mode and args.n_samples is not None and args.output:
        # Single-item CLI mode
        if args.mode == "1kg":
            _build_one_1kg(args.n_samples, args.output)
        else:
            _build_one_synth(args.n_samples, args.output)
        return

    # Full-sweep mode (backward compatible)
    inventory = {}

    logger.info("=== Building 1KG databases ===")
    inventory["1kg"] = build_1kg_databases()

    logger.info("=== Building synthetic databases ===")
    inventory["synth"] = build_synthetic_databases()

    out = RESULTS_DIR / "data_inventory.json"
    out.write_text(json.dumps(inventory, indent=2))
    logger.info("Data inventory saved to %s", out)


if __name__ == "__main__":
    main()
