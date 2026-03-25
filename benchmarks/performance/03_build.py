#!/usr/bin/env python3
"""Experiment 2: Database build performance.

Measures build time (ingest + consolidate + build phases) across different
sample counts and thread counts. Also records disk usage.

Output: results/build_perf.json
"""

import json
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from config import (
    BUILD_SCALES,
    BUILD_THREAD_COUNTS,
    DATA_DIR,
    GENOME_BUILD,
    RESULTS_DIR,
    SEED,
    SYNTH_CHROMS,
    SYNTH_VARIANTS_PER_CHROM,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def _dir_size_bytes(path: Path) -> int:
    """Recursively compute total size of a directory in bytes."""
    total = 0
    for entry in path.rglob("*"):
        if entry.is_file():
            total += entry.stat().st_size
    return total


def _run_build_timed(manifest_path: str, db_dir: str, threads: int) -> dict:
    """Run afquery create-db via subprocess with /usr/bin/time for peak RSS.

    Returns timing and memory metrics.
    """
    if os.path.exists(db_dir):
        shutil.rmtree(db_dir)

    cmd = [
        "/usr/bin/time", "-v",
        sys.executable, "-m", "afquery",
        "create-db",
        "--manifest", manifest_path,
        "--output-dir", db_dir,
        "--genome-build", GENOME_BUILD,
        "--threads", str(threads),
        "--build-threads", str(threads),
    ]

    t0 = time.perf_counter()
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(Path(__file__).resolve().parents[2])
    )
    wall_s = time.perf_counter() - t0

    if result.returncode != 0:
        logger.error("Build failed:\nstdout: %s\nstderr: %s", result.stdout, result.stderr)
        return {"error": result.stderr, "wall_s": wall_s}

    # Parse peak RSS from /usr/bin/time -v output
    peak_rss_kb = 0
    for line in result.stderr.split("\n"):
        if "Maximum resident set size" in line:
            peak_rss_kb = int(line.strip().split()[-1])
            break

    db_size_bytes = _dir_size_bytes(Path(db_dir))

    return {
        "wall_s": round(wall_s, 2),
        "peak_rss_mb": round(peak_rss_kb / 1024, 1),
        "db_size_mb": round(db_size_bytes / (1024 * 1024), 1),
    }


def main():
    ensure_dirs()

    from afquery.preprocess.synth import generate_synthetic_manifest

    build_dir = DATA_DIR / "build_bench"
    all_results = []

    for n_samples in BUILD_SCALES:
        # Generate synthetic data once per scale
        synth_out = build_dir / f"synth_{n_samples}"
        manifest_path = synth_out / "manifest.tsv"

        if not manifest_path.exists():
            logger.info("Generating synthetic data for %d samples...", n_samples)
            if synth_out.exists():
                shutil.rmtree(synth_out)
            generate_synthetic_manifest(
                output_dir=synth_out,
                n_samples=n_samples,
                n_variants_per_chrom=SYNTH_VARIANTS_PER_CHROM,
                chroms=SYNTH_CHROMS,
                seed=SEED,
            )

        # Compute raw VCF size
        vcf_dir = synth_out / "vcfs"
        raw_vcf_size_mb = round(_dir_size_bytes(vcf_dir) / (1024 * 1024), 1)

        for threads in BUILD_THREAD_COUNTS:
            db_dir = build_dir / f"db_{n_samples}_{threads}t"

            logger.info(
                "=== Build: %d samples, %d threads ===", n_samples, threads
            )

            metrics = _run_build_timed(str(manifest_path), str(db_dir), threads)
            metrics["n_samples"] = n_samples
            metrics["threads"] = threads
            metrics["raw_vcf_size_mb"] = raw_vcf_size_mb

            all_results.append(metrics)
            logger.info("  Result: %s", metrics)

            # Clean up DB to save disk space (keep only the last one per scale)
            if threads != BUILD_THREAD_COUNTS[-1]:
                if db_dir.exists():
                    shutil.rmtree(db_dir)

    out = RESULTS_DIR / "build_perf.json"
    out.write_text(json.dumps(all_results, indent=2))
    logger.info("Results saved to %s", out)


if __name__ == "__main__":
    main()
