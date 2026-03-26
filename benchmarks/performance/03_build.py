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

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

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


def _gen_synth(n_samples: int, output_path: Path):
    """Generate synthetic VCFs for one scale and write a sentinel JSON."""
    from afquery.preprocess.synth import generate_synthetic_manifest

    build_dir = DATA_DIR / "build_bench"
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

    vcf_dir = synth_out / "vcfs"
    raw_vcf_size_mb = round(_dir_size_bytes(vcf_dir) / (1024 * 1024), 1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(
            {
                "phase": "gen",
                "n_samples": n_samples,
                "manifest_path": str(manifest_path),
                "raw_vcf_size_mb": raw_vcf_size_mb,
            },
            indent=2,
        )
    )
    logger.info("Gen sentinel saved to %s", output_path)


def _build_one(n_samples: int, threads: int, output_path: Path):
    """Run one (n_samples, threads) build benchmark and write a per-run JSON."""
    build_dir = DATA_DIR / "build_bench"
    manifest_path = build_dir / f"synth_{n_samples}" / "manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Synth manifest not found: {manifest_path}. Run --phase gen first."
        )

    vcf_dir = build_dir / f"synth_{n_samples}" / "vcfs"
    raw_vcf_size_mb = round(_dir_size_bytes(vcf_dir) / (1024 * 1024), 1)

    db_dir = build_dir / f"db_{n_samples}_{threads}t"
    logger.info("=== Build: %d samples, %d threads ===", n_samples, threads)

    metrics = _run_build_timed(str(manifest_path), str(db_dir), threads)
    metrics["n_samples"] = n_samples
    metrics["threads"] = threads
    metrics["raw_vcf_size_mb"] = raw_vcf_size_mb
    logger.info("  Result: %s", metrics)

    # Always clean up DB after measuring (no "last threads" to keep in parallel mode)
    if db_dir.exists():
        shutil.rmtree(db_dir)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(metrics, indent=2))
    logger.info("Result saved to %s", output_path)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["gen", "build"])
    parser.add_argument("--n-samples", type=int)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    ensure_dirs()

    if args.phase == "gen" and args.n_samples is not None and args.output:
        _gen_synth(args.n_samples, args.output)
        return

    if args.phase == "build" and args.n_samples is not None and args.threads is not None and args.output:
        _build_one(args.n_samples, args.threads, args.output)
        return

    # Full-sweep mode (backward compatible)
    from afquery.preprocess.synth import generate_synthetic_manifest

    build_dir = DATA_DIR / "build_bench"
    all_results = []

    for n_samples in BUILD_SCALES:
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

        vcf_dir = synth_out / "vcfs"
        raw_vcf_size_mb = round(_dir_size_bytes(vcf_dir) / (1024 * 1024), 1)

        for threads in BUILD_THREAD_COUNTS:
            db_dir = build_dir / f"db_{n_samples}_{threads}t"

            logger.info("=== Build: %d samples, %d threads ===", n_samples, threads)

            metrics = _run_build_timed(str(manifest_path), str(db_dir), threads)
            metrics["n_samples"] = n_samples
            metrics["threads"] = threads
            metrics["raw_vcf_size_mb"] = raw_vcf_size_mb

            all_results.append(metrics)
            logger.info("  Result: %s", metrics)

            # Always clean up DB after measuring
            if db_dir.exists():
                shutil.rmtree(db_dir)

    out = RESULTS_DIR / "build_perf.json"
    out.write_text(json.dumps(all_results, indent=2))
    logger.info("Results saved to %s", out)


if __name__ == "__main__":
    main()
