#!/usr/bin/env python3
"""Experiment 3: VCF annotation throughput.

Measures variants/second when annotating VCFs of different sizes against the
1KG chr22 database (2504 samples), varying thread count.

Output: results/annotate_throughput.json
"""

import json
import logging
import os
import random
import sys
import tempfile
import time
from pathlib import Path

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from config import (
    ANNOTATE_THREAD_COUNTS,
    ANNOTATE_VARIANT_COUNTS,
    ONEKG_DB_DIR,
    ONEKG_SUBSETS,
    RESULTS_DIR,
    SEED,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def _create_input_vcf(db_path: Path, n_variants: int, seed: int, output_path: Path) -> int:
    """Create a VCF with n_variants sampled from the database.

    Returns actual number of variants written.
    """
    from afquery.benchmark import _find_test_variants

    all_variants = _find_test_variants(db_path, n=n_variants + 500)
    if len(all_variants) < n_variants:
        logger.warning(
            "Only %d variants available, requested %d",
            len(all_variants),
            n_variants,
        )

    rng = random.Random(seed)
    selected = rng.sample(all_variants, min(n_variants, len(all_variants)))
    # Sort by (chrom, pos) for VCF order
    selected.sort(key=lambda v: (v[0], v[1]))

    chroms = sorted(set(v[0] for v in selected))

    with open(output_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        for chrom in chroms:
            f.write(f"##contig=<ID={chrom}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for chrom, pos, ref, alt in selected:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/1\n")

    return len(selected)


def _run_one_annotate(db_path: Path, n_variants: int, threads: int, output_path: Path):
    """Run one (n_variants, threads) annotation benchmark and write a per-run JSON."""
    from afquery.database import Database

    with tempfile.TemporaryDirectory(prefix="afquery_annot_bench_") as tmpdir:
        tmpdir = Path(tmpdir)
        input_vcf = tmpdir / f"input_{n_variants}.vcf"
        actual_n = _create_input_vcf(db_path, n_variants, SEED, input_vcf)
        logger.info("Created input VCF with %d variants", actual_n)

        output_vcf = tmpdir / f"output_{n_variants}_{threads}t.vcf"
        logger.info("=== Annotate: %d variants, %d threads ===", actual_n, threads)

        db = Database(str(db_path))
        t0 = time.perf_counter()
        db.annotate_vcf(
            input_vcf=str(input_vcf),
            output_vcf=str(output_vcf),
            n_workers=threads,
        )
        elapsed_s = time.perf_counter() - t0

    throughput = actual_n / elapsed_s if elapsed_s > 0 else 0
    result = {
        "n_variants": actual_n,
        "threads": threads,
        "wall_s": round(elapsed_s, 3),
        "throughput_vars_per_s": round(throughput, 1),
    }
    logger.info("  Result: %s", result)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2))
    logger.info("Result saved to %s", output_path)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--n-variants", type=int)
    parser.add_argument("--threads", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    ensure_dirs()

    max_subset = max(ONEKG_SUBSETS)
    db_path = ONEKG_DB_DIR / f"1kg_{max_subset}"
    if not (db_path / "manifest.json").exists():
        logger.error("1KG DB not found at %s. Run 01_prepare_data.py first.", db_path)
        return

    if args.n_variants is not None and args.threads is not None and args.output:
        # Single-run CLI mode — no speedup computation (done in collect_annotate.py)
        _run_one_annotate(db_path, args.n_variants, args.threads, args.output)
        return

    # Full-sweep mode (backward compatible)
    all_results = []

    with tempfile.TemporaryDirectory(prefix="afquery_annot_bench_") as tmpdir:
        tmpdir = Path(tmpdir)

        for n_variants in ANNOTATE_VARIANT_COUNTS:
            input_vcf = tmpdir / f"input_{n_variants}.vcf"
            actual_n = _create_input_vcf(db_path, n_variants, SEED, input_vcf)
            logger.info("Created input VCF with %d variants", actual_n)

            for threads in ANNOTATE_THREAD_COUNTS:
                output_vcf = tmpdir / f"output_{n_variants}_{threads}t.vcf"
                logger.info("=== Annotate: %d variants, %d threads ===", actual_n, threads)

                from afquery.database import Database

                db = Database(str(db_path))
                t0 = time.perf_counter()
                db.annotate_vcf(
                    input_vcf=str(input_vcf),
                    output_vcf=str(output_vcf),
                    n_workers=threads,
                )
                elapsed_s = time.perf_counter() - t0

                throughput = actual_n / elapsed_s if elapsed_s > 0 else 0
                result = {
                    "n_variants": actual_n,
                    "threads": threads,
                    "wall_s": round(elapsed_s, 3),
                    "throughput_vars_per_s": round(throughput, 1),
                }
                all_results.append(result)
                logger.info("  Result: %s", result)

                if output_vcf.exists():
                    output_vcf.unlink()

    # Compute speedup ratios
    for n_variants in ANNOTATE_VARIANT_COUNTS:
        baseline = next(
            (r for r in all_results if r["n_variants"] == n_variants and r["threads"] == 1),
            None,
        )
        if baseline:
            for r in all_results:
                if r["n_variants"] == n_variants:
                    r["speedup"] = (
                        round(baseline["wall_s"] / r["wall_s"], 2) if r["wall_s"] > 0 else 0
                    )

    out = RESULTS_DIR / "annotate_throughput.json"
    out.write_text(json.dumps(all_results, indent=2))
    logger.info("Results saved to %s", out)


if __name__ == "__main__":
    main()
