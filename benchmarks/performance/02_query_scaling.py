#!/usr/bin/env python3
"""Experiment 1: Query latency scaling with sample count.

Measures point, region, and batch query latencies across synthetic databases
of increasing sample count (1K to 50K). Tests with and without filters.

Output: results/query_scaling.json
"""

import json
import logging
import sys
from pathlib import Path

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from shared.utils import stats as _stats, time_ms as _time_ms  # noqa: E402
from config import (  # noqa: E402
    QUERY_BATCH_REPS,
    QUERY_COLD_REPS,
    QUERY_REGION_REPS,
    QUERY_WARM_REPS,
    QUERY_WARMUP,
    RESULTS_DIR,
    SYNTH_SCALES,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def _load_inventory() -> dict:
    inv_path = RESULTS_DIR / "data_inventory.json"
    if not inv_path.exists():
        raise FileNotFoundError("Run 01_prepare_data.py first.")
    return json.loads(inv_path.read_text())


def _find_variants(db_path: Path, n: int = 1100):
    """Find test variants from a database."""
    from afquery.benchmark import _find_test_variants

    return _find_test_variants(db_path, n=n)


# Filter configurations to test
FILTER_CONFIGS = {
    "no_filter": {"phenotype": None, "sex": "both", "tech": None},
    "sex_female": {"phenotype": None, "sex": "female", "tech": None},
    "phenotype": {"phenotype": ["E11.9"], "sex": "both", "tech": None},
    "combined": {"phenotype": ["E11.9"], "sex": "female", "tech": None},
}


def bench_one_db(db_path: str, n_samples: int) -> dict:
    """Run all query benchmarks on a single database."""
    from afquery.database import Database

    db = Database(db_path)
    variants = _find_variants(Path(db_path))
    if not variants:
        logger.warning("No variants found in %s", db_path)
        return {"error": "no variants"}

    test_chrom = variants[0][0]
    test_pos = variants[0][1]
    same_chrom = [(p, r, a) for c, p, r, a in variants if c == test_chrom]
    batch_100 = same_chrom[:100]
    batch_1000 = same_chrom[:1000]

    # Find a region with variants (1 Mbp around first variant)
    region_start = max(1, test_pos - 500_000)
    region_end = test_pos + 500_000

    results = {"n_samples": n_samples, "chrom": test_chrom, "filters": {}}

    for filter_name, filt in FILTER_CONFIGS.items():
        logger.info("  Filter: %s", filter_name)
        filt_results = {}

        # --- Point query (cold) ---
        cold_times = []
        for _ in range(QUERY_COLD_REPS):
            # Recreate Database to reset any internal caches
            db_fresh = Database(db_path)
            _, ms = _time_ms(db_fresh.query, chrom=test_chrom, pos=test_pos, **filt)
            cold_times.append(ms)
        filt_results["point_cold"] = _stats(cold_times)

        # --- Point query (warm) ---
        # Warmup
        for _ in range(QUERY_WARMUP):
            db.query(chrom=test_chrom, pos=test_pos, **filt)
        warm_times = []
        for _ in range(QUERY_WARM_REPS):
            _, ms = _time_ms(db.query, chrom=test_chrom, pos=test_pos, **filt)
            warm_times.append(ms)
        filt_results["point_warm"] = _stats(warm_times)

        # --- Region query (1 Mbp) ---
        # Warmup
        for _ in range(QUERY_WARMUP):
            db.query_region(chrom=test_chrom, start=region_start, end=region_end, **filt)
        region_times = []
        for _ in range(QUERY_REGION_REPS):
            _, ms = _time_ms(
                db.query_region,
                chrom=test_chrom,
                start=region_start,
                end=region_end,
                **filt,
            )
            region_times.append(ms)
        filt_results["region_1mbp"] = _stats(region_times)

        # --- Batch 100 ---
        for _ in range(QUERY_WARMUP):
            db.query_batch(chrom=test_chrom, variants=batch_100, **filt)
        batch100_times = []
        for _ in range(QUERY_BATCH_REPS):
            _, ms = _time_ms(db.query_batch, chrom=test_chrom, variants=batch_100, **filt)
            batch100_times.append(ms)
        filt_results["batch_100"] = _stats(batch100_times)

        # --- Batch 1000 ---
        for _ in range(QUERY_WARMUP):
            db.query_batch(chrom=test_chrom, variants=batch_1000, **filt)
        batch1000_times = []
        for _ in range(QUERY_BATCH_REPS):
            _, ms = _time_ms(
                db.query_batch, chrom=test_chrom, variants=batch_1000, **filt
            )
            batch1000_times.append(ms)
        filt_results["batch_1000"] = _stats(batch1000_times)

        results["filters"][filter_name] = filt_results

    return results


def main():
    import argparse

    from config import SYNTH_DB_DIR  # noqa: E402

    parser = argparse.ArgumentParser()
    parser.add_argument("--scale", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    ensure_dirs()

    if args.scale is not None and args.output:
        # Single-scale CLI mode — DB path is deterministic, no inventory needed
        db_path = str(SYNTH_DB_DIR / f"synth_{args.scale}")
        logger.info("=== Benchmarking synth DB: %d samples ===", args.scale)
        result = bench_one_db(db_path, args.scale)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2))
        logger.info("Result saved to %s", args.output)
        return

    # Full-sweep mode (backward compatible)
    inventory = _load_inventory()

    all_results = []
    for n in SYNTH_SCALES:
        db_path = inventory["synth"].get(str(n))
        if not db_path:
            logger.warning("No synth DB for %d samples, skipping", n)
            continue

        logger.info("=== Benchmarking synth DB: %d samples ===", n)
        result = bench_one_db(db_path, n)
        all_results.append(result)

    out = RESULTS_DIR / "query_scaling.json"
    out.write_text(json.dumps(all_results, indent=2))
    logger.info("Results saved to %s", out)


if __name__ == "__main__":
    main()
