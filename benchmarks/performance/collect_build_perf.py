#!/usr/bin/env python3
"""Aggregate per-run build JSONs into build_perf.json.

When multiple repetitions exist per (n_samples, threads), computes
median and IQR for wall_s, peak_rss_mb, and db_size_mb.

Usage:
    python collect_build_perf.py \
        --inputs results/raw/build/build_1000_1t_r1.json ... \
        --output results/build_perf.json
"""

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path


def _summarize(values: list[float]) -> dict:
    """Return median, q1, q3 for a list of values."""
    s = sorted(values)
    n = len(s)
    return {
        "median": statistics.median(s),
        "q1": s[n // 4],
        "q3": s[(3 * n) // 4],
        "n_reps": n,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    all_entries = [json.loads(f.read_text()) for f in args.inputs]
    valid = [r for r in all_entries if "error" not in r]
    if len(valid) < len(all_entries):
        n_skipped = len(all_entries) - len(valid)
        print(f"WARNING: skipped {n_skipped} entry/entries with build errors")

    # Group by (n_samples, threads)
    groups = defaultdict(list)
    for r in valid:
        groups[(r["n_samples"], r["threads"])].append(r)

    results = []
    for (n_samples, threads), entries in sorted(groups.items()):
        wall_s_vals = [e["wall_s"] for e in entries]
        rss_vals = [e["peak_rss_mb"] for e in entries]
        db_vals = [e["db_size_mb"] for e in entries]

        row = {
            "n_samples": n_samples,
            "threads": threads,
            "raw_vcf_size_mb": entries[0].get("raw_vcf_size_mb", 0),
        }

        if len(entries) == 1:
            # Single rep: keep flat format for backward compat
            row["wall_s"] = entries[0]["wall_s"]
            row["peak_rss_mb"] = entries[0]["peak_rss_mb"]
            row["db_size_mb"] = entries[0]["db_size_mb"]
        else:
            ws = _summarize(wall_s_vals)
            row["wall_s"] = ws["median"]
            row["wall_s_q1"] = ws["q1"]
            row["wall_s_q3"] = ws["q3"]
            row["wall_s_n_reps"] = ws["n_reps"]
            row["peak_rss_mb"] = statistics.median(rss_vals)
            row["db_size_mb"] = statistics.median(db_vals)

        results.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2))
    print(f"build_perf.json written to {args.output} ({len(results)} entries)")


if __name__ == "__main__":
    main()
