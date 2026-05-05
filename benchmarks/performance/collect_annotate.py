#!/usr/bin/env python3
"""Aggregate per-run annotation JSONs into annotate_throughput.json.

When multiple repetitions exist per (n_variants, threads), computes
median and IQR for wall_s and throughput. Speedup ratios are computed
on median values.

Usage:
    python collect_annotate.py \
        --inputs results/raw/annotate/annotate_10000_1t_r1.json ... \
        --output results/annotate_throughput.json
"""

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    all_entries = [json.loads(f.read_text()) for f in args.inputs]

    # Group by (n_variants, threads)
    groups = defaultdict(list)
    for r in all_entries:
        groups[(r["n_variants"], r["threads"])].append(r)

    results = []
    for (n_variants, threads), entries in sorted(groups.items()):
        wall_vals = [e["wall_s"] for e in entries]
        tp_vals = [e["throughput_vars_per_s"] for e in entries]

        row = {
            "n_variants": n_variants,
            "threads": threads,
        }

        if len(entries) == 1:
            row["wall_s"] = entries[0]["wall_s"]
            row["throughput_vars_per_s"] = entries[0]["throughput_vars_per_s"]
        else:
            med_wall = statistics.median(wall_vals)
            med_tp = statistics.median(tp_vals)
            s = sorted(wall_vals)
            n = len(s)
            row["wall_s"] = round(med_wall, 3)
            row["wall_s_q1"] = round(s[n // 4], 3)
            row["wall_s_q3"] = round(s[(3 * n) // 4], 3)
            row["wall_s_n_reps"] = n
            row["throughput_vars_per_s"] = round(med_tp, 1)

        results.append(row)

    # Compute speedup ratios per n_variants group (on median wall_s)
    n_variants_values = sorted(set(r["n_variants"] for r in results))
    for n_variants in n_variants_values:
        baseline = next(
            (r for r in results if r["n_variants"] == n_variants and r["threads"] == 1),
            None,
        )
        if baseline:
            for r in results:
                if r["n_variants"] == n_variants:
                    r["speedup"] = (
                        round(baseline["wall_s"] / r["wall_s"], 2)
                        if r["wall_s"] > 0
                        else 0
                    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2))
    print(f"annotate_throughput.json written to {args.output} ({len(results)} entries)")


if __name__ == "__main__":
    main()
