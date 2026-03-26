#!/usr/bin/env python3
"""Aggregate per-run annotation JSONs into annotate_throughput.json.

Computes speedup ratios: for each n_variants group, speedup = wall_s(1t) / wall_s(Nt).

Usage:
    python collect_annotate.py \
        --inputs results/raw/annotate/annotate_10000_1t.json ... \
        --output results/annotate_throughput.json
"""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    results = [json.loads(f.read_text()) for f in args.inputs]
    results.sort(key=lambda r: (r["n_variants"], r["threads"]))

    # Compute speedup ratios per n_variants group
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
