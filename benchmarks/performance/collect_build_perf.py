#!/usr/bin/env python3
"""Aggregate per-run build JSONs into build_perf.json.

Usage:
    python collect_build_perf.py \
        --inputs results/raw/build/build_1000_1t.json ... \
        --output results/build_perf.json
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
    results.sort(key=lambda r: (r["n_samples"], r["threads"]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2))
    print(f"build_perf.json written to {args.output} ({len(results)} entries)")


if __name__ == "__main__":
    main()
