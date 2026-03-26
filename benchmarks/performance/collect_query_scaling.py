#!/usr/bin/env python3
"""Aggregate per-scale query JSONs into query_scaling.json.

Usage:
    python collect_query_scaling.py \
        --inputs results/raw/query/query_1000.json ... \
        --output results/query_scaling.json
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
    results.sort(key=lambda r: r["n_samples"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2))
    print(f"query_scaling.json written to {args.output} ({len(results)} entries)")


if __name__ == "__main__":
    main()
