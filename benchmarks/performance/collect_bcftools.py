#!/usr/bin/env python3
"""Aggregate per-subset bcftools JSONs into bcftools_comparison.json.

Note: concordance.json is produced directly by the concordance rule and
does not go through this script.

Usage:
    python collect_bcftools.py \
        --inputs results/raw/bcftools/bcftools_500.json ... \
        --output results/bcftools_comparison.json
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
    print(f"bcftools_comparison.json written to {args.output} ({len(results)} entries)")


if __name__ == "__main__":
    main()
