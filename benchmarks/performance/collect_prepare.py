#!/usr/bin/env python3
"""Aggregate per-DB build JSONs into data_inventory.json.

Usage:
    python collect_prepare.py \
        --onekg-jsons results/raw/prepare/1kg_500.json ... \
        --synth-jsons results/raw/prepare/synth_1000.json ... \
        --output results/data_inventory.json
"""

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--onekg-jsons", nargs="+", type=Path, required=True)
    parser.add_argument("--synth-jsons", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    inventory = {"1kg": {}, "synth": {}}

    for f in args.onekg_jsons:
        data = json.loads(f.read_text())
        inventory["1kg"][str(data["n_samples"])] = data["db_path"]

    for f in args.synth_jsons:
        data = json.loads(f.read_text())
        inventory["synth"][str(data["n_samples"])] = data["db_path"]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(inventory, indent=2))
    print(f"data_inventory.json written to {args.output}")


if __name__ == "__main__":
    main()
