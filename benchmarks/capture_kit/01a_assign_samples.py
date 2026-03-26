#!/usr/bin/env python3
"""Subsample 1KG panel, assign kit technologies, write assignments.json.

Outputs:
  - manifests/assignments.json  — selected samples + scenario→sample→tech mapping
  - per_sample/ symlinks        — symlinks to per-sample VCFs (for rules that need them)

No subprocess calls: this script is purely Python.
Per-sample VCF splitting and BED masking are handled as separate Snakemake rules.
"""

import json
import logging
import random
import sys
from pathlib import Path

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from config import (
    MANIFEST_DIR,
    N_SAMPLES,
    ONEKG_PANEL,
    ONEKG_VCF_DIR,
    PER_SAMPLE_DIR,
    SCENARIOS,
    SEED,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def read_panel(panel_path: Path) -> list[tuple[str, str]]:
    """Read 1KG panel file, return [(sample_name, sex), ...]."""
    samples = []
    for line in panel_path.read_text().strip().split("\n"):
        if line.startswith("sample"):
            continue
        parts = line.split("\t")
        sample, sex = parts[0], parts[3]
        samples.append((sample, sex))
    return samples


def subsample(all_samples, n, seed):
    """Randomly subsample n samples."""
    rng = random.Random(seed)
    if n >= len(all_samples):
        return list(all_samples)
    return rng.sample(all_samples, n)


def assign_technologies(samples, scenario, seed):
    """Assign kit versions to samples by shuffling then slicing.

    Returns: {sample_name: tech_name} dict.
    """
    rng = random.Random(seed)
    shuffled = list(samples)
    rng.shuffle(shuffled)

    assignments = {}
    offset = 0
    for tech_name, count in scenario.items():
        for sample, _sex in shuffled[offset : offset + count]:
            assignments[sample] = tech_name
        offset += count

    return assignments


def main():
    ensure_dirs()

    # 1. Read panel, subsample
    all_samples = read_panel(ONEKG_PANEL)
    selected = subsample(all_samples, N_SAMPLES, SEED)
    logger.info("Selected %d samples (seed=%d)", len(selected), SEED)

    # 2. Keep only samples that have a VCF (already split by a prior run or rule)
    available = [
        (s, sex) for s, sex in selected
        if (ONEKG_VCF_DIR / f"{s}.vcf.gz").exists()
    ]
    if len(available) < len(selected):
        logger.warning(
            "Only %d/%d samples have VCFs; using available subset",
            len(available), len(selected),
        )
    selected = available if available else selected

    # 3. Symlink per-sample VCFs into our working directory (idempotent)
    for sample, _ in selected:
        src = ONEKG_VCF_DIR / f"{sample}.vcf.gz"
        dst = PER_SAMPLE_DIR / f"{sample}.vcf.gz"
        if src.exists() and not dst.exists():
            dst.symlink_to(src)

    # 4. Build scenario assignments: {scenario: {sample: tech}}
    scenarios = {}
    for scenario_name, scenario_dist in SCENARIOS.items():
        scenarios[scenario_name] = assign_technologies(selected, scenario_dist, SEED)

    # 5. Write assignments.json
    assignments = {
        "samples": [[s, sex] for s, sex in selected],
        "scenarios": scenarios,
    }
    out_path = MANIFEST_DIR / "assignments.json"
    out_path.write_text(json.dumps(assignments, indent=2))
    logger.info("Wrote assignments: %s (%d samples, %d scenarios)",
                out_path, len(selected), len(scenarios))


if __name__ == "__main__":
    main()
