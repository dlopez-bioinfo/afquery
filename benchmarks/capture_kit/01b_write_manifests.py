#!/usr/bin/env python3
"""Write AFQuery manifest TSVs from assignments.json.

Reads manifests/assignments.json (written by 01a_assign_samples.py after all
per-sample VCFs are split and masked) and writes 4 manifest TSV files:
  - manifest_wgs.tsv        — all samples as WGS (unmasked VCFs)
  - manifest_balanced.tsv   — balanced kit scenario
  - manifest_skewed.tsv     — skewed kit scenario
  - manifest_extreme.tsv    — extreme kit scenario

No subprocess calls: pure Python.
"""

import json
import logging
import sys
from pathlib import Path

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from config import (
    MANIFEST_DIR,
    MASKED_DIR,
    PER_SAMPLE_DIR,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def write_manifest(assignments, output_path, is_wgs=False):
    """Write AFQuery manifest TSV.

    assignments: list of (sample_name, sex, tech_name) tuples.
    WGS: vcf from per_sample/ (unmasked). WES: vcf from masked/{tech}/.
    """
    lines = ["sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes"]
    for sample, sex, tech in assignments:
        if is_wgs:
            vcf_path = str(PER_SAMPLE_DIR / f"{sample}.vcf.gz")
            tech_name = "WGS"
        else:
            vcf_path = str(MASKED_DIR / tech / f"{sample}.vcf.gz")
            tech_name = tech
        lines.append(f"{sample}\t{sex}\t{tech_name}\t{vcf_path}\tBENCH")

    output_path.write_text("\n".join(lines) + "\n")
    logger.info("Wrote manifest: %s (%d samples)", output_path, len(assignments))


def main():
    ensure_dirs()

    assignments_path = MANIFEST_DIR / "assignments.json"
    data = json.loads(assignments_path.read_text())

    samples = data["samples"]          # [[name, sex], ...]
    sex_map = {s: sex for s, sex in samples}
    scenarios = data["scenarios"]      # {scenario: {sample: tech}}

    # WGS manifest — all samples, unmasked
    wgs_assignments = [(s, sex, "WGS") for s, sex in samples]
    write_manifest(wgs_assignments, MANIFEST_DIR / "manifest_wgs.tsv", is_wgs=True)

    # WES manifests — one per scenario
    for scenario_name, tech_map in scenarios.items():
        wes_assignments = [
            (sample, sex_map[sample], tech)
            for sample, tech in tech_map.items()
        ]
        write_manifest(
            wes_assignments,
            MANIFEST_DIR / f"manifest_{scenario_name}.tsv",
            is_wgs=False,
        )

    logger.info("Done. %d samples, %d WES scenarios.", len(samples), len(scenarios))


if __name__ == "__main__":
    main()
