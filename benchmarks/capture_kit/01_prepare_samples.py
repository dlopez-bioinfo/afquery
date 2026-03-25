#!/usr/bin/env python3
"""Subsample 1KG, split into per-sample VCFs, mask with capture BEDs, write manifests.

Inputs:
  - 1KG Phase 3 chr22 merged VCF (from 00_download_1kg.sh)
  - 1KG panel file (sample metadata)
  - BED files in beds/masking/ (already prepared)

Outputs:
  - per_sample/   symlinks to per-sample VCFs
  - masked/{tech}/ masked per-sample VCFs (variants outside BED removed)
  - manifests/    4 TSV manifests (1 WGS + 3 WES scenarios)
"""

import logging
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from capture_kit_bench.config import (
    CAPTURE_DIR,
    MASKING_BED_DIR,
    MASKED_DIR,
    MANIFEST_DIR,
    N_SAMPLES,
    ONEKG_DIR,
    ONEKG_MERGED_VCF,
    ONEKG_PANEL,
    ONEKG_VCF_DIR,
    PER_SAMPLE_DIR,
    SCENARIOS,
    SEED,
    THREADS,
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


def ensure_per_sample_vcfs(samples, merged_vcf, vcf_dir, threads):
    """Split per-sample VCFs from merged VCF using parallel bcftools view.

    Reuses existing VCFs — only processes samples that don't have a file yet.
    """
    missing = [s for s, _ in samples if not (vcf_dir / f"{s}.vcf.gz").exists()]
    if not missing:
        logger.info("All %d per-sample VCFs already exist in %s", len(samples), vcf_dir)
        return

    logger.info(
        "Splitting %d missing per-sample VCFs (of %d total)...",
        len(missing), len(samples),
    )

    sample_list_file = vcf_dir / "_samples_to_split.txt"
    sample_list_file.write_text("\n".join(missing) + "\n")

    cmd = (
        f"cat {sample_list_file} | parallel -j {threads} "
        f"'bcftools view -s {{}} --min-ac 1:nref {merged_vcf} "
        f"-Oz -o {vcf_dir}/{{}}.vcf.gz 2>/dev/null'"
    )
    subprocess.run(cmd, shell=True, check=True)
    sample_list_file.unlink()
    logger.info("VCF splitting complete.")


def assign_technologies(samples, scenario, seed):
    """Assign kit versions to samples by shuffling then slicing.

    Returns: list of (sample_name, sex, tech_name) tuples.
    """
    rng = random.Random(seed)
    shuffled = list(samples)
    rng.shuffle(shuffled)

    assignments = []
    offset = 0
    for tech_name, count in scenario.items():
        for sample, sex in shuffled[offset : offset + count]:
            assignments.append((sample, sex, tech_name))
        offset += count

    return assignments


def mask_all_vcfs(assignments, per_sample_dir, masking_bed_dir, masked_dir, threads):
    """Mask VCFs using bedtools intersect via GNU parallel.

    For each sample, removes variants outside its assigned kit's BED regions.
    """
    tasks = []
    for sample, _, tech in assignments:
        output_vcf = masked_dir / tech / f"{sample}.vcf.gz"
        if not output_vcf.exists():
            tasks.append(f"{sample}\t{tech}")

    if not tasks:
        logger.info("All masked VCFs already exist")
        return

    task_file = masked_dir / "_mask_tasks.tsv"
    task_file.write_text("\n".join(tasks) + "\n")
    logger.info("Masking %d VCFs with %d threads...", len(tasks), threads)

    cmd = (
        f"cat {task_file} | parallel -j {threads} --colsep '\\t' "
        f"'bedtools intersect "
        f"-a {per_sample_dir}/{{1}}.vcf.gz "
        f"-b {masking_bed_dir}/{{2}}.bed "
        f"-header -wa | bgzip > {masked_dir}/{{2}}/{{1}}.vcf.gz'"
    )
    subprocess.run(cmd, shell=True, check=True)
    task_file.unlink()
    logger.info("Masking complete.")


def write_manifest(assignments, vcf_base_dir, output_path, is_wgs=False):
    """Write AFQuery manifest TSV.

    WGS: tech_name="WGS", vcf from per_sample/ (no masking).
    WES: tech_name=assigned kit, vcf from masked/{tech}/.
    """
    lines = ["sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes"]
    for sample, sex, tech in assignments:
        if is_wgs:
            vcf_path = str(vcf_base_dir / f"{sample}.vcf.gz")
            tech_name = "WGS"
        else:
            vcf_path = str(vcf_base_dir / tech / f"{sample}.vcf.gz")
            tech_name = tech
        lines.append(f"{sample}\t{sex}\t{tech_name}\t{vcf_path}\tBENCH")

    output_path.write_text("\n".join(lines) + "\n")
    logger.info("Wrote manifest: %s (%d samples)", output_path, len(assignments))


def main():
    ensure_dirs()

    # 1. Read panel, subsample
    all_samples = read_panel(ONEKG_PANEL)
    selected = subsample(all_samples, N_SAMPLES, SEED)
    logger.info("Selected %d samples (seed=%d)", len(selected), SEED)

    # 2. Ensure per-sample VCFs exist (split from merged)
    ensure_per_sample_vcfs(selected, ONEKG_MERGED_VCF, ONEKG_VCF_DIR, THREADS)

    # Keep only samples that actually have a VCF
    available = [
        (s, sex) for s, sex in selected
        if (ONEKG_VCF_DIR / f"{s}.vcf.gz").exists()
    ]
    if len(available) < len(selected):
        logger.warning(
            "Only %d/%d samples have VCFs; using available subset",
            len(available), len(selected),
        )
    selected = available[:N_SAMPLES]

    # 3. Symlink per-sample VCFs into our working directory
    for sample, _ in selected:
        src = ONEKG_VCF_DIR / f"{sample}.vcf.gz"
        dst = PER_SAMPLE_DIR / f"{sample}.vcf.gz"
        if not dst.exists():
            dst.symlink_to(src)

    # 4. For each scenario: assign kits, mask VCFs, write WES manifest
    for scenario_name, scenario_dist in SCENARIOS.items():
        logger.info("=== Scenario: %s ===", scenario_name)
        assignments = assign_technologies(selected, scenario_dist, SEED)
        mask_all_vcfs(
            assignments, PER_SAMPLE_DIR, MASKING_BED_DIR, MASKED_DIR, THREADS,
        )
        write_manifest(
            assignments,
            MASKED_DIR,
            MANIFEST_DIR / f"manifest_{scenario_name}.tsv",
            is_wgs=False,
        )

    # 5. Write WGS manifest (all samples as WGS, unmasked VCFs)
    wgs_assignments = [(s, sex, "WGS") for s, sex in selected]
    write_manifest(
        wgs_assignments,
        PER_SAMPLE_DIR,
        MANIFEST_DIR / "manifest_wgs.tsv",
        is_wgs=True,
    )

    logger.info("Done. %d samples prepared.", len(selected))


if __name__ == "__main__":
    main()
