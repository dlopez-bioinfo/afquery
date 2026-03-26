"""Shared Snakemake rules: download and split 1000 Genomes Phase 3 chr22.

Three rules handle the full pipeline:

  checkpoint download_1kg_data  — wget chr22 VCF + panel (1 job)
  rule      split_1kg_sample    — bcftools view for one sample (2504 independent jobs)
  rule      download_1kg        — write manifest.tsv once all splits are done (1 job)

ONEKG_DIR and ONEKG_MANIFEST are imported from the root Snakefile namespace.
"""

import os
from pathlib import Path

# Create log directory before any rule runs
onstart:
    os.makedirs("logs/download_1kg", exist_ok=True)

# Derive paths from ONEKG_DIR (set by root Snakefile, updated by _set_data_dir()).
# Using Path(str(...)) avoids stale references if ONEKG_DIR is reassigned.
_ONEKG_VCF_DIR    = Path(str(ONEKG_DIR)) / "vcfs"
_ONEKG_MERGED_VCF = (
    Path(str(ONEKG_DIR))
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
_ONEKG_PANEL = Path(str(ONEKG_DIR)) / "integrated_call_samples_v3.20130502.ALL.panel"

_FTP_BASE   = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
_CHR22_VCF  = "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
_PANEL_FILE = "integrated_call_samples_v3.20130502.ALL.panel"


# ---------------------------------------------------------------------------
# Step 1: download merged VCF + panel
# ---------------------------------------------------------------------------

checkpoint download_1kg_data:
    """Download chr22 multi-sample VCF, its index, and the sample panel file."""
    output:
        merged_vcf=str(_ONEKG_MERGED_VCF),
        tbi=str(_ONEKG_MERGED_VCF) + ".tbi",
        panel=str(_ONEKG_PANEL),
    log:
        "logs/download_1kg/download_data.log",
    threads: 1
    resources:
        mem_mb=4_000,
        runtime=120,
    conda: config["conda_env_file"]
    shell:
        """
        set -euo pipefail
        LOG=$(realpath {log}) && mkdir -p "$(dirname "$LOG")" "$(dirname "{output.merged_vcf}")"
        wget -c -q --show-progress \
            -O "{output.merged_vcf}" "{_FTP_BASE}/{_CHR22_VCF}" >> "$LOG" 2>&1
        wget -c -q --show-progress \
            -O "{output.tbi}" "{_FTP_BASE}/{_CHR22_VCF}.tbi" >> "$LOG" 2>&1
        wget -c -q --show-progress \
            -O "{output.panel}" "{_FTP_BASE}/{_PANEL_FILE}" >> "$LOG" 2>&1
        """


# ---------------------------------------------------------------------------
# Step 2: split per-sample VCFs — one independent SLURM job per sample
# ---------------------------------------------------------------------------

rule split_1kg_sample:
    """Extract one per-sample VCF from the merged 1KG chr22 VCF."""
    input:
        merged=str(_ONEKG_MERGED_VCF),
        tbi=str(_ONEKG_MERGED_VCF) + ".tbi",
    output:
        vcf=str(_ONEKG_VCF_DIR / "{sample}.vcf.gz"),
    log:
        "logs/download_1kg/split/{sample}.log",
    threads: 1
    resources:
        mem_mb=4_000,
        runtime=10,
    conda: config["conda_env_file"]
    shell:
        """
        set -euo pipefail
        LOG=$(realpath {log}) && mkdir -p "$(dirname "$LOG")" "$(dirname "{output.vcf}")"
        bcftools view -s {wildcards.sample} --min-ac 1:nref {input.merged} \
            -Oz -o {output.vcf} 2>>"$LOG"
        """


# ---------------------------------------------------------------------------
# Step 3: write manifest.tsv once all splits are ready
# ---------------------------------------------------------------------------

def _all_1kg_vcfs(wildcards):
    """Return paths to all per-sample VCFs; re-evaluated after download_1kg_data."""
    checkpoints.download_1kg_data.get(**wildcards)
    samples = []
    with open(_ONEKG_PANEL) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if parts[0] == "sample":
                continue
            samples.append(parts[0])
    return expand(str(_ONEKG_VCF_DIR / "{sample}.vcf.gz"), sample=samples)


rule download_1kg:
    """Write manifest.tsv after all per-sample VCFs have been split.

    This rule keeps the same name and outputs as before so that downstream
    rules (performance, capture_kit) that depend on manifest.tsv or
    .download_done are unaffected.
    """
    input:
        vcfs=_all_1kg_vcfs,
        panel=str(_ONEKG_PANEL),
    output:
        manifest=str(ONEKG_MANIFEST),
        done=touch(str(ONEKG_DIR / ".download_done")),
    params:
        vcf_dir=str(_ONEKG_VCF_DIR),
    log:
        "logs/download_1kg/manifest.log",
    threads: 1
    resources:
        mem_mb=2_000,
        runtime=10,
    conda: config["conda_env_file"]
    shell:
        """
        set -euo pipefail
        LOG=$(realpath {log}) && mkdir -p "$(dirname "$LOG")"
        printf 'sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n' \
            > {output.manifest}
        tail -n +2 {input.panel} \
        | while IFS=$'\t' read -r sample pop super_pop gender; do
            vcf="{params.vcf_dir}/${{sample}}.vcf.gz"
            [ -f "$vcf" ] || continue
            printf '%s\t%s\t%s\t%s\t%s\n' \
                "$sample" "$gender" "WGS" "$vcf" "$super_pop"
        done >> {output.manifest}
        n=$(tail -n +2 {output.manifest} | wc -l)
        echo "Manifest written: $n samples → {output.manifest}" >> "$LOG"
        """
