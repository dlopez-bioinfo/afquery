"""Shared Snakemake rule: download and split 1000 Genomes Phase 3 chr22."""

import os
from pathlib import Path

# Create log directory before any rule runs
onstart:
    os.makedirs("logs", exist_ok=True)

# ONEKG_DIR and ONEKG_MANIFEST are imported from the root Snakefile via
# shared.config (already set in the root namespace when this file is included).


rule download_1kg:
    """Download 1KG chr22, split per-sample VCFs and write manifest.tsv.

    Uses file-existence guards inside the shell script so the rule is
    idempotent: re-running after a partial failure resumes from where it
    stopped without re-downloading what already exists.
    """
    output:
        manifest=str(ONEKG_MANIFEST),
        done=touch(str(ONEKG_DIR / ".download_done")),
    log:
        "logs/download_1kg.log",
    threads: 32
    resources:
        mem_mb=8_000,
        runtime=180,  # minutes
        slurm_extra="--nodes=1",
    shell:
        """
        mkdir -p $(dirname {log}) && module load BCFtools/1.18-GCC-12.3.0 parallel 2>/dev/null || true
        export ONEKG_DIR="{ONEKG_DIR}"
        export THREADS={threads}
        bash {workflow.basedir}/shared/scripts/download_1kg.sh \
            > {log} 2>&1
        """
