"""Shared configuration for all benchmarks.

Edit DATA_DIR (or set the AFQUERY_BENCH_DATA environment variable) to point
to a directory with enough space for VCFs and databases.
"""

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths — adjust via environment variable or edit directly
# ---------------------------------------------------------------------------
DATA_DIR = Path(
    os.environ.get(
        "AFQUERY_BENCH_DATA",
        "/mnt/lustre/home/dlopez/projects/afquery_bench_data",
    )
)

# 1000 Genomes Phase 3 chr22
ONEKG_DIR = DATA_DIR / "1kg"
ONEKG_VCF_DIR = ONEKG_DIR / "vcfs"
ONEKG_MERGED_VCF = (
    ONEKG_DIR
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
ONEKG_PANEL = ONEKG_DIR / "integrated_call_samples_v3.20130502.ALL.panel"
ONEKG_MANIFEST = ONEKG_DIR / "manifest.tsv"

# ---------------------------------------------------------------------------
# Experiment constants
# ---------------------------------------------------------------------------
SEED = 42
GENOME_BUILD = "GRCh37"
ONEKG_CHROM = "22"  # Phase 3 uses bare "22", not "chr22"
