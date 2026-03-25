"""Centralized configuration for all benchmark scripts.

Edit DATA_DIR to point to a directory with enough space for VCFs and databases.
All other constants control experiment parameters and should be kept fixed for
reproducibility.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths — adjust DATA_DIR to your environment
# ---------------------------------------------------------------------------
BENCHMARKS_DIR = Path(__file__).resolve().parent
DATA_DIR = Path("/mnt/lustre/home/dlopez/projects/afquery_bench_data")

# 1000 Genomes data
ONEKG_DIR = DATA_DIR / "1kg"
ONEKG_VCF_DIR = ONEKG_DIR / "vcfs"
ONEKG_MERGED_VCF = (
    ONEKG_DIR
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
ONEKG_PANEL = ONEKG_DIR / "integrated_call_samples_v3.20130502.ALL.panel"

# Databases built from 1KG data
ONEKG_DB_DIR = DATA_DIR / "1kg_dbs"

# Synthetic data
SYNTH_DIR = DATA_DIR / "synth"
SYNTH_DB_DIR = DATA_DIR / "synth_dbs"

# Results and figures
RESULTS_DIR = BENCHMARKS_DIR / "results"
FIGURES_DIR = BENCHMARKS_DIR / "figures"

# ---------------------------------------------------------------------------
# Experiment parameters
# ---------------------------------------------------------------------------
SEED = 42
GENOME_BUILD = "GRCh37"
ONEKG_CHROM = "22"  # 1KG Phase 3 uses "22", not "chr22"

# Sample subsets for 1KG experiments
ONEKG_SUBSETS = [500, 1_000, 2_504]

# Synthetic scaling experiment
SYNTH_SCALES = [1_000, 5_000, 10_000, 25_000, 50_000]
SYNTH_VARIANTS_PER_CHROM = 10_000
SYNTH_CHROMS = ("chr22",)

# Query scaling (Experiment 1)
QUERY_COLD_REPS = 5
QUERY_WARM_REPS = 50
QUERY_REGION_REPS = 20
QUERY_BATCH_REPS = 20
QUERY_WARMUP = 3  # discarded warmup iterations

# Build scaling (Experiment 2)
BUILD_THREAD_COUNTS = [1, 4, 8, 16, 32]
BUILD_SCALES = [1_000, 5_000, 10_000]

# Annotation throughput (Experiment 3)
ANNOTATE_THREAD_COUNTS = [1, 4, 8, 16, 32]
ANNOTATE_VARIANT_COUNTS = [10_000, 50_000, 100_000]

# bcftools comparison (Experiment 4)
BCFTOOLS_REPS = 10

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def ensure_dirs() -> None:
    """Create all output directories if they don't exist."""
    for d in [
        DATA_DIR,
        ONEKG_DIR,
        ONEKG_VCF_DIR,
        ONEKG_DB_DIR,
        SYNTH_DIR,
        SYNTH_DB_DIR,
        RESULTS_DIR,
        FIGURES_DIR,
    ]:
        d.mkdir(parents=True, exist_ok=True)
