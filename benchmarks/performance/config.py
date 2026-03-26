"""Configuration for the performance benchmark.

Shared constants (DATA_DIR, 1KG paths, SEED, GENOME_BUILD) come from
benchmarks/shared/config.py and are re-exported here for convenience.
"""

import sys
from pathlib import Path

_BENCH_DIR = Path(__file__).resolve().parent.parent  # benchmarks/
sys.path.insert(0, str(_BENCH_DIR))

from shared.config import (  # noqa: E402
    DATA_DIR,
    GENOME_BUILD,
    ONEKG_DIR,
    ONEKG_MANIFEST,
    ONEKG_MERGED_VCF,
    ONEKG_PANEL,
    ONEKG_VCF_DIR,
    ONEKG_CHROM,
    SEED,
)

# ---------------------------------------------------------------------------
# Performance-specific paths
# ---------------------------------------------------------------------------
BENCHMARKS_DIR = Path(__file__).resolve().parent

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
    from shared.utils import ensure_dirs as _ensure
    _ensure(
        DATA_DIR,
        ONEKG_DIR,
        ONEKG_VCF_DIR,
        ONEKG_DB_DIR,
        SYNTH_DIR,
        SYNTH_DB_DIR,
        RESULTS_DIR,
        FIGURES_DIR,
    )
