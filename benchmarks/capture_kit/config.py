"""Configuration for capture kit mixing benchmark.

BED files downloaded from Agilent SureDesign
(https://earray.chem.agilent.com/suredesign/) and pre-filtered to chr22.

Coverage overlap on chr22:
  All 3 kits intersection: 863,486 bases
  Union (any kit):       1,506,718 bases
  Discordant (1-2 kits):   643,232 bases (42.7%)
"""

import importlib.util
import os
from pathlib import Path

# Import shared 1KG config from performance benchmark
perf_config_path = Path(__file__).resolve().parent.parent / "performance" / "config.py"
spec = importlib.util.spec_from_file_location("perf_config", perf_config_path)
perf_config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(perf_config)

DATA_DIR = perf_config.DATA_DIR
GENOME_BUILD = perf_config.GENOME_BUILD
ONEKG_DIR = perf_config.ONEKG_DIR
ONEKG_MERGED_VCF = perf_config.ONEKG_MERGED_VCF
ONEKG_PANEL = perf_config.ONEKG_PANEL
ONEKG_VCF_DIR = perf_config.ONEKG_VCF_DIR
SEED = perf_config.SEED

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CAPTURE_DIR = DATA_DIR / "capture_kit"
BED_DIR = CAPTURE_DIR / "beds"
MASKING_BED_DIR = BED_DIR / "masking"
AFQUERY_BED_DIR = BED_DIR / "afquery"
PER_SAMPLE_DIR = CAPTURE_DIR / "per_sample"
MASKED_DIR = CAPTURE_DIR / "masked"
MANIFEST_DIR = CAPTURE_DIR / "manifests"
DB_DIR = CAPTURE_DIR / "dbs"
RESULTS_DIR = CAPTURE_DIR / "results"
FIGURES_DIR = Path(__file__).resolve().parent / "figures"

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
N_SAMPLES = 1000
THREADS = os.cpu_count() or 1

SCENARIOS = {
    "balanced": {"SureSelect_v5": 334, "SureSelect_v6": 333, "SureSelect_v7": 333},
    "skewed":   {"SureSelect_v5": 600, "SureSelect_v6": 300, "SureSelect_v7": 100},
    "extreme":  {"SureSelect_v5": 800, "SureSelect_v6": 150, "SureSelect_v7":  50},
}

# ACMG thresholds (ClinGen-inspired, two disease models)
# cardiomyopathy: Whiffin et al. 2017 max credible AF ~ 0.04% -> BS1 > 0.1%
# metabolic: rare disease -> BS1 > 0.01%, PM2 = absent (AC=0)
ACMG_THRESHOLDS = {
    "cardiomyopathy": {"BA1": 0.05, "BS1": 0.001, "PM2": 0.0001},
    "metabolic":      {"BA1": 0.05, "BS1": 0.0001, "PM2": 0.0},
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def ensure_dirs():
    """Create all output directories if they don't exist."""
    for d in [CAPTURE_DIR, BED_DIR, MASKING_BED_DIR, AFQUERY_BED_DIR,
              PER_SAMPLE_DIR, MASKED_DIR, MANIFEST_DIR, DB_DIR,
              RESULTS_DIR, FIGURES_DIR]:
        d.mkdir(parents=True, exist_ok=True)
    for tech in SCENARIOS["balanced"]:
        (MASKED_DIR / tech).mkdir(parents=True, exist_ok=True)
