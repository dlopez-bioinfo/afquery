# Capture Kit Mixing Benchmark

Quantify the impact of capture kit mixing on allele frequency-based variant interpretation.

## Motivation

When datasets combine multiple capture kits (e.g., SureSelect v5, v6, v7), per-kit coverage differences cause different samples to have different detection sensitivities. This benchmark measures the resulting classification discordance for ACMG criteria (BA1, BS1, PM2).

## Overview

1. **Sample generation** — Mix 1000 Genomes samples using three Agilent SureSelect kits in three scenarios
2. **Database builds** — One AFQuery database per scenario
3. **ACMG classification** — Test ACMG criteria thresholds across kits
4. **Metrics** — Coverage overlap, discordance rates, impact on interpretation

## Configuration

Edit `config.py`:

```python
SCENARIOS = {
    "balanced": {...},   # Even distribution (334/333/333)
    "skewed":   {...},   # Realistic (600/300/100)
    "extreme":  {...},   # Worst case (800/150/50)
}

ACMG_THRESHOLDS = {
    "cardiomyopathy": {...},
    "metabolic": {...},
}
```

Uses 1000 Genomes Phase 3 data (chr22) — shared with performance benchmarks.

## Quick Start

```bash
# 1. Download BED files (manual—see 01_prepare_samples.py for instructions)

# 2. Prepare masked VCFs per kit
python 01_prepare_samples.py

# 3. Build AFQuery databases
python 02_build_databases.py

# 4. Compute coverage metrics
python 03_compute_metrics.py

# 5. Classify variants by ACMG criteria
python 04_classify_acmg.py

# 6. Generate plots
python 05_plot_figures.py
```

## Data Files

BED files are downloaded manually from Agilent SureDesign (chr22 pre-filtered versions):

- `SureSelect_v5.bed`
- `SureSelect_v6.bed`
- `SureSelect_v7.bed`

Place in `{DATA_DIR}/capture_kit/beds/` before running.

## Coverage Metrics (chr22)

From Agilent SureDesign:

| Metric | Value |
|--------|-------|
| 3-kit intersection | 863,486 bp |
| Union (any kit) | 1,506,718 bp |
| Discordant (1–2 kits) | 643,232 bp (42.7%) |

## Output

- `results/*.csv` — classification discordance, metrics
- `figures/*.pdf` — publication-quality plots

## References

- Whiffin N et al. (2017). Cardiovascular disease risk gene burden in the founding population of the Amish. *Human Molecular Genetics* 26(17):3423–3428. DOI: 10.1093/hmg/ddx265
- Kalia SS et al. (2017). Recommendations for reporting of secondary findings in clinical exome and genome sequencing, 2016 update. *Genetics in Medicine* 19(2):249–255. DOI: 10.1038/gim.2016.190
- Agilent SureDesign: https://earray.chem.agilent.com/suredesign/
