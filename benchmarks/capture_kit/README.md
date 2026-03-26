# Capture Kit Mixing Benchmark

Quantify the impact of capture kit mixing on allele frequency-based variant interpretation.

## Motivation

When datasets combine multiple capture kits (e.g., SureSelect v5, v6, v7), per-kit coverage differences cause different samples to have different detection sensitivities. This benchmark measures the resulting classification discordance for ACMG criteria (BA1, BS1, PM2).

## Overview

1. **Sample generation** — Mix 1000 Genomes samples using three Agilent SureSelect kits in three scenarios
2. **Database builds** — One AFQuery database per scenario (+ WGS ground truth)
3. **Metrics** — NADI (Normalized AN Deviation Index) and AF error vs. WGS truth
4. **ACMG classification** — Misclassification rates per scenario and disease model
5. **Figures** — High-resolution plots

## Running

All steps are orchestrated by Snakemake from the `benchmarks/` root.

```bash
# Run the full capture kit benchmark (HPC)
snakemake --profile benchmarks/profiles/slurm capture_kit_all

# Run locally
snakemake --cores all capture_kit_all

# Dry run
snakemake --profile benchmarks/profiles/slurm --dry-run capture_kit_all
```

The scripts can also be run directly when working from this directory
(1KG data must already exist — run `snakemake download_1kg` first):

```bash
python 01_prepare_samples.py
python 02_build_databases.py
python 03_compute_metrics.py
python 04_classify_acmg.py
python 05_plot_figures.py
```

## BED Files (manual prerequisite)

BED files must be downloaded manually from Agilent SureDesign before running
the benchmark. Two sets are needed, placed in `{DATA_DIR}/capture_kit/beds/`:

- `masking/SureSelect_v5.bed` — bare chromosome names (e.g. `22`), for `bedtools intersect`
- `masking/SureSelect_v6.bed`
- `masking/SureSelect_v7.bed`
- `afquery/SureSelect_v5.bed` — `chr`-prefixed names (e.g. `chr22`), for AFQuery's CaptureIndex
- `afquery/SureSelect_v6.bed`
- `afquery/SureSelect_v7.bed`

Pre-filter BED files to chr22 before placing them.

## Coverage Metrics (chr22)

| Metric | Value |
|--------|-------|
| 3-kit intersection | 863,486 bp |
| Union (any kit) | 1,506,718 bp |
| Discordant (1–2 kits) | 643,232 bp (42.7%) |

## Configuration

Scenarios and ACMG thresholds are in `config.py`:

```python
SCENARIOS = {
    "balanced": {"SureSelect_v5": 334, "SureSelect_v6": 333, "SureSelect_v7": 333},
    "skewed":   {"SureSelect_v5": 600, "SureSelect_v6": 300, "SureSelect_v7": 100},
    "extreme":  {"SureSelect_v5": 800, "SureSelect_v6": 150, "SureSelect_v7":  50},
}

ACMG_THRESHOLDS = {
    "cardiomyopathy": {"BA1": 0.05, "BS1": 0.001,  "PM2": 0.0001},
    "metabolic":      {"BA1": 0.05, "BS1": 0.0001, "PM2": 0.0},
}
```

Shared constants (DATA_DIR, 1KG paths, SEED) come from `../shared/config.py`.

## Output

- `results/merged.parquet` — merged AF comparison table (all scenarios)
- `results/nadi_summary.json` — summary statistics per scenario
- `results/acmg_results.json` — misclassification counts per disease/scenario
- `figures/fig_capkit_*.{pdf,png}` — high-resolution plots

## References

- Whiffin N et al. (2017). Cardiovascular disease risk gene burden in the founding population of the Amish. *Human Molecular Genetics* 26(17):3423–3428. DOI: 10.1093/hmg/ddx265
- Kalia SS et al. (2017). Recommendations for reporting of secondary findings in clinical exome and genome sequencing, 2016 update. *Genetics in Medicine* 19(2):249–255. DOI: 10.1038/gim.2016.190
- Agilent SureDesign: https://earray.chem.agilent.com/suredesign/
