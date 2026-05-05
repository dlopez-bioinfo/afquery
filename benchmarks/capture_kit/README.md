# Capture Kit Mixing Benchmark

Quantify the impact of capture kit mixing on allele frequency-based variant interpretation.

## Motivation

When datasets combine multiple capture kits (e.g., SureSelect v5, v6, v7), per-kit coverage differences cause different samples to have different detection sensitivities. This benchmark measures the resulting classification discordance for ACMG criteria (BA1, BS1, PM2).

## Overview

1. **Sample assignment** (`01a_assign_samples.py`) -- Subsample 1000 Genomes, assign kits per scenario
2. **VCF masking** (Snakemake rule) -- bedtools intersect per (sample, kit) pair
3. **Manifest generation** (`01b_write_manifests.py`) -- Write AFQuery manifests (1 WGS + 3 WES)
4. **Database builds** (`02_build_databases.py`) -- One AFQuery database per scenario (+ WGS ground truth)
5. **Metrics** (`03_compute_metrics.py`) -- NADI (Normalized AN Deviation Index) and AF error vs. WGS truth
6. **ACMG classification** (`04_classify_acmg.py`) -- Misclassification rates per scenario and disease model
7. **Figures** (`05_plot_figures.py`) -- High-resolution plots

## Running

All steps are orchestrated by Snakemake from the `benchmarks/` root.

```bash
# Run the full capture kit benchmark
snakemake --cores 52 capture_kit_all

# Dry run
snakemake --cores 52 --dry-run capture_kit_all

# Smoke test (50 samples, balanced scenario only)
snakemake --cores 52 --config smoke_test=true capture_kit_all
```

## BED Files

BED files are committed to the repository under `benchmarks/capture_kit/beds/`:

```
beds/
├── SureSelect_v5.bed
├── SureSelect_v6.bed
└── SureSelect_v7.bed
```

These are pre-filtered to chr22 with `chr`-prefixed chromosome names. The same BED files are used both for VCF masking (via bedtools, after `awk` adds `chr` prefix to 1KG VCFs) and for AFQuery's CaptureIndex (technology-aware AN computation).

## Coverage Metrics (chr22)

| Metric | Value |
|--------|-------|
| 3-kit intersection | 863,486 bp |
| Union (any kit) | 1,506,718 bp |
| Discordant (1-2 kits) | 643,232 bp (42.7%) |

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

Technology assignments use per-scenario deterministic seeds derived from the global SEED, ensuring independent (non-correlated) kit assignments across scenarios.

## Experiment Design

### The Problem: AN Distortion from Kit Mixing

When a cohort mixes exome capture kits, not all samples cover all positions. Consider a variant at position X covered by SureSelect v5 but not v6. If 600 samples used v5 and 400 used v6, only 600 samples actually have data at position X. The correct AN (allele number) is 1,200 (600 diploid samples x 2), not 2,000 (1000 x 2).

A **naive approach** uses AN = 2N for all positions, inflating AN and deflating AF. AFQuery's **technology-aware approach** knows which kit each sample used, checks the BED file, and excludes uncovered samples from AN.

This benchmark quantifies the magnitude of that distortion and its downstream clinical impact.

### Ground Truth and Scenarios

A WGS (whole-genome sequencing) database serves as ground truth: all 1,000 samples cover every position, so AN is always 2,000 on autosomes. Three WES scenarios test increasing levels of kit imbalance:

| Scenario | v5 | v6 | v7 | Intent |
|----------|----|----|----|----|
| Balanced | 334 | 333 | 333 | Equal distribution; minimal distortion expected |
| Skewed | 600 | 300 | 100 | Realistic imbalance (one dominant kit) |
| Extreme | 800 | 150 | 50 | Pathological case; maximum distortion |

For each scenario, two AF estimates are compared against WGS truth:
- **AF_AFQuery** = AC_wes / AN_wes (technology-aware AN from AFQuery)
- **AF_naive** = AC_wes / AN_wgs (naive AN = 2N, ignoring kit coverage)

### Metrics

**AF error** = |AF_method - AF_wgs| for each variant. Reported as MAE (mean absolute error), median AE, and max error, overall and stratified by number of covering kits (1, 2, or 3).

**NADI** (Normalized AN Deviation Index) = |AN_naive - AN_true| / AN_true. Measures how far the naive AN is from the correct value. NADI_AFQuery is 0 by definition (AFQuery computes AN correctly).

**AN ratio** = AN_naive / AN_AFQuery = AN_wgs / AN_wes. Values > 1 indicate AN inflation. A ratio of 1.5 means the naive approach overcounts alleles by 50%.

### ACMG Classification Analysis

Variants are classified into ACMG evidence categories based on AF thresholds:

| Category | Meaning | Threshold |
|----------|---------|-----------|
| BA1 | Benign (standalone) | AF > 5% |
| BS1 | Benign (strong) | AF > disease-specific threshold |
| PM2 | Pathogenic (moderate) | AF <= threshold or absent |
| Neutral | None of the above | -- |

Two disease models with different stringency are tested:
- **Cardiomyopathy:** BS1 > 0.1%, PM2 <= 0.01% (moderate thresholds)
- **Metabolic:** BS1 > 0.01%, PM2 = absent (AC=0; very strict)

For each (scenario, disease) combination, variants are classified using AF_wgs (truth), AF_AFQuery, and AF_naive. The key metrics are:
- **Discordance rate:** fraction of variants classified differently from truth
- **Toward-pathogenic rate:** fraction of variants where the method predicts a more pathogenic-supporting class than truth (e.g., BS1→neutral, BA1→BS1). These are errors caused by AF deflation — the method underestimates allele frequency.
- **Toward-benign rate:** fraction of variants where the method predicts a more benign-supporting class than truth (e.g., neutral→BS1). These reflect AF overestimation due to sampling variance in low-coverage kit positions.
- **BS1 recall:** fraction of true BS1 variants correctly classified as BS1.
- **False pathogenic (FP):** truth is benign (BA1/BS1), method calls PM2 (clinically dangerous: could trigger unnecessary follow-up)
- **Missed pathogenic (MP):** truth is PM2, method calls benign (clinically dangerous: could miss a pathogenic variant)

---

## Figures

### Figure 1: AF Error Distribution (`fig_capkit_af_error_violin`)

Violin plot (log-scale Y). Three scenario groups along X-axis. Within each group, two violins: AFQuery (blue) and naive (orange). Only non-zero errors are plotted (zero-error variants are excluded since log(0) is undefined).

**How to interpret:** The naive violins should be wider and shifted upward (larger errors), especially in the extreme scenario. AFQuery violins should be narrow and close to the bottom (small errors). If AFQuery's violin is empty or nearly so, it means most variants have zero or near-zero error -- this is the expected result, confirming that technology-aware AN eliminates the distortion.

### Figure 2: AF Error by Kit Coverage (`fig_capkit_error_by_coverage`)

Three subplots (one per scenario), each a grouped bar chart. X-axis: 1 kit, 2 kits, 3 kits. Bars: naive MAE vs. AFQuery MAE.

**How to interpret:** Errors should decrease with more kit coverage. At 3-kit positions (covered by all kits), both methods should show near-zero error because AN_wes = AN_wgs when all kits cover a position. At 1-kit positions, the naive error should be largest because only a fraction of samples are covered but the naive method assumes all are. AFQuery should maintain low error across all coverage levels.

### Figure 3: AF Scatter (`fig_capkit_scatter`)

2x3 grid of scatter plots. Rows: AFQuery (top) and naive (bottom). Columns: balanced, skewed, extreme. Each point is a variant colored by kit coverage (red=1, yellow=2, green=3). Diagonal dashed line = perfect agreement.

**How to interpret:** AFQuery row should show tight clustering around the diagonal regardless of scenario. Naive row should show increasing scatter (points below the diagonal, since AF_naive < AF_wgs due to AN inflation) as kit imbalance increases. Red points (1-kit positions) in the naive row should deviate most from the diagonal.

### Figure 4: ACMG Misclassification (`fig_capkit_acmg`)

Two subplots (cardiomyopathy and metabolic). Each shows grouped bars: naive vs. AFQuery discordance rate (%) per scenario. Red annotations on naive bars show false pathogenic (FP) and missed pathogenic (MP) counts.

**How to interpret:** AFQuery bars should be near zero (correct classification). Naive bars should grow with scenario severity. FP and MP counts are the clinically critical numbers -- even a small number of false pathogenic calls can trigger invasive follow-up. The metabolic model (stricter thresholds) should show more discordance than cardiomyopathy because the PM2 threshold (AC=0) is extremely sensitive to even small AN changes.

### Figure 5: Toward-Pathogenic Misclassification Rate (`fig_capkit_toward_pathogenic`)

Two subplots (cardiomyopathy and metabolic). Grouped bars per scenario showing the percentage of variants misclassified in the toward-pathogenic direction (predicted class is more pathogenic-supporting than truth), for AFQuery vs. Naive.

**How to interpret:** Naive's bars should dominate — AF deflation from AN inflation systematically converts benign evidence (BA1, BS1) into weaker or absent evidence (BS1→neutral, BA1→BS1). AFQuery's bars should be near zero because technology-aware AN eliminates this deflation. The gap between methods widens with kit imbalance (extreme scenario) and with stricter disease thresholds (metabolic). This figure provides the key clinical argument: Naive's errors are not random but systematically in the direction that increases VUS burden and misses benign evidence.

### Figure 6: AN Inflation Histogram (`fig_capkit_an_ratio`)

Three subplots (one per scenario). Histogram of AN_naive / AN_AFQuery per variant. Green dashed line at 1.0 (no inflation). Red dotted line at median.

**How to interpret:** The distribution should be right-skewed (ratios > 1). The median ratio directly measures the "average AN inflation." In the extreme scenario, the median might reach 1.5-2x, meaning the naive approach systematically overcounts alleles by 50-100%. Variants at ratio = 1.0 are at 3-kit positions (no inflation). The spread of the distribution reflects the heterogeneity of kit coverage across the exome.

---

## Output

- `results/merged.parquet` -- merged AF comparison table (all scenarios)
- `results/nadi_summary.json` -- summary statistics per scenario
- `results/acmg_results.json` -- misclassification counts per disease/scenario (includes `toward_pathogenic`, `toward_benign`, `bs1_recall` per method)
- `figures/fig_capkit_*.{pdf,png}` -- high-resolution plots

## Known Limitations

- **Chromosome scope:** Only chr22 is analyzed. Kit coverage overlap patterns on this small autosome (~51 Mb) may not generalize genome-wide.
- **Kit vendor:** Only Agilent SureSelect kits (v5, v6, v7) are tested. Other vendors (Illumina TruSeq, Roche SeqCap) and custom panels may exhibit different overlap characteristics.
- **Scenario design:** All scenarios include all three kits. The case where a kit is entirely absent (e.g., only 2 kits in the cohort) is not modeled.

## References

- Whiffin N et al. (2017). Cardiovascular disease risk gene burden in the founding population of the Amish. *Human Molecular Genetics* 26(17):3423-3428. DOI: 10.1093/hmg/ddx265
- Kalia SS et al. (2017). Recommendations for reporting of secondary findings in clinical exome and genome sequencing, 2016 update. *Genetics in Medicine* 19(2):249-255. DOI: 10.1038/gim.2016.190
- Agilent SureDesign: https://earray.chem.agilent.com/suredesign/
