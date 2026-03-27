# AFQuery Benchmarking Suite

This directory contains benchmarking experiments for AFQuery.

## Benchmark Categories

### 1. [Performance Benchmarks](./performance/README.md)

Core performance characterization across 4 experiments:

- **Experiment 1:** Query latency scaling with sample count
- **Experiment 2:** Build time scaling with parallelism
- **Experiment 3:** VCF annotation throughput
- **Experiment 4:** AFQuery vs. bcftools comparison

Uses real 1000 Genomes data (chr22) and synthetic datasets (1K-50K samples).

### 2. [Capture Kit Benchmark](./capture_kit/README.md)

Capture kit mixing impact on allele frequency classification:

- Sample generation with three Agilent SureSelect kits (v5, v6, v7)
- Three mixing scenarios (balanced, skewed, extreme)
- ACMG classification discordance analysis
- Coverage overlap metrics

## Prerequisites

- micromamba with `snakemake` environment already set up (with Snakemake >= 8 and SLURM executor plugin)
- Benchmark dependencies installed in the `snakemake` environment (see Environment Setup section)
- `/usr/bin/time` (GNU time) for memory profiling -- should be available on most Linux systems
- ~200 GB disk space for all 1KG data and databases

## Environment Setup

All benchmark dependencies are installed in the `snakemake` micromamba environment.

### Install dependencies

The benchmark requires:
- External bioinformatics tools: bcftools >= 1.18, bedtools, bgzip/tabix
- Python packages: pandas, numpy, scipy, matplotlib
- AFQuery and its dependencies: pyroaring, pyarrow, duckdb, cyvcf2, pyranges, click, tqdm

Install them in your `snakemake` environment:

```bash
micromamba run -n snakemake mamba install \
  -c bioconda -c conda-forge \
  bcftools bedtools matplotlib numpy pandas scipy wget \
  && micromamba run -n snakemake \
  pip install -e <path_to_afquery_repo>
```

### Verify

```bash
micromamba run -n snakemake bcftools --version
micromamba run -n snakemake python -c "import afquery; print('OK')"
```

## Running the Benchmarks

All commands below are run from the `benchmarks/` directory. Activate the `snakemake` environment first:

```bash
micromamba activate snakemake
```

Then use Snakemake with `--cores 52` to utilize all local CPU cores:

```bash
# Run everything (both benchmarks)
snakemake --cores 52 all

# Run individual benchmarks
snakemake --cores 52 performance_all
snakemake --cores 52 capture_kit_all

# Download 1KG data only (prerequisite for both)
snakemake --cores 52 download_1kg

# Dry run (preview what will execute)
snakemake --cores 52 --dry-run all

# Smoke test (fast validation with small parameter scales)
snakemake --cores 52 --config smoke_test=true all
```

### Smoke Test

The `smoke_test=true` flag reduces both pipelines to minimal scales for fast validation:

- **Performance:** 1 synthetic scale (1K), 1 1KG subset (500), 2 thread counts (1, 4), 1 rep
- **Capture kit:** 50 samples, balanced scenario only

This is useful for verifying pipeline logic without waiting for full-scale runs.

### Resuming after a failure

Snakemake uses output files to track completed steps. Re-running any command above will automatically skip steps whose outputs already exist and resume from the first incomplete step.

## Directory Structure

```
benchmarks/
├── Snakefile              # Root pipeline (includes both benchmarks)
├── config.yaml            # Global parameters (data_dir, threads, etc.)
├── envs/
│   └── benchmark.yaml     # conda/micromamba environment spec
├── shared/
│   ├── __init__.py
│   ├── config.py          # Common constants: DATA_DIR, 1KG paths, SEED
│   ├── utils.py           # Common helpers: stats, time_ms, save_figure, WONG_COLORS
│   └── rules/
│       └── download_1kg.smk   # Shared Snakemake rules: download + split 1KG
├── performance/
│   ├── Snakefile          # Performance-specific rules
│   ├── config.py          # Performance parameters (scales, reps, thread counts)
│   ├── config_smoke.py    # Minimal config for quick smoke testing
│   ├── 01_prepare_data.py
│   ├── 02_query_scaling.py
│   ├── 03_build.py
│   ├── 04_annotate.py
│   ├── 05_vs_bcftools.py
│   ├── 06_plot.py
│   ├── collect_prepare.py        # Aggregation scripts
│   ├── collect_query_scaling.py
│   ├── collect_build_perf.py
│   ├── collect_annotate.py
│   └── collect_bcftools.py
└── capture_kit/
    ├── Snakefile          # Capture kit-specific rules
    ├── config.py          # Capture kit parameters (scenarios, ACMG thresholds)
    ├── beds/              # Agilent SureSelect BED files (chr22, committed)
    │   ├── SureSelect_v5.bed
    │   ├── SureSelect_v6.bed
    │   └── SureSelect_v7.bed
    ├── 01a_assign_samples.py    # Subsample & assign technologies
    ├── 01b_write_manifests.py   # Write AFQuery manifest TSVs
    ├── 02_build_databases.py
    ├── 03_compute_metrics.py
    ├── 04_classify_acmg.py
    └── 05_plot_figures.py
```

## Configuration

Edit `config.yaml` to set the data directory and global parameters:

```yaml
data_dir: "/path/to/bench_data"   # needs ~200 GB free
build_memory: "8GB"               # DuckDB memory per worker
```

The data directory can also be set via the `AFQUERY_BENCH_DATA` environment variable.

## Known Limitations

- **Chromosome scope:** All benchmarks use only chromosome 22 (1000 Genomes Phase 3). Chr22 is one of the smallest autosomes (~51 Mb) and coverage patterns may not be representative of genome-wide behavior. The paper's methods section should discuss generalizability.
- **Synthetic data for scaling:** Query scaling experiments (1K-50K samples) use synthetic data. Allele frequency spectra and LD patterns differ from real population data.
- **Capture kits:** The capture kit benchmark uses only Agilent SureSelect kits (v5, v6, v7). Other vendors (Illumina, Roche) and custom panels are not evaluated.

## Reproducibility

All random seeds are fixed in `shared/config.py` (SEED = 42).
Re-running with the same configuration produces identical results within timing noise.

## References

- Weber LM et al. (2019). Essential guidelines for computational method benchmarking. *Genome Biology* 20:125. DOI: 10.1186/s13059-019-1738-8
- Auton A et al. (2015). A global reference for human genetic variation. *Nature* 526:68-74. DOI: 10.1038/nature15393
