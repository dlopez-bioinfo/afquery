# AFQuery Benchmarking Suite

This directory contains benchmarking experiments for AFQuery.

## Benchmark Categories

### 1. [Performance Benchmarks](./performance/README.md)

Core performance characterization across 4 experiments:

- **Experiment 1:** Query latency scaling with sample count
- **Experiment 2:** Build time scaling with parallelism
- **Experiment 3:** VCF annotation throughput
- **Experiment 4:** AFQuery vs. bcftools comparison

Uses real 1000 Genomes data (chr22) and synthetic datasets (1K–50K samples).

### 2. [Capture Kit Benchmark](./capture_kit/README.md)

Capture kit mixing impact on allele frequency classification:

- Sample generation with three Agilent SureSelect kits (v5, v6, v7)
- Three mixing scenarios (balanced, skewed, extreme)
- ACMG classification discordance analysis
- Coverage overlap metrics

## Prerequisites

- micromamba (or conda) to create the benchmark environment
- Snakemake ≥ 8 with SLURM executor plugin installed in your base environment:
  ```bash
  pip install snakemake snakemake-executor-plugin-slurm
  ```
- `/usr/bin/time` (GNU time) for memory profiling
- ~200 GB disk space for all 1KG data and databases

## Environment Setup

All benchmark tools and Python dependencies are declared in `envs/benchmark.yaml`.
afquery is installed from PyPI and covers most Python dependencies; the file adds
the external bioinformatics tools (bcftools ≥ 1.18, bedtools, bgzip/tabix)
and plotting libraries (pandas, numpy, matplotlib).

### Create the environment

```bash
micromamba env create -f envs/benchmark.yaml
```

### Verify

```bash
micromamba run -n afquery_bench bcftools --version
micromamba run -n afquery_bench python -c "import afquery; print('OK')"
```

Snakemake activates `afquery_bench` automatically for each job.
To use a different environment spec, change `conda_env_file` in `config.yaml`.

## Running the Benchmarks

All commands below are run from the `benchmarks/` directory.

### On an HPC with SLURM (recommended)

```bash
# Run everything (both benchmarks) across as many nodes as needed
# conda/micromamba support is already configured in profiles/slurm/config.yaml
snakemake --profile profiles/slurm all

# Run only one benchmark
snakemake --profile profiles/slurm performance_all
snakemake --profile profiles/slurm capture_kit_all

# Download 1KG data only (prerequisite for both)
snakemake --profile profiles/slurm download_1kg
```

### Locally (without SLURM)

```bash
snakemake --cores all --software-deployment-method conda --conda-frontend micromamba all
snakemake --cores all --software-deployment-method conda --conda-frontend micromamba performance_all
```

### Dry run (preview what will execute)

```bash
snakemake --profile profiles/slurm --dry-run all
```

### Resuming after a failure

Snakemake uses output files to track completed steps. Re-running any
command above will automatically skip steps whose outputs already exist
and resume from the first incomplete step.

## Directory Structure

```
benchmarks/
├── Snakefile              # Root pipeline (includes both benchmarks)
├── config.yaml            # Global parameters (data_dir, threads, conda_env_file, etc.)
├── envs/
│   └── benchmark.yaml     # conda/micromamba environment spec
├── profiles/
│   └── slurm/
│       └── config.yaml    # SLURM executor settings (resources per rule)
├── shared/
│   ├── config.py          # Common constants: DATA_DIR, 1KG paths, SEED
│   ├── utils.py           # Common helpers: stats, time_ms, save_figure, WONG_COLORS
│   ├── rules/
│   │   └── download_1kg.smk   # Shared Snakemake rules: download + split 1KG
│   └── scripts/
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
│   ├── results/           # JSON timing data
│   └── figures/           # High-resolution PDFs/PNGs
└── capture_kit/
    ├── Snakefile          # Capture kit-specific rules
    ├── config.py          # Capture kit parameters (scenarios, ACMG thresholds)
    ├── 01_prepare_samples.py
    ├── 02_build_databases.py
    ├── 03_compute_metrics.py
    ├── 04_classify_acmg.py
    ├── 05_plot_figures.py
    ├── results/           # Analysis parquet/JSON
    └── figures/           # Plots
```

## Configuration

Edit `config.yaml` to set the data directory and global parameters:

```yaml
data_dir: "/path/to/bench_data"   # needs ~200 GB free
build_memory: "8GB"               # DuckDB memory per worker
```

The data directory can also be set via the `AFQUERY_BENCH_DATA` environment variable.

## Reproducibility

All random seeds are fixed in `shared/config.py` (SEED = 42).
Re-running with the same configuration produces identical results within timing noise.

## References

- Weber LM et al. (2019). Essential guidelines for computational method benchmarking. *Genome Biology* 20:125. DOI: 10.1186/s13059-019-1738-8
- Auton A et al. (2015). A global reference for human genetic variation. *Nature* 526:68-74. DOI: 10.1038/nature15393
