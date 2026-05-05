# Installation Guide for AFQuery Benchmarks

This guide covers installing benchmark dependencies in your existing `afquery_bench` micromamba environment.
See [README.md](./README.md) for the recommended setup using `envs/benchmark.yaml`.

## Option 1: Using conda environment file (Recommended)

The easiest way is to create the pre-configured environment:

```bash
cd benchmarks/
micromamba env create -f envs/benchmark.yaml
micromamba activate afquery_bench
```

Then skip to [Running Benchmarks](#running-benchmarks).

## Option 2: Manual installation

If you prefer to install into an existing environment, install dependencies into your `snakemake` micromamba environment:

### Prerequisites

You must have:
- Snakemake ≥ 8 installed in a micromamba environment named `snakemake`

If not, install it:
```bash
micromamba install -n snakemake snakemake
```

### Install all dependencies at once

```bash
micromamba run -n snakemake mamba install \
  -c bioconda -c conda-forge \
  bcftools bedtools matplotlib numpy pandas scipy wget \
  && micromamba run -n snakemake pip install afquery
```

### Or step-by-step installation

#### 1. Bioinformatics tools

```bash
micromamba run -n snakemake mamba install \
  -c bioconda -c conda-forge \
  bcftools bedtools wget
```

#### 2. Python scientific libraries

```bash
micromamba run -n snakemake mamba install \
  -c conda-forge \
  numpy pandas scipy matplotlib
```

#### 3. AFQuery

```bash
micromamba run -n snakemake pip install afquery
```

## Verification

Test that everything is installed:

```bash
# Test bcftools
micromamba run -n afquery_bench bcftools --version

# Test Python imports
micromamba run -n afquery_bench python -c "
import afquery
import numpy
import pandas
import scipy
import matplotlib
print('All dependencies OK!')
"
```

Expected output:
```
bcftools >=1.17
All dependencies OK!
```

## Running Benchmarks

Once dependencies are installed and the environment is active:

```bash
cd benchmarks/
micromamba activate afquery_bench
snakemake --cores 52 all
```

Or run with conda activation:

```bash
micromamba run -n afquery_bench snakemake --cores 52 all
```

See [README.md](./README.md) for more options (smoke test, dry-run, etc.).

## Troubleshooting

### `afquery` import fails

Verify the installation:
```bash
micromamba run -n afquery_bench pip show afquery
```

If not found, reinstall:
```bash
micromamba run -n afquery_bench pip install --upgrade afquery
```

### Missing Python packages

Check installed packages:
```bash
micromamba run -n afquery_bench pip list | grep -E "numpy|pandas|scipy|matplotlib"
```

If any are missing, re-run the mamba install commands above.

### `bcftools` not found

Verify bcftools is installed:
```bash
micromamba run -n afquery_bench which bcftools
```

If not found, install it:
```bash
micromamba run -n afquery_bench mamba install -c bioconda bcftools
```
