# Installation Guide for AFQuery Benchmarks

This guide covers installing benchmark dependencies in your existing `snakemake` micromamba environment.

## Prerequisites

You must already have:
- Snakemake ≥ 8 installed in a micromamba environment named `snakemake`

If not, install it:
```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba install -n snakemake snakemake
```

## Install Benchmark Dependencies

Run this command to install all required tools and Python packages:

```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake mamba install \
  -c bioconda -c conda-forge \
  bcftools bedtools matplotlib numpy pandas scipy wget \
  && /mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake \
  pip install -e /mnt/lustre/home/dlopez/projects/afquery
```

## Step-by-step (alternative)

If you prefer to install step-by-step:

### 1. Bioinformatics tools

```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake mamba install \
  -c bioconda -c conda-forge \
  bcftools bedtools wget
```

### 2. Python scientific libraries

```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake mamba install \
  -c conda-forge \
  numpy pandas scipy matplotlib
```

### 3. AFQuery and dependencies

```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake \
  pip install -e /mnt/lustre/home/dlopez/projects/afquery
```

## Verification

Test that everything is installed:

```bash
# Test bcftools
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake bcftools --version

# Test Python imports
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake python -c "
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
bcftools 1.18 (using htslib 1.18)
All dependencies OK!
```

## Running Benchmarks

Once dependencies are installed, run benchmarks as documented in `README.md`:

```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake \
  snakemake --cores 52 all
```

## Troubleshooting

### `afquery` import fails

If `import afquery` fails after `pip install -e`:
- Verify the afquery source exists at `/mnt/lustre/home/dlopez/projects/afquery`
- Check that `src/afquery/__init__.py` exists
- Verify the install ran without errors (check for `Successfully installed afquery`)

### Missing Python packages

If any import fails (numpy, pandas, etc.):
```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake pip list | grep -E "numpy|pandas|scipy|matplotlib"
```

Re-run the mamba install commands above if any are missing.

### `bcftools` not found

Verify bcftools installation:
```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake which bcftools
```

If not found, re-run:
```bash
/mnt/lustre/home/dlopez/.local/bin/micromamba run -n snakemake mamba install -c bioconda bcftools
```
