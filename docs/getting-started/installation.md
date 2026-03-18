# Installation

## Requirements

- Python 3.10 or later
- For the build phase (`create-db`): 2–8 GB RAM per worker (see [Performance Tuning](../advanced/performance.md))

---

## pip (recommended)

```bash
pip install afquery
```

To install the latest development version directly from GitHub:

```bash
pip install git+https://github.com/dlopez-bioinfo/afquery.git
```

---

## conda / mamba (bioconda)

```bash
conda install -c bioconda -c conda-forge afquery
```

Or with mamba (faster):

```bash
mamba install -c bioconda -c conda-forge afquery
```

!!! note
    The bioconda package includes all compiled dependencies (cyvcf2, pyarrow, duckdb, pyroaring).

---

## From Source

```bash
git clone https://github.com/dlopez-bioinfo/afquery.git
cd afquery
pip install -e .
```

To also install documentation dependencies:

```bash
pip install -e ".[docs]"
```

---

## Verify Installation

```bash
afquery --version
```

You should see output like:

```
afquery, version 1.2.0
```

---

## Optional: Documentation Dependencies

If you want to build or serve the documentation locally:

```bash
pip install -e ".[docs]"
mkdocs serve
```
