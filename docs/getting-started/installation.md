# Installation

## Requirements

- Python 3.10 or later
- For the build phase (`create-db`): 2–8 GB RAM per worker (see [Performance Tuning](../advanced/performance.md))

---

## pip

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
conda install bioconda::afquery
```

Or with mamba (faster):

```bash
mamba install bioconda::afquery
```


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

## Optional: Documentation Dependencies

If you want to build or serve the documentation locally:

```bash
pip install -e ".[docs]"
mkdocs serve
```
