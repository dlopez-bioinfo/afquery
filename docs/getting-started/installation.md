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

## Docker

Official images are published to [GitHub Container Registry](https://github.com/dlopez-bioinfo/afquery/pkgs/container/afquery) for every release.

```bash
docker pull ghcr.io/dlopez-bioinfo/afquery:latest
```

Run the CLI by mounting a local directory with your database:

```bash
docker run --rm \
  -v /path/to/db:/db \
  ghcr.io/dlopez-bioinfo/afquery:latest \
  query --db /db --locus chr1:925952 --phenotype E11.9
```

!!! note
    Images are built for `linux/amd64`. The `latest` tag always points to the most recent stable release.

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
