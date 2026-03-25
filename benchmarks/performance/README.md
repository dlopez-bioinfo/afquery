# Performance Benchmarks

Characterize AFQuery query latency, build scalability, annotation throughput, and competitive performance.

## Prerequisites

- AFQuery installed in development mode: `pip install -e ".[dev,docs]"`
- `bcftools` available in `$PATH` (tested with v1.17+)
- `/usr/bin/time` (GNU time, not shell builtin)
- Sufficient disk space (~50 GB for all experiments)

## Hardware

Document the hardware when reporting results:

- **CPU:** (e.g., 2× Intel Xeon Gold 6230R, 52 cores total)
- **RAM:** (e.g., 384 GB DDR4)
- **Storage:** (e.g., Lustre parallel filesystem)
- **OS:** (e.g., RHEL 8.5, kernel 4.18)

## Quick Start

```bash
# 1. Edit config.py to set DATA_DIR to a directory with enough space
#    Default: /mnt/lustre/home/dlopez/projects/afquery_bench_data

# 2. Download and prepare 1000 Genomes chr22 data
bash 00_download_1kg.sh

# 3. Build AFQuery databases (1KG subsets + synthetic at all scales)
python 01_prepare_data.py

# 4. Run experiments (each is independent, can run in any order)
python 02_query_scaling.py     # Exp 1: query latency vs. samples
python 03_build.py             # Exp 2: build time vs. threads
python 04_annotate.py          # Exp 3: annotation throughput
python 05_vs_bcftools.py       # Exp 4: AFQuery vs. bcftools

# 5. Generate all figures
python 06_plot.py
```

## Experiments

| # | Script | Question |
|---|--------|----------|
| 1 | `02_query_scaling.py` | How does query latency scale with sample count? |
| 2 | `03_build.py` | How does build time scale with parallelism? |
| 3 | `04_annotate.py` | What is the VCF annotation throughput? |
| 4 | `05_vs_bcftools.py` | How does AFQuery compare to bcftools? |

## Data

Data files are **not** included in the repository. Scripts download and generate everything:

- **1000 Genomes Phase 3 chr22:** Downloaded from EBI FTP (~200 MB compressed)
- **Synthetic data:** Generated via `afquery.preprocess.synth` at scales 1K–50K

## Output

- `results/*.json` — raw timing data (gitignored)
- `figures/*.{pdf,png}` — publication-quality figures (gitignored)

## Reproducibility

All random seeds are fixed in `config.py` (default: 42). Re-running with the same configuration should produce identical results within timing noise.

## Configuration

Edit `config.py` to customize:

```python
DATA_DIR = Path("/mnt/lustre/home/dlopez/projects/afquery_bench_data")

# 1000 Genomes subsets
ONEKG_SUBSETS = [500, 1_000, 2_504]

# Synthetic scaling
SYNTH_SCALES = [1_000, 5_000, 10_000, 25_000, 50_000]

# Experiment parameters
QUERY_WARM_REPS = 50
BUILD_THREAD_COUNTS = [1, 4, 8, 16, 32]
ANNOTATE_THREAD_COUNTS = [1, 4, 8, 16]
BCFTOOLS_REPS = 10
```

## References

- Weber LM et al. (2019). Essential guidelines for computational method benchmarking. *Genome Biology* 20:125. DOI: 10.1186/s13059-019-1738-8
- Auton A et al. (2015). A global reference for human genetic variation. *Nature* 526:68-74. DOI: 10.1038/nature15393
- 1000 Genomes Phase 3 FTP: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
