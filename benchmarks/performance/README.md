# Performance Benchmarks

Characterize AFQuery query latency, build scalability, annotation throughput, and competitive performance.

## Experiments

| # | Script | Question |
|---|--------|----------|
| 1 | `02_query_scaling.py` | How does query latency scale with sample count? |
| 2 | `03_build.py` | How does build time scale with parallelism? |
| 3 | `04_annotate.py` | What is the VCF annotation throughput? |
| 4 | `05_vs_bcftools.py` | How does AFQuery compare to bcftools? |

## Running

All experiments are orchestrated by Snakemake from the `benchmarks/` root.

```bash
# Run the full performance benchmark (HPC)
snakemake --profile benchmarks/profiles/slurm performance_all

# Run locally
snakemake --cores 52 performance_all

# Run a single experiment (e.g. query scaling only)
snakemake --cores 4 benchmarks/performance/results/query_scaling.json

# Dry run
snakemake --profile benchmarks/profiles/slurm --dry-run performance_all
```

The scripts can also be run directly when working from this directory:

```bash
# 1KG data must already be downloaded (snakemake download_1kg, or run from benchmarks/)
python 01_prepare_data.py
python 02_query_scaling.py
python 03_build.py
python 04_annotate.py
python 05_vs_bcftools.py
python 06_plot.py
```

## Configuration

Experiment parameters are in `config.py` (imports shared constants from `../shared/config.py`):

```python
# Synthetic scaling experiment
SYNTH_SCALES = [1_000, 5_000, 10_000, 25_000, 50_000]

# 1KG subsets
ONEKG_SUBSETS = [500, 1_000, 2_504]

# Experiment repetitions
QUERY_WARM_REPS = 50
BUILD_THREAD_COUNTS = [1, 4, 8, 16, 32]
ANNOTATE_THREAD_COUNTS = [1, 4, 8, 16, 32]
BCFTOOLS_REPS = 10
```

For a quick smoke test (tiny scales, few reps):

```bash
# From the performance/ directory
PYTHONPATH=.. python -c "import config_smoke"
```

## Data

Not included in the repository. Downloaded and generated automatically by the pipeline:

- **1000 Genomes Phase 3 chr22:** Downloaded from EBI FTP (~200 MB compressed)
- **Synthetic data:** Generated via `afquery.preprocess.synth` at scales 1K–50K

## Output

- `results/*.json` — raw timing data (gitignored)
- `figures/*.{pdf,png}` — publication-quality figures (gitignored)

## Hardware

Document hardware when reporting results:

- **CPU:** (e.g., 2× Intel Xeon Gold 6230R, 52 cores total)
- **RAM:** (e.g., 370 GB DDR4)
- **Storage:** (e.g., Lustre parallel filesystem)
- **OS:** (e.g., RHEL 8.5, kernel 4.18)

## References

- Weber LM et al. (2019). Essential guidelines for computational method benchmarking. *Genome Biology* 20:125. DOI: 10.1186/s13059-019-1738-8
- Auton A et al. (2015). A global reference for human genetic variation. *Nature* 526:68-74. DOI: 10.1038/nature15393
- 1000 Genomes Phase 3 FTP: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
