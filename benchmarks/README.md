# AFQuery Benchmarking Suite

This directory contains benchmarking experiments for the AFQuery publication.

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

## Prerequisites (Both)

- AFQuery installed in development mode: `pip install -e ".[dev,docs]"`
- `bcftools` available in `$PATH` (v1.17+)
- `/usr/bin/time` (GNU time)
- ~50 GB disk space for all data

## Quick Navigation

| Benchmark | Purpose | Time | Entry Point |
|-----------|---------|------|-------------|
| **Performance** | AFQuery performance characterization | Hours–days | `performance/README.md` |
| **Capture Kit** | Kit mixing impact on variant interpretation | Hours | `capture_kit/README.md` |

## Output Structure

Each benchmark maintains its own results and figures:

```
benchmarks/
├── performance/
│   ├── results/           # JSON timing data (gitignored)
│   └── figures/           # Publication-quality PDFs (gitignored)
└── capture_kit/
    ├── results/           # Analysis CSVs (gitignored)
    └── figures/           # Plots (gitignored)
```

## Reproducibility

All random seeds are fixed in respective `config.py` files (default: 42).
Re-running with the same configuration should produce identical results within timing noise.

## References

- Weber LM et al. (2019). Essential guidelines for computational method benchmarking. *Genome Biology* 20:125. DOI: 10.1186/s13059-019-1738-8
- Auton A et al. (2015). A global reference for human genetic variation. *Nature* 526:68-74. DOI: 10.1038/nature15393
