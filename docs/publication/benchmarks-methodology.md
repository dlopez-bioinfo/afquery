# Benchmarks Methodology

This document describes the reproducible benchmark protocol for AFQuery performance evaluation.

---

## Overview

Benchmarks measure:
1. **Query latency** — point, region, and batch queries at varying cohort sizes
2. **Build throughput** — database construction time at varying CPU counts
3. **Incremental update time** — adding samples to an existing database
4. **Memory footprint** — peak RAM during build and at query time

---

## Synthetic Cohort Generation

Benchmarks use the built-in synthetic data generator:

```bash
afquery benchmark --n-samples N --n-variants M --output benchmark_N_M.json
```

The generator creates:
- N samples with random sex (50% male, 50% female)
- M variants per chromosome with uniformly random positions across chr1–22, chrX, chrY, chrMT
- Random genotypes: carrier rate uniformly drawn from [0.001, 0.1] per variant
- Heterozygous:homozygous ratio fixed at 3:1 (reflecting typical population genetics)
- Random technology assignment: 60% WGS, 40% WES (with synthetic BED files)
- Random phenotype assignment: 10 labels, each sample assigned 0–3 labels

For real-database benchmarks:
```bash
afquery benchmark --db-dir ./production_db/ --output real_db_benchmark.json
```

---

## Query Latency Protocol

### Setup

```bash
# Warm up OS page cache (run once, not timed)
afquery query --db ./db/ --chrom chr1 --region 1-250000000 > /dev/null
```

### Measurements (100 iterations each)

```bash
# Point query — cold (flush OS cache between runs with `sync; echo 3 > /proc/sys/vm/drop_caches`)
for i in $(seq 100); do
  afquery query --db ./db/ --chrom chr1 --pos 925952 --format json > /dev/null
  # Flush cache between iterations for cold measurements
done

# Point query — warm (no cache flush)
for i in $(seq 100); do
  afquery query --db ./db/ --chrom chr1 --pos 925952 --format json > /dev/null
done

# Region query (1 Mbp)
afquery query --db ./db/ --chrom chr1 --region 900000-1900000 --format json > /dev/null

# Batch query (100 variants)
afquery query --db ./db/ --chrom chr1 --from-file 100_variants.tsv --format json > /dev/null
```

Report: median, p95, p99 latency across 100 iterations.

---

## Build Throughput Protocol

### Thread scaling (fixed cohort size)

```bash
for threads in 1 4 8 16 32 52; do
  time afquery create-db \
    --manifest manifest.tsv \
    --output-dir ./db_${threads}/ \
    --genome-build GRCh38 \
    --build-threads ${threads} \
    --build-memory 2GB \
    --force
done
```

Expected scaling table (50K samples, 2500 buckets, GRCh38):

| Cores | Build time | Speedup vs. 1 core |
|-------|-----------|-------------------|
| 1 | ~8 hours | 1× |
| 4 | ~2 hours | ~4× |
| 8 | ~1 hour | ~7× |
| 16 | ~30 min | ~14× |
| 32 | ~18 min | ~24× |
| 52 | ~13 min | ~38× |

*(Fill with measured values. Near-linear scaling up to ~32 cores; beyond that, I/O and SQLite contention limit further gains.)*

### Cohort size scaling (fixed thread count = 52)

```bash
for n_samples in 1000 5000 10000 50000; do
  time afquery create-db \
    --manifest manifest_${n_samples}.tsv \
    --output-dir ./db_${n_samples}/ \
    --genome-build GRCh38 \
    --build-threads 52 \
    --build-memory 2GB
done
```

---

## Incremental Update Protocol

```bash
# Baseline: build with N samples
afquery create-db --manifest baseline_N.tsv --output-dir ./db/ --genome-build GRCh38

# Time: add M new samples
time afquery update-db --db ./db/ --add-samples new_M.tsv --build-threads 52

# Verify: results identical to full rebuild
afquery create-db --manifest combined_N_plus_M.tsv --output-dir ./db_full/ --genome-build GRCh38
# Compare AF for 1000 random variants between ./db/ and ./db_full/
```

---

## Memory Footprint

### Build phase

Peak memory per worker = `--build-memory` setting. Total = `build_threads × build_memory`.

Monitor with:
```bash
/usr/bin/time -v afquery create-db --manifest manifest.tsv --output-dir ./db/ \
  --genome-build GRCh38 --build-threads 16 --build-memory 2GB 2> timing.txt
grep "Maximum resident" timing.txt
```

### Query phase

Query memory is very low. Measure with:
```python
import tracemalloc
from afquery import Database

tracemalloc.start()
db = Database("./db/")
results = db.query("chr1", pos=925952)
current, peak = tracemalloc.get_traced_memory()
print(f"Peak query memory: {peak / 1e6:.1f} MB")
tracemalloc.stop()
```

Typical: <500 MB regardless of cohort size.

---

## Comparison Protocol

To compare AFQuery against alternative approaches for allele frequency computation:

### bcftools stats (VCF scan baseline)

```bash
# Merge all samples into one VCF (required for bcftools)
bcftools merge sample_*.vcf.gz -Oz -o merged.vcf.gz

# Compute AF at specific positions
time bcftools stats --samples - --regions chr1:925952 merged.vcf.gz | grep "AF"
```

### VCFtools --freq

```bash
time vcftools --gzvcf merged.vcf.gz --freq --chr chr1 --from-bp 925952 --to-bp 925952
```

Report build time (preparing the merged VCF) + query time for AFQuery comparison.

---

## Reporting

All benchmark results should report:
- AFQuery version
- Hardware specification (CPU count, RAM, storage type)
- Cohort size (N samples, M variants)
- Genome build
- All query types with median ± SD and p95 latency
- Build time with thread scaling table
