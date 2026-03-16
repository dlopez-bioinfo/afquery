# AFQuery

**Fast, file-based genomic allele frequency queries for large cohorts. No server, no cloud — just files.**

AFQuery stores genotype data as Roaring Bitmaps in Parquet files and answers allele frequency queries in under 100 ms across 10K–50K samples, with flexible filtering by sex, phenotype (ICD codes), and sequencing technology.

---

## Features

- **Sub-100 ms point queries** on 50K-sample cohorts, ~10 ms warm
- **Filter by sex, phenotype (ICD codes), and sequencing technology** — any combination
- **Ploidy-aware AN** for sex chromosomes (chrX, chrY, chrMT)
- **Roaring Bitmap compression** — ~2 bytes/sample/variant typical storage
- **Incremental updates** — add or remove samples without full rebuild
- **VCF annotation** with custom sample subsets
- **Bulk export** disaggregated by sex, technology, or phenotype
- **Zero infrastructure** — purely file-based, no server required

---

## Architecture

```
  Input VCFs (single-sample)
        │
        ▼
  ┌─────────────┐
  │  Ingest      │  cyvcf2 reads VCFs → per-sample genotype rows
  │  (SQLite)    │  stored in SQLite temp DB
  └──────┬──────┘
         │
         ▼
  ┌─────────────┐
  │  Build       │  DuckDB aggregates per 1-Mbp bucket →
  │  (Parquet)   │  Roaring Bitmaps per variant stored as Parquet
  └──────┬──────┘
         │
         ▼
  ┌─────────────────────────────┐
  │  Database on disk           │
  │  variants/chr1/bucket_0/    │  Hive-partitioned Parquet
  │  variants/chr2/bucket_1/    │
  │  capture/<tech>.pkl         │  Interval trees for WES coverage
  │  metadata.sqlite            │  Sample/phenotype/tech metadata
  │  manifest.json              │  Build configuration
  └──────┬──────────────────────┘
         │
         ▼
  ┌─────────────┐
  │  Query       │  Load bitmap → intersect with eligible samples →
  │  Engine      │  compute AC/AN/AF in microseconds
  └─────────────┘
```

---

## Quick Start

```bash
# 1. Install
pip install afquery

# 2. Build a database from your VCFs
afquery create-db \
  --manifest samples.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38

# 3. Query allele frequency
afquery query \
  --db ./db/ \
  --chrom chr1 \
  --pos 123456 \
  --phenotype E11.9 \
  --sex female

# 4. Annotate a VCF
afquery annotate \
  --db ./db/ \
  --input variants.vcf \
  --output annotated.vcf
```

!!! tip "Performance"
    Sub-100 ms cold point queries. ~10 ms warm. Scales to 50K samples with no infrastructure changes.

---

## Next Steps

- [Installation](getting-started/installation.md) — pip, conda, from source
- [Quickstart](getting-started/quickstart.md) — 5-minute end-to-end tutorial
- [Key Concepts](getting-started/concepts.md) — bitmaps, Parquet, manifest
