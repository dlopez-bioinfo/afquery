# Scientific Overview

This document provides a narrative description of AFQuery suitable for use in scientific publications, supplementary methods, and reviewer responses.

---

## Problem Statement

Large-scale genomic cohort studies require fast, flexible computation of allele frequencies across dynamically defined sample subsets. Existing tools for allele frequency computation — including bcftools stats, VCFtools, and GATK — operate on static VCF files and require reprocessing the entire cohort when the sample set or filter criteria change. This design is adequate for one-time analyses but is prohibitive for interactive clinical variant interpretation, where a researcher may need to compute AF over dozens of different subsets (by sex, phenotype, technology, or arbitrary combinations) in a single session.

A second limitation of general-purpose VCF tools is their lack of metadata-aware filtering. Computing AF over a subgroup — for example, male WGS samples without a specific diagnosis — requires pre-selecting samples by external metadata, subsetting the VCF, and then running AF computation. This multi-step process precludes real-time exploratory analysis.

---

## Approach

AFQuery introduces a pre-indexed database architecture that separates the slow step (building the genotype index from raw VCFs) from the fast step (querying allele frequencies on arbitrary subcohorts). The key data structure is the **Roaring Bitmap**, a compressed bitset that represents, for each variant, the set of samples carrying the alternate allele. At query time, computing AC/AN/AF requires only:

1. Loading the relevant bitmaps from Parquet storage
2. Intersecting with the bitmap of eligible samples (determined by sex, metadata, and capture filters)
3. Summing bits (popcount)

This reduces the per-query work to microsecond-scale bitmap operations, achieving sub-100 ms end-to-end latency including Parquet I/O.

---

## Novelty

AFQuery addresses the following methodological gaps not covered by existing tools:

### 1. Dynamic subcohort queries at sub-100 ms latency

Existing tools (bcftools stats, VCFtools --freq) scan VCF files linearly, scaling with file size. AFQuery queries execute in under 100 ms regardless of cohort size (tested up to 50,000 samples and 5M variants per chromosome), because queries access only the bitmaps for the relevant position rather than scanning the full dataset.

### 2. Incremental database updates without reprocessing

When new samples are added to the cohort, AFQuery merges new genotype data into the existing bitmap index without rebuilding from scratch. This enables real-time cohort growth in clinical settings where samples are added continuously.

### 3. Multi-dimensional metadata filtering

AFQuery supports simultaneous filtering by sex, arbitrary metadata codes, and sequencing technology, with both inclusion and exclusion semantics (`^CODE` prefix). Metadata codes are arbitrary strings defined by the researcher — ICD-10 codes, HPO terms, project tags, or any user-defined labels. No controlled vocabulary is required.

### 4. Server-less, portable database format

The AFQuery database is a directory of standard Parquet files with a SQLite metadata database. It requires no server process, can be shared by copying, and can be queried from any machine with AFQuery installed. This is particularly valuable for clinical and research settings where infrastructure deployment is constrained.

### 5. Ploidy-aware AF for sex chromosomes

AFQuery correctly handles chrX (PAR/non-PAR), chrY, and chrMT ploidy rules per sample. Males at chrX non-PAR contribute AN=1; females contribute AN=2. This ensures accurate hemizygous frequency computation for X-linked variant analysis without manual adjustment.

---

## Methods Summary

### Database Construction

1. **Ingest phase**: Each single-sample VCF is parsed by a dedicated worker process using cyvcf2. Per-sample genotype data is written to a temporary SQLite database.

2. **Build phase**: DuckDB aggregates per-sample rows into per-variant bitmap rows. Processing is parallelized across 1-Mbp genomic buckets; all buckets are distributed to a ProcessPoolExecutor pool, enabling near-linear scaling with CPU count. Each bucket is written as a Parquet file with Roaring Bitmap-serialized `het_bitmap`, `hom_bitmap` (and `fail_bitmap` in schema v2) columns.

3. **Metadata storage**: Sample metadata (sex, technology, phenotype codes) and precomputed sample-group bitmaps are stored in SQLite for fast retrieval at query time.

### Query Execution

1. **Filter resolution**: The eligible sample set is computed by ANDing precomputed bitmaps for sex, metadata, and technology filters.

2. **Parquet scan**: DuckDB reads the relevant bucket file(s) and returns rows matching the query position.

3. **Bitmap operations**: For each matching row, the variant's `het_bitmap` and `hom_bitmap` are ANDed with the eligible sample bitmap. AC = popcount(het_bitmap) + 2 × popcount(hom_bitmap). AN is computed per the chromosome-specific ploidy rules.

### Incremental Updates

New samples are ingested and built using the same pipeline. The new variant bitmaps are merged into existing Parquet files using a DuckDB bitwise OR operation. Removed samples are cleared by bitwise AND NOT with a removal mask; compaction rewrites Parquet files to remove dead bits.

---

## Validation

### Consistency

Rebuilding the same cohort from scratch produces bitwise-identical AC/AN/AF results. Incremental updates (build with N samples, add M samples) produce results identical to a full build with N+M samples.

### Ploidy correctness

chrX non-PAR AF for male-only queries equals AC/AN where AN = number of male samples (not 2×). This was verified by comparing manual counts from raw VCFs against AFQuery output.

### Subcohort reproducibility

Phenotype-filtered queries return bitwise-identical results regardless of query order, confirming correct bitmap isolation.

---

## Comparison with Related Tools

| Tool | Architecture | Subcohort queries | Incremental updates | Latency |
|------|-------------|------------------|--------------------|---------|
| AFQuery | Bitmap index + Parquet | Yes, dynamic | Yes, without rebuild | <100 ms |
| bcftools stats | VCF scan | No (requires pre-selection) | No | Minutes |
| VCFtools --freq | VCF scan | No | No | Minutes |
| GATK GenomicsDB | Tile-based DB | Partial (pre-defined strata) | Yes (import) | Seconds |
| Hail | Spark-based | Yes, programmatic | Rebuild required | Seconds–minutes |

AFQuery is optimized for the specific use case of interactive, metadata-driven subcohort AF computation at clinical scale (10K–50K samples). It is not a general-purpose variant database; tools like GATK GenomicsDB or Hail are more appropriate for joint genotyping pipelines.

---

## Citing AFQuery

*(Citation information to be added upon publication.)*
