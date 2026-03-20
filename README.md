# AFQuery

[![CI](https://github.com/dlopez-bioinfo/afquery/actions/workflows/ci.yml/badge.svg)](https://github.com/dlopez-bioinfo/afquery/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/dlopez-bioinfo/afquery/graph/badge.svg)](https://codecov.io/gh/dlopez-bioinfo/afquery)
[![Docs](https://img.shields.io/website?url=https%3A%2F%2Fdlopez-bioinfo.github.io%2Fafquery%2F&label=docs)](https://dlopez-bioinfo.github.io/afquery/)
<br>
[![PyPI](https://img.shields.io/pypi/v/afquery.svg?color=blue)](https://pypi.org/project/afquery/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/afquery.svg)](https://anaconda.org/bioconda/afquery)
[![Python](https://img.shields.io/pypi/pyversions/afquery.svg)](https://pypi.org/project/afquery/)
[![License: MIT](https://img.shields.io/github/license/dlopez-bioinfo/afquery)](https://github.com/dlopez-bioinfo/afquery/blob/master/LICENSE)

AFQuery enables fast allele frequency queries on user-defined subsets of local genomic cohorts, without rescanning VCFs.

AFQuery is a bitmap-indexed engine that efficiently recomputes AC/AN/AF for dynamically defined subcohorts (e.g., by phenotype, sex, or sequencing technology), a common requirement in ACMG/AMP variant classification. It stores per-variant genotype data as Roaring Bitmaps in Parquet files and resolves sample filters into bitmaps that can be intersected in microseconds, enabling sub-100 ms queries on large cohorts. The system accounts for ploidy in sex chromosomes, adjusts AN based on sequencing technology, supports incremental updates, and runs locally using a file-based setup (Parquet + SQLite) without requiring server or cloud infrastructure.

[Full Documentation→](https://dlopez-bioinfo.github.io/afquery/)

## When to use AFQuery

- You need allele frequencies for phenotype or user-defined subcohorts  
- You work with mixed sequencing technologies or capture kits versions (WGS, WES, targeted panels)  
- You require fast, repeated queries without rescanning VCFs  
- You want a local, reproducible workflow without cloud or cluster dependencies  

## Features

- **Dynamic subcohort queries (<100 ms)** — bitmap intersections at query time; no VCF re-scan required  
- **Technology-aware** — avoids bias when mixing WGS, WES, and panels using different BED capture indexes  
- **Ploidy-aware** — correct handling of sex chromosomes (PAR/non-PAR, chrX, chrY)  
- **ACMG-compatible allele counting** — AC/AN/AF computed per standard definitions  
- **Flexible metadata filtering** — arbitrary labels (ICD-10, HPO, custom fields) with inclusion/exclusion rules 
- **Incremental updates** — add or remove samples and update metadata without rebuilding the database
- **VCF annotation** — annotate variants using subcohort-specific frequencies
- **FILTER/call quality tracking** — failed calls (FILTER!=PASS) tracked per variant and reported as N_FAIL
- **Batch and region queries** — query a single locus, a genomic region, or a list of variants from a file
- **Bulk CSV export** — export all variant frequencies with optional disaggregation by sex, technology, or phenotype
- **Audit changelog** — all database operations logged with timestamps and operator notes
- **Database validation** — integrity checks with scripted exit codes
- **Portable and serverless** — file-based system, no infrastructure required

## Performance

- Query latency: <100 ms (tested up to 50,000 samples)  
- Storage: ~2 bytes/sample/variant  
- Scales to millions of variants per chromosome  

## Comparison with Alternative Tools

| | AFQuery | bcftools | GATK GenomicsDB | Hail |
|---|---|---|---|---|
| Technology-aware AN | Yes | No | No | No |
| Metadata filtering | Arbitrary labels | No | No | Custom code |
| Ploidy-aware sex chromosomes | Yes | Manual | No | Manual |
| Dynamic subcohort queries | Yes | No | Limited | Requires code |
| FILTER/call quality tracking | Per variant (N_FAIL) | Manual | No | Manual |
| Incremental updates | Yes | No | Yes | No |
| Infrastructure required | None (file-based) | None | Java/server | Spark cluster |
| Query latency (50K samples) | <100 ms | ~5 min | <1 min | 1–2 min |

## Algorithm Overview

AFQuery pre-indexes per-variant genotype data as [Roaring Bitmaps](https://roaringbitmap.org/) stored in Parquet files. Each variant row holds three bitmaps: heterozygous carriers, homozygous alt carriers, and samples with FILTER!=PASS. Sample metadata (sex, phenotype, technology) is pre-serialized as bitmaps in SQLite.

At query time, the requested sample filter is resolved to a single candidate bitmap via bitmap intersections and differences — taking microseconds regardless of cohort size. For each variant, the candidate bitmap is intersected with the genotype bitmaps to compute AC/AN/AF. AN accounts for WES capture regions (via BED-indexed interval trees) and for ploidy on sex chromosomes (males are haploid on non-PAR chrX and chrY).

## Input Requirements

- VCF files: normalized and consistent with the selected genome build (GRCh37 or GRCh38)  
- Sample metadata: must include sex, sequencing technology, and any fields used for filtering (e.g., phenotype)  
- BED files (optional): define capture regions for each sequencing technology  

## Quick Start

Example workflow from raw VCFs to query, export, and annotation:

```bash
pip install afquery

# Build the database
afquery create-db --manifest samples.tsv --output-dir ./db/ --genome-build GRCh38

# Inspect the database
afquery info --db ./db/

# Query a single position, filtered to a phenotype
afquery query --db ./db/ --locus chr1:925952 --phenotype E11.9 --sex female

# Query a genomic region
afquery query --db ./db/ --region chr1:900000-1000000

# Export all variant frequencies to CSV
afquery dump --db ./db/ --output all_variants.csv

# Annotate a VCF with cohort frequencies
afquery annotate --db ./db/ --input patient.vcf --output annotated.vcf --threads 12

# Add new samples to an existing database
afquery update-db --db ./db/ --add-samples new_samples.tsv
```

## Documentation

- [Full Documentation](https://dlopez-bioinfo.github.io/afquery/)
- [Getting Started](https://dlopez-bioinfo.github.io/afquery/getting-started/)
- [Why local allele frequencies?](https://dlopez-bioinfo.github.io/afquery/getting-started/motivation/)
- [CLI Reference](https://dlopez-bioinfo.github.io/afquery/reference/cli/)

## Citation

If you use AFQuery, please cite:

> AFQuery: fast, metadata-aware allele frequency queries on local genomic cohorts.  
> *(manuscript in preparation)*
