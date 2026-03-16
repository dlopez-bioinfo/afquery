# AFQuery

[![CI](https://github.com/dlopez-bioinfo/afquery/actions/workflows/ci.yml/badge.svg)](https://github.com/dlopez-bioinfo/afquery/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/dlopez-bioinfo/afquery/graph/badge.svg)](https://codecov.io/gh/dlopez-bioinfo/afquery)
[![Docs](https://img.shields.io/website?url=https%3A%2F%2Fdlopez-bioinfo.github.io%2Fafquery%2F&label=docs)](https://dlopez-bioinfo.github.io/afquery/)
<br>
[![PyPI](https://img.shields.io/pypi/v/afquery.svg?color=blue)](https://pypi.org/project/afquery/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/afquery.svg)](https://anaconda.org/bioconda/afquery)
[![Python](https://img.shields.io/pypi/pyversions/afquery.svg)](https://pypi.org/project/afquery/)
[![License: MIT](https://img.shields.io/github/license/dlopez-bioinfo/afquery)](https://github.com/dlopez-bioinfo/afquery/blob/master/LICENSE)

Fast, file-based genomic allele frequency queries for large cohorts (10K–50K samples).
No server, no cloud — just files. Sub-100ms point queries with flexible filtering
by sex, phenotype (ICD codes), and sequencing technology.

## Quick Example

```bash
pip install afquery

# Build database from your VCFs
afquery create-db --manifest samples.tsv --output-dir ./db/ --genome-build GRCh38

# Query allele frequency
afquery query --db ./db/ --chrom chr1 --pos 123456 --phenotype E11.9 --sex female

# Annotate a VCF
afquery annotate --db ./db/ --input variants.vcf --output annotated.vcf
```

**[Full documentation →](https://dlopez-bioinfo.github.io/afquery/)**

## Features

- Sub-100ms point queries on 50K-sample cohorts
- Filter by sex, phenotype (ICD codes), and sequencing technology
- Bitmap-compressed storage (Roaring Bitmaps + Parquet)
- Incremental updates (add/remove samples without full rebuild)
- VCF annotation with custom sample subsets
- Ploidy-aware AN for sex chromosomes (X/Y/MT)
- Zero infrastructure — purely file-based

## Installation

```bash
pip install afquery
# or
conda install -c bioconda -c conda-forge afquery
```

Requires Python 3.10+.

## Quick Start

### 1. Create a Database

First, prepare a manifest TSV with sample metadata:

```tsv
sample_name	sex	tech_name	vcf_path	phenotype_codes
sample_1	male	wgs	vcfs/sample_1.vcf	E11.9,I10
sample_2	female	wes_kit_a	vcfs/sample_2.vcf	E11.9
sample_3	male	wgs	vcfs/sample_3.vcf	I10
```

**Key points about the manifest:**
- `tech_name`: Either `wgs` (case-insensitive) for whole genome, or a custom technology name
- `vcf_path`: Path to the single-sample VCF file (relative to manifest directory, or absolute)
- For WES/exome technologies, capture regions are loaded from `--bed-dir/{tech_name}.bed`
  - Example: `tech_name=wes_kit_a` → loads `beds/wes_kit_a.bed`

Organize your files:
```
project/
├── manifest.tsv
├── vcfs/
│   ├── sample_1.vcf
│   ├── sample_2.vcf
│   └── sample_3.vcf
└── beds/              # Required for non-WGS technologies
    └── wes_kit_a.bed  # BED file for WES kit A
```

Then build your database:

```bash
afquery create-db \
  --manifest manifest.tsv \
  --bed-dir ./beds/ \
  --output-dir ./my_db/ \
  --genome-build GRCh38
```

This creates:
- `my_db/manifest.json` — database metadata
- `my_db/metadata.sqlite` — samples, technologies, phenotype codes, precomputed bitmaps
- `my_db/variants/{chrom}.parquet` — variant data with encoded genotypes
- `my_db/capture/` — capture regions for each technology

### 2. Query Allele Frequencies

```bash
# Point query
afquery query --db my_db --chrom chr1 --pos 1000 --alt G

# Batch query (from file with columns: pos ref alt)
afquery query --db my_db --chrom chr1 --from-file positions.tsv --phenotype E11.9

# Region query
afquery query --db my_db --chrom chr1 --start 1000 --end 10000 --phenotype E11.9 --sex M
```

### 3. Annotate VCF Files

```bash
afquery annotate \
  --db my_db \
  --input input.vcf \
  --output annotated.vcf \
  --phenotype E11.9 \
  --tech WGS
```

Adds `AFQUERY_AC`, `AFQUERY_AN`, `AFQUERY_AF` and genotype fields to your VCF.

## Python API

```python
from afquery import Database

db = Database("/path/to/db")

# Single position query
# Automatically filters samples by: sex + phenotype codes + capture coverage
results = db.query(
    chrom="chr1",
    pos=1000,
    alt="G",
    phenotype=["E11.9"],
    sex="both"
)
for r in results:
    print(f"AC={r.AC}, AN={r.AN}, AF={r.AF}")

# Batch query (multi-variant)
results = db.query_batch(
    "chr1",
    variants=[(1500, "A", "T"), (3500, "G", "C")],
    phenotype=["E11.9"],
)

# Region query (genomic range)
results = db.query_region(
    chrom="chr1",
    start=1000,
    end=10000,
    phenotype=["E11.9", "I10"]
)

# Annotate VCF with allele frequencies
# Note: tech filters annotation to samples of that technology
db.annotate_vcf(
    input_vcf="input.vcf",
    output_vcf="annotated.vcf",
    phenotype=["E11.9"],
    tech=["wgs"]  # Only annotate using WGS samples
)
```

**How samples are filtered in queries**:
- Sex filter: `male`, `female`, or `both`
- phenotype filter: All codes must match
- Capture filter: Automatic—only samples whose tech's BED covers the position

## Database Structure

```
my_db/
├── manifest.json          # Metadata: genome_build, sample_count, schema_version
├── metadata.sqlite        # SQLite: samples, technologies, phenotype codes, bitmaps
├── variants/
│   ├── chr1.parquet
│   ├── chr2.parquet
│   └── ...
└── capture/
    ├── tech_0.pickle      # WGS capture region (always covered)
    └── tech_1.pickle      # WES kit capture region
```

Each variant row contains:
- `pos` — 1-based genomic position
- `ref` — reference allele
- `alt` — alternate allele
- `het_bitmap` — Roaring Bitmap of heterozygous samples
- `hom_bitmap` — Roaring Bitmap of homozygous samples

## How Capture BED Files Are Associated with Samples

Samples are linked to capture regions through their **technology**:

1. **Manifest specifies technology**: Each sample lists `tech_name` (e.g., `wgs`, `wes_kit_a`)
2. **Technology maps to BED file**:
   - **WGS**: No BED file needed (always fully covered)
   - **WES/Custom**: BED file loaded from `{bed_dir}/{tech_name}.bed`
3. **Storage in database**:
   - `metadata.sqlite::technologies` stores tech_id, tech_name, and bed_path
   - `metadata.sqlite::samples` stores sample_id, sample_name, and tech_id (foreign key)
4. **Query-time filtering**: When querying, samples are filtered by:
   - Sex (male/female/both)
   - phenotype diagnosis codes
   - Capture region coverage (via tech's BED file)

**Example**: If you have samples on two exome kits:
```tsv
sample_name	sex	tech_name	vcf_path	phenotype_codes
S001	male	exome_v1	vcfs/S001.vcf	E11.9
S002	female	exome_v1	vcfs/S002.vcf	E11.9
S003	male	exome_v2	vcfs/S003.vcf	I10
```

Then provide:
```
beds/
├── exome_v1.bed    # Coverage for samples S001, S002
└── exome_v2.bed    # Coverage for sample S003
```

At query time, each sample's eligible regions are determined by its tech's BED file.

## Advanced Features

### Incremental Updates

```bash
afquery update-db \
  --db my_db \
  --add-samples new_samples.tsv
```

Adds new samples without rebuilding the entire database.

### Remove Samples and Compact

Remove samples and reclaim disk space:

```bash
afquery update-db --db my_db --remove-samples sample_1,sample_2 --compact
```

### Run Benchmarks

```bash
afquery benchmark --n-samples 5000 --n-variants 100000
```

## Ploidy Rules

AF computation respects chromosome-specific ploidy:

| Chromosome | Formula |
|---|---|
| Autosomes | `AN = 2 × eligible_samples` |
| chrM | `AN = 1 × eligible_samples` |
| chrY | `AN = 1 × eligible_males` |
| chrX (PAR) | `AN = 2 × eligible_samples` |
| chrX (non-PAR) | `AN = 2 × eligible_females + 1 × eligible_males` |

Where `eligible` = samples matching sex, phenotype, and technology capture filters.

## Performance Targets

- **Point query (cold)**: <100 ms
- **Point query (warm)**: ~10 ms
- **Batch 100 positions**: ~200 ms
- **VCF annotation (5K variants)**: ~30 s
- **VCF annotation (5M variants)**: ~30 min

## Command Reference

```
afquery query                Query one position, region, or batch (--from-file)
afquery annotate             Annotate VCF file with AF info fields
afquery dump                 Export allele frequencies for all variants to CSV
afquery info                 Show database metadata and sample list
afquery check                Validate database integrity
afquery create-db            Build database from a VCF manifest
afquery update-db            Add/remove samples or compact the database
afquery version              Show or set the database version label
afquery benchmark            Run performance benchmarks
```

Run `afquery --help` for full options.

## Development

### Running Tests

```bash
# All tests
python3 -m pytest --tb=short -q

# Specific test module
python3 -m pytest tests/test_query.py -v
```

### Key Modules

- `afquery.query` — QueryEngine, point/batch/region queries
- `afquery.annotate` — VCF annotation pipeline
- `afquery.database` — Database wrapper (public API)
- `afquery.preprocess` — Manifest parsing, VCF ingestion, Parquet building
- `afquery.bitmaps` — Roaring Bitmap encoding/decoding
- `afquery.ploidy` — Chromosome-specific ploidy rules
- `afquery.models` — Data classes (QueryResult, ParsedSample, etc.)


## Genome Builds

- **GRCh37** (hg19) — PAR regions: chrX:1-2649520
- **GRCh38** (hg38) — PAR regions: chrX:1-3099677

## Troubleshooting

### ImportError: `cyvcf2`

`cyvcf2` import happens inside worker processes during preprocessing. Do not import at module level.

### DuckDB Temp Files

Uses Parquet format (not Arrow IPC) for compatibility. Set `DUCKDB_TEMP_DIRECTORY` if needed.

## License

(Add your license here)

## Citation

If you use afquery in research, please cite:

```
(Citation format to be determined)
```

---
