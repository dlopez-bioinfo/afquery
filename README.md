# afquery

Genomic allele frequency query engine with bitmap-encoded genotypes. Fast, file-based queries over 10K-50K samples at sub-100ms latency.

## Features

- **Fast point queries**: <100ms cold start, ~10ms warm queries on single positions
- **Batch queries**: Multi-position queries via SQL IN clauses or temporary tables
- **Region queries**: Genomic range queries with automatic partitioning
- **VCF annotation**: Annotate VCF files with computed allele frequencies and sample genotypes
- **Ploidy-aware**: Correct AN/AC computation for autosomes, chrX, chrY, and chrM
- **Bitmap compression**: Roaring Bitmaps for efficient genotype storage
- **In-process**: No server process—queries run locally with DuckDB
- **Incremental updates**: Add new samples to existing databases
- **Multiple genome builds**: Support for GRCh37 and GRCh38

## Installation

```bash
pip install -e /home/dan/projects/afquery/ --break-system-packages
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

Then preprocess your VCF files:

```bash
afquery preprocess \
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

# Batch query (100 positions)
afquery query-batch --db my_db --positions positions.tsv --phenotype E11.9

# Region query
afquery query --db my_db --chrom chr1 --start 1000 --end 10000 --phenotype E11.9 --sex M
```

### 3. Annotate VCF Files

```bash
afquery annotate \
  --db my_db \
  --vcf input.vcf \
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
result = db.query(
    chrom="chr1",
    pos=1000,
    alt="G",
    phenotype_codes=["E11.9"],
    sex="both"
)
print(f"AC={result.ac}, AN={result.an}, AF={result.af}")

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
    phenotype_codes=["E11.9", "I10"]
)

# Annotate VCF with allele frequencies
# Note: tech_name filters annotation to samples of that technology
db.annotate_vcf(
    vcf_path="input.vcf",
    output_path="annotated.vcf",
    phenotype_codes=["E11.9"],
    tech_name="wgs"  # Only annotate using WGS samples
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

### Incremental Updates (add_samples)

```bash
afquery add-samples \
  --db my_db \
  --manifest new_samples.tsv \
  --vcf-dir ./new_vcfs/
```

Adds new samples without rebuilding the entire database.

### Compact Database

Remove samples and reclaim disk space:

```bash
afquery compact --db my_db --samples-to-remove sample_1,sample_2
```

### Run Benchmarks

```bash
afquery benchmark --db my_db --n-queries 1000 --query-type point
```

### Generate Synthetic Data

```bash
afquery synth --output synthetic_db/ --n-samples 5000 --n-variants 100000
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
afquery query                Query single position
afquery query-batch          Batch query multiple positions
afquery annotate             Annotate VCF file
afquery info                 Show database info
afquery preprocess           Build database from VCFs
afquery add-samples          Add new samples to database
afquery compact              Remove samples and reclaim space
afquery synth                Generate synthetic test database
afquery benchmark            Run performance benchmarks
```

Run `afquery --help` for full options.

## Development

### Running Tests

```bash
# All 190 tests
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

### Architecture

See `brain/architecture.md` for detailed system design, data flow, and query algorithm.

## Genome Builds

- **GRCh37** (hg19) — PAR regions: chrX:1-2649520
- **GRCh38** (hg38) — PAR regions: chrX:1-3099677

## Technologies Supported

- **WGS** — Whole genome sequencing (always fully covered, no BED file needed)
  - Manifest: `tech_name = wgs` (case-insensitive)
  - Query-time: All positions in genome considered covered

- **WES** — Whole exome sequencing (coverage defined by BED file)
  - Manifest: `tech_name = wes_kit_a` (or any custom name)
  - Preprocessing: Loads `{bed_dir}/wes_kit_a.bed` (0-based half-open BED format)
  - Query-time: Only positions within BED intervals considered covered

- **Custom** — Any technology with a BED file (e.g., gene panels, targeted sequencing)
  - Manifest: Use any `tech_name`
  - Preprocessing: Loads `{bed_dir}/{tech_name}.bed`
  - Query-time: Respects BED file coverage

## Troubleshooting

### ImportError: `cyvcf2`

`cyvcf2` import happens inside worker processes during preprocessing. Do not import at module level.

### DuckDB Temp Files

Uses Parquet format (not Arrow IPC) for compatibility. Set `DUCKDB_TEMP_DIRECTORY` if needed.

### Sample IDs

Sample IDs are 0-indexed and monotonically increasing. Never reuse removed IDs—use `compact` to reclaim space.

## Contributing

1. Read `brain/project_state.json` for current phase and test count
2. Read `brain/architecture.md` for system design
3. Follow code conventions in `CLAUDE.md`
4. Update `brain/` docs after architectural changes
5. Run tests before submitting

## License

(Add your license here)

## Citation

If you use afquery in research, please cite:

```
(Citation format to be determined)
```

---

**Status**: Phase 5 complete (190 tests passing). Active development.
