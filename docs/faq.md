# FAQ

## Does AFQuery need a server?

No. AFQuery is entirely file-based. The database is a directory of Parquet files and a SQLite metadata file. Queries run in-process using DuckDB and pyroaring. No server, no daemon, no cloud service is required.

---

## Can I query multiple chromosomes at once?

Not with a single `afquery query` call — each call is restricted to one chromosome. For multi-chromosome queries, use a shell loop:

```bash
for chrom in chr1 chr2 chr3; do
  afquery query --db ./db/ --chrom $chrom --pos 925952 --format tsv
done
```

Or use the Python API:

```python
from afquery import Database

db = Database("./db/")
for chrom in ["chr1", "chr2", "chr3"]:
    results = db.query(chrom, pos=925952)
    for r in results:
        print(chrom, r.AC, r.AN, r.AF)
```

---

## What VCF format is supported?

Single-sample VCF files, either uncompressed (`.vcf`) or bgzip-compressed (`.vcf.gz`). The VCF must contain genotype (GT) information.

Multi-sample VCFs are not supported. Split them first with bcftools:

```bash
bcftools +split multi_sample.vcf.gz -Oz -o split_vcfs/
```

---

## Can I use multi-sample VCFs?

Not directly. AFQuery expects one VCF per sample. Split multi-sample VCFs with bcftools before building the database:

```bash
bcftools +split multi_sample.vcf.gz -Oz -o ./split_vcfs/
```

Then create a manifest pointing to the split files.

---

## How do I update phenotype metadata without re-ingesting?

This is not yet supported. Phenotype codes are stored in `metadata.sqlite` and are associated with samples at ingest time. To change phenotype assignments, you would need to rebuild the database.

This feature (in-place metadata update) is planned for a future release.

---

## What is the maximum cohort size?

AFQuery has been tested with up to 50,000 samples. Bitmap operations remain fast at this scale because Roaring Bitmaps are highly compressed for sparse data. At 50K samples, a typical variant bitmap is ~64 KB.

For cohorts larger than 50K, performance should remain sub-second for point queries, but build-phase memory requirements scale with cohort size. See [Performance Tuning](advanced/performance.md).

---

## Is GRCh37/hg19 supported?

Yes. Both GRCh37 and GRCh38 are supported. Specify at database creation:

```bash
afquery create-db --genome-build GRCh37 ...
```

The genome build affects PAR1/PAR2 coordinates on chrX (see [Ploidy & Sex Chromosomes](advanced/ploidy-and-sex-chroms.md)). Chromosome names in your VCFs should match the chosen build (`chr1`/`1` both work — `normalize_chrom()` handles the `chr` prefix).

---

## Can I share the database?

Yes. The database is just files — copy or rsync the entire directory:

```bash
rsync -av ./db/ user@remote:/path/to/db/
```

The database is read-only after creation (queries do not write). Multiple processes can query the same database simultaneously without conflict.

---

## How do I see what phenotype codes are in my database?

```bash
afquery info --db ./db/
```

Or via Python:

```python
from afquery import Database
db = Database("./db/")
print(db.get_all_phenotypes())
```

---

## Does filtering by technology affect WGS samples?

WGS samples are always covered at every position (no BED file). Technology filters work by sample membership, not coverage:
- `--tech wgs` restricts to samples with `tech=wgs` in the manifest
- WES samples at positions outside their capture BED are excluded by coverage, not by the tech filter

---

## Why does AF differ between sex-filtered and unfiltered queries on chrX?

On chrX non-PAR positions, males contribute AN=1 and females contribute AN=2. When you filter by `--sex female`, AN increases (purely diploid denominator) and AF may change. This is correct ploidy-aware behavior. See [Ploidy & Sex Chromosomes](advanced/ploidy-and-sex-chroms.md).

---

## What does FAIL_SAMPLES mean in query output?

`FAIL_SAMPLES` (shown as `FAIL=N` in text output) is the count of eligible samples that had the alt allele called but with `FILTER≠PASS`. These samples are not counted in AC/AN. Available only in schema v2 databases. See [FILTER=PASS Tracking](advanced/filter-pass-tracking.md).
