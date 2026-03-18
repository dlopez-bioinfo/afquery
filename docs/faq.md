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

**This is not supported.** Phenotype codes are assigned at ingest time and stored in `metadata.sqlite`. They cannot be changed without removing and re-adding the sample.

To update a sample's phenotype codes:
1. Remove the sample: `afquery update-db --db ./db/ --remove-samples SAMP_001`
2. Update your manifest with the new codes
3. Re-add the sample: `afquery update-db --db ./db/ --add-samples updated_manifest.tsv`

**Plan your phenotype code assignments carefully before building the database.**

!!! warning "Phenotype codes are permanent"
    Codes are case-sensitive and matched exactly. `E11.9` ≠ `e11.9`. Typos are silently stored and will produce empty results when queried. Run `afquery info --db ./db/` to see all registered codes and verify they match your expectations.

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

---

## Can I use any strings as phenotype codes?

Yes. The `phenotype_codes` column in the manifest accepts arbitrary comma-separated strings. Common choices include ICD-10 codes (`E11.9`), HPO terms (`HP:0001250`), or project-specific tags (`control`, `pilot`). There is no controlled vocabulary — you define the coding scheme for your cohort.

See [Manifest Format](guides/manifest-format.md#phenotype-codes) for details.

---

## What happens if I query a phenotype code that does not exist?

AFQuery returns `AC=0, AN=0, AF=None` — the same result as zero eligible samples. There is no error or warning for unknown codes. Use `afquery info --db ./db/` to list all registered codes before running queries.

---

## Common Pitfalls

### What if AN is very low?

Low AN means the allele frequency estimate is unreliable. For example: AC=1, AN=4, AF=0.25 — this is not a robust 25% frequency estimate.

**Rules of thumb:**
- **AN ≥ 100**: bare minimum for any interpretation
- **AN ≥ 500**: necessary for rare variant filtering
- **AN ≥ 1000**: required for clinical variant interpretation

Always check `AFQUERY_AN` alongside `AFQUERY_AF` in downstream analyses. A variant with high AF but very low AN should be treated with skepticism.

---

### What about population stratification in mixed-ancestry cohorts?

If your cohort is a mix of ancestries, the AF reflects a weighted average across all ancestries. A variant at AF=0.001 globally may be at AF=0.01 in a subpopulation.

**Mitigation:**
- Tag samples by ancestry or population in `phenotype_codes`
- Use stratified queries: `afquery query --phenotype AFR` for African samples, `--phenotype EUR` for European samples, etc.
- Compare AF across subgroups to detect ancestry-specific signals

Without stratification, rare variant interpretation in mixed cohorts can be misleading.

---

### Why does my WES database have inflated AN at off-target positions?

AFQuery silently treats WES as WGS if the BED file is absent during `create-db`. The database will issue a warning but does not error. Result: off-target positions count towards AN as if they were sequenced.

**Always verify after build:**
```bash
afquery check --db ./db/
```

The check command validates that all WES samples have proper BED file coverage and reports any missing BED assignments.

---

### Can I compare AF across schema versions?

No, not directly. v1 databases do not track `FILTER≠PASS` failures; `AFQUERY_N_FAIL` is always `0` on v1. If your source VCFs had common non-PASS genotypes (e.g., high sequencing error rate or low-confidence calls), AC computed from v1 databases may differ from v2.

**Solution:** Rebuild v1 databases with the current version of AFQuery for consistent FAIL tracking across your cohort.
