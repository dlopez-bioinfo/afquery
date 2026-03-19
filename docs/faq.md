# FAQ

## Does AFQuery need a server?

No. AFQuery is entirely file-based. The database is a directory of Parquet files and a SQLite metadata file. Queries run in-process using DuckDB and pyroaring. No server, no daemon, no cloud service is required.

---

## Can I query multiple chromosomes at once?

Yes. The `--from-file` batch mode supports variants across multiple chromosomes in a single call. Place variants from any chromosome in the same TSV file:

```tsv
chr1	925952	G	A
chrX	5000000	A	G
chr2	166845670	C	T
```

```bash
afquery query --db ./db/ --from-file variants.tsv --format tsv
```

Results are returned in the same order as the input file. You can also use the Python API directly:

```python
from afquery import Database

db = Database("./db/")
results = db.query_batch_multi([
    ("chr1", 925952, "G", "A"),
    ("chrX", 5000000, "A", "G"),
])
for r in results:
    print(r.variant.chrom, r.variant.pos, r.AC, r.AN, r.AF)
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

Use `afquery update-db --update-sample` to correct a sample's `sex` or `phenotype_codes` without touching the VCF or the Parquet files. The change is recorded in the changelog and precomputed bitmaps are regenerated automatically.

Update a single sample:

```bash
afquery update-db --db ./db/ --update-sample SAMP_001 --set-phenotype "E11.9,I10"
afquery update-db --db ./db/ --update-sample SAMP_001 --set-sex female
```

Update many samples at once using a TSV file (header: `sample_name`, `field`, `new_value`):

```bash
afquery update-db --db ./db/ --update-samples-file corrections.tsv
```

See [Updating sample metadata](guides/update-database.md#update-sample-metadata) for a full walkthrough.

!!! note "Code format"
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

The genome build affects PAR1/PAR2 coordinates on chrX (see [Ploidy & Special Chromosomes](advanced/ploidy-and-sex-chroms.md)). Chromosome names in your VCFs should match the chosen build (`chr1`/`1` both work — `normalize_chrom()` handles the `chr` prefix).

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

On chrX non-PAR positions, males contribute AN=1 and females contribute AN=2. When you filter by `--sex female`, AN increases (purely diploid denominator) and AF may change. This is correct ploidy-aware behavior. See [Ploidy & Special Chromosomes](advanced/ploidy-and-sex-chroms.md).

---

## What does FAIL_SAMPLES mean in query output?

`FAIL_SAMPLES` (shown as `FAIL=N` in text output) is the count of eligible samples that had the alt allele called but with `FILTER≠PASS`. These samples are not counted in AC/AN. See [FILTER=PASS Tracking](advanced/filter-pass-tracking.md).

---

## Can I use any strings as phenotype codes?

Yes. The `phenotype_codes` column in the manifest accepts arbitrary comma-separated strings. Common choices include ICD-10 codes (`E11.9`), HPO terms (`HP:0001250`), or project-specific tags (`control`, `pilot`). There is no controlled vocabulary — you define the coding scheme for your cohort.

See [Manifest Format](guides/manifest-format.md#phenotype-codes) for details.

---

## What happens if I query a phenotype code that does not exist?

AFQuery emits an `AfqueryWarning` to stderr and returns `AC=0, AN=0, AF=None`. The warning message indicates whether the unknown code was used as an include filter (will match 0 samples) or an exclude filter (has no effect):

```
AfqueryWarning: Phenotype 'CODE' not in database — include will match 0 samples.
AfqueryWarning: Phenotype 'CODE' not in database — exclude has no effect.
```

Use `afquery info --db ./db/` to list all registered codes before running queries. Use `--no-warn` to suppress warnings in batch pipelines.

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

## What are AFQuery's limitations?

AFQuery is purpose-built for fast subcohort AF computation and is not a general-purpose genomic database:

- **Not a joint genotyper**: AFQuery does not perform joint genotyping. Input VCFs should be individually called before ingestion.
- **Not a variant database**: AFQuery stores only genotype-level summaries (bitmaps). Individual sample genotypes cannot be retrieved from the database.
- **No statistical genetics**: AFQuery does not compute Hardy-Weinberg equilibrium, population stratification, or other statistical genetics metrics.
- **Batch queries**: The `--from-file` batch mode supports variants across multiple chromosomes in a single call. Point queries (`--locus`) and region queries (`--region`) target a single position or range; for multi-position multi-chromosome lookups, use `--from-file`.
- **Cohort size limit**: Performance at >100K samples has not been validated. Memory requirements for the build phase scale with cohort size.

---

### Why does `create-db` abort with "BED file not found"?

AFQuery requires a BED capture file for every non-WGS technology. If the file is missing,
the run aborts immediately — before any VCF ingestion begins — with:

```
ManifestError: BED file not found for technology '<tech>': '<path>'
```

**How to fix:**

- Verify the BED file path exists and is readable.
- Ensure the file is named exactly `<tech_name>.bed` and is located in the directory
  passed to `--bed-dir`.
- If the samples are genuinely whole-genome (no capture), change `tech_name` to `WGS`
  in the manifest (no BED file required for WGS).
