# Understanding Output

This page explains what each field in AFQuery output means and how to interpret special cases.

---

## Output Fields

| Field | Type | Description |
|-------|------|-------------|
| **AC** | int | Allele count — number of alt allele copies in eligible samples |
| **AN** | int | Allele number — total alleles examined (2 × diploid samples, adjusted for ploidy) |
| **AF** | float or None | Allele frequency — `AC / AN`. `None` when AN=0 |
| **N_HET** | int | Number of eligible samples heterozygous for the alt allele (GT=0/1) |
| **N_HOM_ALT** | int | Number of eligible samples homozygous for the alt allele (GT=1/1) |
| **N_HOM_REF** | int | Number of eligible samples homozygous reference (GT=0/0) at this position |
| **FAIL_SAMPLES** | int or None | Number of eligible samples with FILTER≠PASS at this position. `None` for v1 databases |

### Relationships

```
AC = N_HET + (2 × N_HOM_ALT)           # on autosomes
AN = 2 × (N_HET + N_HOM_ALT + N_HOM_REF)  # on autosomes
AF = AC / AN                             # None when AN=0
```

On sex chromosomes, these relationships adjust for ploidy. See [Ploidy & Sex Chromosomes](../advanced/ploidy-and-sex-chroms.md).

---

## Output Formats

### Text (default)

```bash
afquery query --db ./db/ --chrom chr1 --pos 925952 --ref G --alt A
```

```
chr1:925952 G>A  AC=3  AN=120  AF=0.0250  N_HET=1  N_HOM_ALT=1  N_HOM_REF=57
```

### TSV

```bash
afquery query --db ./db/ --chrom chr1 --pos 925952 --ref G --alt A --format tsv
```

```
chrom	pos	ref	alt	AC	AN	AF	N_HET	N_HOM_ALT	N_HOM_REF	FAIL_SAMPLES
chr1	925952	G	A	3	120	0.0250	1	1	57	0
```

### JSON

```bash
afquery query --db ./db/ --chrom chr1 --pos 925952 --ref G --alt A --format json
```

```json
{
  "chrom": "chr1",
  "pos": 925952,
  "ref": "G",
  "alt": "A",
  "AC": 3,
  "AN": 120,
  "AF": 0.025,
  "N_HET": 1,
  "N_HOM_ALT": 1,
  "N_HOM_REF": 57,
  "FAIL_SAMPLES": 0
}
```

---

## Special Cases

### AN=0 and AF=None

AN=0 means no eligible samples have coverage at this position. This happens when:

- All eligible samples are WES and the position is outside capture regions
- The phenotype/sex/technology filter excludes all samples
- The chromosome is not in the database

When AN=0, AF cannot be computed (division by zero) and is reported as `None` (text/JSON) or `.` (TSV/VCF).

!!! warning "AN=0 does not mean the variant is absent"
    AN=0 means AFQuery has no data to compute frequency. It is not evidence of rarity. Do not use AN=0 results for PM2 classification.

### AC=0 with High AN

AC=0 with a high AN (e.g., AN=4000) means the variant was **genuinely not observed** in a well-covered cohort. This is meaningful evidence that the variant is rare or absent in your population.

### FAIL_SAMPLES=None

A `None` value for FAIL_SAMPLES indicates the database was built with schema v1 (before FILTER=PASS tracking was added). To get FAIL_SAMPLES data, rebuild the database with the current AFQuery version. See [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md).

### FAIL_SAMPLES > 0

When FAIL_SAMPLES > 0, some eligible samples had a non-PASS filter at this position in their source VCF. A high FAIL_SAMPLES relative to AN may indicate a problematic site (e.g., systematic sequencing artifacts). Consider filtering results where `FAIL_SAMPLES / (AN / 2) > 0.1` (more than 10% of samples fail).

---

## VCF Annotation Fields

When using `afquery annotate`, the following INFO fields are added to each variant:

| INFO field | Description |
|-----------|-------------|
| `AFQUERY_AC` | Allele count |
| `AFQUERY_AN` | Allele number |
| `AFQUERY_AF` | Allele frequency |
| `AFQUERY_N_HET` | Heterozygous sample count |
| `AFQUERY_N_HOM_ALT` | Homozygous alt sample count |
| `AFQUERY_N_HOM_REF` | Homozygous ref sample count |
| `AFQUERY_FAIL` | Fail sample count (schema v2 only) |

These fields can be used directly in downstream filtering with `bcftools filter`:

```bash
# Keep variants rare in cohort with sufficient coverage
bcftools filter -i 'AFQUERY_AF < 0.001 && AFQUERY_AN >= 1000' annotated.vcf.gz
```

---

## Related Pages

- [Key Concepts](concepts.md) — how AC, AN, and AF are computed
- [Query Allele Frequencies](../guides/query.md) — full query CLI reference
- [Debugging Results](../advanced/debugging-results.md) — diagnosing unexpected output
- [Glossary](../reference/glossary.md) — definitions of all terms
