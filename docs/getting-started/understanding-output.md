# Understanding Output

This page explains what each field in AFQuery output means and how to interpret special cases.

---

## Output Fields

| Field | Type | Description |
|-------|------|-------------|
| **AC** | int | Allele count — number of alt allele copies in eligible samples |
| **AN** | int | Allele number — total alleles examined (adjusted for ploidy and eligible samples) |
| **AF** | float | Allele frequency — `AC / AN`. `None` when AN=0 |
| **N_HET** | int | Number of eligible samples heterozygous for the alt allele (GT=0/1) |
| **N_HOM_ALT** | int | Number of eligible samples homozygous for the alt allele (GT=1/1 or GT=1) |
| **N_HOM_REF** | int | Number of eligible samples homozygous reference (GT=0/0 or GT=0) |
| **FAIL_SAMPLES** | int | Number of eligible samples with FILTER≠PASS at this position. |


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

- All eligible samples are WES or panels and the position is outside capture regions
- The phenotype/sex/technology filter excludes all samples
- The chromosome is not in the database


!!! warning "AN=0 does not mean the variant is absent"
    AN=0 means AFQuery has no data to compute frequency. It is not evidence of rarity.

### Warnings

afquery emits a `AfqueryWarning` to stderr when a query may silently return fewer or no results. Common causes:

| Situation | Warning message |
|-----------|----------------|
| Chromosome not in database | `Chromosome 'chrXX' has no data in this database. Available: [...]` |
| Unknown phenotype code (include) | `Phenotype 'CODE' not in database — include will match 0 samples.` |
| Unknown phenotype code (exclude) | `Phenotype 'CODE' not in database — exclude has no effect.` |
| Unknown technology name (include) | `Technology 'NAME' not in database — include will match 0 samples.` |
| Unknown technology name (exclude) | `Technology 'NAME' not in database — exclude has no effect.` |
| Contradictory filters (e.g. include + exclude same code) | `Sample filter produces an empty eligible set — all queries will return AN=0.` |

Use `--no-warn` to suppress these warnings:

```bash
afquery query --db ./my_db/ --chrom chr22 --pos 1000 --no-warn
afquery annotate --db ./my_db/ --input in.vcf --output out.vcf --no-warn
```

### AC=0 with High AN

AC=0 with a high AN (e.g., AN=4000) means the variant was **genuinely not observed** in a well-covered cohort. This is meaningful evidence that the variant is rare or absent in your population.


### FAIL_SAMPLES > 0

When FAIL_SAMPLES > 0, some eligible samples had a non-PASS filter at this position in their source VCF. A high FAIL_SAMPLES relative to AN may indicate a problematic site (e.g., systematic sequencing artifacts). Consider filtering positions with more than 10% of failing (`FAIL_SAMPLES / (AN / 2) > 0.1`).

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

## Next steps

- [Key Concepts](concepts.md) — how AC, AN, and AF are computed
- [Bioinformatics user guide](../guides/create-database.md) — 
- [Clinical user guide](../guides/create-database.md) — 

NOTA_PARA_LLM: Añadir descripción a los enlaces next steps que no tienen descripción.
