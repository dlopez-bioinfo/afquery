# Annotate a VCF

`afquery annotate` adds allele frequency information to an existing VCF file as INFO fields. It queries each variant in the VCF against the database and writes results inline.

---

## Basic Usage


```bash
afquery annotate \
  --db ./db/ \
  --input variants.vcf.gz \
  --output annotated.vcf.gz
```

---

## Added INFO Fields

| Field | Type | Number | Description |
|-------|------|--------|-------------|
| `AFQUERY_AC` | Integer | A (per ALT) | Allele count in eligible samples |
| `AFQUERY_AN` | Integer | 1 (per site) | Allele number (total alleles examined) |
| `AFQUERY_AF` | Float | A (per ALT) | Allele frequency (`AC / AN`) |
| `AFQUERY_N_HET` | Integer | A (per ALT) | Heterozygous sample count |
| `AFQUERY_N_HOM_ALT` | Integer | A (per ALT) | Homozygous alt sample count |
| `AFQUERY_N_HOM_REF` | Integer | A (per ALT) | Homozygous ref sample count |
| `AFQUERY_N_FAIL` | Integer | 1 (per site) | Samples with FILTER≠PASS and alt allele called |

!!! note "Multi-allelic sites"
    Number=A fields have one value per ALT allele (comma-separated for multi-allelic sites). Number=1 fields are shared across all ALT alleles at the same position.

---

## Sample Filtering

Annotate with a specific subgroup:

```bash
afquery annotate \
  --db ./db/ \
  --input variants.vcf \
  --output annotated.vcf \
  --phenotype E11.9 \
  --sex female \
  --tech wgs
```

The INFO fields will reflect AF computed only over the filtered sample set. This allows generating population-specific frequency tracks.

---

## Parallelism

Annotation runs in parallel across variants. By default, all available CPU cores are used:

```bash
afquery annotate \
  --db ./db/ \
  --input variants.vcf \
  --output annotated.vcf \
  --threads 8
```

For large VCFs (100K+ variants), set `--threads` to the number of available cores.

---


## Using Annotated Output

### BCFtools

Filter variants with high AF:

```bash
bcftools filter -i 'AFQUERY_AF > 0.01' annotated.vcf
```

Extract specific fields:

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AFQUERY_AC\t%AFQUERY_AN\t%AFQUERY_AF\n' annotated.vcf
```

### Python (pysam / cyvcf2)

```python
import cyvcf2

vcf = cyvcf2.VCF("annotated.vcf")
for variant in vcf:
    ac = variant.INFO.get("AFQUERY_AC")
    an = variant.INFO.get("AFQUERY_AN")
    af = variant.INFO.get("AFQUERY_AF")
    print(f"{variant.CHROM}:{variant.POS} AC={ac} AN={an} AF={af}")
```

### R (VariantAnnotation)

```r
library(VariantAnnotation)
vcf <- readVcf("annotated.vcf")
info(vcf)$AFQUERY_AF
```

---

## Full Option Reference

See [CLI Reference → annotate](../reference/cli.md#annotate).
