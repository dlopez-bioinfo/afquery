# Sex-Specific Allele Frequency

## Scenario

You are analyzing a candidate X-linked recessive gene. A hemizygous variant in males may be pathogenic, but the same variant in females acts as a heterozygous carrier state. Allele frequency means something different for each sex due to ploidy, and a combined AF obscures the clinically relevant frequencies.

## Why Standard Databases Fall Short

Most population databases report a single AF for chrX variants, combining diploid female calls and haploid male calls. This average is appropriate for some purposes but misleading for X-linked recessive analysis, where the male hemizygous frequency is the clinically relevant metric.

For a variant at chrX non-PAR with AC=10 in 500 males (AN=500) and AC=5 in 1000 females (AN=2000):
- Male hemizygous rate: 10/500 = 2%
- Female carrier rate: 5/2000 = 0.25%
- Combined AF: 15/2500 = 0.6%

The combined AF (0.6%) obscures the fact that 2% of males are hemizygous carriers — a crucial number for X-linked recessive risk assessment.

## Solution with AFQuery

AFQuery applies ploidy-aware AN computation automatically. For chrX non-PAR positions:
- Males contribute AN=1 (haploid)
- Females contribute AN=2 (diploid)

Using `--sex male` returns AF over a purely haploid denominator — the true hemizygous frequency. Using `--sex female` returns AF over the diploid denominator.

## Computing Sex-Specific AF at chrX

### 1. Query full cohort (combined)

```bash
afquery query --db ./db/ --locus chrX:153296777
```

```
chrX:153296777
  REF=C  ALT=T  AC=15  AN=2500  AF=0.0060  N_HET=5  N_HOM_ALT=0
```

AN=2500: 500 males × 1 + 1000 females × 2 = 2500 (mixed ploidy)

### 2. Query males only (hemizygous frequency)

```bash
afquery query --db ./db/ --locus chrX:153296777 --sex male
```

```
chrX:153296777
  REF=C  ALT=T  AC=10  AN=500  AF=0.0200  N_HET=10  N_HOM_ALT=0
```

Hemizygous rate: 2% (N_HET here represents hemizygous males, GT=1)

### 3. Query females only (carrier frequency)

```bash
afquery query --db ./db/ --locus chrX:153296777 --sex female
```

```
chrX:153296777
  REF=C  ALT=T  AC=5  AN=2000  AF=0.0025  N_HET=5  N_HOM_ALT=0
```

Carrier rate: 0.25%

### 4. Python API comparison

```python
from afquery import Database

db = Database("./db/")

male = db.query("chrX", pos=153296777, sex="male")[0]
female = db.query("chrX", pos=153296777, sex="female")[0]
combined = db.query("chrX", pos=153296777)[0]

print(f"Male hemizygous rate:  {male.AF:.4f}  (AN={male.AN})")
print(f"Female carrier rate:   {female.AF:.4f}  (AN={female.AN})")
print(f"Combined AF:           {combined.AF:.4f}  (AN={combined.AN})")
```

### 5. Annotating with sex-stratified AF

```bash
# Annotate using only male samples (for X-linked recessive analysis)
afquery annotate \
  --db ./db/ \
  --input patient.vcf.gz \
  --output annotated_males.vcf.gz \
  --sex male
```

## Biological Interpretation

| Metric | Value | Clinical relevance |
|--------|-------|-------------------|
| Male hemizygous rate | 2% | ACMG BS1 for X-linked recessive (exceeds 0.01 threshold) |
| Female carrier rate | 0.25% | Heterozygous carriers; does not inform X-linked recessive pathogenicity |
| Combined AF | 0.6% | Misleading for X-linked analysis — neither the hemizygous nor carrier rate |

For X-linked recessive disorders, always use `--sex male` when evaluating the hemizygous carrier frequency.

## PAR Regions

The Pseudoautosomal Regions (PAR1/PAR2) on chrX behave like autosomes. At PAR positions, males contribute AN=2 and sex-stratified queries are less meaningful. AFQuery applies the correct PAR/non-PAR boundary automatically.

See [Ploidy & Special Chromosomes](../advanced/ploidy-and-sex-chroms.md) for PAR coordinates.

## Related Features

- [Ploidy & Special Chromosomes](../advanced/ploidy-and-sex-chroms.md) — PAR regions, ploidy rules
- [Sample Filtering](../guides/sample-filtering.md) — sex filter syntax
- [Clinical Prioritization](clinical-prioritization.md) — VCF annotation with sex filter
