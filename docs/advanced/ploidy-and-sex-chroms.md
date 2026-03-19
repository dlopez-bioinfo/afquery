# Ploidy & Special Chromosomes

AFQuery computes ploidy-aware AN for sex chromosomes (chrX, chrY) and the mitochondrial chromosome (chrMT). This ensures that allele frequencies are correct when querying these chromosomes, where the number of alleles per sample differs from the diploid autosomes.

---

## Ploidy Rules

| Chromosome | Female AN contribution | Male AN contribution |
|------------|----------------------|---------------------|
| Autosomes (chr1–22) | 2 | 2 |
| chrX (non-PAR) | 2 | 1 |
| chrX (PAR1, PAR2) | 2 | 2 |
| chrY | 0 | 1 |
| chrMT | 1 | 1 |

For each eligible sample at a given position, AFQuery adds the appropriate ploidy count to AN based on the sample's sex and the chromosome/position.

---

## Pseudoautosomal Regions (PAR)

The pseudoautosomal regions on chrX behave like autosomes — both males and females contribute AN=2. PAR coordinates by genome build:

### GRCh38

| Region | Start | End |
|--------|-------|-----|
| PAR1 | 60,001 | 2,699,520 |
| PAR2 | 154,931,044 | 155,270,560 |

### GRCh37 / hg19

| Region | Start | End |
|--------|-------|-----|
| PAR1 | 60,001 | 2,699,520 |
| PAR2 | 154,931,044 | 155,260,560 |

Positions within PAR1 or PAR2 on chrX are treated as diploid for all samples.

---

## Effect on AF Queries

### chrY

Querying chrY with `--sex female` returns `AN=0` (females have no Y chromosome):

```bash
afquery query --db ./db/ --locus chrY:2787758 --sex female
# chrY:2787758 — no results (AN=0 for all variants)

afquery query --db ./db/ --locus chrY:2787758 --sex male
# chrY:2787758  REF=C  ALT=T  AC=3  AN=856  AF=0.0035
```

### chrX non-PAR

Male samples contribute AN=1, female samples contribute AN=2. This means a cohort of 500 females and 500 males has AN = 500×2 + 500×1 = 1500 at a non-PAR X position.

```bash
afquery query --db ./db/ --locus chrX:100000000
# N_total = 1000 samples, AN = 1500 (not 2000)
```

### chrMT

All samples are haploid at mitochondrial loci:

```bash
afquery query --db ./db/ --locus chrMT:3243
# AN = n_samples (one allele per sample)
```

---

## Genotype Counting

At non-PAR chrX positions:
- A male with GT=`1` contributes AC=1, AN=1
- A female with GT=`0/1` contributes AC=1, AN=2
- A female with GT=`1/1` contributes AC=2, AN=2

N_HET and N_HOM_ALT are counted per sample (not per allele):
- Males at chrX non-PAR (haploid positions) are counted in **N_HOM_ALT** when GT=1, because all alleles at that position are alternate. N_HET is reserved for diploid positions where both reference and alternate alleles are present.
- Females at chrX with GT=`0/1` are counted in N_HET; with GT=`1/1` in N_HOM_ALT.

---

## Sex Filter Interaction

When `--sex female` is used on chrX (non-PAR), AN is purely diploid:
- Each eligible female contributes AN=2
- AF is computed over a fully diploid denominator

When `--sex male` is used on chrX (non-PAR), AN is purely haploid:
- Each eligible male contributes AN=1
- AF reflects the observed allele frequency in haploid male calls

This makes it straightforward to compare X-linked variant frequencies between sexes without manual ploidy adjustment.
