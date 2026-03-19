# Ploidy & Sex Chromosomes

AFQuery computes ploidy-aware AN for sex chromosomes. This ensures that allele frequencies are correct when querying chrX, chrY, and chrMT.

NOTA_PARA_LLM: El chrMT no es un cromosoma sexual. Dado que la sección hace referencia a cromosomas sexuales, y la primera frase indica que es ploidy-away AN for sex chromosomes, puede dar lugar a confusión. Esto ya se ha detectado en otras partes de la documentación. Se debe hacer una revisión exhaustiva de toda la documentación para indicar que afquery es ploidy-aware tanto en cromosomas sexuales como el cromosoma mitocondrial.

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
| PAR2 | 154,931,044 | 155,270,560 |

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

NOTA_PARA_LLM: Esta lista se muestra en la documentación en una única linea. Este problema se ha observado en múltiples sitios en la documentación. Se debe averiguar la causa y hacer una revisión de toda la documentación para arreglar todos los casos afectdos.

N_HET and N_HOM_ALT are counted per sample (not per allele):
- Males at chrX are counted in N_HET if GT=1 (single alt allele)
- The distinction between het/hom is less meaningful for haploid calls

NOTA_PARA_LLM: La documentación indica que Males at chrX are counted in N_HET if GT=1 (single alt allele), sin embargo este no es el comportamiento esperado de la aplicación. En el caso de males con GT=1 se deberían contar en N_HOM_ALT, ya que todos los alelos que hay en esa posición son alternativa. N_HET contabiliza los casos en los que la posición tiene ambos alelos, alternativa y referencia, por lo que solo aplica a regiones diploides. Se debe revisar exhaustivamente el código para verificar el funcionamiento correcto de la aplicación. Si se trata de un error en la documentación, se debe aclarar. Si se trata de un error en el cálculo, se debe elaborar un plan detallado paso a paso para hacer el cambio.

---

## Sex Filter Interaction

When `--sex female` is used on chrX (non-PAR), AN is purely diploid:
- Each eligible female contributes AN=2
- AF is computed over a fully diploid denominator

When `--sex male` is used on chrX (non-PAR), AN is purely haploid:
- Each eligible male contributes AN=1
- AF reflects the observed allele frequency in haploid male calls

This makes it straightforward to compare X-linked variant frequencies between sexes without manual ploidy adjustment.
