# Clinical Variant Prioritization

## Scenario

A rare disease patient has undergone whole-exome sequencing. After standard filtering, 500,000 candidate variants remain. You want to annotate each variant with cohort-specific allele frequency and filter to those that are genuinely rare in your population — removing variants that are common locally but might appear rare in gnomAD.

## Why Standard Databases Fall Short

gnomAD provides an excellent first filter, but:

1. **Population mismatch**: A variant at AF=0.001 in gnomAD may be at AF=0.02 in your local cohort — common locally but appearing rare globally
2. **Coverage gaps**: gnomAD WES has incomplete coverage of some gene panels; your WGS database has full coverage
3. **Local artifacts**: Systematic sequencing artifacts in your pipeline appear as rare variants in gnomAD but are common in your cohort

AFQuery lets you apply cohort-specific AF as an additional filter layer on top of gnomAD, removing locally common variants that standard databases miss.

## Step-by-Step Example

### 1. Annotate patient VCF with cohort AF

```bash
afquery annotate \
  --db ./db/ \
  --input patient.vcf.gz \
  --output patient_annotated.vcf.gz \
  --threads 16
```

This adds to each variant:
- `AFQUERY_AC`: allele count in cohort
- `AFQUERY_AN`: allele number (eligible samples at this position)
- `AFQUERY_AF`: allele frequency
- `AFQUERY_N_HET`, `AFQUERY_N_HOM_ALT`: genotype counts

### 2. Filter for rare variants with reliable AN

```bash
bcftools filter \
  -i 'AFQUERY_AF < 0.001 && AFQUERY_AN >= 1000' \
  patient_annotated.vcf.gz \
  -o patient_rare.vcf.gz
```

The `AFQUERY_AN >= 1000` threshold ensures the AF estimate is based on sufficient data. An AF estimate from AN=10 is meaningless — with AN=10, even a single carrier gives AF=0.1.

### 3. Handle variants not in the cohort

Variants absent from the database have `AFQUERY_AN=0` (no eligible samples, or variant not observed). These require separate treatment:

```bash
# Variants NOT in cohort (AN=0): treat as novel
bcftools filter -i 'AFQUERY_AN == 0' patient_annotated.vcf.gz

# Variants in cohort with sufficient coverage
bcftools filter -i 'AFQUERY_AN >= 1000 && AFQUERY_AF < 0.001' patient_annotated.vcf.gz
```

### 4. Annotation with subgroup AF (optional)

If you want AF relative to a matched control group:

```bash
afquery annotate \
  --db ./db/ \
  --input patient.vcf.gz \
  --output patient_control_af.vcf.gz \
  --phenotype ^rare_disease     # AF in non-rare-disease samples
```

### 5. Full prioritization workflow

```bash
# Step 1: Annotate with cohort AF
afquery annotate --db ./db/ --input patient.vcf.gz --output step1.vcf.gz --threads 16

# Step 2: Keep variants that are rare OR have no cohort coverage
bcftools filter -i '(AFQUERY_AN >= 1000 && AFQUERY_AF < 0.001) || AFQUERY_AN == 0' \
  step1.vcf.gz -o step2.vcf.gz

# Step 3: Flag locally common variants that gnomAD marks as rare
bcftools filter -i 'AFQUERY_AN >= 1000 && AFQUERY_AF >= 0.01' \
  step1.vcf.gz -o locally_common.vcf.gz
```

### 6. Python workflow

```python
import cyvcf2
from afquery import Database

db = Database("./db/")

# Annotate and filter in memory
vcf = cyvcf2.VCF("patient.vcf.gz")
rare_candidates = []

for variant in vcf:
    results = db.query(variant.CHROM, pos=variant.POS, alt=variant.ALT[0])
    if not results:
        rare_candidates.append(variant)  # Not in cohort → novel
        continue
    r = results[0]
    if r.AN >= 1000 and r.AF < 0.001:
        rare_candidates.append(variant)

print(f"Rare candidates: {len(rare_candidates)}")
```

## Biological Interpretation

| Filter | Retained | Removed | Reason |
|--------|----------|---------|--------|
| AFQUERY_AN >= 1000 | 45,000 | 455,000 | Insufficient cohort coverage |
| AFQUERY_AF < 0.001 | 1,200 | 43,800 | Locally common variants |
| Novel (AN=0) | 300 | — | Not observed in cohort |

A typical clinical pipeline retains ~1,500 rare/novel candidates after cohort AF filtering, compared to 500,000 before — a 300× reduction.

**AN threshold guidance:**
- `AN >= 100`: minimum for any AF interpretation
- `AN >= 500`: recommended for rare variant filtering
- `AN >= 1000`: conservative threshold for robust AF estimates

### Mapping to ACMG Criteria

The prioritization workflow above directly supports ACMG/AMP variant classification:

- **BA1 (Stand-alone benign)**: Variants with `AFQUERY_AF > 0.05` and `AFQUERY_AN >= 1000` meet the BA1 threshold for benign classification in your local cohort.
- **PM2 (Supporting pathogenic)**: Variants with `AFQUERY_AC == 0` and `AFQUERY_AN >= 2000` are genuinely absent in your cohort, supporting PM2 application.
- **PS4 (Strong pathogenic)**: Compare AF between affected and unaffected subgroups using `--phenotype` filters to identify enrichment in cases.

For detailed ACMG workflows with worked examples and AN threshold guidance, see [ACMG Criteria (BA1/PM2/PS4)](../clinical/acmg-use-cases.md).

## Related Features

- [Annotate a VCF](../guides/annotate-vcf.md) — full annotation options
- [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md) — using N_FAIL in filtering
- [Population-Specific AF](population-specific-af.md) — local vs. gnomAD comparison
