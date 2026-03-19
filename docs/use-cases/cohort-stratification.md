# Cohort Stratification

## Scenario

Your cohort includes samples from multiple disease groups and multiple sequencing technologies. You want to compare allele frequencies across groups to (a) identify technology-driven artifacts and (b) test whether a variant is enriched in a specific disease group.

## Why Standard Databases Fall Short

General population databases do not allow stratification by disease group or sequencing technology. They provide a single aggregate frequency, hiding heterogeneity that may be scientifically or clinically significant.

## Solution with AFQuery

AFQuery allows computing AF over any combination of metadata filters on the same database. You can compare:
- WGS vs. WES (technology artifact detection)
- Disease group A vs. disease group B (disease enrichment)
- Within-technology disease comparisons (controlling for technology)

All comparisons happen on the same database.

## Comparing WGS vs. WES and Disease Groups

### Setup: Multi-technology, multi-disease cohort

```tsv
sample_name	vcf_path	sex	tech_name	phenotype_codes
SAMP_001	vcfs/SAMP_001.vcf.gz	female	wgs	epilepsy
SAMP_002	vcfs/SAMP_002.vcf.gz	male	wgs	control
SAMP_003	vcfs/SAMP_003.vcf.gz	female	wes_v1	epilepsy
SAMP_004	vcfs/SAMP_004.vcf.gz	male	wes_v1	control
...
```

### 1. Technology comparison (artifact detection)

```bash
# WGS samples
afquery query --db ./db/ --locus chr1:925952 --tech wgs --format tsv

# WES samples
afquery query --db ./db/ --locus chr1:925952 --tech wes_v1 --format tsv
```

If AF differs substantially between technologies (>2×), this may indicate a capture artifact or systematic genotyping difference. Use `afquery check` to verify BED file coverage.

### 2. Disease group comparison

```bash
# Epilepsy group
afquery query --db ./db/ --locus chr1:925952 --phenotype epilepsy --format json

# Control group (explicitly tagged)
afquery query --db ./db/ --locus chr1:925952 --phenotype control --format json

# Or use exclusion-based pseudo-controls
afquery query --db ./db/ --locus chr1:925952 --phenotype ^epilepsy --format json
```

### 3. Within-technology disease comparison

```bash
# Epilepsy WGS only (controls for technology)
afquery query --db ./db/ --locus chr1:925952 \
  --phenotype epilepsy --tech wgs --format tsv

# Control WGS only
afquery query --db ./db/ --locus chr1:925952 \
  --phenotype control --tech wgs --format tsv
```

### 4. Systematic stratification with dump

For genome-wide stratified analysis:

```bash
afquery dump \
  --db ./db/ \
  --by-phenotype epilepsy --by-phenotype control \
  --by-tech \
  --output stratified_af.csv
```

### 5. Python: compute enrichment ratio

```python
from afquery import Database

db = Database("./db/")

disease = db.query("chr1", pos=925952, phenotype=["epilepsy"])
controls = db.query("chr1", pos=925952, phenotype=["control"])

if disease and controls and disease[0].AF and controls[0].AF:
    enrichment = disease[0].AF / controls[0].AF
    print(f"Disease AF: {disease[0].AF:.4f}  (AN={disease[0].AN})")
    print(f"Control AF: {controls[0].AF:.4f}  (AN={controls[0].AN})")
    print(f"Enrichment: {enrichment:.1f}×")
```

## Biological Interpretation

| Comparison | Disease AF | Control AF | Enrichment | Interpretation |
|------------|------------|------------|------------|----------------|
| Epilepsy vs. control (WGS) | 0.012 | 0.003 | 4× | Possible disease association |
| WGS vs. WES | 0.011 | 0.010 | 1.1× | No technology artifact |
| Epilepsy WGS vs. WES | 0.012 | 0.013 | 0.9× | Consistent across technologies |

Low AN values (<100) indicate small group sizes — interpret those frequencies cautiously.

## Related Features

- [Bulk Export](../guides/dump-export.md) — export stratified AF across all variants
- [Pseudo-controls](pseudo-controls.md) — exclusion-based background frequency
- [Technology Integration](technology-integration.md) — WES capture regions and AN
