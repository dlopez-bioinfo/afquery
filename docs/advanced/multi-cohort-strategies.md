# Multi-Cohort Strategies

When your organization manages samples from multiple cohorts — different institutions, studies, or disease programs — you need a strategy for how to organize AFQuery databases. This page covers three common patterns with trade-offs.

---

## Pattern 1: One Database per Cohort

Each cohort gets its own database, queried independently.

```
/databases/
  cardiology_cohort/       ← 5000 samples, cardiology
  neurology_cohort/        ← 3000 samples, neurology
  rare_disease_registry/   ← 8000 samples, mixed rare diseases
```

### When to Use

- Cohorts come from different institutions with separate data governance
- Sample sets have no overlap
- Each cohort has its own VCF pipeline and genome build
- You need to annotate VCFs against a specific cohort only

### Trade-offs

| Advantage | Disadvantage |
|-----------|-------------|
| Simple data governance — each database is self-contained | Cannot compute cross-cohort AF in a single query |
| Independent updates — rebuild one without touching others | Duplicate storage if samples overlap |
| Clear provenance — each database tracks its own manifest | More databases to maintain |

### Cross-Cohort Comparison (Python)

```python
from afquery import Database

dbs = {
    "cardiology": Database("/databases/cardiology_cohort/"),
    "neurology": Database("/databases/neurology_cohort/"),
    "rare_disease": Database("/databases/rare_disease_registry/"),
}

chrom, pos, alt = "chr1", 12345678, "T"

for name, db in dbs.items():
    results = db.query(chrom, pos=pos, alt=alt)
    if results:
        r = results[0]
        print(f"{name:20s}  AC={r.AC:4d}  AN={r.AN:5d}  AF={r.AF:.4f}")
    else:
        print(f"{name:20s}  not observed")
```

---

## Pattern 2: Merged Database with Phenotype Codes

All samples in one database, with phenotype codes distinguishing cohorts.

### Manifest Design

```tsv
sample_id	vcf_path	sex	technology	phenotype_codes
CARD_001	/data/card/001.vcf.gz	female	wgs	cardiology,EUR,control
CARD_002	/data/card/002.vcf.gz	male	wgs	cardiology,EUR,case_HCM
NEURO_001	/data/neuro/001.vcf.gz	female	wes	neurology,AFR,case_epilepsy
RD_001	/data/rd/001.vcf.gz	male	wgs	rare_disease,SAS,case_LQTS
```

Key: use cohort names (`cardiology`, `neurology`, `rare_disease`) as phenotype codes alongside disease-specific and ancestry labels.

### Querying by Cohort

```bash
# AF in cardiology cohort only
afquery query --db ./merged/ --locus chr1:12345678 --ref C --alt T \
  --phenotype cardiology

# AF in everyone except rare disease
afquery query --db ./merged/ --locus chr1:12345678 --ref C --alt T \
  --phenotype ^rare_disease

# AF in European subset across all cohorts
afquery query --db ./merged/ --locus chr1:12345678 --ref C --alt T \
  --phenotype EUR
```

### When to Use

- All cohorts share the same genome build and VCF pipeline
- You want cross-cohort AF queries without scripting
- Data governance allows combining samples
- Phenotype code design can capture all relevant groupings

### Trade-offs

| Advantage | Disadvantage |
|-----------|-------------|
| Single database to maintain | Rebuilding requires all VCFs accessible |
| Cross-cohort queries via phenotype filters | Phenotype code design must be planned upfront |
| One annotation pass covers all cohorts | Adding a new cohort requires `afquery update` |
| Flexible ad-hoc stratification | Larger database, longer rebuild time |

---

## Pattern 3: Tiered Approach

Maintain both per-cohort and merged databases.

```
/databases/
  institutional/
    cardiology_cohort/     ← institutional, restricted access
    neurology_cohort/      ← institutional, restricted access
  shared/
    combined_controls/     ← merged control samples, broader access
```

### When to Use

- Some cohorts have access restrictions that prevent full merging
- You need a shared "reference panel" of controls from multiple sources
- Institutional databases are updated independently on different schedules

### Implementation

1. Each institution maintains its own database
2. Control samples (or a consented subset) are merged into a shared database
3. Clinical queries annotate against both institutional and shared databases

```bash
# Annotate against institutional cohort
afquery annotate --db /databases/institutional/cardiology/ \
  --input patient.vcf.gz --output step1.vcf.gz

# Annotate against shared controls (use a different INFO prefix via Python API)
```

---

## Decision Matrix

| Factor | Pattern 1 (Separate) | Pattern 2 (Merged) | Pattern 3 (Tiered) |
|--------|---------------------|--------------------|--------------------|
| Data governance | Easiest | Requires agreement | Flexible |
| Cross-cohort queries | Scripting required | Built-in | Partial |
| Rebuild cost | Per-cohort only | All samples | Both |
| Storage | Proportional | Slightly less | Higher (duplication) |
| Maintenance complexity | Low per-database | Low (one database) | Higher |
| Best for | Multi-institution | Single organization | Mixed governance |

---

## Next Steps

- [Create a Database](../guides/create-database.md) — database creation options
- [Update a Database](../guides/update-database.md) — adding samples to an existing database
- [Sample Filtering](../guides/sample-filtering.md) — phenotype and technology filters
- [Pipeline Integration](pipeline-integration.md) — using databases in automated pipelines
