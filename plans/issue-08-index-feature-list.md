# Issue 08 — Index Page Feature List

## Summary

The `docs/index.md` feature list is incomplete (feature 6 lacks a description) and missing two features (VCF annotation and changelog). The ordering of features should also be revised to lead with the most impactful capabilities.

## Source notes

| File | Line | Note |
|------|------|------|
| `docs/index.md` | 61 | Feature 6 is present but has no description; add a description explaining how AFQuery handles cohorts with mixed sequencing technologies and correctly computes AN per captured region. |
| `docs/index.md` | 63 | Add **Feature 7**: VCF annotation with custom sample subsets. |
| `docs/index.md` | 65 | Add **Feature 8**: The database includes a changelog to ensure reproducibility. |
| `docs/index.md` | 67 | Reorder features: 1, 3, 6, 5, 7, 8, 2, 4. |

## Scope

Docs only — `docs/index.md`.

## Steps

### 1. Add description to Feature 6

Current feature 6 is about mixed-technology cohorts. Add a concise description explaining that AFQuery correctly computes AN by intersecting each sample's capture BED with the queried position, even when the cohort mixes WGS, WES kits, and gene panels (including different versions of the same kit).

### 2. Add Feature 7 — VCF annotation

```markdown
7. **VCF annotation** with custom sample subsets — annotate a patient VCF with `AFQUERY_AC`, `AFQUERY_AN`, and `AFQUERY_AF` INFO fields computed from any combination of phenotype, sex, and technology filters.
```

### 3. Add Feature 8 — Changelog

```markdown
8. **Audit changelog** — every database operation (sample add, remove, or metadata update) is recorded in a tamper-evident changelog, ensuring result reproducibility.
```

### 4. Reorder features

Apply the reordering 1 → 3 → 6 → 5 → 7 → 8 → 2 → 4 to the feature list in `docs/index.md`.

The rationale for this order: lead with the core capability (allele frequency query), then the key differentiators (ploidy-aware AN, mixed-technology AN), then the clinical use cases (sample filtering), then the operational tools (annotation, changelog), and close with the internals (bitmap engine, database operations).

### 5. Verification

```bash
mkdocs build
mkdocs serve  # visually verify the index page feature list
```
