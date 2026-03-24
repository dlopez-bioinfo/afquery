# FILTER=PASS Tracking

AFQuery tracks variants that are called but fail quality filters (`FILTER≠PASS`). This allows distinguishing between "variant not present" and "variant present but low quality" in your cohort.

---

## Background: VCF FILTER Field

In VCF format, the FILTER column indicates whether a variant call passed quality filters:

- `PASS` or `.` (missing) — the variant passed all filters
- Any other value (e.g., `LowQual`, `VQSRTrancheSNP99.90to100.00`) — the variant failed one or more filters

AFQuery default behavior:

- **PASS-only**: only `FILTER=PASS` variants are counted in AC/AN. This is always enforced.

---

## fail_bitmap

AFQuery stores a third bitmap per variant alongside `het_bitmap` and `hom_bitmap`:

- **`fail_bitmap`** — bit set for each sample that has a non-ref genotype (AC>0) AND `FILTER≠PASS`

This means:

- A sample in `fail_bitmap` was genotyped with the alt allele but the call failed QC
- Such samples are **not** counted in AC/AN (they don't affect AF)
- Their count is exposed as `N_FAIL`

---

## Database Creation

The `fail_bitmap` is always written and always tracks PASS-only ingestion:

```bash
afquery create-db --manifest manifest.tsv --output-dir ./db/ --genome-build GRCh38
```

---

## Querying N_FAIL

### CLI

The `N_FAIL` count is shown automatically in query output:

```bash
afquery query --db ./db/ --locus chr1:925952
```

```
chr1:925952 G>A  AC=142  AN=2742  AF=0.0518  n_eligible=1371  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1231  N_FAIL=7
```

`N_FAIL=7` means 7 eligible samples had the alt allele called but with FILTER≠PASS.

### Python API

```python
results = db.query("chr1", pos=925952)
for r in results:
    print(f"AC={r.AC}  AN={r.AN}  FAIL={r.N_FAIL}")
    if r.N_FAIL > 0:
        print(f"  Warning: {r.N_FAIL} samples have low-quality calls at this site")
```

`N_FAIL` is always an `int` (default `0`).

### Identifying specific FAIL samples

To see which individual samples have FAIL status at a position, use `variant-info`:

```bash
afquery variant-info --db ./db/ --locus chr1:925952
```

Each carrier row shows its `filter` column as `PASS` or `FAIL`, along with sample metadata (technology, phenotype codes). This helps pinpoint whether failures cluster in a specific technology or sample group. See [Variant Info](../guides/variant-info.md) for full options.

---

## VCF Annotation

AFQuery adds an additional INFO field to annotated VCFs:

| Field | Type | Description |
|-------|------|-------------|
| `AFQUERY_N_FAIL` | Integer | Eligible samples with FILTER≠PASS at this variant |

```bash
afquery annotate --db ./db/ --input variants.vcf --output annotated.vcf
```

---

## PASS-Only Enforcement

AF reflects the quality-filtered allele frequency — the frequency of the alt allele among high-quality calls. This is appropriate for most clinical and research use cases. PASS-only ingestion is always enforced.

---

## Next Steps

- [Understanding Output](../getting-started/understanding-output.md) — interpreting N_FAIL in query and annotate results
- [ACMG Criteria](../use-cases/acmg-use-cases.md) — using N_FAIL to assess site quality before applying BA1/PM2
- [Debugging Results](debugging-results.md) — diagnosing unexpectedly high N_FAIL values
