# FILTER=PASS Tracking

AFQuery tracks variants that are called but fail quality filters (`FILTER≠PASS`). This allows distinguishing between "variant not present" and "variant present but low quality" in your cohort.

---

## Background: VCF FILTER Field

In VCF format, the FILTER column indicates whether a variant call passed quality filters:

- `PASS` or `.` (missing) — the variant passed all filters
- Any other value (e.g., `LowQual`, `VQSRTrancheSNP99.90to100.00`) — the variant failed one or more filters

AFQuery's default behavior:

- Only `FILTER=PASS` calls contribute to **AC** (and therefore to **AF**). This is always enforced.
- **AN** is not affected by the filter — it counts every eligible sample at the position, failed calls included.

---

## fail_bitmap

AFQuery stores a third bitmap per variant alongside `het_bitmap` and `hom_bitmap`:

- **`fail_bitmap`** — bit set for each sample whose call at this site has `FILTER≠PASS`. In practice this is two cases: a sample that carries the alt allele but whose call failed a filter, and a sample with a missing genotype (`./.`) at a site that itself failed a filter.

What this means for a sample in `fail_bitmap`:

- Its alt alleles are **not** counted in AC, so it does not raise AF.
- It is still an eligible sample, so it **is** counted in AN. The failed call lowers AF by sitting in the denominator without contributing to the numerator.
- Its count is exposed separately as `N_FAIL`.

---

## Database Creation

The `fail_bitmap` is always written, regardless of build options:

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
chr1:925952 G>A  AC=142  AN=2742  AF=0.0518  n_eligible=1371  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1224  N_FAIL=7
```

`N_FAIL=7` means 7 eligible samples had a call with FILTER≠PASS at this site. They are part of `n_eligible` (and of AN), so the genotype counts still add up: 138 + 2 + 1224 + 7 = 1371.

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

Among the INFO fields `afquery annotate` writes, one reports the failed-call count:

| Field | Type | Description |
|-------|------|-------------|
| `AFQUERY_N_FAIL` | Integer | Eligible samples with `FILTER≠PASS` at this variant |

```bash
afquery annotate --db ./db/ --input variants.vcf --output annotated.vcf
```

See [Understanding Output](../getting-started/understanding-output.md#vcf-annotation-fields) for the full set of `AFQUERY_*` fields.

---

## PASS-Only Enforcement

PASS-only counting is always enforced and cannot be turned off. It applies to the numerator only: a failed call never adds an alt allele to AC, but the sample stays in the eligible set and so remains in AN. The result is a deliberately conservative AF — a failed carrier weighs on the denominator without lifting the numerator, which is the safe direction for clinical and research use.

---

## Next Steps

- [Understanding Output](../getting-started/understanding-output.md) — interpreting N_FAIL in query and annotate results
- [ACMG Criteria](../use-cases/acmg-use-cases.md) — using N_FAIL to assess site quality before applying BA1/PM2
- [Debugging Results](debugging-results.md) — diagnosing unexpectedly high N_FAIL values
