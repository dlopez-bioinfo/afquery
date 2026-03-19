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
- Their count is exposed as `FAIL_SAMPLES` / `N_FAIL`

---

## Database Creation

The `fail_bitmap` is always written and always tracks PASS-only ingestion:

```bash
# PASS-only ingestion (fail_bitmap tracks failed calls)
afquery create-db --manifest manifest.tsv --output-dir ./db/ --genome-build GRCh38
```

Support for ingesting all variants regardless of filter status is planned for a future release. Currently, PASS-only ingestion is always enforced.

---

## Querying FAIL_SAMPLES

### CLI

The `FAIL_SAMPLES` count is shown automatically in query output:

```bash
afquery query --db ./db/ --locus chr1:925952
```

```
chr1:925952
  REF=G  ALT=A  AC=142  AN=2742  AF=0.0518  N_HET=138  N_HOM_ALT=2  FAIL=7
```

`FAIL=7` means 7 eligible samples had the alt allele called but with FILTER≠PASS.

### Python API

```python
results = db.query("chr1", pos=925952)
for r in results:
    print(f"AC={r.AC}  AN={r.AN}  FAIL={r.N_FAIL}")
    if r.N_FAIL > 0:
        print(f"  Warning: {r.N_FAIL} samples have low-quality calls at this site")
```

`N_FAIL` is always an `int` (default `0`).

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

AF reflects the quality-filtered allele frequency — the frequency of the alt allele among high-quality calls. This is appropriate for most clinical and research use cases.

Currently, PASS-only ingestion is always enforced. Support for ingesting all variants regardless of filter status is planned for a future release.
