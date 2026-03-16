# FILTER=PASS Tracking

Schema version 2 introduces tracking of variants that are called but fail quality filters (`FILTER≠PASS`). This allows distinguishing between "variant not present" and "variant present but low quality" in your cohort.

---

## Background: VCF FILTER Field

In VCF format, the FILTER column indicates whether a variant call passed quality filters:
- `PASS` or `.` (missing) — the variant passed all filters
- Any other value (e.g., `LowQual`, `VQSRTrancheSNP99.90to100.00`) — the variant failed one or more filters

AFQuery default behavior:
- **PASS-only** (default): only `FILTER=PASS` variants are counted in AC/AN
- **All filters** (`--include-all-filters`): all variant calls are counted regardless of FILTER status

---

## Schema v2: fail_bitmap

In schema v2, AFQuery stores a third bitmap per variant alongside `het_bitmap` and `hom_bitmap`:

- **`fail_bitmap`** — bit set for each sample that has a non-ref genotype (AC>0) AND `FILTER≠PASS`

This means:
- A sample in `fail_bitmap` was genotyped with the alt allele but the call failed QC
- Such samples are **not** counted in AC/AN (they don't affect AF)
- Their count is exposed as `FAIL_SAMPLES` / `N_FAIL`

---

## Enabling Schema v2

Schema v2 is the default for new databases. The `fail_bitmap` is always written; its content depends on the `--include-all-filters` flag:

```bash
# Default: PASS-only ingestion (fail_bitmap tracks failed calls)
afquery create-db --manifest manifest.tsv --output-dir ./db/ --genome-build GRCh38

# Include all filters: AC/AN counts all calls; fail_bitmap is all-zeros
afquery create-db --manifest manifest.tsv --output-dir ./db/ --genome-build GRCh38 \
  --include-all-filters
```

---

## Querying FAIL_SAMPLES

### CLI

The `FAIL_SAMPLES` count is shown automatically in query output for v2 databases:

```bash
afquery query --db ./db/ --chrom chr1 --pos 925952
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
    if r.N_FAIL is not None and r.N_FAIL > 0:
        print(f"  Warning: {r.N_FAIL} samples have low-quality calls at this site")
```

`N_FAIL` is `0` for v1 databases (no fail tracking), not `None`. Use `QueryResult.N_FAIL` directly.

---

## VCF Annotation

Schema v2 databases add an additional INFO field to annotated VCFs:

| Field | Type | Description |
|-------|------|-------------|
| `AFQUERY_N_FAIL` | Integer | Eligible samples with FILTER≠PASS at this variant |

```bash
afquery annotate --db ./db/ --input variants.vcf --output annotated.vcf
```

v1 databases do not write `AFQUERY_N_FAIL`.

---

## Migration from v1

v1 databases (no `schema_version` in `manifest.json`, or `schema_version < 2.0`) are fully compatible with AFQuery v2 code. They return `N_FAIL=0` in query results and do not write `AFQUERY_N_FAIL` in annotation.

There is no automatic migration from v1 to v2. To gain fail tracking, rebuild the database with `create-db`.

---

## When to Use --include-all-filters

By default (PASS-only), AF reflects the quality-filtered allele frequency — the frequency of the alt allele among high-quality calls. This is the recommended setting for most clinical and research use cases.

Use `--include-all-filters` when:
- You want to count all called genotypes regardless of quality
- You are studying filter performance or variant calling artifacts
- Your upstream pipeline handles filtering separately

With `--include-all-filters`, `fail_bitmap` stores zeros for all variants (since all calls are counted in AC/AN).
