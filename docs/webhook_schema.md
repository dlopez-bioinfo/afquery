# Webhook / API JSON Schema

## Batch Query

### Request

```json
{
  "chrom":     "chr1",
  "positions": [1000, 2000, 3000],
  "icd10":     ["E11.9", "I10"],
  "sex":       "both"
}
```

| Field       | Type            | Required | Notes                                      |
|-------------|-----------------|----------|--------------------------------------------|
| `chrom`     | string          | yes      | Any form accepted: `"1"`, `"chr1"`, `"X"` |
| `positions` | array of int    | yes      | 1-based positions; duplicates are ignored  |
| `icd10`     | array of string | yes      | ICD-10 codes; union semantics              |
| `sex`       | string          | no       | `"both"` (default), `"male"`, `"female"`  |

### Response

Array of variant objects, one per variant found (absent variants are omitted). Sorted by `(pos, alt)`.

```json
[
  {
    "chrom":      "chr1",
    "pos":        1000,
    "ref":        "A",
    "alt":        "T",
    "AC":         5,
    "AN":         100,
    "AF":         0.05,
    "n_eligible": 50
  }
]
```

| Field        | Type   | Notes                                        |
|--------------|--------|----------------------------------------------|
| `chrom`      | string | Canonical form (`chr`-prefixed)              |
| `pos`        | int    | 1-based                                      |
| `ref`        | string | Reference allele                             |
| `alt`        | string | Alternate allele (one object per alt)        |
| `AC`         | int    | Allele count in eligible samples             |
| `AN`         | int    | Allele number (ploidy-aware)                 |
| `AF`         | float  | `AC / AN`                                    |
| `n_eligible` | int    | Number of eligible samples at this position  |

Positions where `AN == 0` (no eligible samples covered) are silently omitted.

---

## Single-Position Query

Same as batch query with a single-element `positions` array.

---

## VCF Annotation

### Request

`POST /annotate`  — multipart form data.

| Part     | Content-Type               | Description                       |
|----------|----------------------------|-----------------------------------|
| `vcf`    | `application/octet-stream` | Input VCF file (plain or gzipped) |
| `params` | `application/json`         | JSON object (see below)           |

`params` object:

```json
{
  "icd10": ["E11.9"],
  "sex":   "both"
}
```

### Response

Annotated VCF file stream (`Content-Type: text/plain`).

Added INFO fields:

| Field      | Number | Type    | Description                          |
|------------|--------|---------|--------------------------------------|
| `AFQUERY_AC`  | A      | Integer | Allele count per ALT in eligible set |
| `AFQUERY_AN`  | 1      | Integer | Allele number (0 if uncovered)       |
| `AFQUERY_AF`  | A      | Float   | Allele frequency per ALT             |

**Notes**:
- `AFQUERY_AN=0` means the position was not covered for the given filters; `AFQUERY_AC` and `AFQUERY_AF` are absent (`.`) in that case.
- `AFQUERY_AC=0`, `AFQUERY_AF=0.0` with `AFQUERY_AN>0` means the position was covered but the specific allele was not observed in the database.
- For multi-allelic sites, `AFQUERY_AC` and `AFQUERY_AF` are comma-separated lists aligned to the `ALT` column.

### Completion summary (trailing header)

After the VCF body, a summary is written to stderr (or returned in a trailing `X-AFQUERY-Stats` response header):

```json
{"n_variants": 1000, "n_annotated": 850, "n_uncovered": 12}
```
