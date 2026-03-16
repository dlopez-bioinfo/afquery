# Key Concepts

## Allele Frequency (AC / AN / AF)

For a given variant (chromosome, position, ref, alt):

- **AC** (Allele Count) — number of copies of the alt allele observed across eligible samples
- **AN** (Allele Number) — total number of alleles examined (2× diploid samples, 1× haploid)
- **AF** (Allele Frequency) — `AC / AN` (None when AN = 0)

A heterozygous carrier contributes AC=1; a homozygous alt carrier contributes AC=2.

### Ploidy

AN depends on the chromosome and sex of eligible samples:

| Chromosome | Female | Male |
|------------|--------|------|
| Autosomes (chr1–22) | 2 | 2 |
| chrX (non-PAR) | 2 | 1 |
| chrX (PAR1/PAR2) | 2 | 2 |
| chrY | 0 | 1 |
| chrMT | 1 | 1 |

See [Ploidy & Sex Chromosomes](../advanced/ploidy-and-sex-chroms.md) for PAR coordinates.

---

## How AFQuery Stores Data

### Per-Variant Bitmaps

For each variant row, AFQuery stores two [Roaring Bitmaps](https://roaringbitmap.org/):

- **`het_bitmap`** — bit set for each sample that is heterozygous (GT=0/1 or 1/0)
- **`hom_bitmap`** — bit set for each sample that is homozygous alt (GT=1/1)
- **`fail_bitmap`** *(schema v2 only)* — bit set for each sample called with FILTER≠PASS

Each sample has a stable integer ID (0-indexed). The bit position in the bitmap equals the sample ID.

### Parquet Storage

Bitmaps are serialized and stored in Parquet files, partitioned by chromosome and 1-Mbp bucket:

```
variants/
  chr1/
    bucket_0/   ← positions 0–999,999
    bucket_1/   ← positions 1,000,000–1,999,999
    ...
  chr2/
    ...
```

Rows within each bucket are sorted by `(pos, alt)`.

### Capture Index (WES)

For whole-exome sequencing (WES) technologies, a BED file defines covered regions. AFQuery builds an interval tree (pickle file) per technology so queries can determine which WES samples are eligible at any given position.

```
capture/
  wes_v1.pkl   ← interval tree for wes_v1 BED
  wes_v2.pkl
```

WGS samples are always eligible (no BED file needed).

---

## The Manifest

The manifest is a TSV file that drives database creation. It maps each sample to its:
- VCF file path
- Sex (`male` / `female`)
- Sequencing technology
- Phenotype codes (ICD-10, comma-separated)

The manifest is parsed into `metadata.sqlite` during `create-db`. The original path is recorded in `manifest.json`.

---

## Sample Filtering Model

Queries can restrict the eligible sample set along three independent dimensions:

| Dimension | Filter | Default |
|-----------|--------|---------|
| **Sex** | `male`, `female`, or `both` | `both` |
| **Phenotype** | Include/exclude ICD codes | all samples |
| **Technology** | Include/exclude tech names | all samples |

Filters **compose with AND** across dimensions: a sample must satisfy all three to be eligible.

Within a dimension, multiple include codes compose with **OR** (a sample matching any code is included).

AN is computed only over eligible samples, so AF naturally reflects the chosen subgroup.

---

## Schema Versions

| Version | Feature |
|---------|---------|
| v1 | het_bitmap + hom_bitmap only |
| v2 | Adds fail_bitmap; AFQUERY_FAIL INFO field in VCF annotation |

Use `--include-all-filters` on `create-db` to ingest all variants regardless of FILTER status (default: PASS-only).

See [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md) for details.
