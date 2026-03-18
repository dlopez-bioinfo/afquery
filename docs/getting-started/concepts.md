# Key Concepts

## Mental Model

**AFQuery = bitmap index over genotypes.**

For each variant, AFQuery stores a compressed bitset recording which samples carry the alt allele. For each query, it intersects the carrier bitset with an eligible-sample bitset, counts bits, and returns AC/AN/AF.

```
Samples → Metadata → Filters → Eligible Bitmap
Variants → Genotypes → het/hom Bitmaps
Query: (eligible bitmap) AND (carrier bitmap) → AC; count(eligible) → AN; AC/AN → AF
```

```mermaid
flowchart LR
    subgraph Input
        S["Samples + Metadata"]
        V["Variants + Genotypes"]
    end

    subgraph Bitmaps
        E["Eligible Bitmap<br/>(sex, phenotype, tech filters)"]
        C["Carrier Bitmap<br/>(het | hom per variant)"]
    end

    subgraph Result
        R["AC = popcount(eligible AND carrier)<br/>AN = popcount(eligible) × ploidy<br/>AF = AC / AN"]
    end

    S --> E
    V --> C
    E --> R
    C --> R

    style E fill:#e3f2fd
    style C fill:#f3e5f5
    style R fill:#e8f5e9
```

---

## Why Cohort-Specific AF Matters

Allele frequency is not a fixed property of a variant — it is a property of a population.

Population databases like gnomAD aggregate data from hundreds of thousands of individuals and are invaluable for identifying common variants. But they are not your cohort. Allele frequencies vary substantially across:

- **Ancestry**: A variant at AF=0.001 in gnomAD Europeans may be at AF=0.01 in a specific Iberian registry — a 10× difference that changes pathogenicity interpretation under ACMG criteria.
- **Sequencing technology**: Coverage differences between WES and WGS can create systematic differences in observed genotype rates.
- **Clinical composition**: A cohort enriched for a specific disease may show elevated frequencies for variants associated with that phenotype.

AFQuery lets you compute allele frequencies on exactly the samples in your hands — and on any dynamically defined subset of them — without rebuilding the database.

---

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

### Worked Example

Consider a cohort of 60 samples (50 diploid autosomes) with a variant at position chr1:925952 G>A:

| Group | Eligible samples | AN | AC | AF |
|-------|-----------------|----|----|-----|
| Full cohort (all 60) | 60 samples | 120 | 3 | 0.025 |
| Females only (35) | 35 samples | 70 | 2 | 0.029 |
| Tagged `E11.9` (20) | 20 samples | 40 | 0 | 0.000 |
| Not tagged `E11.9` (40) | 40 samples | 80 | 3 | 0.038 |

**Key point**: AN is not always `2 × cohort_size`. Eligible samples change per query, so AN reflects your chosen subgroup exactly.

### Visualization

The same variant can show dramatically different allele frequencies across subgroups:

```mermaid
graph TB
    subgraph Full["Full Cohort (60 samples)"]
        F["AC=3, AN=120, AF=0.025"]
    end

    subgraph Females["Females Only (35 samples)"]
        FEM["AC=2, AN=70, AF=0.029"]
    end

    subgraph Disease["Tagged E11.9 (20 samples)"]
        DIS["AC=0, AN=40, AF=0.000"]
    end

    subgraph Controls["Without E11.9 (40 samples)"]
        CTL["AC=3, AN=80, AF=0.038"]
    end

    Full -->|filter| Females
    Full -->|stratify| Disease
    Full -->|exclude| Controls

    style F fill:#e3f2fd
    style FEM fill:#f3e5f5
    style DIS fill:#ffebee
    style CTL fill:#e8f5e9
```

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
- Metadata/phenotype codes (arbitrary strings, comma-separated)

The manifest is parsed into `metadata.sqlite` during `create-db`. The original path is recorded in `manifest.json`.

---

## Sample Filtering Model

Queries can restrict the eligible sample set along three independent dimensions:

| Dimension | Filter | Default |
|-----------|--------|---------|
| **Sex** | `male`, `female`, or `both` | `both` |
| **Phenotype** | Include/exclude metadata codes (arbitrary strings) | all samples |
| **Technology** | Include/exclude tech names | all samples |

Filters **compose with AND** across dimensions: a sample must satisfy all three to be eligible.

Within a dimension, multiple include codes compose with **OR** (a sample matching any code is included).

AN is computed only over eligible samples, so AF naturally reflects the chosen subgroup.

### The Metadata Model

**AFQuery treats phenotype codes as arbitrary string labels.** You can use:

- ICD-10 disease codes (`E11.9`, `G40`)
- HPO phenotype terms (`HP:0001250`)
- GO terms
- Custom project tags (`control`, `rare_disease`, `pilot_cohort`)
- Technology subgroups (`panel_v1`, `panel_v2`)
- Any combination of the above

Multiple labels per sample are supported. There is no validation or controlled vocabulary — you define the ontology for your cohort.

When planning phenotype codes, consider:
- Codes are immutable after ingestion (changing requires removing and re-adding samples)
- Codes are case-sensitive: `E11.9` ≠ `e11.9`
- Trailing or leading spaces cause silent mismatch (always use `E11.9,I10`, never `E11.9, I10`)

---

## Schema Versions

| Version | Feature |
|---------|---------|
| v1 | het_bitmap + hom_bitmap only |
| v2 | Adds fail_bitmap; AFQUERY_FAIL INFO field in VCF annotation |

PASS-only ingestion is always enforced by default. See [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md) for details.
