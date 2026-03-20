# Manifest Format

The manifest is a tab-separated (TSV) file that describes your sample cohort. It is the primary input to `afquery create-db` and `afquery update-db --add-samples`.

---

## Column Specification

| Column | Required | Type | Description |
|--------|----------|------|-------------|
| `sample_name` | Yes | string | Unique identifier for the sample. Must be unique across the entire database. |
| `vcf_path` | Yes | string | Absolute or relative path to the single-sample VCF file (plain or `.gz`). |
| `sex` | Yes | `male` \| `female` | Biological sex. Used for ploidy-aware AN computation on sex chromosomes. |
| `tech_name` | Yes | string | Sequencing technology name (e.g., `wgs`, `wes_v1`, `capture_kit_A`). Case-sensitive. |
| `phenotype_codes` | No | string | Comma-separated phenotype codes for this sample (arbitrary strings, e.g., `E11.9,control`). |

---

## Example

```tsv
sample_name	vcf_path	sex	tech_name	phenotype_codes
SAMP_001	/data/vcfs/SAMP_001.vcf.gz	female	wgs	E11.9,I10
SAMP_002	/data/vcfs/SAMP_002.vcf.gz	male	wgs	E11.9
SAMP_003	/data/vcfs/SAMP_003.vcf.gz	female	wes_v1	I10
SAMP_004	/data/vcfs/SAMP_004.vcf.gz	male	wes_v1
SAMP_005	/data/vcfs/SAMP_005.vcf.gz	female	wgs	E11.9,I10,N18.3
```

!!! note "Header row required"
    The first row must be the header with exactly these column names.

---

## Phenotype Codes

The `phenotype_codes` column accepts **any arbitrary strings** — there is no required ontology. Common choices include:

- ICD-10 disease codes (e.g., `E11.9`, `I10`, `N18.3`)
- HPO phenotype terms (e.g., `HP:0001250`)
- Custom project tags (e.g., `control`, `pilot`, `high_coverage_wgs`)
- Population or batch labels (e.g., `EUR`, `batch_2023`)
- Technology subgroups (e.g., `panel_v1`, `panel_v2`)

Multiple codes per sample can be provided (comma-separated, no spaces)
Samples with no phenotype can leave the `phenotype_codes` column empty


!!! important phenotype filtering
    `phenotype_codes` are stored as-is and matched exactly in queries (case-sensitive)

---

## Technology Names

- Technology names are arbitrary strings — use whatever is meaningful for your cohort
- Common conventions: `wgs`, `wes_v1`, `wes_v2`, `capture_twist_v2`
- **WGS technologies** (no BED file): all positions are always eligible
- **WES technologies** (with BED file): coverage is determined by `<tech>.bed` in `--bed-dir`

!!! important "WES BED files"
    For any technology that is not WGS, you must provide a BED file named `<tech>.bed` in
    the `--bed-dir` directory. If the file is not found, AFQuery aborts with an error before
    any ingestion begins:

    ```
    ManifestError: BED file not found for technology 'WES_kit_A': '/path/to/WES_kit_A.bed'
    ```

    To resolve this, either supply the correct BED file or rename the technology to `WGS`
    in the manifest (no BED file required for WGS).

---

## BED File Format

BED files for WES capture regions must be:

- **0-based, half-open** coordinates (standard BED format)
- Tab-separated with at least 3 columns: `chrom`, `start`, `end`
- Chromosome names matching your VCF style (`chr1` or `1`)
- Named `<tech_name>.bed` (exact match to the `tech_name` column in the manifest)

---

## Common Mistakes

| Mistake | Effect | Fix |
|---------|--------|-----|
| Spaces instead of tabs | Parse error | Use `\t` (Tab key) as separator |
| Duplicate `sample_name` | Ingest error | Each sample must have a unique name |
| Wrong sex value | Silent error (no ploidy adjustment) | Use exactly `male` or `female` |
| Spaces in phenotype (`E11.9, I10`) | Code `" I10"` stored with leading space | Use `E11.9,I10` (no spaces) |
| Relative VCF path from wrong CWD | File not found error | Use absolute paths or run from the correct directory |

---
