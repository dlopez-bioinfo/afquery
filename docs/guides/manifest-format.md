# Manifest Format

The manifest is a tab-separated (TSV) file that describes your sample cohort. It is the primary input to `afquery create-db` and `afquery update-db --add-samples`.

---

## Column Specification

| Column | Required | Type | Description |
|--------|----------|------|-------------|
| `sample_name` | Yes | string | Unique identifier for the sample. Must be unique across the entire database. |
| `vcf_path` | Yes | string | Absolute or relative path to the single-sample VCF file (plain or `.gz`). |
| `sex` | Yes | `male` \| `female` | Biological sex. Used for ploidy-aware AN computation on sex chromosomes. |
| `tech` | Yes | string | Sequencing technology name (e.g., `wgs`, `wes_v1`, `capture_kit_A`). Case-sensitive. |
| `phenotype` | No | string | Comma-separated ICD codes for this sample (e.g., `E11.9,I10`). Empty = no phenotype. |

---

## Example

```tsv
sample_name	vcf_path	sex	tech	phenotype
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

- Use standard **ICD-10 codes** (e.g., `E11.9`, `I10`, `N18.3`)
- Multiple codes per sample: comma-separated, no spaces (`E11.9,I10`)
- Samples with no phenotype: leave the `phenotype` column empty
- Codes are stored as-is and matched exactly in queries (case-sensitive)

---

## Technology Names

- Technology names are arbitrary strings — use whatever is meaningful for your cohort
- Common conventions: `wgs`, `wes_v1`, `wes_v2`, `capture_twist_v2`
- **WGS technologies** (no BED file): all positions are always eligible
- **WES technologies** (with BED file): coverage is determined by `<tech>.bed` in `--bed-dir`

!!! important "WES BED files"
    For any technology that is not WGS, you must provide a BED file named `<tech>.bed` in the `--bed-dir` directory. If no BED file is found, AFQuery treats the technology as WGS (all positions covered).

---

## BED File Format

BED files for WES capture regions must be:

- **0-based, half-open** coordinates (standard BED format)
- Tab-separated with at least 3 columns: `chrom`, `start`, `end`
- Chromosome names matching your VCF style (`chr1` or `1`)
- Named `<tech_name>.bed` (exact match to the `tech` column in the manifest)

Example `wes_v1.bed`:
```
chr1	65419	65433
chr1	925952	926117
chr1	1014541	1015072
```

---

## Common Mistakes

| Mistake | Effect | Fix |
|---------|--------|-----|
| Spaces instead of tabs | Parse error | Use `\t` (Tab key) as separator |
| Duplicate `sample_name` | Ingest error | Each sample must have a unique name |
| Wrong sex value | Silent error (no ploidy adjustment) | Use exactly `male` or `female` |
| Spaces in phenotype (`E11.9, I10`) | Code `" I10"` stored with leading space | Use `E11.9,I10` (no spaces) |
| Relative VCF path from wrong CWD | File not found error | Use absolute paths or run from the correct directory |
| Missing BED for WES tech | Treated as WGS (all positions covered) | Place `<tech>.bed` in `--bed-dir` |

---

## Template

Download a blank manifest template:

```tsv
sample_name	vcf_path	sex	tech	phenotype
```

Or generate one with headers only:

```bash
echo -e "sample_name\tvcf_path\tsex\ttech\tphenotype" > manifest.tsv
```
