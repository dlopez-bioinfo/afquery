# VCF Preprocessing

AFQuery ingests single-sample, normalized VCF files. While AFQuery itself does not perform VCF normalization, the accuracy of your allele frequency estimates depends on the quality and consistency of your input VCFs. This page explains why normalization matters and provides a reference pipeline.

A ready-to-use normalization script is provided at [`resources/normalize_vcf.sh`](https://github.com/dlopez-bioinfo/afquery/blob/master/resources/normalize_vcf.sh).

---

## Why Normalize?

### Left-alignment and decomposition

Multi-allelic variants and complex indels can be represented in multiple equivalent ways in VCF format. Without normalization:

- The same deletion may appear as different ALT alleles depending on the caller
- Multi-allelic sites may not be decomposed into biallelic records
- Duplicate records can inflate AC

`bcftools norm` left-aligns indels against the reference genome and decomposes multi-allelic sites, ensuring consistent representation across samples.

### Ploidy correction for sex chromosomes

Male samples should have haploid genotype calls at chrX non-PAR regions (GT=`1`, not GT=`0/1` or `1/1`). Many variant callers output diploid genotypes for males at chrX by default. Without ploidy correction:

- Males at chrX are counted as diploid → AN is inflated by up to 50%
- AF estimates for X-linked variants are systematically underestimated

`bcftools +fixploidy` corrects genotype ploidy based on a sex file.

### Filtering and stripping

- **Non-PASS genotypes**: Masking them as missing ensures that low-quality calls do not inflate AC (AFQuery can track these as FAIL_SAMPLES in schema v2)
- **Homozygous reference calls**: Removing ref/ref genotypes reduces file size and speeds ingestion; they contribute AC=0 and are not needed
- **INFO fields**: Stripping INFO reduces file size; AFQuery uses only the GT field

---

## Reference Pipeline (bcftools)

The following pipeline processes a single-sample VCF into an AFQuery-ready normalized VCF. It requires [bcftools](https://samtools.github.io/bcftools/) ≥ 1.15.

```bash
# Variables
VCF_ID="sample_001"          # Used for output file naming
VCF="sample_001.vcf.gz"      # Input VCF (single-sample)
REF="reference.fa"           # Reference FASTA (must match VCF genome build)
GENDER="male"                # M/F/male/female
THREADS=4
OUT="${VCF_ID}_norm.vcf.gz"

# Step 1: Create chromosome rename map (only needed once per project)
# Maps chr-prefixed names (chr1, chrX) to non-prefixed (1, X)
# AFQuery's normalize_chrom() handles both styles, but consistency is important
{
  echo "chrM MT"
  for i in {1..22}; do echo "chr${i} ${i}"; done
  echo "chrX X"
  echo "chrY Y"
} > chr_list.txt

# Step 2: Create sex file for ploidy correction
SAMPLE_ID=$(bcftools query -l "${VCF}" | head -n1)
echo "${SAMPLE_ID} M" > "${VCF_ID}_gender.txt"   # or F for female

# Step 3: Run normalization pipeline
bcftools annotate \
    -x INFO,^FORMAT/GT \       # Strip INFO fields; keep only GT
    --force \
    --rename-chrs chr_list.txt \
    "${VCF}" | \
bcftools norm \
    -d exact \                  # Remove exact duplicate records
    -d both \                   # Also remove non-exact duplicates (same pos, diff representation)
    -f "${REF}" \               # Reference FASTA for left-alignment
    --check-ref ws \            # Warn on ref mismatch; skip (not error)
    --targets 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT | \
bcftools +setGT -- \
    -n . \                      # Set to missing
    -i 'FILTER!="PASS"' \       # ...any genotype where FILTER is not PASS
    -t q | \
bcftools +fixploidy -- \
    -s "${VCF_ID}_gender.txt" | \ # Correct chrX/Y ploidy for this sample's sex
bcftools view \
    -e 'GT="0" || GT="0/0"' \  # Remove homozygous reference calls
    --threads "${THREADS}" \
    -Oz \
    -o "${OUT}"

# Index
bcftools index --threads "${THREADS}" "${OUT}"
```

The convenience script `resources/normalize_vcf.sh` wraps this pipeline and can be called as:

```bash
./resources/normalize_vcf.sh sample_001 input.vcf.gz reference.fa male 8
```

---

## Input File Formats

### Sex file

The sex file specifies the sample ID and sex for `bcftools +fixploidy`:

```
sample_001 M
```

- One line per sample (use only the sample ID from the VCF header)
- Sex: `M` (male) or `F` (female)
- Separated by a single space

To get the sample ID from a VCF:

```bash
bcftools query -l sample.vcf.gz
```

### Chromosome rename map

Maps chromosome names in your VCF to the internal names AFQuery expects. If your VCF uses `chr1`-style names, use:

```
chrM MT
chr1 1
chr2 2
...
chr22 22
chrX X
chrY Y
```

If your VCF already uses `1`-style names (no `chr` prefix), you can skip the rename step. AFQuery's `normalize_chrom()` handles both styles transparently.

---

## When Can You Skip Steps?

| Step | When to skip |
|------|--------------|
| Chromosome rename | VCF already uses consistent naming and AFQuery normalizes it correctly |
| Left-align + decompose (`bcftools norm`) | VCF is guaranteed to be already normalized by the variant caller |
| Mask non-PASS genotypes (`setGT`) | All genotypes are PASS, or you want to include all calls regardless of filter |
| Ploidy correction (`+fixploidy`) | All samples are female, or cohort is WGS with correct diploid calls everywhere |
| Remove homozygous ref (`bcftools view`) | You want to preserve ref calls (increases file size; no accuracy impact) |

For quick experiments with a small cohort: steps 1 (strip INFO + rename) + 3 (norm) + 6 (remove ref) are the minimum needed.

---

## Common Pitfalls

**Missing chromosome rename map → unknown contigs**
If your VCF uses `chr1` naming but you skip the rename step and AFQuery does not find the chromosome, queries will return no results. Always verify chromosome names match your database build.

**Wrong sex file → inflated AN at chrX**
If a male sample is processed without ploidy correction, GT=`0/1` (diploid) is stored instead of GT=`1` (haploid). This contributes AN=2 instead of AN=1 for that sample at chrX non-PAR positions — inflating AN by up to 50% for male-heavy cohorts.

**Skipping deduplication → duplicate alleles inflate AC**
If the same variant appears twice in a VCF (common with certain multi-sample calling pipelines after splitting), both records are ingested, doubling AC for that sample. Use `-d exact -d both` with `bcftools norm`.

**Multi-sample VCF input**
The normalization script warns if the input VCF contains more than one sample and uses the first sample only. Always split multi-sample VCFs before normalization:

```bash
bcftools +split multi_sample.vcf.gz -Oz -o ./split_vcfs/
```

**Reference genome mismatch**
`--check-ref ws` (warn and skip) is used instead of `s` (error) to avoid aborting on minor reference mismatches. Review warnings in stderr output. Serious mismatches indicate the wrong reference genome is being used.

---

## Integration with AFQuery

Normalization is entirely external to AFQuery. After normalizing all samples, point the AFQuery manifest to the normalized VCF outputs:

```tsv
sample_name	vcf_path	sex	tech_name	phenotype_codes
sample_001	/data/normalized/sample_001_norm.vcf.gz	male	wgs	E11.9
```

Then proceed to [Create a Database](../guides/create-database.md).
