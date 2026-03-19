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

- **Non-PASS genotypes**: Masking them as missing ensures that low-quality calls do not inflate AC. AFQuery tracks these as FAIL_SAMPLES.
- **Homozygous reference calls**: Removing ref/ref genotypes reduces file size and speeds ingestion; they contribute AC=0 and are not needed
- **INFO fields**: Stripping INFO reduces file size; AFQuery uses only the GT field

NOTA_PARA_LLM: La eliminación del campo INFO también evita que VCFs malformados puedan romper el pipeline
