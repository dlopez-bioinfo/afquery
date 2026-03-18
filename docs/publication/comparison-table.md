# Comparison: AFQuery vs. Alternative Tools

This table summarizes how AFQuery compares to alternative approaches for computing cohort allele frequencies.

---

## Feature Comparison

| Feature | AFQuery | bcftools stats | VCFtools --freq | GATK GenomicsDB | Hail |
|---------|---------|---------------|----------------|----------------|------|
| **Query latency** | <100 ms | Minutes (VCF scan) | Minutes (VCF scan) | Seconds | Seconds–minutes |
| **Dynamic subcohort queries** | Yes | No | No | Partial | Yes (programmatic) |
| **Metadata filtering** | Arbitrary labels | No | No | No | User-defined |
| **Sex-stratified AF** | Yes (auto ploidy) | Manual | Manual | No | Manual |
| **Technology-aware AN** | Yes (BED capture) | No | No | No | No |
| **Incremental updates** | Yes (no rebuild) | N/A | N/A | Yes (import) | Rebuild |
| **Infrastructure required** | None (file-based) | None | None | Java/server | Spark cluster |
| **Input format** | Single-sample VCFs | Any VCF | VCF | gVCF | VCF/BGEN |
| **Output format** | JSON/TSV/VCF annotation | Stats text | Freq file | Merged gVCF | Table/VCF |
| **Max cohort size (tested)** | 50,000 | 100,000+ | 100,000+ | 100,000+ | 1,000,000+ |
| **Bitmap compression** | Yes (Roaring) | No | No | Yes (GenomicsDB) | No |
| **VCF annotation** | Yes | No | No | No | Yes |
| **Python API** | Yes | No | No | No | Yes |

---

## Use Case Match

| Use case | Recommended tool | Reason |
|----------|-----------------|--------|
| Interactive variant interpretation (<1s required) | **AFQuery** | Sub-100 ms dynamic queries |
| Cohort-specific rare disease AF | **AFQuery** | Metadata-aware subcohort queries |
| gnomAD-style static population AF export | bcftools stats or VCFtools | Optimized for bulk, fixed-cohort computation |
| Joint genotyping pipeline (gVCF → VCF) | GATK GenomicsDB | Purpose-built for GVCF consolidation |
| Population-scale GWAS analysis | Hail | Designed for million-sample statistical genetics |
| X-linked hemizygous frequency | **AFQuery** | Automatic ploidy-aware AN |

---

## Performance Comparison (50K samples, chr1, point query)

| Tool | Setup time | Per-query time | Subcohort requery |
|------|-----------|---------------|------------------|
| AFQuery | ~13 min (52 cores) | <100 ms cold, ~10 ms warm | <100 ms (no rebuild) |
| bcftools stats | Merge: ~2h, Index: ~30 min | ~5 min (full scan) | ~5 min (new subset VCF required) |
| VCFtools --freq | Same as bcftools | ~5 min | ~5 min |

*Note: Merge/index times for bcftools/VCFtools are one-time costs, but subcohort requery requires extracting a new VCF subset.*

---

## Limitations of AFQuery

AFQuery is purpose-built for fast subcohort AF computation and is not a general-purpose genomic database:

- **Not a joint genotyper**: AFQuery does not perform joint genotyping. Input VCFs should be individually called before ingestion.
- **Not a variant database**: AFQuery stores only genotype-level summaries (bitmaps). Individual sample genotypes cannot be retrieved from the database.
- **No statistical genetics**: AFQuery does not compute Hardy-Weinberg equilibrium, population stratification, or other statistical genetics metrics.
- **Single-chromosome queries**: Each query call operates on one chromosome. Multi-chromosome batch queries require multiple calls (easily parallelized with the Python API).
- **Cohort size limit**: Performance at >100K samples has not been validated. Memory requirements for the build phase scale with cohort size.

---

## Citation

If comparing AFQuery in a publication, cite the original AFQuery paper *(citation to be added upon publication)*.

For bcftools: Danecek et al., *GigaScience*, 2021 (doi:10.1093/gigascience/giab008).
For VCFtools: Danecek et al., *Bioinformatics*, 2011 (doi:10.1093/bioinformatics/btr330).
For GATK/GenomicsDB: Van der Auwera & O'Connor, *Current Protocols*, 2020.
For Hail: Hail Team, https://hail.is.
