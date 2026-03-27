# Performance Benchmarks

Characterize AFQuery query latency, build scalability, annotation throughput, and competitive performance.

## Experiments

| # | Script | Question |
|---|--------|----------|
| 1 | `02_query_scaling.py` | How does query latency scale with sample count? |
| 2 | `03_build.py` | How does build time scale with parallelism? |
| 3 | `04_annotate.py` | What is the VCF annotation throughput? |
| 4 | `05_vs_bcftools.py` | How does AFQuery compare to bcftools? |

## Running

All experiments are orchestrated by Snakemake from the `benchmarks/` root.

```bash
# Run the full performance benchmark
snakemake --cores 52 performance_all

# Dry run
snakemake --cores 52 --dry-run performance_all
```

## Configuration

Experiment parameters are in `config.py` (imports shared constants from `../shared/config.py`):

```python
# Synthetic scaling experiment
SYNTH_SCALES = [1_000, 5_000, 10_000, 25_000, 50_000]

# 1KG subsets
ONEKG_SUBSETS = [500, 1_000, 2_504]

# Experiment repetitions
QUERY_WARM_REPS = 50
BUILD_THREAD_COUNTS = [1, 4, 8, 16, 32]
BUILD_REPS = 3          # repetitions per (n_samples, threads)
ANNOTATE_THREAD_COUNTS = [1, 4, 8, 16, 32]
ANNOTATE_REPS = 3       # repetitions per (n_variants, threads)
BCFTOOLS_REPS = 10
```

## Experiment Details

### Experiment 1: Query Latency Scaling (`02_query_scaling.py`)

**Purpose:** Determine how AFQuery query latency scales as the number of indexed samples grows from 1,000 to 50,000. This is the central claim of the tool: that bitmap-indexed queries remain fast regardless of cohort size, unlike approaches that re-scan genotypes on every query.

**Setup:** Synthetic databases are built at 5 scales (1K, 5K, 10K, 25K, 50K samples), each containing 10,000 variants on chr22. Five query types are benchmarked against each database:

| Query type | Description | Reps |
|-----------|-------------|------|
| Point (cold) | Single-locus AF with a freshly loaded `Database` | 5 |
| Point (warm) | Single-locus AF with a pre-loaded `Database` (after warmup) | 50 |
| Region (1 Mbp) | All variants in a 1 Mbp window | 20 |
| Batch 100 | 100 variants in a single call | 20 |
| Batch 1000 | 1,000 variants in a single call | 20 |

Each query type is also tested under 4 filter configurations: no filter, sex=female, phenotype=E11.9, and combined (sex + phenotype). The primary figure uses "no filter".

**Figure 1 (`fig1_query_scaling`):** Log-log line plot. X-axis: number of samples. Y-axis: median latency in ms. Each line is a query type. Shaded bands show IQR (Q1-Q3).

**How to interpret:** If the lines are flat or near-flat as sample count increases, it confirms that bitmap operations scale sub-linearly. A steep slope would indicate a performance bottleneck. Cold vs. warm gap reveals the cost of initial Parquet file loading. Batch queries should show amortized cost (1000-variant batch should NOT be 1000x a point query).

---

### Experiment 2: Build Performance (`03_build.py`)

**Purpose:** Measure how long it takes to build an AFQuery database from VCFs, and how that time decreases with more threads. This tells users what hardware investment is needed to maintain the database.

**Setup:** Synthetic VCFs are generated at 3 sample scales (1K, 5K, 10K) with 10,000 variants on chr22. For each scale, the full `create-db` pipeline (ingest + consolidate + build) is timed across 5 thread counts (1, 4, 8, 16, 32). Each (scale, threads) combination is repeated `BUILD_REPS` times (default 3). Peak RSS memory is captured via `/usr/bin/time -v`. The resulting database directory size is recorded.

**Figure 2 (`fig2_build_perf`):** Grouped bar chart. X-axis: thread count. Bars grouped by sample scale (3 colors). Y-axis: wall-clock time in seconds. Error bars show IQR when multiple reps are available.

**How to interpret:** Bars should shrink as threads increase, demonstrating parallelism effectiveness. Diminishing returns at high thread counts indicate I/O or memory bottlenecks. Comparing across scales shows how build cost grows with cohort size.

**Figure 6 (`fig6_disk_footprint`):** Stacked bar chart comparing raw VCF size vs. AFQuery database size per scale. Shows the compression ratio of the bitmap encoding.

**How to interpret:** AFQuery databases should be significantly smaller than raw VCFs because they store only genotype bitmaps (not per-sample FORMAT fields). A smaller database also means faster queries (less I/O).

---

### Experiment 3: Annotation Throughput (`04_annotate.py`)

**Purpose:** Measure how many variants per second AFQuery can annotate in a VCF file, and how throughput scales with thread count. VCF annotation is the primary production workflow: a clinical lab submits a patient VCF and receives it back with AC/AN/AF annotations from the cohort database.

**Setup:** Input VCFs of 3 sizes (10K, 50K, 100K variants) are created by sampling from the largest 1KG database (2,504 samples). Each VCF is annotated using 5 thread counts (1, 4, 8, 16, 32). Each combination is repeated `ANNOTATE_REPS` times (default 3). The collect script computes speedup ratios relative to the single-threaded baseline.

**Figure 3 (`fig3_annotate_throughput`):** Line plot. X-axis: threads. Y-axis: throughput (variants/second). Each line is a VCF size.

**How to interpret:** Lines should increase roughly linearly with thread count up to some saturation point. Sub-linear scaling indicates contention (e.g., Parquet I/O serialization). Different VCF sizes should converge at high thread counts if the per-variant work is constant. The absolute throughput number (e.g., "100K variants/s at 16 threads") is the key metric a user cares about.

---

### Experiment 4: AFQuery vs. bcftools (`05_vs_bcftools.py`)

**Purpose:** Compare AFQuery's performance against bcftools, the de facto standard for VCF manipulation. This establishes whether the bitmap approach offers a meaningful speedup over the traditional "recompute from genotypes" approach.

**Setup:** 1KG databases at 3 subset sizes (500, 1K, 2,504 samples) are compared. For each subset, a multi-sample VCF is created with `bcftools view -S`. Three operations are timed for both tools:

| Operation | bcftools equivalent | AFQuery equivalent |
|-----------|--------------------|--------------------|
| Point query | `bcftools view -r CHR:POS \| +fill-tags` | `db.query(chrom, pos)` |
| Subset query | `bcftools view -r CHR:POS -S females.txt \| +fill-tags` | `db.query(chrom, pos, sex="female")` |
| Full dump | `bcftools +fill-tags \| query -f ...` | `db.dump(chrom=chrom)` |

Both tools are timed cold (AFQuery includes `Database()` init; bcftools includes subprocess spawn). Each operation is repeated `BCFTOOLS_REPS` times (default 10).

Additionally, AF concordance is computed on the largest subset: AFQuery and bcftools AF values for all chr22 variants are compared and an R^2 is calculated.

**Figure 4 (`fig4_bcftools_comparison`):** Grouped bar chart (log-scale Y). Three operation groups, each with a bcftools bar and an AFQuery bar. Value labels show absolute time in ms.

**How to interpret:** AFQuery should be faster for point and subset queries because it reads pre-computed bitmaps (O(1) per variant) rather than scanning all genotypes (O(N) per variant). The speedup ratio is the key claim. For full dump, the gap may narrow because both tools must scan all data. If bcftools is faster for dump, that is expected (streaming VCF is hard to beat for sequential reads).

**Figure 5 (`fig5_concordance`):** Scatter plot. X-axis: bcftools AF. Y-axis: AFQuery AF. A text box shows R^2 and number of common variants.

**How to interpret:** Points should fall tightly along the y=x diagonal. R^2 should be very close to 1.0 (ideally > 0.9999). Deviations indicate either a bug in AF computation or differences in multi-allelic handling. This plot validates correctness, not performance. A small number of `afquery_only` or `bcftools_only` variants is expected due to FILTER handling (AFQuery stores PASS-only by default).

---

## Timing Methodology

**Query scaling (Experiment 1):** In-process timing via `time.perf_counter()`. Cold queries reinstantiate the `Database` object per iteration. Warm queries reuse a pre-loaded `Database` after discarding warmup iterations. Reports median and IQR over 50 warm / 5 cold repetitions.

**Build performance (Experiment 2):** Timed via `/usr/bin/time -v` (wall-clock + peak RSS). Each (n_samples, threads) combination is repeated `BUILD_REPS` times; the collect script reports median and IQR across repetitions.

**Annotation throughput (Experiment 3):** In-process timing. Each (n_variants, threads) combination is repeated `ANNOTATE_REPS` times; median throughput reported.

**AFQuery vs. bcftools (Experiment 4):** Both tools are timed cold -- AFQuery timing includes `Database()` initialization, bcftools timing includes subprocess spawn and VCF parsing. This reflects first-query latency for each tool.

## Data

Not included in the repository. Downloaded and generated automatically by the pipeline:

- **1000 Genomes Phase 3 chr22:** Downloaded from EBI FTP (~200 MB compressed)
- **Synthetic data:** Generated via `afquery.preprocess.synth` at scales 1K-50K

## Output

- `perf_results/*.json` -- aggregated timing data
- `perf_results/raw/` -- per-run JSON files
- `figures/*.{pdf,png}` -- high-resolution figures

## Hardware

Document the hardware used when reporting results:

| Parameter | Value |
|-----------|-------|
| CPU | (model, cores, frequency) |
| RAM | (total) |
| Storage | (type: SSD/HDD/Lustre, model) |
| OS | (distribution, kernel version) |

## References

- Weber LM et al. (2019). Essential guidelines for computational method benchmarking. *Genome Biology* 20:125. DOI: 10.1186/s13059-019-1738-8
- Auton A et al. (2015). A global reference for human genetic variation. *Nature* 526:68-74. DOI: 10.1038/nature15393
- 1000 Genomes Phase 3 FTP: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
