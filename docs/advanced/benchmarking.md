# Benchmarking

`afquery benchmark` measures query performance on synthetic or real data and produces a JSON report.

---

## Quick Benchmark

Run against synthetic data (no real database required):

```bash
afquery benchmark
```

This generates 1,000 synthetic samples and 10,000 variants per chromosome, builds an in-memory database, runs a suite of query types, and writes results to `benchmark_report.json`.

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--n-samples` | `1000` | Number of synthetic samples |
| `--n-variants` | `10000` | Number of variants per chromosome |
| `--output` | `benchmark_report.json` | Output path for JSON report |
| `--db-dir` | None | Run against an existing database instead of synthetic data |

---

## Benchmark Against a Real Database

```bash
afquery benchmark --db-dir ./my_db/ --output my_db_benchmark.json
```

This uses your actual Parquet files and sample metadata, giving realistic performance numbers for your specific cohort size and variant density.

---

## Report Format

The JSON report contains timing results for each query type:

```json
{
  "n_samples": 1000,
  "n_variants": 10000,
  "genome_build": "GRCh38",
  "queries": {
    "point_query_cold_ms": 87.3,
    "point_query_warm_ms": 9.1,
    "region_query_1mbp_ms": 312.4,
    "batch_query_100_ms": 198.7,
    "annotation_100_variants_ms": 523.1
  },
  "timestamp": "2026-03-16T10:00:00Z"
}
```

---

## Interpreting Results

| Metric | Target | Notes |
|--------|--------|-------|
| `point_query_cold_ms` | < 100 ms | First query after cold start; includes Parquet I/O |
| `point_query_warm_ms` | < 20 ms | Subsequent queries; OS page cache active |
| `region_query_1mbp_ms` | < 500 ms | Depends on variant density in region |
| `batch_query_100_ms` | < 300 ms | 100 pre-specified variants on same chromosome |
| `annotation_100_variants_ms` | < 2000 ms | VCF annotation, single-threaded |

If `point_query_cold_ms` exceeds 500 ms, check disk I/O performance. If `point_query_warm_ms` is slow, check available RAM for OS page cache.

---

## Synthetic Data Generation

The benchmark generates:
- `N` samples with random sex (50/50) and random technology assignment
- `M` variants per chromosome with uniformly random positions
- Random genotypes with configurable carrier rates

Synthetic data is written to a temporary directory and cleaned up after the benchmark completes.

---

## Tracking Performance Over Time

Run the benchmark after major updates to detect regressions:

```bash
# Before update
afquery benchmark --db-dir ./db/ --output before.json

# After adding 500 samples
afquery update-db --db ./db/ --add-samples new_batch.tsv
afquery benchmark --db-dir ./db/ --output after.json

# Compare
python3 -c "
import json
before = json.load(open('before.json'))['queries']
after  = json.load(open('after.json'))['queries']
for k in before:
    delta = after[k] - before[k]
    print(f'{k}: {before[k]:.1f} → {after[k]:.1f} ms  ({delta:+.1f})')
"
```

---

## Comparison with Alternative Approaches

For context on AFQuery's design trade-offs, see [Comparison Table](../publication/comparison-table.md).

Brief summary of query latency comparison:

| Tool | 50K samples, chr1, point query |
|------|-------------------------------|
| AFQuery (cold) | <100 ms |
| AFQuery (warm) | ~10 ms |
| bcftools stats (VCF scan) | ~5 minutes |
| VCFtools --freq (VCF scan) | ~5 minutes |

AFQuery achieves this speed advantage because queries access only the bitmaps for a single 1-Mbp bucket (one Parquet file) rather than scanning the entire dataset.

---

## Regression Testing

Run the benchmark after each major update to detect performance regressions:

```bash
# Save baseline before changes
afquery benchmark --db-dir ./db/ --output baseline.json

# After update
afquery benchmark --db-dir ./db/ --output after_update.json

# Compare
python3 - <<'EOF'
import json
b = json.load(open('baseline.json'))['queries']
a = json.load(open('after_update.json'))['queries']
for k in b:
    pct = (a[k] - b[k]) / b[k] * 100
    flag = " ⚠" if pct > 20 else ""
    print(f"{k}: {b[k]:.1f} → {a[k]:.1f} ms  ({pct:+.1f}%){flag}")
EOF
```

A regression of >20% on `point_query_cold_ms` warrants investigation.
