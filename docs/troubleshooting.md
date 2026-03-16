# Troubleshooting

---

## create-db Runs Out of Memory

**Symptom:** `create-db` fails with `MemoryError`, `Out of memory`, or a DuckDB `Not enough memory to store` error.

**Cause:** The per-worker DuckDB memory limit is too low for your cohort size or variant density.

**Fix:**
- Lower the number of parallel build workers: `--build-threads 4`
- Increase per-worker memory: `--build-memory 4GB`
- Both together keep total RAM the same with fewer concurrent workers

```bash
# Instead of default (all CPUs × 2GB)
afquery create-db ... --build-threads 8 --build-memory 4GB
```

See [Performance Tuning](advanced/performance.md) for sizing guidance.

---

## Query Returns AN=0

**Symptom:** A query returns `AC=0, AN=0, AF=None` for a variant you know exists.

**Possible causes:**

1. **All samples excluded by filters** — Check your `--phenotype`, `--sex`, or `--tech` filters. AN=0 means no eligible samples.
2. **WES position outside capture regions** — If all samples use WES and the position is not in any BED file, no samples are eligible.
3. **Chromosome not ingested** — Check `afquery info --db ./db/` to confirm the chromosome is present.
4. **Variant not in database** — The variant may not be present in any sample VCF. Use `--region` to check nearby positions.

---

## VCF Annotation Is Slow

**Symptom:** `afquery annotate` takes much longer than expected.

**Fixes:**
- Increase thread count: `--threads 16`
- Check disk I/O: annotation reads many small Parquet files; SSDs are significantly faster than spinning disks
- For very large VCFs (1M+ variants), annotation time scales linearly with variant count

---

## DuckDB Errors During Build

**Symptom:** Errors like `"Not supported: Writing to Arrow IPC"` or `"Unsupported file format"` during `create-db`.

**Cause:** An older DuckDB version is attempting to write Arrow IPC format for temporary files instead of Parquet.

**Fix:** Upgrade to DuckDB ≥ 0.10:

```bash
pip install --upgrade duckdb
```

AFQuery requires DuckDB to use Parquet for all temporary files. Arrow IPC is not supported.

---

## cyvcf2 ImportError in Workers

**Symptom:** `ImportError: cannot import name 'VCF' from 'cyvcf2'` in worker processes during ingest.

**Cause:** cyvcf2 cannot be pickled for `ProcessPoolExecutor` — it must be imported inside the worker function body.

**Fix:** This is an internal invariant of AFQuery. If you see this error with the latest version, please file a bug report. Do not import cyvcf2 at module level in any code that runs in worker processes.

---

## Wrong Bucket IDs (Silent Bug)

**Symptom:** Queries return no results for known variants, but the database was built without errors.

**Cause:** A DuckDB float-division rounding bug was present in older AFQuery versions. When computing bucket IDs, `CAST(pos / 1000000 AS BIGINT)` rounds floats incorrectly (e.g., position 1,500,000 → bucket 2 instead of 1).

**Fix:** Upgrade to the latest AFQuery version. The fix uses `CAST(pos AS BIGINT) / 1000000` for correct integer division. Rebuild the database after upgrading.

---

## Compact Takes a Long Time

**Symptom:** `afquery update-db --compact` runs for many minutes or hours.

**Cause:** Compaction rewrites every Parquet file in the database. For large databases (many chromosomes × many buckets), this is expected.

**Recommendation:** Run compact during off-hours or overnight. It is safe to interrupt (resume is not supported; re-run to complete).

---

## Sample Not Found in Remove Operation

**Symptom:** `afquery update-db --remove-samples SAMP_001` fails with `"Sample not found"`.

**Cause:** The sample name is case-sensitive and must match exactly what was in the original manifest.

**Fix:** Check the exact sample name:

```bash
afquery info --db ./db/ --samples | grep SAMP
```

---

## afquery check Reports Errors

**Symptom:** `afquery check --db ./db/` exits non-zero and prints error messages.

**Common errors:**

| Error | Fix |
|-------|-----|
| `Missing Parquet for chromosome chr3` | Re-run `create-db` or investigate incomplete build |
| `Manifest mismatch: expected N samples, found M` | Database may be partially updated; re-run `update-db` |
| `Capture file missing for wes_v1` | BED file was not provided at build time; rebuild with `--bed-dir` |

---

## Version Compatibility

AFQuery databases are forward-compatible but not backward-compatible. A database built with AFQuery 1.x can be read by AFQuery 2.x, but not vice versa.

If you get `"Unsupported schema_version"` errors, upgrade AFQuery:

```bash
pip install --upgrade afquery
```
