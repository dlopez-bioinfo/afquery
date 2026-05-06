# Create a Database

`afquery create-db` builds the AFQuery database from a manifest of single-sample VCFs. This is a one-time setup step; incremental updates use `afquery update-db`.

---

## Basic Usage

```bash
afquery create-db \
  --manifest manifest.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38
```

For cohorts with WES/panel samples, provide the BED file directory:

```bash
afquery create-db \
  --manifest manifest.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --bed-dir ./beds/
```

---

## What Happens

1. **Ingest phase** — Each VCF is parsed with cyvcf2. Genotypes and INFO fields are written to a SQLite temporary database, one row per variant per sample.
2. **Build phase** — DuckDB reads the SQLite data, aggregates per 1-Mbp bucket, and writes Roaring Bitmap Parquet files partitioned by chromosome/bucket.
3. **Finalize** — `manifest.json` and `metadata.sqlite` are written to the output directory.

---

## Directory Layout After Creation

```
./db/
├── manifest.json          # Build configuration (genome build, schema version, etc.)
├── metadata.sqlite        # Sample/phenotype/technology/changelog metadata
├── variants/              # Hive-partitioned Parquet files
│   ├── chr1/
│   │   ├── bucket_0/      # Positions 0–999,999
│   │   │   └── data.parquet
│   │   ├── bucket_1/      # Positions 1,000,000–1,999,999
│   │   │   └── data.parquet
│   │   └── ...
│   ├── chr2/
│   └── ...
└── capture/               # Interval trees for WES technologies (pickle files)
    ├── wes_v1.pkl
    └── wes_v2.pkl
```

---

## Memory and Thread Tuning

For large cohorts, tune these options:

```bash
afquery create-db \
  --manifest manifest.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --build-threads 32 \
  --build-memory 4GB
```

| Option | Default | Recommendation |
|--------|---------|----------------|
| `--build-threads` | all CPUs | Set to `min(cpu_count, available_RAM_GB / 2)` |
| `--build-memory` | `2GB` | Increase for dense WGS regions or large cohorts |
| `--threads` | all CPUs | Controls ingest parallelism (VCF parsing) |

The build phase uses one DuckDB process per 1-Mbp bucket. With `--build-threads 32` and `--build-memory 4GB`, peak RAM usage is approximately `32 × 4 = 128 GB`.

---

## Resume Behavior

If `create-db` is interrupted, it resumes automatically from where it left off. Individual bucket Parquet files that were already written are skipped.

To force a complete restart:

```bash
afquery create-db --manifest manifest.tsv --output-dir ./db/ --genome-build GRCh38 --force
```

!!! warning
    `--force` deletes all existing output in `--output-dir`. Use with caution.

---

## FILTER=PASS Behavior

By default, only variants with `FILTER=PASS` (or no FILTER field) are counted in AC/AN. Variants that fail filters are tracked in `fail_bitmap`. PASS-only ingestion is always enforced — there is currently no CLI option to change this behaviour.

See [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md) for details.

---

## Coverage-Evidence Filters (Phase 2)

Four optional flags enable per-sample, quality-aware tracking of which positions
each WES technology actually covered. They are fully opt-in: omit them and the
database behaves exactly as before.

| Flag | Default | Effect |
|------|---------|--------|
| `--min-dp D`     | 0   | Minimum `FORMAT/DP` for a carrier to count as quality evidence. |
| `--min-gq G`     | 0   | Minimum `FORMAT/GQ` for a carrier to count as quality evidence. |
| `--min-qual Q`   | 0.0 | Minimum VCF `QUAL` field for a carrier to count as quality evidence. |
| `--min-covered K`| 0   | Per WES tech, position is "trusted" only if at least K of its carriers pass the quality thresholds. Triggers Phase 2 storage when ≥1. |

When any of these flags is non-zero AFQuery:

1. Reads `FORMAT/DP`, `FORMAT/GQ`, and `QUAL` from each variant call during ingest.
   Use the bundled `resources/normalize_vcf.sh` (which preserves these FORMAT
   fields) or ensure your own preprocessing keeps them.
2. Stores two new columns in each variant Parquet row:
   - `quality_pass_bitmap` — carriers meeting all thresholds.
   - `filtered_bitmap` — non-carrier WES samples whose tech failed the
     `--min-covered` gate. The query layer reports these in `N_NO_COVERAGE`
     instead of `N_HOM_REF`.
3. Bumps `schema_version` from `2.0` to `3.0` and persists the chosen thresholds
   under `coverage_filter` in `manifest.json`.

Example:

```bash
afquery create-db \
  --manifest samples.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --bed-dir ./beds/ \
  --min-dp 30 --min-gq 20 --min-covered 1
```

Phase 2 thresholds are fixed at creation time. `update-db --add-samples` reuses
them and recomputes `quality_pass_bitmap` and `filtered_bitmap` for every
position whose WES tech receives new samples (see
[Update Database](update-database.md)).

See [Coverage Evidence](../advanced/coverage-evidence.md) for the full data model
and the query-time companion flag `--min-quality-evidence`.

---

## Validating the Result

After creation, run:

```bash
afquery check --db ./db/
```

And inspect database metadata:

```bash
afquery info --db ./db/
```

Example `info` output:
```
Database:       ./db/
Schema version: 2.0
Genome build:   GRCh38
DB version:     1.0
Samples:        1371
Technologies:   wgs, wes_v1, wes_v2
Chromosomes:    chr1 ... chrX chrY chrM
```

These commands are available at any time after database creation — not only immediately after `create-db`. Use `afquery check` to verify database integrity (manifest consistency, Parquet file health, capture index presence) and `afquery info` to inspect sample counts, registered technologies, and phenotype codes.

---

## Full Option Reference

See [CLI Reference → create-db](../reference/cli.md#create-db).

---

## Next Steps

- [Manifest Format](manifest-format.md) — TSV manifest column reference and common mistakes
- [Query Allele Frequencies](query.md) — run your first queries against the new database
- [Performance Tuning](../advanced/performance.md) — build thread and memory configuration for large cohorts
