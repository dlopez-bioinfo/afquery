# CLI Reference

All AFQuery commands follow the pattern `afquery <command> [OPTIONS]`.

---

## create-db

Build a new AFQuery database from a manifest of single-sample VCFs.

```
afquery create-db [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--manifest` | TEXT | **required** | Path to TSV manifest file |
| `--output-dir` | TEXT | **required** | Path to output database directory |
| `--genome-build` | `GRCh37`\|`GRCh38` | **required** | Reference genome build |
| `--threads` | INTEGER | all CPUs | Worker threads for the ingest phase (VCF parsing) |
| `--build-threads` | INTEGER | min(cpu_count, n_buckets) | Max parallel workers for the build phase (DuckDB) |
| `--build-memory` | TEXT | `2GB` | DuckDB memory limit per build worker |
| `--tmp-dir` | TEXT | `{output_dir}/.tmp_preprocess` | Temporary directory for intermediate files |
| `--bed-dir` | TEXT | None | Directory containing BED files for WES technologies |
| `--force` | flag | False | Delete any partial results and restart from scratch |
| `--db-version` | TEXT | `1.0` | Version label for this database |
| `-v, --verbose` | flag | False | Verbose output with per-item progress |

---

## query

Query allele frequencies at one or more positions.

```
afquery query [OPTIONS]
```

Exactly one of `--locus`, `--region`, or `--from-file` must be provided.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--locus` | TEXT | None | Single position as `CHROM:POS` (e.g., `chr1:925952`) |
| `--region` | TEXT | None | Genomic range as `CHROM:START-END` (e.g., `chr1:900000-1000000`) |
| `--from-file` | PATH | None | Headerless TSV with columns `chrom pos [ref [alt]]` (batch query; multi-chromosome supported) |
| `--phenotype` | TEXT | None | Phenotype filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--sex` | `male`\|`female`\|`both` | `both` | Restrict to specified sex |
| `--tech` | TEXT | None | Technology filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--ref` | TEXT | None | Filter to specific reference allele (only for `--locus`) |
| `--alt` | TEXT | None | Filter to specific alternate allele (only for `--locus`) |
| `--format` | `text`\|`json`\|`tsv` | `text` | Output format |

---

## annotate

Annotate a VCF file with allele frequency information.

```
afquery annotate [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--input` | TEXT | **required** | Input VCF file (plain or `.gz`) |
| `--output` | TEXT | **required** | Output annotated VCF file |
| `--phenotype` | TEXT | None | Phenotype filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--sex` | `male`\|`female`\|`both` | `both` | Restrict to specified sex |
| `--tech` | TEXT | None | Technology filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--threads` | INTEGER | all CPUs | Number of worker threads for parallel annotation |
| `-v, --verbose` | flag | False | Verbose output with per-item progress |

---

## dump

Export allele frequency data to CSV.

```
afquery dump [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `-o, --output` | TEXT | stdout | Output CSV file path |
| `--chrom` | TEXT | None | Restrict export to this chromosome |
| `--start` | INTEGER | None | 1-based start position (inclusive). Requires `--chrom`. |
| `--end` | INTEGER | None | 1-based end position (inclusive). Requires `--chrom`. |
| `--phenotype` | TEXT | None | Phenotype filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--sex` | `male`\|`female`\|`both` | `both` | Restrict to specified sex |
| `--tech` | TEXT | None | Technology filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--by-sex` | flag | False | Disaggregate output by sex (adds `AC_male`/`AC_female` columns) |
| `--by-tech` | flag | False | Disaggregate output by technology (adds `AC_<tech>` columns) |
| `--by-phenotype` | TEXT | None | Disaggregate by specific phenotype codes. Repeatable. |
| `--all-groups` | flag | False | Disaggregate by all sexes × technologies × phenotypes (Cartesian product) |
| `--threads` | INTEGER | all CPUs | Number of worker threads for parallel export |
| `-v, --verbose` | flag | False | Verbose output with per-item progress |

---

## update-db

Add samples, remove samples, update sample metadata, or compact the database.

```
afquery update-db [OPTIONS]
```

At least one of `--remove-samples`, `--add-samples`, `--compact`, `--update-sample`, or `--update-samples-file` must be provided. Operations execute in order: remove → update-metadata → add → compact.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--remove-samples` | TEXT | None | Sample name(s) to remove. Repeatable; comma-separated or multiple flags. |
| `--add-samples` | PATH | None | Manifest TSV of new samples to add. Repeatable for multiple manifests. |
| `--compact` | flag | False | Remove dead bits from removed samples to reclaim disk space |
| `--update-sample` | TEXT | None | Sample name to update (single-sample metadata mode). Requires `--set-sex` and/or `--set-phenotype`. |
| `--set-sex` | TEXT | None | New sex for `--update-sample`. Options: `male`, `female`. |
| `--set-phenotype` | TEXT | None | New phenotype codes (comma-separated) for `--update-sample`. Replaces all current codes. |
| `--update-samples-file` | PATH | None | TSV file for batch metadata update. Header: `sample_name`, `field`, `new_value`. Mutually exclusive with `--update-sample`. |
| `--operator-note` | TEXT | None | Free-text note appended to each changelog entry for this metadata update. |
| `--threads` | INTEGER | all CPUs | Number of worker threads for parallel processing |
| `--tmp-dir` | TEXT | system temp | Temporary directory for intermediate files |
| `--bed-dir` | TEXT | None | Directory containing BED files for WES technologies |
| `--db-version` | TEXT | auto-increment | New version label after update |
| `-v, --verbose` | flag | False | Verbose output with per-item progress |

---

## info

Display database metadata and statistics.

```
afquery info [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--samples` | flag | False | List all samples with metadata |
| `--changelog` | flag | False | Show full changelog history |
| `--format` | `table`\|`tsv`\|`json` | `table` | Output format |

---

## check

Validate database integrity.

```
afquery check [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |

Exits with code 0 if the database is healthy, non-zero otherwise.

---

## version show

Display the database version label.

```
afquery version show [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |

---

## version set

Set a new version label on the database.

```
afquery version set [OPTIONS] NEW_VERSION
```

| Argument/Option | Type | Description |
|-----------------|------|-------------|
| `NEW_VERSION` | TEXT (positional) | New version label |
| `--db` | TEXT | **required** — Path to database directory |

---

## benchmark

Run performance benchmarks on synthetic or real data.

```
afquery benchmark [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--n-samples` | INTEGER | `1000` | Number of synthetic samples to generate |
| `--n-variants` | INTEGER | `10000` | Number of variants per chromosome |
| `--output` | TEXT | `benchmark_report.json` | Output path for JSON benchmark report |
| `--db-dir` | TEXT | None | Use an existing database instead of generating synthetic data |

NOTA_PARA_LLM: todos los módulos usan --db para proporcionar la ruta de la base de datos excepto el módulo benchmark. Se debe cambiar por --db para mantener la coherencia con el resto de la herramienta.

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (invalid arguments, file not found, database integrity failure, etc.) |
| 2 | Usage error (Click framework: missing required argument or unknown option) |
