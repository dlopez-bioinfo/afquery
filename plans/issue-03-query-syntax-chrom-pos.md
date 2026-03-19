# Issue 03 — Unified chrom:pos Query Syntax

## Summary

The current query interface requires `--chrom` as a separate flag alongside `--pos`, `--region`, and `--from-file`. The desired UX is a single positional-style argument in `chrom:pos` / `chrom:start-end` / `chrom\tpos\tref\talt` format so that `--chrom` is eliminated and multi-chromosome batch queries become possible.

## Source notes

| File | Line | Note |
|------|------|------|
| `docs/guides/query.md` | 15 | Point query: new format `chr1:925952`; remove `--chrom`. |
| `docs/guides/query.md` | 33 | Region query: new format `chr1:900000-1000000`; remove `--chrom`. |
| `docs/guides/query.md` | 47 | Batch query: TSV columns `chrom\tpos\tref\talt` (ref/alt optional); remove `--chrom`; allow multi-chromosome batches. |
| `docs/guides/query.md` | 70 | Output docs missing `N_FAIL`; update once code changes are in. |

## Scope

Code + Docs (high complexity — CLI, `QueryParams`, query engine, annotate, tests, all documentation).

## Steps

### 1. Update `models.py`

- Modify `QueryParams` to accept a `locus` string (`"chr1:925952"`) instead of separate `chrom` / `pos` fields, **or** keep the internal representation but parse `locus` at the CLI boundary.
- For region queries, parse `"chr1:900000-1000000"` into `(chrom, start, end)`.
- For batch queries, each row in the TSV carries its own `chrom` column.

### 2. Update CLI (`src/afquery/cli.py` or equivalent)

- Point query: replace `--chrom + --pos` with a single `--locus chr1:925952` option (or positional argument).
- Region query: replace `--chrom + --region` with `--region chr1:900000-1000000`.
- Batch query: replace `--chrom + --from-file` with `--from-file variants.tsv`; the file now includes a `chrom` column.
- Remove `--chrom` entirely. Add validation that raises a clear error for malformed locus strings.
- Update `--help` text.

### 3. Update `QueryEngine` / `query.py`

- Adjust any internal call sites that pass `chrom` separately.
- Batch execution loop must now group by chromosome (for Parquet file lookup) while preserving input order in the output.

### 4. Update `annotate.py`

- If `annotate` internally calls the query engine with explicit chrom/pos, update those call sites to use the new interface.

### 5. Update TSV batch format

Old format (headerless, 3 cols: `pos ref alt`):
```
925952	G	A
```
New format (4 cols: `chrom pos ref alt`, ref/alt optional):
```
chr1	925952	G	A
chr2	1014541	C	T
```

### 6. Update output: add `N_FAIL`

- Verify `N_FAIL` is present in the `QueryResult` model.
- Add `N_FAIL` to text, TSV, and JSON output formatters.
- Update `docs/guides/query.md` output examples.

### 7. Update all documentation

- `docs/guides/query.md` — rewrite Point, Region, and Batch sections; update output examples.
- `docs/reference/cli.md` — update `query` command option table.
- `CLAUDE.md` / README examples if present.
- Any other doc file with `--chrom` examples.

### 8. Update tests

- `tests/test_query.py` (and similar) — update all calls that pass `chrom` separately.
- Add tests for multi-chromosome batch queries.
- Add tests for malformed locus strings.

### 9. Verification

```bash
pytest --tb=short -q
afquery query --db ./db/ --locus chr1:925952
afquery query --db ./db/ --region chr1:900000-1000000
afquery query --db ./db/ --from-file multi_chrom.tsv
```
