# Issue 06 — Multi-Chromosome Batch Query

## Summary

The current `--from-file` batch mode requires all variants to be on the same chromosome (set via `--chrom`). A new feature is needed to allow a single batch query across multiple chromosomes. This is closely related to issue 03 (unified chrom:pos syntax) and should be implemented as part of that work.

## Source notes

| File | Line | Note |
|------|------|------|
| `docs/faq.md` | 11 | Multi-chromosome queries are requested; once implemented, remove this FAQ entry. |

## Scope

New feature — depends on Issue 03 (chrom:pos syntax). If Issue 03 is implemented first, this issue may be largely resolved; verify before doing separate work.

## Dependency

This issue is tightly coupled to **Issue 03**. The batch TSV format proposed in Issue 03 already includes a `chrom` column, which is the prerequisite for multi-chromosome queries. Implement Issue 03 first, then verify whether this issue is fully resolved.

## Steps

### 1. Verify overlap with Issue 03

After Issue 03 is merged:
- Run `afquery query --from-file multi_chrom.tsv` with a file containing variants on multiple chromosomes.
- If it works, close this issue and remove the FAQ entry.

### 2. If additional work is needed: group by chromosome in batch execution

In `src/afquery/query.py` (batch query path):
- Group input rows by `chrom`.
- For each chromosome group, load the corresponding Parquet file once.
- Collect results and return in original input order.

### 3. Performance consideration

- Opening multiple Parquet files sequentially is acceptable for typical batch sizes.
- If the batch spans many chromosomes, consider pre-sorting to minimise repeated file opens.

### 4. Update documentation

- `docs/faq.md` — remove the FAQ entry "Can I query multiple chromosomes at once?" and replace with a pointer to the batch query guide.
- `docs/guides/query.md` — add a note in the Batch Query section that the file may contain variants from different chromosomes.

### 5. Update tests

- Add a test that runs a batch query spanning at least two chromosomes and verifies correct results for all rows.

### 6. Verification

```bash
pytest --tb=short -q
afquery query --db ./db/ --from-file multi_chrom.tsv --format tsv
```
