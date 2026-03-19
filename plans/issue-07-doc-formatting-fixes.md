# Issue 07 — Documentation Formatting Fixes

## Summary

Multiple formatting issues have been identified across the documentation: illegible Mermaid diagrams (text rendered on one line or at tiny font size), lists rendered on a single line instead of bullet points, and a CLI inconsistency (`--db-dir` vs `--db` in benchmark).

## Source notes

| File | Line | Issue |
|------|------|-------|
| `docs/advanced/performance.md` | 95 | Mermaid "Query Execution Path" diagram renders on one line — illegible. All Mermaid diagrams should be audited. |
| `docs/advanced/performance.md` | 181 | Two profiling examples are identical; one should be removed or differentiated. |
| `docs/guides/update-database.md` | 26 | Mermaid timeline diagram has very small text — hard to read. |
| `docs/guides/sample-filtering.md` | 43 | Mermaid flow diagram has very small font size. |
| `docs/getting-started/quickstart.md` | 28 | Field list rendered on a single line; should be a bullet list. |
| `docs/use-cases/clinical-prioritization.md` | 43 | Variable list rendered on a single line; should be a bullet list. |
| `docs/use-cases/clinical-prioritization.md` | 120 | AN threshold list rendered on a single line; should be a bullet list. |
| `docs/advanced/ploidy-and-sex-chroms.md` | 86 | List rendered on a single line; should be a bullet list. |
| `docs/reference/cli.md` | 208 | `benchmark` module uses `--db-dir` instead of `--db`; fix in both code and docs. |

## Scope

Docs fixes + one small CLI/code change (`--db-dir` → `--db` in benchmark).

## Steps

### 1. Fix Mermaid diagrams

For each illegible diagram, investigate the root cause:
- Diagrams that render on one line usually have all nodes defined without newlines inside the ```` ```mermaid ```` block, or use a graph direction that forces horizontal layout.
- Fix by rewriting the diagram with proper newlines, splitting long labels, or switching direction (e.g. `LR` → `TD`).

Files to fix:
- `docs/advanced/performance.md` line ~95 — "Query Execution Path"
- `docs/guides/update-database.md` line ~26 — timeline diagram
- `docs/guides/sample-filtering.md` line ~43 — filtering flow diagram

After fixing, run `mkdocs serve` and visually inspect each diagram.

### 2. Remove duplicate profiling example

In `docs/advanced/performance.md` around line 181, two profiling code blocks are identical. Either:
- Remove one entirely, or
- Make the second one a distinct example (e.g. profiling a region query vs. a point query).

### 3. Fix single-line lists

The following files contain lists formatted as a single comma-separated or space-separated line that should instead use Markdown bullet syntax (`- item`):

- `docs/getting-started/quickstart.md` line ~28
- `docs/use-cases/clinical-prioritization.md` lines ~43 and ~120
- `docs/advanced/ploidy-and-sex-chroms.md` line ~86

Read each location, identify the broken list, and reformat with proper `- ` prefixes and newlines.

### 4. Fix `--db-dir` → `--db` in benchmark module

**Code change** (`src/afquery/cli.py` or benchmark entry point):
- Find the `benchmark` command's option definition; rename `--db-dir` to `--db`.
- Update internal usage of `db_dir` parameter name if needed.

**Docs change** (`docs/reference/cli.md` line ~206):
- Update the option table: `--db-dir` → `--db`.

### 5. Verification

```bash
mkdocs build  # no warnings about broken links
pytest --tb=short -q  # benchmark tests still pass
afquery benchmark --help  # confirm --db is shown, not --db-dir
```
