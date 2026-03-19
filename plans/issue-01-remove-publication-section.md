# Issue 01 — Remove Publication Section

## Summary

The `docs/publication/` directory contains three pages (benchmarks-methodology, scientific-overview, comparison-table) that should not exist as a standalone section in the user-facing documentation. The content should either be moved into the appropriate existing sections or deleted if it is redundant.

## Source notes

| File | Line | Note |
|------|------|------|
| `docs/publication/benchmarks-methodology.md` | 1 | No publication section or benchmarking methodology page should appear in the docs. |
| `docs/publication/scientific-overview.md` | 1 | No publication section; content not yet moved should be reviewed, merged into other sections or deleted if redundant. |
| `docs/publication/comparison-table.md` | 1 | No comparison-table page; some content already moved, the rest should be merged or deleted. |

## Scope

- Docs only (no code changes required).

## Steps

1. **Review `docs/publication/scientific-overview.md`** — identify paragraphs already present elsewhere in the docs and delete them; for any paragraph not yet covered, decide the best target section (e.g. `docs/index.md`, `docs/getting-started/concepts.md`) and migrate it.
2. **Review `docs/publication/comparison-table.md`** — same process: detect duplicates and drop them; move non-redundant content (e.g. comparison with gnomAD) to a relevant section such as `docs/index.md` or a new `docs/advanced/` page.
3. **Review `docs/publication/benchmarks-methodology.md`** — move any benchmark information that belongs in `docs/advanced/performance.md` or similar; delete the rest.
4. **Delete the three source files** and the `docs/publication/` directory.
5. **Update `mkdocs.yml`** — remove the `publication` nav section.
6. **Check for internal links** pointing to `publication/` anywhere in the docs and redirect or remove them.
7. Run `mkdocs build` to verify no broken links remain.
