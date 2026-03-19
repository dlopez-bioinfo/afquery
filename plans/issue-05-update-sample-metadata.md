# Issue 05 — Update Sample Metadata (new `update-db --update-metadata` subcommand)

## Summary

Sample phenotype codes and other metadata are currently immutable after ingestion. A new feature is needed so that an operator can update a sample's metadata (sex, phenotype codes, technology) without re-ingesting the VCF. Every such change must be logged in the database changelog for reproducibility.

## Source notes

| File | Line | Note |
|------|------|------|
| `docs/getting-started/concepts.md` | 223 | "Codes are immutable after ingestion" — a new feature to update them is planned; changelog entry required. |
| `docs/faq.md` | 61 | FAQ about metadata correction currently says it is not supported; once implemented, remove this FAQ. |
| `docs/reference/cli.md` | 107 | Documents the desired `update-db` extension; design must be confirmed before implementation. |

## Scope

New feature, high complexity — touches SQLite schema, precomputed bitmaps, CLI, models, and documentation.

> **Design note**: Before coding, clarify with the user the exact fields that can be updated, the single-sample vs. batch interface, and what happens to precomputed bitmaps when phenotypes change.

## Steps

### 1. Design the interface (confirm with user before coding)

Proposed CLI:

```bash
# Single sample
afquery update-db --db ./db/ --update-sample SAMPLE_ID --set-phenotype "E11.9,I10" --set-sex female

# Batch (TSV: sample_id  field  new_value)
afquery update-db --db ./db/ --update-samples-file changes.tsv
```

Fields that may be updated: `sex`, `phenotype_codes`, `technology` (tech change requires BED reassignment — confirm scope).

### 2. SQLite schema changes (`metadata.sqlite`)

- Add `sample_updates` table or reuse existing `changelog` table with a new operation type `UPDATE_SAMPLE`.
- Columns: `timestamp`, `sample_id`, `field`, `old_value`, `new_value`, `operator_note` (optional).

### 3. Core logic (`src/afquery/preprocess/__init__.py` or new `src/afquery/update.py`)

- Fetch the current sample row from SQLite.
- Apply the field change.
- **Recompute affected precomputed bitmaps** in `metadata.sqlite`:
  - Sex bitmaps (`male_bitmap`, `female_bitmap`).
  - Phenotype bitmaps (remove sample bit from old codes; add to new codes).
  - Technology bitmaps if technology changes.
- Write a changelog entry.

### 4. CLI (`src/afquery/cli.py` or equivalent)

- Extend `update-db` command with `--update-sample` / `--update-samples-file` options.
- Add input validation: unknown sample IDs, invalid field names, malformed phenotype codes.
- Print a summary of changes and changelog entry IDs on success.

### 5. Update documentation

- `docs/getting-started/concepts.md` — replace "Codes are immutable…" with the new behaviour description.
- `docs/faq.md` — remove the FAQ entry that says metadata cannot be changed; replace with a pointer to the guide.
- `docs/reference/cli.md` — document new `--update-sample` / `--update-samples-file` options under `update-db`.
- Add a new guide page or section: `docs/guides/update-database.md` (extend existing page).

### 6. Update tests

- Unit tests for bitmap recomputation after phenotype change.
- Integration test: update a sample's phenotype, then query with the new phenotype; verify AC/AN change.
- Test for batch update from TSV.
- Test for changelog entry creation.

### 7. Verification

```bash
pytest --tb=short -q
afquery update-db --db ./db/ --update-sample S1 --set-phenotype "E11.9"
afquery info --db ./db/  # confirm changelog shows the update
afquery query --db ./db/ --locus chr1:925952 --phenotype E11.9  # verify AF reflects the change
```
