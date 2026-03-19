# Issue 02 — Drop Schema v1

## Summary

AFQuery currently ships with two database schema versions (v1: het_bitmap + hom_bitmap; v2: adds fail_bitmap). Because the application is still in development and has no production users, schema v1 and all compatibility code should be removed. Only v2 should remain.

## Source notes

| File | Line | Note |
|------|------|------|
| `docs/getting-started/concepts.md` | 234 | Two schema versions exist; since the app is in early development, v1 and all references to it should be removed from code and docs. |
| `docs/advanced/filter-pass-tracking.md` | 5 | Intro says "Schema version 2 introduces…"; once v1 is gone the versioning framing should be dropped entirely. |

## Scope

Code + Docs (high complexity — touches ingest, build, query, tests, and documentation).

## Steps

### 1. Audit all v1 references in the codebase

Search for every occurrence of `schema_version`, `v1`, `schema=1`, `"version": 1`, or similar patterns in `src/afquery/` and `tests/`.

### 2. Remove v1 compatibility in code

Key files expected to change (verify during audit):

- `src/afquery/preprocess/build.py` — remove any branch that emits Parquet without `fail_bitmap`.
- `src/afquery/preprocess/__init__.py` — remove `schema_version` writing logic for v1; hard-code `schema_version = 2` (or rename to a feature-flag constant).
- `src/afquery/query.py` — remove any fallback that handles absence of `fail_bitmap` column.
- `src/afquery/database.py` — remove v1 migration / detection code if present.
- `manifest.json` schema — keep `schema_version` field but only ever write `2`; remove any v1 migration path.

### 3. Remove v1 compatibility in tests

- `tests/conftest.py` — ensure the synthetic database fixture only builds v2 databases.
- Remove any test that specifically exercises v1 behaviour.

### 4. Update documentation

- `docs/getting-started/concepts.md` — replace the schema version table with a single sentence stating that the database format always includes `fail_bitmap`.
- `docs/advanced/filter-pass-tracking.md` — rewrite the introduction; remove "Schema version 2 introduces…" framing; describe fail_bitmap as a standard feature, not a versioned addition.
- Any other doc mentioning v1 or schema migration.

### 5. Verification

```bash
grep -r "schema_version\|v1\|schema=1" src/ tests/ docs/
pytest --tb=short -q
```

Both commands should return no hits for v1 references, and the test suite must pass.
