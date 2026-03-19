# Issue 04 — WES BED File Validation (Hard Error on Missing BED)

## Summary

When a WES sample's BED capture file is not found, AFQuery currently falls back silently to treating the sample as WGS (all positions covered). This is incorrect behaviour: a missing BED file should cause the pipeline to abort with a clear error message. This plan covers adding strict validation and updating the documentation to reflect it.

## Source notes

| File | Line | Note |
|------|------|------|
| `docs/faq.md` | 191 | Silent WGS fallback for missing BED is undesired; the tool must error and stop. |
| `docs/guides/manifest-format.md` | 64 | Docs say "If no BED file is found, AFQuery treats the technology as WGS"; this must be corrected once the behaviour is fixed. |

## Scope

Code + Docs.

## Steps

### 1. Locate the silent-fallback code

Search in `src/afquery/preprocess/regions.py` and `src/afquery/capture.py` for the logic that handles a missing BED path. Likely pattern:

```python
if bed_path is None or not os.path.exists(bed_path):
    # treat as WGS
    return CaptureIndex(wgs=True)
```

### 2. Replace the fallback with a hard error

For `create-db` and `update-db`:
- If a technology is WES and its BED file path is missing from the manifest, raise a `ValueError` / `FileNotFoundError` with a message like:
  ```
  BED file not found for technology 'WES_kit_A': '/path/to/kit_a.bed'.
  Provide a valid BED path in the manifest or use technology type 'WGS'.
  ```
- Perform this check **before** any VCF ingestion begins, so users get an early fail instead of discovering the problem after hours of processing.

### 3. Add a manifest validation step

In `src/afquery/preprocess/manifest.py` (or wherever manifest rows are parsed):
- After parsing all rows, collect all unique WES technologies.
- For each, verify the BED file path exists on disk.
- Accumulate all missing files and report them in a single error (don't stop at the first).

### 4. Update documentation

- `docs/guides/manifest-format.md` — replace the sentence about WGS fallback with a warning that a missing BED causes an abort; document the expected error message.
- `docs/faq.md` — update the FAQ entry to describe the new strict behaviour and how to fix it (supply the BED path).

### 5. Update tests

- Add a test in `tests/` that constructs a manifest with a WES technology pointing to a non-existent BED and asserts that `create-db` raises the expected error.
- Ensure existing tests that intentionally omit BED paths (WGS-only tests) are unaffected.

### 6. Verification

```bash
pytest --tb=short -q
# Manually: create a manifest with a bad BED path and run afquery create-db; expect a clear error.
```
