# Coverage Evidence

`N_HOM_REF` is computed as a residual:
`len(eligible) − N_HET − N_HOM_ALT − N_FAIL`. For WGS samples that residual is
exactly right: every covered sample without a variant call is hom-ref. For WES
samples it is a *best-effort* assumption: the BED capture region tells us a
position *could* be sequenced, but not that it *was* sequenced at adequate
depth in this particular sample. Standard variant-only VCFs do not contain
hom-ref calls, so AFQuery cannot distinguish "true hom-ref" from "no coverage"
for non-carrier WES samples.

Two opt-in mechanisms let users tighten that assumption. Together they expose
a new field, **`N_NO_COVERAGE`**, that holds samples whose hom-ref status is
not trusted under the chosen criteria. The new genotype invariant is:

```
N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE = n_eligible
```

Samples in `N_NO_COVERAGE` remain in `eligible` and `AN` (just like
`N_FAIL`), so allele frequencies stay conservative.

WGS samples are never re-classified as `N_NO_COVERAGE`. Carrier samples
(`het` / `hom` / `fail`) are never affected by these filters — only WES
non-carriers can move between `N_HOM_REF` and `N_NO_COVERAGE`.

---

## Phase 1 — Query-time, evidence-counting

No schema change. AFQuery counts existing carriers per WES tech at each
position and applies a per-tech gate.

| Flag | Counts |
|------|--------|
| `--min-pass K`     | `het ∪ hom` PASS carriers within the tech |
| `--min-observed K` | `het ∪ hom ∪ fail` carriers within the tech |

If the tech falls below either threshold, *all of its non-carrier samples* at
that position move from `N_HOM_REF` to `N_NO_COVERAGE`. When both flags are
set, both must hold (AND). Default `0` ⇒ no filtering, identical to legacy
behaviour.

```bash
afquery query --db ./db/ --locus chr1:925952 --min-pass 1
```

Cost: a few extra bitmap intersections per WES tech per position. Suitable for
existing databases without re-creation.

---

## Phase 2 — Build-time, quality-aware

Phase 2 stores the result of a quality decision so the query layer does no
extra work. It requires a one-time creation with quality thresholds:

| Flag (`create-db`) | Effect |
|--------------------|--------|
| `--min-dp D`     | Minimum `FORMAT/DP` for a carrier to count as quality evidence. |
| `--min-gq G`     | Minimum `FORMAT/GQ` for a carrier to count as quality evidence. |
| `--min-qual Q`   | Minimum `QUAL` field for a carrier to count as quality evidence. |
| `--min-covered K`| Per WES tech, position is "trusted" only if at least K carriers pass the quality thresholds. |

Two new Parquet columns are written per `(chrom, pos, ref, alt)`:

- `quality_pass_bitmap` — carriers that meet `min_dp` AND `min_gq` AND `min_qual`.
- `filtered_bitmap` — non-carrier WES samples whose tech failed the
  `min_covered` gate.

`schema_version` is bumped to `3.0` and the chosen thresholds are recorded in
`manifest.json` under `coverage_filter`. They are immutable; they apply to
samples added later via `update-db --add-samples`.

At query time `filtered_bitmap` is intersected with `eligible` and added to
`N_NO_COVERAGE` automatically — no additional flag needed.

```bash
afquery create-db \
  --manifest samples.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --bed-dir ./beds/ \
  --min-dp 30 --min-gq 20 --min-covered 1
```

### Query-time companion: `--min-quality-evidence K`

Only valid against `schema_version ≥ 3.0` databases. Tightens the build-time
gate: at query time, a WES tech needs ≥K samples in `quality_pass_bitmap`,
otherwise its remaining (non-carrier, not-already-filtered) samples join
`N_NO_COVERAGE`.

```bash
afquery query --db ./db/ --locus chr1:925952 --min-quality-evidence 5
```

Using this flag against an older DB raises a `ValueError` with a clear message.

---

## Combining Phase 1 and Phase 2

The two phases are additive. `N_NO_COVERAGE` is the *union* of:

1. Stored `filtered_bitmap & eligible` (Phase 2, automatic).
2. Phase 1 dynamic filtering driven by `--min-pass` / `--min-observed`.
3. Phase 2 query-time tightening driven by `--min-quality-evidence`.

Carriers are never included; samples cannot be double-counted.

---

## Choosing thresholds

- **Pure genotyping (no quality info available)**: use `--min-pass 1` or
  `--min-observed 1` at query time. No DB rebuild needed. Conservative;
  positions that were probably sequenced but happen to have zero PASS calls
  in your cohort flip to `N_NO_COVERAGE`.
- **Real cohorts with DP/GQ available**: rebuild with
  `--min-dp 20 --min-gq 20 --min-covered 1`. Carriers with low confidence stop
  validating positions. `N_NO_COVERAGE` rises only where the cohort signal is
  truly weak.
- **High-stakes clinical use**: layer on `--min-quality-evidence 3` (or
  similar) at query time to demand multiple independent quality calls per
  tech.

---

## Output channels

- `query` / `query_region` / `query_batch_multi` / `dump` / `annotate` all
  expose `N_NO_COVERAGE` as a first-class field/column, plus per-group variants
  in `dump` (`N_NO_COVERAGE_<label>`).
- `variant_info` lists individual filtered samples with `genotype = "no_coverage"`.
- The annotate INFO field is `AFQUERY_N_NO_COVERAGE` (Number=A, one entry per
  ALT allele).

---

## What it does *not* do

- It does not fabricate hom-ref calls. Samples that lack any VCF record at a
  position remain invisible to AFQuery; the only choice is whether to assume
  hom-ref or fall back to `N_NO_COVERAGE`.
- It does not import per-sample BAM coverage. Cohorts that need true per-sample
  coverage tracking should provide per-sample BEDs (see
  [issue #29](https://github.com/dlopez-bioinfo/afquery/issues/29)).
- It does not change `AC`. Allele counts are still computed only from
  carriers that survive ploidy adjustment.
