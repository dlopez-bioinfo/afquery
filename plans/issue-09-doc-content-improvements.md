# Issue 09 — Documentation Content Improvements

## Summary

A collection of content-level improvements spread across use-case guides, getting-started pages, and reference docs. Each item is a self-contained editorial fix — no code changes required.

## Source notes

| File | Line | Note |
|------|------|------|
| `docs/use-cases/population-specific-af.md` | 22 | "Step-by-Step Example" section is duplicated across all clinical workflow guides; find a way to reduce redundancy. |
| `docs/use-cases/technology-integration.md` | 33 | Same Step-by-Step duplication. |
| `docs/use-cases/cohort-stratification.md` | 22 | Same Step-by-Step duplication. |
| `docs/use-cases/sex-specific-af.md` | 28 | Same Step-by-Step duplication. |
| `docs/use-cases/clinical-prioritization.md` | 13 | Explain population mismatch via genetic drift (verify with literature). |
| `docs/use-cases/clinical-prioritization.md` | 20 | Extend and formalise the local-artifacts explanation. |
| `docs/getting-started/preprocessing.md` | 37 | Note that INFO field stripping also prevents malformed VCFs from breaking the pipeline. |
| `docs/getting-started/understanding-output.md` | 143 | "Next Steps" links are missing descriptions; add a one-sentence description to each link. |
| `docs/guides/create-database.md` | 129 | "Verify the Build" section implies a mandatory check; rewrite to frame `afquery check` and `afquery info` as optional verification commands. |
| `docs/guides/dump-export.md` | 95 | "Disaggregate by Sex", "Disaggregate by Technology", and "Disaggregate by Phenotype" are three separate sections for the same concept; consolidate into one "Disaggregate output" section. |
| `docs/getting-started/concepts.md` | 210 | Replace `ICD-10` with `International Classification of Diseases (ICD) codes` everywhere in docs and code. |
| `docs/getting-started/concepts.md` | 212 | Replace `HPO` with `Human Phenotype Ontology (HPO) terms` (introduce acronym on first use). |
| `docs/getting-started/concepts.md` | 214 | Replace `GO` with `Online Catalog of Human Genes and Genetic Disorders (OMIM) entries`. |
| `docs/advanced/ploidy-and-sex-chroms.md` | 5 | chrMT is not a sex chromosome; update all docs that claim AFQuery is "ploidy-aware for sex chromosomes" to say "ploidy-aware for sex chromosomes and the mitochondrial chromosome". |
| `docs/faq.md` | 155 | Document the warning system: review the code to identify all warning conditions and add a section to the FAQ or a dedicated page. |
| `docs/advanced/debugging-results.md` | 114 | The table mentions `--include-all-filters`; audit whether this flag exists in the code. If it does, document it properly everywhere. If it does not, remove all mentions. |

## Scope

Docs only (except terminology replacements in `ICD-10`/`HPO`/`GO` which may touch Python source strings — verify and update those too).

## Steps

### 1. Deduplicate "Step-by-Step Example" sections

**Options (discuss and choose one):**
- **Option A — Extract to shared page**: Create `docs/guides/step-by-step-query.md` with the canonical example; in each use-case guide, replace the duplicate section with a short intro + `See [Step-by-Step Query Guide](../guides/step-by-step-query.md)`.
- **Option B — Use MkDocs snippets / includes**: Use the `--8<--` include syntax if the MkDocs config already uses `pymdownx.snippets`; define the example once in a `docs/_partials/` file and include it.
- **Option C — Minimal approach**: Reduce each use-case guide's example to 2–3 lines specific to that use case, with a pointer to the full query guide for details.

Files to update: `population-specific-af.md`, `technology-integration.md`, `cohort-stratification.md`, `sex-specific-af.md`, `clinical-prioritization.md`.

### 2. Improve clinical-prioritization explanations

- **Line 13**: Add a paragraph explaining that population AF differences can arise from genetic drift — the stochastic fixation or loss of alleles in isolated subpopulations. Cite or paraphrase evidence from the literature (e.g. gnomAD population-stratification papers).
- **Line 20**: Extend the local-artifacts section to explain sequencing platform artefacts (e.g. strand bias, systematic errors in repetitive regions) and why a local cohort sharing the same sequencing pipeline is the best control.

### 3. Update preprocessing page (line 37)

Add a sentence after the INFO-stripping description: stripping non-standard INFO fields also prevents malformed VCFs (e.g. with duplicate or mis-typed INFO keys) from causing parse errors during ingestion.

### 4. Add descriptions to "Next Steps" links (`understanding-output.md` line 143)

Read the section around line 143 and add a one-sentence description to each bare link, e.g.:
```markdown
- [Sample Filtering](../guides/sample-filtering.md) — restrict queries to phenotype, sex, or technology subsets.
```

### 5. Rewrite "Verify the Build" section (`create-database.md` line 129)

Replace imperative "Verify the Build" framing with: "After the database is created, you can optionally run the following commands to inspect its contents or check integrity: …"

### 6. Consolidate Disaggregate sections (`dump-export.md` line 95)

Merge "Disaggregate by Sex", "Disaggregate by Technology", and "Disaggregate by Phenotype" into a single section: "Disaggregate Output". Show all three flags (`--by-sex`, `--by-tech`, `--by-phenotype`) together with a combined example.

### 7. Terminology replacements (global)

Run a search-and-replace across `docs/` and `src/`:
- `ICD-10` → `International Classification of Diseases (ICD) codes` (on first use in each file; subsequent uses can be `ICD`)
- `HPO` → `Human Phenotype Ontology (HPO) terms` (introduce acronym; subsequent uses can be `HPO`)
- `GO` → `Online Catalog of Human Genes and Genetic Disorders (OMIM) entries` (note: confirm that "GO" in the codebase refers to OMIM and not Gene Ontology before replacing)

### 8. Clarify chrMT ploidy scope

Search all docs for "ploidy-aware for sex chromosomes" or "sex chromosomes" in the context of ploidy. Update each occurrence to include the mitochondrial chromosome, e.g.: "ploidy-aware for sex chromosomes (chrX, chrY) and the mitochondrial chromosome (chrMT)".

### 9. Document the warning system (`faq.md` line 155)

- Read `src/afquery/` to identify all `warnings.warn(...)` or equivalent calls.
- Add a new FAQ entry (or a dedicated `docs/advanced/warnings.md` page) listing each warning, its cause, and the recommended action.

### 10. Audit `--include-all-filters` (`debugging-results.md` line 114)

```bash
grep -r "include.all.filters\|include_all_filters" src/ tests/ docs/
```
- If the flag **exists**: document it properly in `docs/reference/cli.md` and the debugging guide.
- If the flag **does not exist**: remove every mention from the documentation.

### 11. Verification

```bash
mkdocs build  # no broken links
grep -r "ICD-10\b" docs/ src/  # should return no results
grep -r "\bGO\b" docs/ src/    # verify only OMIM-related hits remain
```
