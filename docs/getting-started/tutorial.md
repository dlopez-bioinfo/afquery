# Tutorial: End-to-End Walkthrough

This tutorial walks through every major AFQuery feature using a synthetic demo dataset. By the end, you will have built a database, queried variants, filtered by metadata, annotated a VCF, and exported results.

!!! tip "Prerequisites"
    Make sure AFQuery is installed. See [Installation](installation.md) if needed.

---

## 1. Generate Demo Data

AFQuery ships with a script that creates 10 synthetic VCFs, a manifest, and BED files for two WES technologies:

```bash
python examples/demo/create_demo_data.py
```

This creates `examples/demo/demo_output/` with:

- `vcfs/` — 10 single-sample VCFs (DEMO_001 through DEMO_010)
- `beds/` — capture BED files for `wes_v1` and `wes_v2`
- `manifest.tsv` — sample metadata

The demo cohort has 4 WGS samples, 3 wes_v1 samples, and 3 wes_v2 samples, with phenotype codes `E11.9`, `I10`, and `control`.

---

## 2. Create the Database

```bash
afquery create-db \
  --manifest examples/demo/demo_output/manifest.tsv \
  --output-dir ./demo_db/ \
  --genome-build GRCh38 \
  --bed-dir examples/demo/demo_output/beds/
```

This ingests all VCFs, builds Roaring Bitmap Parquet files, and writes `manifest.json` and `metadata.sqlite`.

---

## 3. Inspect the Database

```bash
afquery info --db ./demo_db/
```

You should see:

- **10 samples** across 3 technologies (wgs, wes_v1, wes_v2)
- **Chromosomes** present in the database
- **Phenotype codes** registered (E11.9, I10, control)

---

## 4. Query a Single Variant

Pick any position from the database and query it:

```bash
afquery query --db ./demo_db/ --locus chr1:925000
```

The output shows AC, AN, AF, and per-genotype counts. If the position is not in the database, you will see `AN=0` — try a region query to find nearby variants:

```bash
afquery query --db ./demo_db/ --region chr1:900000-950000
```

---

## 5. Filter by Sex

Query only female samples:

```bash
afquery query --db ./demo_db/ --region chr1:900000-950000 --sex female
```

Compare AN between `--sex female`, `--sex male`, and the default (both). On autosomes, AN reflects the number of eligible samples times two.

---

## 6. Filter by Phenotype

Query samples tagged with `E11.9`:

```bash
afquery query --db ./demo_db/ --region chr1:900000-950000 --phenotype E11.9
```

Exclude control samples:

```bash
afquery query --db ./demo_db/ --region chr1:900000-950000 --phenotype ^control
```

The `^` prefix means "exclude". See [Sample Filtering](../guides/sample-filtering.md) for the full syntax.

---

## 7. Filter by Technology

Restrict to WGS samples only:

```bash
afquery query --db ./demo_db/ --region chr1:900000-950000 --tech wgs
```

Compare with WES-only queries — note how AN changes based on capture coverage:

```bash
afquery query --db ./demo_db/ --region chr1:900000-950000 --tech wes_v1
```

---

## 8. Combine Filters

All filter dimensions compose with AND:

```bash
afquery query \
  --db ./demo_db/ \
  --region chr1:900000-950000 \
  --sex female \
  --phenotype E11.9 \
  --tech wgs
```

This selects only female WGS samples tagged with E11.9.

---

## 9. Annotate a VCF

Use one of the demo VCFs as input:

```bash
afquery annotate \
  --db ./demo_db/ \
  --input examples/demo/demo_output/vcfs/DEMO_001.vcf.gz \
  --output ./annotated_demo.vcf
```

The output VCF gains `AFQUERY_AC`, `AFQUERY_AN`, `AFQUERY_AF`, and other INFO fields. Inspect the result:

```bash
grep -v "^##" ./annotated_demo.vcf | head -20
```

See [Annotate a VCF](../guides/annotate-vcf.md) for filtering and downstream usage.

---

## 10. Bulk Export with Dump

Export all variant frequencies to CSV:

```bash
afquery dump --db ./demo_db/ --output demo_dump.csv
```

Disaggregate by sex and technology:

```bash
afquery dump --db ./demo_db/ --output demo_dump_stratified.csv --by-sex --by-tech
```

---

## 11. Interpret Results with ACMG Criteria

With the annotated VCF or query results, you can apply ACMG criteria:

- **BA1**: Is AF > 5% in your cohort? → Variant is benign.
- **PM2_Supporting**: Is the variant absent (AC=0) with high AN? → Supporting pathogenic evidence.
- **PS4**: Is the variant enriched in cases vs. controls? Use `--phenotype` and `--phenotype ^disease_code` to compare.

For detailed guidance, see [ACMG Criteria (BA1/PM2/PS4)](../use-cases/acmg-use-cases.md).

---

## Next Steps

- [Key Concepts](concepts.md) — understand how bitmaps, Parquet, and metadata filtering work
- [Manifest Format](../guides/manifest-format.md) — prepare your own cohort manifest
- [Create a Database](../guides/create-database.md) — build a database from your real data
- [ACMG Criteria](../use-cases/acmg-use-cases.md) — apply local AF to variant classification
