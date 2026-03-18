# Query Allele Frequencies

`afquery query` retrieves allele frequencies from the database. Three query modes are available: point, region, and batch.

---

## Point Query

Query a single genomic position:

```bash
afquery query --db ./db/ --chrom chr1 --pos 925952
```

Output (default `text` format):
```
chr1:925952
  REF=G  ALT=A  AC=142  AN=2742  AF=0.0518  N_HET=138  N_HOM_ALT=2
```

Filter to a specific ref/alt allele (useful at multi-allelic sites):

```bash
afquery query --db ./db/ --chrom chr1 --pos 925952 --ref G --alt A
```

---

## Region Query

Query all variants in a genomic range:

```bash
afquery query --db ./db/ --chrom chr1 --region 900000-1000000
```

The range is 1-based, inclusive on both ends.

---

## Batch Query

Query multiple positions at once from a file:

```bash
afquery query --db ./db/ --chrom chr1 --from-file variants.tsv
```

The input file is a **headerless TSV** with columns `pos ref alt`:

```tsv
925952	G	A
1014541	C	T
1020172	A	G
```

---

## Output Formats

### text (default)

Human-readable, one block per variant:

```
chr1:925952
  REF=G  ALT=A  AC=142  AN=2742  AF=0.0518  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1231
```

### tsv

Tab-separated, one row per variant, suitable for downstream processing:

```bash
afquery query --db ./db/ --chrom chr1 --region 900000-1000000 --format tsv
```

```
chrom	pos	ref	alt	AC	AN	AF	N_HET	N_HOM_ALT	N_HOM_REF
chr1	925952	G	A	142	2742	0.0518	138	2	1231
```

### json

JSON array, one object per variant:

```bash
afquery query --db ./db/ --chrom chr1 --pos 925952 --format json
```

```json
[
  {
    "chrom": "chr1",
    "pos": 925952,
    "ref": "G",
    "alt": "A",
    "AC": 142,
    "AN": 2742,
    "AF": 0.05178,
    "N_HET": 138,
    "N_HOM_ALT": 2,
    "N_HOM_REF": 1231
  }
]
```

---

## Sample Filtering

All query modes support the same filter options:

```bash
afquery query \
  --db ./db/ \
  --chrom chr1 \
  --pos 925952 \
  --phenotype E11.9 \
  --sex female \
  --tech wgs
```

Filters compose with AND:
- `--phenotype E11.9 --sex female` = female samples with E11.9

Multiple values for the same filter compose with OR:
- `--phenotype E11.9 --phenotype I10` = samples with E11.9 OR I10

Exclude with `^` prefix:
- `--tech ^wes_v1` = all technologies except `wes_v1`

See [Sample Filtering](sample-filtering.md) for full syntax.

---

## Results When AN=0

If all samples are excluded by your filters, the result will have `AC=0`, `AN=0`, and `AF=None`. This is expected behavior — it means no samples in the selected subgroup were eligible at that position (e.g., all WES samples and the position is not in their capture regions).

---

## Comparing AF Across Subgroups

Run two queries and compare:

```bash
# Diabetic patients
afquery query --db ./db/ --chrom chr1 --pos 925952 --phenotype E11.9 --format json

# Healthy controls (exclude diabetic)
afquery query --db ./db/ --chrom chr1 --pos 925952 --phenotype ^E11.9 --format json
```

For systematic comparison across many variants, consider [Bulk Export](dump-export.md) with `--by-phenotype`.

---

## Full Option Reference

See [CLI Reference → query](../reference/cli.md#query).
