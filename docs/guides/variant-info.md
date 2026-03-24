# Variant Info

`afquery variant-info` returns the list of samples carrying a specific variant, together with each sample's metadata: sex, sequencing technology, phenotype codes, genotype (het/hom), and FILTER status (PASS/FAIL).

This avoids the need to re-query raw VCF files when inspecting individual variant carriers.

---

## Basic usage

```bash
afquery variant-info --db ./db/ --locus chr1:925952
```

By default all samples are queried and results are printed as an aligned text table:

```
sample_id  sample_name  sex     tech       phenotypes    genotype  filter
---------  -----------  ------  ---------  ------------  --------  ------
3          P003         male    WGS        E11.9,J45     het       PASS
17         P017         female  WES_kit_A  E11.9         hom       PASS
42         P042         male    WGS        I10           alt       FAIL
```

---

## Filtering to a specific allele

When multiple alleles exist at the same position, use `--ref` and `--alt` to restrict to one:

```bash
afquery variant-info --db ./db/ --locus chr17:41245466 --ref A --alt T
```

Without `--ref`/`--alt`, carriers for all alleles at the locus are returned and a warning is emitted if more than one allele is found.

---

## Sample filters

`variant-info` accepts the same sample filters as `query`:

```bash
# Only female carriers
afquery variant-info --db ./db/ --locus chr1:925952 --sex female

# Only carriers with phenotype E11.9
afquery variant-info --db ./db/ --locus chr1:925952 --phenotype E11.9

# Exclude phenotype I10
afquery variant-info --db ./db/ --locus chr1:925952 --phenotype ^I10

# Restrict to WGS samples
afquery variant-info --db ./db/ --locus chr1:925952 --tech WGS

# Combine filters
afquery variant-info --db ./db/ --locus chr1:925952 \
  --sex female --phenotype E11.9 --tech WGS,WES_kit_A
```

See [Sample Filtering](sample-filtering.md) for the full filter syntax.

---

## Output formats

### TSV

Machine-readable tab-separated output, suitable for downstream processing:

```bash
afquery variant-info --db ./db/ --locus chr1:925952 --format tsv > carriers.tsv
```

```
sample_id	sample_name	sex	tech	phenotypes	genotype	filter
3	P003	male	WGS	E11.9,J45	het	PASS
17	P017	female	WES_kit_A	E11.9	hom	PASS
42	P042	male	WGS	I10	alt	FAIL
```

### JSON

Structured output with variant metadata and a sample list:

```bash
afquery variant-info --db ./db/ --locus chr1:925952 --format json
```

```json
{
  "variant": {
    "chrom": "chr1",
    "pos": 925952,
    "ref": ".",
    "alt": "."
  },
  "samples": [
    {
      "sample_id": 3,
      "sample_name": "P003",
      "sex": "male",
      "tech": "WGS",
      "phenotypes": ["E11.9", "J45"],
      "genotype": "het",
      "filter": "PASS"
    },
    {
      "sample_id": 42,
      "sample_name": "P042",
      "sex": "male",
      "tech": "WGS",
      "phenotypes": ["I10"],
      "genotype": "alt",
      "filter": "FAIL"
    }
  ]
}
```

When `--ref` and `--alt` are specified, the `variant` block contains the actual alleles. Otherwise, `"."` is used as a placeholder.

---

## Genotype values

| Value | Meaning |
|---|---|
| `het` | Heterozygous carrier, FILTER=PASS |
| `hom` | Homozygous alt carrier, FILTER=PASS |
| `alt` | Non-ref carrier with FILTER≠PASS (ploidy unknown) |

---

## Python API

```python
from afquery import Database

db = Database("./db/")

# All carriers at a locus
carriers = db.variant_info("chr1", 925952)

# With allele and sample filters
carriers = db.variant_info(
    chrom="chr17",
    pos=41245466,
    ref="A",
    alt="T",
    phenotype=["E11.9"],
    sex="female",
    tech=["WGS"],
)

for c in carriers:
    print(c.sample_name, c.genotype, c.filter_pass, c.phenotypes)
```

The `variant_info` function is also available at the top level for one-off use:

```python
from afquery import variant_info

carriers = variant_info("./db/", "chr1", 925952, ref="G", alt="A")
```

---

## All options

| Option | Default | Description |
|---|---|---|
| `--db` | required | Path to database directory |
| `--locus` | required | `CHROM:POS` (e.g. `chr1:925952`) |
| `--ref` | — | Filter to specific reference allele |
| `--alt` | — | Filter to specific alternate allele |
| `--phenotype` | all | Include phenotype (repeatable; `^CODE` excludes) |
| `--sex` | `both` | `male`, `female`, or `both` |
| `--tech` | all | Include technology (repeatable; `^NAME` excludes) |
| `--format` | `text` | `text`, `tsv`, or `json` |
| `--no-warn` | off | Suppress `AfqueryWarning` messages |
