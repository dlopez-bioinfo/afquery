# Bulk Export (dump)

`afquery dump` exports allele frequency data to CSV. It supports filtering by region and disaggregating output by sex, technology, or phenotype group.

---

## Basic Usage

Export all variants to stdout:

```bash
afquery dump --db ./db/
```

Write to a file:

```bash
afquery dump --db ./db/ --output all_variants.csv
```

---

## Filter by Region

Export a single chromosome:

```bash
afquery dump --db ./db/ --chrom chr1 --output chr1.csv
```

Export a specific region:

```bash
afquery dump --db ./db/ --chrom chr1 --start 900000 --end 1000000 --output region.csv
```

Positions are 1-based, inclusive on both ends.

---

## Sample Filtering

Apply the same sample filters as `query`:

```bash
afquery dump --db ./db/ \
  --phenotype E11.9 \
  --sex female \
  --tech wgs \
  --output diabetic_female_wgs.csv
```

---

## Disaggregate by Sex

Add separate AC/AN columns for male and female:

```bash
afquery dump --db ./db/ --by-sex --output by_sex.csv
```

Output columns:
```
chrom,pos,ref,alt,AC,AN,AF,AC_male,AN_male,AF_male,AC_female,AN_female,AF_female
```

---

## Disaggregate by Technology

Add separate columns per sequencing technology:

```bash
afquery dump --db ./db/ --by-tech --output by_tech.csv
```

Output columns include `AC_wgs`, `AN_wgs`, `AF_wgs`, `AC_wes_v1`, etc.

---

## Disaggregate by Phenotype

Add separate columns for specific phenotype groups:

```bash
afquery dump --db ./db/ \
  --by-phenotype E11.9 \
  --by-phenotype I10 \
  --output by_phenotype.csv
```

Output includes `AC_E11.9`, `AN_E11.9`, `AF_E11.9`, `AC_I10`, etc.

NOTA_PARA_LLM: Disaggregate by Sex, Disaggregate by Technology y Disaggregate by Phenotype pueden ser encapsulados dentro de una misma sección, ya que el concepto es el mismo. De esa manera es más clara la ayuda.

---

## All Groups (Cartesian Product)

Disaggregate by all combinations of sex × technology × phenotype simultaneously:

```bash
afquery dump --db ./db/ --all-groups --output all_groups.csv
```

!!! warning
    The number of columns grows as a Cartesian product. Disaggregating by all combinations may produce very wide CSV file.

---

## Parallelism

Use multiple threads for faster export:

```bash
afquery dump --db ./db/ --threads 16 --output all_variants.csv
```

---


## Full Option Reference

See [CLI Reference → dump](../reference/cli.md#dump).
