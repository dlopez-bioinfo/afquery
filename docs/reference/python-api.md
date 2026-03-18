# Python API

AFQuery exposes a clean Python API through the `Database` class. All CLI functionality is available programmatically.

---

## Installation

```python
from afquery import Database
```

---

## Database

### Constructor

```python
db = Database(db_path: str)
```

Opens an AFQuery database at the given path. The manifest and sample metadata are loaded on initialization.

```python
db = Database("./my_db/")
```

---

### query

```python
db.query(
    chrom: str,
    pos: int,
    phenotype: list[str] | None = None,
    sex: str = "both",
    ref: str | None = None,
    alt: str | None = None,
    tech: list[str] | None = None,
) -> list[QueryResult]
```

Query allele frequencies at a single genomic position.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | str | Chromosome name (e.g., `"chr1"`, `"chrX"`) |
| `pos` | int | 1-based genomic position |
| `phenotype` | list[str] \| None | Phenotype filter codes. Use `"^CODE"` prefix to exclude. |
| `sex` | str | `"both"` (default), `"male"`, or `"female"` |
| `ref` | str \| None | Filter to specific reference allele |
| `alt` | str \| None | Filter to specific alternate allele |
| `tech` | list[str] \| None | Technology filter. Use `"^TECH"` prefix to exclude. |

**Returns:** List of `QueryResult` objects, one per variant at the position (sorted by `(pos, alt)`).

**Example:**

```python
results = db.query("chr1", pos=925952)
for r in results:
    print(f"{r.variant.ref}>{r.variant.alt}  AC={r.AC}  AN={r.AN}  AF={r.AF:.4f}")

# With filters
results = db.query(
    chrom="chr1",
    pos=925952,
    phenotype=["E11.9"],
    sex="female",
    tech=["wgs"],
)
```

---

### query_region

```python
db.query_region(
    chrom: str,
    start: int,
    end: int,
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
) -> list[QueryResult]
```

Query allele frequencies for all variants in a genomic range.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | str | Chromosome name |
| `start` | int | 1-based start position (inclusive) |
| `end` | int | 1-based end position (inclusive) |
| `phenotype` | list[str] \| None | Phenotype filter |
| `sex` | str | Sex filter |
| `tech` | list[str] \| None | Technology filter |

**Example:**

```python
results = db.query_region("chr1", start=900000, end=1000000)
print(f"Found {len(results)} variants")
```

---

### query_batch

```python
db.query_batch(
    chrom: str,
    variants: list[tuple[int, str, str]],
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
) -> list[QueryResult]
```

Query allele frequencies for a list of specific variants.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | str | Chromosome name |
| `variants` | list[tuple[int, str, str]] | List of `(pos, ref, alt)` tuples |
| `phenotype` | list[str] \| None | Phenotype filter |
| `sex` | str | Sex filter |
| `tech` | list[str] \| None | Technology filter |

**Example:**

```python
variants = [(925952, "G", "A"), (1014541, "C", "T"), (1020172, "A", "G")]
results = db.query_batch("chr1", variants=variants)
```

---

### annotate_vcf

```python
db.annotate_vcf(
    input_vcf: str,
    output_vcf: str,
    phenotype: list[str] | None = None,
    sex: str = "both",
    n_workers: int | None = None,
    tech: list[str] | None = None,
) -> dict
```

Annotate a VCF file with allele frequency INFO fields.

**Returns:** Stats dict:
```python
{
    "n_variants": int,    # total variants in input VCF
    "n_annotated": int,   # variants with at least one allele found in DB
    "n_uncovered": int,   # variants with no allele found in DB
}
```

---

### dump

```python
db.dump(
    output: str | None = None,
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
    by_sex: bool = False,
    by_tech: bool = False,
    by_phenotype: list[str] | None = None,
    all_groups: bool = False,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    n_workers: int | None = None,
) -> dict
```

Export allele frequency data to CSV. If `output` is None, writes to stdout.

---

### add_samples

```python
db.add_samples(
    manifest_path: str,
    threads: int = 8,
    tmp_dir: str | None = None,
    bed_dir: str | None = None,
    genome_build: str | None = None,
) -> dict
```

Add samples from a manifest TSV. Returns a stats dict.

---

### remove_samples

```python
db.remove_samples(sample_names: list[str]) -> dict
```

Remove samples by name. Returns a stats dict with `n_removed`.

---

### compact

```python
db.compact() -> dict
```

Compact the database to reclaim space from removed samples.

---

### info

```python
db.info() -> dict
```

Return database metadata as a dict.

---

### list_samples

```python
db.list_samples() -> list[dict]
```

Return a list of all samples with their metadata (name, sex, tech, phenotypes).

---

### check

```python
db.check() -> list
```

Validate database integrity. Returns `list[CheckResult]`. Each item has `.severity` (`"error"`, `"warning"`, or `"info"`) and `.message` (str). An empty list means the database is healthy.

---

### changelog

```python
db.changelog(limit: int | None = None) -> list[dict]
```

Return changelog history. Each entry is a dict:
```python
{
    "event_id": int,
    "event_type": str,        # "preprocess", "add_samples", "remove_samples", "compact"
    "event_time": str,        # ISO datetime string
    "sample_names": list[str] | None,
    "notes": str | None,
}
```

---

### set_db_version

```python
db.set_db_version(version: str) -> None
```

Set the database version label.

---

### get_all_phenotypes

```python
db.get_all_phenotypes() -> list[str]
```

Return all distinct phenotype codes present in the database.

---

## QueryResult

```python
@dataclass
class QueryResult:
    variant: VariantKey      # chrom, pos, ref, alt
    AC: int                  # Allele count
    AN: int                  # Allele number
    AF: float | None         # Allele frequency (None if AN == 0)
    n_samples_eligible: int  # Number of eligible samples at this position
    N_HET: int               # Heterozygous count
    N_HOM_ALT: int           # Homozygous alt count
    N_HOM_REF: int           # Homozygous ref count
    N_FAIL: int              # Samples with alt allele called but FILTER≠PASS (schema v2 only; 0 for v1 databases — not None)
```

### VariantKey

```python
@dataclass
class VariantKey:
    chrom: str   # Canonical form: 'chr1', 'chrX', etc.
    pos: int     # 1-based
    ref: str
    alt: str
```

---

## SampleFilter

```python
@dataclass
class SampleFilter:
    phenotype_include: list[str] = []   # Empty = all samples
    phenotype_exclude: list[str] = []
    tech_include: list[str] = []        # Empty = all samples
    tech_exclude: list[str] = []
    sex: str = "both"                   # 'male' | 'female' | 'both'

    @staticmethod
    def parse(
        phenotype_tokens: list[str],
        tech_tokens: list[str],
        sex: str = "both",
    ) -> SampleFilter
```

`SampleFilter.parse` handles the `^` prefix exclusion syntax:

```python
from afquery.models import SampleFilter

sf = SampleFilter.parse(
    phenotype_tokens=["E11.9", "^I10"],
    tech_tokens=["wgs"],
    sex="female",
)
# sf.phenotype_include = ["E11.9"]
# sf.phenotype_exclude = ["I10"]
# sf.tech_include = ["wgs"]
# sf.sex = "female"
```

---

## Full Example

```python
from afquery import Database

db = Database("./my_db/")

# Point query
results = db.query("chr1", pos=925952, phenotype=["E11.9"], sex="female")
for r in results:
    print(f"AF={r.AF:.4f}  AC={r.AC}/{r.AN}  HET={r.N_HET}  HOM={r.N_HOM_ALT}")

# Region query
region_results = db.query_region("chrX", start=154931044, end=155270560)
print(f"Found {len(region_results)} variants in PAR2")

# Batch query
variants = [(925952, "G", "A"), (1014541, "C", "T")]
batch_results = db.query_batch("chr1", variants=variants)

# Database info
meta = db.info()
print(f"Samples: {meta['n_samples']}, Build: {meta['genome_build']}")

# List all phenotypes
phenotypes = db.get_all_phenotypes()
print("Phenotypes:", phenotypes)
```
