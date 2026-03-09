from dataclasses import dataclass, field


@dataclass
class Sample:
    sample_id: int      # 0-indexed; used as bit position in Roaring Bitmaps
    sample_name: str
    sex: str            # 'male' | 'female'
    tech_id: int


@dataclass
class Technology:
    tech_id: int
    tech_name: str
    bed_path: str | None  # None means WGS (always covered)


@dataclass
class VariantKey:
    chrom: str   # canonical form: 'chr1', 'chrX', etc.
    pos: int     # 1-based
    ref: str
    alt: str


@dataclass
class QueryParams:
    chrom: str
    pos: int               # 1-based
    icd10_codes: list[str]
    sex_filter: str = "both"  # 'male' | 'female' | 'both'


@dataclass
class QueryResult:
    variant: VariantKey
    AC: int
    AN: int
    AF: float | None       # None if AN == 0
    n_samples_eligible: int
