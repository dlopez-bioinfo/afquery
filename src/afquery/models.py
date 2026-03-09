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
    phenotype_codes: list[str] = field(default_factory=list)
    sex_filter: str = "both"  # 'male' | 'female' | 'both'
    ref: str | None = None   # optional filter to specific ref allele
    alt: str | None = None   # optional filter to specific alt allele


@dataclass
class QueryResult:
    variant: VariantKey
    AC: int
    AN: int
    AF: float | None       # None if AN == 0
    n_samples_eligible: int
