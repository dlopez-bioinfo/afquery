"""variant_info — retrieve sample carriers for a specific variant.

Thin convenience wrapper around :class:`~afquery.database.Database` for
programmatic use without instantiating the class directly.
"""
from .database import Database
from .models import SampleCarrier


def variant_info(
    db_path: str,
    chrom: str,
    pos: int,
    ref: str | None = None,
    alt: str | None = None,
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
) -> list[SampleCarrier]:
    """Return all samples carrying the given variant, with their metadata.

    Args:
        db_path: Path to the AFQuery database directory.
        chrom: Chromosome (e.g. ``"chr1"`` or ``"1"``).
        pos: Position, 1-based.
        ref: Reference allele filter. If omitted, all alleles at *pos* are returned.
        alt: Alternate allele filter. If omitted, all alleles at *pos* are returned.
        phenotype: Phenotype filter tokens. Use ``"^CODE"`` prefix to exclude.
        sex: ``"both"`` (default), ``"male"``, or ``"female"``.
        tech: Technology filter tokens. Use ``"^TECH"`` prefix to exclude.

    Returns:
        List of :class:`~afquery.models.SampleCarrier` sorted by
        ``sample_id``.  Empty list if the variant is absent or no eligible
        carrier exists.
    """
    db = Database(db_path)
    return db.variant_info(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        phenotype=phenotype,
        sex=sex,
        tech=tech,
    )


__all__ = ["variant_info"]
