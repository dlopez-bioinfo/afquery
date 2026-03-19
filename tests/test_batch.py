"""Tests for QueryEngine.query_batch() and query_region() (Phase 3)."""
import pytest
from afquery.database import Database


def _db(test_db):
    return Database(test_db)


# Known variants from conftest.py VARIANTS fixture
_V1 = (1500, "A", "T")
_V2 = (3500, "G", "C")
_V3 = (5000, "T", "G")


# ---------------------------------------------------------------------------
# query_batch — correctness
# ---------------------------------------------------------------------------

def test_query_batch_matches_single_query(test_db):
    db = _db(test_db)
    single = (
        db.query("chr1", 1500, ["E11.9"])
        + db.query("chr1", 3500, ["E11.9"])
        + db.query("chr1", 5000, ["E11.9"])
    )
    batch = db.query_batch("chr1", [_V1, _V2, _V3], ["E11.9"])

    single_idx = {(r.variant.pos, r.variant.alt): r for r in single}
    batch_idx  = {(r.variant.pos, r.variant.alt): r for r in batch}

    assert set(single_idx) == set(batch_idx)
    for key in single_idx:
        assert single_idx[key].AC == batch_idx[key].AC
        assert single_idx[key].AN == batch_idx[key].AN
        assert abs(single_idx[key].AF - batch_idx[key].AF) < 1e-9


def test_query_batch_deduplicates_variants(test_db):
    db = _db(test_db)
    results = db.query_batch("chr1", [_V1, _V1, _V2], ["E11.9"])
    positions = [r.variant.pos for r in results]
    assert positions.count(1500) == 1


def test_query_batch_sorted_result(test_db):
    db = _db(test_db)
    results = db.query_batch("chr1", [_V3, _V1, _V2], ["E11.9"])
    positions = [r.variant.pos for r in results]
    assert positions == sorted(positions)


def test_query_batch_empty_no_matching_position(test_db):
    db = _db(test_db)
    results = db.query_batch("chr1", [(99999999, "A", "T")], ["E11.9"])
    assert results == []


def test_query_batch_all_uncovered_phenotype(test_db):
    """Unknown Phenotype → eligible bitmap empty → AN=0 → no results."""
    db = _db(test_db)
    results = db.query_batch("chr1", [_V1, _V2], ["ZZUNKNOWN"])
    assert results == []


def test_query_batch_chrom_normalisation(test_db):
    """chr-prefix should be normalised transparently."""
    db = _db(test_db)
    r_prefix = db.query_batch("chr1", [_V1], ["E11.9"])
    r_bare   = db.query_batch("1",    [_V1], ["E11.9"])
    assert len(r_prefix) == len(r_bare)
    if r_prefix:
        assert r_prefix[0].AC == r_bare[0].AC


def test_query_batch_sex_filter(test_db):
    """Female-only batch should give same result as single query with female filter."""
    db = _db(test_db)
    single = db.query("chr1", 1500, ["E11.9"], sex="female")
    batch  = db.query_batch("chr1", [_V1], ["E11.9"], sex="female")
    assert len(single) == len(batch)
    if single:
        assert single[0].AC == batch[0].AC
        assert single[0].AN == batch[0].AN


def test_query_batch_nonexistent_chrom(test_db):
    db = _db(test_db)
    results = db.query_batch("chr99", [(1000, "A", "T")], ["E11.9"])
    assert results == []


# ---------------------------------------------------------------------------
# query_batch — multi-ALT at same position
# ---------------------------------------------------------------------------

def test_query_batch_multi_alt_same_pos(test_db):
    """Two variants at the same position with different ALTs → independent results."""
    db = _db(test_db)
    # Query the known variant plus a fictional one at the same pos; only the
    # real one should come back.
    results = db.query_batch("chr1", [_V1, (1500, "A", "X")], ["E11.9"])
    alts = {r.variant.alt for r in results}
    assert "T" in alts        # the real ALT is returned
    assert "X" not in alts    # fictional ALT is absent


# ---------------------------------------------------------------------------
# query_batch — large batch (temp-table path, ≥ BATCH_IN_THRESHOLD)
# ---------------------------------------------------------------------------

def test_query_batch_large(test_db):
    """≥10 000 unique positions → temp-table path; known variants still returned."""
    db = _db(test_db)
    # Generate a large list; include the real known variants so they are found
    variants = [(i, "A", "T") for i in range(1, 10_001)]
    # Replace entries for known positions with their real (ref, alt)
    variants[1499] = _V1   # pos 1500
    variants[3499] = _V2   # pos 3500
    results = db.query_batch("chr1", variants, ["E11.9"])

    found = {r.variant.pos for r in results}
    assert 1500 in found
    assert 3500 in found


# ---------------------------------------------------------------------------
# query_region — correctness
# ---------------------------------------------------------------------------

def test_query_region_full_range(test_db):
    db = _db(test_db)
    results = db.query_region("chr1", 1, 10_000, ["E11.9"])
    positions = {r.variant.pos for r in results}
    assert 1500 in positions
    assert 3500 in positions
    assert 5000 in positions


def test_query_region_narrow_range(test_db):
    db = _db(test_db)
    results = db.query_region("chr1", 1000, 2000, ["E11.9"])
    positions = {r.variant.pos for r in results}
    assert 1500 in positions
    assert 3500 not in positions
    assert 5000 not in positions


def test_query_region_empty(test_db):
    db = _db(test_db)
    results = db.query_region("chr1", 99_000_000, 99_999_999, ["E11.9"])
    assert results == []


def test_query_region_sorted(test_db):
    db = _db(test_db)
    results = db.query_region("chr1", 1, 10_000, ["E11.9"])
    positions = [r.variant.pos for r in results]
    assert positions == sorted(positions)


def test_query_region_matches_batch(test_db):
    """query_region results should be identical to query_batch over the same variants."""
    db = _db(test_db)
    region = db.query_region("chr1", 1, 10_000, ["E11.9"])
    batch  = db.query_batch("chr1", [_V1, _V2, _V3], ["E11.9"])

    region_idx = {(r.variant.pos, r.variant.alt): r for r in region}
    batch_idx  = {(r.variant.pos, r.variant.alt): r for r in batch}

    assert region_idx.keys() == batch_idx.keys()
    for key in region_idx:
        assert region_idx[key].AC == batch_idx[key].AC
        assert region_idx[key].AN == batch_idx[key].AN


def test_query_region_boundary_inclusive(test_db):
    """start and end should both be inclusive (1-based)."""
    db = _db(test_db)
    # pos=1500 should be included at both start=1500 and end=1500
    r_start = db.query_region("chr1", 1500, 9999, ["E11.9"])
    r_end   = db.query_region("chr1", 1, 1500, ["E11.9"])
    assert any(r.variant.pos == 1500 for r in r_start)
    assert any(r.variant.pos == 1500 for r in r_end)


def test_query_region_all_uncovered(test_db):
    db = _db(test_db)
    results = db.query_region("chr1", 1, 10_000, ["ZZUNKNOWN"])
    assert results == []


def test_query_region_nonexistent_chrom(test_db):
    db = _db(test_db)
    results = db.query_region("chr99", 1, 1_000_000, ["E11.9"])
    assert results == []


# ---------------------------------------------------------------------------
# query_batch_multi — multi-chromosome batch
# ---------------------------------------------------------------------------

def test_query_batch_multi_basic(test_db):
    """Multi-chrom batch returns results from both chromosomes."""
    db = _db(test_db)
    variants = [
        ("chr1", 1500, "A", "T"),
        ("chrX", 5000000, "A", "G"),
    ]
    results = db.query_batch_multi(variants, ["E11.9"])
    chroms = {r.variant.chrom for r in results}
    assert "chr1" in chroms
    assert "chrX" in chroms


def test_query_batch_multi_preserves_order(test_db):
    """Results are returned in input order."""
    db = _db(test_db)
    variants = [
        ("chrX", 5000000, "A", "G"),
        ("chr1", 1500, "A", "T"),
        ("chr1", 3500, "G", "C"),
    ]
    results = db.query_batch_multi(variants)
    assert results[0].variant.chrom == "chrX"
    assert results[1].variant.pos == 1500
    assert results[2].variant.pos == 3500


def test_query_batch_multi_matches_single_queries(test_db):
    """Multi-chrom batch gives same AC/AN as individual queries."""
    db = _db(test_db)
    variants = [
        ("chr1", 1500, "A", "T"),
        ("chrX", 5000000, "A", "G"),
    ]
    results = db.query_batch_multi(variants, ["E11.9"])
    result_map = {(r.variant.chrom, r.variant.pos): r for r in results}

    single_chr1 = db.query("chr1", 1500, ["E11.9"])
    single_chrX = db.query("chrX", 5000000, ["E11.9"])

    assert result_map[("chr1", 1500)].AC == single_chr1[0].AC
    assert result_map[("chr1", 1500)].AN == single_chr1[0].AN
    assert result_map[("chrX", 5000000)].AC == single_chrX[0].AC
    assert result_map[("chrX", 5000000)].AN == single_chrX[0].AN


def test_query_batch_multi_nonexistent_chrom(test_db):
    """Variants on nonexistent chroms are silently skipped."""
    db = _db(test_db)
    variants = [
        ("chr99", 1000, "A", "T"),
        ("chr1", 1500, "A", "T"),
    ]
    results = db.query_batch_multi(variants, ["E11.9"])
    assert all(r.variant.chrom == "chr1" for r in results)


# ---------------------------------------------------------------------------
# CLI parser helpers
# ---------------------------------------------------------------------------

def test_parse_locus_valid():
    from afquery.cli import _parse_locus
    assert _parse_locus("chr1:925952") == ("chr1", 925952)
    assert _parse_locus("1:1500") == ("1", 1500)


def test_parse_locus_invalid():
    import click
    from afquery.cli import _parse_locus
    with pytest.raises(click.UsageError):
        _parse_locus("bad_format")
    with pytest.raises(click.UsageError):
        _parse_locus("chr1:notanint")


def test_parse_region_valid():
    from afquery.cli import _parse_region
    assert _parse_region("chr1:900000-1000000") == ("chr1", 900000, 1000000)
    assert _parse_region("1:1000-2000") == ("1", 1000, 2000)


def test_parse_region_invalid():
    import click
    from afquery.cli import _parse_region
    with pytest.raises(click.UsageError):
        _parse_region("bad_format")
    with pytest.raises(click.UsageError):
        _parse_region("chr1:notarange")


def test_parse_variants_file_4col(tmp_path):
    from afquery.cli import _parse_variants_file
    f = tmp_path / "v.tsv"
    f.write_text("chr1\t1500\tA\tT\nchrX\t5000000\tA\tG\n")
    result = _parse_variants_file(str(f))
    assert result == [("chr1", 1500, "A", "T"), ("chrX", 5000000, "A", "G")]


def test_parse_variants_file_2col(tmp_path):
    from afquery.cli import _parse_variants_file
    f = tmp_path / "v.tsv"
    f.write_text("chr1\t1500\n")
    result = _parse_variants_file(str(f))
    assert result == [("chr1", 1500, "", "")]
