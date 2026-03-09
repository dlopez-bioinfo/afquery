"""Tests for QueryEngine.query_batch() and query_region() (Phase 3)."""
import pytest
from afquery.database import Database


def _db(test_db):
    return Database(test_db)


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
    batch = db.query_batch("chr1", [1500, 3500, 5000], ["E11.9"])

    single_idx = {(r.variant.pos, r.variant.alt): r for r in single}
    batch_idx  = {(r.variant.pos, r.variant.alt): r for r in batch}

    assert set(single_idx) == set(batch_idx)
    for key in single_idx:
        assert single_idx[key].AC == batch_idx[key].AC
        assert single_idx[key].AN == batch_idx[key].AN
        assert abs(single_idx[key].AF - batch_idx[key].AF) < 1e-9


def test_query_batch_deduplicates_positions(test_db):
    db = _db(test_db)
    results = db.query_batch("chr1", [1500, 1500, 3500], ["E11.9"])
    positions = [r.variant.pos for r in results]
    assert positions.count(1500) == 1


def test_query_batch_sorted_result(test_db):
    db = _db(test_db)
    results = db.query_batch("chr1", [5000, 1500, 3500], ["E11.9"])
    positions = [r.variant.pos for r in results]
    assert positions == sorted(positions)


def test_query_batch_empty_no_matching_position(test_db):
    db = _db(test_db)
    results = db.query_batch("chr1", [99999999], ["E11.9"])
    assert results == []


def test_query_batch_all_uncovered_icd10(test_db):
    """Unknown ICD10 → eligible bitmap empty → AN=0 → no results."""
    db = _db(test_db)
    results = db.query_batch("chr1", [1500, 3500], ["ZZUNKNOWN"])
    assert results == []


def test_query_batch_chrom_normalisation(test_db):
    """chr-prefix should be normalised transparently."""
    db = _db(test_db)
    r_prefix = db.query_batch("chr1", [1500], ["E11.9"])
    r_bare   = db.query_batch("1",    [1500], ["E11.9"])
    assert len(r_prefix) == len(r_bare)
    if r_prefix:
        assert r_prefix[0].AC == r_bare[0].AC


def test_query_batch_sex_filter(test_db):
    """Female-only batch should give same result as single query with female filter."""
    db = _db(test_db)
    single = db.query("chr1", 1500, ["E11.9"], sex="female")
    batch  = db.query_batch("chr1", [1500], ["E11.9"], sex="female")
    assert len(single) == len(batch)
    if single:
        assert single[0].AC == batch[0].AC
        assert single[0].AN == batch[0].AN


def test_query_batch_nonexistent_chrom(test_db):
    db = _db(test_db)
    results = db.query_batch("chr99", [1000], ["E11.9"])
    assert results == []


# ---------------------------------------------------------------------------
# query_batch — large batch (temp-table path, ≥ BATCH_IN_THRESHOLD)
# ---------------------------------------------------------------------------

def test_query_batch_large(test_db):
    """≥10 000 unique positions → temp-table path; known variants still returned."""
    db = _db(test_db)
    positions = list(range(1, 10_001))   # 10 000 positions, includes 1500 & 3500
    results = db.query_batch("chr1", positions, ["E11.9"])

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
    """query_region results should be identical to query_batch over the same positions."""
    db = _db(test_db)
    region = db.query_region("chr1", 1, 10_000, ["E11.9"])
    batch  = db.query_batch("chr1", [1500, 3500, 5000], ["E11.9"])

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
