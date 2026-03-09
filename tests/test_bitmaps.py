import pytest
from pyroaring import BitMap
from afquery.bitmaps import (
    serialize,
    deserialize,
    build_sex_bitmaps,
    build_icd10_bitmaps,
    build_tech_bitmaps,
)
from afquery.models import Sample


def _make_samples():
    return [
        Sample(0, "S00", "male",   0),
        Sample(1, "S01", "male",   0),
        Sample(2, "S02", "female", 0),
        Sample(3, "S03", "female", 0),
        Sample(4, "S04", "male",   1),
        Sample(5, "S05", "female", 1),
        Sample(6, "S06", "female", 1),
        Sample(7, "S07", "male",   2),
        Sample(8, "S08", "male",   2),
        Sample(9, "S09", "female", 2),
    ]


# --- serialize / deserialize ---

def test_serialize_deserialize_empty():
    bm = BitMap()
    assert deserialize(serialize(bm)) == bm


def test_serialize_deserialize_identity():
    bm = BitMap([0, 2, 5, 100, 999])
    assert deserialize(serialize(bm)) == bm


def test_serialize_returns_bytes():
    bm = BitMap([1, 2, 3])
    assert isinstance(serialize(bm), bytes)


def test_deserialize_from_bytes():
    bm = BitMap([10, 20, 30])
    data = serialize(bm)
    result = deserialize(bytes(data))
    assert result == bm


# --- bitmap operations ---

def test_intersection():
    a = BitMap([0, 1, 2, 3])
    b = BitMap([2, 3, 4, 5])
    assert (a & b) == BitMap([2, 3])


def test_union():
    a = BitMap([0, 1])
    b = BitMap([2, 3])
    assert (a | b) == BitMap([0, 1, 2, 3])


def test_len():
    bm = BitMap([0, 5, 10, 15])
    assert len(bm) == 4


# --- build_sex_bitmaps ---

def test_build_sex_bitmaps_keys():
    samples = _make_samples()
    result = build_sex_bitmaps(samples)
    assert set(result.keys()) == {"male", "female"}


def test_build_sex_bitmaps_males():
    samples = _make_samples()
    result = build_sex_bitmaps(samples)
    assert result["male"] == BitMap([0, 1, 4, 7, 8])


def test_build_sex_bitmaps_females():
    samples = _make_samples()
    result = build_sex_bitmaps(samples)
    assert result["female"] == BitMap([2, 3, 5, 6, 9])


def test_build_sex_bitmaps_partition():
    samples = _make_samples()
    result = build_sex_bitmaps(samples)
    # Union should equal all sample IDs; intersection should be empty
    assert (result["male"] | result["female"]) == BitMap(range(10))
    assert (result["male"] & result["female"]) == BitMap()


# --- build_icd10_bitmaps ---

def test_build_icd10_bitmaps_single_code():
    data = [(0, "E11.9"), (1, "E11.9"), (2, "E11.9")]
    result = build_icd10_bitmaps(data)
    assert result["E11.9"] == BitMap([0, 1, 2])


def test_build_icd10_bitmaps_multiple_codes():
    data = [(0, "E11.9"), (1, "I10"), (0, "I10")]
    result = build_icd10_bitmaps(data)
    assert result["E11.9"] == BitMap([0])
    assert result["I10"] == BitMap([0, 1])


def test_build_icd10_bitmaps_overlapping():
    # Sample 0 appears in two codes
    data = [(0, "A"), (0, "B"), (1, "A")]
    result = build_icd10_bitmaps(data)
    assert result["A"] == BitMap([0, 1])
    assert result["B"] == BitMap([0])


def test_build_icd10_bitmaps_known_fixture():
    sample_icd10 = [
        (0, "E11.9"), (1, "E11.9"), (2, "E11.9"), (5, "E11.9"), (6, "E11.9"), (7, "E11.9"),
        (1, "I10"),   (3, "I10"),   (5, "I10"),   (8, "I10"),   (9, "I10"),
    ]
    result = build_icd10_bitmaps(sample_icd10)
    assert result["E11.9"] == BitMap([0, 1, 2, 5, 6, 7])
    assert result["I10"]   == BitMap([1, 3, 5, 8, 9])


# --- build_tech_bitmaps ---

def test_build_tech_bitmaps_keys():
    samples = _make_samples()
    result = build_tech_bitmaps(samples)
    assert set(result.keys()) == {0, 1, 2}


def test_build_tech_bitmaps_wgs():
    samples = _make_samples()
    result = build_tech_bitmaps(samples)
    assert result[0] == BitMap([0, 1, 2, 3])


def test_build_tech_bitmaps_wes_a():
    samples = _make_samples()
    result = build_tech_bitmaps(samples)
    assert result[1] == BitMap([4, 5, 6])


def test_build_tech_bitmaps_wes_b():
    samples = _make_samples()
    result = build_tech_bitmaps(samples)
    assert result[2] == BitMap([7, 8, 9])
