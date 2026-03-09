import json
from pathlib import Path

import pytest

from afquery import Database

DATA_DIR = Path(__file__).parent / "data"


def load_expected():
    return json.loads((DATA_DIR / "expected_results.json").read_text())


@pytest.mark.parametrize("case", load_expected())
def test_query_expected(test_db, case):
    db = Database(test_db)
    results = db.query(
        chrom=case["chrom"],
        pos=case["pos"],
        icd10=case["icd10"],
        sex=case["sex"],
    )

    # Find the result matching the expected alt
    match = next((r for r in results if r.variant.alt == case["alt"]), None)
    assert match is not None, (
        f"No result for alt={case['alt']} at {case['chrom']}:{case['pos']}"
    )
    assert match.AC == case["AC"], f"{case['description']}: AC mismatch"
    assert match.AN == case["AN"], f"{case['description']}: AN mismatch"
    assert match.n_samples_eligible == case["n_eligible"], \
        f"{case['description']}: n_eligible mismatch"

    if case["AF"] is None:
        assert match.AF is None
    else:
        assert abs(match.AF - case["AF"]) < 1e-9, f"{case['description']}: AF mismatch"


def test_query_an_zero_returns_empty(test_db):
    # C50={2,4,6}, male filter → only sample 4 (male, WES_A tech_id=1)
    # pos=3500 is in WES_B only (3000-4000) → sample 4 NOT covered → eligible={}
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=3500, icd10=["C50"], sex="male")
    assert results == []


def test_query_no_variant_at_position(test_db):
    # Position that is covered but has no variant in the Parquet
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1000, icd10=["E11.9"], sex="both")
    assert results == []


def test_query_unknown_icd10(test_db):
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1500, icd10=["UNKNOWN"], sex="both")
    assert results == []


def test_query_chrom_normalization(test_db):
    # "1" should work the same as "chr1"
    db = Database(test_db)
    r1 = db.query(chrom="chr1", pos=1500, icd10=["E11.9"], sex="both")
    r2 = db.query(chrom="1",    pos=1500, icd10=["E11.9"], sex="both")
    assert len(r1) == len(r2)
    assert r1[0].AC == r2[0].AC


def test_query_multiple_icd10_union(test_db):
    # E11.9 + I10 should give more eligible than E11.9 alone (larger union)
    db = Database(test_db)
    r_e119 = db.query(chrom="chr1", pos=1500, icd10=["E11.9"],        sex="both")
    r_both = db.query(chrom="chr1", pos=1500, icd10=["E11.9", "I10"], sex="both")
    assert r_both[0].n_samples_eligible >= r_e119[0].n_samples_eligible


def test_query_nonexistent_chrom(test_db):
    db = Database(test_db)
    results = db.query(chrom="chr22", pos=1000, icd10=["E11.9"], sex="both")
    assert results == []


def test_info(test_db):
    db = Database(test_db)
    data = db.info()
    assert data["genome_build"] == "GRCh37"
    assert data["sample_count"] == 10
