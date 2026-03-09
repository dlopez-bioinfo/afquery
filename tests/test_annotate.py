"""Tests for annotate_vcf() (Phase 3).

Expected values for icd10=["E11.9"], sex="both", GRCh37:

  chr1:1500 A>T   — in DB    — eligible={0,1,2,5,6}  AN=10  AC=4   AF=0.4
  chr1:2000 A>G   — absent   — eligible={0,1,2,5,6}  AN=10  AC=0   AF=0.0
  chr1:3500 G>C,T — C in DB, T absent
                             — eligible={0,1,2,7}    AN=8   AC=[1,0] AF=[0.125,0.0]
  chr1:5000 T>G   — in DB    — eligible={0,1,2}      AN=6   AC=2   AF=0.333…

Coverage derivation:
  WGS (tech 0, samples 0-3):   covers all positions
  WES_kit_A (tech 1, samples 4-6): chr1 999<pos<=2000  → covers 1500, 2000
  WES_kit_B (tech 2, samples 7-9): chr1 2999<pos<=4000 → covers 3500
  pos 5000 not covered by any WES kit.
"""
import pytest
import cyvcf2


DATA_VCF = "tests/data/annotate_input.vcf"


def _as_tuple(val):
    """Normalise scalar or tuple cyvcf2 INFO return to always be a tuple (or None)."""
    if val is None:
        return None
    if hasattr(val, "__len__"):
        return tuple(val)
    return (val,)


def _read_annotated(path: str) -> dict[int, dict]:
    """Return {pos: {"AC": tuple|None, "AN": int|None, "AF": tuple|None}} from annotated VCF."""
    result = {}
    for v in cyvcf2.VCF(path):
        result[v.POS] = {
            "AC": _as_tuple(v.INFO.get("AFQUERY_AC")),
            "AN": v.INFO.get("AFQUERY_AN"),
            "AF": _as_tuple(v.INFO.get("AFQUERY_AF")),
            "ALT": list(v.ALT),
        }
    return result


# ---------------------------------------------------------------------------
# Round-trip: known variants get correct AFQUERY_AC / AFQUERY_AN / AFQUERY_AF
# ---------------------------------------------------------------------------

def test_annotate_known_variant_1500(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"], sex="both")

    records = _read_annotated(out)
    r = records[1500]
    assert r["AN"] == 10
    assert tuple(r["AC"]) == (4,)
    assert abs(r["AF"][0] - 0.4) < 1e-6


def test_annotate_known_variant_5000(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"], sex="both")

    records = _read_annotated(out)
    r = records[5000]
    assert r["AN"] == 6
    assert tuple(r["AC"]) == (2,)
    assert abs(r["AF"][0] - 2 / 6) < 1e-6


# ---------------------------------------------------------------------------
# Absent variant at covered position → AC=0, AF=0.0, AN>0
# ---------------------------------------------------------------------------

def test_annotate_absent_covered(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"], sex="both")

    records = _read_annotated(out)
    r = records[2000]
    assert r["AN"] == 10
    assert tuple(r["AC"]) == (0,)
    assert r["AF"][0] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Multi-allelic: one ALT known, one absent
# ---------------------------------------------------------------------------

def test_annotate_multiallelic(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"], sex="both")

    records = _read_annotated(out)
    r = records[3500]
    assert r["ALT"] == ["C", "T"]
    assert r["AN"] == 8
    ac = tuple(r["AC"])
    af = tuple(r["AF"])
    assert ac == (1, 0)
    assert af[0] == pytest.approx(0.125)
    assert af[1] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Uncovered position (AN=0) → AFQUERY_AN=0, AFQUERY_AC and AFQUERY_AF absent/missing
# ---------------------------------------------------------------------------

def test_annotate_uncovered_position(test_db, tmp_path):
    """Unknown ICD10 → no eligible samples → AN=0 everywhere."""
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    db.annotate_vcf(DATA_VCF, out, icd10=["ZZUNKNOWN"], sex="both")

    records = _read_annotated(out)
    for pos, r in records.items():
        assert r["AN"] == 0, f"Expected AN=0 at pos {pos}"
        assert r["AC"] is None, f"Expected AFQUERY_AC absent at pos {pos}"
        assert r["AF"] is None, f"Expected AFQUERY_AF absent at pos {pos}"


# ---------------------------------------------------------------------------
# Stats dict returned by annotate_vcf
# ---------------------------------------------------------------------------

def test_annotate_stats(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    stats = db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"], sex="both")

    assert stats["n_variants"] == 4          # 4 records in input VCF
    assert stats["n_uncovered"] == 0         # all covered with E11.9
    # 1500, 3500 (C allele), 5000 all have at least one DB hit
    assert stats["n_annotated"] == 3


def test_annotate_stats_uncovered(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    stats = db.annotate_vcf(DATA_VCF, out, icd10=["ZZUNKNOWN"], sex="both")

    assert stats["n_variants"] == 4
    assert stats["n_uncovered"] == 4
    assert stats["n_annotated"] == 0


# ---------------------------------------------------------------------------
# Output VCF contains exactly the same number of records as input
# ---------------------------------------------------------------------------

def test_annotate_record_count_preserved(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"], sex="both")

    in_count  = sum(1 for _ in cyvcf2.VCF(DATA_VCF))
    out_count = sum(1 for _ in cyvcf2.VCF(out))
    assert in_count == out_count


# ---------------------------------------------------------------------------
# INFO headers are present in output
# ---------------------------------------------------------------------------

def test_annotate_info_headers_present(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"], sex="both")

    header_str = cyvcf2.VCF(out).raw_header
    assert "AFQUERY_AC" in header_str
    assert "AFQUERY_AN" in header_str
    assert "AFQUERY_AF" in header_str


# ---------------------------------------------------------------------------
# Database.annotate_vcf wrapper
# ---------------------------------------------------------------------------

def test_database_annotate_vcf_wrapper(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    stats = db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"])
    assert stats["n_variants"] == 4


# ---------------------------------------------------------------------------
# Sex filter changes eligible set and AC
# ---------------------------------------------------------------------------

def test_annotate_sex_filter_female(test_db, tmp_path):
    from afquery.database import Database
    db = Database(test_db)
    out = str(tmp_path / "out.vcf")
    db.annotate_vcf(DATA_VCF, out, icd10=["E11.9"], sex="female")

    records = _read_annotated(out)
    r = records[1500]
    # female E11.9 samples: {2,5,6}; covered at 1500 by WGS+WES_kit_A → eligible={2,5,6}
    # het={0,5}∩{2,5,6}={5}, hom={2}∩{2,5,6}={2} → AC=1+2=3, AN=6
    assert r["AN"] == 6
    assert tuple(r["AC"]) == (3,)
    assert abs(r["AF"][0] - 0.5) < 1e-6
