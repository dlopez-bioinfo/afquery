import pytest
from afquery.constants import normalize_chrom, is_autosome, is_sex_chrom, is_mito


def test_normalize_bare_number():
    assert normalize_chrom("1") == "chr1"


def test_normalize_already_prefixed():
    assert normalize_chrom("chr1") == "chr1"


def test_normalize_X():
    assert normalize_chrom("X") == "chrX"


def test_normalize_MT():
    assert normalize_chrom("MT") == "chrM"


def test_normalize_whitespace():
    assert normalize_chrom(" chrX ") == "chrX"


def test_is_autosome_true():
    assert is_autosome("chr1") is True
    assert is_autosome("22") is True


def test_is_autosome_false():
    assert is_autosome("chrX") is False
    assert is_autosome("chrM") is False


def test_is_sex_chrom_true():
    assert is_sex_chrom("chrX") is True
    assert is_sex_chrom("Y") is True


def test_is_sex_chrom_false():
    assert is_sex_chrom("chr1") is False


def test_is_mito_true():
    assert is_mito("chrM") is True
    assert is_mito("MT") is True


def test_is_mito_false():
    assert is_mito("chr1") is False
