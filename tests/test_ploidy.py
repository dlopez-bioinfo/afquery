import pytest
from pyroaring import BitMap
from afquery.ploidy import compute_AN, is_par

# Sample setup: IDs 0-9, males={0,1,4,7,8}, females={2,3,5,6,9}
ALL_10 = BitMap(range(10))
MALES = BitMap([0, 1, 4, 7, 8])
FEMALES = BitMap([2, 3, 5, 6, 9])


def test_is_par_chrX_in_par1_grch37():
    assert is_par("chrX", 100_000, "GRCh37") is True


def test_is_par_chrX_outside_par_grch37():
    assert is_par("chrX", 5_000_000, "GRCh37") is False


def test_is_par_autosome():
    assert is_par("chr1", 1_000_000, "GRCh37") is False


def test_is_par_chrY_in_par1_grch37():
    assert is_par("chrY", 50_000, "GRCh37") is True


def test_is_par_boundary_start():
    # GRCh37 chrX PAR1: 60001-2699520
    assert is_par("chrX", 60_001, "GRCh37") is True


def test_is_par_boundary_end():
    assert is_par("chrX", 2_699_520, "GRCh37") is True


def test_is_par_just_outside_boundary():
    assert is_par("chrX", 60_000, "GRCh37") is False
    assert is_par("chrX", 2_699_521, "GRCh37") is False


# --- compute_AN tests ---

def test_autosome_all_samples():
    AN = compute_AN(ALL_10, MALES, FEMALES, "chr1", 1_000_000, "GRCh37")
    assert AN == 20


def test_autosome_subset():
    subset = BitMap([0, 1, 2, 3, 4])  # 5 samples
    AN = compute_AN(subset, MALES, FEMALES, "chr1", 1_000_000, "GRCh37")
    assert AN == 10


def test_chrX_par_all():
    # PAR1 GRCh37: 60001-2699520
    AN = compute_AN(ALL_10, MALES, FEMALES, "chrX", 100_000, "GRCh37")
    assert AN == 20


def test_chrX_non_par_all():
    # pos 5_000_000 is non-PAR
    # 5 females * 2 + 5 males * 1 = 15
    AN = compute_AN(ALL_10, MALES, FEMALES, "chrX", 5_000_000, "GRCh37")
    assert AN == 15


def test_chrX_non_par_females_only():
    AN = compute_AN(FEMALES, MALES, FEMALES, "chrX", 5_000_000, "GRCh37")
    assert AN == 10  # 5 females * 2


def test_chrX_non_par_males_only():
    AN = compute_AN(MALES, MALES, FEMALES, "chrX", 5_000_000, "GRCh37")
    assert AN == 5  # 5 males * 1


def test_chrY_all():
    # Only males contribute, ploidy=1
    AN = compute_AN(ALL_10, MALES, FEMALES, "chrY", 500_000, "GRCh37")
    assert AN == 5  # 5 males


def test_chrY_females_only():
    AN = compute_AN(FEMALES, MALES, FEMALES, "chrY", 500_000, "GRCh37")
    assert AN == 0


def test_chrM_all():
    AN = compute_AN(ALL_10, MALES, FEMALES, "chrM", 100, "GRCh37")
    assert AN == 10


def test_chrM_subset():
    subset = BitMap([0, 2, 5])
    AN = compute_AN(subset, MALES, FEMALES, "chrM", 100, "GRCh37")
    assert AN == 3


def test_chrom_normalization():
    # "1" should normalize to "chr1"
    AN = compute_AN(ALL_10, MALES, FEMALES, "1", 1_000_000, "GRCh37")
    assert AN == 20


def test_empty_eligible():
    empty = BitMap()
    assert compute_AN(empty, MALES, FEMALES, "chr1", 1_000_000, "GRCh37") == 0
    assert compute_AN(empty, MALES, FEMALES, "chrX", 5_000_000, "GRCh37") == 0
    assert compute_AN(empty, MALES, FEMALES, "chrY", 500_000, "GRCh37") == 0
    assert compute_AN(empty, MALES, FEMALES, "chrM", 100, "GRCh37") == 0


def test_icd10_subset_eligible():
    # Subset: {0,1,2} — 2 males (0,1), 1 female (2)
    subset = BitMap([0, 1, 2])
    AN = compute_AN(subset, MALES, FEMALES, "chrX", 5_000_000, "GRCh37")
    assert AN == 4  # 1 female*2 + 2 males*1
