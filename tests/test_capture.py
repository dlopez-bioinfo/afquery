import pickle
import tempfile
from pathlib import Path

import pytest

from afquery.capture import CaptureIndex

DATA_DIR = Path(__file__).parent / "data"
BED_A = str(DATA_DIR / "beds" / "wes_kit_a.bed")  # chr1:999-2000, chrX:999-2000
BED_B = str(DATA_DIR / "beds" / "wes_kit_b.bed")  # chr1:2999-4000


# --- WGS sentinel ---

def test_wgs_always_covered():
    idx = CaptureIndex.wgs()
    assert idx.covers("chr1", 1) is True
    assert idx.covers("chrX", 5_000_000) is True
    assert idx.covers("chrY", 1) is True
    assert idx.covers("chrM", 100) is True


# --- from_bed: chr1:999-2000 (0-based half-open → 1-based: 1000-2000) ---

def test_inside_region():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 1500) is True


def test_outside_region_before():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 999) is False


def test_outside_region_after():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 2001) is False


def test_exact_start_boundary():
    # BED Start=999 (0-based) → first covered 1-based pos = 1000
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 1000) is True
    assert idx.covers("chr1", 999) is False


def test_exact_end_boundary():
    # BED End=2000 → last covered 1-based pos = 2000
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 2000) is True
    assert idx.covers("chr1", 2001) is False


def test_wrong_chrom():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr2", 1500) is False


def test_chrX_covered():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chrX", 1500) is True


def test_chrX_not_covered_by_bed_b():
    idx = CaptureIndex.from_bed(BED_B)
    assert idx.covers("chrX", 1500) is False


def test_bed_b_covers_3500():
    idx = CaptureIndex.from_bed(BED_B)
    assert idx.covers("chr1", 3500) is True


def test_bed_b_does_not_cover_1500():
    idx = CaptureIndex.from_bed(BED_B)
    assert idx.covers("chr1", 1500) is False


# --- pickle save/load round-trip ---

def test_save_load_wgs(tmp_path):
    idx = CaptureIndex.wgs()
    path = str(tmp_path / "cap.pickle")
    idx.save(path)
    loaded = CaptureIndex.load(path)
    assert loaded.covers("chr1", 9999) is True


def test_save_load_bed(tmp_path):
    idx = CaptureIndex.from_bed(BED_A)
    path = str(tmp_path / "cap.pickle")
    idx.save(path)
    loaded = CaptureIndex.load(path)
    assert loaded.covers("chr1", 1500) is True
    assert loaded.covers("chr1", 999) is False
