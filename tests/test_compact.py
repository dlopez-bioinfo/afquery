import json
import shutil
import sqlite3
from pathlib import Path

import pyarrow.parquet as pq
import pytest
from pyroaring import BitMap

from afquery.bitmaps import deserialize
from afquery.database import Database
from afquery.preprocess.compact import compact_database
from afquery.preprocess.update import remove_samples


# ---------------------------------------------------------------------------
# Fixture: a copy of test_db with sample S07 removed
#   S07 (ID=7) is the only carrier of chr1:3500 (het_ids=[7]).
#   After removal, chr1:3500 has empty bitmaps → compact should drop that row.
# ---------------------------------------------------------------------------

@pytest.fixture
def compact_db(tmp_path, test_db):
    """Copy test_db, remove S07 to create an empty-bitmap row, return db path."""
    db_copy = str(tmp_path / "compact_db")
    shutil.copytree(test_db, db_copy)
    remove_samples(db_copy, ["S07"])
    return db_copy


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_compact_removes_empty_rows(compact_db):
    """Rows where all active samples have 0 genotype should be removed."""
    # Before compact: chr1:3500 exists with empty bitmaps
    tbl_before = pq.read_table(str(Path(compact_db) / "variants" / "chr1.parquet"))
    positions_before = tbl_before["pos"].to_pylist()
    assert 3500 in positions_before, "chr1:3500 should exist before compact"

    compact_database(Path(compact_db))

    tbl_after = pq.read_table(str(Path(compact_db) / "variants" / "chr1.parquet"))
    positions_after = tbl_after["pos"].to_pylist()
    assert 3500 not in positions_after, "chr1:3500 should be removed after compact"


def test_compact_stats(compact_db):
    """compact_database returns expected stats dict with correct keys and types."""
    stats = compact_database(Path(compact_db))

    assert isinstance(stats, dict)
    for key in ("files_rewritten", "rows_processed", "rows_removed", "size_before", "size_after"):
        assert key in stats, f"Missing stat: {key}"
    assert stats["rows_processed"] > 0
    assert stats["rows_removed"] >= 1    # at least the chr1:3500 empty row
    assert stats["files_rewritten"] >= 1


def test_compact_idempotent(compact_db):
    """Running compact twice returns zero rows_removed on second run."""
    compact_database(Path(compact_db))
    stats2 = compact_database(Path(compact_db))
    assert stats2["rows_removed"] == 0
    assert stats2["files_rewritten"] == 0


def test_compact_updates_manifest(compact_db):
    """Compact adds last_compact timestamp to manifest.json."""
    assert "last_compact" not in json.loads(
        (Path(compact_db) / "manifest.json").read_text()
    )
    compact_database(Path(compact_db))
    manifest = json.loads((Path(compact_db) / "manifest.json").read_text())
    assert "last_compact" in manifest
    assert isinstance(manifest["last_compact"], str)


def test_compact_query_results_unchanged(compact_db):
    """Queries return same results before and after compact."""
    db_before = Database(compact_db)
    results_before = db_before.query(
        chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both"
    )

    compact_database(Path(compact_db))

    db_after = Database(compact_db)
    results_after = db_after.query(
        chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both"
    )

    assert len(results_before) == len(results_after)
    for r_b, r_a in zip(results_before, results_after):
        assert r_b.AC == r_a.AC
        assert r_b.AN == r_a.AN
        assert abs(r_b.AF - r_a.AF) < 1e-9


def test_compact_preserves_active_bitmaps(compact_db):
    """After compact, bitmaps only contain active sample IDs."""
    compact_database(Path(compact_db))

    con = sqlite3.connect(str(Path(compact_db) / "metadata.sqlite"))
    active_ids = BitMap(
        [r[0] for r in con.execute("SELECT sample_id FROM samples").fetchall()]
    )
    con.close()

    tbl = pq.read_table(str(Path(compact_db) / "variants" / "chr1.parquet"))
    for i in range(len(tbl)):
        het_bm = deserialize(tbl["het_bitmap"][i].as_py())
        hom_bm = deserialize(tbl["hom_bitmap"][i].as_py())
        assert het_bm.issubset(active_ids), f"Row {i}: het_bm contains inactive sample"
        assert hom_bm.issubset(active_ids), f"Row {i}: hom_bm contains inactive sample"


def test_compact_atomic_write(compact_db, monkeypatch):
    """Compact writes .tmp then renames — tmp file is cleaned up."""
    db_path = Path(compact_db)
    compact_database(db_path)
    # No .tmp files should remain
    tmp_files = list(db_path.rglob("*.parquet.tmp"))
    assert tmp_files == [], f"Leftover tmp files: {tmp_files}"


def test_compact_no_op_when_no_removed_samples(test_db, tmp_path):
    """Compact on a DB with no removed samples rewrites nothing."""
    db_copy = str(tmp_path / "clean_db")
    shutil.copytree(test_db, db_copy)

    stats = compact_database(Path(db_copy))
    assert stats["rows_removed"] == 0
    # May or may not rewrite files depending on bitmap state, but rows_removed=0


def test_database_compact_method(compact_db):
    """Database.compact() delegates to compact_database and reloads."""
    db = Database(compact_db)
    stats = db.compact()
    assert "rows_removed" in stats
    assert stats["rows_removed"] >= 1


def test_compact_partitioned_format(tmp_path, data_dir):
    """Compact also works on partitioned (bucketed) Parquet format."""
    from afquery.preprocess import run_preprocess
    db_path = tmp_path / "partitioned_db"
    db_path.mkdir()

    run_preprocess(
        manifest_path=str(data_dir / "manifest.tsv"),
        output_dir=str(db_path),
        genome_build="GRCh37",
        bed_dir=str(data_dir / "beds"),
        n_threads=2,
    )

    # Remove a sample to create empty-bitmap rows
    remove_samples(str(db_path), ["S07"])

    stats = compact_database(db_path)
    assert stats["rows_processed"] > 0
    # After compact, no .tmp files left
    tmp_files = list(db_path.rglob("*.parquet.tmp"))
    assert tmp_files == []
