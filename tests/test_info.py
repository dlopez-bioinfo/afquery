"""Tests for database introspection, audit trail, and migration."""
import json
import sqlite3
from pathlib import Path

import pytest

from afquery.database import Database
from afquery.preprocess.migrate import migrate_sqlite
from afquery.preprocess import _write_sqlite
from afquery.models import Sample, Technology


# ---------------------------------------------------------------------------
# Migration tests
# ---------------------------------------------------------------------------

class TestMigrate:
    def _make_old_db(self, path: Path) -> Path:
        """Create a SQLite DB with the pre-migration 4-column samples schema."""
        db = path / "metadata.sqlite"
        con = sqlite3.connect(str(db))
        con.executescript("""
            CREATE TABLE samples (
                sample_id   INTEGER PRIMARY KEY,
                sample_name TEXT NOT NULL,
                sex         TEXT NOT NULL,
                tech_id     INTEGER NOT NULL
            );
            CREATE TABLE technologies (
                tech_id   INTEGER PRIMARY KEY,
                tech_name TEXT NOT NULL,
                bed_path  TEXT
            );
            CREATE TABLE sample_phenotype (
                sample_id      INTEGER NOT NULL,
                phenotype_code TEXT NOT NULL,
                PRIMARY KEY (sample_id, phenotype_code)
            );
            CREATE TABLE precomputed_bitmaps (
                bitmap_type TEXT NOT NULL,
                bitmap_key  TEXT NOT NULL,
                bitmap_data BLOB NOT NULL,
                PRIMARY KEY (bitmap_type, bitmap_key)
            );
        """)
        con.execute("INSERT INTO samples VALUES (0, 'S0', 'male', 0)")
        con.commit()
        con.close()
        return db

    def test_migrate_adds_columns_to_old_db(self, tmp_path):
        db = self._make_old_db(tmp_path)
        migrate_sqlite(db)
        con = sqlite3.connect(str(db))
        cols = {row[1] for row in con.execute("PRAGMA table_info(samples)").fetchall()}
        con.close()
        assert "vcf_path" in cols
        assert "ingested_at" in cols

    def test_migrate_creates_changelog_table(self, tmp_path):
        db = self._make_old_db(tmp_path)
        migrate_sqlite(db)
        con = sqlite3.connect(str(db))
        tables = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()}
        con.close()
        assert "changelog" in tables

    def test_migrate_idempotent(self, tmp_path):
        db = self._make_old_db(tmp_path)
        migrate_sqlite(db)
        migrate_sqlite(db)  # second call must not raise

    def test_migrate_preserves_existing_rows(self, tmp_path):
        db = self._make_old_db(tmp_path)
        migrate_sqlite(db)
        con = sqlite3.connect(str(db))
        rows = con.execute("SELECT sample_name FROM samples").fetchall()
        con.close()
        assert rows == [("S0",)]


# ---------------------------------------------------------------------------
# info() tests
# ---------------------------------------------------------------------------

EXPECTED_INFO_KEYS = {
    "db_path", "db_version", "genome_build", "schema_version",
    "created_at", "updated_at", "sample_count",
    "by_sex", "by_tech", "by_phenotype", "changelog_recent",
}


class TestInfo:
    def test_info_returns_all_keys(self, test_db):
        db = Database(test_db)
        data = db.info()
        assert EXPECTED_INFO_KEYS == set(data.keys())

    def test_info_sex_breakdown_sums_to_total(self, test_db):
        db = Database(test_db)
        data = db.info()
        assert sum(data["by_sex"].values()) == data["sample_count"]

    def test_info_tech_breakdown_sums_to_total(self, test_db):
        db = Database(test_db)
        data = db.info()
        assert sum(data["by_tech"].values()) == data["sample_count"]

    def test_info_db_version_unknown_for_old_manifest(self, test_db):
        """Fixture manifest has no db_version → shown as 'unknown'."""
        db = Database(test_db)
        data = db.info()
        # conftest manifest has no db_version key
        assert data["db_version"] == "unknown"

    def test_info_sample_count_matches(self, test_db):
        db = Database(test_db)
        data = db.info()
        assert data["sample_count"] == len(db.list_samples())


# ---------------------------------------------------------------------------
# list_samples() tests
# ---------------------------------------------------------------------------

EXPECTED_SAMPLE_FIELDS = {
    "sample_id", "sample_name", "sex", "tech", "phenotypes", "vcf_path", "ingested_at"
}


class TestListSamples:
    def test_list_samples_fields(self, test_db):
        db = Database(test_db)
        samples = db.list_samples()
        assert len(samples) > 0
        for s in samples:
            assert set(s.keys()) == EXPECTED_SAMPLE_FIELDS

    def test_list_samples_vcf_null_for_fixture(self, test_db):
        """Fixture inserts samples without vcf_path → should be None."""
        db = Database(test_db)
        samples = db.list_samples()
        for s in samples:
            assert s["vcf_path"] is None

    def test_list_samples_ingested_at_null_for_fixture(self, test_db):
        """Fixture inserts samples without ingested_at → should be None."""
        db = Database(test_db)
        samples = db.list_samples()
        for s in samples:
            assert s["ingested_at"] is None

    def test_list_samples_ordered_by_id(self, test_db):
        db = Database(test_db)
        samples = db.list_samples()
        ids = [s["sample_id"] for s in samples]
        assert ids == sorted(ids)


# ---------------------------------------------------------------------------
# changelog() tests
# ---------------------------------------------------------------------------

class TestChangelog:
    def test_changelog_empty_on_fixture_db(self, test_db):
        """Fixture builds DB without calling run_preprocess → no changelog entries."""
        db = Database(test_db)
        assert db.changelog() == []

    def test_changelog_after_write_sqlite(self, tmp_path):
        """Calling _write_sqlite inserts a preprocess changelog entry."""
        (tmp_path / "variants").mkdir()
        (tmp_path / "capture").mkdir()
        samples = [Sample(0, "S0", "male", 0)]
        technologies = [Technology(0, "WGS", None)]
        vcf_paths = ["/path/to/s0.vcf.gz"]
        ingested_at = "2026-01-01T00:00:00+00:00"
        _write_sqlite(str(tmp_path), samples, technologies, [], vcf_paths, ingested_at)

        con = sqlite3.connect(str(tmp_path / "metadata.sqlite"))
        rows = con.execute(
            "SELECT event_type, notes FROM changelog"
        ).fetchall()
        con.close()
        assert len(rows) == 1
        assert rows[0][0] == "preprocess"
        assert "1 samples ingested" in rows[0][1]

    def test_changelog_after_write_sqlite_has_sample_names(self, tmp_path):
        (tmp_path / "variants").mkdir()
        (tmp_path / "capture").mkdir()
        samples = [Sample(0, "Alpha", "male", 0), Sample(1, "Beta", "female", 0)]
        technologies = [Technology(0, "WGS", None)]
        _write_sqlite(str(tmp_path), samples, technologies, [], [None, None], None)

        con = sqlite3.connect(str(tmp_path / "metadata.sqlite"))
        row = con.execute("SELECT sample_names FROM changelog").fetchone()
        con.close()
        names = json.loads(row[0])
        assert "Alpha" in names
        assert "Beta" in names


# ---------------------------------------------------------------------------
# set_db_version() tests
# ---------------------------------------------------------------------------

class TestSetDbVersion:
    def test_set_db_version(self, test_db):
        db = Database(test_db)
        db.set_db_version("2.5")
        db2 = Database(test_db)
        assert db2.info()["db_version"] == "2.5"

    def test_set_db_version_writes_changelog(self, test_db):
        """set_db_version appends a changelog entry with event_type='version_set'."""
        db = Database(test_db)
        db.set_db_version("9.9")
        con = sqlite3.connect(str(Path(test_db) / "metadata.sqlite"))
        rows = con.execute(
            "SELECT event_type, notes FROM changelog WHERE event_type='version_set'"
            " ORDER BY event_id DESC LIMIT 1"
        ).fetchall()
        con.close()
        assert len(rows) == 1
        assert rows[0][0] == "version_set"
        assert "9.9" in rows[0][1]

    def test_set_db_version_updates_manifest_file(self, test_db):
        db = Database(test_db)
        db.set_db_version("99.0")
        manifest = json.loads((Path(test_db) / "manifest.json").read_text())
        assert manifest["db_version"] == "99.0"
        # Reset for other tests
        db.set_db_version("unknown_reset")


# ---------------------------------------------------------------------------
# _bump_version tests
# ---------------------------------------------------------------------------

class TestBumpVersion:
    def test_bump_integer(self):
        from afquery.preprocess.update import _bump_version
        assert _bump_version("1") == "2"

    def test_bump_semver(self):
        from afquery.preprocess.update import _bump_version
        assert _bump_version("1.0") == "1.1"
        assert _bump_version("1.2.3") == "1.2.4"

    def test_bump_non_numeric_suffix(self):
        from afquery.preprocess.update import _bump_version
        assert _bump_version("v2-alpha") == "v2-alpha.1"
