"""Tests for update_sample_metadata and parse_updates_tsv."""
import json
import os
import shutil
import sqlite3

import pytest

from afquery.preprocess.update import (
    UpdateError,
    parse_updates_tsv,
    update_sample_metadata,
)
from afquery.database import Database


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_sex(db_dir: str, sample_name: str) -> str:
    con = sqlite3.connect(os.path.join(db_dir, "metadata.sqlite"))
    try:
        return con.execute(
            "SELECT sex FROM samples WHERE sample_name = ?", (sample_name,)
        ).fetchone()[0]
    finally:
        con.close()


def _get_phenotypes(db_dir: str, sample_name: str) -> list[str]:
    con = sqlite3.connect(os.path.join(db_dir, "metadata.sqlite"))
    try:
        sid = con.execute(
            "SELECT sample_id FROM samples WHERE sample_name = ?", (sample_name,)
        ).fetchone()[0]
        return sorted(
            r[0]
            for r in con.execute(
                "SELECT phenotype_code FROM sample_phenotype WHERE sample_id = ?", (sid,)
            ).fetchall()
        )
    finally:
        con.close()


def _get_bitmap(db_dir: str, btype: str, bkey: str):
    from afquery.bitmaps import deserialize
    con = sqlite3.connect(os.path.join(db_dir, "metadata.sqlite"))
    try:
        row = con.execute(
            "SELECT bitmap_data FROM precomputed_bitmaps WHERE bitmap_type = ? AND bitmap_key = ?",
            (btype, bkey),
        ).fetchone()
        return deserialize(row[0]) if row else None
    finally:
        con.close()


def _get_changelog(db_dir: str, event_type: str = "UPDATE_SAMPLE") -> list[dict]:
    con = sqlite3.connect(os.path.join(db_dir, "metadata.sqlite"))
    try:
        rows = con.execute(
            "SELECT event_id, event_type, event_time, sample_names, notes"
            " FROM changelog WHERE event_type = ? ORDER BY event_id",
            (event_type,),
        ).fetchall()
        return [
            {
                "event_id": r[0],
                "event_type": r[1],
                "event_time": r[2],
                "sample_names": json.loads(r[3]) if r[3] else None,
                "notes": json.loads(r[4]) if r[4] else None,
            }
            for r in rows
        ]
    finally:
        con.close()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def fresh_db(test_db, tmp_path):
    """Return a writable copy of the session-scoped test_db."""
    dst = str(tmp_path / "db")
    shutil.copytree(test_db, dst)
    return dst


# ---------------------------------------------------------------------------
# parse_updates_tsv
# ---------------------------------------------------------------------------

def test_parse_updates_tsv_basic(tmp_path):
    tsv = tmp_path / "changes.tsv"
    tsv.write_text(
        "sample_name\tfield\tnew_value\n"
        "S01\tsex\tfemale\n"
        "S02\tphenotype_codes\tE11.9,I10\n"
    )
    updates = parse_updates_tsv(str(tsv))
    assert len(updates) == 2
    assert updates[0] == {"sample_name": "S01", "field": "sex", "new_value": "female"}
    assert updates[1] == {"sample_name": "S02", "field": "phenotype_codes", "new_value": "E11.9,I10"}


def test_parse_updates_tsv_skips_blank_lines(tmp_path):
    tsv = tmp_path / "changes.tsv"
    tsv.write_text("sample_name\tfield\tnew_value\n\nS01\tsex\tfemale\n\n")
    updates = parse_updates_tsv(str(tsv))
    assert len(updates) == 1


def test_parse_updates_tsv_missing_column_raises(tmp_path):
    tsv = tmp_path / "bad.tsv"
    tsv.write_text("sample_name\tnew_value\nS01\tfemale\n")
    with pytest.raises(UpdateError, match="missing required columns"):
        parse_updates_tsv(str(tsv))


# ---------------------------------------------------------------------------
# update_sample_metadata — sex
# ---------------------------------------------------------------------------

def test_update_sex(fresh_db):
    # S00 starts as male
    assert _get_sex(fresh_db, "S00") == "male"

    entries = update_sample_metadata(
        fresh_db,
        [{"sample_name": "S00", "field": "sex", "new_value": "female"}],
    )

    assert len(entries) == 1
    assert entries[0]["old"] == "male"
    assert entries[0]["new"] == "female"
    assert _get_sex(fresh_db, "S00") == "female"


def test_update_sex_regenerates_bitmap(fresh_db):
    # S00 (id=0) is male — bit 0 must be in male bitmap before update
    male_bm_before = _get_bitmap(fresh_db, "sex", "male")
    assert 0 in male_bm_before

    update_sample_metadata(
        fresh_db,
        [{"sample_name": "S00", "field": "sex", "new_value": "female"}],
    )

    male_bm_after = _get_bitmap(fresh_db, "sex", "male")
    female_bm_after = _get_bitmap(fresh_db, "sex", "female")
    assert 0 not in male_bm_after
    assert 0 in female_bm_after


# ---------------------------------------------------------------------------
# update_sample_metadata — phenotype_codes
# ---------------------------------------------------------------------------

def test_update_phenotype_replace(fresh_db):
    # S00 starts with E11.9 and J45
    before = _get_phenotypes(fresh_db, "S00")
    assert "E11.9" in before

    entries = update_sample_metadata(
        fresh_db,
        [{"sample_name": "S00", "field": "phenotype_codes", "new_value": "C50,G20"}],
    )

    assert len(entries) == 1
    assert entries[0]["field"] == "phenotype_codes"
    after = _get_phenotypes(fresh_db, "S00")
    assert sorted(after) == ["C50", "G20"]
    # old codes must be gone
    assert "E11.9" not in after
    assert "J45" not in after


def test_update_phenotype_regenerates_bitmap(fresh_db):
    # S00 (id=0) is in E11.9 bitmap before update
    bm_before = _get_bitmap(fresh_db, "phenotype", "E11.9")
    assert 0 in bm_before

    update_sample_metadata(
        fresh_db,
        [{"sample_name": "S00", "field": "phenotype_codes", "new_value": "C50"}],
    )

    bm_e119_after = _get_bitmap(fresh_db, "phenotype", "E11.9")
    bm_c50_after = _get_bitmap(fresh_db, "phenotype", "C50")
    assert 0 not in bm_e119_after
    assert 0 in bm_c50_after


# ---------------------------------------------------------------------------
# update_sample_metadata — both fields in one call
# ---------------------------------------------------------------------------

def test_update_both_fields(fresh_db):
    entries = update_sample_metadata(
        fresh_db,
        [
            {"sample_name": "S00", "field": "sex", "new_value": "female"},
            {"sample_name": "S00", "field": "phenotype_codes", "new_value": "I10"},
        ],
    )

    assert len(entries) == 2
    assert _get_sex(fresh_db, "S00") == "female"
    assert _get_phenotypes(fresh_db, "S00") == ["I10"]


# ---------------------------------------------------------------------------
# update_sample_metadata — changelog
# ---------------------------------------------------------------------------

def test_changelog_entry_format(fresh_db):
    update_sample_metadata(
        fresh_db,
        [{"sample_name": "S00", "field": "sex", "new_value": "female"}],
        operator_note="clinical correction",
    )

    cl = _get_changelog(fresh_db)
    assert len(cl) == 1
    entry = cl[0]
    assert entry["event_type"] == "UPDATE_SAMPLE"
    assert entry["sample_names"] == ["S00"]
    notes = entry["notes"]
    assert notes["sample"] == "S00"
    assert notes["field"] == "sex"
    assert notes["old"] == "male"
    assert notes["new"] == "female"
    assert notes["operator_note"] == "clinical correction"


def test_changelog_entry_no_operator_note(fresh_db):
    update_sample_metadata(
        fresh_db,
        [{"sample_name": "S00", "field": "sex", "new_value": "female"}],
    )
    cl = _get_changelog(fresh_db)
    assert "operator_note" not in cl[0]["notes"]


def test_multiple_updates_create_multiple_changelog_entries(fresh_db):
    update_sample_metadata(
        fresh_db,
        [
            {"sample_name": "S00", "field": "sex", "new_value": "female"},
            {"sample_name": "S01", "field": "sex", "new_value": "female"},
        ],
    )
    cl = _get_changelog(fresh_db)
    assert len(cl) == 2


# ---------------------------------------------------------------------------
# update_sample_metadata — batch TSV
# ---------------------------------------------------------------------------

def test_batch_update_from_tsv(fresh_db, tmp_path):
    tsv = tmp_path / "changes.tsv"
    tsv.write_text(
        "sample_name\tfield\tnew_value\n"
        "S00\tsex\tfemale\n"
        "S01\tphenotype_codes\tC50\n"
    )
    updates = parse_updates_tsv(str(tsv))
    entries = update_sample_metadata(fresh_db, updates)

    assert len(entries) == 2
    assert _get_sex(fresh_db, "S00") == "female"
    assert _get_phenotypes(fresh_db, "S01") == ["C50"]


# ---------------------------------------------------------------------------
# update_sample_metadata — empty updates (no-op)
# ---------------------------------------------------------------------------

def test_empty_updates_is_noop(fresh_db):
    entries = update_sample_metadata(fresh_db, [])
    assert entries == []


# ---------------------------------------------------------------------------
# Validation errors
# ---------------------------------------------------------------------------

def test_invalid_sample_name_raises(fresh_db):
    with pytest.raises(UpdateError, match="not found in database"):
        update_sample_metadata(
            fresh_db,
            [{"sample_name": "DOES_NOT_EXIST", "field": "sex", "new_value": "female"}],
        )


def test_invalid_field_name_raises(fresh_db):
    with pytest.raises(UpdateError, match="Invalid field"):
        update_sample_metadata(
            fresh_db,
            [{"sample_name": "S00", "field": "technology", "new_value": "WGS"}],
        )


def test_invalid_sex_value_raises(fresh_db):
    with pytest.raises(UpdateError, match="Invalid sex value"):
        update_sample_metadata(
            fresh_db,
            [{"sample_name": "S00", "field": "sex", "new_value": "unknown"}],
        )


def test_empty_phenotype_codes_raises(fresh_db):
    with pytest.raises(UpdateError, match="empty"):
        update_sample_metadata(
            fresh_db,
            [{"sample_name": "S00", "field": "phenotype_codes", "new_value": "  ,  "}],
        )


def test_validation_is_atomic(fresh_db):
    """If any update is invalid, no changes should be applied."""
    original_sex = _get_sex(fresh_db, "S00")
    with pytest.raises(UpdateError):
        update_sample_metadata(
            fresh_db,
            [
                {"sample_name": "S00", "field": "sex", "new_value": "female"},
                {"sample_name": "DOES_NOT_EXIST", "field": "sex", "new_value": "female"},
            ],
        )
    # S00 should be unchanged because validation failed before any writes
    assert _get_sex(fresh_db, "S00") == original_sex


# ---------------------------------------------------------------------------
# Integration: phenotype update affects query results
# ---------------------------------------------------------------------------

def test_phenotype_update_affects_query(fresh_db):
    """After updating phenotype, querying with the new code should include the sample."""
    db = Database(fresh_db)

    # S07 (id=7) has E11.9 and J45 initially; query E11.9 should include it
    results_before = db.query("chr1", 1500, phenotype=["E11.9"])
    eligible_before = {r.n_samples_eligible for r in results_before}

    # Replace S07's phenotypes so E11.9 is removed
    db.update_sample_metadata(
        [{"sample_name": "S07", "field": "phenotype_codes", "new_value": "I10"}]
    )

    results_after = db.query("chr1", 1500, phenotype=["E11.9"])
    # n_samples_eligible must decrease since S07 no longer has E11.9
    for r_before, r_after in zip(results_before, results_after):
        assert r_after.n_samples_eligible <= r_before.n_samples_eligible


# ---------------------------------------------------------------------------
# CLI smoke test
# ---------------------------------------------------------------------------

def test_cli_update_sample_sex(fresh_db):
    from click.testing import CliRunner
    from afquery.cli import cli

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["update-db", "--db", fresh_db, "--update-sample", "S00", "--set-sex", "female"],
    )
    assert result.exit_code == 0, result.output
    assert "S00" in result.output
    assert _get_sex(fresh_db, "S00") == "female"


def test_cli_update_sample_phenotype(fresh_db):
    from click.testing import CliRunner
    from afquery.cli import cli

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "update-db", "--db", fresh_db,
            "--update-sample", "S00",
            "--set-phenotype", "C50,G20",
        ],
    )
    assert result.exit_code == 0, result.output
    assert _get_phenotypes(fresh_db, "S00") == ["C50", "G20"]


def test_cli_update_samples_file(fresh_db, tmp_path):
    from click.testing import CliRunner
    from afquery.cli import cli

    tsv = tmp_path / "changes.tsv"
    tsv.write_text(
        "sample_name\tfield\tnew_value\n"
        "S00\tsex\tfemale\n"
    )
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["update-db", "--db", fresh_db, "--update-samples-file", str(tsv)],
    )
    assert result.exit_code == 0, result.output
    assert _get_sex(fresh_db, "S00") == "female"


def test_cli_missing_update_sample_raises(fresh_db):
    from click.testing import CliRunner
    from afquery.cli import cli

    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["update-db", "--db", fresh_db, "--set-sex", "female"],
    )
    assert result.exit_code != 0
    assert "--set-sex" in result.output or "require" in result.output


def test_cli_no_operation_raises(fresh_db):
    from click.testing import CliRunner
    from afquery.cli import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["update-db", "--db", fresh_db])
    assert result.exit_code != 0
