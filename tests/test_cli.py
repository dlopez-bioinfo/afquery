import json

import pytest
from click.testing import CliRunner

from afquery.cli import cli


@pytest.fixture
def runner():
    return CliRunner()


# --- afquery info ---

def test_info_shows_genome_build(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db])
    assert result.exit_code == 0
    assert "GRCh37" in result.output


def test_info_shows_sample_count(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db])
    assert result.exit_code == 0
    assert "10" in result.output


def test_info_missing_db(runner):
    result = runner.invoke(cli, ["info"])
    assert result.exit_code != 0


# --- afquery query: text format ---

def test_query_text_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "1", "--pos", "1500", "--icd10", "E11.9",
    ])
    assert result.exit_code == 0
    assert "AC=4" in result.output
    assert "AN=10" in result.output
    assert "AF=0.4000" in result.output


def test_query_text_no_results(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "1", "--pos", "9999", "--icd10", "UNKNOWN",
    ])
    assert result.exit_code == 0
    assert "No variants found" in result.output


# --- afquery query: json format ---

def test_query_json_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "chr1", "--pos", "1500", "--icd10", "E11.9",
        "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) == 1
    assert data[0]["AC"] == 4
    assert data[0]["AN"] == 10
    assert abs(data[0]["AF"] - 0.4) < 1e-9


def test_query_json_empty(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "chr1", "--pos", "9999", "--icd10", "UNKNOWN",
        "--format", "json",
    ])
    assert result.exit_code == 0
    assert json.loads(result.output) == []


# --- afquery query: tsv format ---

def test_query_tsv_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "1", "--pos", "1500", "--icd10", "E11.9",
        "--format", "tsv",
    ])
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert lines[0].startswith("chrom\tpos\tref\talt")
    assert len(lines) == 2  # header + 1 variant
    fields = lines[1].split("\t")
    assert fields[4] == "4"   # AC
    assert fields[5] == "10"  # AN


# --- sex filter ---

def test_query_sex_female(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "chr1", "--pos", "1500", "--icd10", "E11.9",
        "--sex", "female", "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data[0]["AC"] == 3
    assert data[0]["AN"] == 6


# --- multiple ICD10 codes ---

def test_query_multiple_icd10(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "chr1", "--pos", "1500",
        "--icd10", "E11.9", "--icd10", "I10",
        "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) == 1


# --- missing required option ---

def test_query_missing_chrom(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--pos", "1500", "--icd10", "E11.9",
    ])
    assert result.exit_code != 0


def test_query_missing_db(runner):
    result = runner.invoke(cli, [
        "query", "--chrom", "1", "--pos", "1500", "--icd10", "E11.9",
    ])
    assert result.exit_code != 0
