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
        "--chrom", "1", "--pos", "1500", "--phenotype", "E11.9",
    ])
    assert result.exit_code == 0
    assert "AC=4" in result.output
    assert "AN=10" in result.output
    assert "AF=0.4000" in result.output


def test_query_text_no_results(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "1", "--pos", "9999", "--phenotype", "UNKNOWN",
    ])
    assert result.exit_code == 0
    assert "No variants found" in result.output


# --- afquery query: json format ---

def test_query_json_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "chr1", "--pos", "1500", "--phenotype", "E11.9",
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
        "--chrom", "chr1", "--pos", "9999", "--phenotype", "UNKNOWN",
        "--format", "json",
    ])
    assert result.exit_code == 0
    assert json.loads(result.output) == []


# --- afquery query: tsv format ---

def test_query_tsv_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "1", "--pos", "1500", "--phenotype", "E11.9",
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
        "--chrom", "chr1", "--pos", "1500", "--phenotype", "E11.9",
        "--sex", "female", "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data[0]["AC"] == 3
    assert data[0]["AN"] == 6


# --- multiple Phenotype codes ---

def test_query_multiple_phenotype(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--chrom", "chr1", "--pos", "1500",
        "--phenotype", "E11.9", "--phenotype", "I10",
        "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) == 1


# --- missing required option ---

def test_query_missing_chrom(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--pos", "1500", "--phenotype", "E11.9",
    ])
    assert result.exit_code != 0


def test_query_missing_db(runner):
    result = runner.invoke(cli, [
        "query", "--chrom", "1", "--pos", "1500", "--phenotype", "E11.9",
    ])
    assert result.exit_code != 0


# --- afquery annotate ---

def test_annotate_with_phenotype(runner, test_db, tmp_path):
    """Test annotate command with explicit phenotype."""
    import cyvcf2
    from pathlib import Path

    # Create a minimal input VCF
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"

    with open(input_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("1\t1500\t.\tC\tT\t.\t.\t.\n")

    result = runner.invoke(cli, [
        "annotate", "--db", test_db,
        "--input", str(input_vcf),
        "--output", str(output_vcf),
        "--phenotype", "E11.9",
    ])
    assert result.exit_code == 0
    assert "Annotated" in result.output
    assert output_vcf.exists()


def test_annotate_without_phenotype(runner, test_db, tmp_path):
    """Test annotate command without phenotype (should use all phenotypes)."""
    import cyvcf2
    from pathlib import Path

    # Create a minimal input VCF
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"

    with open(input_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("1\t1500\t.\tC\tT\t.\t.\t.\n")

    result = runner.invoke(cli, [
        "annotate", "--db", test_db,
        "--input", str(input_vcf),
        "--output", str(output_vcf),
    ])
    assert result.exit_code == 0
    assert "Annotated" in result.output
    assert output_vcf.exists()

    # Verify that the output VCF has AFQUERY annotations
    vcf = cyvcf2.VCF(str(output_vcf))
    records = list(vcf)
    assert len(records) > 0
    # At pos 1500 with multiple phenotypes, AN should be larger than with single phenotype
    assert records[0].INFO.get("AFQUERY_AN", 0) > 0
    vcf.close()
