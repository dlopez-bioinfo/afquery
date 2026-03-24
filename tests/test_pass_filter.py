"""Tests for FILTER=PASS ingestion (Task 10) and N_FAIL tracking (Task 11)."""
import json
import sqlite3
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from pyroaring import BitMap

from afquery.bitmaps import deserialize, serialize
from afquery.database import Database
from afquery.models import Sample, Technology
from afquery.preprocess import run_preprocess

from conftest import VARIANTS, SAMPLES, SAMPLE_PHENOTYPE, TECHNOLOGIES


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_vcf_with_filter(path: Path, sample_name: str, records: list[tuple]) -> None:
    """Write a minimal VCF. records = [(chrom, pos, ref, alt, gt_str, filter_str), ...]"""
    contigs = sorted({r[0] for r in records}) if records else ["chr1"]
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FILTER=<ID=LowQual,Description="Low quality">\n')
        for contig in contigs:
            f.write(f"##contig=<ID={contig}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )
        for chrom, pos, ref, alt, gt, flt in records:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\tGT\t{gt}\n")


# ---------------------------------------------------------------------------
# Tests: N_FAIL from test_db
# ---------------------------------------------------------------------------

def test_fail_samples_returned(test_db):
    """N_FAIL should always be present and match conftest VARIANTS."""
    db = Database(test_db)
    # chr1:1500 A>T — fail_ids=[0] from conftest
    results = db.query(chrom="chr1", pos=1500)
    assert results, "Expected at least one result"
    r = next((x for x in results if x.variant.ref == "A" and x.variant.alt == "T"), None)
    assert r is not None
    assert r.N_FAIL >= 0, "should return N_FAIL"
    assert r.N_FAIL == 1  # only S00 (id=0) fails


def test_fail_samples_zero_when_no_failures(test_db):
    """Variants without any failures should return N_FAIL=0."""
    db = Database(test_db)
    # chr1:3500 G>C — fail_ids=[] from conftest
    results = db.query(chrom="chr1", pos=3500)
    assert results
    r = results[0]
    assert r.N_FAIL >= 0
    assert r.N_FAIL == 0


def test_fail_samples_chrM(test_db):
    """chrM:100 C>A has fail_ids=[5] from conftest, but S05 is WES_kit_A (no chrM coverage).
    Only WGS samples are eligible at chrM — S05 is not eligible, so N_FAIL=0."""
    db = Database(test_db)
    results = db.query(chrom="chrM", pos=100)
    assert results
    r = next((x for x in results if x.variant.ref == "C" and x.variant.alt == "A"), None)
    assert r is not None
    assert r.N_FAIL >= 0
    # S05 (fail carrier) is WES_kit_A — WES_kit_A BED doesn't cover chrM.
    # Only WGS samples are eligible, none of which fail at this position.
    assert r.N_FAIL == 0


def test_fail_samples_filtered_by_eligible(test_db):
    """N_FAIL should only count samples in the eligible bitmap (tech coverage + filters)."""
    db = Database(test_db)
    # S05 is tech_id=1 (WES_kit_A). chr1:1500 is covered by WES_kit_A (bed has chr1:999-2000).
    # fail_ids for chr1:1500 is [0]. S00 is tech_id=0 (WGS).
    # Filter to tech WES_kit_A — S00 is not in that tech, so N_FAIL should be 0.
    results = db.query(chrom="chr1", pos=1500, tech=["WES_kit_A"])
    if results:
        r = next((x for x in results if x.variant.ref == "A"), None)
        if r:
            assert r.N_FAIL >= 0
            # S00 (fail) is WGS, not WES_kit_A — not eligible under this filter
            assert r.N_FAIL == 0

def test_fail_samples_in_query_region(test_db):
    """query_region also returns N_FAIL."""
    db = Database(test_db)
    results = db.query_region(chrom="chr1", start=1000, end=2000)
    assert results
    for r in results:
        assert r.N_FAIL >= 0, "query_region should return N_FAIL"


def test_fail_samples_in_query_batch(test_db):
    """query_batch also returns N_FAIL."""
    db = Database(test_db)
    results = db.query_batch(chrom="chr1", variants=[(1500, "A", "T"), (3500, "G", "C")])
    assert results
    for r in results:
        assert r.N_FAIL >= 0, "query_batch should return N_FAIL"


# ---------------------------------------------------------------------------
# Tests: run_preprocess with FILTER field
# ---------------------------------------------------------------------------

def _make_vcf_manifest(tmp_path: Path, records_by_sample: dict) -> Path:
    """Create VCFs and manifest for run_preprocess.
    records_by_sample = {sample_name: [(chrom, pos, ref, alt, gt, filter_str), ...]}
    """
    manifest_rows = ["sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes"]
    for sample_name, records in records_by_sample.items():
        vcf_path = tmp_path / f"{sample_name}.vcf"
        _write_vcf_with_filter(vcf_path, sample_name, records)
        manifest_rows.append(f"{sample_name}\tmale\twgs\t{vcf_path}\tE11.9")
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text("\n".join(manifest_rows) + "\n")
    return manifest_path


def test_pass_only_default_in_preprocess(tmp_path):
    """Default preprocess (pass_only=True): non-PASS calls don't count in AC."""
    records_by_sample = {
        "S00": [("chr1", 1000, "A", "T", "0/1", "PASS")],     # passes → AC counts
        "S01": [("chr1", 1000, "A", "T", "0/1", "LowQual")],  # fails → AC does NOT count
    }
    manifest_path = _make_vcf_manifest(tmp_path, records_by_sample)
    db_path = tmp_path / "db_pass"
    db_path.mkdir()
    run_preprocess(
        manifest_path=str(manifest_path),
        output_dir=str(db_path),
        genome_build="GRCh37",
        threads=1,
    )
    db = Database(str(db_path))
    results = db.query(chrom="chr1", pos=1000)
    assert results, "Variant should be present (S01 fails but still tracked in fail_bitmap)"
    r = results[0]
    # AC should only count S00 (PASS) — not S01 (LowQual)
    assert r.AC == 1, f"Expected AC=1 (PASS-only), got AC={r.AC}"
    assert r.N_FAIL == 1, f"Expected N_FAIL=1, got {r.N_FAIL}"


def test_variant_with_only_fail_appears(tmp_path):
    """A variant where all carriers fail FILTER (AC=0) still appears with FAIL>0."""
    records_by_sample = {
        "S00": [("chr1", 2000, "G", "C", "0/1", "LowQual")],  # only carrier, fails
    }
    manifest_path = _make_vcf_manifest(tmp_path, records_by_sample)
    db_path = tmp_path / "db_fail_only"
    db_path.mkdir()
    run_preprocess(
        manifest_path=str(manifest_path),
        output_dir=str(db_path),
        genome_build="GRCh37",
        threads=1,
    )
    db = Database(str(db_path))
    results = db.query(chrom="chr1", pos=2000)
    assert results, "Variant with only fail carriers should still appear"
    r = results[0]
    assert r.AC == 0, f"Expected AC=0 (carrier fails FILTER), got AC={r.AC}"
    assert r.N_FAIL == 1, f"Expected N_FAIL=1, got {r.N_FAIL}"


def test_schema_version_in_manifest(tmp_path):
    """run_preprocess writes schema_version=2.0 in manifest.json."""
    records_by_sample = {
        "S00": [("chr1", 1000, "A", "T", "0/1", "PASS")],
    }
    manifest_path = _make_vcf_manifest(tmp_path, records_by_sample)
    db_path = tmp_path / "db_schema"
    db_path.mkdir()
    run_preprocess(
        manifest_path=str(manifest_path),
        output_dir=str(db_path),
        genome_build="GRCh37",
        threads=1,
    )
    manifest = json.loads((db_path / "manifest.json").read_text())
    assert manifest.get("schema_version") == "2.0"
    assert manifest.get("pass_only_filter") is True
