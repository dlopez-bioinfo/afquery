import json
from pathlib import Path

import pytest

from afquery.benchmark import run_benchmark, run_benchmark_with_synth


def test_run_benchmark_returns_expected_keys(test_db):
    """run_benchmark returns a dict with all required timing keys."""
    results = run_benchmark(Path(test_db))

    assert "error" not in results, f"Unexpected error: {results.get('error')}"
    for key in ("point_query_cold_ms", "point_query_warm_ms", "batch_100_ms", "batch_1000_ms"):
        assert key in results, f"Missing key: {key}"
        assert isinstance(results[key], float)

    assert "targets" in results
    for t_key in ("point_cold_ok", "point_warm_ok", "batch_100_ok"):
        assert t_key in results["targets"]
        assert isinstance(results["targets"][t_key], bool)


def test_run_benchmark_positive_timings(test_db):
    """All timing values should be positive."""
    results = run_benchmark(Path(test_db))
    for key in ("point_query_cold_ms", "point_query_warm_ms", "batch_100_ms", "batch_1000_ms"):
        assert results[key] >= 0, f"{key} should be non-negative"


def test_benchmark_with_synth_small(tmp_path):
    """run_benchmark_with_synth with tiny dataset completes without error."""
    output = str(tmp_path / "report.json")
    results = run_benchmark_with_synth(
        n_samples=10,
        n_variants=100,
        output_report=output,
    )

    assert "error" not in results
    assert Path(output).exists()
    saved = json.loads(Path(output).read_text())
    assert saved["point_query_cold_ms"] == results["point_query_cold_ms"]


def test_benchmark_cli_runs(tmp_path, data_dir):
    """benchmark CLI command with --db-dir runs without error."""
    from click.testing import CliRunner
    from afquery.cli import cli
    from afquery.preprocess import run_preprocess

    # Build a tiny DB to benchmark
    db_path = str(tmp_path / "bench_db")
    run_preprocess(
        manifest_path=str(data_dir / "manifest.tsv"),
        output_dir=db_path,
        genome_build="GRCh37",
        bed_dir=str(data_dir / "beds"),
        n_threads=2,
    )

    output = str(tmp_path / "bench.json")
    runner = CliRunner()
    result = runner.invoke(cli, [
        "benchmark",
        "--db-dir", db_path,
        "--output", output,
    ])

    assert result.exit_code == 0, f"CLI failed: {result.output}"
    assert Path(output).exists()


def test_benchmark_no_variants_returns_error():
    """run_benchmark returns error dict when DB has no variants."""
    import tempfile, sqlite3, json, os
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        db_path = Path(d)
        (db_path / "variants").mkdir()
        (db_path / "capture").mkdir()
        # Minimal manifest.json
        (db_path / "manifest.json").write_text(json.dumps({
            "genome_build": "GRCh37",
            "sample_count": 0,
        }))
        # Minimal SQLite
        con = sqlite3.connect(str(db_path / "metadata.sqlite"))
        con.executescript("""
            CREATE TABLE samples (sample_id INTEGER PRIMARY KEY, sample_name TEXT, sex TEXT, tech_id INTEGER);
            CREATE TABLE technologies (tech_id INTEGER PRIMARY KEY, tech_name TEXT, bed_path TEXT);
            CREATE TABLE sample_icd10 (sample_id INTEGER, icd10_code TEXT, PRIMARY KEY(sample_id, icd10_code));
            CREATE TABLE precomputed_bitmaps (bitmap_type TEXT, bitmap_key TEXT, bitmap_data BLOB, PRIMARY KEY(bitmap_type, bitmap_key));
        """)
        con.close()

        result = run_benchmark(db_path)
        assert "error" in result
