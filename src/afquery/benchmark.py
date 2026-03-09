import json
import sqlite3
import tempfile
import time
from pathlib import Path

import pyarrow.parquet as pq

from .database import Database


def run_benchmark(db_path: Path, n_warmup: int = 3) -> dict:
    """
    Run a battery of queries against db_path and report timing in milliseconds.

    Returns:
        {
          "point_query_cold_ms": float,
          "point_query_warm_ms": float,   # median of n_warmup repeats
          "batch_100_ms": float,
          "batch_1000_ms": float,
          "targets": {
              "point_cold_ok": bool,      # <100ms
              "point_warm_ok": bool,      # <10ms
              "batch_100_ok": bool,       # <200ms
          }
        }
    """
    db_path = Path(db_path)
    db = Database(str(db_path))

    # Discover a valid (chrom, pos) with data
    test_chrom, test_pos = _find_test_variant(db_path)
    if test_chrom is None:
        return {"error": "No variants found in database"}

    # Get an ICD10 code present in the DB
    con = sqlite3.connect(str(db_path / "metadata.sqlite"))
    icd10_row = con.execute(
        "SELECT icd10_code FROM sample_icd10 LIMIT 1"
    ).fetchone()
    con.close()
    icd10_codes = [icd10_row[0]] if icd10_row else ["E11.9"]

    # Cold point query (first call — cold I/O cache)
    t0 = time.perf_counter()
    db.query(chrom=test_chrom, pos=test_pos, icd10=icd10_codes)
    cold_ms = (time.perf_counter() - t0) * 1000

    # Warm point queries (n_warmup repeats, take median)
    warm_times = []
    for _ in range(n_warmup):
        t0 = time.perf_counter()
        db.query(chrom=test_chrom, pos=test_pos, icd10=icd10_codes)
        warm_times.append((time.perf_counter() - t0) * 1000)
    warm_ms = sorted(warm_times)[len(warm_times) // 2]

    # Batch 100 query
    batch_100 = [test_pos + i * 1000 for i in range(100)]
    t0 = time.perf_counter()
    db.query_batch(chrom=test_chrom, positions=batch_100, icd10=icd10_codes)
    batch_100_ms = (time.perf_counter() - t0) * 1000

    # Batch 1000 query
    batch_1000 = [test_pos + i * 1000 for i in range(1000)]
    t0 = time.perf_counter()
    db.query_batch(chrom=test_chrom, positions=batch_1000, icd10=icd10_codes)
    batch_1000_ms = (time.perf_counter() - t0) * 1000

    return {
        "point_query_cold_ms": round(cold_ms, 2),
        "point_query_warm_ms": round(warm_ms, 2),
        "batch_100_ms": round(batch_100_ms, 2),
        "batch_1000_ms": round(batch_1000_ms, 2),
        "targets": {
            "point_cold_ok": cold_ms < 100,
            "point_warm_ok": warm_ms < 10,
            "batch_100_ok": batch_100_ms < 200,
        },
    }


def _find_test_variant(db_path: Path) -> tuple[str | None, int | None]:
    """Return (chrom, pos) of the first variant found in the DB."""
    variants_dir = db_path / "variants"
    if not variants_dir.exists():
        return None, None

    for entry in sorted(variants_dir.iterdir()):
        if entry.suffix == ".parquet":
            tbl = pq.read_table(str(entry), columns=["pos"])
            if len(tbl) > 0:
                return entry.stem, int(tbl["pos"][0].as_py())
        elif entry.is_dir():
            for bucket in sorted(entry.glob("bucket_*.parquet")):
                tbl = pq.read_table(str(bucket), columns=["pos"])
                if len(tbl) > 0:
                    return entry.name, int(tbl["pos"][0].as_py())

    return None, None


def run_benchmark_with_synth(
    n_samples: int = 1000,
    n_variants: int = 10_000,
    output_report: str = "benchmark_report.json",
) -> dict:
    """Build a synthetic DB, run benchmarks, and save a JSON report."""
    from .preprocess import run_preprocess
    from .preprocess.synth import generate_synthetic_manifest

    with tempfile.TemporaryDirectory(prefix="afquery_bench_") as tmp_root:
        tmp_path = Path(tmp_root)
        synth_dir = tmp_path / "synth"
        db_dir = tmp_path / "db"
        db_dir.mkdir()

        manifest_path = generate_synthetic_manifest(
            synth_dir,
            n_samples=n_samples,
            n_variants_per_chrom=n_variants,
            chroms=("chr1",),
        )

        run_preprocess(
            manifest_path=str(manifest_path),
            output_dir=str(db_dir),
            genome_build="GRCh37",
        )

        results = run_benchmark(db_dir)

    with open(output_report, "w") as f:
        json.dump(results, f, indent=2)

    return results
