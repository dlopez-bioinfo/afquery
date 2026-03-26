#!/usr/bin/env python3
"""Experiment 4: AFQuery vs. bcftools comparison.

Compares equivalent operations:
  - Point query (single locus AF)
  - Subset query (AF filtered by sex)
  - Full dump (whole-chromosome AF export)

Also validates AF concordance between both tools.

Output: results/bcftools_comparison.json, results/concordance.json
"""

import csv
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from shared.utils import stats as _stats  # noqa: E402
from config import (
    BCFTOOLS_REPS,
    ONEKG_DB_DIR,
    ONEKG_DIR,
    ONEKG_MERGED_VCF,
    ONEKG_SUBSETS,
    RESULTS_DIR,
    SEED,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def _check_bcftools():
    """Verify bcftools is available."""
    try:
        result = subprocess.run(
            ["bcftools", "--version"], capture_output=True, text=True
        )
        version = result.stdout.split("\n")[0]
        logger.info("Using %s", version)
        return True
    except FileNotFoundError:
        logger.error("bcftools not found in PATH")
        return False


def _create_subset_vcf(merged_vcf: Path, sample_list: list[str], output_vcf: Path):
    """Create a subset multi-sample VCF with only the specified samples."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("\n".join(sample_list) + "\n")
        samples_file = f.name

    try:
        subprocess.run(
            [
                "bcftools", "view",
                "-S", samples_file,
                "--force-samples",
                "-Oz", "-o", str(output_vcf),
                str(merged_vcf),
            ],
            check=True,
            capture_output=True,
        )
        subprocess.run(
            ["bcftools", "index", "-t", str(output_vcf)],
            check=True,
            capture_output=True,
        )
    finally:
        os.unlink(samples_file)


def _get_sample_list(db_path: Path, sex: str = None) -> list[str]:
    """Get sample names from the AFQuery database, optionally filtered by sex."""
    from afquery.database import Database

    db = Database(str(db_path))
    samples = db.list_samples()
    if sex:
        samples = [s for s in samples if s["sex"] == sex]
    return [s["sample_name"] for s in samples]


def _find_test_locus(db_path: Path) -> tuple[str, int]:
    """Find a variant position with data in the database."""
    from afquery.benchmark import _find_test_variants

    variants = _find_test_variants(db_path, n=10)
    if not variants:
        raise RuntimeError("No variants found")
    return variants[0][0], variants[0][1]


def _time_bcftools_point(vcf_path: str, chrom: str, pos: int) -> float:
    """Time a bcftools point query (AF at a single locus). Returns ms."""
    cmd = (
        f"bcftools view -r {chrom}:{pos} {vcf_path} "
        f"| bcftools +fill-tags -- -t AC,AN,AF "
        f"| bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AC\\t%AN\\t%AF\\n'"
    )
    t0 = time.perf_counter()
    subprocess.run(cmd, shell=True, capture_output=True)
    return (time.perf_counter() - t0) * 1000


def _time_bcftools_subset_point(
    vcf_path: str, chrom: str, pos: int, samples_file: str
) -> float:
    """Time a bcftools subset point query (AF for a subset of samples). Returns ms."""
    cmd = (
        f"bcftools view -r {chrom}:{pos} -S {samples_file} {vcf_path} "
        f"| bcftools +fill-tags -- -t AC,AN,AF "
        f"| bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AC\\t%AN\\t%AF\\n'"
    )
    t0 = time.perf_counter()
    subprocess.run(cmd, shell=True, capture_output=True)
    return (time.perf_counter() - t0) * 1000


def _time_bcftools_dump(vcf_path: str, output_path: str) -> float:
    """Time a full-chromosome AF dump with bcftools. Returns ms."""
    cmd = (
        f"bcftools +fill-tags {vcf_path} -- -t AC,AN,AF "
        f"| bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AC\\t%AN\\t%AF\\n' "
        f"> {output_path}"
    )
    t0 = time.perf_counter()
    subprocess.run(cmd, shell=True, capture_output=True)
    return (time.perf_counter() - t0) * 1000


def _time_afquery_point(db_path: str, chrom: str, pos: int) -> float:
    """Time an AFQuery point query. Returns ms."""
    from afquery.database import Database

    db = Database(db_path)
    t0 = time.perf_counter()
    db.query(chrom=chrom, pos=pos)
    return (time.perf_counter() - t0) * 1000


def _time_afquery_subset_point(
    db_path: str, chrom: str, pos: int, sex: str
) -> float:
    """Time an AFQuery subset point query. Returns ms."""
    from afquery.database import Database

    db = Database(db_path)
    t0 = time.perf_counter()
    db.query(chrom=chrom, pos=pos, sex=sex)
    return (time.perf_counter() - t0) * 1000


def _time_afquery_dump(db_path: str, chrom: str, output_path: str) -> float:
    """Time an AFQuery full-chromosome dump. Returns ms."""
    from afquery.database import Database

    db = Database(db_path)
    t0 = time.perf_counter()
    db.dump(output=output_path, chrom=chrom)
    return (time.perf_counter() - t0) * 1000


def _compute_concordance(db_path: Path, vcf_path: str, chrom: str) -> dict:
    """Compare AF values between AFQuery dump and bcftools +fill-tags."""
    with tempfile.TemporaryDirectory(prefix="conc_") as tmpdir:
        tmpdir = Path(tmpdir)

        # AFQuery dump
        afq_out = str(tmpdir / "afquery_dump.csv")
        from afquery.database import Database

        db = Database(str(db_path))
        db.dump(output=afq_out, chrom=chrom)

        # bcftools dump
        bcf_out = str(tmpdir / "bcftools_dump.tsv")
        cmd = (
            f"bcftools +fill-tags {vcf_path} -- -t AC,AN,AF "
            f"| bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AC\\t%AN\\t%AF\\n' "
            f"> {bcf_out}"
        )
        subprocess.run(cmd, shell=True, check=True, capture_output=True)

        # Parse AFQuery results
        afq_data = {}
        with open(afq_out) as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
                afq_data[key] = {
                    "AC": int(row["AC"]),
                    "AN": int(row["AN"]),
                    "AF": float(row["AF"]) if row["AF"] != "" else 0.0,
                }

        # Parse bcftools results
        bcf_data = {}
        with open(bcf_out) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                chrom_val, pos, ref, alt, ac, an, af = parts[:7]
                # Handle multi-allelic (comma-separated AC/AF)
                acs = ac.split(",")
                afs = af.split(",")
                alts = alt.split(",")
                for i, a in enumerate(alts):
                    key = (chrom_val, int(pos), ref, a)
                    bcf_data[key] = {
                        "AC": int(acs[i]) if i < len(acs) else 0,
                        "AN": int(an),
                        "AF": float(afs[i]) if i < len(afs) else 0.0,
                    }

        # Compute concordance
        common_keys = set(afq_data.keys()) & set(bcf_data.keys())
        n_common = len(common_keys)
        n_afq_only = len(set(afq_data.keys()) - set(bcf_data.keys()))
        n_bcf_only = len(set(bcf_data.keys()) - set(afq_data.keys()))

        af_pairs = []
        ac_mismatches = 0
        an_mismatches = 0
        for key in common_keys:
            afq = afq_data[key]
            bcf = bcf_data[key]
            af_pairs.append((afq["AF"], bcf["AF"]))
            if afq["AC"] != bcf["AC"]:
                ac_mismatches += 1
            if afq["AN"] != bcf["AN"]:
                an_mismatches += 1

        # R-squared
        if len(af_pairs) > 1:
            afq_afs = [p[0] for p in af_pairs]
            bcf_afs = [p[1] for p in af_pairs]
            mean_afq = sum(afq_afs) / len(afq_afs)
            mean_bcf = sum(bcf_afs) / len(bcf_afs)
            ss_res = sum((a - b) ** 2 for a, b in af_pairs)
            ss_tot = sum((b - mean_bcf) ** 2 for b in bcf_afs)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
        else:
            r_squared = None

        return {
            "n_common_variants": n_common,
            "n_afquery_only": n_afq_only,
            "n_bcftools_only": n_bcf_only,
            "ac_mismatches": ac_mismatches,
            "an_mismatches": an_mismatches,
            "r_squared": round(r_squared, 6) if r_squared is not None else None,
            "af_pairs_sample": af_pairs[:100],  # first 100 for plotting
        }


def _bench_one_subset(n_samples: int, output_path: Path):
    """Run bcftools comparison for one subset and write a per-run JSON."""
    if not _check_bcftools():
        raise RuntimeError("bcftools not found")
    if not ONEKG_MERGED_VCF.exists():
        raise FileNotFoundError(f"1KG merged VCF not found: {ONEKG_MERGED_VCF}")

    db_path = ONEKG_DB_DIR / f"1kg_{n_samples}"
    if not (db_path / "manifest.json").exists():
        raise FileNotFoundError(f"No 1KG DB for {n_samples} samples")

    logger.info("=== bcftools comparison: %d samples ===", n_samples)

    with tempfile.TemporaryDirectory(prefix="bcftools_bench_") as tmpdir:
        tmpdir = Path(tmpdir)

        sample_list = _get_sample_list(db_path)
        subset_vcf = tmpdir / f"subset_{n_samples}.vcf.gz"
        logger.info("  Creating subset VCF (%d samples)...", len(sample_list))
        _create_subset_vcf(ONEKG_MERGED_VCF, sample_list, subset_vcf)

        chrom, pos = _find_test_locus(db_path)
        logger.info("  Test locus: %s:%d", chrom, pos)

        female_samples = _get_sample_list(db_path, sex="female")
        females_file = tmpdir / f"females_{n_samples}.txt"
        females_file.write_text("\n".join(female_samples) + "\n")

        result = {"n_samples": n_samples, "chrom": chrom}

        logger.info("  Point query benchmark...")
        bcf_times = [_time_bcftools_point(str(subset_vcf), chrom, pos) for _ in range(BCFTOOLS_REPS)]
        afq_times = [_time_afquery_point(str(db_path), chrom, pos) for _ in range(BCFTOOLS_REPS)]
        result["point_query"] = {"bcftools": _stats(bcf_times), "afquery": _stats(afq_times)}

        logger.info("  Subset query benchmark (females)...")
        bcf_sub_times = [
            _time_bcftools_subset_point(str(subset_vcf), chrom, pos, str(females_file))
            for _ in range(BCFTOOLS_REPS)
        ]
        afq_sub_times = [
            _time_afquery_subset_point(str(db_path), chrom, pos, "female")
            for _ in range(BCFTOOLS_REPS)
        ]
        result["subset_query"] = {"bcftools": _stats(bcf_sub_times), "afquery": _stats(afq_sub_times)}

        logger.info("  Full dump benchmark...")
        bcf_dump_out = str(tmpdir / f"bcf_dump_{n_samples}.tsv")
        afq_dump_out = str(tmpdir / f"afq_dump_{n_samples}.csv")
        bcf_dump_times = [
            _time_bcftools_dump(str(subset_vcf), bcf_dump_out)
            for _ in range(min(3, BCFTOOLS_REPS))
        ]
        afq_dump_times = [
            _time_afquery_dump(str(db_path), chrom, afq_dump_out)
            for _ in range(min(3, BCFTOOLS_REPS))
        ]
        result["dump"] = {"bcftools": _stats(bcf_dump_times), "afquery": _stats(afq_dump_times)}

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2))
    logger.info("Result saved to %s", output_path)


def _bench_concordance(n_samples: int, output_path: Path):
    """Compute AF concordance for one subset and write concordance JSON."""
    if not _check_bcftools():
        raise RuntimeError("bcftools not found")
    if not ONEKG_MERGED_VCF.exists():
        raise FileNotFoundError(f"1KG merged VCF not found: {ONEKG_MERGED_VCF}")

    db_path = ONEKG_DB_DIR / f"1kg_{n_samples}"
    if not (db_path / "manifest.json").exists():
        raise FileNotFoundError(f"No 1KG DB for {n_samples} samples")

    logger.info("=== Computing AF concordance (%d samples) ===", n_samples)

    with tempfile.TemporaryDirectory(prefix="bcftools_conc_") as tmpdir:
        tmpdir = Path(tmpdir)
        sample_list = _get_sample_list(db_path)
        subset_vcf = tmpdir / f"subset_{n_samples}.vcf.gz"
        _create_subset_vcf(ONEKG_MERGED_VCF, sample_list, subset_vcf)
        chrom, _ = _find_test_locus(db_path)
        concordance = _compute_concordance(db_path, str(subset_vcf), chrom)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(concordance, indent=2))
    logger.info("Concordance: R² = %s", concordance.get("r_squared"))
    logger.info("Result saved to %s", output_path)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--n-samples", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--concordance-only", action="store_true")
    args = parser.parse_args()

    ensure_dirs()

    if args.n_samples is not None and args.output:
        if args.concordance_only:
            _bench_concordance(args.n_samples, args.output)
        else:
            _bench_one_subset(args.n_samples, args.output)
        return

    # Full-sweep mode (backward compatible)
    if not _check_bcftools():
        return
    if not ONEKG_MERGED_VCF.exists():
        logger.error("1KG merged VCF not found at %s.", ONEKG_MERGED_VCF)
        return

    all_results = []

    with tempfile.TemporaryDirectory(prefix="bcftools_bench_") as tmpdir:
        tmpdir = Path(tmpdir)

        for n_samples in ONEKG_SUBSETS:
            db_path = ONEKG_DB_DIR / f"1kg_{n_samples}"
            if not (db_path / "manifest.json").exists():
                logger.warning("No 1KG DB for %d samples, skipping", n_samples)
                continue

            logger.info("=== bcftools comparison: %d samples ===", n_samples)

            sample_list = _get_sample_list(db_path)
            subset_vcf = tmpdir / f"subset_{n_samples}.vcf.gz"
            logger.info("  Creating subset VCF (%d samples)...", len(sample_list))
            _create_subset_vcf(ONEKG_MERGED_VCF, sample_list, subset_vcf)

            chrom, pos = _find_test_locus(db_path)
            logger.info("  Test locus: %s:%d", chrom, pos)

            female_samples = _get_sample_list(db_path, sex="female")
            females_file = tmpdir / f"females_{n_samples}.txt"
            females_file.write_text("\n".join(female_samples) + "\n")

            result = {"n_samples": n_samples, "chrom": chrom}

            logger.info("  Point query benchmark...")
            bcf_times = [_time_bcftools_point(str(subset_vcf), chrom, pos) for _ in range(BCFTOOLS_REPS)]
            afq_times = [_time_afquery_point(str(db_path), chrom, pos) for _ in range(BCFTOOLS_REPS)]
            result["point_query"] = {"bcftools": _stats(bcf_times), "afquery": _stats(afq_times)}

            logger.info("  Subset query benchmark (females)...")
            bcf_sub_times = [
                _time_bcftools_subset_point(str(subset_vcf), chrom, pos, str(females_file))
                for _ in range(BCFTOOLS_REPS)
            ]
            afq_sub_times = [
                _time_afquery_subset_point(str(db_path), chrom, pos, "female")
                for _ in range(BCFTOOLS_REPS)
            ]
            result["subset_query"] = {"bcftools": _stats(bcf_sub_times), "afquery": _stats(afq_sub_times)}

            logger.info("  Full dump benchmark...")
            bcf_dump_out = str(tmpdir / f"bcf_dump_{n_samples}.tsv")
            afq_dump_out = str(tmpdir / f"afq_dump_{n_samples}.csv")
            bcf_dump_times = [
                _time_bcftools_dump(str(subset_vcf), bcf_dump_out)
                for _ in range(min(3, BCFTOOLS_REPS))
            ]
            afq_dump_times = [
                _time_afquery_dump(str(db_path), chrom, afq_dump_out)
                for _ in range(min(3, BCFTOOLS_REPS))
            ]
            result["dump"] = {"bcftools": _stats(bcf_dump_times), "afquery": _stats(afq_dump_times)}

            all_results.append(result)

        # Concordance on largest subset
        max_subset = max(ONEKG_SUBSETS)
        db_path = ONEKG_DB_DIR / f"1kg_{max_subset}"
        if (db_path / "manifest.json").exists():
            logger.info("=== Computing AF concordance (%d samples) ===", max_subset)
            sample_list = _get_sample_list(db_path)
            subset_vcf = tmpdir / f"subset_{max_subset}.vcf.gz"
            if not subset_vcf.exists():
                _create_subset_vcf(ONEKG_MERGED_VCF, sample_list, subset_vcf)
            chrom, _ = _find_test_locus(db_path)
            concordance = _compute_concordance(db_path, str(subset_vcf), chrom)
            conc_out = RESULTS_DIR / "concordance.json"
            conc_out.write_text(json.dumps(concordance, indent=2))
            logger.info("Concordance: R² = %s", concordance.get("r_squared"))

    out = RESULTS_DIR / "bcftools_comparison.json"
    out.write_text(json.dumps(all_results, indent=2))
    logger.info("Results saved to %s", out)


if __name__ == "__main__":
    main()
