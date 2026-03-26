#!/usr/bin/env python3
"""Dump variants from all DBs, merge, compute NADI and AF error metrics.

Inputs:
  - dbs/db_wgs/                WGS ground truth database
  - dbs/db_{scenario}/         WES databases (one per scenario)
  - beds/afquery/*.bed         for computing per-position kit coverage

Outputs:
  - results/merged.parquet     merged AF comparison table (all scenarios)
  - results/nadi_summary.json  summary statistics per scenario
"""

import json
import logging
import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from config import (
    AFQUERY_BED_DIR,
    DB_DIR,
    RESULTS_DIR,
    SCENARIOS,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)


def dump_db(
    db_path: Path, by_tech: bool = False, include_ac_zero: bool = False,
) -> pd.DataFrame:
    """Dump all variants from an AFQuery database to a DataFrame."""
    from afquery.dump import _build_groups, dump_database
    from afquery.models import SampleFilter
    from afquery.query import QueryEngine

    engine = QueryEngine(str(db_path))
    sf = SampleFilter()

    groups = []
    if by_tech:
        groups = _build_groups(
            engine, sf,
            by_sex=False, by_tech=True, by_phenotype=[], all_groups=False,
        )

    buf = StringIO()
    dump_database(
        engine, buf, sf, groups,
        chrom_filter=None,
        pos_start=None, pos_end=None,
        n_workers=None,
        include_ac_zero=include_ac_zero,
    )
    buf.seek(0)
    return pd.read_csv(buf)


def compute_kit_coverage(positions: pd.Series, bed_dir: Path) -> pd.DataFrame:
    """For each unique position, count how many kits cover it.

    Returns DataFrame with columns: pos, n_kits, kit_pattern.
    """
    from afquery.capture import CaptureIndex

    kits = {}
    for bed_file in sorted(bed_dir.glob("*.bed")):
        tech_name = bed_file.stem
        kits[tech_name] = CaptureIndex.from_bed(str(bed_file))

    unique_positions = sorted(positions.unique())
    rows = []
    for pos in unique_positions:
        covering = [
            name for name, idx in kits.items() if idx.covers("chr22", int(pos))
        ]
        rows.append({
            "pos": pos,
            "n_kits": len(covering),
            "kit_pattern": "+".join(covering) if covering else "none",
        })

    return pd.DataFrame(rows)


def merge_scenario(
    wgs_df: pd.DataFrame,
    wes_df: pd.DataFrame,
    scenario: str,
    coverage_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge WGS and WES dumps, compute AF_naive, errors, NADI."""
    merged = wgs_df[["chrom", "pos", "ref", "alt", "AC", "AN", "AF"]].merge(
        wes_df[["chrom", "pos", "ref", "alt", "AC", "AN", "AF"]],
        on=["chrom", "pos", "ref", "alt"],
        how="outer",
        suffixes=("_wgs", "_wes"),
    )

    # Variants in WGS but not WES: AC_wes=0 (carriers lost to masking)
    merged["AC_wes"] = merged["AC_wes"].fillna(0).astype(int)
    merged["AN_wes"] = merged["AN_wes"].fillna(0).astype(int)
    merged["AC_wgs"] = merged["AC_wgs"].fillna(0).astype(int)
    merged["AN_wgs"] = merged["AN_wgs"].fillna(0).astype(int)

    # Recompute AFs to avoid float precision issues from CSV round-trip
    merged["AF_wgs"] = np.where(
        merged["AN_wgs"] > 0, merged["AC_wgs"] / merged["AN_wgs"], 0.0,
    )
    merged["AF_afquery"] = np.where(
        merged["AN_wes"] > 0, merged["AC_wes"] / merged["AN_wes"], np.nan,
    )
    # AF_naive = AC_masked / AN_wgs  (naive: assume all samples cover everything)
    merged["AF_naive"] = np.where(
        merged["AN_wgs"] > 0, merged["AC_wes"] / merged["AN_wgs"], 0.0,
    )

    # Kit coverage info
    merged = merged.merge(coverage_df, on="pos", how="left")

    # Errors (only for variants covered by at least 1 kit)
    covered = merged["AN_wes"] > 0
    merged["error_afquery"] = np.where(
        covered, np.abs(merged["AF_afquery"] - merged["AF_wgs"]), np.nan,
    )
    merged["error_naive"] = np.where(
        covered, np.abs(merged["AF_naive"] - merged["AF_wgs"]), np.nan,
    )

    # NADI: |AN_method - AN_true| / AN_true  where AN_true = AN_wes (2*K)
    merged["NADI_naive"] = np.where(
        merged["AN_wes"] > 0,
        np.abs(merged["AN_wgs"] - merged["AN_wes"]) / merged["AN_wes"],
        np.nan,
    )
    merged["NADI_afquery"] = 0.0

    # AN ratio for histogram (AN_naive / AN_AFQuery = AN_wgs / AN_wes)
    merged["AN_ratio"] = np.where(
        merged["AN_wes"] > 0, merged["AN_wgs"] / merged["AN_wes"], np.nan,
    )

    merged["scenario"] = scenario
    return merged


def compute_summary(merged: pd.DataFrame) -> dict:
    """Compute summary statistics for one scenario."""
    df = merged[merged["AN_wes"] > 0].copy()

    summary = {
        "n_variants_total": int(len(merged)),
        "n_variants_covered": int(len(df)),
        "n_variants_by_nkits": {
            str(k): int(v)
            for k, v in df["n_kits"].value_counts().sort_index().items()
        },
    }

    for method in ["afquery", "naive"]:
        errs = df[f"error_{method}"].dropna()
        summary[f"{method}_mae"] = float(errs.mean()) if len(errs) else 0.0
        summary[f"{method}_median_ae"] = float(errs.median()) if len(errs) else 0.0
        summary[f"{method}_max_error"] = float(errs.max()) if len(errs) else 0.0
        summary[f"{method}_rmse"] = (
            float(np.sqrt((errs**2).mean())) if len(errs) else 0.0
        )

        # Stratified by n_kits
        summary[f"{method}_by_nkits"] = {}
        for nk in sorted(df["n_kits"].dropna().unique()):
            sub_err = df.loc[df["n_kits"] == nk, f"error_{method}"].dropna()
            summary[f"{method}_by_nkits"][str(int(nk))] = {
                "n": int(len(sub_err)),
                "mae": float(sub_err.mean()) if len(sub_err) else 0.0,
            }

    # NADI summary
    nadi = df["NADI_naive"].dropna()
    summary["nadi_naive_mean"] = float(nadi.mean()) if len(nadi) else 0.0
    summary["nadi_naive_median"] = float(nadi.median()) if len(nadi) else 0.0
    summary["nadi_naive_max"] = float(nadi.max()) if len(nadi) else 0.0

    return summary


# ---- Sanity checks --------------------------------------------------------

def run_sanity_checks(merged: pd.DataFrame, scenario: str):
    """Validate invariants on the merged data."""
    df = merged[merged["AN_wes"] > 0]
    errors = []

    # AC monotonicity: masking cannot create alleles
    violations = (df["AC_wes"] > df["AC_wgs"]).sum()
    if violations > 0:
        errors.append(f"AC_wes > AC_wgs in {violations} rows")

    # AN ordering: WES AN <= WGS AN
    violations = (df["AN_wes"] > df["AN_wgs"]).sum()
    if violations > 0:
        errors.append(f"AN_wes > AN_wgs in {violations} rows")

    # 3-kit convergence: where all kits cover, AN should be equal
    three_kit = df[df["n_kits"] == 3]
    if len(three_kit) > 0:
        an_mismatch = (three_kit["AN_wes"] != three_kit["AN_wgs"]).sum()
        if an_mismatch > 0:
            errors.append(
                f"AN_wes != AN_wgs for {an_mismatch}/{len(three_kit)} "
                f"3-kit variants"
            )

    # AF_AFQuery should be unbiased (mean error centered near 0)
    bias = (df["AF_afquery"] - df["AF_wgs"]).mean()
    if abs(bias) > 0.01:
        errors.append(f"AFQuery bias = {bias:.4f} (expected ~0)")

    # AF_naive should be biased low
    naive_bias = (df["AF_naive"] - df["AF_wgs"]).mean()
    if naive_bias > 0.001:
        errors.append(f"Naive bias = {naive_bias:.4f} (expected < 0)")

    if errors:
        for e in errors:
            logger.error("  SANITY CHECK FAILED [%s]: %s", scenario, e)
    else:
        logger.info("  All sanity checks passed for %s", scenario)
        logger.info(
            "    AFQuery bias=%.6f, naive bias=%.6f", bias, naive_bias,
        )


def main():
    ensure_dirs()

    # 1. Dump WGS ground truth (shared across scenarios)
    logger.info("Dumping WGS ground truth database...")
    wgs_df = dump_db(DB_DIR / "db_wgs", by_tech=False, include_ac_zero=True)
    logger.info("WGS dump: %d variants", len(wgs_df))

    # 2. Compute kit coverage per position (shared — same BEDs for all scenarios)
    logger.info(
        "Computing kit coverage for %d unique positions...",
        wgs_df["pos"].nunique(),
    )
    coverage_df = compute_kit_coverage(wgs_df["pos"], AFQUERY_BED_DIR)
    kit_dist = coverage_df["n_kits"].value_counts().sort_index()
    logger.info("Kit coverage distribution:\n%s", kit_dist.to_string())

    # 3. Process each scenario
    all_merged = []
    summaries = {}

    for scenario_name in SCENARIOS:
        logger.info("=== Processing scenario: %s ===", scenario_name)

        wes_df = dump_db(
            DB_DIR / f"db_{scenario_name}",
            by_tech=True, include_ac_zero=True,
        )
        logger.info("  WES dump: %d variants", len(wes_df))

        merged = merge_scenario(wgs_df, wes_df, scenario_name, coverage_df)
        run_sanity_checks(merged, scenario_name)
        all_merged.append(merged)

        summaries[scenario_name] = compute_summary(merged)
        logger.info(
            "  MAE  naive=%.6f  AFQuery=%.6f",
            summaries[scenario_name]["naive_mae"],
            summaries[scenario_name]["afquery_mae"],
        )

    # 4. Save results
    full_df = pd.concat(all_merged, ignore_index=True)
    full_df.to_parquet(RESULTS_DIR / "merged.parquet", index=False)
    logger.info("Saved merged.parquet (%d rows)", len(full_df))

    (RESULTS_DIR / "nadi_summary.json").write_text(
        json.dumps(summaries, indent=2)
    )
    logger.info("Saved nadi_summary.json")


if __name__ == "__main__":
    main()
