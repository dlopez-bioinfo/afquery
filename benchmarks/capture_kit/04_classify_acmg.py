#!/usr/bin/env python3
"""Apply ACMG-like classification thresholds and compute misclassification rates.

Inputs:
  - results/merged.parquet

Outputs:
  - results/acmg_results.json    misclassification counts and confusion matrices
"""

import json
import logging
import sys
from pathlib import Path

import pandas as pd

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from config import (
    ACMG_THRESHOLDS,
    RESULTS_DIR,
    SCENARIOS,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger(__name__)

CLASSES = ["BA1", "BS1", "PM2", "neutral"]


def classify_variant(af: float, ac: int, thresholds: dict) -> str:
    """Classify variant into ACMG evidence category.

    BA1: AF > 5% (benign standalone)
    BS1: AF > disease-specific threshold (benign strong)
    PM2: AF <= threshold or absent (pathogenic moderate)
    neutral: none of the above
    """
    if af > thresholds["BA1"]:
        return "BA1"
    if af > thresholds["BS1"]:
        return "BS1"
    # PM2: absent (AC=0) when threshold is 0.0, otherwise AF <= threshold
    if thresholds["PM2"] == 0.0:
        if ac == 0:
            return "PM2"
    elif af <= thresholds["PM2"]:
        return "PM2"
    return "neutral"


def compute_confusion(true_classes, pred_classes):
    """Compute confusion matrix as nested dict."""
    matrix = {t: {p: 0 for p in CLASSES} for t in CLASSES}
    for t, p in zip(true_classes, pred_classes):
        matrix[t][p] += 1
    return matrix


def analyze_scenario_disease(df: pd.DataFrame, thresholds: dict) -> dict:
    """Analyze one scenario x disease combination."""
    sub = df[df["AN_wes"] > 0].copy()

    sub["class_truth"] = sub.apply(
        lambda r: classify_variant(r["AF_wgs"], int(r["AC_wgs"]), thresholds),
        axis=1,
    )
    sub["class_afquery"] = sub.apply(
        lambda r: classify_variant(
            r["AF_afquery"], int(r["AC_wes"]), thresholds,
        ),
        axis=1,
    )
    sub["class_naive"] = sub.apply(
        lambda r: classify_variant(r["AF_naive"], int(r["AC_wes"]), thresholds),
        axis=1,
    )

    n = len(sub)
    result = {
        "n_variants": n,
        "truth_distribution": sub["class_truth"].value_counts().to_dict(),
    }

    for method in ["afquery", "naive"]:
        col = f"class_{method}"
        concordant = int((sub[col] == sub["class_truth"]).sum())
        discordant = n - concordant

        # Critical errors with clinical direction
        # Benign (truth) -> PM2 (pred): false pathogenic evidence
        false_pathogenic = int(
            (
                sub["class_truth"].isin(["BA1", "BS1"])
                & (sub[col] == "PM2")
            ).sum()
        )
        # PM2 (truth) -> benign (pred): missed pathogenic evidence
        missed_pathogenic = int(
            (
                (sub["class_truth"] == "PM2")
                & sub[col].isin(["BA1", "BS1"])
            ).sum()
        )

        result[method] = {
            "concordant": concordant,
            "discordant": discordant,
            "discordance_rate": round(discordant / n, 6) if n > 0 else 0,
            "false_pathogenic": false_pathogenic,
            "missed_pathogenic": missed_pathogenic,
            "confusion_matrix": compute_confusion(
                sub["class_truth"], sub[col],
            ),
        }

    return result


def main():
    ensure_dirs()

    merged_path = RESULTS_DIR / "merged.parquet"
    if not merged_path.exists():
        logger.error("merged.parquet not found. Run 03_compute_metrics.py first.")
        sys.exit(1)

    df = pd.read_parquet(merged_path)
    logger.info("Loaded merged.parquet: %d rows", len(df))

    results = {}
    for disease, thresholds in ACMG_THRESHOLDS.items():
        results[disease] = {}
        for scenario in SCENARIOS:
            logger.info("Classifying: %s x %s", disease, scenario)
            scenario_df = df[df["scenario"] == scenario]
            results[disease][scenario] = analyze_scenario_disease(
                scenario_df, thresholds,
            )

            r = results[disease][scenario]
            logger.info(
                "  naive:   %d discordant (%.2f%%), "
                "FP=%d, MP=%d",
                r["naive"]["discordant"],
                r["naive"]["discordance_rate"] * 100,
                r["naive"]["false_pathogenic"],
                r["naive"]["missed_pathogenic"],
            )
            logger.info(
                "  afquery: %d discordant (%.2f%%), "
                "FP=%d, MP=%d",
                r["afquery"]["discordant"],
                r["afquery"]["discordance_rate"] * 100,
                r["afquery"]["false_pathogenic"],
                r["afquery"]["missed_pathogenic"],
            )

    out_path = RESULTS_DIR / "acmg_results.json"
    out_path.write_text(json.dumps(results, indent=2))
    logger.info("Saved %s", out_path)


if __name__ == "__main__":
    main()
