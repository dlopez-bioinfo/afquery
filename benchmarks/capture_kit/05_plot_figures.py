#!/usr/bin/env python3
"""Generate publication-quality figures for capture kit mixing benchmark.

Figures:
  1. AF error distribution (violin): AFQuery vs. naive x 3 scenarios
  2. AF error by kit coverage level (grouped bar): 1/2/3 kits x 3 scenarios
  3. AF scatter: method vs. WGS truth (2x3 grid)
  4. ACMG misclassification (grouped bar): 2 disease models x 3 scenarios
  5. AN inflation histogram: AN_naive/AN_AFQuery distribution x 3 scenarios

Inputs:
  - results/merged.parquet
  - results/acmg_results.json

Outputs:
  - figures/fig_capkit_*.{pdf,png}
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))
sys.path.insert(0, str(_BENCH_DIR.parent / "src"))

from shared.utils import (  # noqa: E402
    WONG_COLORS as COLORS,
    FIG_W_SINGLE as SINGLE_COL_WIDTH,
    FIG_W_DOUBLE as DOUBLE_COL_WIDTH,
    apply_style,
    save_figure,
)
from config import FIGURES_DIR, RESULTS_DIR, SCENARIOS, ensure_dirs  # noqa: E402

apply_style()

SCENARIO_ORDER = ["balanced", "skewed", "extreme"]
SCENARIO_LABELS = {
    "balanced": "Balanced\n(334/333/333)",
    "skewed": "Realistic skew\n(600/300/100)",
    "extreme": "Extreme skew\n(800/150/50)",
}


def _save(fig, name: str):
    save_figure(fig, name, FIGURES_DIR)


def _color_violin(parts, color):
    """Apply color to violin plot parts."""
    for pc in parts["bodies"]:
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
    for key in ("cmins", "cmaxes", "cbars", "cmedians"):
        if key in parts:
            parts[key].set_color(color)


# ---------------------------------------------------------------------------
# Figure 1: AF Error Distribution (violin plot)
# ---------------------------------------------------------------------------
def plot_af_error_violin(df):
    fig, ax = plt.subplots(figsize=(DOUBLE_COL_WIDTH, 3.5))

    centers = []
    for i, scenario in enumerate(SCENARIO_ORDER):
        sub = df[(df["scenario"] == scenario) & (df["AN_wes"] > 0)]

        pos_a = i * 3
        pos_n = i * 3 + 1

        err_afq = sub["error_afquery"].dropna()
        err_afq_nz = err_afq[err_afq > 0]
        if len(err_afq_nz) > 0:
            parts = ax.violinplot(
                [err_afq_nz.values], positions=[pos_a],
                showmedians=True, widths=0.8,
            )
            _color_violin(parts, COLORS["blue"])

        err_naive = sub["error_naive"].dropna()
        err_naive_nz = err_naive[err_naive > 0]
        if len(err_naive_nz) > 0:
            parts = ax.violinplot(
                [err_naive_nz.values], positions=[pos_n],
                showmedians=True, widths=0.8,
            )
            _color_violin(parts, COLORS["orange"])

        centers.append((pos_a + pos_n) / 2)

    ax.set_xticks(centers)
    ax.set_xticklabels([SCENARIO_LABELS[s] for s in SCENARIO_ORDER])
    ax.set_yscale("log")
    ax.set_ylabel("|AF error|  (|AF_method - AF_WGS|)")
    ax.set_title("AF Estimation Error: Technology-Aware vs. Naive AN")

    ax.legend(
        handles=[
            Patch(facecolor=COLORS["blue"], alpha=0.7,
                  label="AFQuery (tech-aware AN)"),
            Patch(facecolor=COLORS["orange"], alpha=0.7,
                  label="Naive (AN = 2N)"),
        ],
        frameon=False, loc="upper left",
    )

    _save(fig, "fig_capkit_af_error_violin")


# ---------------------------------------------------------------------------
# Figure 2: AF Error by Kit Coverage (grouped bar)
# ---------------------------------------------------------------------------
def plot_af_error_by_coverage(df):
    fig, axes = plt.subplots(
        1, 3, figsize=(DOUBLE_COL_WIDTH, 3.0), sharey=True,
    )

    for ax, scenario in zip(axes, SCENARIO_ORDER):
        sub = df[(df["scenario"] == scenario) & (df["AN_wes"] > 0)]
        bar_width = 0.35
        x = np.arange(3)

        mae_naive = []
        mae_afq = []
        for nk in [1, 2, 3]:
            s = sub[sub["n_kits"] == nk]
            mae_naive.append(s["error_naive"].mean() if len(s) > 0 else 0)
            mae_afq.append(s["error_afquery"].mean() if len(s) > 0 else 0)

        ax.bar(
            x - bar_width / 2, mae_naive, bar_width,
            label="Naive", color=COLORS["orange"], edgecolor="white",
        )
        ax.bar(
            x + bar_width / 2, mae_afq, bar_width,
            label="AFQuery", color=COLORS["blue"], edgecolor="white",
        )

        ax.set_xticks(x)
        ax.set_xticklabels(["1 kit", "2 kits", "3 kits"])
        ax.set_title(scenario.capitalize())
        if ax is axes[0]:
            ax.set_ylabel("Mean |AF error|")
            ax.legend(frameon=False, fontsize=7)

    fig.suptitle("AF Error by Capture Kit Coverage", y=1.02, fontsize=11)
    plt.tight_layout()
    _save(fig, "fig_capkit_error_by_coverage")


# ---------------------------------------------------------------------------
# Figure 3: AF Scatter (2x3 grid)
# ---------------------------------------------------------------------------
def plot_af_scatter(df):
    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE_COL_WIDTH, 5.0))
    methods = [("AF_afquery", "AFQuery"), ("AF_naive", "Naive (AN=2N)")]
    kit_colors = {
        1: COLORS["red"], 2: COLORS["yellow"], 3: COLORS["green"],
    }

    for row, (af_col, method_label) in enumerate(methods):
        for col, scenario in enumerate(SCENARIO_ORDER):
            ax = axes[row, col]
            sub = df[(df["scenario"] == scenario) & (df["AN_wes"] > 0)]

            # Plot 1-kit last (on top) to highlight discordant positions
            for nk in [3, 2, 1]:
                s = sub[sub["n_kits"] == nk]
                ax.scatter(
                    s["AF_wgs"], s[af_col], s=2, alpha=0.3,
                    color=kit_colors[nk], label=f"{nk} kit(s)",
                    rasterized=True,
                )

            ax.plot(
                [0, 1], [0, 1], "--", color="black",
                linewidth=0.5, alpha=0.5,
            )
            ax.set_xlim(-0.02, 0.52)
            ax.set_ylim(-0.02, 0.52)
            ax.set_aspect("equal")

            if row == 1:
                ax.set_xlabel("AF (WGS truth)")
            if col == 0:
                ax.set_ylabel(f"AF ({method_label})")
            if row == 0:
                ax.set_title(
                    SCENARIO_LABELS[scenario].replace("\n", " "),
                )

    axes[0, 2].legend(frameon=False, fontsize=7, markerscale=3, loc="lower right")
    plt.tight_layout()
    _save(fig, "fig_capkit_scatter")


# ---------------------------------------------------------------------------
# Figure 4: ACMG Misclassification (grouped bar)
# ---------------------------------------------------------------------------
def plot_acmg(acmg_data):
    fig, axes = plt.subplots(
        1, 2, figsize=(DOUBLE_COL_WIDTH, 3.5), sharey=True,
    )

    disease_titles = {
        "cardiomyopathy": "Cardiomyopathy\n(BS1>0.1%, PM2\u22640.01%)",
        "metabolic": "Metabolic\n(BS1>0.01%, PM2=absent)",
    }

    for ax, disease in zip(axes, ["cardiomyopathy", "metabolic"]):
        x = np.arange(len(SCENARIO_ORDER))
        bar_width = 0.35

        rates_naive = []
        rates_afq = []
        for scenario in SCENARIO_ORDER:
            d = acmg_data[disease][scenario]
            rates_naive.append(100 * d["naive"]["discordance_rate"])
            rates_afq.append(100 * d["afquery"]["discordance_rate"])

        bars_n = ax.bar(
            x - bar_width / 2, rates_naive, bar_width,
            label="Naive", color=COLORS["orange"], edgecolor="white",
        )
        ax.bar(
            x + bar_width / 2, rates_afq, bar_width,
            label="AFQuery", color=COLORS["blue"], edgecolor="white",
        )

        # Annotate critical errors on naive bars
        for i, scenario in enumerate(SCENARIO_ORDER):
            d = acmg_data[disease][scenario]
            n_fp = d["naive"]["false_pathogenic"]
            n_mp = d["naive"]["missed_pathogenic"]
            if n_fp > 0 or n_mp > 0:
                ax.annotate(
                    f"FP:{n_fp}\nMP:{n_mp}",
                    xy=(
                        bars_n[i].get_x() + bars_n[i].get_width() / 2,
                        bars_n[i].get_height(),
                    ),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=6,
                    color=COLORS["red"],
                )

        ax.set_xticks(x)
        ax.set_xticklabels([s.capitalize() for s in SCENARIO_ORDER])
        ax.set_title(disease_titles[disease], fontsize=9)
        if ax is axes[0]:
            ax.set_ylabel("Misclassified variants (%)")
        ax.legend(frameon=False, fontsize=7)

    fig.suptitle(
        "ACMG Classification Errors from AN Distortion",
        y=1.05, fontsize=11,
    )
    plt.tight_layout()
    _save(fig, "fig_capkit_acmg")


# ---------------------------------------------------------------------------
# Figure 5: AN Inflation Histogram
# ---------------------------------------------------------------------------
def plot_an_ratio(df):
    fig, axes = plt.subplots(
        1, 3, figsize=(DOUBLE_COL_WIDTH, 2.8), sharey=True,
    )

    for ax, scenario in zip(axes, SCENARIO_ORDER):
        sub = df[(df["scenario"] == scenario) & (df["AN_wes"] > 0)]
        ratio = sub["AN_ratio"].dropna()
        ratio = ratio[ratio > 0]

        ax.hist(
            ratio, bins=50, color=COLORS["blue"], alpha=0.7,
            edgecolor="white", linewidth=0.3,
        )
        ax.axvline(
            1.0, color=COLORS["green"], linestyle="--", linewidth=1.2,
            label="No inflation",
        )

        median_r = float(ratio.median())
        ax.axvline(
            median_r, color=COLORS["red"], linestyle=":", linewidth=0.8,
            label=f"Median: {median_r:.2f}\u00d7",
        )

        ax.set_xlabel("AN_naive / AN_AFQuery")
        ax.set_title(
            SCENARIO_LABELS[scenario].replace("\n", " "), fontsize=9,
        )
        if ax is axes[0]:
            ax.set_ylabel("Number of variants")
            ax.legend(frameon=False, fontsize=7)

    fig.suptitle("AN Inflation Distribution", y=1.02, fontsize=11)
    plt.tight_layout()
    _save(fig, "fig_capkit_an_ratio")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ensure_dirs()
    print("Generating capture kit benchmark figures...")

    merged_path = RESULTS_DIR / "merged.parquet"
    acmg_path = RESULTS_DIR / "acmg_results.json"

    if not merged_path.exists():
        print(f"ERROR: {merged_path} not found. Run 03_compute_metrics.py first.")
        sys.exit(1)

    df = pd.read_parquet(merged_path)
    print(f"Loaded {len(df)} rows from merged.parquet")

    print("\nFigure 1: AF error violin")
    plot_af_error_violin(df)

    print("\nFigure 2: AF error by kit coverage")
    plot_af_error_by_coverage(df)

    print("\nFigure 3: AF scatter")
    plot_af_scatter(df)

    if acmg_path.exists():
        acmg_data = json.loads(acmg_path.read_text())
        print("\nFigure 4: ACMG misclassification")
        plot_acmg(acmg_data)
    else:
        print("\nSkipping Figure 4: acmg_results.json not found")

    print("\nFigure 5: AN inflation histogram")
    plot_an_ratio(df)

    print(f"\nAll figures saved to {FIGURES_DIR}/")


if __name__ == "__main__":
    main()
