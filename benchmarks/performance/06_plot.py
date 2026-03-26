#!/usr/bin/env python3
"""Generate publication-quality figures from benchmark results.

Reads JSON files from results/ and produces PDF + PNG figures in figures/.

Figures:
  1. Query latency vs. sample count (line plot, log-scale Y)
  2. Build time vs. threads (grouped bar chart)
  3. Annotation throughput vs. threads (line plot)
  4. AFQuery vs. bcftools comparison (grouped bar chart)
  5. AF concordance scatter plot
  6. Disk footprint (stacked bar)
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

_BENCH_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_BENCH_DIR))

from shared.utils import (  # noqa: E402
    WONG_COLORS as COLORS,
    FIG_W_SINGLE as SINGLE_COL_WIDTH,
    FIG_W_DOUBLE as DOUBLE_COL_WIDTH,
    apply_style,
    save_figure,
)
from config import FIGURES_DIR, RESULTS_DIR, ensure_dirs  # noqa: E402

apply_style()


def _save(fig, name: str):
    save_figure(fig, name, FIGURES_DIR)


# ---------------------------------------------------------------------------
# Figure 1: Query latency scaling
# ---------------------------------------------------------------------------
def plot_query_scaling():
    path = RESULTS_DIR / "query_scaling.json"
    if not path.exists():
        print("  Skipping: query_scaling.json not found")
        return

    data = json.loads(path.read_text())
    filter_name = "no_filter"  # primary plot uses no filter

    query_types = [
        ("point_warm", "Point (warm)", COLORS["blue"]),
        ("point_cold", "Point (cold)", COLORS["cyan"]),
        ("batch_100", "Batch 100", COLORS["orange"]),
        ("batch_1000", "Batch 1000", COLORS["red"]),
        ("region_1mbp", "Region 1 Mbp", COLORS["green"]),
    ]

    fig, ax = plt.subplots(figsize=(SINGLE_COL_WIDTH, 3.0))

    for qtype, label, color in query_types:
        xs = []
        ys = []
        y_lo = []
        y_hi = []
        for entry in data:
            if filter_name not in entry.get("filters", {}):
                continue
            fdata = entry["filters"][filter_name]
            if qtype not in fdata:
                continue
            xs.append(entry["n_samples"])
            ys.append(fdata[qtype]["median_ms"])
            y_lo.append(fdata[qtype]["q1_ms"])
            y_hi.append(fdata[qtype]["q3_ms"])

        if xs:
            ax.plot(xs, ys, "o-", label=label, color=color, markersize=4, linewidth=1.5)
            ax.fill_between(xs, y_lo, y_hi, alpha=0.15, color=color)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Latency (ms)")
    ax.set_title("Query Latency Scaling")
    ax.legend(loc="upper left", frameon=False)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

    _save(fig, "fig1_query_scaling")


# ---------------------------------------------------------------------------
# Figure 2: Build time vs. threads
# ---------------------------------------------------------------------------
def plot_build_perf():
    path = RESULTS_DIR / "build_perf.json"
    if not path.exists():
        print("  Skipping: build_perf.json not found")
        return

    data = json.loads(path.read_text())

    # Group by n_samples
    scales = sorted(set(r["n_samples"] for r in data))
    threads = sorted(set(r["threads"] for r in data))

    fig, ax = plt.subplots(figsize=(DOUBLE_COL_WIDTH, 3.5))

    n_groups = len(threads)
    n_bars = len(scales)
    bar_width = 0.8 / n_bars
    scale_colors = [COLORS["blue"], COLORS["orange"], COLORS["green"]]

    for i, scale in enumerate(scales):
        positions = [j + i * bar_width for j in range(n_groups)]
        heights = []
        for t in threads:
            entry = next(
                (r for r in data if r["n_samples"] == scale and r["threads"] == t),
                None,
            )
            heights.append(entry["wall_s"] if entry and "wall_s" in entry else 0)

        color = scale_colors[i % len(scale_colors)]
        ax.bar(
            positions, heights, bar_width,
            label=f"{scale:,} samples",
            color=color, edgecolor="white", linewidth=0.5,
        )

    ax.set_xticks([j + bar_width * (n_bars - 1) / 2 for j in range(n_groups)])
    ax.set_xticklabels([str(t) for t in threads])
    ax.set_xlabel("Build threads")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("Database Build Time")
    ax.legend(frameon=False)

    _save(fig, "fig2_build_perf")


# ---------------------------------------------------------------------------
# Figure 3: Annotation throughput
# ---------------------------------------------------------------------------
def plot_annotate_throughput():
    path = RESULTS_DIR / "annotate_throughput.json"
    if not path.exists():
        print("  Skipping: annotate_throughput.json not found")
        return

    data = json.loads(path.read_text())
    variant_counts = sorted(set(r["n_variants"] for r in data))
    var_colors = [COLORS["blue"], COLORS["orange"], COLORS["green"]]

    fig, ax = plt.subplots(figsize=(SINGLE_COL_WIDTH, 3.0))

    for i, nv in enumerate(variant_counts):
        subset = [r for r in data if r["n_variants"] == nv]
        subset.sort(key=lambda r: r["threads"])
        xs = [r["threads"] for r in subset]
        ys = [r["throughput_vars_per_s"] for r in subset]
        color = var_colors[i % len(var_colors)]
        ax.plot(xs, ys, "o-", label=f"{nv:,} variants", color=color, markersize=4, linewidth=1.5)

    ax.set_xlabel("Threads")
    ax.set_ylabel("Throughput (variants/s)")
    ax.set_title("VCF Annotation Throughput")
    ax.legend(frameon=False)

    _save(fig, "fig3_annotate_throughput")


# ---------------------------------------------------------------------------
# Figure 4: AFQuery vs. bcftools
# ---------------------------------------------------------------------------
def plot_bcftools_comparison():
    path = RESULTS_DIR / "bcftools_comparison.json"
    if not path.exists():
        print("  Skipping: bcftools_comparison.json not found")
        return

    data = json.loads(path.read_text())
    # Use largest sample count
    entry = max(data, key=lambda r: r["n_samples"])

    operations = ["point_query", "subset_query", "dump"]
    op_labels = ["Point query", "Subset query\n(by sex)", "Full dump\n(chr22)"]

    fig, ax = plt.subplots(figsize=(SINGLE_COL_WIDTH, 3.0))

    bar_width = 0.35
    x_pos = range(len(operations))

    bcf_vals = []
    afq_vals = []
    for op in operations:
        if op in entry:
            bcf_vals.append(entry[op]["bcftools"]["median_ms"])
            afq_vals.append(entry[op]["afquery"]["median_ms"])
        else:
            bcf_vals.append(0)
            afq_vals.append(0)

    bars1 = ax.bar(
        [x - bar_width / 2 for x in x_pos], bcf_vals, bar_width,
        label="bcftools", color=COLORS["orange"], edgecolor="white",
    )
    bars2 = ax.bar(
        [x + bar_width / 2 for x in x_pos], afq_vals, bar_width,
        label="AFQuery", color=COLORS["blue"], edgecolor="white",
    )

    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(op_labels)
    ax.set_ylabel("Time (ms)")
    ax.set_yscale("log")
    ax.set_title(f"AFQuery vs. bcftools ({entry['n_samples']:,} samples)")
    ax.legend(frameon=False)

    # Add value labels on bars
    for bar in list(bars1) + list(bars2):
        h = bar.get_height()
        if h > 0:
            ax.annotate(
                f"{h:.0f}",
                xy=(bar.get_x() + bar.get_width() / 2, h),
                xytext=(0, 3), textcoords="offset points",
                ha="center", va="bottom", fontsize=7,
            )

    _save(fig, "fig4_bcftools_comparison")


# ---------------------------------------------------------------------------
# Figure 5: AF concordance
# ---------------------------------------------------------------------------
def plot_concordance():
    path = RESULTS_DIR / "concordance.json"
    if not path.exists():
        print("  Skipping: concordance.json not found")
        return

    data = json.loads(path.read_text())
    pairs = data.get("af_pairs_sample", [])
    r_sq = data.get("r_squared")

    if not pairs:
        print("  Skipping: no AF pairs in concordance data")
        return

    fig, ax = plt.subplots(figsize=(SINGLE_COL_WIDTH, 3.2))

    afq_afs = [p[0] for p in pairs]
    bcf_afs = [p[1] for p in pairs]

    ax.scatter(bcf_afs, afq_afs, s=8, alpha=0.5, color=COLORS["blue"], edgecolors="none")
    ax.plot([0, 1], [0, 1], "--", color=COLORS["black"], linewidth=0.8, alpha=0.5)

    ax.set_xlabel("bcftools AF")
    ax.set_ylabel("AFQuery AF")
    ax.set_title("Allele Frequency Concordance")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal")

    if r_sq is not None:
        ax.text(
            0.05, 0.92, f"R² = {r_sq:.6f}\nn = {data['n_common_variants']:,}",
            transform=ax.transAxes, fontsize=9,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

    _save(fig, "fig5_concordance")


# ---------------------------------------------------------------------------
# Figure 6: Disk footprint
# ---------------------------------------------------------------------------
def plot_disk_footprint():
    path = RESULTS_DIR / "build_perf.json"
    if not path.exists():
        print("  Skipping: build_perf.json not found")
        return

    data = json.loads(path.read_text())

    # Get one entry per scale (last thread count has the DB still on disk)
    scales = sorted(set(r["n_samples"] for r in data))
    max_threads = max(r["threads"] for r in data)

    fig, ax = plt.subplots(figsize=(SINGLE_COL_WIDTH, 3.0))

    bar_width = 0.35
    x_pos = range(len(scales))

    raw_sizes = []
    db_sizes = []
    for scale in scales:
        entry = next(
            (r for r in data if r["n_samples"] == scale and r["threads"] == max_threads),
            None,
        )
        raw_sizes.append(entry.get("raw_vcf_size_mb", 0) if entry else 0)
        db_sizes.append(entry.get("db_size_mb", 0) if entry else 0)

    ax.bar(
        [x - bar_width / 2 for x in x_pos], raw_sizes, bar_width,
        label="Raw VCFs", color=COLORS["orange"], edgecolor="white",
    )
    ax.bar(
        [x + bar_width / 2 for x in x_pos], db_sizes, bar_width,
        label="AFQuery DB", color=COLORS["blue"], edgecolor="white",
    )

    ax.set_xticks(list(x_pos))
    ax.set_xticklabels([f"{s:,}" for s in scales])
    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Size (MB)")
    ax.set_title("Disk Footprint")
    ax.legend(frameon=False)

    _save(fig, "fig6_disk_footprint")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ensure_dirs()
    print("Generating figures...")

    print("\nFigure 1: Query latency scaling")
    plot_query_scaling()

    print("\nFigure 2: Build performance")
    plot_build_perf()

    print("\nFigure 3: Annotation throughput")
    plot_annotate_throughput()

    print("\nFigure 4: bcftools comparison")
    plot_bcftools_comparison()

    print("\nFigure 5: AF concordance")
    plot_concordance()

    print("\nFigure 6: Disk footprint")
    plot_disk_footprint()

    print(f"\nAll figures saved to {FIGURES_DIR}/")


if __name__ == "__main__":
    main()
