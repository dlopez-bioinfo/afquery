"""Shared utilities for all benchmarks."""

import statistics
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Publication style
# ---------------------------------------------------------------------------

# Colorblind-safe palette (Wong, 2011 — Nature Methods 8:441)
WONG_COLORS = {
    "blue":   "#0072B2",
    "orange": "#E69F00",
    "green":  "#009E73",
    "red":    "#D55E00",
    "purple": "#CC79A7",
    "cyan":   "#56B4E9",
    "yellow": "#F0E442",
    "black":  "#000000",
}

FIG_W_SINGLE = 3.5  # inches (journal single column)
FIG_W_DOUBLE = 7.0  # inches (journal double column)

RCPARAMS = {
    "font.size": 10,
    "font.family": "sans-serif",
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
}

# ---------------------------------------------------------------------------
# Timing utilities
# ---------------------------------------------------------------------------


def time_ms(func, *args, **kwargs):
    """Run func(*args, **kwargs) and return (result, elapsed_ms)."""
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    elapsed = (time.perf_counter() - t0) * 1000
    return result, elapsed


def stats(times: list) -> dict:
    """Compute summary statistics for a list of times (ms)."""
    s = sorted(times)
    n = len(s)
    return {
        "median_ms": statistics.median(s),
        "q1_ms": s[n // 4],
        "q3_ms": s[(3 * n) // 4],
        "min_ms": s[0],
        "max_ms": s[-1],
    }


# ---------------------------------------------------------------------------
# Figure utilities
# ---------------------------------------------------------------------------


def apply_style():
    """Apply shared publication style to matplotlib."""
    plt.rcParams.update(RCPARAMS)


def save_figure(fig, name: str, figures_dir: Path):
    """Save figure as PDF and PNG in figures_dir."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures_dir / f"{name}.pdf", format="pdf")
    fig.savefig(figures_dir / f"{name}.png", format="png")
    plt.close(fig)
    print(f"  Saved {name}.pdf + {name}.png")


# ---------------------------------------------------------------------------
# Path utilities
# ---------------------------------------------------------------------------


def ensure_dirs(*paths: Path):
    """Create all given directories (and parents) if they do not exist."""
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)
