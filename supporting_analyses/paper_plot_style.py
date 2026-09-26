"""Shared Matplotlib style for project-generated quantitative figures.

Keep this module beside the analysis data and plotting scripts.  New
Matplotlib figures should import ``apply_paper_style`` and ``COLOURS`` instead
of introducing an incompatible local theme.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

# These must be set before importing Matplotlib for reliable headless runs.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-botryococcus"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLOURS = {
    "blue": "#0072B2",
    "orange": "#D55E00",
    "green": "#009E73",
    "purple": "#7c5aa6",
    "grid": "#d8d8d8",
    "axis": "#333333",
    "identity": "#555555",
}


def apply_paper_style() -> None:
    """Apply the white-background, light-grid manuscript plotting language."""
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.linewidth": 0.8,
            "axes.edgecolor": COLOURS["axis"],
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "legend.fontsize": 8.2,
            "legend.frameon": True,
            "legend.framealpha": 0.94,
            "legend.edgecolor": "#808080",
            "legend.fancybox": False,
            "mathtext.fontset": "dejavusans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def add_paper_grid(axis, *, minor: bool = True) -> None:
    """Add the shared unobtrusive grid to a Matplotlib axis."""
    axis.grid(True, which="major", color=COLOURS["grid"], linewidth=0.65, alpha=0.75)
    if minor:
        axis.grid(True, which="minor", color=COLOURS["grid"], linewidth=0.35, alpha=0.40)
