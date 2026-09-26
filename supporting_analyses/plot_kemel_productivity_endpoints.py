#!/usr/bin/env python3
"""Create a vector, color-independent plot of the two Kemel productivity endpoints."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from paper_plot_style import COLOURS, add_paper_grid, apply_paper_style, plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "Comparison_Literature" / "Comparison_Kemel" / "Kemel2025_comparison_package" / "kemel_2025_productivity_targets.csv"
DEFAULT_OUTPUT = HERE.parent / "Draft_2026-09-07" / "figures"

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    output_dir = parser.parse_args().output_dir
    apply_paper_style()
    with DATA.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    diameters = np.asarray([float(row["mean_diameter_um"]) for row in rows])
    figure, axis = plt.subplots(figsize=(6.25, 4.45), layout="constrained")
    for column, label, color, marker, linestyle in (
        ("biomass_productivity_g_L_d", "Biomass", COLOURS["blue"], "o", "-"),
        ("hydrocarbon_productivity_g_L_d", "Hydrocarbon", COLOURS["orange"], "s", "--"),
    ):
        values = np.asarray([float(row[column]) for row in rows])
        axis.plot(diameters, values, color=color, marker=marker, linestyle=linestyle,
                  linewidth=1.8, markersize=5.5, markeredgecolor="white",
                  markeredgewidth=0.7, label=label, zorder=3)
    axis.set_xlim(130, 450)
    axis.set_ylim(0.014, 0.128)
    axis.set_xticks([150, 200, 250, 300, 350, 400])
    axis.set_yticks(np.arange(0.02, 0.13, 0.02))
    axis.set_xlabel(r"Mean colony diameter ($\mu$m)")
    axis.set_ylabel(r"Volumetric productivity (g L$^{-1}$ d$^{-1}$)")
    axis.set_title("Kemel et al. (2025): size–productivity targets")
    axis.legend(loc="upper right", fontsize=9)
    add_paper_grid(axis, minor=False)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_dir / "fig4_kemel_productivity.pdf", bbox_inches="tight")
    figure.savefig(output_dir / "fig4_kemel_productivity.png", dpi=300, bbox_inches="tight")
    plt.close(figure)

if __name__ == "__main__":
    main()
