#!/usr/bin/env python3
"""Create shared-style figures for the Zhang--Kojima test hierarchy.

Run ``python3 CODES/zhang_hierarchy_point2.py`` first.  This plotting script
uses its frozen-parameter Fig. 5/Fig. 3 results and the separate Fig. 6
diagnostic records; it does not perform any fit itself.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

from paper_plot_style import COLOURS, add_paper_grid, apply_paper_style, plt


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
DEFAULT_FIGURE_DIRECTORY = HERE.parent / "Draft_2026-09-07" / "figures"


def read_rows(filename: str) -> list[dict[str, str]]:
    with (RESULTS / filename).open(newline="") as handle:
        return list(csv.DictReader(handle))


def save_pair(figure, output_directory: Path, stem: str) -> tuple[Path, Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    pdf_path = output_directory / f"{stem}.pdf"
    png_path = output_directory / f"{stem}.png"
    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, bbox_inches="tight")
    plt.close(figure)
    return pdf_path, png_path


def plot_hierarchy(output_directory: Path) -> tuple[Path, Path]:
    """Plot Fig. 5 scalar validation beside the no-refit distribution test."""
    apply_paper_style()
    fig5 = read_rows("zhang_fig5_predictions.csv")
    distribution = read_rows("zhang_fig3_distribution_metrics.csv")
    figure, axes = plt.subplots(1, 2, figsize=(7.1, 3.15), layout="constrained")

    parity = axes[0]
    parity.plot([0.035, 0.25], [0.035, 0.25], "--", color=COLOURS["identity"], linewidth=1.15, label="Identity")
    split_style = {
        "calibration": (COLOURS["blue"], "o", "3-klx calibration"),
        "prediction": (COLOURS["orange"], "D", "10-klx no-refit prediction"),
    }
    for split, (colour, marker, label) in split_style.items():
        selected = [row for row in fig5 if row["split"] == split]
        parity.scatter(
            [float(row["observed_dv_mm"]) for row in selected],
            [float(row["predicted_dv_mm"]) for row in selected],
            s=28,
            marker=marker,
            color=colour,
            edgecolor=COLOURS["axis"],
            linewidth=0.55,
            label=label,
            zorder=3,
        )
    parity.set_xlim(0.035, 0.25)
    parity.set_ylim(0.035, 0.25)
    parity.set_aspect("equal", adjustable="box")
    parity.set_xlabel(r"Observed $D_v$ [mm]")
    parity.set_ylabel(r"Predicted $D_v$ [mm]")
    parity.text(0.03, 0.97, "(a) Fig. 5 transient", transform=parity.transAxes, va="top", fontweight="bold")
    add_paper_grid(parity, minor=False)
    parity.legend(loc="lower right", fontsize=7.1, handlelength=1.8)

    distance_axis = axes[1]
    groups = [
        (3, COLOURS["blue"], "3 klx\nno-refit"),
        (10, COLOURS["orange"], "10 klx\nindependent no-refit"),
    ]
    for xpos, (preculture, colour, _label) in enumerate(groups):
        values = np.asarray(
            [float(row["total_variation_distance"]) for row in distribution if int(row["preculture_irradiance_klx"]) == preculture]
        )
        jitter = np.linspace(-0.09, 0.09, len(values))
        distance_axis.scatter(
            xpos + jitter,
            values,
            s=26,
            color=colour,
            edgecolor=COLOURS["axis"],
            linewidth=0.45,
            alpha=0.90,
            zorder=3,
            label="Individual distribution" if xpos == 0 else None,
        )
        mean_value = float(values.mean())
        distance_axis.hlines(mean_value, xpos - 0.15, xpos + 0.15, color=COLOURS["axis"], linewidth=2.0, zorder=4)
        distance_axis.plot(xpos, mean_value, marker="_", markersize=14, color=COLOURS["axis"], label="Mean" if xpos == 0 else None)
    distance_axis.set_xlim(-0.38, 1.38)
    distance_axis.set_ylim(0.0, 1.0)
    distance_axis.set_xticks([0, 1], [group[2] for group in groups])
    distance_axis.set_ylabel(r"Distribution total variation, $d_{\rm TV}$")
    distance_axis.text(0.03, 0.97, "(b) Fig. 3 distribution prediction", transform=distance_axis.transAxes, va="top", fontweight="bold")
    add_paper_grid(distance_axis, minor=False)
    distance_axis.legend(loc="lower right", fontsize=7.2)

    return save_pair(figure, output_directory, "fig7_zhang_hierarchy_validation")


def plot_fig6_diagnostic(output_directory: Path) -> tuple[Path, Path]:
    """Plot the quasi-steady Fig. 6 diagnostic without pooling it with Fig. 5."""
    apply_paper_style()
    rows = read_rows("zhang_fig6_diagnostic_points.csv")
    figure, axis = plt.subplots(figsize=(6.2, 4.05), layout="constrained")
    groups = [
        (3, COLOURS["blue"], "o", "3-klx"),
        (10, COLOURS["orange"], "D", "10-klx"),
    ]
    for preculture, colour, marker, label in groups:
        selected = [row for row in rows if int(row["preculture_irradiance_klx"]) == preculture]
        axis.scatter(
            [float(row["average_light_intensity_V"]) for row in selected],
            [float(row["observed_dv_mm"]) for row in selected],
            s=28,
            marker=marker,
            color=colour,
            edgecolor=COLOURS["axis"],
            linewidth=0.50,
            alpha=0.88,
            label=f"{label} observations",
            zorder=3,
        )
        by_light: dict[float, list[float]] = defaultdict(list)
        for row in selected:
            by_light[float(row["average_light_intensity_V"])].append(float(row["predicted_dv_mm"]))
        xvalues = sorted(by_light)
        axis.plot(
            xvalues,
            [np.mean(by_light[value]) for value in xvalues],
            color=colour,
            linewidth=1.65,
            label=f"{label} quasi-steady PBE",
            zorder=2,
        )
    axis.set_xlim(0.2, 4.0)
    axis.set_ylim(0.03, 0.26)
    axis.set_xlabel(r"Average bubble-column light signal, $I_{\rm BC}$ [V]")
    axis.set_ylabel(r"Volume-averaged diameter, $D_v$ [mm]")
    add_paper_grid(axis)
    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.19),
        ncols=2,
        fontsize=7.2,
        columnspacing=1.0,
        handlelength=1.8,
    )
    return save_pair(figure, output_directory, "fig8_zhang_fig6_quasisteady_diagnostic")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIRECTORY)
    args = parser.parse_args()
    for output in (plot_hierarchy(args.output_dir), plot_fig6_diagnostic(args.output_dir)):
        print(f"Wrote {output[0]}")
        print(f"Wrote {output[1]}")


if __name__ == "__main__":
    main()
