#!/usr/bin/env python3
"""Create the shared-style figure for the memory and bimodality test.

Run ``python3 CODES/memory_bimodality_point3.py`` first.  This script only
reads its reproducible CSV output; it never fits a distribution or alters the
model calculations.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from paper_plot_style import COLOURS, add_paper_grid, apply_paper_style, plt


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
DEFAULT_FIGURE_DIRECTORY = HERE.parent / "Draft_2026-09-07" / "figures"


def read_rows(filename: str) -> list[dict[str, str]]:
    with (RESULTS / filename).open(newline="") as handle:
        return list(csv.DictReader(handle))


def select_distribution(rows: list[dict[str, str]], model: str, days: float) -> tuple[np.ndarray, np.ndarray]:
    selected = [row for row in rows if row["model"] == model and np.isclose(float(row["days"]), days)]
    selected.sort(key=lambda row: float(row["diameter_um"]))
    return (
        np.asarray([float(row["diameter_um"]) for row in selected]),
        np.asarray([float(row["volume_percent"]) for row in selected]),
    )


def select_timecourse(rows: list[dict[str, str]], model: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    selected = [row for row in rows if row["model"] == model]
    selected.sort(key=lambda row: float(row["days"]))
    return (
        np.asarray([float(row["days"]) for row in selected]),
        np.asarray([float(row["small_volume_fraction"]) for row in selected]),
        np.asarray([float(row["large_residual_fraction"]) for row in selected]),
    )


def metric_row(rows: list[dict[str, str]], model: str) -> dict[str, str]:
    return next(row for row in rows if row["model"] == model)


def smooth_log_binned_curve(values: np.ndarray) -> np.ndarray:
    """Use the same fixed light smoothing as the declared modal diagnostic.

    The CSV files retain unsmoothed section-to-source-bin values.  This helper
    only prevents the display of a finite-volume interpolation saw-tooth as if
    it were a biological mode, and is applied uniformly to all four curves.
    """

    kernel = np.asarray([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0
    result = values.copy()
    for _ in range(2):
        result = np.convolve(np.pad(result, (2, 2), mode="edge"), kernel, mode="valid")
    return result


def save_pair(figure, output_directory: Path, stem: str) -> tuple[Path, Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    pdf_path = output_directory / f"{stem}.pdf"
    png_path = output_directory / f"{stem}.png"
    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, bbox_inches="tight")
    plt.close(figure)
    return pdf_path, png_path


def make_figure(output_directory: Path) -> tuple[Path, Path]:
    apply_paper_style()
    distributions = read_rows("memory_bimodality_distributions.csv")
    timecourse = read_rows("memory_bimodality_timecourse.csv")
    metrics = read_rows("memory_bimodality_metrics.csv")
    sensitivity = read_rows("memory_bimodality_sensitivity.csv")

    # A full-width distribution panel above two supporting diagnostics keeps
    # labels legible after the vector graphic is placed at manuscript-column
    # width.  A three-across strip would reduce the 9-pt shared style below a
    # useful print size.
    figure = plt.figure(figsize=(7.25, 5.10), layout="constrained")
    grid = figure.add_gridspec(2, 2, height_ratios=(1.02, 0.98))
    distribution_axis = figure.add_subplot(grid[0, :])
    trajectory_axis = figure.add_subplot(grid[1, 0])
    selection_axis = figure.add_subplot(grid[1, 1])

    initial_x, initial_y = select_distribution(distributions, "medium_light_initial", 0.0)
    target_x, target_y = select_distribution(distributions, "low_light_reconstruction_mean", 20.0)
    one_x, one_y = select_distribution(distributions, "one_state_collapsed_limit", 20.0)
    structured_x, structured_y = select_distribution(distributions, "structured_memory_biological_kernel", 20.0)

    target_y_raw = target_y.copy()
    initial_y = smooth_log_binned_curve(initial_y)
    target_y = smooth_log_binned_curve(target_y)
    one_y = smooth_log_binned_curve(one_y)
    structured_y = smooth_log_binned_curve(structured_y)
    distribution_axis.plot(initial_x, initial_y, color="#777777", linestyle="--", linewidth=1.25, label="Medium-light start")
    distribution_axis.plot(target_x, target_y, color=COLOURS["axis"], linewidth=2.0, label="LL day-20 reconstruction")
    distribution_axis.plot(one_x, one_y, color=COLOURS["orange"], linewidth=1.65, label="One-state limit")
    distribution_axis.plot(structured_x, structured_y, color=COLOURS["blue"], linewidth=1.85, label="Structured model")
    distribution_axis.set_xscale("log")
    distribution_axis.set_xlim(20.0, 2000.0)
    distribution_axis.set_ylim(0.0, max(target_y.max(), one_y.max(), structured_y.max()) * 1.14)
    distribution_axis.set_xticks([20, 50, 100, 200, 500, 1000, 2000])
    distribution_axis.set_xticklabels(["20", "50", "100", "200", "500", "1000", "2000"])
    distribution_axis.set_xlabel(r"Equivalent diameter [$\mu$m]")
    distribution_axis.set_ylabel("Volume distribution [%]")
    distribution_axis.text(0.03, 0.97, "(a) Same start, day 20", transform=distribution_axis.transAxes, va="top", fontweight="bold")
    add_paper_grid(distribution_axis, minor=False)
    distribution_axis.legend(loc="lower right", fontsize=6.7, handlelength=1.9)

    model_styles = (
        ("one_state_collapsed_limit", COLOURS["orange"], "One-state"),
        ("structured_memory_biological_kernel", COLOURS["blue"], "Structured"),
    )
    for model, colour, label in model_styles:
        time, small, large = select_timecourse(timecourse, model)
        trajectory_axis.plot(time, small, color=colour, marker="o", markersize=3.3, linewidth=1.5, label=f"{label}: small")
        trajectory_axis.plot(time, large, color=colour, marker="s", markersize=3.1, linewidth=1.35, linestyle="--", label=f"{label}: large")
    target_small = float(target_y_raw[target_x < 200.0].sum() / 100.0)
    target_large = float(target_y_raw[target_x >= 300.0].sum() / 100.0)
    trajectory_axis.axhline(target_small, color=COLOURS["axis"], linewidth=0.85, linestyle=":")
    trajectory_axis.axhline(target_large, color=COLOURS["axis"], linewidth=0.85, linestyle=":")
    trajectory_axis.text(19.7, target_small + 0.012, "LL small", ha="right", va="bottom", fontsize=7.0, color=COLOURS["axis"])
    trajectory_axis.text(19.7, target_large + 0.012, "LL large", ha="right", va="bottom", fontsize=7.0, color=COLOURS["axis"])
    trajectory_axis.set_xlim(0.0, 20.0)
    trajectory_axis.set_ylim(0.0, 0.78)
    trajectory_axis.set_xlabel("Days after low-light shift")
    trajectory_axis.set_ylabel("Volume fraction")
    trajectory_axis.text(0.03, 0.97, "(b) Mode and residual", transform=trajectory_axis.transAxes, va="top", fontweight="bold")
    add_paper_grid(trajectory_axis, minor=False)
    trajectory_axis.legend(loc="lower center", ncols=2, fontsize=6.5, columnspacing=0.8, handlelength=1.4)

    for row in sensitivity:
        selection_axis.scatter(
            float(row["small_volume_fraction"]),
            float(row["large_residual_fraction"]),
            s=37,
            marker="o",
            color=COLOURS["green"],
            edgecolor=COLOURS["axis"],
            linewidth=0.40,
            alpha=0.88,
            label="Structured sensitivity" if row is sensitivity[0] else None,
            zorder=2,
        )
    for model, colour, marker, label in (
        ("one_state_collapsed_limit", COLOURS["orange"], "D", "One-state reference"),
        ("structured_memory_biological_kernel", COLOURS["blue"], "o", "Structured reference"),
    ):
        row = metric_row(metrics, model)
        selection_axis.scatter(
            float(row["small_volume_fraction"]),
            float(row["large_residual_fraction"]),
            s=55,
            marker=marker,
            color=colour,
            edgecolor=COLOURS["axis"],
            linewidth=0.60,
            label=label,
            zorder=4,
        )
    selection_axis.scatter(
        target_small,
        target_large,
        s=57,
        marker="*",
        color=COLOURS["purple"],
        edgecolor=COLOURS["axis"],
        linewidth=0.50,
        label="LL reconstruction",
        zorder=5,
    )
    selection_axis.axvline(0.10, color=COLOURS["axis"], linewidth=0.75, linestyle=":")
    selection_axis.axhline(0.10, color=COLOURS["axis"], linewidth=0.75, linestyle=":")
    selection_axis.set_xlim(0.0, 0.72)
    selection_axis.set_ylim(0.0, 0.78)
    selection_axis.set_xlabel(r"Small fraction, $D<200\,\mu$m")
    selection_axis.set_ylabel(r"Large residual, $D\geq300\,\mu$m")
    selection_axis.text(0.03, 0.97, "(c) Fixed closure sensitivity", transform=selection_axis.transAxes, va="top", fontweight="bold")
    add_paper_grid(selection_axis, minor=False)
    selection_axis.legend(loc="lower right", fontsize=6.5, handlelength=1.35)

    return save_pair(figure, output_directory, "fig9_memory_bimodality_model_selection")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIRECTORY)
    args = parser.parse_args()
    pdf_path, png_path = make_figure(args.output_dir)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
