#!/usr/bin/env python3
"""Plot the Point-4 García-Cubero environmental/chemostat test.

Run ``python3 CODES/garcia_cubero_chemostat_point4.py`` first.  This script
only reads the resulting CSV files; it neither refits the state closure nor
changes the declared calibration/prediction split.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from paper_plot_style import COLOURS, add_paper_grid, apply_paper_style, plt
from matplotlib.lines import Line2D


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
DEFAULT_FIGURE_DIRECTORY = HERE.parent / "Draft_2026-09-07" / "figures"

STATE_MODEL = "environmental_state_closure"
NULL_MODEL = "uniform_washout_null"
SPLIT_STYLE = {
    "calibration": ("o", "Calibration subset"),
    "unused_design_check": ("s", "Unused design check"),
    "optimized_prediction": ("*", "Held-out optimized condition"),
}


def read_rows(filename: str) -> list[dict[str, str]]:
    with (RESULTS / filename).open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str) -> float | None:
    return None if value == "" else float(value)


def metric(rows: list[dict[str, str]], model: str, split: str) -> dict[str, str]:
    return next(row for row in rows if row["model"] == model and row["split"] == split)


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
    predictions = read_rows("garcia_chemostat_predictions.csv")
    metrics = read_rows("garcia_chemostat_metrics.csv")
    dilution = read_rows("garcia_chemostat_dilution_response.csv")
    washout = read_rows("garcia_chemostat_uniform_washout.csv")
    invariance = read_rows("garcia_chemostat_washout_invariance.csv")

    # Match the preprint text width closely so the shared 9-pt plotting style
    # remains legible after inclusion at ``\linewidth`` rather than being
    # aggressively downscaled from a landscape-sized source graphic.
    figure = plt.figure(figsize=(5.55, 4.50), layout="constrained")
    grid = figure.add_gridspec(2, 2, height_ratios=(1.04, 0.96))
    parity_axis = figure.add_subplot(grid[0, :])
    dilution_axis = figure.add_subplot(grid[1, 0])
    washout_axis = figure.add_subplot(grid[1, 1])

    state_rows = [row for row in predictions if row["model"] == STATE_MODEL]
    null_rows = {row["condition_id"]: row for row in predictions if row["model"] == NULL_MODEL}
    for row in state_rows:
        marker, _ = SPLIT_STYLE[row["split"]]
        observed = float(row["observed_volume_weighted_diameter_um"])
        predicted = float(row["predicted_volume_weighted_diameter_um"])
        point_size = 70 if marker == "*" else 35
        sd = as_float(row["observed_sd_um"])
        if sd is not None:
            parity_axis.errorbar(
                observed,
                predicted,
                xerr=sd,
                color=COLOURS["blue"],
                linewidth=0.65,
                capsize=1.5,
                zorder=2,
            )
        parity_axis.scatter(
            observed,
            predicted,
            s=point_size,
            marker=marker,
            color=COLOURS["blue"],
            edgecolor="white",
            linewidth=0.55,
            zorder=3,
        )
        null_row = null_rows[row["condition_id"]]
        parity_axis.scatter(
            observed,
            float(null_row["predicted_volume_weighted_diameter_um"]),
            s=27,
            marker="x",
            color=COLOURS["orange"],
            linewidth=1.10,
            zorder=4,
        )

    limits = (35.0, 375.0)
    parity_axis.plot(limits, limits, color=COLOURS["axis"], linewidth=0.85, linestyle=":", zorder=1)
    parity_axis.set_xlim(*limits)
    parity_axis.set_ylim(*limits)
    parity_axis.set_xlabel(r"Observed volume-weighted $D_v$ [$\mu$m]")
    parity_axis.set_ylabel(r"Predicted $D_v$ [$\mu$m]")
    parity_axis.set_title("(a) Fixed calibration/prediction split", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    parity_handles = [
        Line2D([], [], marker="o", linestyle="", color=COLOURS["blue"], markeredgecolor="white", label="Environmental state"),
        Line2D([], [], marker="x", linestyle="", color=COLOURS["orange"], label="Uniform-washout null"),
        # The split is encoded by marker shape on the blue state-closure
        # points; these handles deliberately match those plotted markers.
        Line2D([], [], marker="o", linestyle="", color=COLOURS["blue"], markeredgecolor="white", label="Calibration subset"),
        Line2D([], [], marker="s", linestyle="", color=COLOURS["blue"], markeredgecolor="white", label="Unused design check"),
        Line2D([], [], marker="*", linestyle="", color=COLOURS["blue"], markeredgecolor="white", markersize=8.4, label="Held-out optimization"),
    ]
    # The lower-left corner contains an unused design-check point. Keep the key in
    # the otherwise empty upper-left region so every plotted marker remains
    # visible in the final manuscript figure.
    parity_axis.legend(
        handles=parity_handles,
        loc="upper left",
        ncols=2,
        fontsize=7.2,
        title="Colour = model; marker = data split",
        title_fontsize=7.2,
        columnspacing=0.8,
        handletextpad=0.45,
    )
    add_paper_grid(parity_axis, minor=False)

    for model, colour, linestyle, label in (
        (STATE_MODEL, COLOURS["blue"], "-", r"Environmental state, $q_*$"),
        (NULL_MODEL, COLOURS["orange"], "--", "Uniform washout only"),
    ):
        selected = [row for row in dilution if row["model"] == model]
        selected.sort(key=lambda row: float(row["dilution_d_1"]))
        dilution_values = np.asarray([float(row["dilution_d_1"]) for row in selected])
        diameter_values = np.asarray([float(row["predicted_volume_weighted_diameter_um"]) for row in selected])
        dilution_axis.plot(dilution_values, diameter_values, color=colour, linestyle=linestyle, marker="o", markersize=3.5, linewidth=1.7, label=label)

    observed_rows = [row for row in dilution if row["model"] == STATE_MODEL and row["observed_volume_weighted_diameter_um"] != ""]
    dilution_axis.scatter(
        [float(row["dilution_d_1"]) for row in observed_rows],
        [float(row["observed_volume_weighted_diameter_um"]) for row in observed_rows],
        marker="D",
        s=36,
        color=COLOURS["purple"],
        edgecolor="white",
        linewidth=0.50,
        label="Source-study observation",
        zorder=4,
    )
    dilution_axis.set_xlim(0.085, 0.315)
    dilution_axis.set_xticks([0.10, 0.20, 0.30])
    dilution_axis.set_ylim(75.0, 310.0)
    dilution_axis.set_xlabel(r"Dilution rate [$\mathrm{d}^{-1}$]")
    dilution_axis.set_ylabel(r"Volume-weighted $D_v$ [$\mu$m]")
    dilution_axis.set_title(r"(b) Matched $T=22.5^\circ$C, $I=1800$", loc="left", fontsize=8.6, fontweight="bold", pad=3)
    dilution_axis.text(0.97, 0.08, r"State route: $1.99\times$", transform=dilution_axis.transAxes, ha="right", va="bottom", fontsize=8.0)
    dilution_axis.legend(loc="upper left", fontsize=7.2, handlelength=1.8)
    add_paper_grid(dilution_axis, minor=False)

    for washout_rate, colour, linestyle, label in (
        (0.10, COLOURS["axis"], "-", r"$-0.10n$"),
        (0.30, COLOURS["orange"], "--", r"$-0.30n$"),
    ):
        selected = [row for row in washout if np.isclose(float(row["washout_d_1"]), washout_rate)]
        selected.sort(key=lambda row: float(row["diameter_um"]))
        washout_axis.plot(
            [float(row["diameter_um"]) for row in selected],
            [float(row["normalized_volume_fraction"]) for row in selected],
            color=colour,
            linestyle=linestyle,
            linewidth=1.55,
            label=label,
        )
    maximum_tv = max(float(row["maximum_normalized_distribution_total_variation"]) for row in invariance)
    maximum_relative_diameter_change = max(float(row["relative_volume_weighted_diameter_difference"]) for row in invariance)
    washout_axis.set_xscale("log")
    washout_axis.set_xlim(35.0, 700.0)
    washout_axis.set_xlabel(r"Equivalent diameter [$\mu$m]")
    washout_axis.set_ylabel("Normalised volume fraction")
    washout_axis.set_title("(c) Uniform-washout control", loc="left", fontsize=8.6, fontweight="bold", pad=3)
    washout_axis.text(
        0.97,
        0.08,
        rf"$d_{{\rm TV}}={maximum_tv:.1e}$" + "\n" + rf"relative $\Delta D_v={maximum_relative_diameter_change:.1e}$",
        transform=washout_axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.0,
    )
    washout_axis.legend(loc="upper right", fontsize=7.2, handlelength=1.8)
    add_paper_grid(washout_axis, minor=False)

    return save_pair(figure, output_directory, "fig10_garcia_chemostat_state_test")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIRECTORY)
    args = parser.parse_args()
    pdf_path, png_path = make_figure(args.output_dir)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
