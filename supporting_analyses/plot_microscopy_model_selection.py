#!/usr/bin/env python3
"""Create the shared-style figure for the microscopy-constrained test.

Run ``python3 CODES/microscopy_model_selection_point5.py`` first.  The figure
only visualises its prospective, non-fitted structural-comparison output.
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
THRESHOLD = 3.0


def read_rows(filename: str) -> list[dict[str, str]]:
    with (RESULTS / filename).open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def save_pair(figure, output_directory: Path, stem: str) -> tuple[Path, Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    pdf_path = output_directory / f"{stem}.pdf"
    png_path = output_directory / f"{stem}.png"
    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, bbox_inches="tight")
    plt.close(figure)
    return pdf_path, png_path


def selected_candidate(candidates: list[dict[str, str]], kernel: str, reset_rule: str = "deep_reset") -> str:
    return next(
        row["candidate_id"]
        for row in candidates
        if row["kernel"] == kernel and row["shell_alternative"] == "no_resolved_shell" and row["reset_rule"] == reset_rule
    )


def make_figure(output_directory: Path) -> tuple[Path, Path]:
    apply_paper_style()
    candidates = read_rows("microscopy_model_selection_candidates.csv")
    lineages = read_rows("microscopy_model_selection_lineages.csv")
    summaries = read_rows("microscopy_model_selection_lineage_summary.csv")
    distributions = read_rows("microscopy_model_selection_distributions.csv")
    pairwise = read_rows("microscopy_model_selection_pairwise.csv")

    equal = selected_candidate(candidates, "equal_binary")
    broad = selected_candidate(candidates, "broad_binary")
    retained = selected_candidate(candidates, "retained_branch_satellites")
    deep = selected_candidate(candidates, "equal_binary", "deep_reset")
    partial = selected_candidate(candidates, "equal_binary", "partial_reset")
    kernel_styles = (
        (equal, COLOURS["orange"], "Equal binary"),
        (broad, COLOURS["blue"], "Broad binary"),
        (retained, COLOURS["green"], "Retained branch + satellites"),
    )

    figure = plt.figure(figsize=(7.25, 5.25), layout="constrained")
    grid = figure.add_gridspec(2, 2, height_ratios=(1.0, 1.0))
    ratio_axis = figure.add_subplot(grid[0, 0])
    state_axis = figure.add_subplot(grid[0, 1])
    distribution_axis = figure.add_subplot(grid[1, 0])
    protocol_axis = figure.add_subplot(grid[1, 1])

    ratio_edges = np.linspace(0.05, 1.05, 24)
    for candidate_id, colour, label in kernel_styles:
        ratios = np.asarray(
            [float(row["descendant_parent_diameter_ratio"]) for row in lineages if row["candidate_id"] == candidate_id]
        )
        histogram, edges = np.histogram(ratios, bins=ratio_edges, density=True)
        centers = np.sqrt(edges[:-1] * edges[1:])
        ratio_axis.plot(centers, histogram, color=colour, linewidth=1.75, drawstyle="steps-mid", label=label)
    ratio_axis.axvline(2.0 ** (-1.0 / 3.0), color=COLOURS["axis"], linewidth=0.80, linestyle=":")
    ratio_axis.text(0.795, 0.94, r"$2^{-1/3}$", transform=ratio_axis.get_xaxis_transform(), ha="center", va="top", fontsize=8.0)
    ratio_axis.set_xlim(0.05, 1.02)
    ratio_axis.set_xlabel(r"Descendant/parent diameter, $D_d/D_m$")
    ratio_axis.set_ylabel("Density")
    ratio_axis.set_title("(a) Paired lineage kernel signature", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    ratio_axis.legend(loc="upper left", fontsize=7.2, handlelength=1.65)
    add_paper_grid(ratio_axis, minor=False)

    summary_by_id = {row["candidate_id"]: row for row in summaries}
    time = np.asarray((0.0, 0.5, 1.0))
    for candidate_id, colour, linestyle, label in (
        (deep, COLOURS["purple"], "-", "Deep reset"),
        (partial, COLOURS["blue"], "--", "Partial reset"),
    ):
        row = summary_by_id[candidate_id]
        state = np.asarray([float(row["state_t0"]), float(row["state_t0_5"]), float(row["state_t1"])])
        state_axis.plot(time, state, color=colour, linestyle=linestyle, marker="o", markersize=3.8, linewidth=1.65, label=label)
    state_axis.set_xlim(-0.04, 1.04)
    state_axis.set_ylim(0.0, 1.02)
    state_axis.set_xticks(time)
    state_axis.set_xlabel("Cycles after recorded event")
    state_axis.set_ylabel("Reduced daughter state")
    state_axis.set_title("(b) Time-resolved reset signature", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    state_axis.legend(loc="lower right", fontsize=7.2, handlelength=1.65)
    add_paper_grid(state_axis, minor=False)

    for candidate_id, colour, label in kernel_styles:
        selected = [
            row
            for row in distributions
            if row["candidate_id"] == candidate_id
            and np.isclose(float(row["time_after_recorded_event_pulse_cycles"]), 0.75)
        ]
        selected.sort(key=lambda row: float(row["diameter_um"]))
        distribution_axis.plot(
            [float(row["diameter_um"]) for row in selected],
            [float(row["normalised_number_fraction"]) for row in selected],
            color=colour,
            linewidth=1.70,
            label=label,
        )
    distribution_axis.set_xscale("log")
    distribution_axis.set_xlim(40.0, 1000.0)
    distribution_axis.set_ylim(bottom=0.0)
    distribution_axis.set_xlabel(r"Equivalent diameter [$\mu$m]")
    distribution_axis.set_ylabel("Normalised number fraction")
    distribution_axis.set_title("(c) Distribution after known event pulses", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    distribution_axis.legend(loc="upper left", fontsize=7.2, handlelength=1.65)
    add_paper_grid(distribution_axis, minor=False)

    modalities = (
        ("Time-resolved\ndistributions", "distribution_discriminable", COLOURS["orange"]),
        ("Paired\nlineage", "paired_lineage_discriminable", COLOURS["blue"]),
        ("Shell\nlabel", "shell_label_discriminable", COLOURS["purple"]),
        ("State\ntrajectory", "state_trajectory_discriminable", COLOURS["green"]),
        ("Combined", "combined_discriminable", COLOURS["axis"]),
    )
    total_pairs = len(pairwise)
    counts = np.asarray([sum(as_bool(row[column]) for row in pairwise) for _, column, _ in modalities])
    bars = protocol_axis.bar(
        np.arange(len(modalities)),
        counts,
        width=0.72,
        color=[colour for _, _, colour in modalities],
        edgecolor="white",
        linewidth=0.6,
    )
    for bar, count in zip(bars, counts):
        protocol_axis.text(bar.get_x() + bar.get_width() / 2.0, count + 1.3, f"{count}/{total_pairs}", ha="center", va="bottom", fontsize=8.0)
    protocol_axis.set_xticks(np.arange(len(modalities)))
    protocol_axis.set_xticklabels([label for label, _, _ in modalities], fontsize=8.0)
    protocol_axis.set_ylim(0.0, total_pairs + 8.0)
    protocol_axis.set_ylabel("Candidate pairs separated")
    protocol_axis.set_title("(d) Prospective protocol adequacy", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    protocol_axis.text(
        0.02,
        0.95,
        f"Threshold: {THRESHOLD:.0f} s.u.\nNo rate fit",
        transform=protocol_axis.transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
        color=COLOURS["axis"],
    )
    add_paper_grid(protocol_axis, minor=False)

    return save_pair(figure, output_directory, "fig11_microscopy_model_selection")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIRECTORY)
    args = parser.parse_args()
    pdf_path, png_path = make_figure(args.output_dir)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
