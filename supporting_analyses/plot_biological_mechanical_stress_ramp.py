#!/usr/bin/env python3
"""Plot the shared-style prospective biological-versus-mechanical test.

Run ``python3 CODES/biological_mechanical_stress_ramp_point6.py`` first.  The
figure displays the seeded protocol benchmark only; its generator curves are
shown solely to assess recovery and are not experimental observations.
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


def fitted_or_generator_parameter(metrics: list[dict[str, str]], metric: str, *, generator: bool = False) -> float:
    row = next(row for row in metrics if row["metric"] == metric)
    return float(row["reference_value"] if generator else row["value"])


def save_pair(figure, output_directory: Path, stem: str) -> tuple[Path, Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    pdf_path = output_directory / f"{stem}.pdf"
    png_path = output_directory / f"{stem}.png"
    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, bbox_inches="tight")
    plt.close(figure)
    return pdf_path, png_path


def hydrodynamic_rate(stress: np.ndarray, prefactor: float, exponent: float) -> np.ndarray:
    return prefactor * np.maximum(stress - 1.0, 0.0) ** exponent


def biological_rate(state: np.ndarray, maximum: float, center: float, width: float) -> np.ndarray:
    return maximum / (1.0 + np.exp(-(state - center) / width))


def make_figure(output_directory: Path) -> tuple[Path, Path]:
    apply_paper_style()
    exposure = read_rows("biological_mechanical_stress_ramp_exposure.csv")
    rates = read_rows("biological_mechanical_stress_ramp_rate_predictions.csv")
    kernels = read_rows("biological_mechanical_stress_ramp_kernel_density.csv")
    lognormal = read_rows("biological_mechanical_stress_ramp_lognormal_response.csv")
    metrics = read_rows("biological_mechanical_stress_ramp_metrics.csv")

    hyd_prefactor = fitted_or_generator_parameter(metrics, "beta_hyd_prefactor_d_1")
    hyd_exponent = fitted_or_generator_parameter(metrics, "hyd_exponent")
    hyd_prefactor_generator = fitted_or_generator_parameter(metrics, "beta_hyd_prefactor_d_1", generator=True)
    hyd_exponent_generator = fitted_or_generator_parameter(metrics, "hyd_exponent", generator=True)
    bio_maximum = fitted_or_generator_parameter(metrics, "beta_bio_max_d_1")
    bio_center = fitted_or_generator_parameter(metrics, "bio_maturation_center")
    bio_width = fitted_or_generator_parameter(metrics, "bio_maturation_width")
    bio_maximum_generator = fitted_or_generator_parameter(metrics, "beta_bio_max_d_1", generator=True)
    bio_center_generator = fitted_or_generator_parameter(metrics, "bio_maturation_center", generator=True)
    bio_width_generator = fitted_or_generator_parameter(metrics, "bio_maturation_width", generator=True)

    figure = plt.figure(figsize=(7.35, 6.15), layout="constrained")
    grid = figure.add_gridspec(2, 3)
    protocol_axis = figure.add_subplot(grid[0, 0])
    hyd_axis = figure.add_subplot(grid[0, 1])
    bio_axis = figure.add_subplot(grid[0, 2])
    kernel_axis = figure.add_subplot(grid[1, 0])
    center_axis = figure.add_subplot(grid[1, 1])
    width_axis = figure.add_subplot(grid[1, 2])

    ramp_styles = {
        "subthreshold_train": (COLOURS["green"], "Subthreshold train", "-"),
        "mild_up_train": (COLOURS["blue"], "Mild up train", "-"),
        "strong_up_train": (COLOURS["orange"], "Strong up train", "-"),
        "mid_up_test": (COLOURS["purple"], "Mid up test", "--"),
        "down_test": (COLOURS["axis"], "Down test", "--"),
    }
    for ramp_id, (colour, label, linestyle) in ramp_styles.items():
        selected = [row for row in exposure if row["ramp_id"] == ramp_id]
        selected.sort(key=lambda row: float(row["time_d"]))
        protocol_axis.plot(
            [float(row["time_d"]) for row in selected],
            [float(row["normalized_stress_tau_over_sigma_ref"]) for row in selected],
            color=colour,
            linestyle=linestyle,
            linewidth=1.35,
            label=label,
        )
    shared = [row for row in exposure if row["ramp_id"] == "subthreshold_train"]
    shared.sort(key=lambda row: float(row["time_d"]))
    maturity_axis = protocol_axis.twinx()
    maturity_axis.plot(
        [float(row["time_d"]) for row in shared],
        [float(row["maturation_state"]) for row in shared],
        color="#111111",
        linestyle=":",
        linewidth=1.35,
        label="Shared maturity",
    )
    protocol_axis.axhline(1.0, color=COLOURS["axis"], linestyle=":", linewidth=0.8)
    protocol_axis.set_xlim(0.0, 8.0)
    protocol_axis.set_ylim(0.55, 1.90)
    maturity_axis.set_ylim(0.0, 1.0)
    protocol_axis.set_xlabel("Time [d]")
    protocol_axis.set_ylabel(r"Stress ratio, $\tau_h/\sigma_{\rm ref}$")
    maturity_axis.set_ylabel(r"Maturation, $a$")
    protocol_axis.set_title("(a) Common-history stress programs", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    handles, labels = protocol_axis.get_legend_handles_labels()
    maturity_handles, maturity_labels = maturity_axis.get_legend_handles_labels()
    protocol_axis.legend(
        handles + maturity_handles,
        labels + maturity_labels,
        loc="upper left",
        fontsize=7.4,
        handlelength=1.45,
        borderpad=0.32,
        labelspacing=0.25,
    )
    add_paper_grid(protocol_axis, minor=False)

    rate_rows = [row for row in rates if int(float(row["at_risk_lineages"])) >= 100]
    for split, marker, face, label in (
        ("train", "o", COLOURS["blue"], "Training intervals"),
        ("test", "o", "white", "Held-out intervals"),
    ):
        selected = [row for row in rate_rows if row["split"] == split]
        hyd_axis.scatter(
            [float(row["normalized_stress_tau_over_sigma_ref"]) for row in selected],
            [float(row["observed_hydrodynamic_rate_d_1"]) for row in selected],
            s=9.5,
            marker=marker,
            facecolor=face,
            edgecolor=COLOURS["blue"],
            linewidth=0.45,
            alpha=0.62,
            label=label,
        )
    stress_grid = np.linspace(0.65, 1.85, 301)
    hyd_axis.plot(
        stress_grid,
        hydrodynamic_rate(stress_grid, hyd_prefactor, hyd_exponent),
        color=COLOURS["blue"],
        linewidth=1.75,
        label=r"Fitted $\beta_{\rm hyd}$",
    )
    hyd_axis.plot(
        stress_grid,
        hydrodynamic_rate(stress_grid, hyd_prefactor_generator, hyd_exponent_generator),
        color=COLOURS["axis"],
        linewidth=1.00,
        linestyle=":",
        label="Benchmark generator",
    )
    hyd_axis.axvline(1.0, color=COLOURS["axis"], linestyle=":", linewidth=0.75)
    hyd_axis.set_xlim(0.62, 1.88)
    hyd_axis.set_ylim(bottom=-0.01)
    hyd_axis.set_xlabel(r"Stress ratio, $\tau_h/\sigma_{\rm ref}$")
    hyd_axis.set_ylabel(r"Hydrodynamic rate [d$^{-1}$]")
    hyd_axis.set_title(r"(b) Stress-selected $\beta_{\rm hyd}$", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    hyd_axis.legend(loc="upper left", fontsize=7.2, handlelength=1.6, borderpad=0.4, labelspacing=0.4)
    add_paper_grid(hyd_axis, minor=False)

    for split, marker, face, label in (
        ("train", "o", COLOURS["orange"], "Training intervals"),
        ("test", "o", "white", "Held-out intervals"),
    ):
        selected = [row for row in rate_rows if row["split"] == split]
        bio_axis.scatter(
            [float(row["maturation_state"]) for row in selected],
            [float(row["observed_biological_rate_d_1"]) for row in selected],
            s=9.5,
            marker=marker,
            facecolor=face,
            edgecolor=COLOURS["orange"],
            linewidth=0.45,
            alpha=0.62,
            label=label,
        )
    state_grid = np.linspace(0.08, 0.94, 301)
    bio_axis.plot(
        state_grid,
        biological_rate(state_grid, bio_maximum, bio_center, bio_width),
        color=COLOURS["orange"],
        linewidth=1.75,
        label=r"Fitted $\beta_{\rm bio}$",
    )
    bio_axis.plot(
        state_grid,
        biological_rate(state_grid, bio_maximum_generator, bio_center_generator, bio_width_generator),
        color=COLOURS["axis"],
        linewidth=1.00,
        linestyle=":",
        label="Benchmark generator",
    )
    bio_axis.set_xlim(0.07, 0.95)
    bio_axis.set_ylim(bottom=-0.01)
    bio_axis.set_xlabel(r"Maturation state, $a$")
    bio_axis.set_ylabel(r"Biological rate [d$^{-1}$]")
    bio_axis.set_title(r"(c) Maturation-gated $\beta_{\rm bio}$", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    bio_axis.legend(loc="upper left", fontsize=7.2, handlelength=1.6, borderpad=0.4, labelspacing=0.4)
    add_paper_grid(bio_axis, minor=False)

    for channel, colour, label in (
        ("hydrodynamic", COLOURS["blue"], r"$b_{\rm hyd}$"),
        ("biological", COLOURS["orange"], r"$b_{\rm bio}$"),
    ):
        selected = [row for row in kernels if row["channel"] == channel]
        selected.sort(key=lambda row: float(row["material_fraction_of_parent"]))
        fractions = [float(row["material_fraction_of_parent"]) for row in selected]
        kernel_axis.plot(
            fractions,
            [float(row["fitted_daughter_kernel_density"]) for row in selected],
            color=colour,
            linewidth=1.65,
            label=f"Fitted {label}",
            linestyle="-" if channel == "hydrodynamic" else "--",
        )
        kernel_axis.plot(
            fractions,
            [float(row["generator_daughter_kernel_density"]) for row in selected],
            color=colour,
            linewidth=0.95,
            linestyle=":" if channel == "hydrodynamic" else "-.",
            label=f"Generator {label}",
        )
    kernel_axis.set_xlim(0.0, 1.0)
    kernel_axis.set_ylim(bottom=0.0)
    kernel_axis.set_xlabel(r"Daughter material fraction, $z=x_d/x_m$")
    kernel_axis.set_ylabel("Expected daughter density")
    kernel_axis.set_title("(d) Separately observed kernels", loc="left", fontsize=9.2, fontweight="bold", pad=10)
    kernel_axis.legend(loc="upper center", fontsize=7.2, handlelength=1.5, ncol=2, borderpad=0.4, labelspacing=0.4)
    add_paper_grid(kernel_axis, minor=False)

    response_stress = np.asarray([float(row["normalized_stress_tau_over_sigma_ref"]) for row in lognormal])
    fitted_median = np.asarray([float(row["fitted_median_diameter_um"]) for row in lognormal])
    generator_median = np.asarray([float(row["generator_median_diameter_um"]) for row in lognormal])
    fitted_width = np.asarray([float(row["fitted_log_standard_deviation_s"]) for row in lognormal])
    generator_width = np.asarray([float(row["generator_log_standard_deviation_s"]) for row in lognormal])
    center_axis.plot(response_stress, fitted_median, color=COLOURS["purple"], linewidth=1.75, label="Recovered event fit")
    center_axis.plot(response_stress, generator_median, color=COLOURS["axis"], linewidth=1.00, linestyle=":", label="Benchmark generator")
    center_axis.axvline(1.0, color=COLOURS["axis"], linestyle=":", linewidth=0.75)
    center_axis.set_xlim(0.62, 1.88)
    center_axis.set_xlabel(r"Stress ratio, $\tau_h/\sigma_{\rm ref}$")
    center_axis.set_ylabel(r"Log-normal median [$\mu$m]")
    center_axis.set_title(r"(e) Predicted center, $\exp(m)$", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    center_axis.text(0.04, 0.05, r"Fixed $a=0.86$", transform=center_axis.transAxes, fontsize=8.0)
    center_axis.legend(loc="upper right", fontsize=7.2, handlelength=1.6, borderpad=0.4, labelspacing=0.4)
    add_paper_grid(center_axis, minor=False)

    width_axis.plot(response_stress, fitted_width, color=COLOURS["green"], linewidth=1.75, label="Recovered event fit")
    width_axis.plot(response_stress, generator_width, color=COLOURS["axis"], linewidth=1.00, linestyle=":", label="Benchmark generator")
    width_axis.axvline(1.0, color=COLOURS["axis"], linestyle=":", linewidth=0.75)
    width_axis.set_xlim(0.62, 1.88)
    width_lower = min(float(np.min(fitted_width)), float(np.min(generator_width))) - 0.006
    width_upper = max(float(np.max(fitted_width)), float(np.max(generator_width))) + 0.006
    width_axis.set_ylim(width_lower, width_upper)
    width_axis.set_xlabel(r"Stress ratio, $\tau_h/\sigma_{\rm ref}$")
    width_axis.set_ylabel(r"Log-normal width, $s$")
    width_axis.set_title("(f) Predicted width", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    width_axis.text(0.04, 0.05, r"Fixed $a=0.86$", transform=width_axis.transAxes, fontsize=8.0)
    width_axis.legend(loc="upper right", fontsize=7.2, handlelength=1.6, borderpad=0.4, labelspacing=0.4)
    add_paper_grid(width_axis, minor=False)

    return save_pair(figure, output_directory, "fig12_biological_mechanical_stress_ramp")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIRECTORY)
    args = parser.parse_args()
    pdf_path, png_path = make_figure(args.output_dir)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
