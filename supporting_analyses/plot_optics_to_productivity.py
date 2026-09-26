#!/usr/bin/env python3
"""Plot the source-bound Point-7 optics-to-productivity test.

Run python3 CODES/optics_to_productivity_point7.py first. The conditional
optical curves in panel (a) are not fitted observations: they show why the
unreported optical diameter/state pairing prevents a unique a_abs(D, a)
relation. Panels (c) and (d) display the no-refit shared-yield invariant and
the product-allocation contrast that follows from the reported endpoints.
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
    """Read one Point-7 result table."""

    with (RESULTS / filename).open(newline="") as handle:
        return list(csv.DictReader(handle))


def source_metric(rows: list[dict[str, str]], metric: str) -> float:
    """Read a stable source-audit scalar."""

    return float(next(row["value"] for row in rows if row["metric"] == metric))


def ratio_metric(rows: list[dict[str, str]], quantity: str) -> float:
    """Read a stable ratio-test scalar."""

    return float(next(row["value"] for row in rows if row["quantity"] == quantity))


def save_pair(figure, output_directory: Path, stem: str) -> tuple[Path, Path]:
    """Write vector manuscript art plus a high-resolution preview."""

    output_directory.mkdir(parents=True, exist_ok=True)
    pdf_path = output_directory / f"{stem}.pdf"
    png_path = output_directory / f"{stem}.png"
    figure.savefig(pdf_path, bbox_inches="tight")
    figure.savefig(png_path, bbox_inches="tight")
    plt.close(figure)
    return pdf_path, png_path


def make_figure(output_directory: Path) -> tuple[Path, Path]:
    """Create the four-panel Point-7 evidence and identifiability figure."""

    apply_paper_style()
    source = read_rows("optics_to_productivity_source_audit.csv")
    curves = read_rows("optics_to_productivity_conditional_curves.csv")
    exponents = read_rows("optics_to_productivity_effective_exponent.csv")
    ratios = read_rows("optics_to_productivity_ratio_test.csv")

    intact_absorption = source_metric(source, "intact_optical_mass_absorption_cross_section")
    disrupted_absorption = source_metric(source, "disrupted_optical_mass_absorption_cross_section")
    medium_diameter = source_metric(source, "medium_productivity_culture_mean_diameter")
    large_diameter = source_metric(source, "large_productivity_culture_mean_diameter")
    biomass_ratio = ratio_metric(ratios, "observed_medium_to_large_biomass_productivity_ratio")
    hydrocarbon_ratio = ratio_metric(ratios, "observed_medium_to_large_hydrocarbon_productivity_ratio")
    allocation_contrast = ratio_metric(ratios, "required_medium_to_large_hydrocarbon_allocation_contrast")
    conditional_gamma = ratio_metric(
        ratios, "conditional_same_state_optical_exponent_if_145_and_434_are_paired"
    )

    figure = plt.figure(figsize=(7.35, 5.65), layout="constrained")
    grid = figure.add_gridspec(2, 2)
    optical_axis = figure.add_subplot(grid[0, 0])
    exponent_axis = figure.add_subplot(grid[0, 1])
    ratio_axis = figure.add_subplot(grid[1, 0])
    allocation_axis = figure.add_subplot(grid[1, 1])

    curve_colours = (COLOURS["blue"], COLOURS["orange"], COLOURS["purple"])
    curve_styles = ("-", "--", "-.")
    hypothetical_diameters = sorted({float(row["hypothetical_disrupted_diameter_um"]) for row in curves})
    for colour, linestyle, hypothetical_diameter in zip(curve_colours, curve_styles, hypothetical_diameters):
        selected = [
            row for row in curves if np.isclose(float(row["hypothetical_disrupted_diameter_um"]), hypothetical_diameter)
        ]
        selected.sort(key=lambda row: float(row["diameter_um"]))
        optical_axis.plot(
            [float(row["diameter_um"]) for row in selected],
            [float(row["conditional_mass_absorption_cross_section_m2_kg"]) for row in selected],
            color=colour,
            linestyle=linestyle,
            linewidth=1.55,
            label=rf"$D_s={hypothetical_diameter:.0f}\,\mu$m",
        )

    optical_axis.axhline(intact_absorption, color=COLOURS["axis"], linestyle=":", linewidth=0.75)
    optical_axis.axhline(disrupted_absorption, color=COLOURS["axis"], linestyle=":", linewidth=0.75)
    optical_axis.axvline(medium_diameter, color=COLOURS["orange"], linestyle="--", linewidth=0.9)
    optical_axis.axvline(large_diameter, color=COLOURS["blue"], linestyle="-.", linewidth=0.9)
    optical_axis.text(
        medium_diameter,
        4.2,
        "145",
        color=COLOURS["orange"],
        ha="center",
        va="bottom",
        fontsize=8.0,
    )
    optical_axis.text(
        large_diameter,
        4.2,
        "434",
        color=COLOURS["blue"],
        ha="center",
        va="bottom",
        fontsize=8.0,
    )
    optical_axis.set_xscale("log")
    optical_axis.set_xlim(38.0, 500.0)
    optical_axis.set_ylim(0.0, 88.0)
    optical_axis.set_xlabel(r"Diameter used only for conditional sensitivity [$\mu$m]")
    optical_axis.set_ylabel(r"$a_{\rm abs}$ [m$^2$ kg$^{-1}$]")
    optical_axis.set_title("(a) Conditional optical curves, not fits", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    optical_axis.legend(loc="lower left", fontsize=7.6, handlelength=1.65, borderpad=0.4, labelspacing=0.4)
    add_paper_grid(optical_axis, minor=False)

    exponent_diameter = np.asarray([float(row["hypothetical_disrupted_diameter_um"]) for row in exponents])
    exponent_values = np.asarray([float(row["conditional_effective_exponent_gamma"]) for row in exponents])
    exponent_axis.plot(exponent_diameter, exponent_values, color=COLOURS["purple"], linewidth=1.65)
    exponent_axis.axvline(medium_diameter, color=COLOURS["orange"], linestyle="--", linewidth=0.75)
    exponent_axis.scatter(
        [medium_diameter],
        [conditional_gamma],
        s=24,
        color=COLOURS["orange"],
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
    )
    exponent_axis.text(
        0.03,
        0.95,
        rf"If $D_s=145\,\mu$m: $\gamma={conditional_gamma:.2f}$"
        "\nconditional only",
        transform=exponent_axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        color=COLOURS["axis"],
    )
    exponent_axis.set_xscale("log")
    exponent_axis.set_xlim(18.0, 430.0)
    exponent_axis.set_ylim(0.0, max(6.0, float(np.nanpercentile(exponent_values, 99.0)) * 1.08))
    exponent_axis.set_xlabel(r"Unreported disrupted diameter, $D_s$ [$\mu$m]")
    exponent_axis.set_ylabel(r"Conditional $\gamma$ in $a_{\rm abs}\propto D^{-\gamma}$")
    exponent_axis.set_title("(b) One contrast cannot identify $D$ dependence", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    add_paper_grid(exponent_axis, minor=False)

    diagonal = np.linspace(1.0, 2.02, 100)
    ratio_axis.plot(
        diagonal,
        diagonal,
        color=COLOURS["axis"],
        linestyle=":",
        linewidth=1.0,
        label=r"All common-yield optics predictions, $R_{\rm HC}=R_X$",
    )
    ratio_axis.scatter(
        [biomass_ratio],
        [hydrocarbon_ratio],
        s=48,
        color=COLOURS["green"],
        edgecolor="white",
        linewidth=0.65,
        label="Reported productivity target",
        zorder=4,
    )
    ratio_axis.annotate(
        "",
        xy=(biomass_ratio, hydrocarbon_ratio),
        xytext=(biomass_ratio, biomass_ratio),
        arrowprops={"arrowstyle": "<->", "color": COLOURS["green"], "lw": 1.0},
    )
    ratio_axis.text(
        biomass_ratio + 0.035,
        0.5 * (biomass_ratio + hydrocarbon_ratio),
        rf"$\Delta R={hydrocarbon_ratio - biomass_ratio:.2f}$",
        color=COLOURS["green"],
        fontsize=7.2,
        va="center",
    )
    ratio_axis.set_xlim(0.96, 2.02)
    ratio_axis.set_ylim(0.96, 2.02)
    ratio_axis.set_aspect("equal", adjustable="box")
    ratio_axis.set_xlabel(r"Biomass ratio, $R_X=P_{X,m}/P_{X,l}$")
    ratio_axis.set_ylabel(r"Hydrocarbon ratio, $R_{\rm HC}$")
    ratio_axis.set_title("(c) Fixed-yield optical null is rejected", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    add_paper_grid(ratio_axis, minor=False)

    medium_share = source_metric(source, "medium_hydrocarbon_to_biomass_productivity_fraction")
    large_share = source_metric(source, "large_hydrocarbon_to_biomass_productivity_fraction")
    bars = allocation_axis.bar(
        (0, 1),
        (large_share, medium_share),
        width=0.58,
        color=(COLOURS["blue"], COLOURS["green"]),
        edgecolor="white",
        linewidth=0.65,
    )
    for bar, value in zip(bars, (large_share, medium_share)):
        allocation_axis.text(
            bar.get_x() + bar.get_width() / 2.0,
            value + 0.010,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=7.2,
        )
    allocation_axis.annotate(
        "",
        xy=(1.0, medium_share + 0.002),
        xytext=(0.0, large_share + 0.002),
        arrowprops={"arrowstyle": "<->", "color": COLOURS["purple"], "lw": 1.0},
    )
    allocation_axis.set_xticks((0, 1))
    allocation_axis.set_xticklabels(("Large\n434 um", "Medium\n145 um"))
    allocation_axis.set_xlim(-0.62, 1.62)
    allocation_axis.set_ylim(0.0, max(medium_share, large_share) + 0.11)
    allocation_axis.set_ylabel(r"HC/biomass productivity fraction")
    allocation_axis.set_title("(d) Required allocation contrast", loc="left", fontsize=9.2, fontweight="bold", pad=4)
    allocation_axis.text(
        0.03,
        0.94,
        "Requires state-dependent\nhydrocarbon allocation",
        transform=allocation_axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        color=COLOURS["axis"],
    )
    allocation_axis.text(
        0.97,
        0.94,
        "Medium/large = $" + f"{allocation_contrast:.2f}" + r"\times$",
        transform=allocation_axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.0,
        color=COLOURS["purple"],
    )
    add_paper_grid(allocation_axis, minor=False)

    return save_pair(figure, output_directory, "fig13_optics_to_productivity_test")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIRECTORY)
    args = parser.parse_args()
    pdf_path, png_path = make_figure(args.output_dir)
    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
