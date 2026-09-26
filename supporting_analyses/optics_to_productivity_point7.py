#!/usr/bin/env python3
"""Source-bound optics-to-productivity test for the Kemel et al. endpoints.

The accessible Kemel comparison package contains two unpaired evidence blocks:

1. intact and high-pressure-disrupted colonies have mass absorption
   cross-sections of 8.0 and 79 m2 kg-1, respectively, but the diameters and
   reduced physiological/ECM states of those optical samples are not reported;
2. medium (145 um) and large (434 um) cultures have biomass and hydrocarbon
   productivity endpoints, but their mass absorption cross-sections are not
   reported.

It would therefore be invalid to fit a continuous empirical a_abs(D, a)
relation from these scalar records. This script makes that identification
boundary quantitative, then applies a no-refit structural test to a
shared-optical-yield null. With colony distributions/dynamics held fixed, a
common absorbed-photon factor must multiply biomass and hydrocarbon
productivity equally. The null consequently predicts equal medium/large
ratios; it can be rejected without fitting an optical exponent or a
productivity coefficient.

For transparency, the script also reports a conditional sensitivity in which
the optical endpoints are hypothetically co-registered with the 145 and
434 um productivity states at common reduced state. That calculation is
deliberately marked as conditional throughout: it is not a source-derived
optical fit and does not turn missing paired measurements into data.

Run from the workspace root:

    python3 CODES/optics_to_productivity_point7.py
    python3 CODES/plot_optics_to_productivity.py
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
RESULTS = HERE / "results"
KEMEL_PACKAGE = PROJECT_ROOT / "Comparison_Literature" / "Comparison_Kemel" / "Kemel2025_comparison_package"
OPTICAL_TARGETS = KEMEL_PACKAGE / "kemel_2025_optical_targets.csv"
PRODUCTIVITY_TARGETS = KEMEL_PACKAGE / "kemel_2025_productivity_targets.csv"
SOURCE_README = KEMEL_PACKAGE / "README.txt"

# These values are not observations. They enumerate three visibly different
# ways in which the unreported disrupted optical diameter could be assigned
# when drawing a conditional same-state power law. The 145-um value is
# included because it is the medium-culture productivity diameter, not because
# the source associates it with the disrupted optical sample.
ILLUSTRATIVE_DISRUPTED_DIAMETERS_UM = (50.0, 145.0, 300.0)
EFFECTIVE_EXPONENT_DIAMETER_GRID_UM = np.geomspace(20.0, 420.0, 401)

# A minimally informative follow-up design for the log-linear optical
# diagnostic: log a_abs = c - gamma log(D/D0) + eta a. The values are design
# locations, not new Kemel measurements.
RECOMMENDED_DIAMETER_BINS_UM = (75.0, 145.0, 250.0, 434.0)
RECOMMENDED_STATE_STRATA = ("low measured state", "high measured state")


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a CSV table while preserving the source text fields."""

    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    """Write a nonempty, consistently ordered CSV result table."""

    rows = list(rows)
    if not rows:
        raise ValueError(f"cannot write an empty result table: {path.name}")
    fieldnames = list(rows[0])
    if any(list(row) != fieldnames for row in rows):
        raise ValueError(f"inconsistent fields while writing {path.name}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def source_inputs() -> tuple[float, float, dict[str, str], dict[str, str]]:
    """Return the reported optical endpoints and medium/large cultures."""

    optical = read_rows(OPTICAL_TARGETS)
    productivity = read_rows(PRODUCTIVITY_TARGETS)
    if len(optical) != 2 or len(productivity) != 2:
        raise ValueError("the Kemel package must contain exactly two optical and two productivity endpoint rows")

    optical_by_state = {row["state"]: row for row in optical}
    intact_key = "large-colony state before disruption"
    disrupted_key = "smaller-colony state after high-pressure homogenization"
    if set(optical_by_state) != {intact_key, disrupted_key}:
        raise ValueError("unexpected optical endpoint labels in the Kemel package")
    intact = float(optical_by_state[intact_key]["mass_absorption_cross_section_m2_kg"])
    disrupted = float(optical_by_state[disrupted_key]["mass_absorption_cross_section_m2_kg"])
    if not (intact > 0.0 and disrupted > intact):
        raise ValueError("Kemel optical endpoints must be positive and increase after disruption")

    productivity_by_culture = {row["culture"]: row for row in productivity}
    if set(productivity_by_culture) != {"medium colonies", "large colonies"}:
        raise ValueError("unexpected productivity endpoint labels in the Kemel package")
    medium = productivity_by_culture["medium colonies"]
    large = productivity_by_culture["large colonies"]
    if float(medium["mean_diameter_um"]) >= float(large["mean_diameter_um"]):
        raise ValueError("Kemel medium culture must be smaller than large culture")
    return intact, disrupted, medium, large


def conditional_power_law_exponent(
    absorption_ratio: float, reference_diameter_um: float, hypothetical_disrupted_diameter_um: float
) -> float:
    """Return gamma for a conditional a_abs proportional to D^-gamma relation."""

    if not (0.0 < hypothetical_disrupted_diameter_um < reference_diameter_um):
        raise ValueError("the hypothetical disrupted diameter must be strictly below the reference diameter")
    return math.log(absorption_ratio) / math.log(reference_diameter_um / hypothetical_disrupted_diameter_um)


def saturation_scale_for_target_ratio(
    target_ratio: float, absorption_large: float, absorption_small: float
) -> float:
    """Invert Phi(a)=a/(K+a) for one conditionally paired target ratio.

    This is only a conditional diagnostic. The source does not pair absorption
    values with the productivity cultures, and the resulting scale is not
    fitted or interpreted as a physical reactor parameter.
    """

    absorption_ratio = absorption_small / absorption_large
    if not (1.0 < target_ratio < absorption_ratio):
        raise ValueError("target ratio must lie between one and the conditional absorption contrast")
    return absorption_small * (target_ratio - 1.0) / (absorption_ratio - target_ratio)


def make_source_audit_rows(
    absorption_large: float,
    absorption_small: float,
    medium: dict[str, str],
    large: dict[str, str],
) -> list[dict[str, object]]:
    """Create source-supported quantities and explicitly label their pairing."""

    biomass_ratio = float(medium["biomass_productivity_g_L_d"]) / float(large["biomass_productivity_g_L_d"])
    hydrocarbon_ratio = float(medium["hydrocarbon_productivity_g_L_d"]) / float(
        large["hydrocarbon_productivity_g_L_d"]
    )
    medium_hydrocarbon_share = float(medium["hydrocarbon_productivity_g_L_d"]) / float(
        medium["biomass_productivity_g_L_d"]
    )
    large_hydrocarbon_share = float(large["hydrocarbon_productivity_g_L_d"]) / float(
        large["biomass_productivity_g_L_d"]
    )
    allocation_contrast = medium_hydrocarbon_share / large_hydrocarbon_share

    return [
        {
            "metric": "intact_optical_mass_absorption_cross_section",
            "value": absorption_large,
            "units": "m2 kg-1",
            "evidence_block": "optical",
            "source_pairing_status": "diameter and reduced state not reported with this optical endpoint",
            "interpretation": "reported intact large-colony optical endpoint",
        },
        {
            "metric": "disrupted_optical_mass_absorption_cross_section",
            "value": absorption_small,
            "units": "m2 kg-1",
            "evidence_block": "optical",
            "source_pairing_status": "diameter and reduced state not reported with this optical endpoint",
            "interpretation": "reported high-pressure-disrupted optical endpoint",
        },
        {
            "metric": "disrupted_to_intact_optical_contrast",
            "value": absorption_small / absorption_large,
            "units": "dimensionless",
            "evidence_block": "optical",
            "source_pairing_status": "reported optical contrast; not paired to the productivity cultures",
            "interpretation": "size/architecture sensitivity exists, but its D and a dependence is not identified",
        },
        {
            "metric": "medium_productivity_culture_mean_diameter",
            "value": float(medium["mean_diameter_um"]),
            "units": "um",
            "evidence_block": "productivity",
            "source_pairing_status": "no mass-absorption value reported for this culture endpoint",
            "interpretation": "reported medium-culture productivity diameter",
        },
        {
            "metric": "large_productivity_culture_mean_diameter",
            "value": float(large["mean_diameter_um"]),
            "units": "um",
            "evidence_block": "productivity",
            "source_pairing_status": "no mass-absorption value reported for this culture endpoint",
            "interpretation": "reported large-culture productivity diameter",
        },
        {
            "metric": "medium_to_large_biomass_productivity_ratio",
            "value": biomass_ratio,
            "units": "dimensionless",
            "evidence_block": "productivity",
            "source_pairing_status": "large productivity inferred from a rounded reported gain",
            "interpretation": "reported/inferred process endpoint",
        },
        {
            "metric": "medium_to_large_hydrocarbon_productivity_ratio",
            "value": hydrocarbon_ratio,
            "units": "dimensionless",
            "evidence_block": "productivity",
            "source_pairing_status": "large productivity inferred from a rounded reported gain",
            "interpretation": "reported/inferred process endpoint",
        },
        {
            "metric": "medium_hydrocarbon_to_biomass_productivity_fraction",
            "value": medium_hydrocarbon_share,
            "units": "g hydrocarbon per g biomass productivity",
            "evidence_block": "productivity",
            "source_pairing_status": "reported medium-culture productivity values",
            "interpretation": "product allocation observable",
        },
        {
            "metric": "large_hydrocarbon_to_biomass_productivity_fraction",
            "value": large_hydrocarbon_share,
            "units": "g hydrocarbon per g biomass productivity",
            "evidence_block": "productivity",
            "source_pairing_status": "large productivity inferred from rounded gains",
            "interpretation": "product allocation observable",
        },
        {
            "metric": "medium_to_large_hydrocarbon_allocation_contrast",
            "value": allocation_contrast,
            "units": "dimensionless",
            "evidence_block": "derived no-refit invariant",
            "source_pairing_status": "does not require optical co-registration",
            "interpretation": "minimum state/product-specific contrast required beyond a common optical multiplier",
        },
    ]


def make_conditional_optical_rows(
    absorption_large: float,
    absorption_small: float,
    reference_diameter_um: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Enumerate conditional curves and the effective-exponent continuum."""

    ratio = absorption_small / absorption_large
    curve_rows: list[dict[str, object]] = []
    for hypothetical_diameter in ILLUSTRATIVE_DISRUPTED_DIAMETERS_UM:
        exponent = conditional_power_law_exponent(ratio, reference_diameter_um, hypothetical_diameter)
        diameter_grid = np.geomspace(hypothetical_diameter, reference_diameter_um, 181)
        absorption_grid = absorption_large * (diameter_grid / reference_diameter_um) ** (-exponent)
        for diameter, absorption in zip(diameter_grid, absorption_grid):
            curve_rows.append(
                {
                    "curve_type": "conditional_same_state_power_law_not_a_fit",
                    "hypothetical_disrupted_diameter_um": hypothetical_diameter,
                    "reference_diameter_um": reference_diameter_um,
                    "conditional_effective_exponent_gamma": exponent,
                    "diameter_um": float(diameter),
                    "conditional_mass_absorption_cross_section_m2_kg": float(absorption),
                }
            )

    exponent_rows: list[dict[str, object]] = []
    for hypothetical_diameter in EFFECTIVE_EXPONENT_DIAMETER_GRID_UM:
        exponent_rows.append(
            {
                "curve_type": "conditional_same_state_power_law_not_a_fit",
                "hypothetical_disrupted_diameter_um": float(hypothetical_diameter),
                "reference_diameter_um": reference_diameter_um,
                "conditional_effective_exponent_gamma": conditional_power_law_exponent(
                    ratio, reference_diameter_um, float(hypothetical_diameter)
                ),
            }
        )
    return curve_rows, exponent_rows


def make_ratio_rows(
    absorption_large: float,
    absorption_small: float,
    medium: dict[str, str],
    large: dict[str, str],
) -> list[dict[str, object]]:
    """Report the no-refit invariant and conditional transfer sensitivity."""

    biomass_ratio = float(medium["biomass_productivity_g_L_d"]) / float(large["biomass_productivity_g_L_d"])
    hydrocarbon_ratio = float(medium["hydrocarbon_productivity_g_L_d"]) / float(
        large["hydrocarbon_productivity_g_L_d"]
    )
    absorption_ratio = absorption_small / absorption_large
    common_best_ratio = math.sqrt(biomass_ratio * hydrocarbon_ratio)
    scale_biomass = saturation_scale_for_target_ratio(biomass_ratio, absorption_large, absorption_small)
    scale_hydrocarbon = saturation_scale_for_target_ratio(hydrocarbon_ratio, absorption_large, absorption_small)
    scale_common = saturation_scale_for_target_ratio(common_best_ratio, absorption_large, absorption_small)
    allocation_contrast = hydrocarbon_ratio / biomass_ratio

    rows = [
        {
            "quantity": "observed_medium_to_large_biomass_productivity_ratio",
            "value": biomass_ratio,
            "units": "dimensionless",
            "test_status": "source target",
            "interpretation": "reported/inferred Kemel endpoint",
        },
        {
            "quantity": "observed_medium_to_large_hydrocarbon_productivity_ratio",
            "value": hydrocarbon_ratio,
            "units": "dimensionless",
            "test_status": "source target",
            "interpretation": "reported/inferred Kemel endpoint",
        },
        {
            "quantity": "shared_optical_yield_invariant_ratio_difference",
            "value": hydrocarbon_ratio - biomass_ratio,
            "units": "dimensionless",
            "test_status": "no-refit structural rejection",
            "interpretation": "a common absorbed-photon multiplier predicts zero; observed values differ",
        },
        {
            "quantity": "required_medium_to_large_hydrocarbon_allocation_contrast",
            "value": allocation_contrast,
            "units": "dimensionless",
            "test_status": "no-refit structural implication",
            "interpretation": "equals (P_HC/P_X)_medium divided by (P_HC/P_X)_large",
        },
        {
            "quantity": "conditional_linear_shared_optical_prediction",
            "value": absorption_ratio,
            "units": "dimensionless",
            "test_status": "conditional sensitivity, not an empirical prediction",
            "interpretation": "if the 8-to-79 optical contrast were co-registered to the two productivity states and photon conversion were linear, both output ratios would equal this value",
        },
        {
            "quantity": "conditional_same_state_optical_exponent_if_145_and_434_are_paired",
            "value": conditional_power_law_exponent(
                absorption_ratio, float(large["mean_diameter_um"]), float(medium["mean_diameter_um"])
            ),
            "units": "dimensionless",
            "test_status": "conditional sensitivity, not an empirical fit",
            "interpretation": "requires unreported optical diameters to equal the productivity diameters and no state contrast",
        },
        {
            "quantity": "conditional_shared_saturation_scale_matching_biomass_ratio",
            "value": scale_biomass,
            "units": "m2 kg-1",
            "test_status": "conditional diagnostic, not a fitted reactor parameter",
            "interpretation": "a shared Phi(a)=a/(K+a) would need this K to match biomass alone",
        },
        {
            "quantity": "conditional_shared_saturation_scale_matching_hydrocarbon_ratio",
            "value": scale_hydrocarbon,
            "units": "m2 kg-1",
            "test_status": "conditional diagnostic, not a fitted reactor parameter",
            "interpretation": "the same shared transfer would need a different K to match hydrocarbon alone",
        },
        {
            "quantity": "conditional_log_symmetric_common_ratio",
            "value": common_best_ratio,
            "units": "dimensionless",
            "test_status": "conditional diagnostic",
            "interpretation": "geometric compromise for a shared output ratio; not fitted to a productivity model",
        },
        {
            "quantity": "conditional_shared_saturation_scale_for_log_symmetric_common_ratio",
            "value": scale_common,
            "units": "m2 kg-1",
            "test_status": "conditional diagnostic, not a fitted reactor parameter",
            "interpretation": "shared-transfer scale at the geometric compromise",
        },
    ]

    # This equality is used as a direct numerical guard on the ratio
    # factorisation, independent of any a_abs curve.
    if not math.isclose(
        allocation_contrast,
        (float(medium["hydrocarbon_productivity_g_L_d"]) / float(medium["biomass_productivity_g_L_d"]))
        / (float(large["hydrocarbon_productivity_g_L_d"]) / float(large["biomass_productivity_g_L_d"])),
        rel_tol=0.0,
        abs_tol=1.0e-12,
    ):
        raise RuntimeError("hydrocarbon allocation factorisation is inconsistent")
    if not (scale_hydrocarbon > scale_biomass > 0.0):
        raise RuntimeError("conditional saturation diagnostic should increase with target ratio")
    return rows


def make_measurement_plan_rows() -> list[dict[str, object]]:
    """State the minimum co-registered test needed to complete identification."""

    rows: list[dict[str, object]] = []
    for state_stratum in RECOMMENDED_STATE_STRATA:
        for diameter in RECOMMENDED_DIAMETER_BINS_UM:
            rows.append(
                {
                    "measurement_block": "optical calibration",
                    "diameter_bin_center_um": diameter,
                    "measured_state_stratum": state_stratum,
                    "required_co_registered_readouts": "diameter distribution, reduced-state marker, dry mass, a_abs or MRPA/scattering measurement",
                    "fit_or_holdout_role": "fit a_abs(D,a) only; keep colony-dynamics parameters frozen",
                    "status": "prospective minimum factorial design, not source data",
                }
            )
    rows.extend(
        [
            {
                "measurement_block": "reactor prediction",
                "diameter_bin_center_um": 145.0,
                "measured_state_stratum": "culture-specific measured state",
                "required_co_registered_readouts": "incident/internal light, biomass productivity, hydrocarbon productivity, PBE distribution",
                "fit_or_holdout_role": "reserved medium-culture productivity target after optical calibration",
                "status": "prospective no-refit validation",
            },
            {
                "measurement_block": "reactor prediction",
                "diameter_bin_center_um": 434.0,
                "measured_state_stratum": "culture-specific measured state",
                "required_co_registered_readouts": "incident/internal light, biomass productivity, hydrocarbon productivity, PBE distribution",
                "fit_or_holdout_role": "reserved large-culture productivity target after optical calibration",
                "status": "prospective no-refit validation",
            },
        ]
    )
    return rows


def metric_value(rows: list[dict[str, object]], quantity: str) -> float:
    """Find a scalar result by its stable quantity identifier."""

    return float(next(row["value"] for row in rows if row["quantity"] == quantity))


def verify_analysis(
    source_rows: list[dict[str, object]],
    curve_rows: list[dict[str, object]],
    exponent_rows: list[dict[str, object]],
    ratio_rows: list[dict[str, object]],
    plan_rows: list[dict[str, object]],
) -> None:
    """Guard the source boundary and the no-refit invariant."""

    source_readme = SOURCE_README.read_text().lower()
    if "mean diameters corresponding to the 8.0 and 79" not in source_readme:
        raise RuntimeError("the source package no longer documents the missing optical diameter boundary")
    if sum("not reported" in str(row["source_pairing_status"]) for row in source_rows[:2]) != 2:
        raise RuntimeError("optical endpoints were accidentally treated as paired size/state observations")
    if not all(row["curve_type"] == "conditional_same_state_power_law_not_a_fit" for row in curve_rows + exponent_rows):
        raise RuntimeError("a conditional optical sensitivity lost its explicit boundary label")
    biomass_ratio = metric_value(ratio_rows, "observed_medium_to_large_biomass_productivity_ratio")
    hydrocarbon_ratio = metric_value(ratio_rows, "observed_medium_to_large_hydrocarbon_productivity_ratio")
    allocation = metric_value(ratio_rows, "required_medium_to_large_hydrocarbon_allocation_contrast")
    if not (hydrocarbon_ratio > biomass_ratio > 1.0 and allocation > 1.0):
        raise RuntimeError("the source productivity ratios must reject a shared-yield optical invariant")
    if len(plan_rows) != len(RECOMMENDED_DIAMETER_BINS_UM) * len(RECOMMENDED_STATE_STRATA) + 2:
        raise RuntimeError("the prospective co-registered optical design is incomplete")


def write_summary(
    source_rows: list[dict[str, object]],
    ratio_rows: list[dict[str, object]],
    plan_rows: list[dict[str, object]],
) -> None:
    """Write a human-readable interpretation boundary beside the CSV outputs."""

    optical_contrast = next(
        float(row["value"]) for row in source_rows if row["metric"] == "disrupted_to_intact_optical_contrast"
    )
    biomass_ratio = metric_value(ratio_rows, "observed_medium_to_large_biomass_productivity_ratio")
    hydrocarbon_ratio = metric_value(ratio_rows, "observed_medium_to_large_hydrocarbon_productivity_ratio")
    allocation = metric_value(ratio_rows, "required_medium_to_large_hydrocarbon_allocation_contrast")
    conditional_gamma = metric_value(
        ratio_rows, "conditional_same_state_optical_exponent_if_145_and_434_are_paired"
    )
    conditional_biomass_scale = metric_value(
        ratio_rows, "conditional_shared_saturation_scale_matching_biomass_ratio"
    )
    conditional_hydrocarbon_scale = metric_value(
        ratio_rows, "conditional_shared_saturation_scale_matching_hydrocarbon_ratio"
    )

    text = f"""# Optics-to-productivity source-bound test (Point 7)

## Source-supported inputs

- The accessible Kemel package reports an intact-to-disrupted mass-absorption
  contrast of {optical_contrast:.3f} (8.0 to 79 m2 kg-1).
- It separately reports medium/large biomass and hydrocarbon productivity
  ratios of {biomass_ratio:.3f} and {hydrocarbon_ratio:.3f}. The large-culture
  values are derived from rounded reported percentage gains.
- Neither optical endpoint is paired with a diameter or reduced state, and
  neither productivity endpoint is paired with a mass-absorption measurement.
  The source therefore does not identify an empirical continuous a_abs(D, a)
  relation.

## Fixed-dynamics, no-refit structural test

With the colony distribution/dynamics fixed, an optics-only common-yield
submodel has P_X = Y_X Phi and P_HC = Y_HC Phi. It necessarily predicts
R_X = R_HC irrespective of the unidentifiable a_abs curve or the detailed
photon-transfer function. The two reported ratios differ. Equivalently, the
hydrocarbon-per-biomass productivity fraction is {allocation:.3f} times larger
for the medium than for the large culture. A state- or product-specific
hydrocarbon-allocation response is consequently required beyond a common
optical multiplier. This is a structural rejection, not a fit of an
allocation coefficient.

## Conditional sensitivity only

If, contrary to the available record, 8.0 and 79 m2 kg-1 were assigned to the
434 and 145 um productivity states at the same reduced state, the implied
power-law exponent would be gamma = {conditional_gamma:.3f}. A shared
Phi(a)=a/(K+a) transfer would need K={conditional_biomass_scale:.3f} m2 kg-1
to match biomass alone and K={conditional_hydrocarbon_scale:.3f} m2 kg-1 to
match hydrocarbon alone. These incompatible conditional values illustrate the
same fixed-yield failure; they are not fitted physical parameters and must not
be read as source-derived optical measurements.

## Measurement consequence

The proposed optical calibration grid has
{len(RECOMMENDED_DIAMETER_BINS_UM)} diameter bins by
{len(RECOMMENDED_STATE_STRATA)} measured state strata
({len(plan_rows) - 2} co-registered optical cells), followed by the 145 and
434 um productivity cultures as reserved no-refit reactor targets. Each
optical row must contain diameter distribution, state marker, dry mass, and
an absorption/MRPA or scattering measurement. The biological PBE parameters
remain frozen while only the optical relation is calibrated.

## Interpretation boundary

The present calculation does not establish a strain-specific a_abs(D, a)
curve, a reactor saturation constant, or a hydrocarbon-allocation law. It
does establish which missing pairing prevents that calibration and which
shared-yield optical null is incompatible with the two reported productivity
ratios.
"""
    (RESULTS / "optics_to_productivity_summary.md").write_text(text)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    absorption_large, absorption_small, medium, large = source_inputs()
    source_rows = make_source_audit_rows(absorption_large, absorption_small, medium, large)
    curve_rows, exponent_rows = make_conditional_optical_rows(
        absorption_large, absorption_small, float(large["mean_diameter_um"])
    )
    ratio_rows = make_ratio_rows(absorption_large, absorption_small, medium, large)
    plan_rows = make_measurement_plan_rows()
    verify_analysis(source_rows, curve_rows, exponent_rows, ratio_rows, plan_rows)

    write_csv(RESULTS / "optics_to_productivity_source_audit.csv", source_rows)
    write_csv(RESULTS / "optics_to_productivity_conditional_curves.csv", curve_rows)
    write_csv(RESULTS / "optics_to_productivity_effective_exponent.csv", exponent_rows)
    write_csv(RESULTS / "optics_to_productivity_ratio_test.csv", ratio_rows)
    write_csv(RESULTS / "optics_to_productivity_measurement_plan.csv", plan_rows)
    write_summary(source_rows, ratio_rows, plan_rows)

    print(f"Wrote {RESULTS / 'optics_to_productivity_summary.md'}")
    print(
        "Optical contrast={:.3f}; biomass ratio={:.3f}; hydrocarbon ratio={:.3f}; "
        "required allocation contrast={:.3f}".format(
            absorption_small / absorption_large,
            metric_value(ratio_rows, "observed_medium_to_large_biomass_productivity_ratio"),
            metric_value(ratio_rows, "observed_medium_to_large_hydrocarbon_productivity_ratio"),
            metric_value(ratio_rows, "required_medium_to_large_hydrocarbon_allocation_contrast"),
        )
    )


if __name__ == "__main__":
    main()
