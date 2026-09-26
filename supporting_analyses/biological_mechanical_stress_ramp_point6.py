#!/usr/bin/env python3
"""Prospective biological-versus-mechanical stress-ramp identification test.

This is an in-silico protocol benchmark, not a fit to a published
``Botryococcus braunii`` stress experiment.  The project literature supports
separate biological and hydrodynamic event channels, but it does not contain
the controlled stress-ramp, paired mother--daughter dataset required to
estimate their rates and kernels.  The benchmark therefore creates a
reproducible event-level dataset from declared contrasts and asks whether the
proposed measurement design recovers them without refitting on held-out ramps.

Every stress program receives the *same measured maturation history*.  Stress
is reported as tau_h / sigma_ref[a(t)], where the state-conditioned reference
cohesive stress is assumed to have been measured independently along that
common history.  Below-threshold
reference ramps identify the maturation-gated biological hazard; the resulting
biological fit is frozen while above-threshold training ramps identify the
hydrodynamic hazard.  Event branch topology (a retained mother plus new branch
versus two free fragments) and a shell label provide the channel readout used
to fit the two daughter kernels separately.  A held-out up-ramp and a
held-out down-ramp test the resulting rates and kernel-derived local
log-normal response.

The declared parameter values, parent cohort, and local OU coefficients are
measurement-design contrasts only.  They must not be interpreted as measured
rates, cohesive stresses, daughter kernels, or log-normal parameters for a
real strain.  A real implementation would replace the synthetic rows with
tracked events and retain the same split-before-fit logic.

Run from the workspace root:

    python3 CODES/biological_mechanical_stress_ramp_point6.py
    python3 CODES/plot_biological_mechanical_stress_ramp.py
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, expit
from scipy.stats import beta as beta_distribution


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


# ---------------------------------------------------------------------------
# Declared protocol and generator contrasts.  None is inferred from a source
# experiment.  The stress coordinate is normalized by an independently
# measured, state-conditioned cohesive reference so that the onset
# r = tau_h/sigma_ref[a(t)] = 1 is not another freely tuned stress threshold
# in this prospective exercise.
# ---------------------------------------------------------------------------
RANDOM_SEED = 20260908
N_LINEAGES_PER_RAMP = 3000
DURATION_D = 8.0
DT_D = 0.05
PARENT_MEDIAN_DIAMETER_UM = 340.0
PARENT_LOG_STANDARD_DEVIATION = 0.17

TRUE_BIOLOGICAL_MAX_RATE_D_1 = 0.185
TRUE_BIOLOGICAL_MATURATION_CENTER = 0.60
TRUE_BIOLOGICAL_MATURATION_WIDTH = 0.070
TRUE_HYDRODYNAMIC_PREFactor_D_1 = 0.480
TRUE_HYDRODYNAMIC_EXPONENT = 1.55

# Hydrodynamic events have two free, material-conserving fragments.  The
# smaller-fragment material fraction is 0.5 U, U ~ Beta(alpha, beta).  A
# biological event has a retained mother branch and a new branch, with a
# shell-labelled material fraction removed; the new-branch fraction is
# (1-shell) U.  These are deliberately distinct, observable event maps.
TRUE_HYDRODYNAMIC_KERNEL_ALPHA = 5.0
TRUE_HYDRODYNAMIC_KERNEL_BETA = 3.0
TRUE_BIOLOGICAL_KERNEL_ALPHA = 3.5
TRUE_BIOLOGICAL_KERNEL_BETA = 9.0
TRUE_BIOLOGICAL_SHELL_LOSS_FRACTION = 0.08

# Fixed local OU diagnostic closure.  These coefficients are not fitted to
# synthetic histograms: the event-rate and kernel fits alone drive the stress
# response displayed in the companion figure.
OU_REFERENCE_LOG_MEDIAN = math.log(PARENT_MEDIAN_DIAMETER_UM)
OU_BASE_RESTORING_RATE_D_1 = 0.270
OU_CONTINUOUS_DIFFUSION_D_1 = 0.018
OU_HYDRODYNAMIC_LINEAR_RESTORING_COEFFICIENT = 0.050
OU_HYDRODYNAMIC_QUADRATIC_RESTORING_COEFFICIENT_D = 0.85
OU_REFERENCE_MATURATION_STATE = 0.86


@dataclass(frozen=True)
class Ramp:
    """One controlled stress program applied to a common maturation history."""

    identifier: str
    split: str
    stress_start: float
    stress_end: float


RAMPS = (
    Ramp("subthreshold_train", "train", 0.65, 0.90),
    Ramp("mild_up_train", "train", 0.65, 1.35),
    Ramp("strong_up_train", "train", 0.65, 1.80),
    Ramp("mid_up_test", "test", 0.65, 1.55),
    Ramp("down_test", "test", 1.80, 0.65),
)


@dataclass(frozen=True)
class BiologicalRateFit:
    maximum_rate_d_1: float
    maturation_center: float
    maturation_width: float
    objective: float


@dataclass(frozen=True)
class HydrodynamicRateFit:
    prefactor_d_1: float
    exponent: float
    objective: float


@dataclass(frozen=True)
class BetaKernelFit:
    alpha: float
    beta: float
    objective: float


def maturation_history(time_d: float | np.ndarray) -> float | np.ndarray:
    """Shared, recorded maturation trajectory for every stress program."""

    return 0.10 + 0.82 * expit((np.asarray(time_d) - 2.60) / 0.45)


def stress_ratio(ramp: Ramp, time_d: float | np.ndarray) -> float | np.ndarray:
    """Linearly ramp normalized stress between days 1 and 7."""

    fraction = np.clip((np.asarray(time_d) - 1.0) / 6.0, 0.0, 1.0)
    return ramp.stress_start + (ramp.stress_end - ramp.stress_start) * fraction


def biological_rate(
    maturation_state: float | np.ndarray,
    maximum_rate_d_1: float,
    maturation_center: float,
    maturation_width: float,
) -> float | np.ndarray:
    """Maturation-gated biological daughter-production hazard."""

    return maximum_rate_d_1 * expit((np.asarray(maturation_state) - maturation_center) / maturation_width)


def hydrodynamic_rate(
    normalized_stress: float | np.ndarray,
    prefactor_d_1: float,
    exponent: float,
) -> float | np.ndarray:
    """Excess-stress hydrodynamic breakup hazard at fixed physiology."""

    excess = np.maximum(np.asarray(normalized_stress) - 1.0, 0.0)
    return prefactor_d_1 * np.power(excess, exponent)


def probability_components(hyd_rate: float, bio_rate: float, dt_d: float) -> tuple[float, float, float]:
    """Exact competing-risk event probabilities over one discrete interval."""

    total_rate = hyd_rate + bio_rate
    no_event = math.exp(-total_rate * dt_d)
    event_probability = 1.0 - no_event
    if total_rate == 0.0:
        return 0.0, 0.0, 1.0
    return (
        event_probability * hyd_rate / total_rate,
        event_probability * bio_rate / total_rate,
        no_event,
    )


def generate_synthetic_protocol() -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Generate at-risk intervals, event records, and branch observations.

    Each lineage leaves the risk set after its first recorded separation event.
    This makes the at-risk count explicit and permits a cause-specific survival
    fit rather than treating a post-event descendant as a new independent
    mother.  The event branch topology is the prospective structural readout:
    ``retained_mother`` classifies biological events without using a fitted
    rate label.
    """

    rng = np.random.default_rng(RANDOM_SEED)
    exposure_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    branch_rows: list[dict[str, object]] = []
    event_identifier = 0

    times = np.arange(0.0, DURATION_D, DT_D)
    for ramp in RAMPS:
        at_risk = N_LINEAGES_PER_RAMP
        for interval_index, time_left_d in enumerate(times):
            if at_risk == 0:
                break
            time_mid_d = float(time_left_d + 0.5 * DT_D)
            state = float(maturation_history(time_mid_d))
            stress = float(stress_ratio(ramp, time_mid_d))
            beta_bio = float(
                biological_rate(
                    state,
                    TRUE_BIOLOGICAL_MAX_RATE_D_1,
                    TRUE_BIOLOGICAL_MATURATION_CENTER,
                    TRUE_BIOLOGICAL_MATURATION_WIDTH,
                )
            )
            beta_hyd = float(
                hydrodynamic_rate(stress, TRUE_HYDRODYNAMIC_PREFactor_D_1, TRUE_HYDRODYNAMIC_EXPONENT)
            )
            p_hyd, p_bio, _ = probability_components(beta_hyd, beta_bio, DT_D)
            number_hyd = int(rng.binomial(at_risk, p_hyd))
            remaining_after_hyd = at_risk - number_hyd
            conditional_bio_probability = p_bio / (1.0 - p_hyd) if p_hyd < 1.0 else 0.0
            number_bio = int(rng.binomial(remaining_after_hyd, conditional_bio_probability))

            exposure_rows.append(
                {
                    "ramp_id": ramp.identifier,
                    "split": ramp.split,
                    "interval_index": interval_index,
                    "time_d": time_mid_d,
                    "interval_d": DT_D,
                    "normalized_stress_tau_over_sigma_ref": stress,
                    "maturation_state": state,
                    "at_risk_lineages": at_risk,
                    "hydrodynamic_events": number_hyd,
                    "biological_events": number_bio,
                    "true_hydrodynamic_rate_d_1": beta_hyd,
                    "true_biological_rate_d_1": beta_bio,
                }
            )

            for _ in range(number_hyd):
                event_identifier += 1
                parent_diameter = float(
                    math.exp(rng.normal(math.log(PARENT_MEDIAN_DIAMETER_UM), PARENT_LOG_STANDARD_DEVIATION))
                )
                scaled_small_fraction = float(rng.beta(TRUE_HYDRODYNAMIC_KERNEL_ALPHA, TRUE_HYDRODYNAMIC_KERNEL_BETA))
                small_fraction = 0.5 * scaled_small_fraction
                large_fraction = 1.0 - small_fraction
                event_rows.append(
                    {
                        "event_id": event_identifier,
                        "ramp_id": ramp.identifier,
                        "split": ramp.split,
                        "time_d": time_mid_d,
                        "normalized_stress_tau_over_sigma_ref": stress,
                        "maturation_state": state,
                        "observed_channel": "hydrodynamic",
                        "retained_mother_observed": False,
                        "shell_loss_fraction_observed": 0.0,
                        "viable_descendant_count": 2,
                        "parent_diameter_um": parent_diameter,
                    }
                )
                for branch_role, fraction in (("small_fragment", small_fraction), ("large_fragment", large_fraction)):
                    branch_rows.append(
                        {
                            "event_id": event_identifier,
                            "ramp_id": ramp.identifier,
                            "split": ramp.split,
                            "time_d": time_mid_d,
                            "observed_channel": "hydrodynamic",
                            "branch_role": branch_role,
                            "parent_diameter_um": parent_diameter,
                            "material_fraction_of_parent": fraction,
                            "descendant_diameter_um": parent_diameter * fraction ** (1.0 / 3.0),
                            "descendant_parent_diameter_ratio": fraction ** (1.0 / 3.0),
                        }
                    )

            for _ in range(number_bio):
                event_identifier += 1
                parent_diameter = float(
                    math.exp(rng.normal(math.log(PARENT_MEDIAN_DIAMETER_UM), PARENT_LOG_STANDARD_DEVIATION))
                )
                available_fraction = 1.0 - TRUE_BIOLOGICAL_SHELL_LOSS_FRACTION
                scaled_new_fraction = float(rng.beta(TRUE_BIOLOGICAL_KERNEL_ALPHA, TRUE_BIOLOGICAL_KERNEL_BETA))
                new_fraction = available_fraction * scaled_new_fraction
                retained_fraction = available_fraction - new_fraction
                event_rows.append(
                    {
                        "event_id": event_identifier,
                        "ramp_id": ramp.identifier,
                        "split": ramp.split,
                        "time_d": time_mid_d,
                        "normalized_stress_tau_over_sigma_ref": stress,
                        "maturation_state": state,
                        "observed_channel": "biological",
                        "retained_mother_observed": True,
                        "shell_loss_fraction_observed": TRUE_BIOLOGICAL_SHELL_LOSS_FRACTION,
                        "viable_descendant_count": 2,
                        "parent_diameter_um": parent_diameter,
                    }
                )
                for branch_role, fraction in (("retained_mother", retained_fraction), ("new_daughter", new_fraction)):
                    branch_rows.append(
                        {
                            "event_id": event_identifier,
                            "ramp_id": ramp.identifier,
                            "split": ramp.split,
                            "time_d": time_mid_d,
                            "observed_channel": "biological",
                            "branch_role": branch_role,
                            "parent_diameter_um": parent_diameter,
                            "material_fraction_of_parent": fraction,
                            "descendant_diameter_um": parent_diameter * fraction ** (1.0 / 3.0),
                            "descendant_parent_diameter_ratio": fraction ** (1.0 / 3.0),
                        }
                    )
            at_risk -= number_hyd + number_bio

    return exposure_rows, event_rows, branch_rows


def fit_biological_rate(exposure_rows: list[dict[str, object]]) -> BiologicalRateFit:
    """Fit beta_bio only on the wholly subthreshold reference ramp."""

    reference = [row for row in exposure_rows if row["ramp_id"] == "subthreshold_train"]
    if any(int(row["hydrodynamic_events"]) != 0 for row in reference):
        raise RuntimeError("subthreshold reference unexpectedly contains hydrodynamic events")
    state = np.asarray([float(row["maturation_state"]) for row in reference])
    at_risk = np.asarray([float(row["at_risk_lineages"]) for row in reference])
    events = np.asarray([float(row["biological_events"]) for row in reference])

    def objective(unconstrained: np.ndarray) -> float:
        maximum_rate = math.exp(float(unconstrained[0]))
        center = float(unconstrained[1])
        width = math.exp(float(unconstrained[2]))
        rate = np.asarray(biological_rate(state, maximum_rate, center, width), dtype=float)
        probability = np.clip(-np.expm1(-rate * DT_D), 1.0e-15, 1.0 - 1.0e-15)
        return float(-np.sum(events * np.log(probability) + (at_risk - events) * np.log1p(-probability)))

    result = minimize(
        objective,
        x0=np.asarray((math.log(0.16), 0.56, math.log(0.09))),
        method="L-BFGS-B",
        bounds=((math.log(0.02), math.log(1.0)), (0.35, 0.85), (math.log(0.02), math.log(0.25))),
    )
    if not result.success:
        raise RuntimeError(f"biological-rate fit failed: {result.message}")
    return BiologicalRateFit(
        maximum_rate_d_1=math.exp(float(result.x[0])),
        maturation_center=float(result.x[1]),
        maturation_width=math.exp(float(result.x[2])),
        objective=float(result.fun),
    )


def fit_hydrodynamic_rate(exposure_rows: list[dict[str, object]], biological_fit: BiologicalRateFit) -> HydrodynamicRateFit:
    """Fit beta_hyd on training ramps with the biological hazard frozen.

    The likelihood is the exact three-outcome competing-risk likelihood for
    hydrodynamic event, biological event, or no event.  It therefore respects
    depletion of the at-risk cohort by the already fitted biological channel.
    """

    training = [row for row in exposure_rows if row["split"] == "train"]
    stress = np.asarray([float(row["normalized_stress_tau_over_sigma_ref"]) for row in training])
    maturation = np.asarray([float(row["maturation_state"]) for row in training])
    at_risk = np.asarray([float(row["at_risk_lineages"]) for row in training])
    observed_hyd = np.asarray([float(row["hydrodynamic_events"]) for row in training])
    observed_bio = np.asarray([float(row["biological_events"]) for row in training])
    observed_none = at_risk - observed_hyd - observed_bio
    frozen_bio_rate = np.asarray(
        biological_rate(
            maturation,
            biological_fit.maximum_rate_d_1,
            biological_fit.maturation_center,
            biological_fit.maturation_width,
        ),
        dtype=float,
    )

    def objective(unconstrained: np.ndarray) -> float:
        prefactor = math.exp(float(unconstrained[0]))
        exponent = math.exp(float(unconstrained[1]))
        hyd_rate = np.asarray(hydrodynamic_rate(stress, prefactor, exponent), dtype=float)
        total_rate = hyd_rate + frozen_bio_rate
        p_none = np.exp(-total_rate * DT_D)
        p_event = 1.0 - p_none
        p_hyd = np.where(total_rate > 0.0, p_event * hyd_rate / total_rate, 0.0)
        p_bio = np.where(total_rate > 0.0, p_event * frozen_bio_rate / total_rate, 0.0)
        epsilon = 1.0e-300
        return float(
            -np.sum(
                observed_hyd * np.log(np.maximum(p_hyd, epsilon))
                + observed_bio * np.log(np.maximum(p_bio, epsilon))
                + observed_none * np.log(np.maximum(p_none, epsilon))
            )
        )

    result = minimize(
        objective,
        x0=np.asarray((math.log(0.40), math.log(1.30))),
        method="L-BFGS-B",
        bounds=((math.log(0.02), math.log(2.0)), (math.log(0.25), math.log(4.0))),
    )
    if not result.success:
        raise RuntimeError(f"hydrodynamic-rate fit failed: {result.message}")
    return HydrodynamicRateFit(
        prefactor_d_1=math.exp(float(result.x[0])),
        exponent=math.exp(float(result.x[1])),
        objective=float(result.fun),
    )


def fit_beta_kernel(values: np.ndarray, label: str) -> BetaKernelFit:
    """Maximum-likelihood fit for a dimensionless observed kernel coordinate."""

    clipped = np.clip(np.asarray(values, dtype=float), 1.0e-8, 1.0 - 1.0e-8)
    if clipped.size < 20:
        raise RuntimeError(f"too few {label} kernel observations for a stable fit")

    def objective(unconstrained: np.ndarray) -> float:
        alpha = math.exp(float(unconstrained[0]))
        beta_value = math.exp(float(unconstrained[1]))
        log_density = (
            (alpha - 1.0) * np.log(clipped)
            + (beta_value - 1.0) * np.log1p(-clipped)
            - betaln(alpha, beta_value)
        )
        return float(-np.sum(log_density))

    mean = float(np.mean(clipped))
    variance = float(np.var(clipped, ddof=1))
    concentration = max(mean * (1.0 - mean) / max(variance, 1.0e-5) - 1.0, 2.0)
    initial_alpha = max(mean * concentration, 0.25)
    initial_beta = max((1.0 - mean) * concentration, 0.25)
    result = minimize(
        objective,
        x0=np.log((initial_alpha, initial_beta)),
        method="L-BFGS-B",
        bounds=((math.log(0.15), math.log(60.0)), (math.log(0.15), math.log(60.0))),
    )
    if not result.success:
        raise RuntimeError(f"{label} kernel fit failed: {result.message}")
    return BetaKernelFit(alpha=math.exp(float(result.x[0])), beta=math.exp(float(result.x[1])), objective=float(result.fun))


def fit_kernels(branch_rows: list[dict[str, object]], event_rows: list[dict[str, object]]) -> tuple[BetaKernelFit, BetaKernelFit, float]:
    """Fit b_hyd and b_bio from their observed branch maps separately."""

    training_event_ids = {int(row["event_id"]) for row in event_rows if row["split"] == "train"}
    hyd_scaled_small = np.asarray(
        [
            2.0 * float(row["material_fraction_of_parent"])
            for row in branch_rows
            if int(row["event_id"]) in training_event_ids and row["branch_role"] == "small_fragment"
        ]
    )
    biological_new = np.asarray(
        [
            float(row["material_fraction_of_parent"])
            for row in branch_rows
            if int(row["event_id"]) in training_event_ids and row["branch_role"] == "new_daughter"
        ]
    )
    shell_values = np.asarray(
        [
            float(row["shell_loss_fraction_observed"])
            for row in event_rows
            if row["split"] == "train" and row["observed_channel"] == "biological"
        ]
    )
    if shell_values.size == 0:
        raise RuntimeError("no biological shell labels were available for the kernel fit")
    observed_shell_loss = float(np.mean(shell_values))
    biological_scaled_new = biological_new / (1.0 - observed_shell_loss)
    return (
        fit_beta_kernel(hyd_scaled_small, "hydrodynamic"),
        fit_beta_kernel(biological_scaled_new, "biological"),
        observed_shell_loss,
    )


def beta_kernel_density(
    material_fraction: np.ndarray,
    fit: BetaKernelFit,
    channel: str,
    shell_loss_fraction: float,
) -> np.ndarray:
    """Expected daughter-number density b(z) for a two-branch event map."""

    z = np.asarray(material_fraction, dtype=float)
    density = np.zeros_like(z)
    if channel == "hydrodynamic":
        small_mask = (z > 0.0) & (z < 0.5)
        large_mask = (z > 0.5) & (z < 1.0)
        density[small_mask] = beta_distribution.pdf(2.0 * z[small_mask], fit.alpha, fit.beta)
        density[large_mask] = beta_distribution.pdf(2.0 * (1.0 - z[large_mask]), fit.alpha, fit.beta)
        return density
    if channel == "biological":
        viable = 1.0 - shell_loss_fraction
        mask = (z > 0.0) & (z < viable)
        scaled_new = z[mask] / viable
        scaled_retained = (viable - z[mask]) / viable
        density[mask] = 0.5 * (
            beta_distribution.pdf(scaled_new, fit.alpha, fit.beta)
            + beta_distribution.pdf(scaled_retained, fit.alpha, fit.beta)
        ) / viable
        return density
    raise ValueError(f"unknown channel {channel!r}")


def kernel_log_jump_moments(fit: BetaKernelFit, channel: str, shell_loss_fraction: float) -> tuple[float, float]:
    """Return E[delta] and E[delta^2] over viable daughter branches."""

    u = np.linspace(1.0e-5, 1.0 - 1.0e-5, 20001)
    density = beta_distribution.pdf(u, fit.alpha, fit.beta)
    density /= np.trapezoid(density, u)
    if channel == "hydrodynamic":
        first_fraction = 0.5 * u
        second_fraction = 1.0 - first_fraction
    elif channel == "biological":
        viable = 1.0 - shell_loss_fraction
        first_fraction = viable * u
        second_fraction = viable - first_fraction
    else:
        raise ValueError(f"unknown channel {channel!r}")
    first_jump = np.log(first_fraction) / 3.0
    second_jump = np.log(second_fraction) / 3.0
    mean = float(np.trapezoid(0.5 * (first_jump + second_jump) * density, u))
    second_moment = float(np.trapezoid(0.5 * (first_jump**2 + second_jump**2) * density, u))
    return mean, second_moment


def local_lognormal_response(
    hyd_fit: HydrodynamicRateFit,
    bio_fit: BiologicalRateFit,
    hyd_kernel: BetaKernelFit,
    bio_kernel: BetaKernelFit,
    shell_loss: float,
    stress_values: np.ndarray,
) -> list[dict[str, object]]:
    """Calculate Appendix-A local OU center and width from fitted event maps."""

    hyd_mean_jump, hyd_second_jump = kernel_log_jump_moments(hyd_kernel, "hydrodynamic", shell_loss)
    bio_mean_jump, bio_second_jump = kernel_log_jump_moments(bio_kernel, "biological", shell_loss)
    hyd_true = BetaKernelFit(TRUE_HYDRODYNAMIC_KERNEL_ALPHA, TRUE_HYDRODYNAMIC_KERNEL_BETA, math.nan)
    bio_true = BetaKernelFit(TRUE_BIOLOGICAL_KERNEL_ALPHA, TRUE_BIOLOGICAL_KERNEL_BETA, math.nan)
    hyd_mean_true, hyd_second_true = kernel_log_jump_moments(hyd_true, "hydrodynamic", 0.0)
    bio_mean_true, bio_second_true = kernel_log_jump_moments(
        bio_true, "biological", TRUE_BIOLOGICAL_SHELL_LOSS_FRACTION
    )

    fitted_bio_rate = float(
        biological_rate(
            OU_REFERENCE_MATURATION_STATE,
            bio_fit.maximum_rate_d_1,
            bio_fit.maturation_center,
            bio_fit.maturation_width,
        )
    )
    true_bio_rate = float(
        biological_rate(
            OU_REFERENCE_MATURATION_STATE,
            TRUE_BIOLOGICAL_MAX_RATE_D_1,
            TRUE_BIOLOGICAL_MATURATION_CENTER,
            TRUE_BIOLOGICAL_MATURATION_WIDTH,
        )
    )

    def solve_one(hyd_rate_value: float, bio_rate_value: float, mean_hyd: float, second_hyd: float, mean_bio: float, second_bio: float) -> tuple[float, float, float]:
        kappa = (
            OU_BASE_RESTORING_RATE_D_1
            + OU_HYDRODYNAMIC_LINEAR_RESTORING_COEFFICIENT * hyd_rate_value
            + OU_HYDRODYNAMIC_QUADRATIC_RESTORING_COEFFICIENT_D * hyd_rate_value**2
        )
        drift_offset = bio_rate_value * mean_bio + hyd_rate_value * mean_hyd
        center = OU_REFERENCE_LOG_MEDIAN + drift_offset / kappa
        diffusion_coefficient_twice = (
            2.0 * OU_CONTINUOUS_DIFFUSION_D_1
            + bio_rate_value * second_bio
            + hyd_rate_value * second_hyd
        )
        log_standard_deviation = math.sqrt(diffusion_coefficient_twice / (2.0 * kappa))
        return center, log_standard_deviation, kappa

    rows: list[dict[str, object]] = []
    for stress in stress_values:
        fitted_hyd_rate = float(hydrodynamic_rate(stress, hyd_fit.prefactor_d_1, hyd_fit.exponent))
        true_hyd_rate = float(
            hydrodynamic_rate(stress, TRUE_HYDRODYNAMIC_PREFactor_D_1, TRUE_HYDRODYNAMIC_EXPONENT)
        )
        fitted_center, fitted_width, fitted_kappa = solve_one(
            fitted_hyd_rate,
            fitted_bio_rate,
            hyd_mean_jump,
            hyd_second_jump,
            bio_mean_jump,
            bio_second_jump,
        )
        true_center, true_width, true_kappa = solve_one(
            true_hyd_rate,
            true_bio_rate,
            hyd_mean_true,
            hyd_second_true,
            bio_mean_true,
            bio_second_true,
        )
        rows.append(
            {
                "normalized_stress_tau_over_sigma_ref": float(stress),
                "reference_maturation_state": OU_REFERENCE_MATURATION_STATE,
                "fitted_hydrodynamic_rate_d_1": fitted_hyd_rate,
                "fitted_biological_rate_d_1": fitted_bio_rate,
                "fitted_log_center_m": fitted_center,
                "fitted_median_diameter_um": math.exp(fitted_center),
                "fitted_log_standard_deviation_s": fitted_width,
                "fitted_geometric_standard_deviation": math.exp(fitted_width),
                "fitted_restoring_rate_kappa_d_1": fitted_kappa,
                "generator_log_center_m": true_center,
                "generator_median_diameter_um": math.exp(true_center),
                "generator_log_standard_deviation_s": true_width,
                "generator_geometric_standard_deviation": math.exp(true_width),
                "generator_restoring_rate_kappa_d_1": true_kappa,
            }
        )
    return rows


def make_prediction_rows(
    exposure_rows: list[dict[str, object]], biological_fit: BiologicalRateFit, hydrodynamic_fit: HydrodynamicRateFit
) -> list[dict[str, object]]:
    """Attach no-refit fitted hazards to every training and held-out interval."""

    records: list[dict[str, object]] = []
    for row in exposure_rows:
        state = float(row["maturation_state"])
        stress = float(row["normalized_stress_tau_over_sigma_ref"])
        at_risk = float(row["at_risk_lineages"])
        interval = float(row["interval_d"])
        fitted_bio = float(
            biological_rate(
                state,
                biological_fit.maximum_rate_d_1,
                biological_fit.maturation_center,
                biological_fit.maturation_width,
            )
        )
        fitted_hyd = float(hydrodynamic_rate(stress, hydrodynamic_fit.prefactor_d_1, hydrodynamic_fit.exponent))
        records.append(
            {
                **row,
                "observed_hydrodynamic_rate_d_1": float(row["hydrodynamic_events"]) / (at_risk * interval),
                "observed_biological_rate_d_1": float(row["biological_events"]) / (at_risk * interval),
                "fitted_hydrodynamic_rate_d_1": fitted_hyd,
                "fitted_biological_rate_d_1": fitted_bio,
            }
        )
    return records


def metric_rows(
    prediction_rows: list[dict[str, object]],
    hyd_fit: HydrodynamicRateFit,
    bio_fit: BiologicalRateFit,
    hyd_kernel: BetaKernelFit,
    bio_kernel: BetaKernelFit,
    shell_loss: float,
    lognormal_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Compute compact recovery and held-out diagnostics for the manuscript."""

    held_out = [row for row in prediction_rows if row["split"] == "test"]

    def rmse(column_fitted: str, column_true: str) -> float:
        difference = np.asarray([float(row[column_fitted]) - float(row[column_true]) for row in held_out])
        return float(math.sqrt(np.mean(difference**2)))

    grid = np.linspace(1.0e-4, 0.9999, 5000)
    hyd_true_density = beta_kernel_density(
        grid,
        BetaKernelFit(TRUE_HYDRODYNAMIC_KERNEL_ALPHA, TRUE_HYDRODYNAMIC_KERNEL_BETA, math.nan),
        "hydrodynamic",
        0.0,
    )
    hyd_fitted_density = beta_kernel_density(grid, hyd_kernel, "hydrodynamic", 0.0)
    bio_true_density = beta_kernel_density(
        grid,
        BetaKernelFit(TRUE_BIOLOGICAL_KERNEL_ALPHA, TRUE_BIOLOGICAL_KERNEL_BETA, math.nan),
        "biological",
        TRUE_BIOLOGICAL_SHELL_LOSS_FRACTION,
    )
    bio_fitted_density = beta_kernel_density(grid, bio_kernel, "biological", shell_loss)
    hyd_kernel_tv = 0.5 * float(np.trapezoid(np.abs(hyd_true_density - hyd_fitted_density), grid))
    bio_kernel_tv = 0.5 * float(np.trapezoid(np.abs(bio_true_density - bio_fitted_density), grid))

    median_errors = np.asarray(
        [float(row["fitted_median_diameter_um"]) - float(row["generator_median_diameter_um"]) for row in lognormal_rows]
    )
    width_errors = np.asarray(
        [
            float(row["fitted_log_standard_deviation_s"])
            - float(row["generator_log_standard_deviation_s"])
            for row in lognormal_rows
        ]
    )
    parameter_records = (
        ("biological_rate", "beta_bio_max_d_1", TRUE_BIOLOGICAL_MAX_RATE_D_1, bio_fit.maximum_rate_d_1),
        ("biological_rate", "bio_maturation_center", TRUE_BIOLOGICAL_MATURATION_CENTER, bio_fit.maturation_center),
        ("biological_rate", "bio_maturation_width", TRUE_BIOLOGICAL_MATURATION_WIDTH, bio_fit.maturation_width),
        ("hydrodynamic_rate", "beta_hyd_prefactor_d_1", TRUE_HYDRODYNAMIC_PREFactor_D_1, hyd_fit.prefactor_d_1),
        ("hydrodynamic_rate", "hyd_exponent", TRUE_HYDRODYNAMIC_EXPONENT, hyd_fit.exponent),
        ("hydrodynamic_kernel", "scaled_small_fragment_alpha", TRUE_HYDRODYNAMIC_KERNEL_ALPHA, hyd_kernel.alpha),
        ("hydrodynamic_kernel", "scaled_small_fragment_beta", TRUE_HYDRODYNAMIC_KERNEL_BETA, hyd_kernel.beta),
        ("biological_kernel", "scaled_new_branch_alpha", TRUE_BIOLOGICAL_KERNEL_ALPHA, bio_kernel.alpha),
        ("biological_kernel", "scaled_new_branch_beta", TRUE_BIOLOGICAL_KERNEL_BETA, bio_kernel.beta),
        ("biological_kernel", "shell_loss_fraction", TRUE_BIOLOGICAL_SHELL_LOSS_FRACTION, shell_loss),
    )
    rows: list[dict[str, object]] = [
        {
            "metric_group": group,
            "metric": name,
            "value": fitted,
            "reference_value": true,
            "absolute_difference": abs(fitted - true),
            "definition": "Declared generator versus independently fitted benchmark value.",
        }
        for group, name, true, fitted in parameter_records
    ]
    rows.extend(
        (
            {
                "metric_group": "held_out_rate_prediction",
                "metric": "hydrodynamic_rate_rmse_d_1",
                "value": rmse("fitted_hydrodynamic_rate_d_1", "true_hydrodynamic_rate_d_1"),
                "reference_value": 0.0,
                "absolute_difference": math.nan,
                "definition": "No-refit fitted-versus-generator hydrodynamic hazard across both held-out ramps.",
            },
            {
                "metric_group": "held_out_rate_prediction",
                "metric": "biological_rate_rmse_d_1",
                "value": rmse("fitted_biological_rate_d_1", "true_biological_rate_d_1"),
                "reference_value": 0.0,
                "absolute_difference": math.nan,
                "definition": "No-refit fitted-versus-generator biological hazard across both held-out ramps.",
            },
            {
                "metric_group": "kernel_recovery",
                "metric": "hydrodynamic_kernel_total_variation",
                "value": hyd_kernel_tv,
                "reference_value": 0.0,
                "absolute_difference": math.nan,
                "definition": "Total variation between fitted and declared two-fragment b_hyd.",
            },
            {
                "metric_group": "kernel_recovery",
                "metric": "biological_kernel_total_variation",
                "value": bio_kernel_tv,
                "reference_value": 0.0,
                "absolute_difference": math.nan,
                "definition": "Total variation between fitted and declared retained-mother/new-branch b_bio.",
            },
            {
                "metric_group": "lognormal_projection",
                "metric": "maximum_median_diameter_error_um",
                "value": float(np.max(np.abs(median_errors))),
                "reference_value": 0.0,
                "absolute_difference": math.nan,
                "definition": "Maximum recovered-versus-generator local OU median-diameter difference over the stress grid.",
            },
            {
                "metric_group": "lognormal_projection",
                "metric": "maximum_log_standard_deviation_error",
                "value": float(np.max(np.abs(width_errors))),
                "reference_value": 0.0,
                "absolute_difference": math.nan,
                "definition": "Maximum recovered-versus-generator local OU log-width difference over the stress grid.",
            },
        )
    )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write a non-empty list of homogeneous records as a CSV file."""

    if not rows:
        raise RuntimeError(f"refusing to write empty output {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metric_value(metrics: list[dict[str, object]], name: str) -> float:
    return float(next(row["value"] for row in metrics if row["metric"] == name))


def write_summary(
    exposure_rows: list[dict[str, object]],
    event_rows: list[dict[str, object]],
    metrics: list[dict[str, object]],
    biological_fit: BiologicalRateFit,
    hydrodynamic_fit: HydrodynamicRateFit,
    hyd_kernel: BetaKernelFit,
    bio_kernel: BetaKernelFit,
    shell_loss: float,
    lognormal_rows: list[dict[str, object]],
) -> None:
    """Document the design, results, and interpretation boundary in Markdown."""

    training_events = sum(1 for row in event_rows if row["split"] == "train")
    test_events = sum(1 for row in event_rows if row["split"] == "test")
    training_hyd = sum(1 for row in event_rows if row["split"] == "train" and row["observed_channel"] == "hydrodynamic")
    training_bio = sum(1 for row in event_rows if row["split"] == "train" and row["observed_channel"] == "biological")
    low_stress = lognormal_rows[0]
    high_stress = lognormal_rows[-1]
    text = f"""# Biological-versus-mechanical stress-ramp identification test (Point 6)

## Protocol

- This is a seeded prospective synthetic benchmark, not a fit to a published
  stress-ramp experiment.  It contains {len(exposure_rows)} at-risk intervals
  and {len(event_rows)} recorded separation events ({training_events} training,
  {test_events} held out).
- Five controlled ramps share one deterministic recorded maturation history.
  Stress is normalized as `tau_h / sigma_ref[a(t)]`, with the threshold fixed
  at one by an independently measured, state-conditioned reference cohesive
  stress.  Thus physiology is not changed while stress is varied.
- The subthreshold training ramp identifies only the maturation-gated
  biological hazard.  That fit is frozen before fitting the excess-stress
  hydrodynamic hazard with the other training ramps.  The held-out mid-up and
  down ramps receive no parameter refit.
- A retained mother plus new branch and a structural shell label identify a
  biological event; two free material-conserving fragments identify a
  hydrodynamic event.  The paired branch maps fit `b_bio` and `b_hyd`
  separately.  In a real experiment these are required observations, not
  assumptions available from a terminal size histogram.

## Results

- The training event set contains {training_bio} biological and {training_hyd}
  hydrodynamic events.  The fitted biological parameters are
  `beta_bio,max = {biological_fit.maximum_rate_d_1:.4f} d^-1`,
  `a_c = {biological_fit.maturation_center:.4f}`, and
  `Delta a = {biological_fit.maturation_width:.4f}`.  The fitted hydrodynamic
  parameters are `k_h = {hydrodynamic_fit.prefactor_d_1:.4f} d^-1` and
  `q = {hydrodynamic_fit.exponent:.4f}`.
- Across both held-out ramps, no-refit fitted-versus-generator rate RMSE is
  {metric_value(metrics, 'hydrodynamic_rate_rmse_d_1'):.4f} d^-1 for the
  hydrodynamic channel and {metric_value(metrics, 'biological_rate_rmse_d_1'):.4f}
  d^-1 for the biological channel.  The fitted two-branch kernel total
  variations are {metric_value(metrics, 'hydrodynamic_kernel_total_variation'):.4f}
  for `b_hyd` and {metric_value(metrics, 'biological_kernel_total_variation'):.4f}
  for `b_bio`; the directly labelled shell fraction is {shell_loss:.3f}.
- At the fixed reference maturity `a = {OU_REFERENCE_MATURATION_STATE:.2f}`, the
  event-derived local OU projection changes its fitted median from
  {float(low_stress['fitted_median_diameter_um']):.1f} um at
  `tau_h/sigma_ref = {float(low_stress['normalized_stress_tau_over_sigma_ref']):.2f}`
  to {float(high_stress['fitted_median_diameter_um']):.1f} um at
  `tau_h/sigma_ref = {float(high_stress['normalized_stress_tau_over_sigma_ref']):.2f}`.
  The corresponding fitted log-standard-deviation changes from
  {float(low_stress['fitted_log_standard_deviation_s']):.3f} to
  {float(high_stress['fitted_log_standard_deviation_s']):.3f}.  This is a
  declared local-closure prediction, not a measured distribution response.

## Interpretation boundary

The benchmark establishes protocol adequacy under its declared contrasts: a
fixed-history stress design plus event-resolved branch topology can separate
the two hazards and their daughter kernels, then propagate them into the
Appendix-A log-normal diagnostic without fitting a terminal histogram.  It
does not demonstrate that these numerical rates, kernel shapes, shell loss,
or OU coefficients apply to a real strain.  Without a channel-resolving
lineage/shell readout, aggregate event counts or size distributions would
instead require an explicitly validated latent competing-risk model and would
not justify a separate `beta_bio`, `beta_hyd`, `b_bio`, and `b_hyd` fit.
"""
    (RESULTS / "biological_mechanical_stress_ramp_summary.md").write_text(text)


def verify_protocol(
    exposure_rows: list[dict[str, object]], event_rows: list[dict[str, object]], branch_rows: list[dict[str, object]]
) -> None:
    """Fail loudly if a future edit breaks the stated experimental contrast."""

    if any(
        int(row["hydrodynamic_events"]) > 0 and float(row["normalized_stress_tau_over_sigma_ref"]) <= 1.0
        for row in exposure_rows
    ):
        raise RuntimeError("below-threshold hydrodynamic events violate the declared excess-stress generator")
    state_by_time: dict[float, set[float]] = {}
    for row in exposure_rows:
        state_by_time.setdefault(float(row["time_d"]), set()).add(round(float(row["maturation_state"]), 12))
    if any(len(values) != 1 for values in state_by_time.values()):
        raise RuntimeError("stress programs do not share a fixed physiological history")
    event_by_id = {int(row["event_id"]): row for row in event_rows}
    material_by_event: dict[int, float] = {}
    for row in branch_rows:
        identifier = int(row["event_id"])
        material_by_event[identifier] = material_by_event.get(identifier, 0.0) + float(row["material_fraction_of_parent"])
    for identifier, event in event_by_id.items():
        expected = 1.0 - float(event["shell_loss_fraction_observed"])
        if not math.isclose(material_by_event[identifier], expected, rel_tol=0.0, abs_tol=2.0e-12):
            raise RuntimeError(f"event {identifier} violates declared branch material accounting")
        is_bio = str(event["observed_channel"]) == "biological"
        if bool(event["retained_mother_observed"]) != is_bio:
            raise RuntimeError(f"event {identifier} has an inconsistent structural channel readout")


def main() -> None:
    exposure_rows, event_rows, branch_rows = generate_synthetic_protocol()
    verify_protocol(exposure_rows, event_rows, branch_rows)
    biological_fit = fit_biological_rate(exposure_rows)
    hydrodynamic_fit = fit_hydrodynamic_rate(exposure_rows, biological_fit)
    hyd_kernel, bio_kernel, shell_loss = fit_kernels(branch_rows, event_rows)
    prediction_rows = make_prediction_rows(exposure_rows, biological_fit, hydrodynamic_fit)
    material_grid = np.linspace(1.0e-4, 0.9999, 500)
    kernel_density_rows: list[dict[str, object]] = []
    for channel, fitted, truth, fitted_shell, true_shell in (
        (
            "hydrodynamic",
            hyd_kernel,
            BetaKernelFit(TRUE_HYDRODYNAMIC_KERNEL_ALPHA, TRUE_HYDRODYNAMIC_KERNEL_BETA, math.nan),
            0.0,
            0.0,
        ),
        (
            "biological",
            bio_kernel,
            BetaKernelFit(TRUE_BIOLOGICAL_KERNEL_ALPHA, TRUE_BIOLOGICAL_KERNEL_BETA, math.nan),
            shell_loss,
            TRUE_BIOLOGICAL_SHELL_LOSS_FRACTION,
        ),
    ):
        fitted_density = beta_kernel_density(material_grid, fitted, channel, fitted_shell)
        generator_density = beta_kernel_density(material_grid, truth, channel, true_shell)
        for fraction, fitted_value, generator_value in zip(material_grid, fitted_density, generator_density):
            kernel_density_rows.append(
                {
                    "channel": channel,
                    "material_fraction_of_parent": float(fraction),
                    "fitted_daughter_kernel_density": float(fitted_value),
                    "generator_daughter_kernel_density": float(generator_value),
                }
            )
    lognormal_rows = local_lognormal_response(
        hydrodynamic_fit,
        biological_fit,
        hyd_kernel,
        bio_kernel,
        shell_loss,
        np.linspace(0.65, 1.85, 121),
    )
    metrics = metric_rows(
        prediction_rows,
        hydrodynamic_fit,
        biological_fit,
        hyd_kernel,
        bio_kernel,
        shell_loss,
        lognormal_rows,
    )

    write_csv(RESULTS / "biological_mechanical_stress_ramp_exposure.csv", exposure_rows)
    write_csv(RESULTS / "biological_mechanical_stress_ramp_events.csv", event_rows)
    write_csv(RESULTS / "biological_mechanical_stress_ramp_branches.csv", branch_rows)
    write_csv(RESULTS / "biological_mechanical_stress_ramp_rate_predictions.csv", prediction_rows)
    write_csv(RESULTS / "biological_mechanical_stress_ramp_kernel_density.csv", kernel_density_rows)
    write_csv(RESULTS / "biological_mechanical_stress_ramp_lognormal_response.csv", lognormal_rows)
    write_csv(RESULTS / "biological_mechanical_stress_ramp_metrics.csv", metrics)
    write_summary(
        exposure_rows,
        event_rows,
        metrics,
        biological_fit,
        hydrodynamic_fit,
        hyd_kernel,
        bio_kernel,
        shell_loss,
        lognormal_rows,
    )
    print(f"Wrote {RESULTS / 'biological_mechanical_stress_ramp_summary.md'}")
    print(
        "Recovered beta_bio,max={:.4f} d^-1, k_h={:.4f} d^-1, and q={:.4f}".format(
            biological_fit.maximum_rate_d_1,
            hydrodynamic_fit.prefactor_d_1,
            hydrodynamic_fit.exponent,
        )
    )


if __name__ == "__main__":
    main()
