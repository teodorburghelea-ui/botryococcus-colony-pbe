#!/usr/bin/env python3
"""Environmental and chemostat model-selection test for García-Cubero et al.

This is intentionally a small, transparent steady-state closure rather than a
new response-surface fit.  It asks the limited question that the García-Cubero
light--temperature--dilution design can answer: can a reduced environmental
state account for a size response to dilution that a uniform washout term
cannot produce?

The fit is restricted, before calculation, to nine Table-2 observations:
three centre replicates plus one high/low contrast for each of temperature,
light, and dilution.  The four published colony-size validation conditions
from Tables 3--4 are kept wholly out of fitting.  The other four nonzero
Table-2 rows are reported as an unused design check, not silently folded into
calibration.  The two zero-biomass/washout Table-2 rows have no defined colony
diameter and are excluded.

At a culture steady state the reduced state is

    q_* = theta_0 + theta_T T~ + theta_I I~ + theta_D D~,
    D_v = D_ref exp(q_*),

where T~, I~, and D~ are centred, scaled environmental coordinates and the
unobserved temporal form would be tau_q dq/dt = q_* - q.  The unbounded q is
an effective log-size-selection coordinate, not the bounded PBE coordinate a.
Only q_* is identifiable from steady measurements, so tau_q is deliberately
not fitted.
The dilution coefficient is an *environmental-state* route: in a full PBE it
must be represented in a state-dependent biological selection or retention
operator.  It is not inserted into the uniform -D_dil n washout term.

The nested uniform-washout null omits theta_D and is also refitted to the same
calibration observations.  A separate finite-volume calculation, with no
inlet or other state/size operator, verifies that applying only
n_i -> exp(-D_dil t) n_i leaves normalized volume distributions and their
volume-weighted diameter invariant to numerical precision.

Requirements: Python 3 with NumPy.  Run from the workspace root:

    python3 CODES/garcia_cubero_chemostat_point4.py

Then make the manuscript figure with:

    python3 CODES/plot_garcia_chemostat.py
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
RESULTS = HERE / "results"
SOURCE_DESIGN = (
    PROJECT_ROOT
    / "Comparison_Literature"
    / "Comparison_With_Garcia_Cubero"
    / "GarciaCubero2021_comparison_package"
    / "garcia_cubero_2021_table2_transcribed.csv"
)
OPTIMIZED_TARGETS = HERE / "data" / "garcia_cubero_2021_optimized_size_validation.csv"

# This is a declared, small calibration set: three Table-2 centre replicates
# plus paired contrasts for temperature (2/8), light (3/11), and dilution
# (5/7).  It spans the design without treating every available datum as a fit.
CALIBRATION_EXPERIMENTS = frozenset({1, 2, 3, 5, 6, 7, 8, 11, 15})

D_REF_UM = 150.0
ENVIRONMENT_CENTER = np.asarray((22.5, 1200.0, 0.20), dtype=float)
ENVIRONMENT_SCALE = np.asarray((7.5, 600.0, 0.10), dtype=float)
WASHOUT_RATES_D_INV = (0.10, 0.20, 0.30)
WASHOUT_DURATION_D = 4.0
PREDICTION_T_CRIT_DF5 = 2.570582


@dataclass(frozen=True)
class Observation:
    """One volume-weighted colony-size observation and its declared role."""

    condition_id: str
    source: str
    split: str
    temperature_C: float
    light_umol_m2_s: float
    dilution_d_1: float
    colony_size_um: float
    colony_size_sd_um: float | None


@dataclass(frozen=True)
class FittedClosure:
    """A nested reduced environmental-state closure fitted in log diameter."""

    name: str
    include_dilution_state: bool
    theta: np.ndarray

    @property
    def parameter_names(self) -> tuple[str, ...]:
        if self.include_dilution_state:
            return ("theta_0", "theta_temperature", "theta_light", "theta_dilution_state")
        return ("theta_0", "theta_temperature", "theta_light")

    def state_equilibrium(self, observation: Observation) -> float:
        coordinates = standardised_environment(observation)
        if self.include_dilution_state:
            return float(np.dot(np.r_[1.0, coordinates], self.theta))
        return float(np.dot(np.r_[1.0, coordinates[:2]], self.theta))

    def predict_diameter_um(self, observation: Observation) -> float:
        return D_REF_UM * math.exp(self.state_equilibrium(observation))


def standardised_environment(observation: Observation) -> np.ndarray:
    """Return centred dimensionless T, incident-light, and dilution inputs."""

    values = np.asarray(
        (observation.temperature_C, observation.light_umol_m2_s, observation.dilution_d_1),
        dtype=float,
    )
    return (values - ENVIRONMENT_CENTER) / ENVIRONMENT_SCALE


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_table2_design() -> list[Observation]:
    """Read nonzero steady cultures and attach the declared calibration split."""

    observations: list[Observation] = []
    for row in read_csv(SOURCE_DESIGN):
        experiment = int(float(row["experiment"]))
        diameter = float(row["colony_size_um"])
        # Rows 12 and 14 are washouts with zero biomass and therefore no
        # colony-size response to model.  They are not zero-diameter colonies.
        if diameter <= 0.0:
            continue
        split = "calibration" if experiment in CALIBRATION_EXPERIMENTS else "unused_design_check"
        observations.append(
            Observation(
                condition_id=f"table2_experiment_{experiment}",
                source="García-Cubero et al. (2021), Table 2",
                split=split,
                temperature_C=float(row["temperature_C"]),
                light_umol_m2_s=float(row["light_umol_m2_s"]),
                dilution_d_1=float(row["dilution_d_1"]),
                colony_size_um=diameter,
                colony_size_sd_um=None,
            )
        )
    expected = sorted(CALIBRATION_EXPERIMENTS)
    found = sorted(int(item.condition_id.rsplit("_", 1)[1]) for item in observations if item.split == "calibration")
    if found != expected:
        raise ValueError(f"declared Table-2 calibration subset changed: expected {expected}, found {found}")
    return observations


def read_optimized_targets() -> list[Observation]:
    """Read published optimizer validation observations kept outside fitting."""

    observations: list[Observation] = []
    for row in read_csv(OPTIMIZED_TARGETS):
        observations.append(
            Observation(
                condition_id=row["condition_id"],
                source=f"García-Cubero et al. (2021), {row['source_table']}",
                split="optimized_prediction",
                temperature_C=float(row["temperature_C"]),
                light_umol_m2_s=float(row["light_umol_m2_s"]),
                dilution_d_1=float(row["dilution_d_1"]),
                colony_size_um=float(row["colony_size_um"]),
                colony_size_sd_um=float(row["colony_size_sd_um"]),
            )
        )
    if len(observations) != 4:
        raise ValueError("the point-4 optimized prediction set must contain four published size targets")
    return observations


def matrix_for(observations: Iterable[Observation], include_dilution_state: bool) -> np.ndarray:
    rows: list[np.ndarray] = []
    for observation in observations:
        coordinate = standardised_environment(observation)
        if include_dilution_state:
            rows.append(np.r_[1.0, coordinate])
        else:
            rows.append(np.r_[1.0, coordinate[:2]])
    return np.asarray(rows, dtype=float)


def fit_closure(observations: list[Observation], *, include_dilution_state: bool, name: str) -> FittedClosure:
    """Fit the only identifiable steady-state quantities, log(Dv / D_ref)."""

    if len(observations) < (4 if include_dilution_state else 3):
        raise ValueError("insufficient declared calibration observations")
    response = np.log(np.asarray([item.colony_size_um / D_REF_UM for item in observations], dtype=float))
    theta, _, rank, _ = np.linalg.lstsq(matrix_for(observations, include_dilution_state), response, rcond=None)
    expected_rank = len(theta)
    if rank != expected_rank:
        raise ValueError("declared environmental calibration matrix is rank deficient")
    return FittedClosure(name=name, include_dilution_state=include_dilution_state, theta=theta)


def fit_constant_mean(observations: list[Observation]) -> FittedClosure:
    """Fit an arithmetic constant-diameter baseline on the calibration rows."""

    if not observations:
        raise ValueError("constant-mean baseline requires observations")
    return FittedClosure(
        name="constant_mean_baseline",
        include_dilution_state=False,
        theta=np.asarray((math.log(float(np.mean([item.colony_size_um for item in observations])) / D_REF_UM), 0.0, 0.0)),
    )


def design_matrix_with_constant(observations: list[Observation], include_dilution_state: bool) -> np.ndarray:
    """Design matrix for the log-diameter OLS diagnostics."""

    if include_dilution_state:
        return matrix_for(observations, True)
    return np.ones((len(observations), 1), dtype=float)


def fit_log_ols(observations: list[Observation], include_dilution_state: bool) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Return coefficients, design matrix, residual variance and inverse cross-product."""

    X = design_matrix_with_constant(observations, include_dilution_state)
    y = np.log(np.asarray([item.colony_size_um / D_REF_UM for item in observations], dtype=float))
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != X.shape[1]:
        raise ValueError("diagnostic design matrix is rank deficient")
    residual = y - X @ beta
    dof = len(observations) - X.shape[1]
    if dof <= 0:
        raise ValueError("diagnostic residual degrees of freedom must be positive")
    sigma2 = float(np.sum(residual**2) / dof)
    xtx_inv = np.linalg.inv(X.T @ X)
    return beta, X, sigma2, xtx_inv


def prediction_interval_records(model_name: str, observations: list[Observation], targets: list[Observation], include_dilution_state: bool) -> list[dict[str, float | str]]:
    """95% log-scale predictive intervals, including residual prediction variance."""

    beta, _, sigma2, xtx_inv = fit_log_ols(observations, include_dilution_state)
    records: list[dict[str, float | str]] = []
    for target in targets:
        xrow = design_matrix_with_constant([target], include_dilution_state)[0]
        log_prediction = math.log(D_REF_UM) + float(xrow @ beta)
        leverage = float(xrow @ xtx_inv @ xrow)
        se_prediction = math.sqrt(sigma2 * (1.0 + leverage))
        center = math.exp(log_prediction)
        records.append({
            "model": model_name,
            "condition_id": target.condition_id,
            "predicted_um": center,
            "prediction_interval_low_um": math.exp(log_prediction - PREDICTION_T_CRIT_DF5 * se_prediction),
            "prediction_interval_high_um": math.exp(log_prediction + PREDICTION_T_CRIT_DF5 * se_prediction),
            "design_leverage": leverage,
            "residual_sd_log": math.sqrt(sigma2),
            "interval_level": 0.95,
        })
    return records


def leave_one_out_records(observations: list[Observation], include_dilution_state: bool, model_name: str) -> list[dict[str, float | str]]:
    """Leave-one-condition-out predictions for the prespecified nine-row design."""

    records: list[dict[str, float | str]] = []
    for index, held_out in enumerate(observations):
        training = observations[:index] + observations[index + 1:]
        if model_name == "constant_mean_baseline":
            model = fit_constant_mean(training)
        else:
            model = fit_closure(training, include_dilution_state=include_dilution_state, name=model_name)
        predicted = model.predict_diameter_um(held_out)
        records.append({
            "model": model_name,
            "condition_id": held_out.condition_id,
            "observed_um": held_out.colony_size_um,
            "predicted_um": predicted,
            "residual_um": predicted - held_out.colony_size_um,
        })
    return records


def influence_records(observations: list[Observation], include_dilution_state: bool, model_name: str) -> list[dict[str, float | str]]:
    """Hat leverage and Cook-style influence diagnostics in log-diameter space."""

    beta, X, sigma2, xtx_inv = fit_log_ols(observations, include_dilution_state)
    y = np.log(np.asarray([item.colony_size_um / D_REF_UM for item in observations], dtype=float))
    residual = y - X @ beta
    p = X.shape[1]
    records: list[dict[str, float | str]] = []
    for item, row, err in zip(observations, X, residual):
        leverage = float(row @ xtx_inv @ row)
        cooks = float((err**2 / (p * sigma2)) * (leverage / max((1.0 - leverage) ** 2, 1.0e-15)))
        records.append({
            "model": model_name,
            "condition_id": item.condition_id,
            "leverage": leverage,
            "studentized_log_residual": float(err / math.sqrt(sigma2 * max(1.0 - leverage, 1.0e-15))),
            "cooks_distance": cooks,
        })
    return records


def metric_record(model: FittedClosure, split: str, observations: list[Observation]) -> dict[str, float | int | str]:
    observed = np.asarray([item.colony_size_um for item in observations], dtype=float)
    predicted = np.asarray([model.predict_diameter_um(item) for item in observations], dtype=float)
    residual = predicted - observed
    centered_sum_squares = float(np.sum((observed - observed.mean()) ** 2))
    r_squared = float("nan") if centered_sum_squares == 0.0 else 1.0 - float(np.sum(residual**2)) / centered_sum_squares
    log_residual = np.log(predicted / observed)
    return {
        "model": model.name,
        "split": split,
        "n_observations": len(observations),
        "n_fitted_parameters": len(model.theta),
        "rmse_um": float(np.sqrt(np.mean(residual**2))),
        "mae_um": float(np.mean(np.abs(residual))),
        "mean_bias_um": float(np.mean(residual)),
        "r_squared": r_squared,
        "log_rss": float(np.sum(log_residual**2)),
    }


def prediction_records(models: list[FittedClosure], observations: list[Observation]) -> list[dict[str, float | str]]:
    records: list[dict[str, float | str]] = []
    for model in models:
        for observation in observations:
            predicted = model.predict_diameter_um(observation)
            records.append(
                {
                    "model": model.name,
                    "split": observation.split,
                    "condition_id": observation.condition_id,
                    "source": observation.source,
                    "temperature_C": observation.temperature_C,
                    "light_umol_m2_s": observation.light_umol_m2_s,
                    "dilution_d_1": observation.dilution_d_1,
                    "observed_volume_weighted_diameter_um": observation.colony_size_um,
                    "observed_sd_um": "" if observation.colony_size_sd_um is None else observation.colony_size_sd_um,
                    "predicted_volume_weighted_diameter_um": predicted,
                    "residual_um": predicted - observation.colony_size_um,
                    "state_equilibrium": model.state_equilibrium(observation),
                }
            )
    return records


def named_observation(temperature_C: float, light_umol_m2_s: float, dilution_d_1: float) -> Observation:
    """Construct a no-target condition for the declared dilution contrast."""

    return Observation(
        condition_id="state_dilution_contrast",
        source="declared point-4 numerical contrast",
        split="numerical_contrast",
        temperature_C=temperature_C,
        light_umol_m2_s=light_umol_m2_s,
        dilution_d_1=dilution_d_1,
        colony_size_um=float("nan"),
        colony_size_sd_um=None,
    )


def dilution_response_records(models: list[FittedClosure], design_observations: list[Observation]) -> list[dict[str, float | str]]:
    """Report the matched 22.5 C, 1800-light dilution contrast used in Table 2."""

    records: list[dict[str, float | str]] = []
    for model in models:
        for dilution in WASHOUT_RATES_D_INV:
            condition = named_observation(22.5, 1800.0, dilution)
            matched = next(
                (
                    item
                    for item in design_observations
                    if np.isclose(item.temperature_C, condition.temperature_C)
                    and np.isclose(item.light_umol_m2_s, condition.light_umol_m2_s)
                    and np.isclose(item.dilution_d_1, condition.dilution_d_1)
                ),
                None,
            )
            records.append(
                {
                    "model": model.name,
                    "temperature_C": condition.temperature_C,
                    "light_umol_m2_s": condition.light_umol_m2_s,
                    "dilution_d_1": dilution,
                    "predicted_volume_weighted_diameter_um": model.predict_diameter_um(condition),
                    "state_equilibrium": model.state_equilibrium(condition),
                    "observed_volume_weighted_diameter_um": "" if matched is None else matched.colony_size_um,
                    "observation_id": "" if matched is None else matched.condition_id,
                }
            )
    return records


def uniform_washout_control() -> tuple[list[dict[str, float]], list[dict[str, float | str]]]:
    """Numerically verify normalized-distribution invariance under -D n alone."""

    edges = np.geomspace(20.0, 1000.0, 241)
    diameter = np.sqrt(edges[:-1] * edges[1:])
    # The precise initial log-normal is immaterial: uniform washout multiplies
    # every size section by the same positive scalar.  It is deliberately
    # broad enough for the displayed overlay to make that invariance visible.
    number_initial = np.exp(-0.5 * ((np.log(diameter) - np.log(170.0)) / 0.38) ** 2)
    distributions: dict[float, np.ndarray] = {}
    diameters: dict[float, float] = {}
    records: list[dict[str, float]] = []
    for washout in WASHOUT_RATES_D_INV:
        number = number_initial * math.exp(-washout * WASHOUT_DURATION_D)
        normalized_volume = number * diameter**3
        normalized_volume /= normalized_volume.sum()
        volume_weighted_diameter = float(np.sum(number * diameter**4) / np.sum(number * diameter**3))
        distributions[washout] = normalized_volume
        diameters[washout] = volume_weighted_diameter
        for section_diameter, fraction in zip(diameter, normalized_volume):
            records.append(
                {
                    "washout_d_1": washout,
                    "duration_d": WASHOUT_DURATION_D,
                    "diameter_um": float(section_diameter),
                    "normalized_volume_fraction": float(fraction),
                    "volume_weighted_diameter_um": volume_weighted_diameter,
                }
            )

    reference_rate = WASHOUT_RATES_D_INV[0]
    metrics: list[dict[str, float | str]] = []
    for washout in WASHOUT_RATES_D_INV[1:]:
        total_variation = 0.5 * float(np.sum(np.abs(distributions[reference_rate] - distributions[washout])))
        diameter_difference = abs(diameters[reference_rate] - diameters[washout])
        metrics.append(
            {
                "comparison": f"uniform_washout_{reference_rate:.2f}_vs_{washout:.2f}_d_1",
                "maximum_normalized_distribution_total_variation": total_variation,
                "absolute_volume_weighted_diameter_difference_um": diameter_difference,
                "relative_volume_weighted_diameter_difference": diameter_difference / diameters[reference_rate],
            }
        )
    if max(float(item["maximum_normalized_distribution_total_variation"]) for item in metrics) > 1.0e-12:
        raise RuntimeError("uniform washout changed the normalized size distribution")
    if max(float(item["relative_volume_weighted_diameter_difference"]) for item in metrics) > 1.0e-12:
        raise RuntimeError("uniform washout changed the volume-weighted diameter")
    return records, metrics


def write_csv(path: Path, records: list[dict]) -> None:
    if not records:
        raise ValueError(f"refusing to write empty output {path.name}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def parameter_records(models: list[FittedClosure]) -> list[dict[str, float | str]]:
    records: list[dict[str, float | str]] = []
    for model in models:
        records.append(
            {
                "model": model.name,
                "parameter": "D_ref",
                "value": D_REF_UM,
                "units": "um",
                    "meaning": "fixed diameter reference in D_v = D_ref exp(q_*)",
            }
        )
        for name, value in zip(model.parameter_names, model.theta):
            records.append(
                {
                    "model": model.name,
                    "parameter": name,
                    "value": float(value),
                    "units": "dimensionless",
                    "meaning": "coefficient in the steady environmental state q_*",
                }
            )
    return records


def get_metric(metrics: list[dict[str, float | int | str]], model: str, split: str) -> dict[str, float | int | str]:
    return next(item for item in metrics if item["model"] == model and item["split"] == split)


def write_extended_diagnostics(
    metrics: list[dict[str, float | int | str]],
    loo_state: list[dict[str, float | str]],
    loo_null: list[dict[str, float | str]],
    loo_mean: list[dict[str, float | str]],
    influence_state: list[dict[str, float | str]],
    prediction_intervals: list[dict[str, float | str]],
) -> None:
    """Write baseline, LOOCV, influence and predictive-interval summaries."""

    def loo_rmse(records: list[dict[str, float | str]]) -> float:
        return math.sqrt(float(np.mean([float(item["residual_um"]) ** 2 for item in records])))

    state_pi = [item for item in prediction_intervals if item["model"] == "environmental_state_closure"]
    max_influence = max(influence_state, key=lambda item: float(item["cooks_distance"]))
    max_leverage = max(influence_state, key=lambda item: float(item["leverage"]))
    baseline = get_metric(metrics, "constant_mean_baseline", "calibration")
    lines = [
        "# Continuous-culture comparator diagnostics",
        "",
        "The nine-row calibration subset was prespecified from the process design before model fitting or inspection of model performance: three centre replicates plus paired high/low contrasts for temperature (2/8), incident light (3/11), and dilution (5/7). The four remaining nonzero Table-2 rows and all four optimized-condition targets were withheld.",
        "",
        f"The constant arithmetic-mean baseline (calibration mean {float(np.mean([item.colony_size_um for item in read_table2_design() if item.split == 'calibration'])):.1f} micrometres) has calibration RMSE {float(baseline['rmse_um']):.1f} micrometres, unused-design-check RMSE {float(get_metric(metrics, 'constant_mean_baseline', 'unused_design_check')['rmse_um']):.1f} micrometres, and optimized-condition RMSE {float(get_metric(metrics, 'constant_mean_baseline', 'optimized_prediction')['rmse_um']):.1f} micrometres.",
        f"Leave-one-condition-out calibration RMSE is {loo_rmse(loo_state):.1f} micrometres for the environmental-state closure, {loo_rmse(loo_null):.1f} micrometres for the dilution-free null, and {loo_rmse(loo_mean):.1f} micrometres for the constant-mean baseline.",
        f"Environmental-state calibration leverage ranges from {min(float(item['leverage']) for item in influence_state):.3f} to {max(float(item['leverage']) for item in influence_state):.3f}; the largest leverage is {max_leverage['condition_id']} and the largest Cook-style influence is {max_influence['condition_id']} (D={float(max_influence['cooks_distance']):.3f}). These diagnostics flag influential design points; they do not establish independent biological replicates.",
        "The 95% predictive intervals below are log-scale OLS intervals including residual variance and are diagnostic, not confidence intervals for a universal response surface.",
    ]
    for item in state_pi:
        lines.append(f"- {item['condition_id']}: {float(item['predicted_um']):.1f} micrometres [{float(item['prediction_interval_low_um']):.1f}, {float(item['prediction_interval_high_um']):.1f}].")
    (RESULTS / "garcia_chemostat_extended_diagnostics.md").write_text("\n".join(lines) + "\n")


def write_summary(
    metrics: list[dict[str, float | int | str]],
    closures: list[FittedClosure],
    dilution_records: list[dict[str, float | str]],
    invariance_metrics: list[dict[str, float | str]],
) -> None:
    state_name = "environmental_state_closure"
    null_name = "uniform_washout_null"
    state_train = get_metric(metrics, state_name, "calibration")
    state_optimized = get_metric(metrics, state_name, "optimized_prediction")
    state_unused = get_metric(metrics, state_name, "unused_design_check")
    null_train = get_metric(metrics, null_name, "calibration")
    null_optimized = get_metric(metrics, null_name, "optimized_prediction")
    state = next(item for item in closures if item.name == state_name)
    state_curve = [
        item
        for item in dilution_records
        if item["model"] == state_name
        and (np.isclose(float(item["dilution_d_1"]), 0.10) or np.isclose(float(item["dilution_d_1"]), 0.30))
    ]
    state_curve.sort(key=lambda item: float(item["dilution_d_1"]))
    dilution_ratio = float(state_curve[-1]["predicted_volume_weighted_diameter_um"]) / float(
        state_curve[0]["predicted_volume_weighted_diameter_um"]
    )
    max_tv = max(float(item["maximum_normalized_distribution_total_variation"]) for item in invariance_metrics)
    max_relative_dv = max(float(item["relative_volume_weighted_diameter_difference"]) for item in invariance_metrics)

    text = f"""# García-Cubero environmental and chemostat test (Point 4)

## Protocol

- Calibration is fixed before numerical fitting: nine nonzero Table-2 observations
  (three centre replicates and the temperature, light, and dilution contrasts
  2/8, 3/11, and 5/7).
- Four remaining nonzero Table-2 observations are an unused design check.
- All four published optimization validations from Tables 3--4 are reserved
  for prediction. They do not enter either fit.
- The reduced environmental state has four fitted equilibrium coefficients:
  `q_* = theta_0 + theta_T T~ + theta_I I~ + theta_D D~` and
  `D_v = 150 exp(q_*)` micrometres. Its unobserved relaxation time is not fit
  from steady-state data.
- The nested washout-only null is refit with `theta_D = 0`. A separate
  sectional calculation, with no inlet or other state/size operator, applies
  only `n_i -> exp(-D t)n_i`.

## Results

- Environmental-state closure: calibration RMSE = {float(state_train['rmse_um']):.1f} um;
  unused-design-check RMSE = {float(state_unused['rmse_um']):.1f} um; held-out optimized-condition
  RMSE = {float(state_optimized['rmse_um']):.1f} um.
- Uniform-washout null: calibration RMSE = {float(null_train['rmse_um']):.1f} um;
  held-out optimized-condition RMSE = {float(null_optimized['rmse_um']):.1f} um.
- Fitted environmental-state coefficients are theta_0 = {state.theta[0]:.4f},
  theta_T = {state.theta[1]:.4f}, theta_I = {state.theta[2]:.4f}, and
  theta_D = {state.theta[3]:.4f}. The declared 22.5 C, 1800 umol m^-2 s^-1
  dilution contrast changes the state-closure prediction by a factor of
  {dilution_ratio:.2f} from D = 0.10 to 0.30 d^-1; the washout-only null is
  exactly diameter-invariant by construction.
- The independent sectional washout control gives maximum normalized-volume
  total variation {max_tv:.3e} and relative volume-weighted-diameter change
  {max_relative_dv:.3e} across 0.10--0.30 d^-1.

## Interpretation boundary

The positive dilution response is routed explicitly through `q_*`, which in a
full PBE must mean state-dependent biological size selection, state-dependent
retention, or both. It is not evidence that uniform washout changes a
normalized distribution. The small design and the weak unused cold-condition
check do not identify a unique temperature law, relaxation time, or retention
mechanism. This is a structural/calibration test, not a replacement for a
time-resolved chemostat PBE fit.
"""
    (RESULTS / "garcia_chemostat_summary.md").write_text(text)


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    design = read_table2_design()
    optimized = read_optimized_targets()
    calibration = [item for item in design if item.split == "calibration"]
    unused_design = [item for item in design if item.split == "unused_design_check"]
    if len(calibration) != 9 or len(unused_design) != 4:
        raise ValueError("unexpected García-Cubero Table-2 split")

    environmental_state = fit_closure(
        calibration,
        include_dilution_state=True,
        name="environmental_state_closure",
    )
    washout_null = fit_closure(
        calibration,
        include_dilution_state=False,
        name="uniform_washout_null",
    )
    constant_mean = fit_constant_mean(calibration)
    models = [environmental_state, washout_null, constant_mean]
    all_observations = design + optimized
    metrics: list[dict[str, float | int | str]] = []
    for model in models:
        for split, subset in (
            ("calibration", calibration),
            ("unused_design_check", unused_design),
            ("all_nonzero_table2", design),
            ("optimized_prediction", optimized),
        ):
            metrics.append(metric_record(model, split, subset))

    dilution_records = dilution_response_records([environmental_state, washout_null], design)
    washout_records, invariance_metrics = uniform_washout_control()
    loo_state = leave_one_out_records(calibration, True, environmental_state.name)
    loo_null = leave_one_out_records(calibration, False, washout_null.name)
    loo_mean = leave_one_out_records(calibration, False, constant_mean.name)
    influence_state = influence_records(calibration, True, environmental_state.name)
    influence_null = influence_records(calibration, False, washout_null.name)
    prediction_intervals = (
        prediction_interval_records(environmental_state.name, calibration, optimized, True)
        + prediction_interval_records(washout_null.name, calibration, optimized, False)
        + prediction_interval_records(constant_mean.name, calibration, optimized, False)
    )
    write_csv(RESULTS / "garcia_chemostat_predictions.csv", prediction_records(models, all_observations))
    write_csv(RESULTS / "garcia_chemostat_metrics.csv", metrics)
    write_csv(RESULTS / "garcia_chemostat_parameters.csv", parameter_records(models))
    write_csv(RESULTS / "garcia_chemostat_dilution_response.csv", dilution_records)
    write_csv(RESULTS / "garcia_chemostat_uniform_washout.csv", washout_records)
    write_csv(RESULTS / "garcia_chemostat_washout_invariance.csv", invariance_metrics)
    write_csv(RESULTS / "garcia_chemostat_leave_one_condition_out.csv", loo_state + loo_null + loo_mean)
    write_csv(RESULTS / "garcia_chemostat_influence_diagnostics.csv", influence_state + influence_null)
    write_csv(RESULTS / "garcia_chemostat_prediction_intervals.csv", prediction_intervals)
    write_summary(metrics, models, dilution_records, invariance_metrics)
    write_extended_diagnostics(metrics, loo_state, loo_null, loo_mean, influence_state, prediction_intervals)

    state_optimized = get_metric(metrics, "environmental_state_closure", "optimized_prediction")
    null_optimized = get_metric(metrics, "uniform_washout_null", "optimized_prediction")
    print("Point-4 García-Cubero environmental/chemostat test complete")
    print(f"  state closure held-out optimized RMSE: {float(state_optimized['rmse_um']):.2f} um")
    print(f"  washout-only held-out optimized RMSE: {float(null_optimized['rmse_um']):.2f} um")
    print(f"  wrote reproducible outputs under {RESULTS}")


if __name__ == "__main__":
    main()
