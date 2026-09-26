#!/usr/bin/env python3
"""Reproduce the Zhang--Kojima calibration/prediction hierarchy.

The script deliberately keeps the Fig. 5 transient calibration distinct from
the Fig. 6 quasi-steady diagnostic.  Only three light--separation parameters
are fitted, and only to the 3-klx Fig. 5 trajectories.  The 10-klx Fig. 5
trajectories and all non-initial Fig. 3 histograms are evaluated without
refitting.  Source digitizations remain in the literature-comparison package;
this script writes all derived, reproducible test results to CODES/results.

Requirements: Python 3, NumPy, SciPy.  No biological parameter is inferred
from the Fig. 3 distributions.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
RESULTS = HERE / "results"
SOURCE_DATA = (
    PROJECT_ROOT
    / "Comparison_Literature"
    / "ComparisonWithLiterature_complete_package"
    / "ComparisonWithLiterature"
    / "data"
)
FIG6_SOURCE = (
    PROJECT_ROOT
    / "Comparison_Literature"
    / "Comparison_Zhang_Koshima_complete_package"
    / "ZK98_Fig6_PBE_predictions.csv"
)


@dataclass(frozen=True)
class FitParameters:
    """The only adjustable terms in the minimal Fig. 5 submodel."""

    dc_min: float
    dc_max: float
    kc: float


@dataclass(frozen=True)
class FixedParameters:
    """Predeclared baseline kinetics retained fixed during the fit."""

    gmax: float = 0.055
    kg: float = 0.70
    beta_max: float = 0.80
    sharpness: float = 8.0


FIXED = FixedParameters()
FIT_PRECULTURE = 3
PREDICTION_PRECULTURE = 10
VOLUMES = (5, 25, 50, 100)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def as_int(row: dict[str, str], key: str) -> int:
    return int(float(row[key]))


class LightHistory:
    """Source-specific Zhang--Kojima light reconstruction from Table 1."""

    def __init__(self, rows: list[dict[str, str]]) -> None:
        self.series: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
        grouped: dict[tuple[int, int], list[tuple[float, float]]] = defaultdict(list)
        for row in rows:
            key = (
                as_int(row, "preculture_irradiance_klx"),
                as_int(row, "lighted_volume_percent"),
            )
            grouped[key].append((as_float(row, "days"), as_float(row, "dry_weight_kg_m3")))
        for key, values in grouped.items():
            values.sort()
            self.series[key] = (
                np.asarray([value[0] for value in values]),
                np.asarray([value[1] for value in values]),
            )

    def dry_weight(self, time: float, preculture: int, volume: int) -> float:
        times, biomass = self.series[(preculture, volume)]
        if time <= times[0]:
            return float(biomass[0])
        if time >= times[-1]:
            slope = (biomass[-1] - biomass[-2]) / (times[-1] - times[-2])
            return float(biomass[-1] + slope * (time - times[-1]))
        return float(np.interp(time, times, biomass))

    def average_light(self, time: float, preculture: int, volume: int) -> float:
        i100 = 6.05 * math.exp(-0.290 * self.dry_weight(time, preculture, volume))
        return volume / 100.0 * i100


class DiameterPBE:
    """Minimal one-coordinate growth--equal-volume-breakup PBE.

    The implementation follows the pre-existing baseline comparison: a
    log-diameter finite-volume grid, positive SSPRK2 updates, material-aware
    equal-volume daughter interpolation, and a prescribed biomass-derived
    light history.  It is deliberately not the structured Eq. (3) model.
    """

    def __init__(self, n_sections: int = 120, dmin: float = 0.015, dmax: float = 0.80) -> None:
        self.edges = np.exp(np.linspace(math.log(dmin), math.log(dmax), n_sections + 1))
        self.diameter = np.sqrt(self.edges[:-1] * self.edges[1:])
        self.width = np.diff(self.edges)
        self.left, self.right, self.wleft, self.wright = self._daughter_map()

    def _daughter_map(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        n = len(self.diameter)
        left = np.zeros(n, dtype=int)
        right = np.zeros(n, dtype=int)
        wleft = np.zeros(n)
        wright = np.zeros(n)
        for index, parent in enumerate(self.diameter):
            daughter = parent / 2.0 ** (1.0 / 3.0)
            if daughter <= self.diameter[0]:
                left[index] = right[index] = 0
                wleft[index] = parent**3 / self.diameter[0] ** 3
                continue
            bin_index = int(np.clip(np.searchsorted(self.diameter, daughter, side="right") - 1, 0, n - 2))
            dl, dr = self.diameter[bin_index], self.diameter[bin_index + 1]
            wr = (parent**3 - 2.0 * dl**3) / (dr**3 - dl**3)
            left[index], right[index] = bin_index, bin_index + 1
            wleft[index], wright[index] = max(2.0 - wr, 0.0), max(wr, 0.0)
        return left, right, wleft, wright

    def volume_diameter(self, population: np.ndarray) -> float:
        return float((np.dot(self.diameter**3, population) / np.sum(population)) ** (1.0 / 3.0))

    def lognormal_initial(self, target_dv: float, sigma: float = 0.12) -> np.ndarray:
        median = target_dv
        population = np.empty_like(self.diameter)
        for _ in range(3):
            population = (
                np.exp(-0.5 * ((np.log(self.diameter) - math.log(median)) / sigma) ** 2)
                / (self.diameter * sigma * math.sqrt(2.0 * math.pi))
                * self.width
            )
            population /= population.sum()
            median *= target_dv / self.volume_diameter(population)
        return population

    def coefficients(self, light: float, parameters: FitParameters) -> tuple[np.ndarray, np.ndarray]:
        growth_faces = FIXED.gmax * light / (FIXED.kg + light) * self.edges
        critical = parameters.dc_min + (parameters.dc_max - parameters.dc_min) * light / (parameters.kc + light)
        z = (self.diameter / critical) ** FIXED.sharpness
        breakup = FIXED.beta_max * z / (1.0 + z)
        return growth_faces, breakup

    def rhs(self, population: np.ndarray, light: float, parameters: FitParameters) -> np.ndarray:
        growth_faces, breakup = self.coefficients(light, parameters)
        flux = np.zeros(len(population) + 1)
        flux[1:] = growth_faces[1:] * population / self.width
        change = flux[:-1] - flux[1:]
        loss = breakup * population
        change -= loss
        np.add.at(change, self.left, self.wleft * loss)
        np.add.at(change, self.right, self.wright * loss)
        return change

    def timestep(self, light: float, parameters: FitParameters) -> float:
        growth_faces, breakup = self.coefficients(light, parameters)
        growth_limit = 0.35 * np.min(self.width / growth_faces[1:])
        breakup_limit = 0.35 / np.max(breakup)
        return float(min(0.20, growth_limit, breakup_limit))

    def simulate(
        self,
        history: LightHistory,
        preculture: int,
        volume: int,
        times: list[float],
        parameters: FitParameters,
        initial: np.ndarray,
    ) -> dict[float, np.ndarray]:
        population = initial.copy()
        current = times[0]
        output = {current: population.copy()}
        for target in times[1:]:
            while current < target - 1e-12:
                light_1 = history.average_light(current, preculture, volume)
                dt = min(self.timestep(light_1, parameters), target - current)
                derivative_1 = self.rhs(population, light_1, parameters)
                stage = np.maximum(population + dt * derivative_1, 0.0)
                light_2 = history.average_light(current + dt, preculture, volume)
                derivative_2 = self.rhs(stage, light_2, parameters)
                population = np.maximum(0.5 * population + 0.5 * (stage + dt * derivative_2), 0.0)
                current += dt
            output[target] = population.copy()
        return output


def group_fig5(rows: list[dict[str, str]]) -> dict[tuple[int, int], list[dict[str, str]]]:
    grouped: dict[tuple[int, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(as_int(row, "preculture_irradiance_klx"), as_int(row, "lighted_volume_percent"))].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: as_float(row, "days"))
    return grouped


def admissible(values: np.ndarray) -> bool:
    return (
        0.02 <= values[0] <= 0.18
        and 0.08 <= values[1] <= 1.00
        and 0.03 <= values[2] <= 12.0
        and values[1] > values[0]
    )


def fit_fig5(
    model: DiameterPBE,
    history: LightHistory,
    grouped: dict[tuple[int, int], list[dict[str, str]]],
) -> FitParameters:
    """Fit only the 3-klx Fig. 5 transient trajectories."""

    def objective(values: np.ndarray) -> float:
        if not admissible(values):
            return 1.0e3
        parameters = FitParameters(*map(float, values))
        errors: list[float] = []
        for volume in VOLUMES:
            series = grouped[(FIT_PRECULTURE, volume)]
            times = [as_float(row, "days") for row in series]
            observations = np.asarray([as_float(row, "volume_averaged_diameter_mm") for row in series])
            snapshots = model.simulate(
                history,
                FIT_PRECULTURE,
                volume,
                times,
                parameters,
                model.lognormal_initial(float(observations[0])),
            )
            predictions = np.asarray([model.volume_diameter(snapshots[time]) for time in times])
            errors.extend((predictions - observations) ** 2)
        return float(np.mean(errors))

    result = minimize(
        objective,
        x0=np.asarray([0.040, 0.90, 10.0]),
        method="Nelder-Mead",
        options={"maxiter": 350, "xatol": 1e-8, "fatol": 1e-11},
    )
    if not result.success:
        raise RuntimeError(f"Fig. 5 optimization did not converge: {result.message}")
    if not admissible(result.x):
        raise RuntimeError("Fig. 5 optimizer returned inadmissible parameters")
    return FitParameters(*map(float, result.x))


def regression_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = predicted - observed
    denominator = np.sum((observed - observed.mean()) ** 2)
    return {
        "n_points": float(len(observed)),
        "rmse_mm": float(np.sqrt(np.mean(residual**2))),
        "mae_mm": float(np.mean(np.abs(residual))),
        "mape_percent": float(100.0 * np.mean(np.abs(residual / observed))),
        "r2": float(1.0 - np.sum(residual**2) / denominator),
    }


def fig5_predictions(
    model: DiameterPBE,
    history: LightHistory,
    grouped: dict[tuple[int, int], list[dict[str, str]]],
    parameters: FitParameters,
) -> tuple[list[dict[str, float | int | str]], list[dict[str, float | int | str]]]:
    rows: list[dict[str, float | int | str]] = []
    summaries: list[dict[str, float | int | str]] = []
    for preculture in (FIT_PRECULTURE, PREDICTION_PRECULTURE):
        observed_all: list[float] = []
        predicted_all: list[float] = []
        for volume in VOLUMES:
            series = grouped[(preculture, volume)]
            times = [as_float(row, "days") for row in series]
            observations = np.asarray([as_float(row, "volume_averaged_diameter_mm") for row in series])
            snapshots = model.simulate(
                history,
                preculture,
                volume,
                times,
                parameters,
                model.lognormal_initial(float(observations[0])),
            )
            predictions = np.asarray([model.volume_diameter(snapshots[time]) for time in times])
            for source, prediction in zip(series, predictions):
                observed = as_float(source, "volume_averaged_diameter_mm")
                rows.append(
                    {
                        "split": "calibration" if preculture == FIT_PRECULTURE else "prediction",
                        "panel": source["panel"],
                        "preculture_irradiance_klx": preculture,
                        "lighted_volume_percent": volume,
                        "days": as_float(source, "days"),
                        "observed_dv_mm": observed,
                        "predicted_dv_mm": float(prediction),
                        "residual_mm": float(prediction - observed),
                    }
                )
            observed_all.extend(observations)
            predicted_all.extend(predictions)
        metrics = regression_metrics(np.asarray(observed_all), np.asarray(predicted_all))
        summaries.append(
            {
                "split": "calibration_3klx" if preculture == FIT_PRECULTURE else "prediction_10klx",
                "preculture_irradiance_klx": preculture,
                **metrics,
            }
        )
    return rows, summaries


def group_fig3(rows: list[dict[str, str]]) -> dict[tuple[int, int, float], list[dict[str, str]]]:
    grouped: dict[tuple[int, int, float], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (
            as_int(row, "preculture_irradiance_klx"),
            as_int(row, "lighted_volume_percent"),
            as_float(row, "days"),
        )
        grouped[key].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: as_float(row, "diameter_bin_center_mm"))
    return grouped


def histogram_initial(model: DiameterPBE, source: list[dict[str, str]], sigma: float = 0.012) -> np.ndarray:
    population = np.zeros_like(model.diameter)
    for row in source:
        weight = as_float(row, "frequency_percent") / 100.0
        if weight <= 0:
            continue
        center = as_float(row, "diameter_bin_center_mm")
        kernel = np.exp(-0.5 * ((model.diameter - center) / sigma) ** 2)
        kernel /= kernel.sum()
        population += weight * kernel
    return population / population.sum()


def sectional_histogram(model: DiameterPBE, population: np.ndarray, centers: np.ndarray) -> np.ndarray:
    values = np.zeros(len(centers))
    for index, center in enumerate(centers):
        lower, upper = index * 0.05, (index + 1) * 0.05
        values[index] = population[(model.diameter >= lower) & (model.diameter < upper)].sum()
    return 100.0 * values / values.sum()


def distribution_metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    difference = predicted - observed
    return {
        "distribution_rmse_percent": float(np.sqrt(np.mean(difference**2))),
        "total_variation_distance": float(0.5 * np.sum(np.abs(difference)) / 100.0),
        "wasserstein_distance_mm": float(np.sum(np.abs(np.cumsum(observed / 100.0) - np.cumsum(predicted / 100.0))) * 0.05),
    }


def fig3_predictions(
    model: DiameterPBE,
    history: LightHistory,
    grouped: dict[tuple[int, int, float], list[dict[str, str]]],
    parameters: FitParameters,
) -> tuple[list[dict[str, float | int | str]], list[dict[str, float | int | str]], list[dict[str, float | int | str]]]:
    details: list[dict[str, float | int | str]] = []
    metrics_rows: list[dict[str, float | int | str]] = []
    summary_rows: list[dict[str, float | int | str]] = []
    for preculture in (FIT_PRECULTURE, PREDICTION_PRECULTURE):
        all_metrics: list[dict[str, float]] = []
        for volume in VOLUMES:
            available_times = sorted(time for pc, vl, time in grouped if pc == preculture and vl == volume)
            initial_time = available_times[0]
            initial = histogram_initial(model, grouped[(preculture, volume, initial_time)])
            snapshots = model.simulate(history, preculture, volume, available_times, parameters, initial)
            for time in available_times:
                source = grouped[(preculture, volume, time)]
                centers = np.asarray([as_float(row, "diameter_bin_center_mm") for row in source])
                observed = np.asarray([as_float(row, "frequency_percent") for row in source])
                predicted = sectional_histogram(model, snapshots[time], centers)
                for center, experimental, calculated in zip(centers, observed, predicted):
                    details.append(
                        {
                            "preculture_irradiance_klx": preculture,
                            "lighted_volume_percent": volume,
                            "days": time,
                            "diameter_bin_center_mm": float(center),
                            "experimental_frequency_percent": float(experimental),
                            "predicted_frequency_percent": float(calculated),
                            "residual_percent": float(calculated - experimental),
                        }
                    )
                if time != initial_time:
                    measures = distribution_metrics(observed, predicted)
                    metrics_rows.append(
                        {
                            "split": "no_refit_same_preculture" if preculture == FIT_PRECULTURE else "no_refit_independent_preculture",
                            "preculture_irradiance_klx": preculture,
                            "lighted_volume_percent": volume,
                            "days": time,
                            **measures,
                        }
                    )
                    all_metrics.append(measures)
        summary_rows.append(
            {
                "split": "no_refit_3klx_distributions" if preculture == FIT_PRECULTURE else "no_refit_10klx_distributions",
                "preculture_irradiance_klx": preculture,
                "n_distributions": len(all_metrics),
                "mean_distribution_rmse_percent": float(np.mean([row["distribution_rmse_percent"] for row in all_metrics])),
                "mean_total_variation_distance": float(np.mean([row["total_variation_distance"] for row in all_metrics])),
                "mean_wasserstein_distance_mm": float(np.mean([row["wasserstein_distance_mm"] for row in all_metrics])),
            }
        )
    return details, metrics_rows, summary_rows


def fig6_diagnostic() -> tuple[list[dict[str, float | int | str]], list[dict[str, float | int | str]]]:
    """Re-tabulate the separately fitted Fig. 6 quasi-steady diagnostic."""
    source = read_csv(FIG6_SOURCE)
    details: list[dict[str, float | int | str]] = []
    summaries: list[dict[str, float | int | str]] = []
    for preculture in (FIT_PRECULTURE, PREDICTION_PRECULTURE):
        selected = [row for row in source if as_int(row, "preculture_irradiance_klx") == preculture]
        observed = np.asarray([as_float(row, "volume_averaged_diameter_mm") for row in selected])
        predicted = np.asarray([as_float(row, "predicted_Dv_mm") for row in selected])
        for row, calculation in zip(selected, predicted):
            details.append(
                {
                    "split": "fig6_diagnostic_calibration" if preculture == FIT_PRECULTURE else "fig6_diagnostic_prediction",
                    "preculture_irradiance_klx": preculture,
                    "lighted_volume_percent": as_int(row, "lighted_volume_percent"),
                    "average_light_intensity_V": as_float(row, "average_light_intensity_V"),
                    "observed_dv_mm": as_float(row, "volume_averaged_diameter_mm"),
                    "predicted_dv_mm": float(calculation),
                    "residual_mm": float(calculation - as_float(row, "volume_averaged_diameter_mm")),
                }
            )
        summaries.append(
            {
                "split": "fig6_quasisteady_3klx" if preculture == FIT_PRECULTURE else "fig6_quasisteady_10klx",
                "preculture_irradiance_klx": preculture,
                **regression_metrics(observed, predicted),
            }
        )
    return details, summaries


def write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write to {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    path: Path,
    parameters: FitParameters,
    fig5: list[dict[str, float | int | str]],
    fig3: list[dict[str, float | int | str]],
    fig6: list[dict[str, float | int | str]],
) -> None:
    calibration, prediction = fig5
    distributions_3, distributions_10 = fig3
    fig6_3, fig6_10 = fig6
    with path.open("w") as handle:
        handle.write("# Zhang--Kojima hierarchy results\n\n")
        handle.write("The Fig. 5 fit uses only 3-klx transient observations and adjusts only ")
        handle.write("Dc_min, Dc_max, and Kc; all other baseline kinetics are fixed.\n\n")
        handle.write(
            f"Fitted values: Dc_min={parameters.dc_min:.8f} mm; "
            f"Dc_max={parameters.dc_max:.8f} mm; Kc={parameters.kc:.8f} V.\n\n"
        )
        handle.write(
            f"Fig. 5 3-klx calibration: RMSE={calibration['rmse_mm']:.5f} mm, "
            f"MAE={calibration['mae_mm']:.5f} mm, R2={calibration['r2']:.3f}.\n"
        )
        handle.write(
            f"Fig. 5 10-klx no-refit prediction: RMSE={prediction['rmse_mm']:.5f} mm, "
            f"MAE={prediction['mae_mm']:.5f} mm, R2={prediction['r2']:.3f}.\n\n"
        )
        handle.write(
            f"No-refit Fig. 3 distributions, 3 klx: mean TV={distributions_3['mean_total_variation_distance']:.3f}, "
            f"mean W1={distributions_3['mean_wasserstein_distance_mm']:.4f} mm.\n"
        )
        handle.write(
            f"No-refit Fig. 3 distributions, 10 klx: mean TV={distributions_10['mean_total_variation_distance']:.3f}, "
            f"mean W1={distributions_10['mean_wasserstein_distance_mm']:.4f} mm.\n\n"
        )
        handle.write(
            "Fig. 6 remains a separate quasi-steady light--size diagnostic, with its own fit; it is not pooled with Fig. 5.\n"
        )
        handle.write(
            f"Fig. 6 3-klx diagnostic: RMSE={fig6_3['rmse_mm']:.5f} mm, R2={fig6_3['r2']:.3f}.\n"
        )
        handle.write(
            f"Fig. 6 10-klx diagnostic: RMSE={fig6_10['rmse_mm']:.5f} mm, R2={fig6_10['r2']:.3f}.\n"
        )


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    fig5_source = read_csv(SOURCE_DATA / "zhang_kojima_1998_fig5_digitized.csv")
    table1_source = read_csv(SOURCE_DATA / "zhang_kojima_1998_table1_transcribed.csv")
    fig3_source = read_csv(SOURCE_DATA / "zhang_kojima_1998_fig3_digitized.csv")
    history = LightHistory(table1_source)
    model = DiameterPBE()

    parameters = fit_fig5(model, history, group_fig5(fig5_source))
    fig5_rows, fig5_summary = fig5_predictions(model, history, group_fig5(fig5_source), parameters)
    fig3_rows, fig3_metrics, fig3_summary = fig3_predictions(model, history, group_fig3(fig3_source), parameters)
    fig6_rows, fig6_summary = fig6_diagnostic()

    write_csv(RESULTS / "zhang_fig5_predictions.csv", fig5_rows)
    write_csv(RESULTS / "zhang_fig5_metrics.csv", fig5_summary)
    write_csv(RESULTS / "zhang_fig3_distribution_predictions.csv", fig3_rows)
    write_csv(RESULTS / "zhang_fig3_distribution_metrics.csv", fig3_metrics)
    write_csv(RESULTS / "zhang_fig3_distribution_summary.csv", fig3_summary)
    write_csv(RESULTS / "zhang_fig6_diagnostic_points.csv", fig6_rows)
    write_csv(RESULTS / "zhang_fig6_diagnostic_metrics.csv", fig6_summary)
    write_summary(RESULTS / "zhang_hierarchy_summary.md", parameters, fig5_summary, fig3_summary, fig6_summary)
    print(f"Wrote Zhang hierarchy results to {RESULTS}")


if __name__ == "__main__":
    main()
