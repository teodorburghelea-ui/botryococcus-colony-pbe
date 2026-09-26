#!/usr/bin/env python3
"""Run the van den Berg low-light memory and bimodality model-selection test.

Both calculations begin with exactly the same reconstructed medium-light
(``Start``) *volume* distribution of van den Berg et al.  That distribution is
converted to a number distribution only for PBE initialization.  The day-20
low-light curves are never supplied to either time integrator: they are read
only after the calculation to report diagnostic distances and topology.

The one-state calculation is the state-collapsed limit of the test closure. It
contains the common, weak equal-volume mechanical-breakup operator, but no
finite-time acclimation state and hence no transition-gated biological
daughter-production pulse.  The structured calculation adds (i) relaxation of
an initially medium-light state after a low-light step and (ii) one
material-conserving biological event kernel.  That kernel retains a large
mother remnant and places the transferred material in a fixed, broad small
daughter-colony kernel.  It is a daughter kernel, not a fitted mixture of the
day-20 distribution.

Continuous material transport, shell loss, and growth calibration are set to
zero in both variants.  The test deliberately isolates the topology generated
by finite-time state and a biological daughter channel; it is not a biomass or
mass-growth fit.

There is deliberately no optimizer in this file.  The fixed closure is a
structural feasibility/model-selection test, not parameter identification.
The final sensitivity rows vary only the acclimation time and material fraction
of the same biological event; all use the same initial distribution, daughter
kernel, and low-light forcing.

Requirements: Python 3 with NumPy and SciPy.  Run from the workspace root:

    python3 CODES/memory_bimodality_point3.py

Then create the shared-style manuscript figure with:

    python3 CODES/plot_memory_bimodality.py
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
RESULTS = HERE / "results"
SOURCE_DATA = (
    PROJECT_ROOT
    / "Comparison_Literature"
    / "Comparison_Van_Der_Berg"
    / "VanDenBerg2019_comparison_package"
    / "van_den_berg_2019_fig2_approx_digitized.csv"
)


@dataclass(frozen=True)
class TestParameters:
    """Fixed, illustrative closure values for the structural test.

    Units are mm and d.  The low-light experiment is used only to test the
    topology at 20 d, so these values are intentionally *not* inferred from
    its day-20 distribution.  ``daughter_ratio`` is applied to the observed
    medium-light modal diameter, not to individual parents; it therefore
    represents a biologically selected daughter scale instead of a fitted
    initial mixture component.  Continuous material transport and shell loss
    are intentionally zero in this topology-only calculation.
    """

    n_diameter: int = 180
    n_state: int = 64
    d_min_mm: float = 0.005
    d_max_mm: float = 2.5
    initial_state: float = 0.94
    initial_log_width: float = 0.060
    end_time_d: float = 20.0
    time_step_d: float = 0.025
    mechanical_rate_d_inv: float = 0.020
    size_gate_center_mm: float = 0.250
    size_gate_width_mm: float = 0.035
    state_gate_low: float = 0.22
    state_gate_high: float = 0.72
    state_gate_width: float = 0.035
    biological_rate_d_inv: float = 0.80
    acclimation_time_d: float = 3.2
    bud_material_fraction: float = 0.10
    daughter_ratio: float = 0.18
    daughter_log_width: float = 0.28
    daughter_reset_state: float = 0.04


@dataclass(frozen=True)
class SourceCurves:
    """The source bin centers and volume distributions in percent."""

    diameter_mm: np.ndarray
    start_volume_percent: np.ndarray
    low_light_volume_percent: np.ndarray


@dataclass(frozen=True)
class Simulation:
    """Snapshots and conservation diagnostics from one forward calculation."""

    snapshots: dict[float, np.ndarray]
    material_relative_drift: float
    min_population: float


SMALL_WINDOW_MM = (0.045, 0.180)
LARGE_WINDOW_MM = (0.250, 0.850)
SMALL_FRACTION_CUTOFF_MM = 0.200
LARGE_RESIDUAL_CUTOFF_MM = 0.300
MIN_MODAL_FRACTION = 0.10
SNAPSHOT_TIMES = (0.0, 2.0, 4.0, 6.0, 10.0, 20.0)
SENSITIVITY_TIMES = (2.4, 3.2, 4.0)
SENSITIVITY_FRACTIONS = (0.06, 0.10)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_source_curves(path: Path) -> SourceCurves:
    """Read the common start and mean low-light reconstruction.

    The two low-light biological replicates are averaged only for the final
    diagnostic target.  They are not used to initialise or advance a model.
    """

    rows = read_csv(path)

    def curve(condition: str) -> tuple[np.ndarray, np.ndarray]:
        selected = [row for row in rows if row["condition"] == condition]
        selected.sort(key=lambda row: float(row["diameter_um"]))
        diameter = np.asarray([float(row["diameter_um"]) / 1000.0 for row in selected])
        volume = np.asarray([float(row["volume_percent"]) for row in selected])
        return diameter, volume

    diameter, start = curve("Start")
    diameter_ll1, low_light_1 = curve("LL1")
    diameter_ll2, low_light_2 = curve("LL2")
    if not (np.allclose(diameter, diameter_ll1) and np.allclose(diameter, diameter_ll2)):
        raise ValueError("the reconstructed van den Berg curves do not share bin centers")
    return SourceCurves(diameter, start, 0.5 * (low_light_1 + low_light_2))


def logistic(value: np.ndarray) -> np.ndarray:
    """Numerically stable logistic function for the smooth rate gates."""

    positive = value >= 0.0
    result = np.empty_like(value, dtype=float)
    result[positive] = 1.0 / (1.0 + np.exp(-value[positive]))
    exponential = np.exp(value[~positive])
    result[~positive] = exponential / (1.0 + exponential)
    return result


def sectional_map(
    diameter: np.ndarray,
    target_diameter: np.ndarray,
    expected_number: np.ndarray,
) -> csr_matrix:
    """Map a source section to target sections with exact D-cubed moment.

    Each map column represents one source diameter.  The output weights sum to
    ``expected_number`` and preserve the material moment represented by
    ``expected_number * target_diameter**3``.  Every target used by the test is
    checked to be inside the represented domain; no material is silently lost
    at a numerical boundary.
    """

    n_sections = len(diameter)
    if np.any(target_diameter < diameter[0]) or np.any(target_diameter > diameter[-1]):
        raise ValueError("daughter target lies outside the material grid")

    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for source, (target, count) in enumerate(zip(target_diameter, expected_number)):
        if target <= diameter[0]:
            left = right = 0
            left_weight, right_weight = count, 0.0
        elif target >= diameter[-1]:
            left = right = n_sections - 1
            left_weight, right_weight = count, 0.0
        else:
            left = int(np.searchsorted(diameter, target, side="right") - 1)
            right = left + 1
            right_weight = count * (target**3 - diameter[left] ** 3) / (
                diameter[right] ** 3 - diameter[left] ** 3
            )
            left_weight = count - right_weight
        rows.extend((left, right))
        columns.extend((source, source))
        values.extend((left_weight, right_weight))
    return csr_matrix((values, (rows, columns)), shape=(n_sections, n_sections))


def state_transport_map(state: np.ndarray, relaxation_time: float, dt: float) -> csr_matrix:
    """Exact-in-time relaxation map for da/dt=-a/tau, deposited positively."""

    target = state * math.exp(-dt / relaxation_time)
    n_state = len(state)
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for source, value in enumerate(target):
        if value <= state[0]:
            left = right = 0
            left_weight, right_weight = 1.0, 0.0
        elif value >= state[-1]:
            left = right = n_state - 1
            left_weight, right_weight = 1.0, 0.0
        else:
            left = int(np.searchsorted(state, value) - 1)
            right = left + 1
            right_weight = (value - state[left]) / (state[right] - state[left])
            left_weight = 1.0 - right_weight
        rows.extend((left, right))
        columns.extend((source, source))
        values.extend((left_weight, right_weight))
    return csr_matrix((values, (rows, columns)), shape=(n_state, n_state))


class MemoryBimodalityTest:
    """A shared initial condition and operators for the nested comparison."""

    def __init__(self, source: SourceCurves, parameters: TestParameters) -> None:
        self.source = source
        self.p = parameters
        edges = np.exp(
            np.linspace(math.log(parameters.d_min_mm), math.log(parameters.d_max_mm), parameters.n_diameter + 1)
        )
        self.diameter = np.sqrt(edges[:-1] * edges[1:])
        self.material = self.diameter**3
        self.state = (np.arange(parameters.n_state, dtype=float) + 0.5) / parameters.n_state
        self.initial_structured = self._initial_distribution()
        equal_volume_daughter = self.diameter / 2.0 ** (1.0 / 3.0)
        # A daughter of the lowest represented section would lie outside the
        # grid.  That section has an exactly zero mechanical rate below, so
        # the harmless clipped map is never applied there; this avoids an
        # unresolved boundary event masquerading as material conservation.
        self.hydrodynamic_map = sectional_map(
            self.diameter,
            np.maximum(equal_volume_daughter, self.diameter[0]),
            np.full(parameters.n_diameter, 2.0),
        )
        self.size_gate = logistic((self.diameter - parameters.size_gate_center_mm) / parameters.size_gate_width_mm)
        self.hydrodynamic_rate = parameters.mechanical_rate_d_inv * self.size_gate
        self.hydrodynamic_rate[equal_volume_daughter < self.diameter[0]] = 0.0

    def _initial_distribution(self) -> np.ndarray:
        """Map the sole medium-light start curve to one high-state population."""

        number_weight = self.source.start_volume_percent / self.source.diameter_mm**3
        number_weight /= number_weight.sum()
        population = np.zeros((self.p.n_diameter, self.p.n_state))
        initial_state_index = int(np.argmin(np.abs(self.state - self.p.initial_state)))
        for diameter, weight in zip(self.source.diameter_mm, number_weight):
            kernel = np.exp(-0.5 * ((np.log(self.diameter) - math.log(diameter)) / self.p.initial_log_width) ** 2)
            kernel /= kernel.sum()
            population[:, initial_state_index] += weight * kernel
        return population

    def biological_operators(
        self,
        relaxation_time: float,
        material_fraction: float,
    ) -> tuple[csr_matrix, np.ndarray, csr_matrix, np.ndarray, int]:
        """Build the fixed state-gated mother-plus-daughter biological closure."""

        if not 0.0 < material_fraction < 1.0:
            raise ValueError("the biological material fraction must lie between zero and one")
        mother_target = self.diameter * (1.0 - material_fraction) ** (1.0 / 3.0)
        mother_map = sectional_map(
            self.diameter,
            np.maximum(mother_target, self.diameter[0]),
            np.ones_like(self.diameter),
        )

        # A single broad daughter *material* kernel, centered at the reported
        # five-to-six-fold smaller scale.  Its amplitude is supplied by the
        # conserved material fraction above, not by fitting target mixture
        # weights.  Dividing by D^3 converts material weights to expected
        # daughter counts.
        medium_light_mode = self.source.diameter_mm[int(np.argmax(self.source.start_volume_percent))]
        daughter_diameter = self.p.daughter_ratio * medium_light_mode
        daughter_material_kernel = np.exp(
            -0.5 * ((np.log(self.diameter) - math.log(daughter_diameter)) / self.p.daughter_log_width) ** 2
        )
        daughter_material_kernel /= daughter_material_kernel.sum()
        daughter_number_per_material = material_fraction * daughter_material_kernel / self.material
        if not np.isclose(np.dot(self.material, daughter_number_per_material), material_fraction):
            raise RuntimeError("biological daughter kernel does not conserve its assigned material fraction")

        state_gate = logistic((self.state - self.p.state_gate_low) / self.p.state_gate_width) * logistic(
            (self.p.state_gate_high - self.state) / self.p.state_gate_width
        )
        biological_rate = self.p.biological_rate_d_inv * self.size_gate[:, None] * state_gate[None, :]
        biological_rate[mother_target < self.diameter[0], :] = 0.0
        state_map = state_transport_map(self.state, relaxation_time, self.p.time_step_d)
        daughter_state_index = int(np.argmin(np.abs(self.state - self.p.daughter_reset_state)))
        return mother_map, daughter_number_per_material, state_map, biological_rate, daughter_state_index

    def simulate_one_state(self) -> Simulation:
        """Advance the state-collapsed one-coordinate limit for 20 d."""

        population = self.initial_structured.sum(axis=1).copy()
        initial_material = float(np.dot(self.material, population))
        snapshots: dict[float, np.ndarray] = {}
        snapshot_steps = {int(round(time / self.p.time_step_d)): time for time in SNAPSHOT_TIMES}
        n_steps = int(round(self.p.end_time_d / self.p.time_step_d))
        minimum = float(population.min())
        for step in range(n_steps + 1):
            if step in snapshot_steps:
                snapshots[snapshot_steps[step]] = population.copy()
            if step == n_steps:
                break
            loss = self.p.time_step_d * self.hydrodynamic_rate * population
            population = population - loss + self.hydrodynamic_map @ loss
            minimum = min(minimum, float(population.min()))
        final_material = float(np.dot(self.material, population))
        return Simulation(snapshots, (final_material - initial_material) / initial_material, minimum)

    def simulate_structured(self, relaxation_time: float, material_fraction: float) -> Simulation:
        """Advance the finite-time-state, biological daughter-production model."""

        mother_map, daughter_per_material, state_map, biological_rate, daughter_state = self.biological_operators(
            relaxation_time, material_fraction
        )
        population = self.initial_structured.copy()
        initial_material = float(np.dot(self.material, population.sum(axis=1)))
        snapshots: dict[float, np.ndarray] = {}
        snapshot_steps = {int(round(time / self.p.time_step_d)): time for time in SNAPSHOT_TIMES}
        n_steps = int(round(self.p.end_time_d / self.p.time_step_d))
        minimum = float(population.min())
        for step in range(n_steps + 1):
            if step in snapshot_steps:
                snapshots[snapshot_steps[step]] = population.copy()
            if step == n_steps:
                break

            # Low-light step: an initially medium-light state moves through a
            # finite transition window.  The map is number- and material-
            # preserving because it acts only in a.
            population = np.asarray((state_map @ population.T).T)

            # The common weak mechanical operator is deliberately identical
            # to that in the one-state limit.
            hydrodynamic_loss = self.p.time_step_d * self.hydrodynamic_rate[:, None] * population
            population = population - hydrodynamic_loss + self.hydrodynamic_map @ hydrodynamic_loss

            # A biological event retains one large mother branch and creates
            # a reset small-daughter branch.  The mother map retains
            # (1-f) of every affected parent's material; the daughter kernel
            # receives exactly f, so the full event is material conserving.
            biological_loss = self.p.time_step_d * biological_rate * population
            population = population - biological_loss + mother_map @ biological_loss
            daughter_material = float(np.dot(self.material, biological_loss.sum(axis=1)))
            population[:, daughter_state] += daughter_material * daughter_per_material
            minimum = min(minimum, float(population.min()))

        final_material = float(np.dot(self.material, population.sum(axis=1)))
        return Simulation(snapshots, (final_material - initial_material) / initial_material, minimum)

    def volume_curve(self, population: np.ndarray) -> np.ndarray:
        """Transform a PBE number population to the source volume bins."""

        marginal = population.sum(axis=1) if population.ndim == 2 else population
        material_density = self.material * marginal
        log_centers = np.log(self.source.diameter_mm)
        bin_edges = np.concatenate(
            (
                [log_centers[0] - 0.5 * (log_centers[1] - log_centers[0])],
                0.5 * (log_centers[:-1] + log_centers[1:]),
                [log_centers[-1] + 0.5 * (log_centers[-1] - log_centers[-2])],
            )
        )
        log_grid = np.log(self.diameter)
        values = np.asarray(
            [material_density[(log_grid >= lower) & (log_grid < upper)].sum() for lower, upper in zip(bin_edges[:-1], bin_edges[1:])]
        )
        if values.sum() <= 0.0:
            raise RuntimeError("cannot normalize an empty predicted volume distribution")
        return 100.0 * values / values.sum()


def smooth_for_mode_detection(values: np.ndarray) -> np.ndarray:
    """Apply a fixed light log-bin smoother only to identify broad modes."""

    kernel = np.asarray([1.0, 4.0, 6.0, 4.0, 1.0]) / 16.0
    result = values.copy()
    for _ in range(2):
        result = np.convolve(np.pad(result, (2, 2), mode="edge"), kernel, mode="valid")
    return result


def mode_summary(diameter: np.ndarray, values: np.ndarray) -> dict[str, float | int | bool]:
    """Return declared topology criteria without fitting any distribution."""

    smoothed = smooth_for_mode_detection(values)
    local_maxima = np.flatnonzero((smoothed[1:-1] > smoothed[:-2]) & (smoothed[1:-1] >= smoothed[2:])) + 1
    modal_floor = 0.10 * float(smoothed.max())
    local_maxima = np.asarray([index for index in local_maxima if smoothed[index] >= modal_floor], dtype=int)

    def peak_in(window: tuple[float, float]) -> float:
        candidates = [index for index in local_maxima if window[0] <= diameter[index] <= window[1]]
        if not candidates:
            return float("nan")
        return float(diameter[max(candidates, key=lambda index: smoothed[index])])

    small_peak = peak_in(SMALL_WINDOW_MM)
    large_peak = peak_in(LARGE_WINDOW_MM)
    small_fraction = float(values[diameter < SMALL_FRACTION_CUTOFF_MM].sum() / 100.0)
    large_fraction = float(values[diameter >= LARGE_RESIDUAL_CUTOFF_MM].sum() / 100.0)
    separated = bool(
        np.isfinite(small_peak)
        and np.isfinite(large_peak)
        and large_peak / small_peak >= 2.5
    )
    topology_pass = bool(
        separated and small_fraction >= MIN_MODAL_FRACTION and large_fraction >= MIN_MODAL_FRACTION
    )
    return {
        "n_prominent_modes": int(len(local_maxima)),
        "small_mode_peak_um": 1000.0 * small_peak,
        "large_mode_peak_um": 1000.0 * large_peak,
        "small_volume_fraction": small_fraction,
        "large_residual_fraction": large_fraction,
        "topology_pass": topology_pass,
    }


def distribution_metrics(observed: np.ndarray, predicted: np.ndarray, diameter_mm: np.ndarray) -> dict[str, float]:
    """Distribution diagnostics; they are reported, never optimized here."""

    difference = predicted - observed
    return {
        "distribution_rmse_volume_percent": float(np.sqrt(np.mean(difference**2))),
        "total_variation_distance": float(0.5 * np.sum(np.abs(difference)) / 100.0),
        "wasserstein_distance_mm": float(
            np.sum(np.abs(np.cumsum(difference / 100.0)[:-1]) * np.diff(diameter_mm))
        ),
    }


def write_csv(path: Path, rows: list[dict[str, float | int | str | bool]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write to {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def distribution_rows(
    model: str,
    days: float,
    diameter_mm: np.ndarray,
    volume_percent: np.ndarray,
) -> list[dict[str, float | str]]:
    return [
        {
            "model": model,
            "days": days,
            "diameter_um": 1000.0 * float(diameter),
            "volume_percent": float(value),
        }
        for diameter, value in zip(diameter_mm, volume_percent)
    ]


def summary_row(
    model: str,
    predicted: np.ndarray,
    target: np.ndarray,
    diameter: np.ndarray,
    simulation: Simulation,
) -> dict[str, float | int | str | bool]:
    return {
        "model": model,
        **distribution_metrics(target, predicted, diameter),
        **mode_summary(diameter, predicted),
        "material_relative_drift": simulation.material_relative_drift,
        "minimum_cell_population": simulation.min_population,
    }


def write_summary(
    path: Path,
    source: SourceCurves,
    parameters: TestParameters,
    target_topology: dict[str, float | int | bool],
    one_state: dict[str, float | int | str | bool],
    structured: dict[str, float | int | str | bool],
    sensitivity: list[dict[str, float | int | str | bool]],
) -> None:
    pass_count = sum(bool(row["topology_pass"]) for row in sensitivity)
    medium_light_mode = source.diameter_mm[int(np.argmax(source.start_volume_percent))]
    daughter_scale = parameters.daughter_ratio * medium_light_mode
    with path.open("w") as handle:
        handle.write("# Low-light memory and bimodality model-selection test\n\n")
        handle.write(
            "Both PBE variants start from the same reconstructed medium-light volume distribution, "
            "converted once to a number distribution. The mean of the two day-20 low-light curves "
            "is read only after forward simulation for diagnostics; this script contains no optimizer "
            "and fits no mixture weights. Continuous material transport and shell loss are zero in "
            "both variants, so this is a topology/model-selection test rather than a biomass-growth fit.\n\n"
        )
        handle.write("## Fixed reference closure\n\n")
        handle.write(
            f"The state-collapsed one-coordinate limit has a common weak equal-volume mechanical rate "
            f"of {parameters.mechanical_rate_d_inv:.3f} d^-1. The structured calculation uses "
            f"tau_a={parameters.acclimation_time_d:.2f} d, a transition gate from "
            f"a={parameters.state_gate_low:.2f} to {parameters.state_gate_high:.2f}, and a "
            f"material-conserving mother-plus-daughter event. Each event transfers "
            f"f={parameters.bud_material_fraction:.2f} of affected parent material to one broad "
            f"daughter material kernel centered at {1000.0 * daughter_scale:.1f} um "
            f"({parameters.daughter_ratio:.2f} times the {1000.0 * medium_light_mode:.1f}-um "
            "medium-light modal diameter), while retaining the mother branch.\n\n"
        )
        handle.write("## Declared topology criterion\n\n")
        handle.write(
            "A pass requires a broad local maximum in 45--180 um, another in 250--850 um, "
            "a diameter ratio of at least 2.5, and at least 10% of volume in both the <200-um "
            "small-colony and >=300-um large-residual regions. The low-light reconstruction has "
            f"small and large mode locations of {target_topology['small_mode_peak_um']:.1f} and "
            f"{target_topology['large_mode_peak_um']:.1f} um, with fractions "
            f"{target_topology['small_volume_fraction']:.3f} and "
            f"{target_topology['large_residual_fraction']:.3f}.\n\n"
        )
        for row in (one_state, structured):
            handle.write(
                f"{row['model']}: TV={row['total_variation_distance']:.3f}; "
                f"W1={row['wasserstein_distance_mm']:.3f} mm; small fraction="
                f"{row['small_volume_fraction']:.3f}; large residual="
                f"{row['large_residual_fraction']:.3f}; small-mode peak="
                f"{row['small_mode_peak_um']:.1f} um; large-mode peak="
                f"{row['large_mode_peak_um']:.1f} um; topology pass="
                f"{row['topology_pass']}; material drift={row['material_relative_drift']:.2e}.\n"
            )
        handle.write("\n")
        handle.write(
            f"Fixed sensitivity set: {pass_count}/{len(sensitivity)} structured calculations pass the "
            "same topology criterion for tau_a in {2.4, 3.2, 4.0} d and transferred material fraction "
            "in {0.06, 0.10}. These are closure sensitivities, not refits.\n\n"
        )
        handle.write("## Interpretation\n\n")
        handle.write(
            "The state-collapsed one-coordinate limit retains a large population but has no mechanism "
            "to create a separate small-colony mode after the low-light shift. The structured closure "
            "creates a transient small daughter branch while retaining a large mother branch, and does "
            "so without inserting a day-20 mixture into the initial condition. This discriminates the "
            "model structure; it does not identify the acclimation time, daughter number, shell loss, "
            "or biological event rate. The later reported collapse to a low-light unimodal distribution "
            "is outside this 20-d test and remains an independent constraint.\n"
        )


def verify_simulation(simulation: Simulation, label: str) -> None:
    """Fail loudly if the declared conservative forward test is not conservative."""

    if abs(simulation.material_relative_drift) > 5.0e-11:
        raise RuntimeError(f"{label} has unacceptable material drift: {simulation.material_relative_drift:.3e}")
    if simulation.min_population < -5.0e-13:
        raise RuntimeError(f"{label} has a negative population: {simulation.min_population:.3e}")


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    parameters = TestParameters()
    source = read_source_curves(SOURCE_DATA)
    test = MemoryBimodalityTest(source, parameters)

    one_simulation = test.simulate_one_state()
    structured_simulation = test.simulate_structured(parameters.acclimation_time_d, parameters.bud_material_fraction)
    verify_simulation(one_simulation, "one-state calculation")
    verify_simulation(structured_simulation, "structured reference calculation")
    one_final = test.volume_curve(one_simulation.snapshots[parameters.end_time_d])
    structured_final = test.volume_curve(structured_simulation.snapshots[parameters.end_time_d])
    target_topology = mode_summary(source.diameter_mm, source.low_light_volume_percent)
    one_metrics = summary_row("one_state_collapsed_limit", one_final, source.low_light_volume_percent, source.diameter_mm, one_simulation)
    structured_metrics = summary_row(
        "structured_memory_biological_kernel",
        structured_final,
        source.low_light_volume_percent,
        source.diameter_mm,
        structured_simulation,
    )

    distributions: list[dict[str, float | str]] = []
    distributions.extend(distribution_rows("medium_light_initial", 0.0, source.diameter_mm, source.start_volume_percent))
    distributions.extend(
        distribution_rows("low_light_reconstruction_mean", parameters.end_time_d, source.diameter_mm, source.low_light_volume_percent)
    )
    timecourse: list[dict[str, float | int | str | bool]] = []
    for model, simulation in (
        ("one_state_collapsed_limit", one_simulation),
        ("structured_memory_biological_kernel", structured_simulation),
    ):
        for time in SNAPSHOT_TIMES:
            distribution = test.volume_curve(simulation.snapshots[time])
            distributions.extend(distribution_rows(model, time, source.diameter_mm, distribution))
            timecourse.append({"model": model, "days": time, **mode_summary(source.diameter_mm, distribution)})

    sensitivity: list[dict[str, float | int | str | bool]] = []
    for relaxation_time in SENSITIVITY_TIMES:
        for material_fraction in SENSITIVITY_FRACTIONS:
            simulation = test.simulate_structured(relaxation_time, material_fraction)
            verify_simulation(simulation, f"structured sensitivity tau={relaxation_time}, f={material_fraction}")
            prediction = test.volume_curve(simulation.snapshots[parameters.end_time_d])
            sensitivity.append(
                {
                    "model": "structured_sensitivity",
                    "acclimation_time_d": relaxation_time,
                    "bud_material_fraction": material_fraction,
                    **distribution_metrics(source.low_light_volume_percent, prediction, source.diameter_mm),
                    **mode_summary(source.diameter_mm, prediction),
                    "material_relative_drift": simulation.material_relative_drift,
                }
            )

    write_csv(RESULTS / "memory_bimodality_distributions.csv", distributions)
    write_csv(RESULTS / "memory_bimodality_timecourse.csv", timecourse)
    write_csv(RESULTS / "memory_bimodality_metrics.csv", [one_metrics, structured_metrics])
    write_csv(RESULTS / "memory_bimodality_sensitivity.csv", sensitivity)
    write_summary(
        RESULTS / "memory_bimodality_summary.md",
        source,
        parameters,
        target_topology,
        one_metrics,
        structured_metrics,
        sensitivity,
    )
    print(f"Wrote memory/bimodality results to {RESULTS}")
    print(
        "Topology pass: "
        f"one-state={one_metrics['topology_pass']}; structured={structured_metrics['topology_pass']}; "
        f"structured sensitivity={sum(bool(row['topology_pass']) for row in sensitivity)}/{len(sensitivity)}"
    )


if __name__ == "__main__":
    main()
