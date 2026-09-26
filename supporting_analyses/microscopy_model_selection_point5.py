#!/usr/bin/env python3
"""Prospective microscopy-constrained biological-event model-selection test.

Uno et al. (2015) and Suzuki et al. (2013) establish structural facts about
race-B *Botryococcus braunii*: sheath/ECM renewal accompanies division,
old structural material can be released as a shell, extracellular material is
produced on a finite cell-cycle sequence, and immature daughters have a
different lipid/ECM state.  They do *not* report tracked colony-level event
frequencies, daughter/parent size-ratio distributions, or time-resolved
colony-size distributions.  It would therefore be invalid to fit a biological
rate or a daughter-kernel parameter from those microscopy papers.

This script turns that boundary into a reproducible prospective discrimination
test.  It enumerates twelve structural closures:

* three currently unresolved biological mappings (equal binary, broad binary,
  and a retained-branch/satellite representation of a hierarchical colony);
* two event-level shell alternatives (no resolved shell versus a released
  structural-material branch); and
* two new-branch reset rules (deep versus partial reset of a reduced
  lipid/ECM state).

The numerical values specify *measurement contrasts*, not estimates from the
microscopy.  Every candidate is evaluated on the same synthetic parent cohort
and the same known, two-pulse event schedule.  Consequently no beta_bio,
event hazard, growth rate, or relaxation rate is estimated or compared.  The
event schedule is conditioned on, just as it would be when events are first
located in a time-resolved movie.

The test asks which planned observations distinguish the closures: (i) a
normalised size distribution after recorded event pulses, (ii) paired
mother--descendant sizes and descendant count, (iii) a shell-material label,
and (iv) a daughter state-marker trajectory at 0, 0.5, and 1 cycle after the
event.  Pairwise standardized distances use declared design resolution only;
they are not likelihoods for the Uno or Suzuki observations.  Thus a complete
separation of candidates means that the *proposed experiment* is informative,
not that an alternative has been selected from existing microscopy.

Requirements: Python 3 with NumPy.  Run from the workspace root:

    python3 CODES/microscopy_model_selection_point5.py
    python3 CODES/plot_microscopy_model_selection.py
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
RESULTS = HERE / "results"
UNO_CONSTRAINTS = (
    PROJECT_ROOT
    / "Comparison_Literature"
    / "Comparison_Uno_Et_Al"
    / "Uno2015_comparison_package"
    / "uno_2015_model_constraints.csv"
)
SUZUKI_CONSTRAINTS = (
    PROJECT_ROOT
    / "Comparison_Literature"
    / "Comparison_Suzuki"
    / "Suzuki2013_comparison_package"
    / "suzuki_2013_model_constraints.csv"
)


# These values are protocol contrasts used only to test observability.  They
# are deliberately not read from, or calibrated to, the microscopy papers.
RANDOM_SEED = 20260907
N_TRACKED_EVENTS = 720
N_POPULATION_COLONIES = 3600
EVENT_PULSE_COUNTS = (0, 360, 720)
EVENT_PULSE_TIMES_CYCLES = (0.0, 0.25, 0.75)
PARENT_MEDIAN_DIAMETER_UM = 350.0
PARENT_LOG_STANDARD_DEVIATION = 0.22
SHELL_LOSS_FRACTION = 0.12
DEEP_RESET_STATE = 0.15
PARTIAL_RESET_STATE = 0.55
RETAINED_BRANCH_MATERIAL_FRACTION = 0.72
RETAINED_BRANCH_INITIAL_STATE = 0.85
STATE_RECOVERY_TIME_CYCLES = 0.75
STATE_SAMPLE_TIMES_CYCLES = (0.0, 0.5, 1.0)
DIAMETER_EDGES_UM = np.geomspace(25.0, 1200.0, 25)
RATIO_EDGES = np.linspace(0.05, 1.05, 17)

# Declared prospective measurement resolution.  It is used only to express a
# pairwise signal in standardised units; no biological parameter is fitted.
EFFECTIVE_DISTRIBUTION_SAMPLE_SIZE = 600
MINIMUM_DISTRIBUTION_STANDARD_ERROR = 0.007
EFFECTIVE_LINEAGE_SAMPLE_SIZE = 240
MINIMUM_LINEAGE_STANDARD_ERROR = 0.015
DESCENDANT_COUNT_STANDARD_ERROR = 0.08
SHELL_LABEL_STANDARD_ERROR = 0.025
STATE_MARKER_STANDARD_ERROR = 0.06
DISCRIMINATION_THRESHOLD = 3.0


@dataclass(frozen=True)
class Candidate:
    """One structurally admissible, deliberately non-fitted event closure."""

    kernel: str
    shell_alternative: str
    shell_loss_fraction: float
    reset_rule: str
    new_branch_initial_state: float

    @property
    def identifier(self) -> str:
        kernel_code = {
            "equal_binary": "E",
            "broad_binary": "B",
            "retained_branch_satellites": "R",
        }[self.kernel]
        shell_code = "0" if self.shell_loss_fraction == 0.0 else "S"
        reset_code = "D" if self.reset_rule == "deep_reset" else "P"
        return f"{kernel_code}{shell_code}{reset_code}"

    @property
    def kernel_label(self) -> str:
        return {
            "equal_binary": "Equal binary",
            "broad_binary": "Broad binary",
            "retained_branch_satellites": "Retained branch + satellites",
        }[self.kernel]

    @property
    def reset_label(self) -> str:
        return "Deep reset" if self.reset_rule == "deep_reset" else "Partial reset"


@dataclass(frozen=True)
class Descendant:
    """A viable branch produced in one event-conditioned lineage record."""

    event_id: int
    branch: str
    is_new_branch: bool
    parent_diameter_um: float
    material_fraction: float
    descendant_diameter_um: float
    initial_state: float


@dataclass(frozen=True)
class CandidateResult:
    """All model-predicted observables for one closure."""

    candidate: Candidate
    descendants_by_event: tuple[tuple[Descendant, ...], ...]
    distributions: dict[float, np.ndarray]
    number_distributions: dict[float, np.ndarray]
    viable_material_fraction: dict[float, float]
    lineage_vector: np.ndarray
    distribution_vector: np.ndarray
    shell_vector: np.ndarray
    state_vector: np.ndarray
    summary: dict[str, float]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def validate_source_constraints() -> list[dict[str, str]]:
    """Read and verify the source statements used to constrain the closure.

    The explicit checks make it difficult for a future edit to silently turn
    a microscopy structural result into a fictitious rate or kernel fit.
    """

    uno = read_rows(UNO_CONSTRAINTS)
    suzuki = read_rows(SUZUKI_CONSTRAINTS)
    uno_observations = {row["experimental_observation"] for row in uno}
    suzuki_observations = {row["observation"] for row in suzuki}
    required_uno = {"Colony matrix construction", "Shell release", "Daughter number and size ratio", "Event frequency"}
    required_suzuki = {"Colony hierarchy", "Secretion timing", "Persistent lipid sheets", "Lipid-body cycle"}
    if not required_uno.issubset(uno_observations):
        raise ValueError("Uno constraint table no longer contains the required structural observations")
    if not required_suzuki.issubset(suzuki_observations):
        raise ValueError("Suzuki constraint table no longer contains the required structural observations")

    selected = [
        {
            "source": "Uno et al. (2015)",
            "microscopy_constraint": "Sheath construction begins in mitosis and completes during daughter maturation.",
            "closure_consequence": "Condition biological events on a maturation/cell-cycle record; do not use an instantaneous size-only event law.",
            "identification_boundary": "No colony-level event frequency is reported; beta_bio is not fitted here.",
        },
        {
            "source": "Uno et al. (2015)",
            "microscopy_constraint": "Old retaining-wall/thin-layer material is released as shell-like structures.",
            "closure_consequence": "Compare an event with no resolved shell and an explicit structural-material-loss branch.",
            "identification_boundary": "The event-level loss fraction is not reported; the 0.12 protocol contrast is not a fitted value.",
        },
        {
            "source": "Uno et al. (2015)",
            "microscopy_constraint": "No statistical daughter number or daughter/parent size-ratio distribution is reported.",
            "closure_consequence": "Keep equal, broad, and retained-branch daughter maps as unresolved candidate kernels.",
            "identification_boundary": "No daughter-kernel parameter is fitted from microscopy.",
        },
        {
            "source": "Suzuki et al. (2013)",
            "microscopy_constraint": "Extracellular secretion follows division and persistent lipid sheets link related cells.",
            "closure_consequence": "Treat biological descendants as state-structured rather than as products of instantaneous fragmentation.",
            "identification_boundary": "The microscopic sequence does not supply a colony-level kinetic rate.",
        },
        {
            "source": "Suzuki et al. (2013)",
            "microscopy_constraint": "Immature daughters contain few small lipid bodies, which increase during maturation.",
            "closure_consequence": "Compare deep and partial reduced-state resets with a time-resolved daughter marker.",
            "identification_boundary": "The reduced state values 0.15 and 0.55 are protocol labels, not fitted microscopy measurements.",
        },
    ]
    return selected


def make_candidates() -> tuple[Candidate, ...]:
    """Return the 3 x 2 x 2 structural candidate set without a rate axis."""

    candidates: list[Candidate] = []
    for kernel in ("equal_binary", "broad_binary", "retained_branch_satellites"):
        for shell_alternative, shell_loss in (("no_resolved_shell", 0.0), ("released_shell", SHELL_LOSS_FRACTION)):
            for reset_rule, reset_state in (("deep_reset", DEEP_RESET_STATE), ("partial_reset", PARTIAL_RESET_STATE)):
                candidates.append(
                    Candidate(
                        kernel=kernel,
                        shell_alternative=shell_alternative,
                        shell_loss_fraction=shell_loss,
                        reset_rule=reset_rule,
                        new_branch_initial_state=reset_state,
                    )
                )
    return tuple(candidates)


def make_common_protocol_inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate a common event-conditioned cohort and latent kernel draws."""

    rng = np.random.default_rng(RANDOM_SEED)
    parents = np.exp(
        rng.normal(math.log(PARENT_MEDIAN_DIAMETER_UM), PARENT_LOG_STANDARD_DEVIATION, size=N_POPULATION_COLONIES)
    )
    broad_binary_fraction = rng.beta(3.0, 3.0, size=N_TRACKED_EVENTS)
    satellite_weights = rng.dirichlet((7.0, 4.0, 2.0), size=N_TRACKED_EVENTS)
    if np.any(parents <= 0.0) or np.any(broad_binary_fraction <= 0.0) or np.any(broad_binary_fraction >= 1.0):
        raise RuntimeError("invalid common protocol draw")
    return parents, broad_binary_fraction, satellite_weights


def event_descendants(
    candidate: Candidate,
    event_id: int,
    parent_diameter_um: float,
    broad_fraction: float,
    satellite_weights: np.ndarray,
) -> tuple[Descendant, ...]:
    """Apply one structural mapping to a tracked event; no event rate appears.

    Each material fraction is relative to the parent immediately before its
    observed event.  The sum of viable fractions is 1 - shell loss exactly.
    """

    viable = 1.0 - candidate.shell_loss_fraction
    branches: list[tuple[str, bool, float, float]]
    if candidate.kernel == "equal_binary":
        branches = [
            ("daughter", True, 0.5 * viable, candidate.new_branch_initial_state),
            ("daughter", True, 0.5 * viable, candidate.new_branch_initial_state),
        ]
    elif candidate.kernel == "broad_binary":
        branches = [
            ("daughter", True, broad_fraction * viable, candidate.new_branch_initial_state),
            ("daughter", True, (1.0 - broad_fraction) * viable, candidate.new_branch_initial_state),
        ]
    elif candidate.kernel == "retained_branch_satellites":
        retained = RETAINED_BRANCH_MATERIAL_FRACTION * viable
        satellite_total = (1.0 - RETAINED_BRANCH_MATERIAL_FRACTION) * viable
        branches = [("retained_branch", False, retained, RETAINED_BRANCH_INITIAL_STATE)]
        branches.extend(
            ("satellite_daughter", True, satellite_total * float(weight), candidate.new_branch_initial_state)
            for weight in satellite_weights
        )
    else:
        raise ValueError(f"unknown candidate kernel: {candidate.kernel}")

    fractions = np.asarray([branch[2] for branch in branches])
    if not np.isclose(fractions.sum(), viable, rtol=0.0, atol=2.0e-15):
        raise RuntimeError("event mapping does not obey declared shell material balance")
    if np.any(fractions <= 0.0):
        raise RuntimeError("event mapping generated a non-positive viable branch")

    return tuple(
        Descendant(
            event_id=event_id,
            branch=branch,
            is_new_branch=is_new_branch,
            parent_diameter_um=parent_diameter_um,
            material_fraction=material_fraction,
            descendant_diameter_um=parent_diameter_um * material_fraction ** (1.0 / 3.0),
            initial_state=initial_state,
        )
        for branch, is_new_branch, material_fraction, initial_state in branches
    )


def state_recovery(initial_state: np.ndarray, time_cycles: float) -> np.ndarray:
    """A common illustrative readout map, not a fitted relaxation law."""

    return 1.0 - (1.0 - initial_state) * math.exp(-time_cycles / STATE_RECOVERY_TIME_CYCLES)


def normalised_distribution(
    parents: np.ndarray,
    descendants_by_event: tuple[tuple[Descendant, ...], ...],
    event_count: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Map a known event pulse into number- and volume-normalised histograms.

    The pulse count is conditioned on rather than generated by a biological
    hazard.  This makes the distribution check a consequence of the structural
    mapping, not an implicit rate fit.
    """

    diameters: list[float] = list(parents[event_count:])
    diameters.extend(
        descendant.descendant_diameter_um
        for event in descendants_by_event[:event_count]
        for descendant in event
    )
    diameters_array = np.asarray(diameters)
    if diameters_array.min() < DIAMETER_EDGES_UM[0] or diameters_array.max() > DIAMETER_EDGES_UM[-1]:
        raise RuntimeError("protocol distribution exceeds the declared measured diameter range")
    weights = diameters_array**3
    volume_histogram, _ = np.histogram(diameters_array, bins=DIAMETER_EDGES_UM, weights=weights)
    number_histogram, _ = np.histogram(diameters_array, bins=DIAMETER_EDGES_UM)
    if volume_histogram.sum() <= 0.0 or number_histogram.sum() <= 0.0:
        raise RuntimeError("empty volume distribution")
    initial_material = float(np.sum(parents**3))
    return (
        volume_histogram / volume_histogram.sum(),
        number_histogram.astype(float) / number_histogram.sum(),
        float(weights.sum() / initial_material),
    )


def standardised_rms_difference(first: np.ndarray, second: np.ndarray, standard_error: np.ndarray | float) -> float:
    """Root-mean-square difference in declared measurement-resolution units."""

    errors = np.asarray(standard_error, dtype=float)
    if np.any(errors <= 0.0):
        raise ValueError("a prospective measurement error must be positive")
    return float(np.sqrt(np.mean(((first - second) / errors) ** 2)))


def distribution_distance(first: np.ndarray, second: np.ndarray) -> float:
    """Combine independently read number/volume histogram distances.

    The vector contains one 24-bin normalisation per recorded pulse and
    weighting. Averaging all 96 bins together would dilute a coherent signal
    in one planned histogram merely because other planned histograms agree.
    We calculate one standardized RMS signal per histogram and combine the
    independent readouts in quadrature.
    """

    n_bins = len(DIAMETER_EDGES_UM) - 1
    if len(first) != len(second) or len(first) % n_bins != 0:
        raise ValueError("unexpected time-resolved distribution vector shape")
    distances: list[float] = []
    for start in range(0, len(first), n_bins):
        first_block = first[start : start + n_bins]
        second_block = second[start : start + n_bins]
        mean_probability = 0.5 * (first_block + second_block)
        standard_error = np.maximum(
            np.sqrt(mean_probability * np.maximum(1.0 - mean_probability, 0.0) / EFFECTIVE_DISTRIBUTION_SAMPLE_SIZE),
            MINIMUM_DISTRIBUTION_STANDARD_ERROR,
        )
        distances.append(standardised_rms_difference(first_block, second_block, standard_error))
    return float(np.linalg.norm(distances))


def lineage_distance(first: np.ndarray, second: np.ndarray) -> float:
    """Compare ratio histogram and mean viable-descendant count."""

    ratio_first, count_first = first[:-1], first[-1:]
    ratio_second, count_second = second[:-1], second[-1:]
    mean_probability = 0.5 * (ratio_first + ratio_second)
    ratio_error = np.maximum(
        np.sqrt(mean_probability * np.maximum(1.0 - mean_probability, 0.0) / EFFECTIVE_LINEAGE_SAMPLE_SIZE),
        MINIMUM_LINEAGE_STANDARD_ERROR,
    )
    ratio_distance = standardised_rms_difference(ratio_first, ratio_second, ratio_error)
    count_distance = standardised_rms_difference(count_first, count_second, DESCENDANT_COUNT_STANDARD_ERROR)
    return float(math.hypot(ratio_distance, count_distance))


def candidate_result(
    candidate: Candidate,
    parents: np.ndarray,
    broad_fractions: np.ndarray,
    satellite_weights: np.ndarray,
) -> CandidateResult:
    """Calculate all planned observables for one non-fitted closure."""

    descendants_by_event = tuple(
        event_descendants(candidate, event_id, float(parents[event_id]), float(broad_fractions[event_id]), satellite_weights[event_id])
        for event_id in range(N_TRACKED_EVENTS)
    )
    distributions: dict[float, np.ndarray] = {}
    number_distributions: dict[float, np.ndarray] = {}
    viable_material_fraction: dict[float, float] = {}
    for event_count, time_cycles in zip(EVENT_PULSE_COUNTS, EVENT_PULSE_TIMES_CYCLES):
        volume_distribution, number_distribution, material_fraction = normalised_distribution(
            parents, descendants_by_event, event_count
        )
        distributions[time_cycles] = volume_distribution
        number_distributions[time_cycles] = number_distribution
        viable_material_fraction[time_cycles] = material_fraction

    descendants = tuple(descendant for event in descendants_by_event for descendant in event)
    ratios = np.asarray([descendant.descendant_diameter_um / descendant.parent_diameter_um for descendant in descendants])
    ratio_histogram, _ = np.histogram(ratios, bins=RATIO_EDGES)
    ratio_vector = ratio_histogram.astype(float) / ratio_histogram.sum()
    mean_descendant_count = float(np.mean([len(event) for event in descendants_by_event]))
    lineage_vector = np.r_[ratio_vector, mean_descendant_count]

    new_branch_states = np.asarray([descendant.initial_state for descendant in descendants if descendant.is_new_branch])
    state_vector = np.asarray(
        [float(np.mean(state_recovery(new_branch_states, time_cycles))) for time_cycles in STATE_SAMPLE_TIMES_CYCLES]
    )
    shell_vector = np.asarray([candidate.shell_loss_fraction])
    # The planned distribution readout includes both number and volume
    # normalisations after each recorded pulse.  A size histogram alone cannot
    # identify a reset state, but the two normalisations preserve information
    # about descendant number and material redistribution.
    late_time = EVENT_PULSE_TIMES_CYCLES[-1]
    distribution_vector = np.concatenate(
        [
            value
            for time_cycles in EVENT_PULSE_TIMES_CYCLES[1:]
            for value in (number_distributions[time_cycles], distributions[time_cycles])
        ]
    )
    large_branch_ratios = np.asarray(
        [max(descendant.descendant_diameter_um / descendant.parent_diameter_um for descendant in event) for event in descendants_by_event]
    )
    summary = {
        "mean_descendant_count": mean_descendant_count,
        "mean_viable_material_fraction": float(np.mean([sum(descendant.material_fraction for descendant in event) for event in descendants_by_event])),
        "daughter_ratio_q10": float(np.quantile(ratios, 0.10)),
        "daughter_ratio_q50": float(np.quantile(ratios, 0.50)),
        "daughter_ratio_q90": float(np.quantile(ratios, 0.90)),
        "largest_branch_ratio_q50": float(np.quantile(large_branch_ratios, 0.50)),
        "shell_label_fraction": candidate.shell_loss_fraction,
        "state_t0": float(state_vector[0]),
        "state_t0_5": float(state_vector[1]),
        "state_t1": float(state_vector[2]),
        "post_pulse_viable_material_fraction": viable_material_fraction[late_time],
    }
    return CandidateResult(
        candidate=candidate,
        descendants_by_event=descendants_by_event,
        distributions=distributions,
        number_distributions=number_distributions,
        viable_material_fraction=viable_material_fraction,
        lineage_vector=lineage_vector,
        distribution_vector=distribution_vector,
        shell_vector=shell_vector,
        state_vector=state_vector,
        summary=summary,
    )


def pairwise_rows(results: tuple[CandidateResult, ...]) -> list[dict[str, object]]:
    """Report modality-specific and combined separation for every candidate pair."""

    rows: list[dict[str, object]] = []
    for first, second in combinations(results, 2):
        distribution = distribution_distance(first.distribution_vector, second.distribution_vector)
        lineage = lineage_distance(first.lineage_vector, second.lineage_vector)
        shell = standardised_rms_difference(first.shell_vector, second.shell_vector, SHELL_LABEL_STANDARD_ERROR)
        state = standardised_rms_difference(first.state_vector, second.state_vector, STATE_MARKER_STANDARD_ERROR)
        combined = math.sqrt(distribution**2 + lineage**2 + shell**2 + state**2)
        rows.append(
            {
                "candidate_a": first.candidate.identifier,
                "candidate_b": second.candidate.identifier,
                "distribution_standardized_distance": distribution,
                "paired_lineage_standardized_distance": lineage,
                "shell_label_standardized_distance": shell,
                "state_trajectory_standardized_distance": state,
                "combined_standardized_distance": combined,
                "distribution_discriminable": distribution >= DISCRIMINATION_THRESHOLD,
                "paired_lineage_discriminable": lineage >= DISCRIMINATION_THRESHOLD,
                "shell_label_discriminable": shell >= DISCRIMINATION_THRESHOLD,
                "state_trajectory_discriminable": state >= DISCRIMINATION_THRESHOLD,
                "combined_discriminable": combined >= DISCRIMINATION_THRESHOLD,
            }
        )
    return rows


def write_csv(filename: str, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    with (RESULTS / filename).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(
    constraints: list[dict[str, str]],
    results: tuple[CandidateResult, ...],
    pairwise: list[dict[str, object]],
) -> None:
    """Write transparent source constraints, observables, and test results."""

    write_csv(
        "microscopy_model_selection_constraints.csv",
        ["source", "microscopy_constraint", "closure_consequence", "identification_boundary"],
        constraints,
    )
    candidate_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    lineage_rows: list[dict[str, object]] = []
    distribution_rows: list[dict[str, object]] = []
    for result in results:
        candidate = result.candidate
        candidate_rows.append(
            {
                "candidate_id": candidate.identifier,
                "kernel": candidate.kernel,
                "kernel_label": candidate.kernel_label,
                "shell_alternative": candidate.shell_alternative,
                "shell_loss_fraction_protocol": candidate.shell_loss_fraction,
                "reset_rule": candidate.reset_rule,
                "new_branch_initial_state_protocol": candidate.new_branch_initial_state,
                "rate_fitted_from_microscopy": "no",
                "interpretation": "Structural candidate only; numerical contrast values are not microscopy fits.",
            }
        )
        summary_rows.append({"candidate_id": candidate.identifier, **result.summary})
        for event in result.descendants_by_event:
            for descendant in event:
                lineage_rows.append(
                    {
                        "candidate_id": candidate.identifier,
                        "event_id": descendant.event_id,
                        "branch": descendant.branch,
                        "is_new_branch": descendant.is_new_branch,
                        "parent_diameter_um": descendant.parent_diameter_um,
                        "descendant_diameter_um": descendant.descendant_diameter_um,
                        "descendant_parent_diameter_ratio": descendant.descendant_diameter_um / descendant.parent_diameter_um,
                        "material_fraction_of_parent": descendant.material_fraction,
                        "initial_reduced_state": descendant.initial_state,
                    }
                )
        for time_cycles, distribution in result.distributions.items():
            for index, (left, right, fraction) in enumerate(zip(DIAMETER_EDGES_UM[:-1], DIAMETER_EDGES_UM[1:], distribution)):
                distribution_rows.append(
                    {
                        "candidate_id": candidate.identifier,
                        "time_after_recorded_event_pulse_cycles": time_cycles,
                        "diameter_um": math.sqrt(left * right),
                        "normalised_number_fraction": result.number_distributions[time_cycles][index],
                        "normalised_volume_fraction": fraction,
                        "viable_material_fraction_of_initial_population": result.viable_material_fraction[time_cycles],
                    }
                )
    write_csv(
        "microscopy_model_selection_candidates.csv",
        [
            "candidate_id",
            "kernel",
            "kernel_label",
            "shell_alternative",
            "shell_loss_fraction_protocol",
            "reset_rule",
            "new_branch_initial_state_protocol",
            "rate_fitted_from_microscopy",
            "interpretation",
        ],
        candidate_rows,
    )
    write_csv(
        "microscopy_model_selection_lineages.csv",
        [
            "candidate_id",
            "event_id",
            "branch",
            "is_new_branch",
            "parent_diameter_um",
            "descendant_diameter_um",
            "descendant_parent_diameter_ratio",
            "material_fraction_of_parent",
            "initial_reduced_state",
        ],
        lineage_rows,
    )
    write_csv(
        "microscopy_model_selection_lineage_summary.csv",
        [
            "candidate_id",
            "mean_descendant_count",
            "mean_viable_material_fraction",
            "daughter_ratio_q10",
            "daughter_ratio_q50",
            "daughter_ratio_q90",
            "largest_branch_ratio_q50",
            "shell_label_fraction",
            "state_t0",
            "state_t0_5",
            "state_t1",
            "post_pulse_viable_material_fraction",
        ],
        summary_rows,
    )
    write_csv(
        "microscopy_model_selection_distributions.csv",
        [
            "candidate_id",
            "time_after_recorded_event_pulse_cycles",
            "diameter_um",
            "normalised_number_fraction",
            "normalised_volume_fraction",
            "viable_material_fraction_of_initial_population",
        ],
        distribution_rows,
    )
    write_csv(
        "microscopy_model_selection_pairwise.csv",
        [
            "candidate_a",
            "candidate_b",
            "distribution_standardized_distance",
            "paired_lineage_standardized_distance",
            "shell_label_standardized_distance",
            "state_trajectory_standardized_distance",
            "combined_standardized_distance",
            "distribution_discriminable",
            "paired_lineage_discriminable",
            "shell_label_discriminable",
            "state_trajectory_discriminable",
            "combined_discriminable",
        ],
        pairwise,
    )

    total_pairs = len(pairwise)
    modality_columns = (
        ("Time-resolved distributions", "distribution_discriminable"),
        ("Paired lineage only", "paired_lineage_discriminable"),
        ("Shell label only", "shell_label_discriminable"),
        ("State trajectory only", "state_trajectory_discriminable"),
        ("All planned readouts", "combined_discriminable"),
    )
    counts = {label: sum(bool(row[column]) for row in pairwise) for label, column in modality_columns}
    with (RESULTS / "microscopy_model_selection_summary.md").open("w") as handle:
        handle.write("# Microscopy-constrained prospective model-selection test\n\n")
        handle.write("This is an event-conditioned in-silico protocol design, not a fit to Uno or Suzuki microscopy. ")
        handle.write("No biological rate, event hazard, kernel parameter, shell-loss fraction, or reset time was estimated.\n\n")
        handle.write("## Structural candidate set\n\n")
        handle.write("Twelve closures combine three unresolved biological kernels, two shell alternatives, and two daughter reset rules. ")
        handle.write("All use the same 720 synthetic tracked events and known two-pulse schedule, so candidate differences are structural mappings rather than rate effects.\n\n")
        handle.write("## Declared prospective discrimination criterion\n\n")
        handle.write(f"A pair is called distinguishable only when its standardized signal is at least {DISCRIMINATION_THRESHOLD:.1f} under the declared measurement resolution. ")
        handle.write("The result is a protocol-adequacy statement, not an empirical model ranking.\n\n")
        for label, _ in modality_columns:
            handle.write(f"- {label}: {counts[label]}/{total_pairs} candidate pairs.\n")
        handle.write("\nTime-resolved size distributions cannot see a reset rule that leaves the immediate material map unchanged. ")
        handle.write("Paired mother--descendant sizes resolve kernel topology, a shell label resolves structural loss, and a time-resolved daughter marker resolves reset. ")
        handle.write("Existing microscopy supplies none of these statistical event datasets, so it does not select one closure.\n")


def main() -> None:
    constraints = validate_source_constraints()
    candidates = make_candidates()
    parents, broad_fractions, satellite_weights = make_common_protocol_inputs()
    results = tuple(
        candidate_result(candidate, parents, broad_fractions, satellite_weights)
        for candidate in candidates
    )
    pairwise = pairwise_rows(results)
    if len(pairwise) != math.comb(len(candidates), 2):
        raise RuntimeError("incomplete candidate-pair report")
    write_outputs(constraints, results, pairwise)
    total_pairs = len(pairwise)
    combined_pairs = sum(bool(row["combined_discriminable"]) for row in pairwise)
    print("Microscopy-constrained prospective model-selection test completed.")
    print(f"Wrote 12 structural candidates and {total_pairs} pairwise comparisons to {RESULTS}.")
    print(f"All planned readouts distinguish {combined_pairs}/{total_pairs} pairs at the declared protocol resolution.")
    print("No Uno/Suzuki rate, event frequency, daughter-kernel parameter, or shell-loss fraction was fitted.")


if __name__ == "__main__":
    main()
