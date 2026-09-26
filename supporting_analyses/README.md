# Reproducible numerical-integrity and Zhang--Kojima checks

`structured_pbe_integrity.jl` is a dependency-free Julia implementation of the
finite-volume/sectional form of Eq. (3) in the manuscript.  It uses positive
first-order upwinding for transport in material and internal-state space,
conservative sectional deposition for daughter formation and aggregation, and
an explicit Courant--Friedrichs--Lewy (CFL) restriction for non-negative cell
populations.

Run it from this directory with:

```bash
julia structured_pbe_integrity.jl
```

It creates `results/` containing:

- `moment_balance.csv` — independent first-material-moment checks for each
  conservative process, shell loss, chemostat exchange, and the complete
  Eq. (3);
- `diameter_transform_convergence.csv` — refinement of the material-to-
  diameter transformation; and
- `daughter_kernel_convergence.csv` — refinement for the biological and
  hydrodynamic daughter kernels, each compared with an `Nx=512` reference.

The numerical values are fixed illustrative values chosen solely to exercise
the discretisation.  They are not fitted to a biological experiment.

`plot_numerical_integrity.py`, stored beside the numerical data and solver,
uses Matplotlib to convert the convergence CSV files into the vector figure
`fig6_numerical_integrity_convergence.pdf` used in Appendix B. Run it from the
workspace root after the Julia verification:

```bash
python3 CODES/plot_numerical_integrity.py
```

The script also writes a high-resolution PNG preview.

## Zhang--Kojima calibration/prediction hierarchy

`zhang_hierarchy_point2.py` reproduces the additional Zhang--Kojima tests
used in the manuscript.  It fits the deliberately minimal, one-coordinate
growth--equal-volume-separation submodel only to the 3-klx Fig.~5 transient
time series.  The only fitted quantities are $D_{c,\min}$, $D_{c,\max}$, and
$K_c$; continuous-growth and breakup-rate prefactors and the selectivity
sharpness are fixed before fitting.  The 10-klx Fig.~5 trajectories and every
non-initial Fig.~3 size histogram are then evaluated with those parameters
frozen.  The script reports scalar errors, total-variation distances, and
one-dimensional Wasserstein distances in `results/`.

The Fig.~6 analysis is intentionally not part of that transient calibration.
Its pre-existing, separately fitted quasi-steady calculation is re-tabulated
only as a light--size diagnostic, with its own calibration/prediction metrics.

From the workspace root, run:

```bash
python3 CODES/zhang_hierarchy_point2.py
python3 CODES/plot_zhang_hierarchy.py
```

These scripts require Python 3 with NumPy, SciPy, and Matplotlib.  They write
the result tables under `CODES/results/` and the vector manuscript figures
`fig7_zhang_hierarchy_validation.pdf` and
`fig8_zhang_fig6_quasisteady_diagnostic.pdf` under
`Draft_2026-09-07/figures/` (with PNG previews alongside them).

## Low-light memory and bimodality model-selection test

`memory_bimodality_point3.py` is a constrained forward test using the
approximately reconstructed van den Berg et al. medium-light start and
day-20 low-light volume distributions. Both models start from the same
medium-light distribution. The state-collapsed one-coordinate limit contains
the shared weak mechanical equal-volume operator; the structured calculation
adds a finite-time acclimation state and one material-conserving,
mother-retaining biological daughter kernel. The day-20 low-light curve is
read only after simulation to calculate distances and a declared topology
criterion. Continuous material transport and shell loss are zero in both
variants to isolate the topology; this is not a biomass-growth fit. There is
no optimiser, and no initial or final mixture is fitted.

Run it from the workspace root with:

```bash
python3 CODES/memory_bimodality_point3.py
python3 CODES/plot_memory_bimodality.py
```

The analysis writes these reproducible products under `CODES/results/`:

- `memory_bimodality_distributions.csv` — source, one-state, and structured
  source-bin volume distributions at the retained time points;
- `memory_bimodality_timecourse.csv` — small-mode and large-residual fractions
  through the 20-d low-light shift;
- `memory_bimodality_metrics.csv` — day-20 distances, mode topology, and
  material-conservation diagnostics for the nested comparison;
- `memory_bimodality_sensitivity.csv` — a fixed six-case closure sensitivity
  set, not a refit; and
- `memory_bimodality_summary.md` — the protocol, fixed values, and the
  interpretation boundary.

`plot_memory_bimodality.py` creates the common-style vector figure
`fig9_memory_bimodality_model_selection.pdf` and a PNG preview under
`Draft_2026-09-07/figures/`. It uses a fixed, uniform light smoothing only for
the visual display of source-bin curves; raw binned values remain in the CSV.

## García-Cubero environmental and chemostat test

`garcia_cubero_chemostat_point4.py` implements the deliberately minimal
steady-state environmental-state test used for García-Cubero et al. It fits a
four-coefficient log-size-selection state only to a declared nine-observation
Table-2 subset: the three centre replicates and paired contrasts for
temperature, light, and dilution. The four published optimization-validation
conditions from Tables 3--4, transcribed in
`data/garcia_cubero_2021_optimized_size_validation.csv`, are never used for
fitting. The nested null refits the same closure without a dilution-mediated
state route.

The script also independently applies only the homogeneous sectional loss
`n_i -> exp(-D t)n_i`; it verifies that normalized volume distributions and
their volume-weighted diameter are invariant. Thus any nonzero dilution
response in the state closure is explicitly routed through state-dependent
biology or retention, not uniform washout. This test does not identify a
temperature law, relaxation time, or retention mechanism.

From the workspace root, run:

```bash
python3 CODES/garcia_cubero_chemostat_point4.py
python3 CODES/plot_garcia_chemostat.py
```

The scripts write the fitted parameters, all calibration/check/prediction
records, metrics, dilution contrast, uniform-washout control, and a protocol
summary under `CODES/results/`. The common-style vector figure
`fig10_garcia_chemostat_state_test.pdf` and PNG preview are written under
`Draft_2026-09-07/figures/`.

## Microscopy-constrained biological-event model selection

`microscopy_model_selection_point5.py` uses the Uno and Suzuki constraint
tables as structural evidence only. It validates that the source packages
document sheath renewal/release, cell-cycle-coupled secretion, immature
daughters, and the absence of event-rate and daughter-ratio statistics. It
then compares 12 deliberately non-fitted biological-event closures: equal,
broad, and retained-branch/satellite kernels crossed with no-resolved-shell or
released-shell alternatives and deep or partial daughter-state resets.

The calculation is a prospective, event-conditioned protocol design rather
than a fit to microscopy. All closures receive the same synthetic parent cohort
and known two-pulse event schedule; no `beta_bio`, event frequency, daughter
kernel parameter, shell-loss fraction, or reset time is inferred. It reports
how planned time-resolved number/volume distributions, paired
mother--descendant sizes, shell labels, and daughter-state trajectories
separate the alternatives at a declared measurement resolution.

From the workspace root, run:

```bash
python3 CODES/microscopy_model_selection_point5.py
python3 CODES/plot_microscopy_model_selection.py
```

The results include source-derived constraint records, the 12 candidate
definitions, raw event-conditioned lineage records, distribution projections,
pairwise modality-specific discrimination metrics, and a protocol summary
under `CODES/results/`. The shared-style vector figure
`fig11_microscopy_model_selection.pdf` and PNG preview are written under
`Draft_2026-09-07/figures/`. A complete pairwise separation is a statement
that the proposed joint experiment is informative; it does not select a model
from the existing microscopy.

## Biological-versus-mechanical stress-ramp identification test

`biological_mechanical_stress_ramp_point6.py` is a prospective, seeded
event-level benchmark for the experiment required to identify the two
separation channels. It is not a fit to a published stress-ramp dataset. Five
controlled normalized-stress programs share one recorded maturation history;
the subthreshold reference first fits the maturation-gated biological hazard,
then freezes it while the excess-stress hydrodynamic hazard is fitted on the
remaining training ramps. A held-out up-ramp and a held-out down-ramp receive
no refit.

The proposed imaging readout records branch topology and a shell label. A
retained mother plus new daughter identifies a biological event, whereas two
free, material-conserving fragments identify a hydrodynamic event. Those
event-resolved observations fit the two daughter kernels separately. The
recovered rates and kernels then drive the local Appendix-A log-normal center
and width diagnostic at fixed maturation; no terminal histogram is fitted.

From the workspace root, run:

```bash
python3 CODES/biological_mechanical_stress_ramp_point6.py
python3 CODES/plot_biological_mechanical_stress_ramp.py
```

The scripts write reproducible at-risk, event, branch, rate-prediction,
kernel-density, log-normal-response, parameter-recovery, and summary products
under `CODES/results/`. The common-style vector figure
`fig12_biological_mechanical_stress_ramp.pdf` and PNG preview are written
under `Draft_2026-09-07/figures/`. The known benchmark generator is retained
only to quantify protocol recovery; its numerical values are not biological
estimates.

## Source-bound optics-to-productivity test

The optics_to_productivity_point7.py analysis evaluates the currently
accessible Kemel et al. endpoint tables without creating a fictitious
size-resolved optical dataset. The two reported optical values (8.0 and
79 m2 kg-1) have no paired diameter or reduced-state records, while the
145- and 434-um productivity cultures have no paired mass-absorption
measurements. The script therefore reports a source-supported identification
boundary rather than fitting a continuous a_abs(D,a) curve.

With the colony-dynamics solution and yield terms held fixed, the
common-yield optical null has an exact no-refit invariant: its medium/large
biomass and hydrocarbon productivity ratios must be equal. The Kemel targets
are 1.34 and 1.74, so the null is structurally rejected. The derived
hydrocarbon-to-biomass productivity fraction is 1.299 times higher for the
medium than the large culture. This establishes the need for a state- or
product-specific hydrocarbon-allocation route, but does not fit one.

The code also writes clearly labelled conditional same-state optical curves:
they demonstrate how the unreported disrupted diameter would control an
apparent optical exponent, and must not be interpreted as observations or
fits. A small co-registered optical measurement grid is included as the
minimum prospective design: four diameter bins by two measured state strata,
followed by the two Kemel productivity cultures as reserved no-refit targets.

From the workspace root, run:

    python3 CODES/optics_to_productivity_point7.py
    python3 CODES/plot_optics_to_productivity.py

The source audit, conditional sensitivity records, ratio test, proposed
measurement grid, and interpretation boundary are written under CODES/results
with an optics_to_productivity_ prefix. The common-style vector figure
fig13_optics_to_productivity_test.pdf and a PNG preview are written under
Draft_2026-09-07/figures/.

## Shared plot style

`paper_plot_style.py` is the common white-background, light-grid,
blue/orange/green Matplotlib style.  Both project-generated quantitative
figure scripts import it, so their typography, colour palette, grids, and
vector-output settings remain consistent.
