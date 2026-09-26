# Matched biological-production ablation

> **Update 25 Sep 2026.** Inputs are now point re-extractions of van den Berg et al. (2019) Fig. 2 (`data/vdb_reextraction_2026-09-25/`). They replace the hand-shaped smooth curves, which are kept as `data/source_distributions_archived_smooth_reconstruction.csv`. All results were regenerated. Outputs from before 25 Sep are kept in `results_archive_2026-09-22/`. Numbers quoted further below in this README (the table, the LL example, the convergence values) describe that archived run.
>
> Current mean excluded-condition TV:
>
> | Daughter law | Size only | Instantaneous | Finite time |
> |---|---:|---:|---:|
> | Retained satellites | 0.217 | 0.206 | 0.220 |
> | Asymmetric binary | 0.176 | 0.229 | 0.227 |
> | Broad binary | 0.172 | 0.191 | 0.193 |
>
> Persistence gives a mean of 0.159, and 5 of 27 predictions beat it.
>
> Maximum TV differences in the convergence checks: half timestep 6.5e-6, grid ×2 0.0096, grid ×4 0.0144, expanded domain 6.5e-9.

This package replaces the confounded low-light mechanism comparison in the acceptance draft. All three comparators retain biological production. It tests whether adding a finite relaxation time improves excluded-condition predictions under explicitly matched mechanisms.

**Result:** the completed comparison does not demonstrate a consistent predictive benefit from finite-time relaxation. After equal-budget refinement, mean excluded-condition TV is:

| Daughter law | Size only | Instantaneous | Finite time |
|---|---:|---:|---:|
| Retained satellites | 0.221 | 0.234 | 0.231 |
| Asymmetric binary | 0.205 | 0.206 | 0.214 |
| Broad binary | 0.213 | 0.206 | 0.228 |

Smaller is better. Each number averages three fits, each excluding a different light condition. No daughter law was selected using test errors. These are retrospective predictions of smooth reconstructed distributions, not original-data statistical validation.

## What is matched

- The same observed starting volume distribution, converted to numbers, and the same volume-bin observation operator. Exterior predicted mass remains in the score.
- Retained-material growth `G_x = 0.07 x d^-1`, with matched refits at 0, 0.06 and 0.08 d^-1 as sensitivities.
- The same size gate, weak equal-material mechanical breakup, zero aggregation and zero shell loss.
- The same parent-relative biological daughter law within each comparison. All event laws conserve material and keep every daughter at or below its own parent size.
- Equal training information: fit two light conditions, predict the excluded third, preserving LL/HL replicate identities and equal weight per light condition.
- Equal search allowance: 512 common-sequence global candidates plus 96 training-only local candidates per model/law/fold. The models have one, two and three adjustable coefficients; equal effort does not imply equal complexity or a proven global optimum.

The three biological laws are retained satellites, 10:90 binary partitions, and a broad distribution of conservative binary partitions. The first implies approximately 24.98 expected descendants including the retained branch. That is a declared many-satellite hypothesis, not a count inferred from microscopy.

The finite-time response uses inherited daughter state and a point-mass initial state. The resulting structured solution has an exact marginal reduction with a common deterministic state trajectory. This isolates relaxation from reset and heterogeneity; **it does not test all full two-coordinate PBE mechanisms**. A one-coordinate model with a prescribed time-dependent biological rate could reproduce that marginal. The state variable and its light-equilibrium map are hypotheses, not measured physiology.

## Why local refinement was added

The initial protocol was fixed before the new searches. Prefix checks at 128, 256 and 512 candidates changed model rankings, so every primary fit received the same additional 96 bounded coordinate-poll evaluations. This extension uses training error only and is recorded in `results/refinement_protocol.json`.

For the satellite-law LL fold, the initial instantaneous/finite-time prediction errors were 0.431/0.286. Improving their training fits with equal effort changed prediction errors to 0.392/0.381. This is a useful result: the apparently large benefit was sensitive to nearly equivalent training parameter choices. The size-only model gives 0.257 for that fold. The low-light ordering is unchanged on the finest tested grid. Near-optimal global candidate ranges are diagnostic sensitivity sets, not confidence intervals.

## Verification

Five scientific test groups pass: event material/number/support identities; common initialization and observation; instantaneous and neutral-response limits; nonzero biological production in every comparator; and independent analytic no-event growth transport.

For the final refined parameters, maximum distribution differences from the main discretization are:

| Check | Maximum TV difference |
|---|---:|
| Half timestep | 0.0000381 |
| Double spatial resolution | 0.00966 |
| Quadruple spatial resolution | 0.01431 |
| Expanded domain | 0.00000000651 |

The event-matrix material residual is below 3.5e-16. The primary time-integrated material error relative to analytic exponential growth is approximately 1.14e-5; it decreases under timestep refinement and is not reported as a machine-precision time-integration result. The largest material fraction in sections with boundary-suppressed biological events among the primary final fits is 7.24e-6. Spatial checks freeze fitted coefficients; they do not certify a reoptimized continuum fit. Growth-only TV error against an independent analytical translated-bin solution decreases from 0.0274 to 0.00347 under refinement.

## Reproduce

From this directory, with Python, NumPy, SciPy and Matplotlib installed:

```bash
python3 test_ablation.py
python3 run_ablation.py
python3 polish_fits.py
python3 report_ablation.py
python3 verify_outputs.py
```

`run_ablation.py` reads the archived protocol and data, runs the global searches and growth sensitivities, and records diagnostics. `polish_fits.py` adds equal-budget refinement and repeats numerical checks on the refined fits. `report_ablation.py` generates the final primary table/figure, frozen-coefficient ablation and search/growth sensitivity reports. The report script also regenerates the parameter table from the full-precision fitted CSV. `verify_outputs.py` independently recomputes all archived replicate and fold scores and checks search budgets, normalization and input hashes.

To restore/reapply the manuscript changes from the saved pre-ablation source, run `python3 update_manuscript.py`. This deliberately rebuilds the acceptance source from that saved snapshot; use it only when later editorial changes do not need preserving. Compile from the draft root using `latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript_acceptance_revision.tex`.

## Files to inspect

- `results/polished_fits.csv`: final fitted coefficients, training/test errors and exact candidate counts.
- `results/polished_replicate_metrics.csv`: individual replicate TV and Wasserstein scores.
- `results/matched_ablation.pdf`: final distributions and excluded-condition prediction errors.
- `results/polished_convergence.csv`: checks on those actual final predictions.
- `results/frozen_coefficient_ablation.csv`: finite-time vs instantaneous comparison with the finite-time fitted rate coefficients held fixed.
- `results/candidate_scores.csv`, `results/refinement_trace.csv`: complete global and local search records.
- `results/search_budget_check.csv`, `results/near_optimal_candidates.csv`: initial-search sensitivity and weakly separated fits.
- `results/growth_sensitivity.pdf`: separately refitted 256-candidate growth scenarios; the 0.07 case uses the same 256-candidate prefix, not the refined primary fits.
- `methods.tex`: full manuscript methods and scope.
- `protocol.json`, `results/refinement_protocol.json`, manifests: declared choices, extension and provenance.

## Source and scope

The input CSVs are copies of the existing local van den Berg reconstruction package. That package describes its distributions as smooth log-normal/mixture reconstructions matched to visible peaks, widths and heights. We did not redigitize or replace those targets. The [primary paper](https://doi.org/10.1104/pp.18.01499) was checked for the split-culture design and light conditions. Its longer low-light experiment is not treated as a continuation of the same tracked population or used as an invented time-series holdout.

This addresses the requested matched-comparator analysis. Independent time-resolved data, daughter-state-reset tests and empirical physiological-state identification remain separate scientific tasks. No newly inferred physiological parameter, confidence interval or experimentally validated superiority is claimed.
