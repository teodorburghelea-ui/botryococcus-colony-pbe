# Matched Zhang time-series revision

> **Update 25 Sep 2026.** The package now reads the re-digitized Fig. 3 (`data/zhang_kojima_1998_fig3_redigitized.csv`; method in `data/redigitization_2026-09-25/README.md`). All 90 fits were rerun: 10 scenarios, including the new class-convention scenario `bins_minus_0p011`.
>
> - **Transfer set.** 9 genuine 10-klx histograms. The three non-source day-25 panels are gone.
> - **Fig. 5.** It now serves as a *paired* scalar check, because Figs. 3 and 5 are the same runs. The "separate source datasets" statement below is superseded; see `results/digitization_audit.md`.
> - **Counting noise.** New script `counting_noise.py` adds a multinomial counting-noise floor.
> - **Legacy figure.** The legacy Fig. 5 figure is no longer generated.
> - **Superseded outputs.** The old package and results are frozen in `../zhang_revision_archive_2026-09-22/`.
>
> Primary transfer TV:
>
> | Model | Transfer TV |
> |---|---|
> | Instantaneous | 0.366–0.397 |
> | Size-only | 0.438–0.461 |
> | Persistence | 0.553 |
> | Growth only | 0.535 |
>
> The η coefficient is at its upper bound (3) in all primary light-response fits.
>
> Run `python3 counting_noise.py` after `report.py`.

This package implements priority 1 of the second systematic manuscript review: matched models on the existing Zhang number-distribution trajectories, no-fit persistence and analytical growth-only benchmarks, reconstruction/closure sensitivity, numerical verification, corrected legacy scores, and manuscript integration.

## Scientific scope

The primary training objective uses the eight non-initial Fig. 3 histograms of the 3-klx preculture, giving equal weight to the four lighted-volume ratios and their two later observations. Twelve non-initial histograms from the 10-klx history are predicted with frozen coefficients. These are retrospective source-graph reconstructions, not newly collected data or a blind validation experiment.

Three rate variants (size-only, instantaneous light, finite-time light) are compared under three conservative binary daughter laws. All retain growth and biological and mechanical separation. Each of nine fits receives 512 global plus 128 local, training-only evaluations. Each of five reconstruction perturbations and three growth/initial-state alternatives receives exactly the same fitting allowance: **81 fits total**. The daughter laws and scenarios are not selected using transfer scores.

All input choices and sensitivity transformations were recorded in `protocol.json` before the new searches. The nested fitted simpler model occupies one of the richer model's global candidate slots and uses only training information. Equal search budgets do not prove global optimality or equal complexity.

The state has homogeneous initialization and exact inheritance. This tests a delayed light response, not a heterogeneous two-coordinate population, daughter-state resets, or a measured ECM state. In particular, the light-response rate is a different hypothesis from the main manuscript's increasing maturity gate.

## Source audit

- The locally available original Zhang–Kojima paper was read, and its Fig. 3 page was visually inspected. Source: **Zhang and Kojima (1998), Journal of Fermentation and Bioengineering 86, 573–576**, DOI [10.1016/S0922-338X(99)80009-9](https://doi.org/10.1016/S0922-338X(99)80009-9).
- Fig. 3 measures number frequencies of microscopic Feret diameters. Fig. 5 reports the cubic mean, D30. Neither is the volume-weighted D43 or the van den Berg volume distribution.
- Exact local copies of the archived input CSVs and their source paths/hashes are in `data/`. No new redigitization is claimed. The old package describes approximate manual extraction from 600-dpi rasterized plots, with frequencies normalized to 100%.
- **Figure 3 and Figure 5 are separate source datasets.** Figure 3 provides run-specific number histograms across lighted-volume ratios, while Figure 5 provides two scalar trajectories with explicitly different initial-size cases. The paper does not provide run-level linkage between these panels, so the scalar traces are retained as diagnostics rather than paired validation targets. `data/figure_observable_registry.csv`, `data/digitization_coordinates.csv`, `results/source_moment_audit.csv`, and `results/source_audit.md` record the observable definitions, axis calibration, and provenance decision. Exact cross-figure validation requires the original numerical records or author confirmation.
- Fig. 5 provides a supplementary scalar check of the new histogram-trained models, not another independent experiment. The original Fig. 5 fit is retained separately as a legacy diagnostic; its updated scores exclude supplied initial diameters.
- Light is prescribed from the source's sparse biomass table and detector-voltage attenuation law. Thus even transfer predictions are conditional on measured biomass-derived light, not autonomous forecasts of reactor operation.
- Source sampling/medium replacement is acknowledged in the methods. No unmeasured selective colony retention, biomass-to-colony-growth identity, or stress-to-rate relationship is inferred.

## Reproduce

Python 3, NumPy, SciPy and Matplotlib are required. No network access or additional data download is needed.

From this directory:

```bash
OPENBLAS_NUM_THREADS=1 python3 test_analysis.py
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario primary
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario shift_left
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario shift_right
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario broaden
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario tail_up
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario tail_down
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario growth_half
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario growth_one_and_half
OPENBLAS_NUM_THREADS=1 python3 analysis.py --scenario initial_equilibrium
OPENBLAS_NUM_THREADS=1 python3 check_numerics.py
OPENBLAS_NUM_THREADS=1 python3 report.py
OPENBLAS_NUM_THREADS=1 python3 verify_results.py
python3 stage_manuscript.py
```

The analysis skips a completed fit when its final JSON exists. To recompute after changing code, inputs or protocol, first preserve/move the corresponding results directory; do not silently combine outputs from different scientific protocols. Run models in their declared order because richer models use their fitted simpler predecessor as a training-only seed.

`stage_manuscript.py` prepares a manuscript in this directory from the saved pre-revision source. It does not overwrite the live manuscript. Reapplying it after subsequent editorial changes requires merging those changes rather than blindly copying the staged file.

Compile from the draft root, using the stage source and a writable output directory, or compile the installed source with `latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript_acceptance_revision.tex`. The source includes this package's methods and generated tables directly.

## Outputs

- `results/summary.csv`: every scenario, law and model, including benchmarks and supplementary scalar RMSE.
- `results/fits.csv`: full-precision fitted coefficients and training/transfer scores.
- `results/<scenario>/*_trace.csv`: all 640 evaluated candidates and training losses.
- `results/<scenario>/*_curves.csv`, `*_metrics.csv`, `*_scalars.csv`: complete distributions, individual errors, and scalar observations/predictions. Initial entries are retained for auditing, but excluded from aggregate scores.
- `results/benchmark_metrics.csv`: persistence and exact growth-only distribution errors, under each scenario.
- `results/sensitivity_contrasts.csv`: paired finite-minus-instantaneous and model-minus-persistence transfer contrasts.
- `results/numerical_checks.csv`: frozen fitted-parameter timestep, grid and domain checks on the actual solver.
- `results/analytic_growth_checks.csv`: independent analytical transport comparison under refinement. These checks use an enlarged upper domain to separate transport accuracy from upper-boundary suppression.
- `results/verification.json`, `results/manifest.json`: independently recomputed output checks and file hashes.
- `results/matched_zhang.pdf`: every transfer observation time and all three laws.
- `results/distributions_<law>.pdf`: all twelve excluded-history histograms, without choosing favorable panels or laws.
- `results/reconstruction_sensitivity.pdf`: paired sensitivity contrasts.
- `results/legacy_zhang_corrected.pdf`: regenerated legacy display with initial scalar points removed.
- `methods.tex`, generated tables, and staged manuscript: the concrete paper changes.

The existing light-acclimation report script in the parent directory additionally includes persistence in its table, figure, and summary. Its fitted simulations are unchanged. Its outputs can be regenerated using `python3 report_ablation.py`, then checked using `python3 verify_outputs.py` from the parent directory.

## Interpretation

Reconstruction alternatives are deliberately specified shape stress tests, not data-derived error bars or confidence intervals. A score reduction should not be interpreted as a significance test. The four operating conditions and multiple times are not interchangeable with independent biological replicates.

Numerical checks retain coefficients fitted on the primary mesh. They establish sensitivity of the actual predictions, not convergence of a reoptimized continuum fit. Differences between nearly identical instantaneous and zero-relaxation solutions can arise from finite parameter searches. The strongest conclusions are benchmark-relative transfer and the absence of an established incremental relaxation benefit, not physiological parameter identification.

The original manuscript source is preserved as `manuscript_before_zhang_revision.tex`. Earlier parent-directory update scripts rebuild older snapshots and would overwrite this revision; use the current source and this package for subsequent work.
