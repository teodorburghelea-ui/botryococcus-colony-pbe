# Pre-declared test of the release-size closure on the light-acclimation data

Written 2026-09-27T10:53:56+02:00, BEFORE any fit.

The release-size closure was identified on the bubble-column data only. The van den Berg et al. (2019) distributions were not used to develop it. They were inspected in earlier, unrelated analyses.

## Closure

- **Separation rate:** beta_bio = b × logistic[(D − D_c(I)) / (0.14 D_c(I))], with D_c(I) = 0.25 mm × (I / 150 µmol m⁻² s⁻¹)^k.
- **Why these values:** at I = 150 this reproduces the existing gate (centre 0.25 mm, width 0.035 mm).
- **Fitted coefficients:** b and k (bounds 10^-2.5 ≤ b ≤ 1 d⁻¹, |k| ≤ 1.5). This is the same count as the light-dependent-rate closure (b, η).
- **Unchanged from the existing ablation:**
  - growth 0.07 d⁻¹ (light-independent);
  - hydrodynamic term and daughter laws (retained satellites, asymmetric binary, broad binary);
  - volume-weighted observation;
  - 20-d horizon;
  - leave-one-light-condition-out folds;
  - 512 Sobol global candidates (seed 20260914) plus 96 training-only coordinate polls.
- **Evaluated two ways:** dynamically (from the common Start distribution) and quasi-steadily (stationary shape at the condition's light).

## Criteria (reported as found)

1. **Persistence.** Mean held-out TV over the 9 law × fold combinations below persistence (0.159).
2. **Existing closures.** Mean held-out TV below the size-only (0.188) and light-dependent-rate (0.208) closures.

Both evaluation modes are reported. The mode is not selected by held-out error; the quasi-steady mode is declared primary, as in the main analysis.

---
## Results (recorded 27 Sep 2026, after all runs)

Criteria failed.

| Evaluation | Mean held-out TV (9 law × fold) | Beats persistence (0.159) | Criterion 1 | Criterion 2 |
|---|---|---|---|---|
| Quasi-steady (declared primary) | 0.393 | 0 of 9 predictions | FAIL | FAIL |
| Dynamic | 0.193 | 2 of 9 predictions | FAIL | FAIL (size-only 0.188) |

- Rate closure: 0.208.
- Fitted k is positive (0.0–1.5) in 17 of 18 fits; release size increases with light, as in the bubble column.
- Interpretation: the 20-d light-acclimation distributions remain close to the starting distribution and far from the light-set equilibrium. The quasi-steady assumption fails for this system.
- Output files: `results/release_size/{qs,dyn}_*.json`.
