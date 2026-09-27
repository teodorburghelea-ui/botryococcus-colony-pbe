# Confirmation protocol for the quasi-steady light-dependent threshold closure

Written 2026-09-27T09:21:12+02:00, BEFORE running any of the tests below.

The closure was proposed after inspecting transfer-set errors of the earlier closures. The primary transfer scores (TV 0.257–0.280) are therefore **not** blind. These checks are declared in advance to test whether the result is robust and light-specific.

## Model (frozen)

- **Separation rate:** beta_bio = b S(D; D_c(I)).
- **Light-dependent threshold:** D_c(I) = D_c0 (I / 2 V)^k.
- **Hydrodynamic term:** 0.02 S(D; D_c0).
- **Growth:** 3 g_D a_eq x.
- **Predicted distribution:** the stationary shape at the current light I(t).
- **Bounds:** 10^-2.5 ≤ b ≤ 1 d^-1; 0.035 ≤ D_c0 ≤ 0.5 mm; -1.5 ≤ k ≤ 1.5.
- **Fitting:** same budget and split as the primary analysis.

## Tests and success criteria

1. **Extraction robustness** (6 alternative extractions × 3 laws).
   Success if the law-averaged transfer TV is at or below that of the static log-normal in at least 5 of 7 scenarios.
2. **Light specificity** (9 derangements of the separation-threshold light × 3 laws; growth keeps the correct light).
   Success if the correct light beats every scrambled assignment for each law.
3. **Independent culture system** (Khatri et al. 2014; frozen primary coefficients; daily-mean light; growth from the biomass balance; no refit).
   Report the predicted 12.5 %/5 % size ratio and the share of the observed log trend. Direction must be reproduced.
   This is not blind: the observed trend is known.
4. **Sensitivity:** the upper bound on b raised to 10 d^-1 (primary scenario, 3 laws). Report only.

The outcome is reported as found, whether it confirms the closure or not.
