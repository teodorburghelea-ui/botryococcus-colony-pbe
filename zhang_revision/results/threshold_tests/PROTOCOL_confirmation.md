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

---
## Results (recorded 27 Sep 2026, after all runs completed)

1. **Extraction robustness: PASS (6 of 7).**
   - Law-averaged transfer TV of the threshold model ranges over 0.255–0.292, against 0.262–0.298 for the static model.
   - The only loss is shift_right: 0.284 against 0.275.
2. **Light specificity: PASS (27 of 27).**
   - The correct light beats every scrambled assignment for every law.

   | Law | Correct light | Best scrambled | Scrambled range |
   |---|---|---|---|
   | Equal binary | 0.257 | 0.355 | 0.355–0.463 |
   | Asymmetric binary | 0.280 | 0.391 | 0.391–0.505 |
   | Broad binary | 0.274 | 0.385 | 0.385–0.468 |

3. **Khatri, no refit: PASS (direction reproduced).**
   - Primary scenario: ratio 1.96–2.81 against 3.69 observed, i.e. 52–79 % of the log trend.
   - RMSE 0.074–0.093 mm against 0.090 mm for the flat line.
   - Steady states exist for every law.
   - The result remains sensitive to the uncalibrated volts-to-photon-flux conversion. Output: `khatri_validation/semicontinuous/threshold_qs_khatri.json`.
4. **Wider b bound (up to 10 d^-1): robust.**
   - Transfer is 0.255, 0.273 and 0.269 (primary 0.257, 0.280 and 0.274).
   - k is stable at 0.37–0.49.
   - b remains above 1 only for the asymmetric law (1.84).
