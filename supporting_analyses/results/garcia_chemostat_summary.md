# García-Cubero environmental and chemostat test (Point 4)

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

- Environmental-state closure: calibration RMSE = 25.9 um;
  unused-design-check RMSE = 44.7 um; held-out optimized-condition
  RMSE = 23.2 um.
- Uniform-washout null: calibration RMSE = 47.2 um;
  held-out optimized-condition RMSE = 66.4 um.
- Fitted environmental-state coefficients are theta_0 = -0.1270,
  theta_T = 0.3216, theta_I = 0.2818, and
  theta_D = 0.3433. The declared 22.5 C, 1800 umol m^-2 s^-1
  dilution contrast changes the state-closure prediction by a factor of
  1.99 from D = 0.10 to 0.30 d^-1; the washout-only null is
  exactly diameter-invariant by construction.
- The independent sectional washout control gives maximum normalized-volume
  total variation 7.449e-17 and relative volume-weighted-diameter change
  2.018e-16 across 0.10--0.30 d^-1.

## Interpretation boundary

The positive dilution response is routed explicitly through `q_*`, which in a
full PBE must mean state-dependent biological size selection, state-dependent
retention, or both. It is not evidence that uniform washout changes a
normalized distribution. The small design and the weak unused cold-condition
check do not identify a unique temperature law, relaxation time, or retention
mechanism. This is a structural/calibration test, not a replacement for a
time-resolved chemostat PBE fit.
