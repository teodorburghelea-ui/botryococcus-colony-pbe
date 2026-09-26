# Biological-versus-mechanical stress-ramp identification test (Point 6)

## Protocol

- This is a seeded prospective synthetic benchmark, not a fit to a published
  stress-ramp experiment.  It contains 800 at-risk intervals
  and 11121 recorded separation events (6346 training,
  4775 held out).
- Five controlled ramps share one deterministic recorded maturation history.
  Stress is normalized as `tau_h / sigma_ref[a(t)]`, with the threshold fixed
  at one by an independently measured, state-conditioned reference cohesive
  stress.  Thus physiology is not changed while stress is varied.
- The subthreshold training ramp identifies only the maturation-gated
  biological hazard.  That fit is frozen before fitting the excess-stress
  hydrodynamic hazard with the other training ramps.  The held-out mid-up and
  down ramps receive no parameter refit.
- A retained mother plus new branch and a structural shell label identify a
  biological event; two free material-conserving fragments identify a
  hydrodynamic event.  The paired branch maps fit `b_bio` and `b_hyd`
  separately.  In a real experiment these are required observations, not
  assumptions available from a terminal size histogram.

## Results

- The training event set contains 5087 biological and 1259
  hydrodynamic events.  The fitted biological parameters are
  `beta_bio,max = 0.1803 d^-1`,
  `a_c = 0.5960`, and
  `Delta a = 0.0629`.  The fitted hydrodynamic
  parameters are `k_h = 0.4504 d^-1` and
  `q = 1.5203`.
- Across both held-out ramps, no-refit fitted-versus-generator rate RMSE is
  0.0069 d^-1 for the
  hydrodynamic channel and 0.0029
  d^-1 for the biological channel.  The fitted two-branch kernel total
  variations are 0.0125
  for `b_hyd` and 0.0027
  for `b_bio`; the directly labelled shell fraction is 0.080.
- At the fixed reference maturity `a = 0.86`, the
  event-derived local OU projection changes its fitted median from
  276.3 um at
  `tau_h/sigma_ref = 0.65`
  to 232.6 um at
  `tau_h/sigma_ref = 1.85`.
  The corresponding fitted log-standard-deviation changes from
  0.338 to
  0.348.  This is a
  declared local-closure prediction, not a measured distribution response.

## Interpretation boundary

The benchmark establishes protocol adequacy under its declared contrasts: a
fixed-history stress design plus event-resolved branch topology can separate
the two hazards and their daughter kernels, then propagate them into the
Appendix-A log-normal diagnostic without fitting a terminal histogram.  It
does not demonstrate that these numerical rates, kernel shapes, shell loss,
or OU coefficients apply to a real strain.  Without a channel-resolving
lineage/shell readout, aggregate event counts or size distributions would
instead require an explicitly validated latent competing-risk model and would
not justify a separate `beta_bio`, `beta_hyd`, `b_bio`, and `b_hyd` fit.
