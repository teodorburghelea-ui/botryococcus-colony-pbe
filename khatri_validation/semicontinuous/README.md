# Khatri et al. (2014): no-refit semi-continuous simulation with frozen Zhang closures (26 Sep 2026)

This package replaces the retired `run_validation.py`. It integrates the matched Zhang PBE (`../../zhang_revision/analysis.py`, `Solver`) for the four Khatri shake-flask cultures (5, 7.5, 10 and 12.5 % daily replacement). It uses the frozen primary coefficients of each daughter law. **Nothing is fitted to Khatri data.**

## Inputs (`khatri_inputs.csv`, from `../source/Khatri_2014.pdf`)
- **Steady dry weight (Fig. 5):** 23.0, 18.1, 13.3 and 9.3 g/L. The text gives 22.3 g/L at 5 %.
- **OD550 per gDW/L (Fig. 6, diamonds):** 0.181, 0.090, 0.084 and 0.077. The text gives 0.18 and 0.076.
- **Mean colony diameter ± SD (Fig. 6, bars):** 0.093, 0.179, 0.224 and 0.343 mm, from image analysis of 40–150 aggregates.
- **Culture conditions:**
  - 50 mL of culture in 500-mL flasks, fed every other day.
  - Light was 225 µmol m⁻² s⁻¹, rising to 325 µmol m⁻² s⁻¹ from day 36.
  - Photoperiod was 15 h light / 9 h dark.
  - Every culture started from 15 g/L and was steady by about day 20.
  - Colony size was measured over days 50–64.

## Chain and declared assumptions (`run_semicontinuous.py`)
1. **Biomass.** X(t) rises linearly from 15 g/L to X_ss over 20 days, then stays constant.
2. **Flask light.** The volume-averaged PAR is I₀(1−e^{−τ})/τ, with τ = ln10 · (OD/gDW) · X · L and a liquid depth L = 0.58 cm (50 mL over an 86.6 cm² base).
3. **Conversion to the Zhang light signal.** 6.05 V (Zhang's empty-column average at 10 klx) is taken to equal 200 µmol m⁻² s⁻¹ (halogen light, about 20 µmol m⁻² s⁻¹ per klx). **Zhang gives no volts-to-PAR calibration**, so this conversion is varied by a factor of 2 in each direction.
4. **Light over the day.** The primary run applies the 15/24 daily mean of the light. An explicit on/off cycle is run as a sensitivity case.
5. **Growth.** Primary: the colony material growth rate equals the specific growth rate that the steady biomass balance requires, μ = −ln(1−2f)/2. Alternative: the frozen Zhang light-driven growth law.
6. **Separation.** b·S(D)·exp[η(½−a_eq)] + 0.02·S(D), with a_eq = I/(I+0.7). Size-only runs set η = 0. The finite-time model equals the instantaneous one because τ = 0 was fitted.
7. **Removal.** Uniform withdrawal does not change the normalized distribution, so the replacement schedule acts only through μ and the light.
8. **Inoculum.** Unknown. The Zhang 3-klx day-0 histogram is used, with the 10-klx one as a sensitivity case.
9. **Observable.** Number-mean diameter of colonies ≥ 0.02 mm, averaged over days 50–64. Thresholds of 0 and 0.05 mm are sensitivity cases.

`tau_scan.py` keeps the frozen instantaneous coefficients and varies only the relaxation time τ (0–8 d) under the explicit light/dark cycle. This is a scan, not a fit. The script also records whether a steady colony size exists: the maximum separation rate must exceed the growth rate.

## Results (`summary.json`, `semicontinuous_results.csv`, `tau_scan_results.csv`, `khatri_semicontinuous.pdf`)

The observed diameter ratio between 12.5 % and 5 % replacement is 3.69. A flat line at the observed mean has RMSE 0.090 mm.

**Direction of the trend.** Colony size increases with replacement rate in every run: all scenarios, laws and models.

**Primary scenario (daily-mean light):**

| Closure | Predicted ratio | Share of observed log trend | RMSE (mm) |
|---|---|---|---|
| Instantaneous light | 2.16–2.66 | 59–75 % | 0.104–0.113 |
| Size only | 1.46–1.83 | 29–46 % | 0.135–0.180 |

- The size-only trend comes purely from the growth route: faster steady growth (μ = D) against a size-gated separation rate.
- **Absolute sizes are under-predicted** at high replacement: 0.17 mm against 0.34 mm observed. Neither closure beats a flat line at the observed mean in RMSE.

**Explicit light/dark cycle.**
- Instantaneous light response gives a ratio of only 1.21–1.24, *below* size-only.
- A relaxing state (τ from 0.1 to 8 d) raises the ratio only to 1.3–1.57.
- The light-response advantage in the primary run therefore depends on the assumption that the state responds to the *daily-mean light* rather than tracking the light/dark cycle. The Zhang data, taken under continuous light, cannot settle this.

**Optical feedback.**
- Holding OD/gDW at 0.181 removes the measured optical change. The instantaneous ratio then falls from 2.16–2.66 to 1.85–2.08.
- The measured size-dependent optical change is part of the predicted effect.

**No steady state at 12.5 %.** For the equal and broad binary laws the maximum separation rate is below the required growth rate at 12.5 %. Colonies therefore keep enlarging after day 64: the day 136–150 means rise to 0.21–0.22 mm. The day 50–64 values are the model state at the measurement window, not a steady state.

## Reproduce
```
OPENBLAS_NUM_THREADS=1 python3 run_semicontinuous.py
OPENBLAS_NUM_THREADS=1 python3 tau_scan.py
python3 plot_semicontinuous.py
```
The runs take about 15 minutes on one core.
