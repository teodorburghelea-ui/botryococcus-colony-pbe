# Pre-declared external tests of the light–size exponent (written 28 Sep 2026, before any extraction or fit)

**Exponent being tested.** From the bubble-column analysis:
- the empirical law gives D ∝ I^0.41 (geometric-mean number diameter);
- the release-size closure gives D_c ∝ I^k, with k = 0.37–0.50.

Neither new dataset was used to derive these values. Both have been read (Kemel et al. 2025 full text; García-Cubero et al. 2021 Table 2). The tests are therefore not blind to the data, but no quantity below is tuned to them.

## Test A: Kemel et al. (2025), Fig. 8

- **Setup.** Torus PBR, BOT-22, incident PFD 100 vs 440 µmol m⁻² s⁻¹. The source states that biomass concentrations were similar.
- **Extraction.** Median and quartiles of colony diameter at each PFD, read from the boxplots by pixel calibration of the y-axis.
- **Predicted ratio.** At equal biomass the internal light scales at least as the incident light, so the median ratio should be 4.4^k:
  - 1.73–2.10 for k = 0.37–0.50;
  - 1.84 for 0.41.
- **Direction of the bias.** Larger colonies at 440 absorb less per unit biomass, which raises internal light. The observed ratio is therefore expected at or above the prediction.
- **Pass criterion.** The observed median ratio lies within the range spanned by k = 0.3–0.6 (1.56–2.43).

## Test B: García-Cubero et al. (2021), Table 2 (13 non-washout steady states)

- **Light available per unit biomass.** Estimated by Beer–Lambert in the 14-mm flat panel: I_avg = I0 (1 − e^−τ)/τ, with τ = a X L, where X is the measured biomass and L = 0.014 m.
- **Absorption coefficient.** a = 8 m² kg⁻¹ is the primary value: native large BOT-22 colonies (Kemel et al. 2025). The sensitivity cases are 4 and 16 m² kg⁻¹.
- **Model.** ln Dv = c + 0.41 ln I_avg. Only c is fitted, as the mean. Dv is the source's volume-weighted diameter.
- **Comparator.** A constant (c only), in leave-one-out RMSE of ln Dv.
- **Pass criterion.** The fixed-exponent model has a lower leave-one-out RMSE than the constant for the primary a.
- **Also reported, not used for pass/fail:**
  - the freely fitted exponent and its correlation;
  - the same test excluding the 15 °C runs (temperature is not in the model).
- **Caveats.**
  - Temperature varies (15–30 °C) and is not represented.
  - The light regime is 12 h:12 h.
  - Size is volume-weighted, whereas the exponent comes from number-weighted data.

---
## Results (recorded 28 Sep 2026, after the tests)

### Test A (Kemel Fig. 8): PASS
- **Medians:** 319 µm at 100 and 699 µm at 440 µmol m⁻² s⁻¹.
- **Median ratio:** 2.20, within the pass range 1.56–2.43. It is slightly above the k = 0.37–0.50 range (1.73–2.10), as expected from the attenuation bias.
- **Implied exponent on incident light:** 0.53.
- **Quartile ratios:** 1.86 (Q1) and 2.68 (Q3).
- **Source:** values read by pixel calibration from the 300-dpi render (`external/kemel_fig8_*`).

### Test B (García-Cubero, all 13 runs, a = 8 m² kg⁻¹): FAIL
- **Leave-one-out RMSE of ln Dv:** constant 0.467; fixed exponent 0.41 0.487.
- **Free slope:** 0.04 (r = 0.03).

### Reported only (not part of pass/fail)
Excluding the three 15 °C runs, the fixed exponent beats the constant (LOO 0.262 against 0.347). The free slope is 0.62 (r = 0.70); with a = 16 m² kg⁻¹ it is 0.67 (r = 0.85). Temperature, which is not in the model, dominates the cold runs.
