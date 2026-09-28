# Pre-declared cross-genus test of the light–size exponent (Microcystis)

Written 28 Sep 2026. Only the abstracts had been read when this protocol was written; the full texts had not been obtained and no value had been extracted.

## Sources

- **Xu et al. 2026**, *Environ. Pollut.*, doi:10.1016/j.envpol.2026.128720. Natural Microcystis communities at five light levels (0–216 µmol m⁻² s⁻¹) for 54 days. This is the primary source.
- **Xu et al. 2023**, *Water Res.*, doi:10.1016/j.watres.2023.119839. Several light intensities. This is the secondary source.

## Prediction

- Where light supports positive growth, colony size scales with light as D ∝ I^k, with k in 0.4–0.6. This is the *B. braunii* range: 0.41 bubble column; 0.53 torus photobioreactor; 0.47–0.67 warm continuous cultures.
- Only the empirical exponent is tested. The release-size mechanism is not, because Microcystis colonies form by division within an EPS matrix and fragment by erosion.

## Extraction and analysis (fixed now)

**Size values.**
- Take the colony-size statistic reported by the source; record whether it is a median, mean or D50, and whether it is number- or volume-weighted.
- Use each light level at the last sampling time, and also at the time of maximum between-treatment contrast.

**Light.**
- Use the incident light as reported.
- If attenuation information (biomass and depth) is given, also compute the depth-averaged light.

**Which light levels enter the fit.**
- Exclude levels with zero or negative growth (including 0 µmol m⁻² s⁻¹): the law is stated only where light supports growth.
- Exclude the highest level if the source reports photoinhibited or suppressed growth there.

**Fit.**
- Fit the slope of ln D against ln I across the remaining levels (at least 3 levels are required).
- Report the slope with a 95 % interval from the reported error bars, if available; otherwise from the residuals.

## Criteria

| Outcome | Condition |
|---|---|
| **Pass** | the 95 % interval of the slope overlaps 0.4–0.6 and excludes 0 |
| **Fail** | the interval excludes 0.4–0.6 |
| **Inconclusive** | the interval includes both 0 and 0.4–0.6, or fewer than 3 usable levels |

Time scale (reported only, not pass/fail): if size time series are given, fit a first-order approach of ln D to its final value for each light level and report the time constant.
