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

---
## Results (recorded 28 Sep 2026, after extraction and fit; `results/external/microcystis_exponent_test.json`)

The protocol wording ("95 % interval from the reported error bars") ignores lack of fit. With three or four light levels that are not on a line, it gives unrealistically narrow intervals. All three interval variants are therefore reported.

### Xu 2026 (primary): day 54, levels 18/54/108 µmol m⁻² s⁻¹

- **Excluded:** 216 µmol m⁻² s⁻¹ (early growth suppressed) and the dark control.
- **Slope:** 0.62 weighted (0.56 OLS). Pairwise exponents: 0.99 (18→54) and −0.23 (54→108). Reduced χ² = 362.

| Interval method | 95 % CI | Verdict |
|---|---|---|
| Error bars only (protocol wording) | [0.56, 0.67] | PASS |
| Scaled by lack of fit | [−5.8, 7.0] | INCONCLUSIVE |
| OLS | [−3.6, 4.7] | INCONCLUSIVE |

### Xu 2023 (secondary): day 42, all four illuminated levels

- **Slope:** 0.26. Pairwise exponents: 0.23, 0.51 and −0.08. Reduced χ² = 12.

| Interval method | 95 % CI | Verdict |
|---|---|---|
| Error bars only | [0.22, 0.30] | FAIL |
| Scaled by lack of fit | [−0.04, 0.55] | INCONCLUSIVE |
| OLS | [−0.04, 0.53] | INCONCLUSIVE |

## Overall verdict

Overall the result is **INCONCLUSIVE**. The literal protocol gives opposite verdicts for the two studies (primary PASS, secondary FAIL), and both become inconclusive once lack of fit is accounted for.

Qualitatively, Microcystis colony size increases with light at low to moderate irradiance (pairwise exponents 0.23–0.99) and saturates or declines above about 108 µmol m⁻² s⁻¹. The 2026 size trajectories are non-monotonic in time and track EPS phases, so the day-54 values are not steady states.

This is not evidence for a common cross-genus exponent. It is compatible with a light–size increase within a window of moderate light.
