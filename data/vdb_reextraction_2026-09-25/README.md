# Point re-extraction of van den Berg et al. (2019), Fig. 2 — 25 September 2026

Replaces, for all analyses, the smooth log-normal/mixture curves that were hand-shaped to match the figure (kept as `../source_distributions_archived_smooth_reconstruction.csv`).

**Source.** van den Berg, T. E., Chukhutsina, V. U., van Amerongen, H., Croce, R., van Oort, B. (2019) *Plant Physiol.* 179:1132–1143, doi:10.1104/pp.18.01499, Fig. 2 ("20 Days of Photoacclimation"). The image is archived here as `VandenBerg2019_Fig2_PMC.jpg`: the open-access PubMed Central figure file PMC6393799 `PP_201801499DR1_f2.jpg`, 810 × 509 px. The full-text PDF could not be downloaded automatically (bot protection on PMC, OUP and HAL); archive it by hand if needed.

**Output.** `van_den_berg_2019_fig2_points.csv` has one row per series × size class.
- The series are Start (day 0), ML, LL1, LL2, HL1 and HL2 (day 20).
- The 54 log-spaced diameters from 1 to about 1,760 µm are the same grid as the archived reconstruction. The archived 55th point at about 2,030 µm lies outside the plotted frame.
- `volume_percent` is renormalized to 100 % for each series. `volume_percent_raw` is the value read from the plot, and `note` records whether the value was traced or interpolated.
- `../source_distributions.csv`, the input to `run_ablation.py`, is a copy of these values.

## Method (`extract_vdb_fig2.py`)
1. **Axis calibration from tick marks.** On the x-axis, 1, 10, 100 and 1000 µm lie at image columns 99, 314, 529 and 743.5 (215 px per decade). On the y-axis, 0 % lies at row 391 and 14 % at row 15.5 (26.8 px per %).
2. **Colour classification.** Each pixel is assigned to the nearest legend colour. Grey error bars and unsaturated pixels are rejected. The legend, the panel title and the axis line are masked.
3. **Tracing.** For each series and image column, the median row of the matching pixels is taken. That trace is then read at each size class over a ±2-px window.
4. **Hidden segments.** Where a series is hidden under a later-drawn curve, it is linearly interpolated in log D between the nearest visible columns, and the row is flagged. This happens mostly on the near-baseline shoulder below about 50 µm, and for HL1 under HL2.
5. **QC plots.** `qc_overlay.png` shows every extracted point on the source image. `qc_vs_archived.png` compares the new points with the archived reconstruction.

## Accuracy and limitations
- **Reading precision.** About ±0.1 % volume (±3 px) for visible trace segments. Values in the sub-1 % shoulder at 5–50 µm are uncertain by about ±0.3 % because the six curves overlap there.
- **Resolution.** The source image is low resolution (810 px wide), and the laser-diffraction size classes are finer than the 54-point grid. The raw sums (82–88 %) are therefore below 100 % and are renormalized.
- **Replicates.** Plotted replicate error bars are not extracted, so the replicate spread is not used.

## Difference from the archived reconstruction

The archived curves are shifted to larger sizes and have heavier upper tails than the published figure. This is the same failure mode found in the archived Zhang digitization.

| Series | TV (new vs archived) | Median D₅₀, new / archived (µm) |
|---|---|---|
| Start | 0.156 | 321 / 425 |
| ML | 0.155 | 321 / 369 |
| LL1 | 0.093 | 279 / 279 |
| LL2 | 0.156 | 183 / 242 |
| HL1 | 0.237 | 321 / 489 |
| HL2 | 0.271 | 369 / 489 |

All results from `run_ablation.py` / `polish_fits.py` / `report_ablation.py` were regenerated from this file on 25 Sep 2026. The earlier outputs are in `../../results_archive_2026-09-22/`.
