# Re-digitization of Zhang & Kojima (1998), Fig. 3 — 25 September 2026

Replaces, for all future analyses, `../zhang_kojima_1998_fig3_digitized.csv`. The old file is kept unchanged for provenance.

**Main file:** `zhang_kojima_1998_fig3_redigitized.csv`. It has one row per histogram × diameter marker (25 histograms × 12 markers) with these columns:

| Column | Content |
|---|---|
| `figure_panel`, `preculture_irradiance_klx`, `lighted_volume_percent`, `days` | Which histogram |
| `marker` | Symbol used in the source |
| `diameter_marker_mm` | Plotted x position of the marker: 0, 0.05, …, 0.55 mm |
| `frequency_percent_raw` | Frequency as plotted |
| `frequency_percent_normalized` | Rescaled so each histogram sums to 100 % |
| `sum_raw_percent` | Sum of the plotted values for that histogram |
| `note` | How each value was obtained |

## Source and what it contains
- **Source file:** Zhang, K., Kojima, E. (1998) *J. Ferment. Bioeng.* 86:573–576, Fig. 3 (p. 574). Local file: `Paper_Draft/PaperPlan.pdf`, MD5 `446956b19a41f2ac9031ee130b539046`. Note the misleading file name.
- **Histograms present:** 25.
  - 3-klx preculture: days 0, 4 and 16 at V_L = 5/25/50/100 %.
  - 10-klx preculture: days 0, 5 and 18 at all four V_L; day 25 **only at V_L = 100 %**.
- **Removed:** the archived day-25 histograms at V_L = 5/25/50 % do not exist in the figure.

## Method (scripted equivalent of WebPlotDigitizer point extraction)
1. **Native image, no resampling.** `extract_page.py` rebuilds the 300-dpi bilevel scan of page 574 from the PDF's embedded CCITT strips (`pdfimages`).
2. **Axis calibration** (`axis_calibration.csv`).
   - For every panel, the four frame-line centres were measured from row and column projection profiles.
   - Frame edges correspond to x = 0 and 0.6 mm and to y = 0 and 100 %. Each panel is calibrated separately to absorb a small scan skew of a few px.
   - Checks: the 50 % tick lies at the frame midpoint (±2 px), and markers fall on the 0.05-mm grid (≈32.5 px spacing; see overlays).
3. **Marker detection** (`detect.py`).
   - Shape templates (open circle, open square, open triangle, filled diamond; 24 px, 3-px stroke) are scored only in narrow columns (±5 px) around each 0.05-mm grid position.
   - Frequency is taken from the fitted marker-centre height. Triangles get a 1.5-px centring correction.
4. **Assembly** (`assemble.py`, `finalize.py`).
   - Where two marker types fire at the same position, the better-scoring one is kept.
   - A series with no marker above the baseline at a grid position is set to 0 %.
5. **Occluded markers** (`refine.py`). Four markers hidden behind others were measured by local template fits (stroke scores 0.87–0.98):
   - 3a-B day 0 at 0.10 mm;
   - 3a-D day 4 at 0.10 mm;
   - 3b-D day 5 at 0.10 and 0.20 mm.
   One value is **inferred**, not measured: 3a-D day 0 at 0.20 mm. That circle is completely hidden inside a same-size square, so it was set to 30.2 %, the value in the identical day-0 histograms of panels A–C.
6. **Small values** (`small_values.py`).
   - Near the baseline, the thick frame line biases heights by up to about 2 %.
   - Small values were therefore measured relative to each series' own baseline marker height.
   - A correction was applied only where the histogram was short and the raised marker was confirmed visually: 3a-C day 4 at 0.25 and 0.30 mm, and 3b-C day 5 at 0.30 mm.
7. **Quality control:** `qc/final_overlay_a.png` and `qc/final_overlay_b.png` plot every extracted value on the scan. All 25 histograms were inspected.

## Validation
- **Internal consistency.**
  - Plotted histograms sum to 93.7–104.8 % (median 99.2 %).
  - The day-0 histograms are identical across the four V_L panels of each preculture, as expected for a split inoculum:
    - 3 klx: 27.4 / 47.1 / 30.2 % at 0.10 / 0.15 / 0.20 mm;
    - 10 klx: 39 / 60 % at 0.05 / 0.10 mm.
  - The source's own day-0 values sum to about 105 %, so they are normalized before use.
- **Independent cross-check against Fig. 5** (`qc/fig5_crosscheck.txt`).
  - D₃,₀ = (Σ n d³ / Σ n)^{1/3} from the re-digitized histograms reproduces the separately digitized Fig. 5 values for all 25 runs: RMSE 0.0115 mm, r = 0.996.
  - The archived digitization gives RMSE 0.0463 mm and r = 0.963.
  - The source text reports initial sizes of 0.15 and 0.07 mm; the re-digitization gives 0.160 and 0.087 mm, against 0.17–0.20 and 0.12 mm in the archive.
- **Bin convention.**
  - A common shift of marker positions by −0.011 mm reduces the Fig. 5 RMSE to 0.0055 mm, which is comparable to the Fig. 5 reading precision.
  - The source does not define its size classes. Use the plotted marker positions as the primary convention and a −0.011 mm (or −0.025 mm upper-edge) shift as a sensitivity case.
- **Difference from the archive.** Over the 25 genuine histograms, mean TV = 0.177 (max 0.392, at 10 klx / 100 % / day 18). See `qc/comparison_redigitized_vs_archived.png`: the archived curves are systematically broader and shifted right.

## Consequences for the manuscript and analysis
1. **Figs. 3 and 5 are the same runs.** The archive's statement that they are "separate source datasets" with different initial-condition designs is not supported: the moments agree to within reading precision. Fig. 5 therefore provides paired scalar observations at additional days (8, 9, 12, 13, 21) and can be used as a legitimate extra validation target.
2. **All Zhang results must be recomputed from this file.** This covers the fits, transfer TV, benchmarks, envelope, figures and tables. The earlier corrected numbers (0.376–0.392 etc.) are based on the old digitization and are superseded.
3. **Uncertainty to state in the methods.**
   - About ±1 percentage point for clearly visible markers: marker-centre precision of about ±4 px ≈ ±1 %.
   - About ±2 percentage points for values near the baseline.
   - Class boundaries are unknown (see the bin convention above).
   - Counting noise: the source measured more than 100 colonies per histogram.

## Reproduce
```
python3 scripts/extract_page.py ../../../../../../PaperPlan.pdf   # or any path to the source PDF
python3 scripts/detect.py && python3 scripts/assemble.py && python3 scripts/refine.py
python3 scripts/small_values.py && python3 scripts/finalize.py
python3 scripts/compare.py && python3 scripts/plot_compare.py
```
Requirements: Python 3 with numpy, Pillow and matplotlib, plus poppler-utils. The page image is regenerated locally and is not redistributed.
