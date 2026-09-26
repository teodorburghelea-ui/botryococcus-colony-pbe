# Zhang–Kojima digitization and observable audit (revised 25 September 2026)

## Correction to the earlier audit (22 September and before)

The earlier version of this audit, `methods.tex` and the CSV notes stated that the original figures had been checked and that Fig. 3 contained day-25 histograms for every lighted-volume ratio of the 10-klx preculture. **That was wrong.** In Zhang and Kojima (1998) Fig. 3b (p. 574) the 25-day series (◆) appears only in panel D (V_L = 100 %). The archived file `data/zhang_kojima_1998_fig3_digitized.csv` contains three day-25 histograms (V_L = 5, 25, 50 %) that do not exist in the source. The earlier audit also concluded that Figs. 3 and 5 were separate datasets; that conclusion was also wrong (see below).

All results produced from the archived file are superseded. They are kept, unchanged, in `../zhang_revision_archive_2026-09-22/` for provenance only.

## Current input

`data/zhang_kojima_1998_fig3_redigitized.csv`, produced by scripted marker extraction on the native 300-dpi scan of p. 574 (method, overlays and validation in `data/redigitization_2026-09-25/README.md`). It contains the 25 histograms present in the source:

- 3-klx preculture: days 0, 4 and 16 at V_L = 5, 25, 50 and 100 %;
- 10-klx preculture: days 0, 5 and 18 at all four V_L; day 25 at V_L = 100 % only.

Uncertainty: about ±1 percentage point for clearly visible markers and ±2 points near the baseline; one occluded value (3a-D day 0 at 0.20 mm) is inferred from the identical day-0 histograms of panels A–C. Class boundaries are not given in the source; the primary convention takes plotted marker positions as class representatives, with a −0.011 mm edge shift as a declared sensitivity case.

## Figures 3 and 5 are the same runs

The cubic mean D₃,₀ = (Σ n d³ / Σ n)^{1/3} computed from each re-digitized histogram reproduces the separately digitized Fig. 5 value of the same (preculture, V_L, day) for all 25 runs: RMSE 0.0115 mm, r = 0.996 (0.0055 mm after the −0.011 mm class shift). The archived histograms gave RMSE 0.046 mm. The earlier "invalid pairing" finding was therefore an artefact of the archived digitization, which was broadened and shifted by +0.025 mm. Fig. 5 now serves as a paired scalar validation of the same cultures, including days without a histogram (3 klx: 8, 12, 21; 10 klx: 9, 13, 25).

The Fig. 5 observable is the cubic volume-averaged diameter

\[
D_v=\left[\sum_i n_iD_i^3/\sum_i n_i\right]^{1/3}.
\]

No original numerical records exist in the local package or public searches; all values remain graph extractions.
