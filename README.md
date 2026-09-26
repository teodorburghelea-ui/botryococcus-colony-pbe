# Population-balance tests of colony growth, separation and light access in *Botryococcus braunii*: reproducibility archive

This archive holds the code, data read from published figures, fitted parameters, search traces, per-histogram predictions and numerical checks for the article:

> T. Burghelea, *Testing population-balance descriptions of colony growth, separation and light access in* Botryococcus braunii (submitted to *Algal Research*).

Nothing here is new experimental data. Every quantitative target was extracted from a published figure or table, and each folder's README documents how.

**Not redistributed.** Publisher PDFs, figure images and page renders are copyrighted and are not included. The extraction scripts rebuild the needed images from a reader's own copy of each paper (the sources are listed below).

## Where each result comes from

| Manuscript item | Folder | Entry point |
|---|---|---|
| Zhang & Kojima (1998) bubble-column transfer test: Zhang comparison table, extraction-sensitivity envelope, Fig. `matched_zhang`, Suppl. S2 | `zhang_revision/` | `analysis.py`, `report.py`, `check_numerics.py`, `counting_noise.py`, `verify_results.py` |
| Re-digitization of Zhang & Kojima Fig. 3 (25 histograms) | `zhang_revision/data/redigitization_2026-09-25/` | `scripts/*.py` (see its README) |
| van den Berg et al. (2019) light-acclimation ablation: prediction table, Fig. `matched_ablation`, Suppl. S1 | top level (`run_ablation.py` …) and `data/`, `results/` | `test_ablation.py`, `run_ablation.py`, `polish_fits.py`, `report_ablation.py`, `verify_outputs.py` |
| Point re-extraction of van den Berg Fig. 2 | `data/vdb_reextraction_2026-09-25/` | `extract_vdb_fig2.py` (see its README) |
| Khatri et al. (2014) semi-continuous simulation, Suppl. S4 | `khatri_validation/semicontinuous/` | `run_semicontinuous.py`, `tau_scan.py`, `plot_semicontinuous.py` |
| García-Cubero et al. (2021) chemostat regressions (`fig10_garcia_chemostat_state_test`) | `supporting_analyses/` | `garcia_cubero_chemostat_point4.py`, `plot_garcia_chemostat.py` |
| Kemel et al. (2025) optical/productivity endpoints (`fig4_kemel_productivity`, `fig13_optics_to_productivity_test`) | `supporting_analyses/` | `plot_kemel_productivity_endpoints.py`, `optics_to_productivity_point7.py`, `plot_optics_to_productivity.py` |
| Microscopy-constrained benchmark and stress ramp (`fig11_microscopy_model_selection`, `fig12_biological_mechanical_stress_ramp`) | `supporting_analyses/` | `microscopy_model_selection_point5.py`, `biological_mechanical_stress_ramp_point6.py` and their `plot_*.py` |
| Structured-PBE numerical integrity (`fig6_numerical_integrity_convergence`) | `supporting_analyses/` | `structured_pbe_integrity.jl` (Julia), `plot_numerical_integrity.py` |
| Final figure files | `figures/` | `fig1_colony_architecture` is a drawn schematic with no generating script |

More notes:
- The **light-acclimation README** is `README_light_acclimation.md`, and the methods text for that analysis is `methods.tex`.
- **Superseded analyses.** `supporting_analyses/zhang_hierarchy_point2.py` and `memory_bimodality_point3.py` belong to earlier drafts and are no longer used by the manuscript. `zhang_revision/data/zhang_kojima_1998_fig3_digitized.csv` is superseded as well; see `data/SUPERSEDED_*.txt`.
- **Figure output folder.** The `plot_*.py` scripts in `supporting_analyses/` write, by default, to a `Draft_2026-09-07/figures/` folder of the original project. Pass `--output-dir <folder>` to write elsewhere.

## Requirements

- Python ≥ 3.10 with NumPy, SciPy and Matplotlib, plus Pillow for the extraction scripts.
- Julia (standard library only), needed only for `structured_pbe_integrity.jl`.
- poppler-utils (`pdfimages`), needed only to rebuild the Zhang page scan.

**Run times on one core:**

| Step | Time |
|---|---|
| Zhang: all 10 scenarios (90 fits) | about 1 h |
| Light-acclimation ablation | about 8 min |
| Khatri simulation | about 15 min |

Set `OPENBLAS_NUM_THREADS=1` for bitwise-stable results.

## Sources of the extracted data

| Source | Figures used |
|---|---|
| Zhang, K., Kojima, E. (1998) *J. Ferment. Bioeng.* 86:573–576. doi:10.1016/S0922-338X(99)80009-9 | Figs. 3 and 5, Table 1 |
| van den Berg, T.E. et al. (2019) *Plant Physiol.* 179:1132–1143. doi:10.1104/pp.18.01499 | Fig. 2 (figure file at PubMed Central, PMC6393799 `PP_201801499DR1_f2.jpg`) |
| Khatri, W. et al. (2014) *Biotechnol. Bioeng.* 111:493–503. doi:10.1002/bit.25126 | Figs. 5 and 6 |
| García-Cubero, R. et al. (2021) *Bioresour. Technol.* 340:125653. doi:10.1016/j.biortech.2021.125653 | Tables 2–4 |
| Kemel, S. et al. (2025) *Algal Res.* 90:104176. doi:10.1016/j.algal.2025.104176 | abstract values |

## Licence

- **Code:** MIT (`LICENSE`).
- **Data and results:** extracted values, CSV/JSON outputs and generated figures are under CC-BY-4.0 (`LICENSE-DATA`).
- **Underlying measurements:** these belong to the cited authors; cite the original articles when reusing them.

## How to cite

See `CITATION.cff`. The archive DOI will be added on deposit.
