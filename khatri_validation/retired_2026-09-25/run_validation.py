"""Frozen-parameter Khatri Figure 6 validation.

The Zhang-calibrated PBE has no dilution-rate-dependent coefficient or inlet
distribution. Uniform semi-continuous withdrawal therefore multiplies every
population bin by the same factor and cannot change a normalized diameter
distribution. This script documents that no-selection prediction and compares
it with the digitized Khatri trend without refitting.
"""
from pathlib import Path
import csv, json, sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "zhang_revision"))
from analysis import Data, Solver

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader((HERE / "data/khatri_2014_fig6_digitized.csv").open()))
obs = np.array([float(r["colony_diameter_mm"]) for r in rows])
fractions = np.array([float(r["daily_replacement_fraction"]) for r in rows])

# A declared reference initial distribution: 3-klx, V_L=100% Zhang histogram.
# It is used only to expose the frozen model's dilution invariance; it is not
# claimed to be Khatri's unreported inlet distribution.
data = Data("primary")
solver = Solver("asymmetric_binary")
n0 = solver.initialize(data.hist[3, 100, 0.0])
reference_d3 = float(np.cbrt((solver.x @ n0) / n0.sum()))

# Uniform withdrawal over one daily replacement interval. The normalized
# distribution and D3 are exactly invariant under n -> (1-f)n.
pred = np.full(len(rows), reference_d3)
residual = pred - obs
rmse = float(np.sqrt(np.mean(residual**2)))
mae = float(np.mean(np.abs(residual)))
r = float(np.corrcoef(fractions, obs)[0, 1])
slope = float(np.polyfit(fractions, obs, 1)[0])

out = HERE / "results"
out.mkdir(exist_ok=True)
with (out / "khatri_frozen_predictions.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["dilution_label_percent", "daily_replacement_fraction", "observed_colony_diameter_mm", "frozen_prediction_mm", "residual_mm"])
    for row, p, e in zip(rows, pred, residual):
        w.writerow([row["dilution_label_percent"], row["daily_replacement_fraction"], row["colony_diameter_mm"], f"{p:.9f}", f"{e:.9f}"])

closures = {
    "asymmetric_binary_instantaneous": "primary/asymmetric_binary_instantaneous.json",
    "equal_binary_instantaneous": "primary/equal_binary_instantaneous.json",
    "broad_binary_instantaneous": "primary/broad_binary_instantaneous.json",
}
closure_rows = []
for name, rel in closures.items():
    params = json.loads((ROOT / "zhang_revision/results" / rel).read_text())["parameters"]
    closure_rows.append({"closure": name, "parameters": params, "frozen_prediction_D3_mm": reference_d3, "rmse_mm": rmse, "mae_mm": mae})
with (out / "khatri_frozen_closure_metrics.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["closure", "parameters", "frozen_prediction_D3_mm", "rmse_mm", "mae_mm"])
    w.writeheader()
    for row in closure_rows:
        row["parameters"] = json.dumps(row["parameters"])
        w.writerow(row)

summary = {
    "source": "Khatri et al. 2014 Figure 6",
    "frozen_parameters": "Zhang primary asymmetric-binary instantaneous closure",
    "reference_initialization": "Zhang 3-klx V_L=100% day-zero histogram; diagnostic only",
    "prediction_definition": "uniform semi-continuous withdrawal with no D-dependent rate or inlet-state closure",
    "reference_predicted_D3_mm": reference_d3,
    "observed_range_mm": [float(obs.min()), float(obs.max())],
    "rmse_mm": rmse,
    "mae_mm": mae,
    "observed_diameter_slope_mm_per_fraction": slope,
    "observed_diameter_fraction_correlation": r,
    "closures_tested": [r["closure"] for r in closure_rows],
    "conclusion": "Frozen Zhang closure cannot reproduce the dilution-dependent colony-size increase under its stated uniform-withdrawal assumptions; adding a D-dependent selection, breakup, retention, or inlet-distribution mechanism would require new calibration or independent data.",
}
(out / "khatri_frozen_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
