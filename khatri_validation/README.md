# Khatri et al. (2014) comparison — status 25 September 2026

**The "frozen-closure prediction" formerly in this package has been retired.** `run_validation.py` did not integrate the PBE: it returned the cubic-mean diameter of a Zhang day-0 histogram (0.2014 mm) for every replacement fraction, and it omitted the biomass → light feedback that the Zhang closure contains. Its code and outputs are kept, unused, in `retired_2026-09-25/` for provenance. They must not be cited as a model prediction.

**Update 26 Sep 2026.** A proper no-refit semi-continuous simulation now exists in `semicontinuous/` (see its README). In it, measured biomass and OD per unit biomass set the flask light, which drives the frozen Zhang closures. The manuscript reports it in the Khatri paragraph and in Supplementary Section S4. It reproduces the direction of the size trend in every run but under-predicts the magnitude, and it does not discriminate the separation closures.

The digitized Fig. 6 bars (`data/`) and the source PDF (`source/`) remain valid.

---

Earlier README (superseded wording retained below for provenance):

# Khatri et al. (2014) validation

This package contains an out-of-sample structural-boundary test of the frozen Zhang-calibrated PBE. The source PDF is stored in `source/Khatri_2014.pdf`. Figure 6 reports colony diameter as horizontal bars and optical density as dashed diamonds; only the bars were digitized for this validation. Each endpoint has an estimated +/-0.002 mm pixel/digitization uncertainty; plotted biological spread widths are stored separately.

The comparison distinguishes daily replacement fraction from continuous dilution rate. Khatri's experiment is semi-continuous, with 5--15% daily replacement; the equivalent continuous dilution rates are $D=-\log(1-f)$ per day, while the paper's Figure 6 labels the replacement fractions in percent.

Frozen primary Zhang parameters are archived in `../zhang_revision/results/primary/`:

- asymmetric binary, instantaneous: `[0.23876300074369664, 0.10873658535099374, 3.0, 0.0]`
- equal binary, instantaneous: `[0.18480533325198165, 0.06553193701003514, 2.4609375, 0.0]`
- broad binary, instantaneous: `[0.19917355582498678, 0.08027126339050926, 2.671875, 0.0]`
