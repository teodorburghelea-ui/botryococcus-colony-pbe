# Low-light memory and bimodality model-selection test

Both PBE variants start from the same reconstructed medium-light volume distribution, converted once to a number distribution. The mean of the two day-20 low-light curves is read only after forward simulation for diagnostics; this script contains no optimizer and fits no mixture weights. Continuous material transport and shell loss are zero in both variants, so this is a topology/model-selection test rather than a biomass-growth fit.

## Fixed reference closure

The state-collapsed one-coordinate limit has a common weak equal-volume mechanical rate of 0.020 d^-1. The structured calculation uses tau_a=3.20 d, a transition gate from a=0.22 to 0.72, and a material-conserving mother-plus-daughter event. Each event transfers f=0.10 of affected parent material to one broad daughter material kernel centered at 76.5 um (0.18 times the 425.2-um medium-light modal diameter), while retaining the mother branch.

## Declared topology criterion

A pass requires a broad local maximum in 45--180 um, another in 250--850 um, a diameter ratio of at least 2.5, and at least 10% of volume in both the <200-um small-colony and >=300-um large-residual regions. The low-light reconstruction has small and large mode locations of 90.4 and 425.2 um, with fractions 0.406 and 0.472.

one_state_collapsed_limit: TV=0.246; W1=0.104 mm; small fraction=0.173; large residual=0.639; small-mode peak=nan um; large-mode peak=320.9 um; topology pass=False; material drift=0.00e+00.
structured_memory_biological_kernel: TV=0.114; W1=0.049 mm; small fraction=0.378; large residual=0.436; small-mode peak=78.5 um; large-mode peak=278.7 um; topology pass=True; material drift=0.00e+00.

Fixed sensitivity set: 5/6 structured calculations pass the same topology criterion for tau_a in {2.4, 3.2, 4.0} d and transferred material fraction in {0.06, 0.10}. These are closure sensitivities, not refits.

## Interpretation

The state-collapsed one-coordinate limit retains a large population but has no mechanism to create a separate small-colony mode after the low-light shift. The structured closure creates a transient small daughter branch while retaining a large mother branch, and does so without inserting a day-20 mixture into the initial condition. This discriminates the model structure; it does not identify the acclimation time, daughter number, shell loss, or biological event rate. The later reported collapse to a low-light unimodal distribution is outside this 20-d test and remains an independent constraint.
