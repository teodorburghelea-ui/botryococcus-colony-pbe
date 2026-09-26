# Frozen Zhang-parameter validation against Khatri et al. (2014)

## Source and digitization

The complete Khatri et al. PDF is archived at `source/Khatri_2014.pdf`. Figure 6 is printed page 499. Its horizontal bars encode mean colony diameter on the upper 0.0--0.5 mm axis; dashed diamonds encode optical density on the lower axis. Four bar endpoints were digitized at the displayed 5, 7.5, 10 and 12.5% daily-replacement conditions. The paper describes a 15% condition as close to washout, but Figure 6 does not show a corresponding diameter bar, so it was not assigned a fabricated value. Each endpoint carries an estimated +/-0.002 mm digitization uncertainty from the 300-dpi rendering; the plotted 0.049--0.111 mm spread widths are retained separately and are not treated as digitization errors. The digitized data and pixel calibration are in `data/khatri_2014_fig6_digitized.csv`.

The paper's operation is semi-continuous. For comparison with a continuous dilution rate, the daily replacement fraction $f$ was converted as $D=-\log(1-f)$ per day. This conversion is reported separately from the source's percentage labels.

## Frozen test

The Zhang-calibrated primary instantaneous closures were not refitted. The frozen PBE contains no dilution-dependent rate coefficient, retention law, or inlet distribution. Uniform withdrawal therefore multiplies all bins by the same factor and leaves the normalized diameter distribution unchanged. The reference diameter was computed from the declared Zhang 3-klx, 100%-lighted-volume day-zero histogram solely to expose this structural prediction; it is not claimed as Khatri's inlet distribution.

All three archived primary closures give the same uniform-withdrawal prediction, $D_3=0.2014$ mm. The observed digitized endpoints increase from 0.093+/-0.002 to 0.343+/-0.002 mm. Because Khatri's exact diameter weighting convention is not documented sufficiently to establish a strict correspondence with the model's $D_3$, pointwise RMSE and MAE are retained only in the machine-readable audit as descriptive scale diagnostics and are not interpreted as validation scores. The trend is used only as a qualitative structural boundary test.

## Interpretation

The frozen Zhang closure does not reproduce the Khatri dilution-dependent increase. This is a model-structure result, not evidence that the Zhang parameters are numerically wrong: uniform dilution alone cannot change a normalized distribution. To test Khatri mechanistically, the PBE needs an independently justified dilution-dependent selection, breakup, retention, or inlet-colony distribution. Such an extension must be calibrated or tested independently and must not be presented as a frozen-parameter prediction.
