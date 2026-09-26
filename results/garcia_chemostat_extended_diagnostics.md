# Continuous-culture comparator diagnostics

The nine-row calibration subset was prespecified from the process design before model fitting or inspection of model performance: three centre replicates plus paired high/low contrasts for temperature (2/8), incident light (3/11), and dilution (5/7). The four remaining nonzero Table-2 rows and all four optimized-condition targets were withheld.

The constant arithmetic-mean baseline (calibration mean 172.9 micrometres) has calibration RMSE 59.9 micrometres, unused-design-check RMSE 74.6 micrometres, and optimized-condition RMSE 79.2 micrometres.
Leave-one-condition-out calibration RMSE is 57.0 micrometres for the environmental-state closure, 71.0 micrometres for the dilution-free null, and 67.3 micrometres for the constant-mean baseline.
Environmental-state calibration leverage ranges from 0.167 to 0.667; the largest leverage is table2_experiment_8 and the largest Cook-style influence is table2_experiment_8 (D=1.560). These diagnostics flag influential design points; they do not establish independent biological replicates.
The 95% predictive intervals below are log-scale OLS intervals including residual variance and are diagnostic, not confidence intervals for a universal response surface.
- single_colony_size_maximum: 340.5 micrometres [160.4, 723.1].
- triple_response_optimum: 157.0 micrometres [79.8, 308.5].
- biomass_colony_size_optimum: 196.6 micrometres [96.4, 401.0].
- hydrocarbon_colony_size_optimum: 199.4 micrometres [107.4, 370.1].
