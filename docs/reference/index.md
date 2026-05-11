# Package index

## Data Preparation

Prepare and diagnose calibration curve precision data

- [`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md)
  : Prepare the Precision Index from Calibration Curve Output
- [`diagnose_cv()`](https://immunoplex.github.io/curveRweights/reference/diagnose_cv.md)
  : Diagnose Precision Index Variation for Scale Estimation

## Model Fitting

Fit the Bayesian location-scale model

- [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
  : Fit a Bayesian Location-Scale Model with Saturated Location and
  Shared Scale
- [`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md)
  : Fit Saturated Weight Models Across Multiple Groups

## Weight Computation

Compute and apply precision weights

- [`compute_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/compute_saturated_weights.md)
  : Compute Precision Weights from Estimated Scale Parameters
- [`apply_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md)
  : Apply Previously Estimated Scale Parameters to New Data

## Diagnostics

Interpret and diagnose weight estimates

- [`interpret_beta1()`](https://immunoplex.github.io/curveRweights/reference/interpret_beta1.md)
  : Interpret the Estimated Beta1 Value
- [`weight_diagnostics()`](https://immunoplex.github.io/curveRweights/reference/weight_diagnostics.md)
  : Compute Summary Diagnostics for Precision Weights

## Data

Example datasets

- [`example_assay`](https://immunoplex.github.io/curveRweights/reference/example_assay.md)
  : Example Luminex Multiplex Immunoassay Data
