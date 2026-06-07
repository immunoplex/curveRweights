# Package index

## curveR Ecosystem Entry Points

The recommended API for users coming from curveRfreq or curveRbayes.
[`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
converts a `calibration_result(_multiplate)` into the standardised input
frame.
[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
estimates phi and beta1.
[`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
applies fitted weights to new data.
[`join_weights()`](https://immunoplex.github.io/curveRweights/reference/join_weights.md)
attaches weights back onto a sample data frame.

- [`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
  : Build precision-weight input from a curveR calibration result
- [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
  : Estimate precision weights from a curveR calibration result
- [`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
  : Apply fitted precision weights to new observations
- [`join_weights()`](https://immunoplex.github.io/curveRweights/reference/join_weights.md)
  : Join estimated weights back onto a data frame

## Precision Index Preparation

Functions that compute and diagnose the precision index used as the
scale predictor.
[`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md)
builds yi, cv_i, and log_cv from concentration and precision columns.
[`diagnose_cv()`](https://immunoplex.github.io/curveRweights/reference/diagnose_cv.md)
checks whether the precision index has enough variation to identify
beta1.
[`interpret_beta1()`](https://immunoplex.github.io/curveRweights/reference/interpret_beta1.md)
classifies the estimated exponent into human-readable regimes.

- [`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md)
  : Prepare the Precision Index from Calibration Curve Output
- [`diagnose_cv()`](https://immunoplex.github.io/curveRweights/reference/diagnose_cv.md)
  : Diagnose Precision Index Variation for Scale Estimation
- [`interpret_beta1()`](https://immunoplex.github.io/curveRweights/reference/interpret_beta1.md)
  : Interpret the Estimated Beta1 Value

## Weight Computation

Deterministic weight computation from estimated scale parameters and
summary diagnostics for the resulting weight distribution.

- [`compute_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/compute_saturated_weights.md)
  : Compute Precision Weights from Estimated Scale Parameters
- [`weight_diagnostics()`](https://immunoplex.github.io/curveRweights/reference/weight_diagnostics.md)
  : Compute Summary Diagnostics for Precision Weights

## Core Fitting Engine

Lower-level functions that fit the joint Bayesian location-scale model
via brms.
[`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
fits one group;
[`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md)
loops over multiple groups.
[`apply_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md)
applies a saved scale_table to new data without refitting.

- [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
  : Fit a Bayesian Location-Scale Model with Saturated Location and
  Shared Scale
- [`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md)
  : Fit Saturated Weight Models Across Multiple Groups
- [`apply_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md)
  : Apply Previously Estimated Scale Parameters to New Data

## Example Dataset

Anonymised Luminex multiplex immunoassay data for testing and
documentation.

- [`example_assay`](https://immunoplex.github.io/curveRweights/reference/example_assay.md)
  : Example Luminex Multiplex Immunoassay Data
