# Estimate precision weights from a curveR calibration result

The curveR-ecosystem entry point. Estimates the power-law relationship
between calibration-curve precision and residual variance, \\\sigma_i =
\phi\\ se_i^{\beta_1}\\, \\w_i = 1/\sigma_i^2\\, using a joint Bayesian
location-scale model with a saturated cell-means location (see
[`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)).
The default scale predictor is the uncapped `se_concentration` from the
calibration curve, so \\\phi = 1\\ means "se_concentration is a
calibrated residual SD".

## Usage

``` r
fit_precision_weights(
  data,
  design = NULL,
  scale_predictor = c("se", "pcov"),
  plate_col = NULL,
  ...
)
```

## Arguments

- data:

  Either the output of
  [`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
  (class `weight_data`), or a `calibration_result(_multiplate)` / plain
  `data.frame`, in which case `design` is required and
  [`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
  is called internally.

- design:

  Character vector of design-group columns (required unless `data` is
  already a `weight_data`).

- scale_predictor:

  `"se"` (default; uncapped log10-scale SD) or `"pcov"` (legacy, lossy
  because of the cap).

- plate_col:

  Plate grouping column for the location random intercept. Defaults to
  `"curve_id"` when present, else `NULL`.

- ...:

  MCMC controls and priors passed to
  [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
  (e.g. `iter`, `warmup`, `chains`, `cores`, `adapt_delta`, `seed`).

## Value

A `precision_weights` object: a list with `$estimates` (phi, beta1, CIs,
interpretation, diagnostics), `$weights` (data frame keyed by
`obs_id`/`sampleid`/`curve_id` with `se`, `pcov`, `sigma`, `w`,
`w_norm`), `$design`, `$scale` (conc_scale + scale_predictor), and
`$fit` (the brms fit).

## See also

[`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md),
[`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md),
[`join_weights()`](https://immunoplex.github.io/curveRweights/reference/join_weights.md),
[`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)

## Examples

``` r
# \donttest{
# mp <- curveRfreq::fit_calibration_freq_multiplate(
#         standards = std, samples = samples_with_design, ...)
# wd <- as_weight_data(mp, design = c("timeperiod", "cohort_arm"))
# pw <- fit_precision_weights(wd, iter = 1000, warmup = 500, chains = 2)
# pw$estimates$phi; pw$estimates$beta1
# }
```
