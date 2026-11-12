# Estimate precision weights from a curveR calibration result

The curveR-ecosystem entry point. Two models are available.

## Usage

``` r
fit_precision_weights(
  data,
  design = NULL,
  scale_predictor = c("se", "pcov"),
  plate_col = NULL,
  patient_col = NULL,
  model = c("multiplicative", "additive"),
  use_censoring = TRUE,
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

  Multiplicative model only: `"se"` (default; uncapped log10-scale SD)
  or `"pcov"` (legacy, lossy because of the cap). The additive model
  always uses `se`.

- plate_col:

  Plate grouping column for the location random intercept. Defaults to
  `"curve_id"` when present, else `NULL`.

- patient_col:

  Subject column for the subject random intercept (additive model).
  Defaults to `"patientid"` when present, else `NULL`.

- model:

  `"multiplicative"` (default, 0.2.x behaviour) or `"additive"`
  (recommended; see Description).

- use_censoring:

  Additive model: use the censoring columns (`censor_type`,
  `censor_lower`, `censor_upper`) when the data carry them. Default
  `TRUE`. If `dilution` is also present the response is placed on the
  final-concentration scale, which is the scale of curveRcore's bounds.

- ...:

  MCMC controls and priors passed to
  [`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md)
  or
  [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
  (e.g. `iter`, `warmup`, `chains`, `cores`, `adapt_delta`, `seed`;
  `fix_beta1` for the additive model).

## Value

A `precision_weights` object: a list with `$estimates`, `$weights`,
`$design`, `$scale`, and `$fit`.

For the additive model `$estimates` holds `tau`, `tau2`, `kappa`,
`beta1`, their intervals, the overall `reliability`, and diagnostics;
`$weights` (keyed by `obs_id`/`sampleid`/`curve_id`) has `se`, `pcov`,
`v_meas`, `tau2`, `sigma` (total), `w`, `reliability`,
`reliability_lo/hi`, `w_norm`, `censor_type` and `bound`. See
[`compute_additive_weights()`](https://immunoplex.github.io/curveRweights/reference/compute_additive_weights.md)
for what each column is for. For the multiplicative model `$estimates`
holds `phi`, `beta1` and `$weights` has `se`, `pcov`, `sigma`, `w`,
`w_norm`.

## Details

`model = "multiplicative"` (the default, unchanged from 0.2.x so that
existing callers such as the i-spi-compute worker keep working) fits
\\\sigma_i = \phi\\ se_i^{\beta_1}\\
([`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)).
It attributes *all* variance to the measurement, so very precisely
measured samples can receive extreme weights.

`model = "additive"` (recommended for new analyses) fits the
measurement-error model of
[`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md):
\\\sigma_i^2 = \tau^2 + \kappa^2 se_i^{2\beta_1}\\, which separates true
between-subject variance (\\\tau^2\\) from measurement error and yields
weights \\w_i = 1/\sigma_i^2\\ that are bounded above by \\1/\tau^2\\, a
per-observation `reliability`, and a censored likelihood for samples
outside the quantifiable range. The additive result has a different
`$estimates` shape (`tau`, `kappa`, `reliability` instead of `phi`);
callers must opt in knowingly.

## See also

[`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md),
[`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md),
[`join_weights()`](https://immunoplex.github.io/curveRweights/reference/join_weights.md),
[`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md),
[`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)

## Examples

``` r
# \donttest{
# mp <- curveRfreq::fit_calibration_freq_multiplate(
#         standards = std, samples = samples_with_design, ...)
# wd <- as_weight_data(mp, design = c("timeperiod", "cohort_arm"))
# pw <- fit_precision_weights(wd, iter = 1000, warmup = 500, chains = 2)
# pw$estimates$tau2; pw$estimates$reliability
# }
```
