# Fit a Bayesian Location-Scale Model with Saturated Location and Shared Scale

Fits a brms Gaussian location-scale model where the location uses
saturated cell means (one coefficient per level of `cell_col`) to absorb
ALL systematic variation in the response, and the scale estimates the
power-law relationship between a per-observation precision index and
residual variance, shared across all cells.

## Usage

``` r
fit_saturated_weight(
  df,
  cell_col = "cell",
  concentration_col = "predicted_concentration",
  se_col = "se_concentration",
  pcov_col = "pcov",
  plate_col = "plate",
  predictor = c("auto", "se", "pcov", "se_over_conc"),
  response_is_log10 = FALSE,
  prior_gamma0 = brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma"),
  prior_gamma1 = brms::set_prior("normal(1, 0.5)", class = "b", dpar = "sigma"),
  prior_location = brms::set_prior("normal(0, 2)", class = "b"),
  prior_plate_sd = brms::set_prior("normal(0, 0.5)", class = "sd"),
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  adapt_delta = 0.95,
  seed = 42
)
```

## Arguments

- df:

  Data frame with observation-level data.

- cell_col:

  Character: name of the saturated cell-means factor column. Create this
  externally, e.g.
  `df$cell <- interaction(df$Arm, df$Timeperiod, drop = TRUE)`, or let
  [`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
  build it for you.

- concentration_col:

  Character: predicted concentration column name.

- se_col:

  Character: SE of concentration column name.

- pcov_col:

  Character: posterior CV column name. `NULL` disables it.

- plate_col:

  Character: plate column name for random intercept. `NULL` = no plate
  random effect.

- predictor:

  Scale predictor, passed to
  [`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md).
  One of `"auto"` (default; legacy behaviour), `"se"` (recommended for
  curveR ecosystem input), `"pcov"`, `"se_over_conc"`.

- response_is_log10:

  Logical, passed to
  [`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md):
  is `concentration_col` already on the log10 scale (as in a curveR
  `calibration_result`)? Default `FALSE`.

- prior_gamma0, prior_gamma1, prior_location, prior_plate_sd:

  [`brms::prior`](https://paulbuerkner.com/brms/reference/set_prior.html)
  objects.

- iter, warmup, chains, cores, adapt_delta, seed:

  MCMC controls.

## Value

Named list (see original documentation): `fit`, `phi`, `beta1`,
`gamma_0`, `gamma_1`, credible intervals, `interpretation`,
`effective_se_power`, `data`, `diagnostics`, `cv_diagnostics`,
`formula`, `priors_used`.

## Details

The location model uses `yi ~ 0 + cell [+ (1|plate)]` (saturated cell
means). The scale model uses `log(sigma) = gamma_0 + gamma_1 * log(cv)`.

Default prior on gamma_0 is Normal(0, 1), centering phi = exp(gamma_0)
at 1. Default prior on gamma_1 is Normal(1, 0.5), centering beta1 at the
delta-method prediction of 1, with a 95 percent interval from 0 to 2.

## See also

[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
for the curveR-ecosystem entry point that wraps this estimator,
[`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md),
[`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md),
[`diagnose_cv()`](https://immunoplex.github.io/curveRweights/reference/diagnose_cv.md),
[`weight_diagnostics()`](https://immunoplex.github.io/curveRweights/reference/weight_diagnostics.md).
