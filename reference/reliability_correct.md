# Regression calibration: shrink measured values toward their cell mean

Replaces each measurement \\y_i\\ by its best prediction of the true
value, \\\mu\_{cell} + \lambda_i (y_i - \mu\_{cell})\\. Use this when
the measured concentration is a *predictor* (correlate-of-protection
regressions, antigen-versus-antigen correlations).

## Usage

``` r
reliability_correct(y, cell, reliability)
```

## Arguments

- y:

  Numeric vector of measurements (log10 scale).

- cell:

  Grouping vector (the design cell, e.g. arm by time point); the
  shrinkage target is the reliability-weighted mean of each cell.

- reliability:

  Numeric vector of per-observation reliabilities in (0, 1\].

## Value

Numeric vector of calibrated values, same length as `y`.

## In plain terms

Weights do not repair measurement error in an explanatory variable.
Noise in X makes a regression slope too small (attenuation, or
"regression dilution"). A reliability \\\lambda\\ of 0.8 means roughly
80\\ see is real and the slope shrinks by that factor. Replacing X by
the shrunken value above removes that bias to first order. Also consider
[`brms::me()`](https://paulbuerkner.com/brms/reference/me.html) for a
fully model-based approach.

## Examples

``` r
set.seed(1)
cell <- rep(c("a", "b"), each = 50)
rel  <- runif(100, 0.4, 0.95)
y    <- rnorm(100, ifelse(cell == "a", 1, 2))
head(reliability_correct(y, cell, rel))
#> [1] 1.2933132037 0.6960542454 1.2915845496 0.0009299955 1.8139752264
#> [6] 2.7884154050
```
