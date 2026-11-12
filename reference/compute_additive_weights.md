# Compute additive-model weights, reliability and measurement variance

Given the estimated signal variance \\\tau^2\\ and measurement-error
scale \\\kappa\\ (and optionally exponent \\\beta_1\\), computes for
every row the quantities a downstream analysis needs.

## Usage

``` r
compute_additive_weights(
  df,
  tau2,
  kappa = 1,
  beta1 = 1,
  se_col = "se_meas",
  censor_col = NULL,
  draws = NULL,
  level = 0.95
)
```

## Arguments

- df:

  Data frame with the measurement SE column.

- tau2:

  Signal variance \\\tau^2\\ (point estimate).

- kappa:

  Measurement-error scale \\\kappa\\. Default 1 (the delta-method value:
  `se` is the true measurement SD).

- beta1:

  Exponent on `se`. Default 1.

- se_col:

  Name of the measurement SE column. Default `"se_meas"`.

- censor_col:

  Name of the censoring-type column (optional).

- draws:

  Optional data frame of posterior draws with columns `tau2`, `kappa`
  and `beta1`; used only for the reliability interval.

- level:

  Credible level for the reliability interval.

## Value

`df` with the columns above added.

## In plain terms

Every assay number is "true value + biology + measurement error". The
calibration curve tells us the measurement error of each sample (`se`);
the model estimates how much true between-subject variation (\\\tau^2\\)
sits on top of it. Because the two *add*, the total variance of a sample
is \\\tau^2 + v_i\\ and its inverse-variance weight \\w_i = 1/(\tau^2 +
v_i)\\ can never exceed \\1/\tau^2\\: a perfectly measured sample counts
at most as much as the biology allows.

The **reliability** \\\lambda_i = \tau^2/(\tau^2 + v_i)\\ is the share
of a sample's variance that is signal rather than measurement error. It
runs from 0 (pure noise) to 1 (perfectly measured) and reads like an
intraclass correlation.

## Output columns

- v_meas:

  Measurement variance \\v_i = \kappa^2\\ se_i^{2\beta_1}\\.

- tau2:

  Signal variance (the same for every row).

- sigma_total:

  \\\sqrt{\tau^2 + v_i}\\.

- w:

  Raw inverse-variance weight \\1/\sigma\_{total}^2\\. Use for
  estimators where the absolute variance matters (e.g. `metafor`'s `vi`
  is `v_meas`, not `1/w`).

- reliability:

  \\\tau^2 w\\, in (0, 1\].

- w_norm:

  `w` divided by its mean over the eligible rows (quantified rows), so
  the mean is 1. Use for software that treats weights as case counts,
  such as `brms`'s [`weights()`](https://rdrr.io/r/stats/weights.html).

- reliability_lo, reliability_hi:

  Posterior interval of the reliability, present when `draws` is
  supplied.

- censor_type, bound:

  Carried through when present. Censored rows get `w_norm = NA`:
  represent them by their censoring bound, not a weight.

## See also

[`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md),
[`normalize_weights()`](https://immunoplex.github.io/curveRweights/reference/normalize_weights.md)

## Examples

``` r
d <- data.frame(se_meas = c(0.03, 0.1, 0.3, 1))
compute_additive_weights(d, tau2 = 0.16, kappa = 1)
#>   se_meas v_meas tau2 sigma_total        w reliability    w_norm
#> 1    0.03  9e-04 0.16   0.4011234 6.215040   0.9944065 1.4658579
#> 2    0.10  1e-02 0.16   0.4123106 5.882353   0.9411765 1.3873914
#> 3    0.30  9e-02 0.16   0.5000000 4.000000   0.6400000 0.9434261
#> 4    1.00  1e+00 0.16   1.0770330 0.862069   0.1379310 0.2033246
```
