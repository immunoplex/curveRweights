# Fit the additive measurement-error model (signal + measurement variance)

Fits a Bayesian Gaussian model in which each observation's variance is
the *sum* of true between-subject variance and measurement variance,
\$\$\sigma_i^2 = \tau^2 + \kappa^2\\ se_i^{2\beta_1},\$\$ with a
censored likelihood for samples outside the quantifiable range. It is
the default estimator of
[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
from curveRweights 0.3.0.

## Usage

``` r
fit_additive_weight(
  df,
  cell_col = "cell",
  concentration_col = "predicted_concentration",
  se_col = "se_concentration",
  plate_col = "plate",
  patient_col = NULL,
  sample_col = NULL,
  censor_col = NULL,
  lower_col = "censor_lower",
  upper_col = "censor_upper",
  se_bound_col = NULL,
  dilution_col = NULL,
  response_is_log10 = FALSE,
  fix_beta1 = TRUE,
  prior_location = brms::set_prior("normal(0, 2)", class = "b"),
  prior_ltau = brms::set_prior("normal(-1.2, 1)", class = "b", nlpar = "ltau"),
  prior_lkappa = brms::set_prior("normal(0, 0.5)", class = "b", nlpar = "lkappa"),
  prior_beta1 = brms::set_prior("normal(1, 0.25)", class = "b", nlpar = "beta1"),
  prior_group_sd = brms::set_prior("normal(0, 0.5)", class = "sd"),
  iter = 2000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  adapt_delta = 0.95,
  seed = 42,
  backend = NULL,
  n_draws_weights = 1000
)
```

## Arguments

- df:

  Data frame with observation-level data.

- cell_col:

  Name of the saturated cell factor column.

- concentration_col:

  Predicted concentration column.

- se_col:

  Measurement SE column (log10 scale; the curveR `se_concentration`).

- plate_col:

  Plate (or curve) grouping column; `NULL` for none.

- patient_col:

  Subject grouping column; `NULL` for none. Strongly recommended when
  subjects are measured more than once.

- sample_col:

  Optional column identifying a physical sample measured several times
  (dilutions/replicates).

- censor_col:

  Name of the censoring-type column (`"none"`, `"left"`, `"right"`,
  `"interval"`). `NULL` (default) treats every row as quantified.

- lower_col, upper_col:

  Natural-scale censoring bounds.

- se_bound_col:

  Optional measurement SE at the censoring limit; see
  [`prepare_censoring()`](https://immunoplex.github.io/curveRweights/reference/prepare_censoring.md).

- dilution_col:

  Optional dilution-factor column. If supplied, the response is
  converted to the final-concentration scale (`yi + log10(dilution)`),
  which is the scale curveRcore's censoring bounds are on. Supply it
  whenever you use curveRcore censoring bounds and dilutions differ.

- response_is_log10:

  Is `concentration_col` already log10?

- fix_beta1:

  Fix \\\beta_1 = 1\\ (default) or estimate it.

- prior_location, prior_ltau, prior_lkappa, prior_beta1, prior_group_sd:

  [`brms::prior`](https://paulbuerkner.com/brms/reference/set_prior.html)
  objects.

- iter, warmup, chains, cores, adapt_delta, seed:

  MCMC controls.

- backend:

  Passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html);
  `NULL` uses the brms default.

- n_draws_weights:

  Number of posterior draws used for the reliability interval (thinned
  from the fit).

## Value

A list with `fit`, `tau`, `tau2`, `kappa`, `beta1`, their credible
intervals (`tau_CI`, `kappa_CI`, `beta1_CI`), overall `reliability` and
`reliability_CI` (\\\tau^2/(\tau^2 + \bar v)\\ over quantified rows),
`data` (the full input with the weight columns of
[`compute_additive_weights()`](https://immunoplex.github.io/curveRweights/reference/compute_additive_weights.md)),
`diagnostics`, `formula` and `priors_used`.

## In plain terms

An assay value varies between subjects for two different reasons.
Subjects truly differ (biology; this is the signal, variance
\\\tau^2\\), and the calibration curve measures each sample with some
error (variance \\v_i\\, known to us through `se_i`). The two add. The
older multiplicative model of curveRweights 0.2.x instead made the
*whole* spread proportional to `se`, which implies that a precisely
measured sample has little biological variation, so it could be given a
weight many times larger than its neighbours. In the additive model the
weight is at most \\1/\tau^2\\; no sample can count for more than a
perfectly measured one.

\\\kappa\\ is a calibration factor: \\\kappa = 1\\ means the reported
`se` is exactly the measurement SD (the delta-method value). \\\beta_1\\
is fixed at 1 by default; \\\kappa\\ and \\\tau\\ are separated only
through the spread of `se` across samples, so freeing \\\beta_1\\ as
well makes the model much harder to identify. Use `fix_beta1 = FALSE`
only with plenty of data over a wide range of `se`.

## Model

Location: `yi ~ 0 + cell + (1|patient) + (1|plate) [+ (1|sample)]`,
where `cell` is the saturated design factor (e.g. arm by time point).
The patient term absorbs correlation between repeated samples of the
same subject (otherwise counted as noise); the optional sample term does
the same for several dilutions or replicates of one physical sample.

Scale: \\\sigma_i = \sqrt{\tau^2 + \kappa^2 se_i^{2\beta_1}}\\ via a
non-linear `sigma` formula. Rows with `censor_type` `left`, `right` or
`interval` contribute a cumulative-probability term at their censoring
limit (see
[`prepare_censoring()`](https://immunoplex.github.io/curveRweights/reference/prepare_censoring.md));
quantified rows contribute a normal density.

Default priors: location `normal(0, 2)`; \\\log\tau \sim N(\log 0.3,
1)\\; \\\log\kappa \sim N(0, 0.5)\\ (centred on the delta-method value);
\\\beta_1 \sim N(1, 0.25)\\ when free; group SDs `normal(0, 0.5)`.

## See also

[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md),
[`compute_additive_weights()`](https://immunoplex.github.io/curveRweights/reference/compute_additive_weights.md),
[`prepare_censoring()`](https://immunoplex.github.io/curveRweights/reference/prepare_censoring.md)
