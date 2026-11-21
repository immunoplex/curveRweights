# Changelog

## curveRweights 0.3.1

### Additive model: options

- `fix_kappa` (in
  [`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md),
  passed on by
  [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)):
  fix kappa = 1, i.e. take the calibration curve’s `se` as the
  measurement SD and estimate only tau. For data whose `se` barely
  varies, where kappa is weakly identified. The identifiability note
  says that kappa was fixed.
- [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
  accepts `patient_col = NA` to leave out the subject random intercept.
- The additive result’s `$estimates$options` records the options
  actually applied (`use_censoring`, `patient_col`, `plate_col`,
  `fix_beta1`, `fix_kappa`). The i-spi-compute worker stores them in
  `calib_weights_fit.fit_options` (curveRcore \>= 0.7.1).

### The same diagnostics for both models

Both fits are Bayesian location-scale models sampled with Stan, so both
now report the same entries in `$diagnostics` (and
`$estimates$diagnostics` from
[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)):

- `rhat_max`, `n_divergent` (already there) and new `ess_min`, the
  lowest bulk or tail effective sample size over the sampled parameters
  (via the posterior package, now in Imports).
- `ident_sd_ratio`: posterior SD / prior SD of the key scale parameter.
  For the additive model it is the existing kappa ratio (now read from
  the `prior_lkappa` actually used rather than assumed 0.5); for the
  multiplicative model it is new, for beta1 (NA when the slope is
  dropped).
- `resid_sd`, `resid_sd_ratio`: a variance-model check. The scaled
  residuals (y - fitted mean) / sigma of the quantified rows should have
  SD 1 everywhere; `resid_sd` is their robust SD, `resid_sd_ratio` the
  robust SD in the lowest-se fifth over that in the highest-se fifth
  (above 1: low-se samples are given too much weight).

The weights themselves are unchanged: the diagnostics are computed after
the fit (the multiplicative weights are identical to 0.3.1 before this
change for the same seed).

### Bug fix

- [`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md)
  no longer fails when the model has no group term at all (no patient
  term and a single plate):
  [`brms::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) errors
  on such a model and is now only called when there are group terms.

### Known limitation

- `fit_additive_weight(sample_col = )` (several readings of one sample)
  is not recommended. With the sample intercept in the location, the
  sample term absorbs the real between-sample variation and tau is left
  to describe only the extra scatter between readings of one sample; tau
  then tends to zero, the ceiling 1/tau^2 on the weights is lost, and
  sampling is poor (high Rhat, treedepth warnings in a two-dilution
  simulation). A formulation in which the sample term *is* the signal
  (residual = measurement error only) is planned.
  [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
  does not pass it on.

### Documentation

- Measurement-error vignette: with the subject random intercept in the
  model, tau is the variation *within* a person and reliability is a
  within-person reliability; reliability belongs to the calibration
  method that produced the `se` (Bayesian and frequentist `se` can split
  the variance very differently); how to read kappa far from 1;
  `fix_kappa`; a new section on checking a fit with the diagnostics both
  models share.

### Tests

- Shape-only tests silence Stan’s low-ESS warnings (`quiet_ess()`
  helper); the censored-likelihood bias test runs longer chains. The
  suite runs without warnings.

## curveRweights 0.3.0

### Additive measurement-error model (new, opt-in)

- [`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md)
  fits a model in which each observation’s variance is the *sum* of true
  between-subject variance and measurement variance, sigma^2 = tau^2 +
  kappa^2 \* se^(2 \* beta1). The weight 1/sigma^2 can never exceed
  1/tau^2, so a very precisely measured sample can no longer dominate
  (the multiplicative model of 0.2.x attributes all variance to the
  measurement). beta1 is fixed at 1 by default (`fix_beta1 = TRUE`)
  because tau and kappa are separated only through the spread of `se`; a
  message reports when kappa is weakly identified.
- Samples outside the quantifiable range enter through a censored
  likelihood (left, right or interval) instead of being treated as exact
  points or weighted:
  [`prepare_censoring()`](https://immunoplex.github.io/curveRweights/reference/prepare_censoring.md)
  builds the censoring flags and the measurement SE at the limit from
  curveRcore’s `censor_type` / `censor_lower` / `censor_upper`. A
  `(1 | patient)` term absorbs repeated measures of one subject;
  optional `sample_col` handles dilutions or replicates of one sample;
  optional `dilution_col` puts the response on the same
  (final-concentration) scale as curveRcore’s bounds.
- [`compute_additive_weights()`](https://immunoplex.github.io/curveRweights/reference/compute_additive_weights.md)
  returns the weight table a downstream analysis needs: `v_meas`,
  `tau2`, `sigma_total`, `w`, `reliability` (tau^2 / (tau^2 + v), with a
  posterior interval), `w_norm` (mean 1 over quantified rows),
  `censor_type`, `bound`.
- New helpers:
  [`normalize_weights()`](https://immunoplex.github.io/curveRweights/reference/normalize_weights.md)
  (re-normalize to mean 1 within an analysis subset),
  [`meas_variance()`](https://immunoplex.github.io/curveRweights/reference/meas_variance.md)
  (the variance to pass to `metafor::rma.mv(V = )`),
  [`reliability_correct()`](https://immunoplex.github.io/curveRweights/reference/reliability_correct.md)
  (regression calibration) and
  [`attenuation_correct()`](https://immunoplex.github.io/curveRweights/reference/attenuation_correct.md).

### Documentation (0.3.0)

- New vignette `measurement-error`: the additive model, reliability,
  censoring, attenuation and a glossary, written for readers without a
  measurement-error background.
- New vignette `downstream-analysis`: recipes by question
  ([`metafor::rma.mv`](https://wviechtb.github.io/metafor/reference/rma.mv.html)
  with
  [`meas_variance()`](https://immunoplex.github.io/curveRweights/reference/meas_variance.md),
  censored `brms`, regression calibration, rank tests with censoring).
  The `method = "FE"` and `svyranktest()` advice in
  `precision-weighting` is replaced.
- `AI_USE.md`: how generative AI was used in this project.

### Existing behaviour is unchanged

- [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
  gains `model = c("multiplicative", "additive")`, `patient_col` and
  `use_censoring`. The default stays `"multiplicative"`, so existing
  callers (including the i-spi-compute worker) get exactly the 0.2.x
  result and shape. The `precision_weights` object now records
  `$estimates$model`.
- The defaults of
  [`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md),
  [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
  and
  [`apply_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md)
  (`predictor = "auto"`) are unchanged;
  [`apply_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md)
  gains a `predictor` argument.
  [`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
  and
  [`join_weights()`](https://immunoplex.github.io/curveRweights/reference/join_weights.md)
  understand both models.
- [`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
  gains `carry` and now carries `patientid`, `dilution` and the
  censoring columns through when present.
- `test-backward-compat.R` pins these defaults and the result shape the
  worker reads.

## curveRweights 0.2.2

- [`vignette("precision-weighting")`](https://immunoplex.github.io/curveRweights/articles/precision-weighting.md)
  gains a “Weights and censored samples” section: why weighting
  out-of-range samples biases location estimates, how to apply weights
  only to quantified samples and represent the rest by their censoring
  bounds (`classify_censoring_multiplate()` in curveRcore), and a note
  that censoring-aware weight estimation is planned, not yet
  implemented.

## curveRweights 0.2.0

### curveR ecosystem integration

- [`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
  — adapter that converts a `calibration_result` or
  `calibration_result_multiplate` (from curveRfreq or curveRbayes) into
  the standardised `weight_data` frame consumed by the estimator.
  Dispatches via curveRcore’s `tidy_samples()` / `tidy_grid()` so the
  adapter never reaches into object internals. Also accepts a plain data
  frame for the legacy / foreign-data path.
- [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
  — curveR-ecosystem entry point. Wraps
  [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
  with ecosystem defaults (`scale_predictor = "se"`, automatic
  `curve_id` plate grouping, `weight_data` input). Returns a
  `precision_weights` S3 object with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods.
- [`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
  — applies fitted (phi, beta1) to new observations or a precision grid
  without refitting. Accepts a `calibration_result(_multiplate)`
  directly (routed through `as_weight_data(source = "grid")`).
- [`join_weights()`](https://immunoplex.github.io/curveRweights/reference/join_weights.md)
  — left-joins the estimated `w`, `w_norm`, and `sigma` columns from a
  `precision_weights` fit onto a user data frame by key.

### Scale predictor overhaul

- [`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md)
  gains a `predictor` argument (`"se"`, `"pcov"`, `"se_over_conc"`,
  `"auto"`) and `response_is_log10` flag. The recommended predictor is
  now `"se"` (uncapped `se_concentration` on the log10 scale), which
  preserves the full precision gradient. `"auto"` reproduces the legacy
  pcov-first behaviour.
- [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
  gains `predictor` and `response_is_log10` passthrough arguments.

### Existing functions (from 0.1.0, now refined)

- [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
  — core: joint Bayesian location-scale model with saturated cell means
  and shared power-law scale.
- [`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md)
  — per-group wrapper (e.g. antigen × source) with combined
  `scale_table` output.
- [`apply_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md)
  — stage-2 weight application from a saved `scale_table` without
  refitting.
- [`compute_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/compute_saturated_weights.md)
  — deterministic weight computation from (gamma_0, gamma_1).
- [`weight_diagnostics()`](https://immunoplex.github.io/curveRweights/reference/weight_diagnostics.md)
  — n_eff, weight_ratio, Gini coefficient.
- [`diagnose_cv()`](https://immunoplex.github.io/curveRweights/reference/diagnose_cv.md)
  — checks whether the precision index has sufficient variation for
  beta1 estimation.
- [`interpret_beta1()`](https://immunoplex.github.io/curveRweights/reference/interpret_beta1.md)
  — classifies beta1 into precision-weighting regimes.
- `example_assay` dataset: 48,224 observations from a Luminex multiplex
  immunoassay (11 antigens, 10 features, 150 subjects, 4 timepoints, 15
  plates).

## curveRweights 0.2.1

### Verified compatible with curveRcore 0.3.0 (mask-aware preprocessing)

- No code changes required. Weighting operates purely on the standard
  points passed to the fit, which are the *included* subset
  (worker-filtered). Verified no references to preprocessing, blanks, or
  the database anywhere in `R/`; weights cannot be influenced by masked
  points.

### Bug fix: `predict_weights()` on a `calibration_result(_multiplate)` always failed

- `as_weight_data(source = "grid")` required every `design` column (e.g.
  `timeperiod`, `cohort_arm`) to be present in the extracted table, but
  [`curveRcore::tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.html)
  returns a per-curve concentration *profile* (keyed by
  `curve_id`/concentration), not a per-design-cell table — it never
  carries the original `samples` design columns. Since
  [`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
  always calls
  `as_weight_data(newdata, design = object$design, source = "grid")`
  internally, calling it with a `calibration_result(_multiplate)` (its
  documented, intended usage) errored unconditionally with
  `"design column(s) not found"`.
- [`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
  now only requires/validates `design` and builds the `.cell`
  saturated-cell factor (with its within-cell-replication checks) for
  `source = "samples"` — the fitting table consumed by
  [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md).
  For `source = "grid"`, `design` is intersected with whatever columns
  are actually present (currently none, from `tidy_grid()`) and the
  cell/replication checks are skipped, since
  [`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
  only ever reads `se`/`concentration`/`pcov` from the profile.
