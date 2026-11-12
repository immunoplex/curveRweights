# curveRweights 0.3.0

## Additive measurement-error model (new, opt-in)

* `fit_additive_weight()` fits a model in which each observation's variance is the
  *sum* of true between-subject variance and measurement variance,
  sigma^2 = tau^2 + kappa^2 * se^(2 * beta1). The weight 1/sigma^2 can never exceed
  1/tau^2, so a very precisely measured sample can no longer dominate (the
  multiplicative model of 0.2.x attributes all variance to the measurement). beta1 is
  fixed at 1 by default (`fix_beta1 = TRUE`) because tau and kappa are separated only
  through the spread of `se`; a message reports when kappa is weakly identified.
* Samples outside the quantifiable range enter through a censored likelihood
  (left, right or interval) instead of being treated as exact points or weighted:
  `prepare_censoring()` builds the censoring flags and the measurement SE at the limit
  from curveRcore's `censor_type` / `censor_lower` / `censor_upper`. A
  `(1 | patient)` term absorbs repeated measures of one subject; optional
  `sample_col` handles dilutions or replicates of one sample; optional `dilution_col`
  puts the response on the same (final-concentration) scale as curveRcore's bounds.
* `compute_additive_weights()` returns the weight table a downstream analysis needs:
  `v_meas`, `tau2`, `sigma_total`, `w`, `reliability` (tau^2 / (tau^2 + v), with a
  posterior interval), `w_norm` (mean 1 over quantified rows), `censor_type`, `bound`.
* New helpers: `normalize_weights()` (re-normalize to mean 1 within an analysis subset),
  `meas_variance()` (the variance to pass to `metafor::rma.mv(V = )`),
  `reliability_correct()` (regression calibration) and `attenuation_correct()`.

## Existing behaviour is unchanged

* `fit_precision_weights()` gains `model = c("multiplicative", "additive")`,
  `patient_col` and `use_censoring`. The default stays `"multiplicative"`, so existing
  callers (including the i-spi-compute worker) get exactly the 0.2.x result and shape.
  The `precision_weights` object now records `$estimates$model`.
* The defaults of `prepare_cv()`, `fit_saturated_weight()` and `apply_saturated_weights()`
  (`predictor = "auto"`) are unchanged; `apply_saturated_weights()` gains a `predictor`
  argument. `predict_weights()` and `join_weights()` understand both models.
* `as_weight_data()` gains `carry` and now carries `patientid`, `dilution` and the
  censoring columns through when present.
* `test-backward-compat.R` pins these defaults and the result shape the worker reads.

# curveRweights 0.2.2

* `vignette("precision-weighting")` gains a "Weights and censored samples"
  section: why weighting out-of-range samples biases location estimates, how to
  apply weights only to quantified samples and represent the rest by their
  censoring bounds (`classify_censoring_multiplate()` in curveRcore), and a note
  that censoring-aware weight estimation is planned, not yet implemented.

# curveRweights 0.2.0

## curveR ecosystem integration

* `as_weight_data()` — adapter that converts a `calibration_result` or
  `calibration_result_multiplate` (from curveRfreq or curveRbayes) into
  the standardised `weight_data` frame consumed by the estimator.
  Dispatches via curveRcore's `tidy_samples()` / `tidy_grid()` so the
  adapter never reaches into object internals. Also accepts a plain
  data frame for the legacy / foreign-data path.
* `fit_precision_weights()` — curveR-ecosystem entry point. Wraps
  `fit_saturated_weight()` with ecosystem defaults (`scale_predictor =
  "se"`, automatic `curve_id` plate grouping, `weight_data` input).
  Returns a `precision_weights` S3 object with `print()` and `summary()`
  methods.
* `predict_weights()` — applies fitted (phi, beta1) to new observations
  or a precision grid without refitting. Accepts a
  `calibration_result(_multiplate)` directly (routed through
  `as_weight_data(source = "grid")`).
* `join_weights()` — left-joins the estimated `w`, `w_norm`, and `sigma`
  columns from a `precision_weights` fit onto a user data frame by key.

## Scale predictor overhaul

* `prepare_cv()` gains a `predictor` argument (`"se"`, `"pcov"`,
  `"se_over_conc"`, `"auto"`) and `response_is_log10` flag. The
  recommended predictor is now `"se"` (uncapped `se_concentration`
  on the log10 scale), which preserves the full precision gradient.
  `"auto"` reproduces the legacy pcov-first behaviour.
* `fit_saturated_weight()` gains `predictor` and `response_is_log10`
  passthrough arguments.

## Existing functions (from 0.1.0, now refined)

* `fit_saturated_weight()` — core: joint Bayesian location-scale model
  with saturated cell means and shared power-law scale.
* `fit_saturated_weight_batch()` — per-group wrapper (e.g.
  antigen × source) with combined `scale_table` output.
* `apply_saturated_weights()` — stage-2 weight application from a saved
  `scale_table` without refitting.
* `compute_saturated_weights()` — deterministic weight computation from
  (gamma_0, gamma_1).
* `weight_diagnostics()` — n_eff, weight_ratio, Gini coefficient.
* `diagnose_cv()` — checks whether the precision index has sufficient
  variation for beta1 estimation.
* `interpret_beta1()` — classifies beta1 into precision-weighting
  regimes.
* `example_assay` dataset: 48,224 observations from a Luminex multiplex
  immunoassay (11 antigens, 10 features, 150 subjects, 4 timepoints,
  15 plates).
  
# curveRweights 0.2.1

## Verified compatible with curveRcore 0.3.0 (mask-aware preprocessing)

* No code changes required. Weighting operates purely on the standard points
  passed to the fit, which are the *included* subset (worker-filtered). Verified
  no references to preprocessing, blanks, or the database anywhere in `R/`;
  weights cannot be influenced by masked points.

## Bug fix: `predict_weights()` on a `calibration_result(_multiplate)` always failed

* `as_weight_data(source = "grid")` required every `design` column (e.g.
  `timeperiod`, `cohort_arm`) to be present in the extracted table, but
  `curveRcore::tidy_grid()` returns a per-curve concentration *profile*
  (keyed by `curve_id`/concentration), not a per-design-cell table — it never
  carries the original `samples` design columns. Since `predict_weights()`
  always calls `as_weight_data(newdata, design = object$design, source =
  "grid")` internally, calling it with a `calibration_result(_multiplate)`
  (its documented, intended usage) errored unconditionally with
  `"design column(s) not found"`.
* `as_weight_data()` now only requires/validates `design` and builds the
  `.cell` saturated-cell factor (with its within-cell-replication checks)
  for `source = "samples"` — the fitting table consumed by
  `fit_precision_weights()`. For `source = "grid"`, `design` is intersected
  with whatever columns are actually present (currently none, from
  `tidy_grid()`) and the cell/replication checks are skipped, since
  `predict_weights()` only ever reads `se`/`concentration`/`pcov` from the
  profile.
