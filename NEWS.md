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
