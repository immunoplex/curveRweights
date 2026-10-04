# Changelog

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
