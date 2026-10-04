# Build precision-weight input from a curveR calibration result

Converts a curveRcore `calibration_result` or
`calibration_result_multiplate` (the object returned by
[`curveRfreq::fit_calibration_freq_multiplate()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq_multiplate.html)
or
[`curveRbayes::fit_calibration_bayes()`](https://immunoplex.github.io/curveRbayes/reference/fit_calibration_bayes.html))
into the tidy, standardized data frame consumed by
[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md).
A plain data frame is also accepted (validation/passthrough), which is
the escape hatch for the original foreign-data workflow.

## Usage

``` r
as_weight_data(
  x,
  design,
  source = c("samples", "grid"),
  conc_scale = c("log10", "natural"),
  include_plate = TRUE,
  plate_in_cell = FALSE,
  drop_oor = FALSE,
  ...
)

# S3 method for class 'calibration_result_multiplate'
as_weight_data(
  x,
  design,
  source = c("samples", "grid"),
  conc_scale = c("log10", "natural"),
  include_plate = TRUE,
  plate_in_cell = FALSE,
  drop_oor = FALSE,
  ...
)

# S3 method for class 'calibration_result'
as_weight_data(
  x,
  design,
  source = c("samples", "grid"),
  conc_scale = c("log10", "natural"),
  include_plate = TRUE,
  plate_in_cell = FALSE,
  drop_oor = FALSE,
  ...
)

# S3 method for class 'data.frame'
as_weight_data(
  x,
  design,
  source = c("samples", "grid"),
  conc_scale = c("log10", "natural"),
  include_plate = TRUE,
  plate_in_cell = FALSE,
  drop_oor = FALSE,
  ...
)
```

## Arguments

- x:

  A `calibration_result`, `calibration_result_multiplate`, or
  `data.frame`.

- design:

  Character vector of design-group column names that define the
  saturated cells (e.g. `c("timeperiod", "cohort_arm")`). Required, and
  must be present on the input, when `source = "samples"`: for a
  `calibration_result(_multiplate)` these are carried through from the
  original `samples` data frame passed at fit time — if they are absent,
  supply them on the fitting input first, or use the `data.frame` method
  with a pre-joined frame. Ignored (silently intersected with the
  available columns) when `source = "grid"`: see `source` below.

- source:

  `"samples"` (default) extracts per-sample predictions, validates that
  every `design` column is present, and builds the `.cell`
  saturated-cell factor used by
  [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md).
  `"grid"` extracts the precision grid (used by
  [`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
  to build a continuous weight profile):
  [`curveRcore::tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.html)
  returns a per-curve concentration profile, not a per-design-cell
  table, so it never carries the original `samples` design columns.
  `design` and the `.cell`/within-cell-replication checks therefore do
  not apply to `source = "grid"` — the profile only needs
  `se`/`concentration`/`pcov`, which
  [`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
  consumes directly.

- conc_scale:

  Location-response scale: `"log10"` (default, uses
  `predicted_concentration`) or `"natural"` (uses `final_concentration`
  when present, else `10^predicted_concentration`).

- include_plate:

  Logical: attach `curve_id` as a column and make it available as a
  plate grouping. Default `TRUE`.

- plate_in_cell:

  Logical: fold `curve_id` into the cell factor alongside `design`.
  Default `FALSE` — plate variation is left to be absorbed by `phi` /
  the plate random effect (see the plan, §5).

- drop_oor:

  Logical: drop rows with `pcov_pass == FALSE`? Default `FALSE` (keep
  but warn) — not discarding out-of-range observations is the point of
  the package.

- ...:

  Passed to methods.

## Value

A data frame (class `weight_data`) with at least: `obs_id`, `sampleid`
(if present), `curve_id` (if `include_plate`), `concentration` (location
response on `conc_scale`), `predicted_concentration` (log10, for the
estimator), `se` (= `se_concentration`, the canonical uncapped scale
predictor), `pcov` (reference only), `pcov_pass`, and, when
`source = "samples"` (or any `design` columns are present on a `"grid"`
input), the `design` columns and `.cell`. Carries attributes
`conc_scale`, `is_log_independent`, `design` (resolved to the columns
actually present), `plate_in_cell`.

## Details

Extraction from the S3 object is delegated to
[`curveRcore::tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.html)
/
[`curveRcore::tidy_grid()`](https://immunoplex.github.io/curveRcore/reference/tidy_grid.html)
so that the data contract has a single owner; this adapter never reaches
into object internals.

## See also

[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md),
[`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md),
[`curveRcore::tidy_samples()`](https://immunoplex.github.io/curveRcore/reference/tidy_samples.html)
