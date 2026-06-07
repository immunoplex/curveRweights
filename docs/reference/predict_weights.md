# Apply fitted precision weights to new observations

Uses the estimated `phi` and `beta1` from a `precision_weights` fit to
compute \\\sigma = \phi\\ s^{\beta_1}\\ and \\w = 1/\sigma^2\\ for new
observations — typically a precision grid (to draw a continuous weight
profile) or held-out samples — without refitting.

## Usage

``` r
predict_weights(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `precision_weights` object from
  [`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md).

- newdata:

  One of: a `weight_data` frame, a plain data frame with the predictor
  column, or a `calibration_result(_multiplate)` (which is run through
  `as_weight_data(source = "grid")` using the fit's design). `NULL`
  (default) re-applies to the fit's own data via the stored weights.

- ...:

  Passed to
  [`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
  when `newdata` is an ecosystem object.

## Value

A data frame with the predictor, `sigma`, `w`, and `w_norm` columns
added.

## See also

[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
