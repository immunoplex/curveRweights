# Prepare censored observations for the additive weight model

Turns the censoring labels produced by curveRcore
(`classify_censoring_multiplate()`: `censor_type`, `censor_lower`,
`censor_upper`) into the columns the additive model's censored
likelihood needs.

## Usage

``` r
prepare_censoring(
  df,
  censor_col = "censor_type",
  lower_col = "censor_lower",
  upper_col = "censor_upper",
  yi_col = "yi",
  se_col = "cv_i",
  se_bound_col = NULL,
  curve_col = "curve_id"
)
```

## Arguments

- df:

  Data frame containing `yi` (log10 response) and `se_col`.

- censor_col:

  Column holding `"none"`, `"left"`, `"interval"`, `"right"` (or `NA`).
  `NULL` means no censoring: every row is treated as quantified.

- lower_col, upper_col:

  Natural-scale censoring bounds.

- yi_col:

  Name of the log10 response column.

- se_col:

  Name of the measurement SE column (log10 scale).

- se_bound_col:

  Optional column with the measurement SE at the censoring limit. `NULL`
  (default) uses the nearest quantified sample of the same curve.

- curve_col:

  Column identifying the curve (used to find the nearest quantified
  sample). Ignored if absent.

## Value

`df` with added columns:

- cens_flag:

  brms censoring code: 0 quantified, -1 left, 1 right, 2 interval. `NA`
  if the row is censored but its bound is missing (such rows cannot be
  used).

- yi:

  Overwritten for censored rows by the (log10) censoring limit: the
  upper bound for left, the lower bound for right and interval.

- yi2:

  Upper bound for interval rows, otherwise equal to `yi`.

- se_meas:

  Measurement SE used by the model: the row's own SE for quantified
  rows, the SE at the limit for censored rows.

- bound:

  log10 censoring limit (`NA` for quantified rows).

## In plain terms

A sample outside the quantifiable range of the standard curve is not a
noisy number; it is a sample we only know to be *below* (or *above*)
some limit. Treating its extrapolated back-calculated value as an
ordinary measurement is wrong, and so is giving it a weight. Instead the
model uses the statement "the true value is below the lower limit"
(left-censored), "above the upper limit" (right-censored) or "between
the two" (interval). This function builds, for every row, the censoring
flag and the limit that statement refers to.

For a censored row the model also needs to know how precisely the limit
itself is determined. That is the measurement SD at the limit, not the
(extrapolated, usually huge) SE of the row's own back-calculated value.
By default it is taken from the quantified sample of the same curve
whose value is closest to the limit; supply `se_bound_col` to use a
better source (for example the SE read off the curve's precision grid at
the limit).

## Scales

`censor_lower`/`censor_upper` are on the natural concentration scale and
must be on the same scale as the response `yi` (see the `dilution_col`
argument of
[`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md)).
They are converted to log10 here.

## See also

[`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md)

## Examples

``` r
d <- data.frame(
  curve_id = "c1",
  yi = c(0.1, 0.8, 1.2, 1.6, NA),
  cv_i = c(NA, 0.12, 0.08, 0.15, NA),
  censor_type = c("left", "none", "none", "none", "right"),
  censor_lower = c(NA, NA, NA, NA, 70),
  censor_upper = c(5, NA, NA, NA, NA))
prepare_censoring(d, se_col = "cv_i")
#>   curve_id       yi cv_i censor_type censor_lower censor_upper      yi2
#> 1       c1 0.698970   NA        left           NA            5 0.698970
#> 2       c1 0.800000 0.12        none           NA           NA 0.800000
#> 3       c1 1.200000 0.08        none           NA           NA 1.200000
#> 4       c1 1.600000 0.15        none           NA           NA 1.600000
#> 5       c1 1.845098   NA       right           70           NA 1.845098
#>   cens_flag se_meas    bound
#> 1        -1    0.12 0.698970
#> 2         0    0.12       NA
#> 3         0    0.08       NA
#> 4         0    0.15       NA
#> 5         1    0.15 1.845098
```
