# Measurement variance of each observation

Returns the measurement-error variance \\v_i\\ (not the total variance).
This is the right input for `metafor::rma.mv(V = )` and `rma(vi = )`:
those functions add the between-subject variance \\\tau^2\\ themselves
and estimate it by REML, so passing the total variance would count the
signal twice.

## Usage

``` r
meas_variance(x, quantified_only = TRUE)
```

## Arguments

- x:

  A `precision_weights` object (additive model) or a data frame with a
  `v_meas` column.

- quantified_only:

  Logical: return `NA` for censored rows? Default `TRUE`.

## Value

Numeric vector, one value per row of the weight table.

## Examples

``` r
d <- data.frame(v_meas = c(0.01, 0.2), censor_type = c("none", "left"))
meas_variance(d)
#> [1] 0.01   NA
```
