# Join estimated weights back onto a data frame

Left-joins the weight columns (`w`, `w_norm`, `sigma` and, for the
additive model, `v_meas`, `tau2`, `reliability` and its interval,
`censor_type`, `bound`) from a `precision_weights` fit onto a user data
frame (e.g. the original `$samples`) by a key.

## Usage

``` r
join_weights(samples_df, object, by = "sampleid")
```

## Arguments

- samples_df:

  Data frame to receive the weights.

- object:

  A `precision_weights` object.

- by:

  Join key column name(s). Default `"sampleid"`.

## Value

`samples_df` with the weight columns added. If you subset the result to
an analysis group, re-normalize with
[`normalize_weights()`](https://immunoplex.github.io/curveRweights/reference/normalize_weights.md)
before using `w_norm` as a case weight.

## See also

[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
