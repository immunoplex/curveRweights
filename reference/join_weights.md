# Join estimated weights back onto a data frame

Left-joins the `w`, `w_norm`, and `sigma` columns from a
`precision_weights` fit onto a user data frame (e.g. the original
`$samples`) by a key.

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

`samples_df` with `w`, `w_norm`, `sigma` added.

## See also

[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
