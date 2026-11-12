# Normalize weights to mean 1 within an analysis subset

Rescales a weight column so its mean is 1 over the rows you will
actually analyse (optionally within groups).

## Usage

``` r
normalize_weights(df, w_col = "w", by = NULL, out_col = "w_norm")
```

## Arguments

- df:

  Data frame (already subset to the analysis rows).

- w_col:

  Weight column to rescale. Default `"w"`.

- by:

  Optional character vector of grouping columns; each group is
  normalized separately.

- out_col:

  Name of the new column. Default `"w_norm"`.

## Value

`df` with `out_col` (mean 1 over finite positive weights).

## Why this matters

Software that treats weights as case counts (`brms`'s `| weights(w)`,
`lm(weights = )`) assumes the weights sum to the number of observations.
Normalizing over the full fitting frame and then filtering to a subset
(one arm, one time point, quantified rows only) breaks that: the
subset's mean is no longer 1, and uncertainty is over- or under-stated.
Subset first, then call this function.

Rescaling never changes estimates from weighted least squares, but it
does change standard errors in software that reads weights as counts. Do
not use it to make weights fall between 0 and 1; for a bounded quantity
use `reliability`.

## Examples

``` r
d <- data.frame(g = c("a", "a", "b", "b"), w = c(1, 3, 10, 30))
normalize_weights(d, by = "g")
#>   g  w w_norm
#> 1 a  1    0.5
#> 2 a  3    1.5
#> 3 b 10    0.5
#> 4 b 30    1.5
```
