# Compute Summary Diagnostics for Precision Weights

Summarises the distribution and effective information content of a set
of precision weights.

## Usage

``` r
weight_diagnostics(w)
```

## Arguments

- w:

  Numeric vector of weights (e.g., `w_saturated` or `w_saturated_norm`).

## Value

Named list:

- n_obs:

  Total observations (including NA).

- n_valid:

  Observations with finite, positive weights.

- n_eff:

  Effective sample size: \\\[\sum w_i\]^2 / \sum w_i^2\\. Equals n_valid
  when all weights are equal; decreases as weights become more
  heterogeneous.

- eff_ratio:

  n_eff / n_valid. Ranges from 0 to 1; 1 = uniform.

- weight_ratio:

  max(w) / min(w) among valid weights.

- gini:

  Gini coefficient of weights. 0 = perfectly uniform, approaching 1 =
  highly concentrated.

## Examples

``` r
w <- c(1.5, 1.2, 0.8, 0.3, 0.1)
weight_diagnostics(w)
#> $n_obs
#> [1] 5
#> 
#> $n_valid
#> [1] 5
#> 
#> $n_eff
#> [1] 3.4
#> 
#> $eff_ratio
#> [1] 0.687
#> 
#> $weight_ratio
#> [1] 15
#> 
#> $gini
#> [1] 0.379
#> 
```
