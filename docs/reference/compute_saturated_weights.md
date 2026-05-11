# Compute Precision Weights from Estimated Scale Parameters

Given the estimated log-scale intercept (gamma_0) and slope (gamma_1),
computes the per-observation residual SD (sigma_i), raw precision weight
(w_saturated), and mean-normalised weight (w_saturated_norm).

## Usage

``` r
compute_saturated_weights(df, gamma_0, gamma_1 = 0)
```

## Arguments

- df:

  Data frame with a `log_cv` column (from
  [`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md)).

- gamma_0:

  Numeric: estimated intercept of the log-sigma model (log(phi)).

- gamma_1:

  Numeric: estimated slope of the log-sigma model (beta1). Set to 0 for
  intercept-only (uniform weights). Default `0`.

## Value

The input data frame with additional columns:

- sigma_i:

  Estimated residual SD for each observation.

- w_saturated:

  Raw precision weight: `1 / sigma_i^2`.

- w_saturated_norm:

  Mean-normalised weight: `w_saturated / mean(w_saturated)`. The shape
  is identical to w_saturated; normalisation to mean = 1 is applied only
  for cross-group comparability.

## Details

The transformation is: \$\$\sigma_i = \exp(\gamma_0 + \gamma_1 \cdot
\log(\mathrm{cv}\_i))\$\$ \$\$w_i = 1 / \sigma_i^2\$\$

## Examples

``` r
data(example_assay)
dat_sub <- example_assay[example_assay$antigen == "prn" &
                         example_assay$feature == "IgG1", ]
d <- prepare_cv(dat_sub, pcov_col = "pcov")

# Apply hypothetical scale estimates
d <- compute_saturated_weights(d, gamma_0 = 0.5, gamma_1 = 1.2)
summary(d$w_saturated_norm)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
#> 0.0005172 0.5825621 1.1687081 1.0000000 1.4015256 1.7504343         6 
weight_diagnostics(d$w_saturated)
#> $n_obs
#> [1] 512
#> 
#> $n_valid
#> [1] 506
#> 
#> $n_eff
#> [1] 402.3
#> 
#> $eff_ratio
#> [1] 0.795
#> 
#> $weight_ratio
#> [1] 3384
#> 
#> $gini
#> [1] 0.284
#> 
```
