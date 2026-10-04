# Fit Saturated Weight Models Across Multiple Groups

Loops over groups defined by `group_vars` (e.g., antigen x source) and
fits one
[`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
model per group. Returns a combined data frame with weights and a
summary table of scale estimates.

## Usage

``` r
fit_saturated_weight_batch(
  datg,
  group_vars = c("antigen", "source"),
  cell_col = "cell",
  ...
)
```

## Arguments

- datg:

  Data frame containing all groups.

- group_vars:

  Character vector: column names defining groups. Each unique
  combination gets its own model.

- cell_col:

  Character: name of the saturated cell-means factor (created
  externally).

- ...:

  Additional arguments passed to
  [`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md),
  such as `pcov_col`, `plate_col`, `iter`, `warmup`, `chains`, `cores`.

## Value

Named list:

- data:

  Full data frame with `w_saturated` and `w_saturated_norm` columns
  added (all groups combined).

- scale_table:

  A [tibble](https://tibble.tidyverse.org/reference/tibble-package.html)
  with one row per group containing: group key columns, `phi`, `beta1`,
  credible intervals, `interpretation`, `n_fit`, `n_eff`,
  `weight_ratio`.

- fits:

  Named list of `brmsfit` objects indexed by group label.

- diagnostics:

  Named list of per-group diagnostic lists.

## Details

Each group gets its own (phi, beta1) because the pcov-to-variance
relationship may differ across antigens (different 4PL curve shapes) and
standard curve sources (different concentration ranges).

## See also

[`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
for the per-group fitting function,
[`apply_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md)
for applying a saved `scale_table` to new data without re-fitting.

## Examples

``` r
# \donttest{
data(example_assay)

# Select IgG1 for pertussis antigens
dat_igg1 <- example_assay[example_assay$feature == "IgG1" &
                          example_assay$antigen %in% c("pt", "fha", "prn"), ]
dat_igg1$cell <- interaction(dat_igg1$group_a, dat_igg1$group_b, drop = TRUE)

# Fit across antigens (reduced iterations for speed)
batch <- fit_saturated_weight_batch(
  datg       = dat_igg1,
  group_vars = c("antigen"),
  cell_col   = "cell",
  pcov_col   = "pcov",
  plate_col  = "plate",
  iter = 1000, warmup = 500, chains = 2, cores = 2
)
#> fit_saturated_weight_batch: 3 groups
#> 
#> [1/3] antigen=prn
#> fit_saturated_weight: 506 of 512 observations usable (6 removed); 8 cell levels
#>   cv: OK: sd(log_cv) = 0.487; beta1 identifiable from 506 observations
#>   location: yi ~ 0 + cell + (1 | plate)
#>   scale:    sigma ~ log_cv
#>   fitting brms model (1000 iter, 2 chains)...
#> Error in .fun(model_code = .x1) : 
#>   Boost not found; call install.packages('BH')
#> Warning: Group [antigen=prn] failed: fit_saturated_weight: brms::brm() failed: Boost not found; call install.packages('BH')
#> 
#> [2/3] antigen=pt
#> fit_saturated_weight: 511 of 512 observations usable (1 removed); 8 cell levels
#>   cv: OK: sd(log_cv) = 0.605; beta1 identifiable from 511 observations
#>   location: yi ~ 0 + cell + (1 | plate)
#>   scale:    sigma ~ log_cv
#>   fitting brms model (1000 iter, 2 chains)...
#> Error in .fun(model_code = .x1) : 
#>   Boost not found; call install.packages('BH')
#> Warning: Group [antigen=pt] failed: fit_saturated_weight: brms::brm() failed: Boost not found; call install.packages('BH')
#> 
#> [3/3] antigen=fha
#> fit_saturated_weight: 511 of 512 observations usable (1 removed); 8 cell levels
#>   cv: OK: sd(log_cv) = 0.546; beta1 identifiable from 511 observations
#>   location: yi ~ 0 + cell + (1 | plate)
#>   scale:    sigma ~ log_cv
#>   fitting brms model (1000 iter, 2 chains)...
#> Error in .fun(model_code = .x1) : 
#>   Boost not found; call install.packages('BH')
#> Warning: Group [antigen=fha] failed: fit_saturated_weight: brms::brm() failed: Boost not found; call install.packages('BH')

# Scale summary: one row per antigen
batch$scale_table
#>   antigen phi beta1 phi_lo phi_hi beta1_lo beta1_hi interpretation n_fit n_eff
#> 1     prn  NA    NA     NA     NA       NA       NA         failed     0    NA
#> 2      pt  NA    NA     NA     NA       NA       NA         failed     0    NA
#> 3     fha  NA    NA     NA     NA       NA       NA         failed     0    NA
#>   weight_ratio
#> 1           NA
#> 2           NA
#> 3           NA

# Weighted data for one comparison
dat_comparison <- batch$data[batch$data$group_a %in% c("vaccine_a", "vaccine_b") &
                             batch$data$group_b == "timepoint_3", ]
nrow(dat_comparison)
#> [1] 450
# }
```
