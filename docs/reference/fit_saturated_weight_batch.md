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
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#>   sigma fixef rows: sigma_Intercept, cellvaccine_a.timepoint_1, cellvaccine_b.timepoint_1, cellvaccine_a.timepoint_2, cellvaccine_b.timepoint_2, cellvaccine_a.timepoint_3, cellvaccine_b.timepoint_3, cellvaccine_a.timepoint_4, cellvaccine_b.timepoint_4, sigma_log_cv
#>   sigma fixef cols: Estimate, Est.Error, Q2.5, Q97.5
#>   phi = 2.885  [2.03, 4.07]
#>   beta1 = 1.024  [0.829, 1.2]
#>   interpretation: moderate precision weighting
#>   n_eff = 415.1 of 506 (ratio = 0.82)
#>   weight_ratio = 1027  gini = 0.259
#> 
#> [2/3] antigen=pt
#> fit_saturated_weight: 511 of 512 observations usable (1 removed); 8 cell levels
#>   cv: OK: sd(log_cv) = 0.605; beta1 identifiable from 511 observations
#>   location: yi ~ 0 + cell + (1 | plate)
#>   scale:    sigma ~ log_cv
#>   fitting brms model (1000 iter, 2 chains)...
#>   sigma fixef rows: sigma_Intercept, cellvaccine_a.timepoint_1, cellvaccine_b.timepoint_1, cellvaccine_a.timepoint_2, cellvaccine_b.timepoint_2, cellvaccine_a.timepoint_3, cellvaccine_b.timepoint_3, cellvaccine_a.timepoint_4, cellvaccine_b.timepoint_4, sigma_log_cv
#>   sigma fixef cols: Estimate, Est.Error, Q2.5, Q97.5
#>   phi = 2.234  [1.5, 3.55]
#>   beta1 = 0.7743  [0.634, 0.935]
#>   interpretation: compressed precision weighting
#>   n_eff = 411.6 of 511 (ratio = 0.805)
#>   weight_ratio = 859.6  gini = 0.277
#> 
#> [3/3] antigen=fha
#> fit_saturated_weight: 511 of 512 observations usable (1 removed); 8 cell levels
#>   cv: OK: sd(log_cv) = 0.546; beta1 identifiable from 511 observations
#>   location: yi ~ 0 + cell + (1 | plate)
#>   scale:    sigma ~ log_cv
#>   fitting brms model (1000 iter, 2 chains)...
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>   sigma fixef rows: sigma_Intercept, cellvaccine_a.timepoint_1, cellvaccine_b.timepoint_1, cellvaccine_a.timepoint_2, cellvaccine_b.timepoint_2, cellvaccine_a.timepoint_3, cellvaccine_b.timepoint_3, cellvaccine_a.timepoint_4, cellvaccine_b.timepoint_4, sigma_log_cv
#>   sigma fixef cols: Estimate, Est.Error, Q2.5, Q97.5
#>   phi = 3.228  [2.26, 4.61]
#>   beta1 = 1.034  [0.904, 1.16]
#>   interpretation: calibrated (pcov ~ residual SD)
#>   n_eff = 288.3 of 511 (ratio = 0.564)
#>   weight_ratio = 24060  gini = 0.469

# Scale summary: one row per antigen
batch$scale_table
#>   antigen      phi     beta1   phi_lo   phi_hi  beta1_lo  beta1_hi
#> 1     prn 2.884783 1.0238494 2.030100 4.074802 0.8291515 1.2042617
#> 2      pt 2.234074 0.7743324 1.502670 3.549606 0.6339866 0.9351635
#> 3     fha 3.227635 1.0338441 2.264336 4.611407 0.9043599 1.1608369
#>                    interpretation n_fit n_eff weight_ratio
#> 1    moderate precision weighting   506 415.1       1027.0
#> 2  compressed precision weighting   511 411.6        859.6
#> 3 calibrated (pcov ~ residual SD)   511 288.3      24060.0

# Weighted data for one comparison
dat_comparison <- batch$data[batch$data$group_a %in% c("vaccine_a", "vaccine_b") &
                             batch$data$group_b == "timepoint_3", ]
nrow(dat_comparison)
#> [1] 450
# }
```
