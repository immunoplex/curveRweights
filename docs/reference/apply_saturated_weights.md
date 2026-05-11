# Apply Previously Estimated Scale Parameters to New Data

When you have scale estimates from a previous fit (via
[`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
or
[`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md))
and want to apply them to a new or subsetted data frame without
re-fitting the brms model.

## Usage

``` r
apply_saturated_weights(
  df,
  scale_table,
  group_vars = c("antigen", "source"),
  concentration_col = "predicted_concentration",
  se_col = "se_concentration",
  pcov_col = "pcov"
)
```

## Arguments

- df:

  Data frame to receive weights.

- scale_table:

  A `scale_table` from
  [`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md),
  or a one-row tibble/data.frame with at least `phi` and `beta1` (or
  `gamma_0` and `gamma_1`).

- group_vars:

  Character vector: column names to match between `df` and `scale_table`
  for group-specific scale parameters. `NULL` = apply a single global
  (phi, beta1) to all rows.

- concentration_col:

  Character: predicted concentration column name.

- se_col:

  Character: SE of concentration column name.

- pcov_col:

  Character: posterior CV column name.

## Value

The input data frame with added columns: `yi`, `cv_i`, `log_cv`,
`cv_source`, `sigma_i`, `w_saturated`, `w_saturated_norm`.

## Details

This is the Stage 2 function in the recommended two-stage workflow: fit
once on the full dataset (all arms, all timepoints), then apply weights
to any comparison subset.

## See also

[`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md)
for producing the `scale_table`,
[`compute_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/compute_saturated_weights.md)
for the underlying weight computation.

## Examples

``` r
# \donttest{
data(example_assay)

# Fit on one group first
dat_prn <- example_assay[example_assay$antigen == "prn" &
                         example_assay$feature == "IgG1", ]
dat_prn$cell <- interaction(dat_prn$group_a, dat_prn$group_b, drop = TRUE)
sw <- fit_saturated_weight(dat_prn, cell_col = "cell", pcov_col = "pcov",
                           plate_col = "plate",
                           iter = 1000, warmup = 500, chains = 2, cores = 2)
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

# Build a scale_table manually (or use batch$scale_table)
st <- data.frame(antigen = "prn", phi = sw$phi, beta1 = sw$beta1)

# Apply to new data without re-fitting
dat_new <- example_assay[example_assay$antigen == "prn" &
                         example_assay$feature == "IgG1" &
                         example_assay$group_b == "timepoint_3", ]
dat_weighted <- apply_saturated_weights(
  df          = dat_new,
  scale_table = st,
  group_vars  = "antigen",
  pcov_col    = "pcov"
)
summary(dat_weighted$w_saturated_norm)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#>  0.0302  0.8877  1.0943  1.0000  1.2029  1.3718       2 
# }
```
