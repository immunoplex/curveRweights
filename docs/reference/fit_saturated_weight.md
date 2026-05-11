# Fit a Bayesian Location-Scale Model with Saturated Location and Shared Scale

Fits a brms Gaussian location-scale model where the location uses
saturated cell means (one coefficient per level of `cell_col`) to absorb
ALL systematic variation in the response, and the scale estimates the
power-law relationship between the calibration curve's posterior CV
(pcov) and residual variance, shared across all cells.

## Usage

``` r
fit_saturated_weight(
  df,
  cell_col = "cell",
  concentration_col = "predicted_concentration",
  se_col = "se_concentration",
  pcov_col = "pcov",
  plate_col = "plate",
  prior_gamma0 = brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma"),
  prior_gamma1 = brms::set_prior("normal(1, 0.5)", class = "b", dpar = "sigma"),
  prior_location = brms::set_prior("normal(0, 2)", class = "b"),
  prior_plate_sd = brms::set_prior("normal(0, 0.5)", class = "sd"),
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  adapt_delta = 0.95,
  seed = 42
)
```

## Arguments

- df:

  Data frame with observation-level data.

- cell_col:

  Character: name of the saturated cell-means factor column. Create this
  externally, e.g.:
  `df$cell <- interaction(df$Arm, df$Timeperiod, drop = TRUE)`.

- concentration_col:

  Character: predicted concentration column name.

- se_col:

  Character: SE of concentration column name.

- pcov_col:

  Character: posterior CV column name. `NULL` = compute from se/conc.

- plate_col:

  Character: plate column name for random intercept. `NULL` = no plate
  random effect.

- prior_gamma0:

  A
  [`brms::prior`](https://paulbuerkner.com/brms/reference/set_prior.html)
  object for gamma_0 (log phi).

- prior_gamma1:

  A
  [`brms::prior`](https://paulbuerkner.com/brms/reference/set_prior.html)
  object for gamma_1 (beta1).

- prior_location:

  A
  [`brms::prior`](https://paulbuerkner.com/brms/reference/set_prior.html)
  object for location (cell mean) coefficients.

- prior_plate_sd:

  A
  [`brms::prior`](https://paulbuerkner.com/brms/reference/set_prior.html)
  object for plate RE standard deviation.

- iter:

  Integer: MCMC iterations per chain.

- warmup:

  Integer: warmup iterations per chain.

- chains:

  Integer: number of MCMC chains.

- cores:

  Integer: parallel cores for chain computation.

- adapt_delta:

  Numeric: target HMC acceptance probability (0-1).

- seed:

  Integer: random seed for reproducibility.

## Value

Named list with elements: `fit` (brmsfit object), `phi` (baseline
scaling factor), `beta1` (precision exponent), `gamma_0` and `gamma_1`
(log-scale estimates), `gamma_0_CI` and `gamma_1_CI` (95 percent
credible intervals), `phi_CI` (back-transformed CI for phi),
`interpretation` (precision regime label), `effective_se_power` (2 \*
beta1), `data` (input data with weights added: `yi`, `cv_i`, `log_cv`,
`sigma_i`, `w_saturated`, `w_saturated_norm`), `diagnostics`
(convergence and weight diagnostics), `cv_diagnostics` (output of
`diagnose_cv`), `formula` (location and scale formula strings),
`priors_used` (priors passed to brms).

## Details

The location model uses `yi ~ 0 + cell [+ (1|plate)]` (saturated cell
means). The scale model uses `log(sigma) = gamma_0 + gamma_1 * log(cv)`.

Default prior on gamma_0 is Normal(0, 1), centering phi = exp(gamma_0)
at 1. Default prior on gamma_1 is Normal(1, 0.5), centering beta1 at the
delta-method prediction of 1, with a 95 percent interval from 0 to 2.

## See also

[`fit_saturated_weight_batch`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md),
[`apply_saturated_weights`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md),
[`prepare_cv`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md),
[`diagnose_cv`](https://immunoplex.github.io/curveRweights/reference/diagnose_cv.md),
[`weight_diagnostics`](https://immunoplex.github.io/curveRweights/reference/weight_diagnostics.md)

## Examples

``` r
# \donttest{
data(example_assay)

# Select one antigen/feature group and create the cell factor
dat_prn <- example_assay[example_assay$antigen == "prn" &
                         example_assay$feature == "IgG1", ]
dat_prn$cell <- interaction(dat_prn$group_a, dat_prn$group_b, drop = TRUE)

# Fit (reduced iterations for speed; use 4000/1000 for real analysis)
sw <- fit_saturated_weight(
  df        = dat_prn,
  cell_col  = "cell",
  pcov_col  = "pcov",
  plate_col = "plate",
  iter = 1000, warmup = 500, chains = 2, cores = 2
)
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

cat("phi =", sw$phi, " beta1 =", sw$beta1, "\n")
#> phi = 2.884783  beta1 = 1.023849 
cat(sw$interpretation, "\n")
#> moderate precision weighting 

# Weights are on sw$data
summary(sw$data$w_saturated_norm)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#> 0.001607 0.645223 1.168663 1.000000 1.364593 1.649590        6 
weight_diagnostics(sw$data$w_saturated)
#> $n_obs
#> [1] 512
#> 
#> $n_valid
#> [1] 506
#> 
#> $n_eff
#> [1] 415.1
#> 
#> $eff_ratio
#> [1] 0.82
#> 
#> $weight_ratio
#> [1] 1027
#> 
#> $gini
#> [1] 0.259
#> 
# }
```
