# Diagnose Precision Index Variation for Scale Estimation

The scale submodel `log(sigma) ~ log_cv` requires meaningful spread in
log_cv across observations. If all observations have nearly identical
pcov (e.g., all in the well-determined midrange), beta1 is not
identifiable and the model falls back to intercept-only sigma (uniform
weights).

## Usage

``` r
diagnose_cv(df)
```

## Arguments

- df:

  Data frame with `cv_i` and `log_cv` columns, typically from
  [`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md).

## Value

Named list:

- n_finite:

  Number of observations with finite cv_i \> 0.

- cv_min, cv_median, cv_max:

  Summary statistics of cv_i.

- log_cv_sd:

  Standard deviation of log_cv.

- log_cv_range:

  Range (max - min) of log_cv.

- use_cv_slope:

  Logical: `TRUE` if `sd(log_cv) >= 0.05`.

- message:

  Human-readable summary of the diagnosis.

## Details

The threshold is `sd(log_cv) >= 0.05`. Below this, there is not enough
contrast in measurement precision across the concentration range to
distinguish differential weighting from uniform weighting.

## Examples

``` r
data(example_assay)
dat_sub <- example_assay[example_assay$antigen == "prn" &
                         example_assay$feature == "IgG1", ]
d <- prepare_cv(dat_sub, pcov_col = "pcov")
cv_diag <- diagnose_cv(d)
cat(cv_diag$message, "\n")
#> OK: sd(log_cv) = 0.487; beta1 identifiable from 506 observations 
cat("sd(log_cv) =", cv_diag$log_cv_sd, "\n")
#> sd(log_cv) = 0.4872 
cat("use slope?", cv_diag$use_cv_slope, "\n")
#> use slope? TRUE 
```
