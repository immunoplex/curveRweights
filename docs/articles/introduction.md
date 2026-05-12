# Introduction to curveRweights

![](../reference/figures/logo.png)

## The problem: unequal precision across the calibration curve

Immunoassay measurements are determined via a calibration curve. Each
observation has a **posterior coefficient of variation (pcov)** that
quantifies how precisely the curve determines that concentration.
Observations near the limits of quantification have large pcov values;
midrange observations have small pcov.

Traditional approaches handle this with a **binary gate**: include if
between LLOQ and ULOQ, exclude otherwise. This discards information and
treats all included observations as equally precise.

`curveRweights` replaces this with **continuous, data-estimated
precision weights** via a Bayesian location-scale model.

## Setup

``` r
library(curveRweights)
library(dplyr)
#> Warning: package 'dplyr' was built under R version 4.5.2
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
data(example_assay)
```

The `example_assay` dataset contains 48,224 observations from a Luminex
multiplex immunoassay: 11 antigens, 10 features, 150 subjects, 4
timepoints, across 15 plates.

``` r
str(example_assay)
#> 'data.frame':    48224 obs. of  11 variables:
#>  $ plate                  : chr  "plate_1" "plate_1" "plate_1" "plate_1" ...
#>  $ nominal_sample_dilution: chr  "dil_6" "dil_6" "dil_6" "dil_6" ...
#>  $ patientid              : chr  "S013" "S013" "S013" "S013" ...
#>  $ group_a                : chr  "vaccine_b" "vaccine_b" "vaccine_b" "vaccine_b" ...
#>  $ group_b                : chr  "timepoint_3" "timepoint_3" "timepoint_3" "timepoint_3" ...
#>  $ antigen                : chr  "act" "dt" "ipv1" "ipv2" ...
#>  $ feature                : chr  "IgG1" "IgG1" "IgG1" "IgG1" ...
#>  $ mfi                    : num  43 1034 101 85 1169 ...
#>  $ predicted_concentration: num  1.2 39.58 8.74 1.33 11.9 ...
#>  $ se_concentration       : num  0.1546 3.5301 0.5288 0.0633 1.1052 ...
#>  $ pcov                   : num  0.1767 0.0789 0.102 0.0405 0.1409 ...
```

## Step 1: Prepare the precision index

[`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md)
computes the outcome variable (`yi = log10(concentration)`) and the
precision index (`cv_i`), preferring the stored pcov from the
calibration curve over the simpler `se / concentration` ratio.

``` r
# Focus on IgG1 responses to pertussis antigens
dat_igg1 <- example_assay |>
  filter(feature == "IgG1", antigen %in% c("pt", "fha", "prn"))

dat_igg1 <- prepare_cv(dat_igg1, pcov_col = "pcov")
head(dat_igg1[, c("antigen", "patientid", "yi", "cv_i", "log_cv", "cv_source")])
#>   antigen patientid         yi       cv_i    log_cv cv_source
#> 1     prn      S013  0.9707595 0.11265055 -2.183465      pcov
#> 2      pt      S013 -0.1981600 0.04432566 -3.116192      pcov
#> 3     fha      S013 -0.1086106 0.07826967 -2.547595      pcov
#> 4     prn      S013  0.4622342 0.12126340 -2.109790      pcov
#> 5      pt      S013 -0.3222110 0.04972253 -3.001297      pcov
#> 6     fha      S013 -0.3375811 0.08815158 -2.428697      pcov
```

## Step 2: Diagnose pcov variation

Before fitting, check whether pcov has enough variation to identify the
precision exponent (beta1). If `sd(log_cv) < 0.05`, the scale model
falls back to uniform weights.

``` r
# Check one antigen
d_prn <- dat_igg1 |> filter(antigen == "prn")
cv_diag <- diagnose_cv(d_prn)
cat(cv_diag$message, "\n")
#> OK: sd(log_cv) = 0.487; beta1 identifiable from 506 observations
cat("cv range:", cv_diag$cv_min, "to", cv_diag$cv_max, "\n")
#> cv range: 0.1056 to 3.12
cat("sd(log_cv):", cv_diag$log_cv_sd, "\n")
#> sd(log_cv): 0.4872
```

## Step 3: Create the saturated cell-means factor

The key design principle is that the cell factor is created
**externally**. This gives you full control over what constitutes a
“cell” in your experiment.

For this dataset, cells are defined by the combination of vaccine group
(`group_a`) and timepoint (`group_b`):

``` r
dat_igg1$cell <- interaction(dat_igg1$group_a, dat_igg1$group_b, drop = TRUE)
table(dat_igg1$cell)
#> 
#> vaccine_a.timepoint_1 vaccine_b.timepoint_1 vaccine_a.timepoint_2 
#>                   195                   216                   123 
#> vaccine_b.timepoint_2 vaccine_a.timepoint_3 vaccine_b.timepoint_3 
#>                   126                   222                   228 
#> vaccine_a.timepoint_4 vaccine_b.timepoint_4 
#>                   201                   225
```

Each cell gets its own mean in the location model (`yi ~ 0 + cell`). The
scale model (`log(sigma) ~ log_cv`) is shared across all cells,
estimated from ~1,500 observations instead of ~150.

## Step 4: Fit the model

[`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md)
fits one model per group (here, per antigen). Each group gets its own
phi and beta1 because different antigens may have different
pcov-to-variance relationships.

``` r
batch <- fit_saturated_weight_batch(
  datg       = dat_igg1,
  group_vars = "antigen",
  cell_col   = "cell",
  pcov_col   = "pcov",
  plate_col  = "plate",
  iter       = 4000,
  warmup     = 1000,
  chains     = 4,
  cores      = 4
)
```

``` r
batch$scale_table
```

A typical `scale_table` for these three antigens might look like:

| antigen | phi  | beta1 | interpretation                 | n_fit | n_eff | weight_ratio |
|---------|------|-------|--------------------------------|-------|-------|--------------|
| fha     | 0.97 | 0.78  | compressed precision weighting | 506   | 440   | 4.2          |
| prn     | 1.84 | 1.95  | amplified precision weighting  | 506   | 390   | 18.7         |
| pt      | 1.21 | 1.03  | moderate precision weighting   | 506   | 460   | 5.8          |

Key observations:

- **fha** has phi near 1 and beta1 \< 1: the calibration curve is well
  calibrated and pcov slightly overstates the precision gradient.
- **prn** has the largest phi (excess variance) and beta1 near 2: pcov
  strongly drives weight differentiation, and observations with high
  pcov are heavily downweighted.
- **pt** is close to the theoretical prediction (phi ~ 1, beta1 ~ 1).

## Step 5: Extract weights for a specific comparison

After fitting, the weights are attached to every row of the data. Filter
to any subset for downstream analysis:

``` r
# Compare vaccine_b vs vaccine_a at timepoint_3
dat_comparison <- batch$data |>
  filter(group_a %in% c("vaccine_a", "vaccine_b"),
         group_b == "timepoint_3")

# Weights are ready
summary(dat_comparison$w_saturated_norm)
```

## Step 6: Downstream arm-effect testing with metafor

The weights can be used directly in
[`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html).
Because phi already captures residual heterogeneity (the tau-squared
analog), use `method = "FE"` to avoid double-counting:

``` r
library(metafor)

# For one antigen at timepoint_3
dat_prn_t3 <- batch$data |>
  filter(antigen == "prn", group_b == "timepoint_3",
         group_a %in% c("vaccine_a", "vaccine_b"))

# z-score the outcome within this comparison
dat_prn_t3$conc_z <- scale(dat_prn_t3$predicted_concentration)

# vi = sigma_i^2 = 1 / w_saturated
dat_prn_t3$vi <- 1 / dat_prn_t3$w_saturated

# Set reference level
dat_prn_t3$group_a <- factor(dat_prn_t3$group_a,
                              levels = c("vaccine_a", "vaccine_b"))

# Fixed-effects meta-regression (tau2 = 0 by design)
fit_rma <- rma(
  yi   = conc_z,
  vi   = vi,
  mods = ~ group_a,
  data = dat_prn_t3,
  method = "FE"
)

summary(fit_rma)
```

The key output is the `group_avaccine_b` coefficient: the estimated
difference in z-scored concentration between vaccine_b and vaccine_a,
weighted by calibration-curve precision.

## Step 7: Survey-weighted rank test (nonparametric alternative)

For a nonparametric arm comparison that respects the precision weights,
use
[`survey::svyranktest()`](https://rdrr.io/pkg/survey/man/svyranktest.html):

``` r
library(survey)

dat_prn_t3$conc_z <- scale(dat_prn_t3$predicted_concentration)

des <- svydesign(ids = ~1, weights = ~w_saturated, data = dat_prn_t3)
rt  <- svyranktest(conc_z ~ group_a, design = des)
rt
```

## Step 8: Apply saved weights to new data

If you receive new data from the same assay platform, apply the
previously estimated scale parameters without re-fitting:

``` r
# Save the scale_table from the batch fit
saveRDS(batch$scale_table, "scale_table_IgG1.rds")

# Later, with new data:
scale_table <- readRDS("scale_table_IgG1.rds")

dat_new_weighted <- apply_saturated_weights(
  df          = dat_new,
  scale_table = scale_table,
  group_vars  = "antigen",
  pcov_col    = "pcov"
)
```

## Interpreting phi and beta1

The two scale parameters have direct scientific meaning:

**phi = exp(gamma_0)** is the baseline scaling factor.

- phi = 1: pcov is a perfectly calibrated proxy for residual SD
- phi \> 1: there is excess variance beyond what pcov predicts
  (biological scatter, plate effects, matrix interference)
- phi \< 1: pcov overstates the true residual variation (unusual)

**beta1 = gamma_1** is the precision exponent.

- beta1 = 1: the theoretical delta-method prediction; sigma is
  proportional to pcov
- beta1 \> 1: amplified weighting; high-pcov observations are
  downweighted more aggressively than theory predicts
- beta1 \< 1: compressed weighting; the pcov gradient is shallower than
  theory predicts
- beta1 near 0: pcov carries almost no information about precision;
  weights are near-uniform

Use
[`interpret_beta1()`](https://immunoplex.github.io/curveRweights/reference/interpret_beta1.md)
for a human-readable classification:

``` r
interpret_beta1(0.98, 0.85, 1.12)
#> [1] "calibrated (pcov ~ residual SD)"
interpret_beta1(2.01, 1.47, 2.57)
#> [1] "amplified precision weighting"
interpret_beta1(0.56, 0.30, 0.82)
#> [1] "compressed precision weighting"
```

## Weight diagnostics

[`weight_diagnostics()`](https://immunoplex.github.io/curveRweights/reference/weight_diagnostics.md)
provides a quick summary of how much the weights differentiate
observations:

``` r
# Example with known weights
w <- c(rep(1.5, 50), rep(0.8, 30), rep(0.2, 20))
weight_diagnostics(w)
#> $n_obs
#> [1] 100
#> 
#> $n_valid
#> [1] 100
#> 
#> $n_eff
#> [1] 80.1
#> 
#> $eff_ratio
#> [1] 0.801
#> 
#> $weight_ratio
#> [1] 7.5
#> 
#> $gini
#> [1] 0.263
```

Key metrics:

- **n_eff**: effective sample size; equals n when all weights are
  uniform, decreases as weights become more heterogeneous
- **eff_ratio**: n_eff / n; 1.0 = uniform, lower = more differentiation
- **weight_ratio**: max / min weight; large values indicate aggressive
  downweighting of some observations
- **gini**: Gini coefficient; 0 = uniform, approaching 1 = concentrated

## Prior sensitivity

The default priors center beta1 at 1 with SD = 0.5. For antigens where
beta1 is far from 1 (like prn with beta1 near 2), the prior exerts
meaningful shrinkage. Test sensitivity by comparing results with wider
priors:

``` r
# Default prior
sw_default <- fit_saturated_weight(dat_prn, cell_col = "cell",
                                    pcov_col = "pcov", plate_col = "plate")

# Wider prior
sw_wide <- fit_saturated_weight(
  dat_prn, cell_col = "cell", pcov_col = "pcov", plate_col = "plate",
  prior_gamma1 = brms::set_prior("normal(1, 1.0)", class = "b", dpar = "sigma")
)

cat("Default beta1:", sw_default$beta1, "\n")
cat("Wide    beta1:", sw_wide$beta1, "\n")
```

If beta1 shifts materially, the data are not fully informative and the
prior is doing real work. The wider prior gives weights that better
reflect the data; the narrower prior gives more conservative
(closer-to-uniform) weights.
