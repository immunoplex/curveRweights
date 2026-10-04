# Precision Weighting with curveRweights

``` r
library(curveRweights)
library(dplyr)
```

------------------------------------------------------------------------

## The problem: unequal precision across the calibration curve

Immunoassay measurements are determined via a calibration curve — a
nonlinear sigmoidal function relating concentration to assay response.
Each back-calculated concentration has a standard error
(`se_concentration`) and a percent coefficient of variation (`pcov`)
that quantify how precisely the curve determines that concentration.
Observations near the limits of quantification are imprecise (the curve
is flat there); midrange observations are precise (the curve is steep).

Traditional approaches handle this with a **binary quantification
gate**: observations between LLOQ and ULOQ are included with equal
weight, and everything else is excluded. This has two problems: it
discards information (out-of-range observations still carry signal) and
treats all included observations as equally precise (they are not).

curveRweights replaces the binary gate with **continuous, data-estimated
precision weights** that retain every observation and down-weight
imprecise ones smoothly.

------------------------------------------------------------------------

## The model

The scale submodel is a power-law relationship between the
calibration-curve precision index and the residual standard deviation:

``` math

\sigma_i = \phi \cdot se_i^{\beta_1}, \qquad w_i = \frac{1}{\sigma_i^2}
```

On the log scale this is linear:

``` math

\log(\sigma_i) = \gamma_0 + \gamma_1 \cdot \log(se_i)
```

where $`\phi = \exp(\gamma_0)`$ and $`\beta_1 = \gamma_1`$.

The scale model is estimated jointly with a **saturated cell-means
location model** — one coefficient per experimental condition — so that
the scale parameters are identified purely from within-cell residual
structure and are not contaminated by misspecified location effects.

### Interpreting phi and beta1

**phi = exp(gamma_0)** is the baseline scaling factor:

- phi = 1: the precision index is a perfectly calibrated residual SD.
- phi \> 1: there is excess variance beyond what the curve predicts —
  biological scatter, plate effects, matrix interference. Analogous to
  positive tau-squared in classical meta-analysis.
- phi \< 1: the precision index overstates the true residual variation
  (unusual).

**beta1 = gamma_1** is the precision exponent:

- beta1 = 1: the delta-method prediction — sigma is proportional to the
  precision index.
- beta1 \> 1: amplified weighting — high-se observations are
  downweighted more aggressively than theory predicts.
- beta1 \< 1: compressed weighting — the precision gradient is shallower
  than theory predicts.
- beta1 near 0: the precision index carries almost no information;
  weights are near-uniform.

[`interpret_beta1()`](https://immunoplex.github.io/curveRweights/reference/interpret_beta1.md)
gives a human-readable classification:

``` r
interpret_beta1(0.98, 0.85, 1.12)
#> [1] "calibrated (pcov ~ residual SD)"
interpret_beta1(2.01, 1.47, 2.57)
#> [1] "amplified precision weighting"
interpret_beta1(0.56, 0.30, 0.82)
#> [1] "compressed precision weighting"
```

### Relationship to meta-analytic tau-squared

In classical random-effects meta-regression, total per-observation
variance is $`\tau^2 + v_i`$ where $`v_i`$ is measurement error and
$`\tau^2`$ is residual heterogeneity. In curveRweights, $`\phi^2`$ plays
the role of $`\tau^2`$ (global inflation) and $`se_i^{2\beta_1}`$ plays
the role of $`v_i`$ (per-observation error). The key difference is that
phi scales multiplicatively rather than adding a fixed floor, and beta1
is data-estimated rather than fixed.

When using these weights in
[`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html),
set `method = "FE"` ($`\tau^2 = 0`$) because heterogeneity is already
captured by phi.

------------------------------------------------------------------------

## Ecosystem orientation

curveRweights has two API surfaces that produce the same weights.

**The curveR ecosystem API** (recommended) accepts a
`calibration_result_multiplate` from curveRfreq or curveRbayes directly:

[`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
→
[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
→
[`predict_weights()`](https://immunoplex.github.io/curveRweights/reference/predict_weights.md)
/
[`join_weights()`](https://immunoplex.github.io/curveRweights/reference/join_weights.md)

**The standalone API** works with any data frame that has a
concentration column and a precision index (se or pcov):

[`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md)
→
[`diagnose_cv()`](https://immunoplex.github.io/curveRweights/reference/diagnose_cv.md)
→
[`fit_saturated_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight.md)
/
[`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md)
→
[`apply_saturated_weights()`](https://immunoplex.github.io/curveRweights/reference/apply_saturated_weights.md)

Both paths call the same underlying brms model; the ecosystem API
handles extraction, column naming, and cell-factor construction
automatically.

------------------------------------------------------------------------

## Worked example: curveR ecosystem path

### Design columns must ride along

The single most important setup step: the design factors that define
your comparison cells (e.g. `timeperiod`, `cohort_arm`) must be
**columns on the `samples` data frame you pass to the fitter**.
curveRfreq and curveRbayes carry the original sample columns through
into `$samples` — they are not invented later.

``` r
library(curveRweights)

# samples_with_design has the usual assay columns
# PLUS timeperiod and cohort_arm
mp <- curveRfreq::fit_calibration_freq_multiplate(
  standards    = std_preprocessed,
  samples      = samples_with_design,
  response_var = "mfi",
  is_log_response    = TRUE,
  is_log_independent = TRUE,
  std_curve_conc     = 30
)
```

### Build the weight-input frame

[`as_weight_data()`](https://immunoplex.github.io/curveRweights/reference/as_weight_data.md)
extracts `se_concentration` (the scale predictor), `pcov` (reference),
the location response, the design columns, and the saturated `.cell`
factor. It validates that design columns are present and that cells have
within-cell replication.

``` r
wd <- as_weight_data(mp, design = c("timeperiod", "cohort_arm"))
head(wd[, c("sampleid", "curve_id", "concentration", "se", "pcov", ".cell")])
```

### Estimate weights

``` r
pw <- fit_precision_weights(wd)
pw
#> <precision_weights>
#>   phi   = 1.12  [0.98, 1.28]
#>   beta1 = 0.94  [0.71, 1.18]
#>   regime: moderate precision weighting
```

### Use the weights

Join them back onto the sample table, or draw a continuous weight
profile across the precision grid:

``` r
# Join weights onto samples by sampleid
samples_w <- join_weights(curveRcore::tidy_samples(mp), pw)

# Continuous weight profile across the precision grid
profile <- predict_weights(pw, newdata = mp)
plot(profile$cv_i, profile$w, log = "xy",
     xlab = "se_concentration", ylab = "weight")
```

Out-of-range samples are retained with small-but-nonzero weights rather
than being discarded:

``` r
table(in_range = samples_w$pcov_pass, has_weight = is.finite(samples_w$w))
```

------------------------------------------------------------------------

## Worked example: standalone path

This path works with any data frame — including data produced outside
the curveR ecosystem — as long as it has a concentration column and a
precision index.

### The example dataset

`example_assay` ships with curveRweights and contains 48,224
observations from a Luminex multiplex immunoassay: 11 antigens, 10
features, 150 subjects, 4 timepoints, across 15 plates.

``` r
data(example_assay)
str(example_assay, give.attr = FALSE)
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

### Step 1: Prepare the precision index

[`prepare_cv()`](https://immunoplex.github.io/curveRweights/reference/prepare_cv.md)
computes the location response (`yi = log10(concentration)`) and the
precision index (`cv_i`). The `predictor` argument controls which column
is used:

``` r
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

### Step 2: Diagnose precision-index variation

Before fitting, check whether the precision index has enough spread to
identify beta1. If `sd(log_cv) < 0.05`, the model falls back to
intercept-only sigma (uniform weights).

``` r
d_prn <- dat_igg1 |> filter(antigen == "prn")
cv_diag <- diagnose_cv(d_prn)
cat(cv_diag$message, "\n")
#> OK: sd(log_cv) = 0.487; beta1 identifiable from 506 observations
cat("cv range:", cv_diag$cv_min, "to", cv_diag$cv_max, "\n")
#> cv range: 0.1056 to 3.12
cat("sd(log_cv):", cv_diag$log_cv_sd, "\n")
#> sd(log_cv): 0.4872
```

### Step 3: Create the saturated cell-means factor

The cell factor is created **externally** — this gives you full control
over what constitutes a “cell” in your experiment. For this dataset,
cells are the combination of vaccine group and timepoint:

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

Each cell gets its own coefficient in the location model
(`yi ~ 0 + cell`). The scale model (`log(sigma) ~ log_cv`) is shared
across all cells, so it is estimated from the full dataset rather than
per-cell subsets.

#### Why saturated cell means?

Any misspecification in the location model (e.g. assuming additivity
when there is an interaction) leaves systematic patterns in the
residuals that contaminate the scale estimates. The saturated model
absorbs ALL systematic location variation by construction, so residuals
are the cleanest possible input to the scale model.

### Step 4: Fit the model

[`fit_saturated_weight_batch()`](https://immunoplex.github.io/curveRweights/reference/fit_saturated_weight_batch.md)
fits one model per group. Each group gets its own (phi, beta1) because
different antigens may have different precision-to-variance
relationships.

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

A typical `scale_table` for pertussis antigens:

| antigen | phi  | beta1 | interpretation                 | n_fit | n_eff | weight_ratio |
|---------|------|-------|--------------------------------|-------|-------|--------------|
| fha     | 0.97 | 0.78  | compressed precision weighting | 506   | 440   | 4.2          |
| prn     | 1.84 | 1.95  | amplified precision weighting  | 506   | 390   | 18.7         |
| pt      | 1.21 | 1.03  | moderate precision weighting   | 506   | 460   | 5.8          |

Observations:

- **fha** has phi near 1 and beta1 \< 1: the calibration curve is well
  calibrated and pcov slightly overstates the precision gradient.
- **prn** has the largest phi (excess variance) and beta1 near 2: pcov
  strongly drives weight differentiation.
- **pt** is close to the delta-method prediction (phi ~ 1, beta1 ~ 1).

### Step 5: Extract weights for a comparison

After fitting, weights are attached to every row. Filter to any subset
for downstream analysis:

``` r
dat_comparison <- batch$data |>
  filter(group_a %in% c("vaccine_a", "vaccine_b"),
         group_b == "timepoint_3")

summary(dat_comparison$w_saturated_norm)
```

------------------------------------------------------------------------

## Downstream analysis

### Weighted meta-regression with metafor

The weights can be used directly in
[`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html).
Because phi already captures residual heterogeneity, use `method = "FE"`
to avoid double-counting:

``` r
library(metafor)

dat_prn_t3 <- batch$data |>
  filter(antigen == "prn", group_b == "timepoint_3",
         group_a %in% c("vaccine_a", "vaccine_b"))

dat_prn_t3$conc_z  <- scale(dat_prn_t3$predicted_concentration)
dat_prn_t3$vi      <- 1 / dat_prn_t3$w_saturated
dat_prn_t3$group_a <- factor(dat_prn_t3$group_a,
                              levels = c("vaccine_a", "vaccine_b"))

fit_rma <- rma(yi = conc_z, vi = vi, mods = ~ group_a,
               data = dat_prn_t3, method = "FE")
summary(fit_rma)
```

The `group_avaccine_b` coefficient is the precision-weighted difference
in z-scored concentration between vaccine_b and vaccine_a.

### Survey-weighted rank test

For a nonparametric comparison that respects the precision weights:

``` r
library(survey)

des <- svydesign(ids = ~1, weights = ~w_saturated, data = dat_prn_t3)
rt  <- svyranktest(conc_z ~ group_a, design = des)
rt
```

------------------------------------------------------------------------

## Applying saved weights to new data

If you receive new data from the same assay platform, apply the
previously estimated scale parameters without refitting:

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

------------------------------------------------------------------------

## Weight diagnostics

[`weight_diagnostics()`](https://immunoplex.github.io/curveRweights/reference/weight_diagnostics.md)
summarises the distribution and effective information content of a
weight vector:

``` r
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

- **n_eff**: effective sample size — equals n when all weights are
  uniform, decreases as weights become more heterogeneous.
- **eff_ratio**: n_eff / n_valid — 1.0 = uniform, lower = more
  differentiation.
- **weight_ratio**: max / min weight — large values indicate aggressive
  downweighting of some observations.
- **gini**: Gini coefficient — 0 = uniform, approaching 1 =
  concentrated.

------------------------------------------------------------------------

## Prior sensitivity

The default priors are:

- gamma_0 ~ Normal(0, 1): phi centered at 1, 95% interval \[0.14, 7.4\]
- gamma_1 ~ Normal(1, 0.5): beta1 centered at the delta-method
  prediction of 1, 95% interval \[0, 2\]

For groups where beta1 is far from 1 (e.g. prn with beta1 near 2), the
prior exerts meaningful shrinkage. Test sensitivity by comparing results
with wider priors:

``` r
# Default prior
sw_default <- fit_saturated_weight(
  dat_prn, cell_col = "cell", pcov_col = "pcov", plate_col = "plate"
)

# Wider prior
sw_wide <- fit_saturated_weight(
  dat_prn, cell_col = "cell", pcov_col = "pcov", plate_col = "plate",
  prior_gamma1 = brms::set_prior("normal(1, 1.0)", class = "b",
                                  dpar = "sigma")
)

cat("Default beta1:", sw_default$beta1, "\n")
cat("Wide    beta1:", sw_wide$beta1, "\n")
```

If beta1 shifts materially, the data are not fully informative and the
prior is doing real work. The wider prior gives weights that better
reflect the data; the narrower prior gives more conservative
(closer-to-uniform) weights.

------------------------------------------------------------------------

## Session information

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 26100)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.1         curveRweights_0.2.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6          tensorA_0.36.2.1      QuickJSR_1.9.2       
#>  [4] xfun_0.57             bslib_0.11.0          ggplot2_4.0.3        
#>  [7] htmlwidgets_1.6.4     inline_0.3.21         lattice_0.22-7       
#> [10] vctrs_0.7.3           tools_4.5.2           generics_0.1.4       
#> [13] curl_7.1.0            stats4_4.5.2          parallel_4.5.2       
#> [16] sandwich_3.1-1        tibble_3.3.1          pkgconfig_2.0.3      
#> [19] brms_2.23.0           Matrix_1.7-4          checkmate_2.3.4      
#> [22] RColorBrewer_1.1-3    S7_0.2.2              desc_1.4.3           
#> [25] distributional_0.8.1  RcppParallel_5.1.11-2 lifecycle_1.0.5      
#> [28] compiler_4.5.2        farver_2.1.2          stringr_1.6.0        
#> [31] textshaping_1.0.5     Brobdingnag_1.2-9     codetools_0.2-20     
#> [34] htmltools_0.5.9       sass_0.4.10           bayesplot_1.15.0     
#> [37] yaml_2.3.12           pillar_1.11.1         pkgdown_2.2.0        
#> [40] jquerylib_0.1.4       MASS_7.3-65           cachem_1.1.0         
#> [43] StanHeaders_2.32.10   bridgesampling_1.2-1  abind_1.4-8          
#> [46] multcomp_1.4-30       nlme_3.1-168          rstan_2.32.7         
#> [49] posterior_1.7.0       tidyselect_1.2.1      digest_0.6.39        
#> [52] mvtnorm_1.3-5         stringi_1.8.7         splines_4.5.2        
#> [55] fastmap_1.2.0         grid_4.5.2            cli_3.6.6            
#> [58] magrittr_2.0.5        loo_2.9.0.9000        pkgbuild_1.4.8       
#> [61] survival_3.8-3        TH.data_1.1-5         scales_1.4.0         
#> [64] backports_1.5.1       estimability_1.5.1    rmarkdown_2.31       
#> [67] matrixStats_1.5.0     emmeans_2.0.2         otel_0.2.0           
#> [70] gridExtra_2.3         ragg_1.5.1            zoo_1.8-15           
#> [73] coda_0.19-4.1         evaluate_1.0.5        knitr_1.51           
#> [76] V8_8.0.1              rstantools_2.6.0      rlang_1.2.0          
#> [79] Rcpp_1.1.1-1.1        xtable_1.8-8          glue_1.8.1           
#> [82] jsonlite_2.0.0        R6_2.6.1              systemfonts_1.3.2    
#> [85] fs_2.1.0
```
