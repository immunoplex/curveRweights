# curveRweights

Bayesian Precision Weighting from Calibration Curve Uncertainty

## The Problem

Immunoassay measurements are determined via a calibration curve (typically a four-parameter logistic, 4PL).  Each observation has a posterior coefficient of variation (pcov) that quantifies how precisely the curve determines that concentration.  Observations near the limits of quantification (LLOQ, ULOQ) have large pcov values because the curve is flat there; midrange observations have small pcov because the curve is steep.

Traditional approaches handle this imprecision with a **binary gate**: observations between LLOQ and ULOQ are included with equal weight, and all others are excluded.  This discards information (excluded observations still carry signal) and treats all included observations as equally precise (they are not).

## The Solution

`curveRweights` estimates the **continuous, power-law relationship** between pcov and residual variance using a joint Bayesian location-scale model:

```
  sigma_i = phi * pcov_i^beta1
  w_i     = 1 / sigma_i^2
```

where:

- **phi** (baseline scaling): captures variance inflation beyond what pcov alone predicts --- biological scatter, plate effects, matrix interference.  When phi = 1, pcov is a perfectly calibrated proxy for residual SD.  Values > 1 indicate excess variance (analogous to positive tau-squared in classical meta-analysis).

- **beta1** (precision exponent): captures how steeply pcov drives weight differentiation.  beta1 = 1 is the theoretical delta-method prediction (sigma proportional to pcov).  beta1 > 1 means the assay's precision gradient is steeper than theory predicts (amplified weighting); beta1 < 1 means it's shallower (compressed weighting).

The estimation uses a **saturated cell-means location model** (one coefficient per experimental condition) so that the scale parameters are identified purely from within-cell residual structure and are not contaminated by misspecified location effects.

## Relationship to Meta-Analytic Tau-Squared

In a classical random-effects meta-regression (DerSimonian & Laird, 1986), total per-observation variance is decomposed as `tau2 + vi` where `vi` is measurement error and `tau2` is residual heterogeneity.  In `curveRweights`:

- `phi^2` plays the role of **tau2**: global variance inflation beyond measurement error
- `pcov_i^(2*beta1)` plays the role of **vi**: per-observation measurement-error variance

The key difference is that phi scales the entire sigma curve multiplicatively rather than adding a fixed floor, and beta1 is data-estimated rather than fixed.  When using these weights in `metafor::rma()`, set `method = "FE"` (tau2 = 0) because heterogeneity is already captured by phi.

## Installation

```r
# Install from GitHub (once published):
# remotes::install_github("yourusername/curveRweights")

# Or install from local source during development:
# In RStudio: Build > Install
# Or from the command line:
devtools::install()
```

### Dependencies

`curveRweights` requires:

- **R >= 4.1.0**
- **brms >= 2.19.0** (Bayesian regression via Stan)
- **dplyr >= 1.1.0**
- **tibble**

Optional (for downstream analysis):

- **metafor** (for `rma()` meta-regression)
- **survey** (for `svyranktest()` weighted rank tests)

## Quick Start

```r
library(curveRweights)

# Step 1: Define the saturated cell-means factor EXTERNALLY
# This gives you full control over what "cells" means in your experiment.
dat$cell <- interaction(dat$Arm, dat$timeperiod, drop = TRUE)

# Step 2: Fit per (antigen, source) group
batch <- fit_saturated_weight_batch(
  datg       = dat,
  group_vars = c("antigen", "source"),
  cell_col   = "cell",
  pcov_col   = "pcov",
  plate_col  = "plate",
  iter       = 4000,
  warmup     = 1000,
  chains     = 4,
  cores      = 4
)

# Step 3: Inspect scale estimates
batch$scale_table

# Step 4: Use weights for any specific comparison
dat_comparison <- batch$data |>
  dplyr::filter(Arm %in% c("TT", "TdaP"), timeperiod == "post3rd")

# dat_comparison$w_saturated is ready for metafor::rma() or survey::svyranktest()
```

## Core Functions

| Function | Purpose |
|---|---|
| `prepare_cv()` | Compute yi, cv_i, log_cv from concentration + pcov |
| `diagnose_cv()` | Check whether pcov has enough variation for beta1 estimation |
| `interpret_beta1()` | Classify estimated beta1 into precision-weighting regimes |
| `compute_saturated_weights()` | Deterministic weight computation from (gamma_0, gamma_1) |
| `weight_diagnostics()` | n_eff, weight_ratio, Gini coefficient summary |
| `fit_saturated_weight()` | Core: fit the joint location-scale brms model |
| `fit_saturated_weight_batch()` | Fit across multiple groups (antigen x source) |
| `apply_saturated_weights()` | Apply saved scale estimates to new data without re-fitting |

## Why Saturated Cell Means?

The location model uses `yi ~ 0 + cell` (no intercept; one coefficient per cell) rather than a structured model like `yi ~ Arm * Timeperiod`.  This is because:

1. **Any misspecification in the location model** (e.g., assuming additivity when there's an interaction) would leave systematic patterns in the residuals that contaminate the scale estimates.

2. **The saturated model absorbs ALL systematic location variation by construction**, so residuals are the cleanest possible input to the scale model.

3. **The cell factor is created externally**, giving you full control.  You can define it to match your experimental design exactly.

## Why the Cell Factor Is External

The `cell_col` parameter takes a pre-existing factor column rather than `Arm` and `Timeperiod` separately because:

- Different experiments have different designs (2 arms x 3 times, 4 arms x 1 time, dose-response x batch, etc.)
- The function doesn't need to know what the cells represent --- it only needs one coefficient per unique experimental condition
- You can verify the cell factor has the right levels before passing it to a long MCMC run

## Prior Sensitivity

The default priors are:

- **gamma_0 ~ Normal(0, 1)**: phi centered at 1, 95% interval [0.14, 7.4]
- **gamma_1 ~ Normal(1, 0.5)**: beta1 centered at theory-predicted 1, 95% interval [0, 2]

For groups where beta1 is far from 1 (e.g., > 2), the prior exerts meaningful shrinkage.  The `prior_gamma1` parameter allows you to pass a custom prior for sensitivity analysis:

```r
# Wider prior for groups with large beta1
sw <- fit_saturated_weight(
  df = dat,
  cell_col = "cell",
  prior_gamma1 = brms::prior(normal(1, 1.0), class = "b", dpar = "sigma")
)
```

## References

- Burkner P-C (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. *Journal of Statistical Software*, 80(1), 1-28. DOI: 10.18637/jss.v080.i01

- Burkner P-C (2018). Advanced Bayesian Multilevel Modeling with the R Package brms. *The R Journal*, 10(1), 395-411. DOI: 10.32614/RJ-2018-017

- Carpenter B, et al. (2017). Stan: A Probabilistic Programming Language. *Journal of Statistical Software*, 76(1), 1-32. DOI: 10.18637/jss.v076.i01

- DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. *Controlled Clinical Trials*, 7(3), 177-188. PMID: 3802833.

- Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a meta-analysis. *Statistics in Medicine*, 21(11), 1539-1558. PMID: 12111919.

- Viechtbauer W (2010). Conducting Meta-Analyses in R with the metafor Package. *Journal of Statistical Software*, 36(3), 1-48. DOI: 10.18637/jss.v036.i03

## License

AGPL-3.0
