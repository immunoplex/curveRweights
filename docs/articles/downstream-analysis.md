# Using the weights downstream: which analysis for which question

Read
[`vignette("measurement-error")`](https://immunoplex.github.io/curveRweights/articles/measurement-error.md)
first for what $`\tau^2`$, $`v_i`$, the weight and the reliability mean.
This page is organised by **what you want to do**.

## Pick the recipe

| You want to… | Noisy variable is the… | Use |
|----|----|----|
| Compare groups or time points | outcome | [meta-regression](#metafor) with `v_meas` (recommended), or [brms](#brms) |
| Compare groups when some samples are below/above the curve | outcome, censored | [censored model](#brms) or [rank test](#rank) |
| Relate antibody to infection, or one antigen to another | predictor | [regression calibration](#predictor) |
| Report how correlated two assays *really* are | both | [disattenuation](#predictor) |
| Compare two groups without a model | outcome | [rank test](#rank) |

Three rules apply to every recipe.

1.  **Subset first, then normalise.** Weights are scaled to mean 1 over
    the rows they were computed on. After you filter to one arm or one
    time point, call
    [`normalize_weights()`](https://immunoplex.github.io/curveRweights/reference/normalize_weights.md)
    on the rows you actually analyse.
2.  **Do not weight censored samples.** Represent them by their limit
    (`censor_type`, `bound`), not by `w`.
3.  **Pass the right variance to the right function.** Software that
    estimates $`\tau^2`$ itself (`metafor`) needs the *measurement*
    variance `v_meas`. Software that expects a weight (`lm(weights =)`,
    `brms`) needs `w_norm`. Passing the total variance to `metafor`
    would count $`\tau^2`$ twice.

## A small data set

The examples below simulate what the weight table looks like for 150
people in two arms, each sampled at two visits, so everything runs
without a database. In practice the table comes from I-SPI (Weights tab,
Data tab download) or from `fit_precision_weights(model = "additive")`
followed by
[`join_weights()`](https://immunoplex.github.io/curveRweights/reference/join_weights.md).

``` r
set.seed(11)
n_pat <- 150
d <- data.frame(patient = rep(seq_len(n_pat), each = 2),
                arm     = rep(c("A", "B"), each = n_pat))
d$obs <- seq_len(nrow(d))
tau <- 0.4
b_pat <- rnorm(n_pat, 0, tau)
true  <- 1 + 0.30 * (d$arm == "B") + b_pat[d$patient]
d$se_meas <- exp(rnorm(nrow(d), log(0.2), 0.8))        # curve SE per sample
d$yi      <- true + rnorm(nrow(d), 0, d$se_meas)       # observed log10 value
d$censor_type <- "none"

d <- compute_additive_weights(d, tau2 = tau^2 + 0.0, kappa = 1,
                              censor_col = "censor_type")
head(d[, c("patient", "arm", "yi", "se_meas", "v_meas", "w", "w_norm", "reliability")], 4)
#>   patient arm        yi    se_meas      v_meas        w    w_norm reliability
#> 1       1   A 0.6873140 0.05578947 0.003112464 6.130739 1.3644421   0.9809183
#> 2       1   A 1.3692189 0.34212506 0.117049557 3.609463 0.8033130   0.5775140
#> 3       2   A 0.6939739 0.28115219 0.079046557 4.183286 0.9310217   0.6693257
#> 4       2   A 1.6282186 0.34723252 0.120570420 3.564168 0.7932323   0.5702668
```

## Recipe 1: meta-regression with `metafor`

[`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html)
is the cleanest tool when the outcome is noisy and each observation has
a known measurement variance. It estimates the extra variance (between
samples $`\tau^2`$; between patients) by REML and applies the correct
weights, so you supply only the measurement variance.

``` r
library(metafor)
#> Warning: package 'metafor' was built under R version 4.5.3
#> Loading required package: Matrix
#> Loading required package: metadat
#> Warning: package 'metadat' was built under R version 4.5.3
#> Loading required package: numDeriv
#> 
#> Loading the 'metafor' package (version 4.8-0). For an
#> introduction to the package please type: help(metafor)
fit <- rma.mv(yi = yi, V = meas_variance(d), mods = ~ arm,
              random = list(~ 1 | patient, ~ 1 | obs),
              data = d, method = "REML")
round(coef(summary(fit)), 3)
#>         estimate    se   zval pval ci.lb ci.ub
#> intrcpt    0.964 0.044 21.899    0 0.878 1.050
#> armB       0.336 0.062  5.393    0 0.214 0.458
```

- `V = meas_variance(d)` is $`v_i`$ (measurement only).
- `~ 1 | obs` (one level per row) is the between-sample variance, the
  model’s $`\tau^2`$; `~ 1 | patient` is a patient-level variance for
  repeated measures.
- The coefficient `armB` is the group difference on the log10 scale.
  Multiply $`10^{\text{estimate}}`$ for a fold change.

Compare with ordinary least squares ignoring the weights:

``` r
fit_ols <- lm(yi ~ arm, data = d)
round(coef(summary(fit_ols)), 3)
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)    0.964      0.043  22.627        0
#> armB           0.350      0.060   5.813        0
```

Both estimate the same difference; the weighted analysis is generally
the more precise one, and it does not let a few very well-measured
samples drive the answer because their weight is capped at $`1/\tau^2`$.

**Do not** use `method = "FE"` with these weights. Older documentation
advised that for the multiplicative model; it assumes there is no
between-sample variance beyond the assay error, which is rarely true for
people.

## Recipe 2: `brms`, censored outcomes, repeated measures

`brms` fits the full model and lets censored samples contribute through
their limits. The same structure the weights come from:

``` r
library(brms)
d_q <- d[d$censor_type == "none", ]          # quantified rows get weights
d_q$w_norm <- normalize_weights(d_q)$w_norm  # re-normalise on the analysis subset

# 1. quantified-only, weighted:
f1 <- brm(yi | weights(w_norm) ~ arm + (1 | patient), data = d_q)

# 2. all rows, with censoring and measurement-error SDs (preferred):
#    d_cens is the output of prepare_censoring(); cens_flag is 0 quantified,
#    -1 left, 1 right, 2 interval. Not run here (needs Stan).
f2 <- brm(bf(yi | cens(cens_flag, yi2) ~ arm + (1 | patient),
             nlf(sigma ~ log(sqrt(exp(2 * ltau) + se_meas^2))),
             ltau ~ 1),
          data = d_cens)
```

Model 2 matches the additive likelihood: total SD
$`\sqrt{\tau^2 + se^2}`$, with censored rows contributing probabilities.
Use it when more than a few percent of your samples are out of range.

## Recipe 3: when the noisy variable is a predictor

Weights cannot fix this. Suppose you regress an immune outcome on
antibody level: noise in the antibody level pulls the slope towards
zero.

``` r
set.seed(5)
n  <- 3000
x  <- rnorm(n, 0, 0.4)                              # true level (log10)
se <- exp(rnorm(n, log(0.25), 0.5))
x_obs <- x + rnorm(n, 0, se)                        # measured level
y  <- 1.0 * x + rnorm(n, 0, 0.3)                    # true slope = 1

p <- data.frame(se_meas = se, x_obs = x_obs, cell = "all")
p <- compute_additive_weights(p, tau2 = 0.4^2, kappa = 1)
x_cal <- reliability_correct(p$x_obs, p$cell, p$reliability)

c(true      = 1,
  naive     = unname(coef(lm(y ~ x_obs))[2]),
  corrected = unname(coef(lm(y ~ x_cal))[2]))
#>      true     naive corrected 
#> 1.0000000 0.6014975 0.9842421
```

The naive slope is biased toward 0 by roughly the mean reliability;
regression calibration
([`reliability_correct()`](https://immunoplex.github.io/curveRweights/reference/reliability_correct.md))
recovers the true value. For publication, note two caveats: standard
errors from the second stage ignore the uncertainty in the calibration
step (bootstrap it, or use
[`brms::me()`](https://paulbuerkner.com/brms/reference/me.html) for a
single joint model), and the correction assumes the measurement error is
independent of the outcome.

For a **correlation** between two assays (for example IgG1 to PT versus
IgG1 to FHA) use
[`attenuation_correct()`](https://immunoplex.github.io/curveRweights/reference/attenuation_correct.md):

``` r
attenuation_correct(r = 0.42, rel_x = p$reliability, rel_y = p$reliability)
#> [1] 0.6127692
```

Report the observed correlation and the corrected one together; a
corrected correlation above 1 in absolute value means the reliabilities
are too optimistic or the sample is small, and the function warns.

## Recipe 4: rank tests with censored data

If you prefer not to assume a distribution, the censoring-aware rank
tests come from survival analysis, where censored observations are the
norm. For **left-censored** assay data, flip the sign so the censoring
becomes right censoring, then use a log-rank-type test with Peto-Peto
(Prentice-Wilcoxon) weights, `rho = 1`, a close relative of the Gehan
test that weights early differences more:

``` r
library(survival)
set.seed(8)
r <- data.frame(arm = rep(c("A", "B"), each = 120))
r$y     <- rnorm(240, 1 + 0.4 * (r$arm == "B"), 0.5)
lod     <- 0.7                                       # lower limit, log10
r$cens  <- r$y < lod                                  # below the limit
r$y_use <- ifelse(r$cens, lod, r$y)                   # replace by the limit
M <- max(r$y_use) + 1
r$t     <- M - r$y_use                                # flip: left becomes right
r$event <- as.integer(!r$cens)                        # 1 = quantified (observed)
survdiff(Surv(t, event) ~ arm, data = r, rho = 1)
#> Call:
#> survdiff(formula = Surv(t, event) ~ arm, data = r, rho = 1)
#> 
#>         N Observed Expected (O-E)^2/E (O-E)^2/V
#> arm=A 120     42.7     69.8      10.5      38.3
#> arm=B 120     72.9     45.8      16.0      38.3
#> 
#>  Chisq= 38.3  on 1 degrees of freedom, p= 6e-10
```

This replaces the older advice of
[`survey::svyranktest()`](https://rdrr.io/pkg/survey/man/svyranktest.html)
with weights: that function treats the weights as sampling weights and
ignores censoring entirely. A rank test does not use `w` at all; it
makes no use of the measurement-error information, which makes it a
useful **sensitivity check** for the model-based results above.

## Checklist before you publish

- State the model used (`model` column) and the curveRweights version
  (`weights_pkg_version`).
- State how many samples were censored (`n_censored`) and how they were
  handled.
- Report group reliabilities; if they are below roughly 0.5 say that
  individual values are not reliable.
- Run one sensitivity analysis that uses no weights, and one that uses a
  different censoring treatment (rank test or Tobit); conclusions should
  agree.
- Do not mix weights from different runs in one analysis: the stored
  weight set for a sample is replaced by the most recent run.

## References

- Viechtbauer W. Conducting meta-analyses in R with the metafor package.
  *J Stat Softw* 2010;36(3).
- Bürkner P-C. brms: An R package for Bayesian multilevel models using
  Stan. *J Stat Softw* 2017;80(1).
- Peto R, Peto J. Asymptotically efficient rank invariant test
  procedures. *J R Stat Soc A* 1972;135:185-207.
- Rosner B, Spiegelman D, Willett WC. Correction of logistic regression
  relative risk estimates and confidence intervals for measurement
  error. *Am J Epidemiol* 1990;132:734-45. (Regression calibration.)
- Hughes MD. Analysis and design issues for studies using censored
  biomarker measurements. *Biometrics* 2000;56:445-52.

See also the calibration pages: curveRfreq
(<https://immunoplex.github.io/curveRfreq/>) and curveRbayes
(<https://immunoplex.github.io/curveRbayes/>).
