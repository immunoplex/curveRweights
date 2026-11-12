# Measurement error in assay data: what the weights mean

**Who this is for.** You know linear models, random effects and
likelihoods, and you know your assay. You have not necessarily met
*measurement-error models*, the branch of statistics that separates “the
sample really differs” from “the assay read it imprecisely”. This page
explains the one idea you need, how curveRweights uses it, and which
numbers it hands you. A companion page,
[`vignette("downstream-analysis")`](https://immunoplex.github.io/curveRweights/articles/downstream-analysis.md),
shows what to do with those numbers.

## The one idea: every reading is signal plus noise

Take one serum sample measured on a bead assay. The number you report,
the back-calculated log10 concentration $`y_i`$, is built from three
things:

``` math

y_i \;=\; \underbrace{\mu_{\text{group}}}_{\text{the effect you study}}
\;+\; \underbrace{b_i}_{\text{real biological spread}}
\;+\; \underbrace{e_i}_{\text{assay measurement error}}
```

- $`\mu_{\text{group}}`$ is the group mean, for example a vaccine arm at
  one time point. Estimating differences in $`\mu`$ is usually the
  scientific question.
- $`b_i`$ is **true biological variation** between people (or between
  samples): some people simply respond more strongly. It has variance
  $`\tau^2`$. More assay precision never removes it.
- $`e_i`$ is **measurement error**. The calibration curve tells us how
  large it is for *this* sample: the standard error `se` of the
  back-calculated concentration, so
  $`\text{Var}(e_i) = v_i = \kappa^2\,se_i^{2\beta_1}`$. Near the top
  and bottom of the standard curve the curve is flat, a small change in
  the signal means a large change in concentration, and $`v_i`$ is
  large. In the middle $`v_i`$ is small.

Because $`b_i`$ and $`e_i`$ are independent, their variances **add**:

``` math

\sigma_i^2 \;=\; \tau^2 + v_i \;=\; \tau^2 + \kappa^2\, se_i^{2\beta_1}.
```

This is the *additive* model. It is the same structure as a
random-effects meta-analysis, where each study has a within-study
variance $`v_i`$ (known) and the studies differ by a between-study
variance $`\tau^2`$ (estimated).

### Why we weight at all

A weighted analysis gives precisely measured samples more say. The
efficient weight is the inverse of the *total* variance:

``` math

w_i = \frac{1}{\tau^2 + v_i}.
```

Look at what this does at the extremes:

- A hopelessly noisy sample ($`v_i \to \infty`$) gets weight near 0.
  Good.
- A perfectly measured sample ($`v_i \to 0`$) gets weight $`1/\tau^2`$,
  not infinity. Good: a perfect measurement of *one person* still tells
  you only about that person, and people differ by $`\tau^2`$.

So the **largest weight ratio you can ever see is bounded**:
$`w_{max}/w_{min} \le (\tau^2 + v_{max})/\tau^2`$.

## Additive versus multiplicative

Versions of curveRweights before 0.3.0 used a *multiplicative* model,
$`\sigma_i = \phi\, se_i^{\beta_1}`$, so
$`\sigma_i^2 = \phi^2 se_i^{2\beta_1}`$. Here **all** the variance is
attributed to the measurement. A sample with a very small $`se`$ then
gets a very large weight, with no ceiling, even when real
person-to-person variation (the $`\tau^2`$ above) would limit what that
sample can tell you. The additive model puts a floor under $`\sigma_i`$
and therefore a ceiling on $`w_i`$.

The simulation below uses the same $`se_i`$ values for both. The true
model is additive with $`\tau = 0.4`$; we ask each approach what weights
it implies.

``` r
set.seed(42)
se  <- exp(rnorm(2000, log(0.15), 0.9))          # a realistic spread of SEs
tau <- 0.4

w_add  <- 1 / (tau^2 + se^2)                      # additive: ceiling 1/tau^2 = 6.25
w_mult <- 1 / (1.1^2 * se^2)                      # multiplicative, phi = 1.1

ratio <- function(w) max(w) / min(w)
c(additive = ratio(w_add), multiplicative = ratio(w_mult))
#>       additive multiplicative 
#>       90.15707   274173.78492
c(`additive ceiling 1/tau^2` = 1 / tau^2,
  `multiplicative max weight` = max(w_mult))
#>  additive ceiling 1/tau^2 multiplicative max weight 
#>                      6.25                  15878.96
```

With the SE spread used here the multiplicative weights span a ratio of
hundreds of thousands, against about 90 for the additive model, and its
largest weight (about 16,000) is thousands of times the additive ceiling
of 6.25. In a real analysis such a sample effectively *is* the group
mean. That is the practical reason the additive model is the one to
prefer when you analyse **differences between people or groups** and
when between-subject variation is not negligible.

Which model do I get?
[`fit_precision_weights()`](https://immunoplex.github.io/curveRweights/reference/fit_precision_weights.md)
still defaults to `model = "multiplicative"` so that existing pipelines
and results are unchanged; ask for the new one with
`model = "additive"`. In I-SPI the choice is recorded with every weight
set (column `model`).

| Situation | Prefer |
|----|----|
| Comparing groups or time points; people differ a lot | additive |
| Assay noise dominates and people barely differ ($`\tau \approx 0`$) | either; they agree |
| You must reproduce an earlier analysis | the model you used then |

## Reliability: how much of what you see is real

The additive model gives each sample a **reliability**

``` math

\lambda_i \;=\; \frac{\tau^2}{\tau^2 + v_i} \;\in\; (0, 1].
```

It is the fraction of that sample’s variance that is genuine biological
signal. If you know the *intraclass correlation* (ICC) from
repeated-measures work or from the Cronbach’s $`\alpha`$ of a
questionnaire, this is the same idea applied sample by sample: 0.9 means
“mostly signal”, 0.3 means “mostly assay noise”.

``` r
d <- data.frame(se_meas = c(0.03, 0.10, 0.30, 1.00))
out <- compute_additive_weights(d, tau2 = 0.16, kappa = 1)
out[, c("se_meas", "v_meas", "w", "reliability", "w_norm")]
#>   se_meas v_meas        w reliability    w_norm
#> 1    0.03  9e-04 6.215040   0.9944065 1.4658579
#> 2    0.10  1e-02 5.882353   0.9411765 1.3873914
#> 3    0.30  9e-02 4.000000   0.6400000 0.9434261
#> 4    1.00  1e+00 0.862069   0.1379310 0.2033246
```

Reliability is the number to **show** to collaborators. A weight depends
on an arbitrary scale; a reliability does not. It also explains
regression dilution, below.

### Weights fix efficiency, not bias in predictors

Weights improve the precision of an estimate when the noisy variable is
the **outcome**. When the noisy variable is a **predictor** (for example
antibody level used to predict infection, or one antigen’s response
correlated with another’s), weights do *not* repair it. Noise in a
predictor flattens the fitted slope toward zero by roughly the factor
$`\lambda`$ (*regression dilution* or *attenuation*). Two remedies are
included:

``` r
# a correlation of 0.42 between two noisy assays with mean reliabilities 0.85, 0.70
attenuation_correct(0.42, rel_x = 0.85, rel_y = 0.70)
#> [1] 0.5444911
```

[`reliability_correct()`](https://immunoplex.github.io/curveRweights/reference/reliability_correct.md)
replaces each measured value by its best prediction of the true value
(shrinking it toward its cell mean by $`\lambda_i`$), which removes the
first-order bias of a regression slope. See
[`vignette("downstream-analysis")`](https://immunoplex.github.io/curveRweights/articles/downstream-analysis.md)
for a worked example.

## Samples outside the standard curve: censoring

A sample whose signal falls below the lowest standard is not “a noisy
concentration”. It is a sample we know only to be *below some limit*.
Giving such a value an ordinary weight is wrong: its precision is poor
*because* the concentration is extreme, so inverse-variance weights
down-weight the very samples that say “low”, and group means drift
towards the middle of the curve.

The additive model treats these samples as **censored**:

- left-censored: the true value is below the lower limit;
- right-censored: above the upper limit;
- interval: between two limits.

Quantified samples contribute the usual Normal density. A censored
sample contributes the *probability* of being beyond its limit, so it
informs the group mean and $`\tau`$ without being treated as a precise
point. In the weight table these rows carry `censor_type` and `bound`,
and their `w_norm` is `NA` on purpose: represent them by the bound, not
by a weight.

``` r
set.seed(3)
n   <- 20000
mu  <- c(low = 0.6, high = 2.0)
g   <- rep(names(mu), each = n)
tru <- rnorm(2 * n, rep(mu, each = n), 0.4)
se_fun <- function(y) 0.05 + 0.35 * plogis((0.5 - y) / 0.15) +
                              0.35 * plogis((y - 3.0) / 0.15)
obs <- tru + rnorm(2 * n, 0, se_fun(tru))
w   <- 1 / (0.4^2 + se_fun(obs)^2)               # additive weights, own-SE
rbind(truth      = mu,
      unweighted = tapply(obs, g, mean)[names(mu)],
      weighted   = tapply(obs * w, g, sum)[names(mu)] / tapply(w, g, sum)[names(mu)])
#>                  low     high
#> truth      0.6000000 2.000000
#> unweighted 0.5964011 2.002382
#> weighted   0.6928282 1.995671
```

Even with the additive ceiling the weighted mean of the low group (about
0.69) sits well above the truth (0.60), while the unweighted mean (0.60)
does not; the high group is barely affected because most of it lies
where precision is good. That is the reason the package treats
out-of-range samples by their limits instead of weights.
[`prepare_censoring()`](https://immunoplex.github.io/curveRweights/reference/prepare_censoring.md)
builds the flags from curveRcore’s `censor_type`, `censor_lower` and
`censor_upper`.

## Repeated measures and dilutions

If the same patient contributes several samples, their $`b_i`$ are not
independent.
[`fit_additive_weight()`](https://immunoplex.github.io/curveRweights/reference/fit_additive_weight.md)
accepts a patient column and adds a patient random intercept, the
familiar `(1 | patientid)`, so the repeated measures are not
over-counted. Dilutions or replicates of the *same sample* are handled
with `sample_col`.

Note on units: the weights are on the **log10 concentration** scale, the
scale the model is fitted on. When censoring is used, the model puts the
response on the *final concentration* scale (it adds
$`\log_{10}(\text{dilution})`$), the same scale as the censoring limits
curveRcore reports.

## What the stored columns mean

When I-SPI or the i-spi-compute worker saves weights, each sample row
(table `calib_weights`) holds:

| Column | Meaning |
|----|----|
| `model` | `multiplicative` or `additive` |
| `se` | the calibration curve’s SE for that sample (measurement SD, log10) |
| `sigma` | total SD $`\sigma_i`$ |
| `w` | raw inverse-variance weight $`1/\sigma_i^2`$ |
| `w_norm` | `w` divided by its mean over quantified rows (mean 1). `NA` for censored rows in the additive model |
| `v_meas` | measurement variance $`v_i`$ (additive only). **Use this** for `metafor` |
| `tau2` | between-sample variance $`\tau^2`$ (additive only) |
| `reliability`, `_lo`, `_hi` | $`\lambda_i`$ and a 95% interval (additive only) |
| `censor_type`, `bound` | none / left / interval / right, and the limit (log10) |
| `weights_pkg_version` | the curveRweights version that produced the row |

Per curve group (`calib_weights_fit`) you also get $`\tau`$, $`\kappa`$,
the group mean reliability, `n_censored`, whether $`\beta_1`$ was fixed,
and an *identifiability* note. Column definitions live in code too:
[`curveRcore::weights_contract()`](https://immunoplex.github.io/curveRcore/reference/weights_contract.html).

## How good is $`\tau`$? Identifiability

$`\tau`$ (signal) and $`\kappa`$ (the scale of the measurement error)
can be separated only because the $`se_i`$ differ from sample to sample:
when $`se`$ is large, the total variance is mostly measurement; when
small, mostly $`\tau^2`$. If all $`se_i`$ are nearly equal, the data
cannot tell the two apart. By default $`\beta_1`$ is fixed at 1 (the
theoretical value from the delta method) and $`\kappa`$ is estimated.
The fit reports a message when $`\kappa`$ is weakly identified; treat
$`\tau`$ and the reliabilities from such a group as provisional, or fix
$`\kappa = 1`$ and assume `se` is the true measurement SD.

## Glossary

| Term | Plain meaning |
|----|----|
| **Measurement error** | Difference between what the assay reports and the true value in that sample. |
| **$`se`$ (standard error of concentration)** | The curve’s own estimate of the measurement SD for one sample. |
| **Signal variance $`\tau^2`$** | True person-to-person (sample-to-sample) variance in the quantity of interest. |
| **Total variance $`\sigma^2`$** | $`\tau^2 + v_i`$. |
| **Weight $`w`$** | $`1/\sigma^2`$; precisely measured samples count more, but never more than $`1/\tau^2`$. |
| **Reliability $`\lambda`$** | $`\tau^2/\sigma^2`$; fraction of the variance that is signal. Like an ICC. |
| **Attenuation / regression dilution** | A noisy predictor makes the fitted slope too flat by about $`\lambda`$. |
| **Censoring** | We know only that the value is below/above a limit. |
| **Left / right / interval** | Below the lower limit / above the upper limit / between two limits. |
| **Identifiability** | Whether the data contain enough information to tell two quantities apart. |
| **REML** | Restricted maximum likelihood; the standard estimator for variance components. |
| **Delta method** | A first-order approximation that converts the curve’s uncertainty into an SE for the concentration. |

## Where next

- [`vignette("downstream-analysis")`](https://immunoplex.github.io/curveRweights/articles/downstream-analysis.md):
  which analysis to run with these numbers.
- [`vignette("precision-weighting")`](https://immunoplex.github.io/curveRweights/articles/precision-weighting.md):
  the multiplicative workflow in detail.
- curveRfreq and curveRbayes: how the calibration curves, standard
  errors and censoring labels are produced
  (<https://immunoplex.github.io/curveRfreq/>,
  <https://immunoplex.github.io/curveRbayes/>).
- The page-by-page table of stored columns:
  [`curveRcore::weights_contract()`](https://immunoplex.github.io/curveRcore/reference/weights_contract.html).

## Further reading

- Carroll RJ, Ruppert D, Stefanski LA, Crainiceanu CM. *Measurement
  Error in Nonlinear Models*, 2nd ed. Chapman & Hall/CRC, 2006. The
  standard textbook.
- Viechtbauer W. Conducting meta-analyses in R with the metafor package.
  *J Stat Softw* 2010;36(3).
- Hughes MD. Analysis and design issues for studies using censored
  biomarker measurements. *Biometrics* 2000;56:445-52. The case for
  handling values below detection limits by likelihood.
- Frost C, Thompson SG. Correcting for regression dilution bias:
  comparison of methods for a single predictor variable. *J R Stat Soc
  A* 2000;163:173-89.
