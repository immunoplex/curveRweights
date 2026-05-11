# curveRweights: Bayesian Precision Weighting from Calibration Curve Uncertainty

Immunoassay observations are measured via a calibration curve (typically
a nonlinear, sigmoidal function of concentration and assay response
following a logistic, log-logistic or Gompertz form). Each observation
has a posterior coefficient of variation (pcov) that quantifies how
precisely the curve determines that concentration. `curveRweights`
estimates the power-law relationship between pcov and residual variance
using a joint Bayesian location-scale model:

\$\$\sigma_i = \phi \cdot \mathrm{pcov}\_i^{\beta_1}\$\$ \$\$w_i = 1 /
\sigma_i^2\$\$

where:

- phi:

  Baseline scaling factor (exp of the log-scale intercept). Values \> 1
  indicate excess variance beyond what pcov predicts — biological
  scatter, plate effects, or matrix interference.

- beta1:

  Precision exponent. beta1 = 1 means pcov is a direct proxy for
  residual SD (the delta-method prediction). beta1 \> 1 means amplified
  weighting; beta1 \< 1 means compressed weighting.

## Relationship to meta-analytic tau-squared

In a classical random-effects meta-regression (DerSimonian & Laird,
1986), total per-observation variance is decomposed as \\\tau^2 + v_i\\
where \\v_i\\ is measurement error and \\\tau^2\\ is residual
heterogeneity. In the curveRweights model, \\\sigma_i^2\\ plays the role
of \\(\tau^2 + v_i)\\ simultaneously: \\\phi^2\\ captures what
\\\tau^2\\ captured globally, and \\\mathrm{pcov}\_i^{2\beta_1}\\
captures what \\v_i\\ captured locally. Consequently, when using these
weights in
[`rma`](https://wviechtb.github.io/metafor/reference/rma.uni.html), set
`method = "FE"` (tau2 fixed at 0) to avoid double-counting
heterogeneity.

## Independence of scale and location estimates

The Fisher information for the scale parameters depends only on the
distribution of log(pcov), not on the experimental condition labels or
the location-model coefficients. The scale estimates are therefore
independent of the treatment effect conditional on pcov. This
independence holds because pcov is a calibration-curve property (how
well the 4PL is determined at that concentration), not a biological
outcome of the treatment. It should be verified empirically by checking
that pcov distributions are balanced across experimental conditions.

## Typical workflow

    library(curveRweights)

    # 1. Define the saturated cell-means factor externally
    dat$cell <- interaction(dat$Arm, dat$timeperiod, drop = TRUE)

    # 2. Fit per (antigen, source) group
    batch <- fit_saturated_weight_batch(
      datg       = dat,
      group_vars = c("antigen", "source"),
      cell_col   = "cell",
      pcov_col   = "pcov",
      plate_col  = "plate"
    )

    # 3. Use weights for any comparison
    dat_sub <- batch$data |>
      dplyr::filter(Arm %in% c("TT", "TdaP"), timeperiod == "post3rd")
    # dat_sub$w_saturated is ready for rma() or svyranktest()

## References

Burkner P-C (2017). brms: An R Package for Bayesian Multilevel Models
Using Stan. *Journal of Statistical Software*, 80(1), 1–28.
[doi:10.18637/jss.v080.i01](https://doi.org/10.18637/jss.v080.i01)

Burkner P-C (2018). Advanced Bayesian Multilevel Modeling with the R
Package brms. *The R Journal*, 10(1), 395–411.
[doi:10.32614/RJ-2018-017](https://doi.org/10.32614/RJ-2018-017)

Carpenter B, Gelman A, Hoffman MD, et al. (2017). Stan: A Probabilistic
Programming Language. *Journal of Statistical Software*, 76(1), 1–32.
[doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)

DerSimonian R, Laird N (1986). Meta-analysis in clinical trials.
*Controlled Clinical Trials*, 7(3), 177–188. PMID: 3802833.

Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
meta-analysis. *Statistics in Medicine*, 21(11), 1539–1558. PMID:
12111919.

Viechtbauer W (2010). Conducting Meta-Analyses in R with the metafor
Package. *Journal of Statistical Software*, 36(3), 1–48.
[doi:10.18637/jss.v036.i03](https://doi.org/10.18637/jss.v036.i03)

## See also

Useful links:

- <https://github.com/immunoplex/curveRweights>

- Report bugs at <https://github.com/immunoplex/curveRweights/issues>

## Author

**Maintainer**: Michael Scot Zens <michael.s.zens@dartmouth.edu>
