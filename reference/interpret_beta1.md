# Interpret the Estimated Beta1 Value

Classifies the pcov-to-variance power-law exponent into interpretable
regimes based on the estimated beta1 and its credible interval.

## Usage

``` r
interpret_beta1(beta1, beta1_lo = NA_real_, beta1_hi = NA_real_)
```

## Arguments

- beta1:

  Numeric: estimated beta1 (gamma_1 from the log-scale model).

- beta1_lo:

  Numeric: lower bound of 95\\ Default `NA_real_`.

- beta1_hi:

  Numeric: upper bound of 95\\ Default `NA_real_`.

## Value

Character string describing the precision weighting regime:

- "not identified (uniform precision)":

  beta1 is NA

- "near-uniform weights":

  beta1 \< 0.2

- "compressed precision weighting":

  0.2 \<= beta1 \< 0.8

- "calibrated (pcov ~ residual SD)":

  CI contained within 0.8 to 1.2

- "amplified precision weighting":

  beta1 \> 1.2

- "moderate precision weighting":

  all other cases

## Details

The theoretical prediction from the delta method applied to the 4PL
calibration curve is beta1 = 1 (pcov is a direct proxy for residual SD
on the log10-concentration scale). Departures from 1 indicate that the
assay's measurement precision maps to residual variance with a different
power than theory predicts.

## Examples

``` r
interpret_beta1(0.98, 0.85, 1.12)
#> [1] "calibrated (pcov ~ residual SD)"
interpret_beta1(2.01, 1.47, 2.57)
#> [1] "amplified precision weighting"
interpret_beta1(0.56, 0.30, 0.82)
#> [1] "compressed precision weighting"
```
