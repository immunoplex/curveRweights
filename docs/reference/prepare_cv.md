# Prepare the Precision Index from Calibration Curve Output

Computes the precision index (cv_i) and its log transform (log_cv) used
as the predictor in the scale submodel. Uses the stored posterior CV
(pcov) from the calibration curve when available; falls back to
se_concentration / predicted_concentration otherwise.

## Usage

``` r
prepare_cv(
  df,
  concentration_col = "predicted_concentration",
  se_col = "se_concentration",
  pcov_col = "pcov"
)
```

## Arguments

- df:

  Data frame with observation-level data.

- concentration_col:

  Character: name of the predicted concentration column.

- se_col:

  Character: name of the SE of concentration column.

- pcov_col:

  Character: name of the posterior CV column. Set to `NULL` to force
  computation from `se_col / concentration_col`.

## Value

The input data frame with additional columns:

- yi:

  log10(predicted_concentration). `NA` when concentration is non-finite
  or non-positive.

- cv_i:

  Precision index: pcov when available and finite, else se/conc. `NA`
  when neither source is usable.

- log_cv:

  log(cv_i). `NA` when cv_i is non-finite or non-positive.

- cv_source:

  Character indicating which source was used for each row: `"pcov"` or
  `"se_over_conc"`.

## Details

The pcov is preferred because it correctly captures the non-Gaussian
posterior near the LLOQ and ULOQ where the standard calibration curve is
flat and the delta-method se/conc approximation breaks down.

## Examples

``` r
data(example_assay)
dat_sub <- example_assay[example_assay$antigen == "prn" &
                         example_assay$feature == "IgG1", ]
d <- prepare_cv(dat_sub, pcov_col = "pcov")
head(d[, c("yi", "cv_i", "log_cv", "cv_source")])
#>             yi      cv_i     log_cv cv_source
#> 11   0.9707595 0.1126505 -2.1834647      pcov
#> 107  0.4622342 0.1212634 -2.1097902      pcov
#> 202 -1.0828663 1.3090398  0.2692939      pcov
#> 269  0.6952331 0.1134851 -2.1760839      pcov
#> 364  1.0235120 0.1133832 -2.1769824      pcov
#> 460  0.3932103 0.1253101 -2.0769637      pcov
```
