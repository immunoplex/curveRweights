# Prepare the Precision Index from Calibration Curve Output

Computes the location response (`yi`), the precision index (`cv_i`) used
as the predictor in the scale submodel, and its log transform
(`log_cv`).

## Usage

``` r
prepare_cv(
  df,
  concentration_col = "predicted_concentration",
  se_col = "se_concentration",
  pcov_col = "pcov",
  predictor = c("auto", "se", "pcov", "se_over_conc"),
  response_is_log10 = FALSE
)
```

## Arguments

- df:

  Data frame with observation-level data.

- concentration_col:

  Character: predicted concentration column.

- se_col:

  Character: SE of concentration column.

- pcov_col:

  Character: posterior CV column. `NULL` disables the `pcov` source.

- predictor:

  One of `"auto"`, `"se"`, `"pcov"`, `"se_over_conc"`. See the *Scale
  predictor* section. Default `"auto"` reproduces the historical
  behaviour.

- response_is_log10:

  Logical: is `concentration_col` already on the log10 scale? Default
  `FALSE`.

## Value

The input data frame with added columns `yi`, `cv_i`, `log_cv`,
`cv_source`.

## Scale predictor (`predictor`)

The scale submodel is \\\log(\sigma_i) = \gamma_0 + \gamma_1
\log(cv_i)\\, i.e. \\\sigma_i = \phi\\ cv_i^{\beta_1}\\. The quantity
placed in `cv_i` therefore defines what \\\phi = 1\\ means.

- `"se"`:

  Use `se_concentration` directly. In the curveR ecosystem this is the
  delta-method SD of the back-calculated concentration on the log10
  scale, i.e. the residual SD of `yi`. This is the recommended
  predictor: it is uncapped, so it preserves the full precision
  gradient, and \\\phi = 1\\ means "se_concentration is a calibrated
  residual SD".

- `"pcov"`:

  Use the (capped) posterior CV. Legacy / foreign-data behaviour. Lossy
  because `pcov` is censored at `cv_x_max`; equivalent to `"se"` only up
  to the constant \\\ln(10)\cdot 100\\ (absorbed into \\\phi\\) and only
  where `pcov` is not capped.

- `"se_over_conc"`:

  Use `se_col / concentration_col` (a natural-scale CV). Useful only
  when no calibration-curve SD/pcov is available.

- `"auto"`:

  Backward-compatible default: prefer `pcov` when present and finite,
  else `se_over_conc` (the original behaviour of this function).

## Response scale (`response_is_log10`)

`yi` is the location response and must be on the log10-concentration
scale. If `concentration_col` already holds log10 concentration (as in a
curveR `calibration_result` with `is_log_independent = TRUE`), set
`response_is_log10 = TRUE` so it is used as-is. If it holds a
natural-scale concentration (the foreign-data convention), leave `FALSE`
so `yi = log10(conc)`.

## Examples

``` r
data(example_assay)
dat_sub <- example_assay[example_assay$antigen == "prn" &
                         example_assay$feature == "IgG1", ]
d <- prepare_cv(dat_sub, pcov_col = "pcov")              # legacy: auto
d2 <- prepare_cv(dat_sub, predictor = "se")              # se as predictor
head(d[, c("yi", "cv_i", "log_cv", "cv_source")])
#>             yi      cv_i     log_cv cv_source
#> 11   0.9707595 0.1126505 -2.1834647      pcov
#> 107  0.4622342 0.1212634 -2.1097902      pcov
#> 202 -1.0828663 1.3090398  0.2692939      pcov
#> 269  0.6952331 0.1134851 -2.1760839      pcov
#> 364  1.0235120 0.1133832 -2.1769824      pcov
#> 460  0.3932103 0.1253101 -2.0769637      pcov
```
