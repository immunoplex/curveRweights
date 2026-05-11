# Example Luminex Multiplex Immunoassay Data

Anonymised observation-level data from a Luminex multiplex immunoassay
measuring antibody responses across two vaccine groups, four timepoints,
eleven antigens, and ten immunological features. Each observation
includes the predicted concentration from a semi-supervised, best
fitting calibration curve, its standard error, and the posterior
coefficient of variation (pcov) that quantifies calibration curve
uncertainty. All calculations were computed in I-SPI using the bayesian
regression approach (https://immunoplex.org/i-spi-docs) .

## Usage

``` r
example_assay
```

## Format

A data frame with 48,224 rows and 11 columns:

- plate:

  Character. Plate identifier (15 plates: plate_1 to plate_15).

- nominal_sample_dilution:

  Character. Anonymised dilution factor (8 levels: dil_1 to dil_8).

- patientid:

  Character. Anonymised subject identifier (150 subjects: S001 to S150).

- group_a:

  Character. Vaccine group (2 levels: vaccine_a, vaccine_b). This is the
  treatment arm used as the independent variable in arm-effect
  comparisons.

- group_b:

  Character. Timepoint (4 levels: timepoint_1 to timepoint_4, in
  chronological order).

- antigen:

  Character. Target antigen (11 levels: act, dt, fha, fim, ipv1, ipv2,
  ipv3, pentamer, prn, pt, tt).

- feature:

  Character. Immunological readout (10 levels: ADCD, ADCP, ADNP, FcgR2a,
  FcgR3b, IgG1, IgG2, IgG3, IgG4, Total_IgG).

- mfi:

  Numeric. Median fluorescence intensity (raw assay readout).

- predicted_concentration:

  Numeric. Concentration predicted from the 4PL calibration curve
  (arbitrary units).

- se_concentration:

  Numeric. Standard error of the predicted concentration from the
  calibration curve posterior.

- pcov:

  Numeric. Posterior coefficient of variation from the calibration
  curve: se_concentration / predicted_concentration for midrange
  samples, but correctly capturing the non-Gaussian posterior near LLOQ
  and ULOQ. This is the precision index used by `curveRweights` to
  estimate observation-level weights.

## Source

Anonymised from a maternal pertussis vaccine immunogenicity study.

## Details

The data were anonymised by replacing patient IDs with random codes
(S001 to S150, shuffled), sample dilutions with generic labels (dil_1 to
dil_8), and timepoint names with sequential labels (timepoint_1 to
timepoint_4). All measurement values (mfi, predicted_concentration,
se_concentration, pcov) are unchanged.

Typical usage selects one (antigen, feature) combination for analysis.
For example, IgG1 responses to pertussis antigens (pt, fha, prn) are a
natural starting point:

    dat_sub <- example_assay |>
      dplyr::filter(antigen == "prn", feature == "IgG1")

## Examples

``` r
data(example_assay)
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
table(example_assay$antigen, example_assay$feature)
#>           
#>            ADCD ADCP ADNP FcgR2a FcgR3b IgG1 IgG2 IgG3 IgG4 Total_IgG
#>   act       510    0    0    512    512  512  512  512  512       446
#>   dt        510  502  502    512    512  512  512  512  512       446
#>   fha       510    0    0    511    512  512  512  512  512       446
#>   fim       510    0    0    512    512  512  512  512  512       446
#>   ipv1      510    0    0    511    512  512  417  512  512       446
#>   ipv2      510    0    0    512    512  512  512  512  512       446
#>   ipv3      510    0    0    512    512  512  512  512  512       446
#>   pentamer  510    0    0    512    512  512  512  512  512       446
#>   prn       510  502  502    511    512  512  512  512  512       446
#>   pt        510  502  502    511    512  512  512  512  512       446
#>   tt        510  502  502    511    512  512  512  512  512       446
```
