# Correct a correlation for attenuation by measurement error

Divides an observed correlation by the square root of the product of the
two variables' mean reliabilities: \\r / \sqrt{\bar\lambda_x
\bar\lambda_y}\\.

## Usage

``` r
attenuation_correct(r, rel_x, rel_y)
```

## Arguments

- r:

  Observed correlation.

- rel_x, rel_y:

  Mean reliability of each variable (a number in (0, 1\], or a vector of
  per-observation reliabilities whose mean is used).

## Value

The disattenuated correlation. Values beyond \\\pm 1\\ indicate the
reliabilities are over-estimated or the sample is small; they are
returned as computed, with a warning.

## Examples

``` r
attenuation_correct(0.42, rel_x = 0.85, rel_y = 0.7)
#> [1] 0.5444911
```
