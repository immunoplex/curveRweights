#' Interpret the Estimated Beta1 Value
#'
#' Classifies the pcov-to-variance power-law exponent into interpretable
#' regimes based on the estimated beta1 and its credible interval.
#'
#' @param beta1 Numeric: estimated beta1 (gamma_1 from the log-scale model).
#' @param beta1_lo Numeric: lower bound of 95\% credible interval.
#'   Default `NA_real_`.
#' @param beta1_hi Numeric: upper bound of 95\% credible interval.
#'   Default `NA_real_`.
#'
#' @return Character string describing the precision weighting regime:
#' \describe{
#'   \item{"not identified (uniform precision)"}{beta1 is NA}
#'   \item{"near-uniform weights"}{beta1 < 0.2}
#'   \item{"compressed precision weighting"}{0.2 <= beta1 < 0.8}
#'   \item{"calibrated (pcov ~ residual SD)"}{CI contained within 0.8 to 1.2}
#'   \item{"amplified precision weighting"}{beta1 > 1.2}
#'   \item{"moderate precision weighting"}{all other cases}
#' }
#'
#' @details
#' The theoretical prediction from the delta method applied to the 4PL
#' calibration curve is beta1 = 1 (pcov is a direct proxy for residual SD
#' on the log10-concentration scale).  Departures from 1 indicate that the
#' assay's measurement precision maps to residual variance with a different
#' power than theory predicts.
#'
#' @examples
#' interpret_beta1(0.98, 0.85, 1.12)
#' interpret_beta1(2.01, 1.47, 2.57)
#' interpret_beta1(0.56, 0.30, 0.82)
#'
#' @export
interpret_beta1 <- function(beta1, beta1_lo = NA_real_, beta1_hi = NA_real_) {
  if (is.na(beta1)) return("not identified (uniform precision)")
  dplyr::case_when(
    !is.na(beta1_lo) & beta1_lo > 0.8 &
      !is.na(beta1_hi) & beta1_hi < 1.2 ~
      "calibrated (pcov ~ residual SD)",
    beta1 < 0.2    ~ "near-uniform weights",
    beta1 < 0.8    ~ "compressed precision weighting",
    beta1 > 1.2    ~ "amplified precision weighting",
    TRUE           ~ "moderate precision weighting"
  )
}
