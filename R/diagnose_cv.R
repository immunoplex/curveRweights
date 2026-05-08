#' Diagnose Precision Index Variation for Scale Estimation
#'
#' The scale submodel `log(sigma) ~ log_cv` requires meaningful spread in
#' log_cv across observations.  If all observations have nearly identical pcov
#' (e.g., all in the well-determined midrange), beta1 is not identifiable and
#' the model falls back to intercept-only sigma (uniform weights).
#'
#' The threshold is `sd(log_cv) >= 0.05`.  Below this, there is not enough
#' contrast in measurement precision across the concentration range to
#' distinguish differential weighting from uniform weighting.
#'
#' @param df Data frame with `cv_i` and `log_cv` columns, typically from
#'   [prepare_cv()].
#'
#' @return Named list:
#' \describe{
#'   \item{n_finite}{Number of observations with finite cv_i > 0.}
#'   \item{cv_min, cv_median, cv_max}{Summary statistics of cv_i.}
#'   \item{log_cv_sd}{Standard deviation of log_cv.}
#'   \item{log_cv_range}{Range (max - min) of log_cv.}
#'   \item{use_cv_slope}{Logical: `TRUE` if `sd(log_cv) >= 0.05`.}
#'   \item{message}{Human-readable summary of the diagnosis.}
#' }
#'
#' @examples
#' data(example_assay)
#' dat_sub <- example_assay[example_assay$antigen == "prn" &
#'                          example_assay$feature == "IgG1", ]
#' d <- prepare_cv(dat_sub, pcov_col = "pcov")
#' cv_diag <- diagnose_cv(d)
#' cat(cv_diag$message, "\n")
#' cat("sd(log_cv) =", cv_diag$log_cv_sd, "\n")
#' cat("use slope?", cv_diag$use_cv_slope, "\n")
#'
#' @importFrom stats sd median
#' @export
diagnose_cv <- function(df) {
  ok   <- is.finite(df$cv_i) & df$cv_i > 0 & is.finite(df$log_cv)
  cv   <- df$cv_i[ok]
  lcv  <- df$log_cv[ok]
  n_ok <- sum(ok)

  if (n_ok < 5) {
    return(list(
      n_finite     = n_ok,
      cv_min       = NA_real_, cv_median = NA_real_, cv_max = NA_real_,
      log_cv_sd    = NA_real_, log_cv_range = NA_real_,
      use_cv_slope = FALSE,
      message      = paste("Only", n_ok, "finite cv_i values; need >= 5")
    ))
  }

  lcv_sd  <- stats::sd(lcv)
  lcv_rng <- diff(range(lcv))
  use_sl  <- lcv_sd >= 0.05

  list(
    n_finite     = n_ok,
    cv_min       = signif(min(cv), 4),
    cv_median    = signif(stats::median(cv), 4),
    cv_max       = signif(max(cv), 4),
    log_cv_sd    = signif(lcv_sd, 4),
    log_cv_range = signif(lcv_rng, 4),
    use_cv_slope = use_sl,
    message      = if (use_sl)
      paste0("OK: sd(log_cv) = ", signif(lcv_sd, 3),
             "; beta1 identifiable from ", n_ok, " observations")
    else
      paste0("WARNING: sd(log_cv) = ", signif(lcv_sd, 3),
             " < 0.05; beta1 not identifiable. ",
             "Falling back to intercept-only sigma (uniform weights)")
  )
}
