#' Compute Precision Weights from Estimated Scale Parameters
#'
#' Given the estimated log-scale intercept (gamma_0) and slope (gamma_1),
#' computes the per-observation residual SD (sigma_i), raw precision weight
#' (w_saturated), and mean-normalised weight (w_saturated_norm).
#'
#' The transformation is:
#' \deqn{\sigma_i = \exp(\gamma_0 + \gamma_1 \cdot \log(\mathrm{cv}_i))}{sigma_i = exp(gamma_0 + gamma_1 * log(cv_i))}
#' \deqn{w_i = 1 / \sigma_i^2}{w_i = 1 / sigma_i^2}
#'
#' @param df Data frame with a `log_cv` column (from [prepare_cv()]).
#' @param gamma_0 Numeric: estimated intercept of the log-sigma model
#'   (log(phi)).
#' @param gamma_1 Numeric: estimated slope of the log-sigma model (beta1).
#'   Set to 0 for intercept-only (uniform weights).  Default `0`.
#'
#' @return The input data frame with additional columns:
#' \describe{
#'   \item{sigma_i}{Estimated residual SD for each observation.}
#'   \item{w_saturated}{Raw precision weight: `1 / sigma_i^2`.}
#'   \item{w_saturated_norm}{Mean-normalised weight: `w_saturated / mean(w_saturated)`.
#'     The shape is identical to w_saturated; normalisation to mean = 1 is
#'     applied only for cross-group comparability.}
#' }
#'
#' @examples
#' data(example_assay)
#' dat_sub <- example_assay[example_assay$antigen == "prn" &
#'                          example_assay$feature == "IgG1", ]
#' d <- prepare_cv(dat_sub, pcov_col = "pcov")
#'
#' # Apply hypothetical scale estimates
#' d <- compute_saturated_weights(d, gamma_0 = 0.5, gamma_1 = 1.2)
#' summary(d$w_saturated_norm)
#' weight_diagnostics(d$w_saturated)
#'
#' @importFrom dplyr mutate if_else
#' @importFrom rlang .data
#' @export
compute_saturated_weights <- function(df, gamma_0, gamma_1 = 0) {

  out <- df |>
    dplyr::mutate(
      sigma_i = dplyr::if_else(
        is.finite(.data$log_cv),
        exp(gamma_0 + gamma_1 * .data$log_cv),
        NA_real_
      ),
      w_saturated = dplyr::if_else(
        is.finite(.data$sigma_i) & .data$sigma_i > 0,
        1 / .data$sigma_i^2,
        NA_real_
      )
    )

  w_valid <- out$w_saturated[is.finite(out$w_saturated) & out$w_saturated > 0]
  w_mean  <- if (length(w_valid) > 0) mean(w_valid) else NA_real_

  out <- out |>
    dplyr::mutate(
      w_saturated_norm = dplyr::if_else(
        is.finite(.data$w_saturated) & .data$w_saturated > 0 & is.finite(w_mean),
        .data$w_saturated / w_mean,
        NA_real_
      )
    )

  out
}


#' Compute Summary Diagnostics for Precision Weights
#'
#' Summarises the distribution and effective information content of a set of
#' precision weights.
#'
#' @param w Numeric vector of weights (e.g., `w_saturated` or
#'   `w_saturated_norm`).
#'
#' @return Named list:
#' \describe{
#'   \item{n_obs}{Total observations (including NA).}
#'   \item{n_valid}{Observations with finite, positive weights.}
#'   \item{n_eff}{Effective sample size:
#'     \eqn{[\sum w_i]^2 / \sum w_i^2}{[sum(w)]^2 / sum(w^2)}.
#'     Equals n_valid when all weights are equal; decreases as weights
#'     become more heterogeneous.}
#'   \item{eff_ratio}{n_eff / n_valid.  Ranges from 0 to 1; 1 = uniform.}
#'   \item{weight_ratio}{max(w) / min(w) among valid weights.}
#'   \item{gini}{Gini coefficient of weights.  0 = perfectly uniform,
#'     approaching 1 = highly concentrated.}
#' }
#'
#' @examples
#' w <- c(1.5, 1.2, 0.8, 0.3, 0.1)
#' weight_diagnostics(w)
#'
#' @export
weight_diagnostics <- function(w) {
  n_obs   <- length(w)
  ok      <- is.finite(w) & w > 0
  w_valid <- w[ok]
  n_valid <- length(w_valid)

  if (n_valid < 2)
    return(list(n_obs = n_obs, n_valid = n_valid,
                n_eff = NA_real_, eff_ratio = NA_real_,
                weight_ratio = NA_real_, gini = NA_real_))

  n_eff <- sum(w_valid)^2 / sum(w_valid^2)
  wr    <- max(w_valid) / min(w_valid)

  ws   <- sort(w_valid)
  n    <- length(ws)
  gini <- sum((2 * seq_len(n) - n - 1) * ws) / (n * sum(ws))

  list(
    n_obs        = n_obs,
    n_valid      = n_valid,
    n_eff        = round(n_eff, 1),
    eff_ratio    = round(n_eff / n_valid, 3),
    weight_ratio = signif(wr, 4),
    gini         = round(gini, 3)
  )
}
