#' Apply fitted precision weights to new observations
#'
#' Uses the estimated \code{phi} and \code{beta1} from a
#' \code{precision_weights} fit to compute \eqn{\sigma = \phi\, s^{\beta_1}} and
#' \eqn{w = 1/\sigma^2} for new observations --- typically a precision grid (to
#' draw a continuous weight profile) or held-out samples --- without refitting.
#'
#' @param object A \code{precision_weights} object from
#'   [fit_precision_weights()].
#' @param newdata One of: a \code{weight_data} frame, a plain data frame with
#'   the predictor column, or a \code{calibration_result(_multiplate)} (which is
#'   run through \code{as_weight_data(source = "grid")} using the fit's design).
#'   \code{NULL} (default) re-applies to the fit's own data via the stored
#'   weights.
#' @param ... Passed to [as_weight_data()] when \code{newdata} is an ecosystem
#'   object.
#'
#' @return A data frame with the predictor, \code{sigma}, \code{w}, and
#'   \code{w_norm} columns added.
#' @seealso [fit_precision_weights()]
#' @export
predict_weights <- function(object, newdata = NULL, ...) {
  stopifnot(inherits(object, "precision_weights"))
  gamma_0 <- object$estimates$gamma_0
  gamma_1 <- object$estimates$gamma_1
  predictor <- object$scale$scale_predictor
  resp_log10 <- object$scale$response_is_log10

  if (is.null(newdata)) return(object$weights)

  if (inherits(newdata, c("calibration_result",
                          "calibration_result_multiplate"))) {
    newdata <- as_weight_data(newdata, design = object$design,
                              source = "grid", ...)
  }
  nd <- as.data.frame(newdata)

  # ensure the predictor column the estimator expects exists
  conc_col <- if ("predicted_concentration" %in% names(nd))
    "predicted_concentration" else "concentration"
  se_col   <- if ("se" %in% names(nd)) "se" else "se_concentration"

  d <- prepare_cv(nd,
                  concentration_col = conc_col,
                  se_col            = se_col,
                  pcov_col          = if ("pcov" %in% names(nd)) "pcov" else NULL,
                  predictor         = predictor,
                  response_is_log10 = resp_log10)
  d <- compute_saturated_weights(d, gamma_0 = gamma_0, gamma_1 = gamma_1)
  d$sigma  <- d$sigma_i
  d$w      <- d$w_saturated
  d$w_norm <- d$w_saturated_norm
  d
}


#' Join estimated weights back onto a data frame
#'
#' Left-joins the \code{w}, \code{w_norm}, and \code{sigma} columns from a
#' \code{precision_weights} fit onto a user data frame (e.g. the original
#' \code{$samples}) by a key.
#'
#' @param samples_df Data frame to receive the weights.
#' @param object A \code{precision_weights} object.
#' @param by Join key column name(s). Default \code{"sampleid"}.
#' @return \code{samples_df} with \code{w}, \code{w_norm}, \code{sigma} added.
#' @seealso [fit_precision_weights()]
#' @importFrom dplyr left_join all_of select any_of
#' @export
join_weights <- function(samples_df, object, by = "sampleid") {
  stopifnot(inherits(object, "precision_weights"))
  if (!all(by %in% names(object$weights)))
    stop("join_weights: key column(s) not in the weights table: ",
         paste(setdiff(by, names(object$weights)), collapse = ", "))
  if (!all(by %in% names(samples_df)))
    stop("join_weights: key column(s) not in samples_df: ",
         paste(setdiff(by, names(samples_df)), collapse = ", "))
  w <- dplyr::select(object$weights,
                     dplyr::all_of(by),
                     dplyr::any_of(c("w", "w_norm", "sigma")))
  dplyr::left_join(samples_df, w, by = by)
}
