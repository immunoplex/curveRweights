#' Estimate precision weights from a curveR calibration result
#'
#' The curveR-ecosystem entry point. Estimates the power-law relationship
#' between calibration-curve precision and residual variance,
#' \eqn{\sigma_i = \phi\, se_i^{\beta_1}}, \eqn{w_i = 1/\sigma_i^2}, using a
#' joint Bayesian location-scale model with a saturated cell-means location
#' (see [fit_saturated_weight()]). The default scale predictor is the uncapped
#' \code{se_concentration} from the calibration curve, so \eqn{\phi = 1} means
#' "se_concentration is a calibrated residual SD".
#'
#' @param data Either the output of [as_weight_data()] (class
#'   \code{weight_data}), or a \code{calibration_result(_multiplate)} / plain
#'   \code{data.frame}, in which case \code{design} is required and
#'   [as_weight_data()] is called internally.
#' @param design Character vector of design-group columns (required unless
#'   \code{data} is already a \code{weight_data}).
#' @param scale_predictor \code{"se"} (default; uncapped log10-scale SD) or
#'   \code{"pcov"} (legacy, lossy because of the cap).
#' @param plate_col Plate grouping column for the location random intercept.
#'   Defaults to \code{"curve_id"} when present, else \code{NULL}.
#' @param ... MCMC controls and priors passed to [fit_saturated_weight()]
#'   (e.g. \code{iter}, \code{warmup}, \code{chains}, \code{cores},
#'   \code{adapt_delta}, \code{seed}).
#'
#' @return A \code{precision_weights} object: a list with \code{$estimates}
#'   (phi, beta1, CIs, interpretation, diagnostics), \code{$weights} (data frame
#'   keyed by \code{obs_id}/\code{sampleid}/\code{curve_id} with \code{se},
#'   \code{pcov}, \code{sigma}, \code{w}, \code{w_norm}), \code{$design},
#'   \code{$scale} (conc_scale + scale_predictor), and \code{$fit} (the brms
#'   fit).
#'
#' @examplesIf requireNamespace("curveRfreq", quietly = TRUE) && requireNamespace("brms", quietly = TRUE)
#' \donttest{
#' # mp <- curveRfreq::fit_calibration_freq_multiplate(
#' #         standards = std, samples = samples_with_design, ...)
#' # wd <- as_weight_data(mp, design = c("timeperiod", "cohort_arm"))
#' # pw <- fit_precision_weights(wd, iter = 1000, warmup = 500, chains = 2)
#' # pw$estimates$phi; pw$estimates$beta1
#' }
#'
#' @seealso [as_weight_data()], [predict_weights()], [join_weights()],
#'   [fit_saturated_weight()]
#' @export
fit_precision_weights <- function(data, design = NULL,
                                  scale_predictor = c("se", "pcov"),
                                  plate_col = NULL, ...) {

  scale_predictor <- match.arg(scale_predictor)

  if (!inherits(data, "weight_data")) {
    if (is.null(design))
      stop("fit_precision_weights: `design` is required when `data` is not ",
           "already a weight_data object (from as_weight_data()).")
    # adapter defaults only; for non-default adapter options pre-build with
    # as_weight_data(). `...` is reserved for MCMC/prior controls below.
    data <- as_weight_data(data, design = design)
  }
  design     <- attr(data, "design")
  conc_scale <- attr(data, "conc_scale")
  is_log_ind <- isTRUE(attr(data, "is_log_independent"))

  if (is.null(plate_col))
    plate_col <- if ("curve_id" %in% names(data)) "curve_id" else NULL

  # The adapter's `predicted_concentration` is log10 when is_log_independent;
  # tell the estimator not to re-log it.
  response_is_log10 <- isTRUE(is_log_ind)

  # Guard: "se" predictor is a clean log10-scale residual SD only when the
  # response is on the log10 scale. Warn otherwise.
  if (scale_predictor == "se" && !response_is_log10)
    warning("fit_precision_weights: scale_predictor = 'se' assumes ",
            "se_concentration is the SD of the log10 response, which holds ",
            "when is_log_independent = TRUE. Your input is not log10; the ",
            "interpretation of phi may not be 'calibrated SD'. Consider ",
            "scale_predictor = 'pcov'.")

  # The estimator builds a formula `yi ~ 0 + cell`; use a non-dotted name to
  # avoid any formula-parsing ambiguity with the leading dot in `.cell`.
  dfe <- as.data.frame(data)
  dfe$cell <- dfe$.cell

  sw <- fit_saturated_weight(
    df                = dfe,
    cell_col          = "cell",
    concentration_col = "predicted_concentration",
    se_col            = "se",
    pcov_col          = if ("pcov" %in% names(dfe)) "pcov" else NULL,
    plate_col         = plate_col,
    predictor         = scale_predictor,
    response_is_log10 = response_is_log10,
    ...
  )

  d <- sw$data
  weights <- data.frame(
    obs_id   = data$obs_id,
    stringsAsFactors = FALSE
  )
  if ("sampleid" %in% names(data)) weights$sampleid <- data$sampleid
  if ("curve_id" %in% names(data)) weights$curve_id <- data$curve_id
  weights$se     <- data$se
  weights$pcov   <- if ("pcov" %in% names(data)) data$pcov else NA_real_
  weights$sigma  <- d$sigma_i
  weights$w      <- d$w_saturated
  weights$w_norm <- d$w_saturated_norm

  structure(
    list(
      estimates = list(
        phi            = sw$phi,
        beta1          = sw$beta1,
        phi_CI         = sw$phi_CI,
        beta1_CI       = sw$gamma_1_CI,
        gamma_0        = sw$gamma_0,
        gamma_1        = sw$gamma_1,
        interpretation = sw$interpretation,
        diagnostics    = sw$diagnostics,
        cv_diagnostics = sw$cv_diagnostics,
        formula        = sw$formula
      ),
      weights = weights,
      design  = design,
      scale   = list(conc_scale = conc_scale,
                     scale_predictor = scale_predictor,
                     response_is_log10 = response_is_log10),
      fit     = sw$fit
    ),
    class = "precision_weights"
  )
}


#' @export
print.precision_weights <- function(x, ...) {
  cat("<precision_weights>\n")
  cat(sprintf("  phi   = %.4g  [%.3g, %.3g]\n",
              x$estimates$phi, x$estimates$phi_CI[["lo"]],
              x$estimates$phi_CI[["hi"]]))
  b_lo <- x$estimates$beta1_CI[["lo"]]; b_hi <- x$estimates$beta1_CI[["hi"]]
  cat(sprintf("  beta1 = %.4g%s\n", x$estimates$beta1,
              if (!is.na(b_lo)) sprintf("  [%.3g, %.3g]", b_lo, b_hi) else ""))
  cat("  regime:", x$estimates$interpretation, "\n")
  cat(sprintf("  predictor = %s on a %s response; %d design cell factor(s): %s\n",
              x$scale$scale_predictor, x$scale$conc_scale,
              length(x$design), paste(x$design, collapse = ", ")))
  wv <- x$weights$w[is.finite(x$weights$w)]
  cat(sprintf("  weights: %d obs, n_eff = %s\n",
              nrow(x$weights), x$estimates$diagnostics$weight$n_eff))
  invisible(x)
}

#' @export
summary.precision_weights <- function(object, ...) {
  cat("Precision-weight model (curveRweights)\n")
  cat("--------------------------------------\n")
  print(object)
  cat("\nScale submodel: log(sigma) = log(phi) + beta1 * log(",
      object$scale$scale_predictor, ")\n", sep = "")
  cat("Location:", object$estimates$formula$location, "\n")
  cat("Scale:   ", object$estimates$formula$scale, "\n")
  cat("\nConvergence: Rhat_max =", object$estimates$diagnostics$rhat_max,
      " divergences =", object$estimates$diagnostics$n_divergent, "\n")
  cat("CV diagnosis:", object$estimates$cv_diagnostics$message, "\n")
  invisible(object)
}
