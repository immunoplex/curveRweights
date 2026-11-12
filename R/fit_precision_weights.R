#' Estimate precision weights from a curveR calibration result
#'
#' The curveR-ecosystem entry point. Two models are available.
#'
#' \code{model = "multiplicative"} (the default, unchanged from 0.2.x so that
#' existing callers such as the i-spi-compute worker keep working) fits
#' \eqn{\sigma_i = \phi\, se_i^{\beta_1}} ([fit_saturated_weight()]). It
#' attributes \emph{all} variance to the measurement, so very precisely measured
#' samples can receive extreme weights.
#'
#' \code{model = "additive"} (recommended for new analyses) fits the
#' measurement-error model of [fit_additive_weight()]:
#' \eqn{\sigma_i^2 = \tau^2 + \kappa^2 se_i^{2\beta_1}}, which separates true
#' between-subject variance (\eqn{\tau^2}) from measurement error and yields
#' weights \eqn{w_i = 1/\sigma_i^2} that are bounded above by \eqn{1/\tau^2}, a
#' per-observation \code{reliability}, and a censored likelihood for samples
#' outside the quantifiable range. The additive result has a different
#' \code{$estimates} shape (\code{tau}, \code{kappa}, \code{reliability} instead
#' of \code{phi}); callers must opt in knowingly.
#'
#' @param data Either the output of [as_weight_data()] (class
#'   \code{weight_data}), or a \code{calibration_result(_multiplate)} / plain
#'   \code{data.frame}, in which case \code{design} is required and
#'   [as_weight_data()] is called internally.
#' @param design Character vector of design-group columns (required unless
#'   \code{data} is already a \code{weight_data}).
#' @param scale_predictor Multiplicative model only: \code{"se"} (default;
#'   uncapped log10-scale SD) or \code{"pcov"} (legacy, lossy because of the
#'   cap). The additive model always uses \code{se}.
#' @param plate_col Plate grouping column for the location random intercept.
#'   Defaults to \code{"curve_id"} when present, else \code{NULL}.
#' @param patient_col Subject column for the subject random intercept (additive
#'   model). Defaults to \code{"patientid"} when present, else \code{NULL}.
#' @param model \code{"multiplicative"} (default, 0.2.x behaviour) or
#'   \code{"additive"} (recommended; see Description).
#' @param use_censoring Additive model: use the censoring columns
#'   (\code{censor_type}, \code{censor_lower}, \code{censor_upper}) when the
#'   data carry them. Default \code{TRUE}. If \code{dilution} is also present
#'   the response is placed on the final-concentration scale, which is the scale
#'   of curveRcore's bounds.
#' @param ... MCMC controls and priors passed to [fit_additive_weight()] or
#'   [fit_saturated_weight()] (e.g. \code{iter}, \code{warmup}, \code{chains},
#'   \code{cores}, \code{adapt_delta}, \code{seed}; \code{fix_beta1} for the
#'   additive model).
#'
#' @return A \code{precision_weights} object: a list with \code{$estimates},
#'   \code{$weights}, \code{$design}, \code{$scale}, and \code{$fit}.
#'
#'   For the additive model \code{$estimates} holds \code{tau}, \code{tau2},
#'   \code{kappa}, \code{beta1}, their intervals, the overall
#'   \code{reliability}, and diagnostics; \code{$weights} (keyed by
#'   \code{obs_id}/\code{sampleid}/\code{curve_id}) has \code{se}, \code{pcov},
#'   \code{v_meas}, \code{tau2}, \code{sigma} (total), \code{w},
#'   \code{reliability}, \code{reliability_lo/hi}, \code{w_norm},
#'   \code{censor_type} and \code{bound}. See [compute_additive_weights()] for
#'   what each column is for. For the multiplicative model \code{$estimates}
#'   holds \code{phi}, \code{beta1} and \code{$weights} has \code{se},
#'   \code{pcov}, \code{sigma}, \code{w}, \code{w_norm}.
#'
#' @examplesIf requireNamespace("curveRfreq", quietly = TRUE) && requireNamespace("brms", quietly = TRUE)
#' \donttest{
#' # mp <- curveRfreq::fit_calibration_freq_multiplate(
#' #         standards = std, samples = samples_with_design, ...)
#' # wd <- as_weight_data(mp, design = c("timeperiod", "cohort_arm"))
#' # pw <- fit_precision_weights(wd, iter = 1000, warmup = 500, chains = 2)
#' # pw$estimates$tau2; pw$estimates$reliability
#' }
#'
#' @seealso [as_weight_data()], [predict_weights()], [join_weights()],
#'   [fit_additive_weight()], [fit_saturated_weight()]
#' @export
fit_precision_weights <- function(data, design = NULL,
                                  scale_predictor = c("se", "pcov"),
                                  plate_col = NULL, patient_col = NULL,
                                  model = c("multiplicative", "additive"),
                                  use_censoring = TRUE, ...) {

  scale_predictor <- match.arg(scale_predictor)
  model           <- match.arg(model)

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
  if (is.null(patient_col))
    patient_col <- if ("patientid" %in% names(data)) "patientid" else NULL

  # The adapter's `predicted_concentration` is log10 when is_log_independent;
  # tell the estimator not to re-log it.
  response_is_log10 <- isTRUE(is_log_ind)

  # The estimator builds a formula `yi ~ 0 + cell`; use a non-dotted name to
  # avoid any formula-parsing ambiguity with the leading dot in `.cell`.
  dfe <- as.data.frame(data)
  dfe$cell <- dfe$.cell

  if (model == "multiplicative") {
    if (scale_predictor == "se" && !response_is_log10)
      warning("fit_precision_weights: scale_predictor = 'se' assumes ",
              "se_concentration is the SD of the log10 response, which holds ",
              "when is_log_independent = TRUE. Your input is not log10; the ",
              "interpretation of phi may not be 'calibrated SD'. Consider ",
              "scale_predictor = 'pcov'.")
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
    return(.as_precision_weights_mult(sw, data, design, conc_scale,
                                      scale_predictor, response_is_log10))
  }

  # ---- additive -------------------------------------------------------------
  if (scale_predictor != "se")
    warning("fit_precision_weights: the additive model always uses `se`; ",
            "scale_predictor is ignored.")
  cens_col <- if (use_censoring && "censor_type" %in% names(dfe) &&
                  any(dfe$censor_type %in% c("left", "right", "interval"))) "censor_type"
  dil_col  <- if (!is.null(cens_col) && "dilution" %in% names(dfe)) "dilution"
  if (!is.null(cens_col))
    message("fit_precision_weights: using censoring columns",
            if (!is.null(dil_col)) " (response on the final-concentration scale via 'dilution')")

  sw <- fit_additive_weight(
    df                = dfe,
    cell_col          = "cell",
    concentration_col = "predicted_concentration",
    se_col            = "se",
    plate_col         = plate_col,
    patient_col       = patient_col,
    censor_col        = cens_col,
    dilution_col      = dil_col,
    response_is_log10 = response_is_log10,
    ...
  )
  d <- sw$data
  weights <- data.frame(obs_id = data$obs_id, stringsAsFactors = FALSE)
  if ("sampleid" %in% names(data)) weights$sampleid <- data$sampleid
  if ("curve_id" %in% names(data)) weights$curve_id <- data$curve_id
  weights$se          <- data$se
  weights$pcov        <- if ("pcov" %in% names(data)) data$pcov else NA_real_
  weights$v_meas      <- d$v_meas
  weights$tau2        <- d$tau2
  weights$sigma       <- d$sigma_total
  weights$w           <- d$w
  weights$reliability <- d$reliability
  weights$reliability_lo <- d$reliability_lo
  weights$reliability_hi <- d$reliability_hi
  weights$w_norm      <- d$w_norm
  weights$censor_type <- if (!is.null(cens_col)) dfe$censor_type else NA_character_
  weights$bound       <- d$bound

  structure(
    list(
      estimates = list(
        model          = "additive",
        tau            = sw$tau,  tau2 = sw$tau2,
        kappa          = sw$kappa, beta1 = sw$beta1,
        tau_CI         = sw$tau_CI, kappa_CI = sw$kappa_CI,
        beta1_CI       = sw$beta1_CI,
        reliability    = sw$reliability,
        reliability_CI = sw$reliability_CI,
        interpretation = sw$diagnostics$identifiability,
        diagnostics    = sw$diagnostics,
        formula        = sw$formula
      ),
      weights = weights,
      design  = design,
      scale   = list(conc_scale = conc_scale, scale_predictor = "se",
                     response_is_log10 = response_is_log10),
      fit     = sw$fit
    ),
    class = "precision_weights"
  )
}


.as_precision_weights_mult <- function(sw, data, design, conc_scale,
                                       scale_predictor, response_is_log10) {
  d <- sw$data
  weights <- data.frame(obs_id = data$obs_id, stringsAsFactors = FALSE)
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
        model          = "multiplicative",
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

.pw_model <- function(x) x$estimates$model %||% "multiplicative"
`%||%` <- function(a, b) if (is.null(a)) b else a


#' @export
print.precision_weights <- function(x, ...) {
  e <- x$estimates
  cat("<precision_weights>  model:", .pw_model(x), "\n")
  if (identical(.pw_model(x), "additive")) {
    cat(sprintf("  tau   = %.3g  [%.3g, %.3g]   (signal SD, log10)\n",
                e$tau, e$tau_CI[["lo"]], e$tau_CI[["hi"]]))
    cat(sprintf("  kappa = %.3g  [%.3g, %.3g]   (1 = se is the measurement SD)\n",
                e$kappa, e$kappa_CI[["lo"]], e$kappa_CI[["hi"]]))
    if (!is.na(e$beta1_CI[["lo"]]))
      cat(sprintf("  beta1 = %.3g  [%.3g, %.3g]\n", e$beta1,
                  e$beta1_CI[["lo"]], e$beta1_CI[["hi"]]))
    cat(sprintf("  overall reliability = %.3f  [%.3f, %.3f]\n",
                e$reliability, e$reliability_CI[["lo"]], e$reliability_CI[["hi"]]))
    cat("  ", e$interpretation, "\n", sep = "")
  } else {
    cat(sprintf("  phi   = %.4g  [%.3g, %.3g]\n",
                e$phi, e$phi_CI[["lo"]], e$phi_CI[["hi"]]))
    b_lo <- e$beta1_CI[["lo"]]; b_hi <- e$beta1_CI[["hi"]]
    cat(sprintf("  beta1 = %.4g%s\n", e$beta1,
                if (!is.na(b_lo)) sprintf("  [%.3g, %.3g]", b_lo, b_hi) else ""))
    cat("  regime:", e$interpretation, "\n")
  }
  cat(sprintf("  predictor = %s on a %s response; %d design cell factor(s): %s\n",
              x$scale$scale_predictor, x$scale$conc_scale,
              length(x$design), paste(x$design, collapse = ", ")))
  cat(sprintf("  weights: %d obs, n_eff = %s\n",
              nrow(x$weights), e$diagnostics$weight$n_eff))
  invisible(x)
}

#' @export
summary.precision_weights <- function(object, ...) {
  cat("Precision-weight model (curveRweights)\n")
  cat("--------------------------------------\n")
  print(object)
  e <- object$estimates
  if (identical(.pw_model(object), "additive")) {
    cat("\nVariance: sigma^2 = tau^2 + kappa^2 * se^(2 * beta1)   (additive)\n")
    cat("Weight:   w = 1/sigma^2 <= 1/tau^2;  reliability = tau^2 * w\n")
    cen <- e$diagnostics$n_by_censoring
    if (!is.null(cen))
      cat("Rows used:", paste(names(cen), unlist(cen), sep = " = ", collapse = ", "), "\n")
  } else {
    cat("\nScale submodel: log(sigma) = log(phi) + beta1 * log(",
        object$scale$scale_predictor, ")   (multiplicative)\n", sep = "")
  }
  cat("Location:", e$formula$location, "\n")
  cat("Scale:   ", e$formula$scale, "\n")
  cat("\nConvergence: Rhat_max =", e$diagnostics$rhat_max,
      " divergences =", e$diagnostics$n_divergent, "\n")
  if (!is.null(e$cv_diagnostics))
    cat("CV diagnosis:", e$cv_diagnostics$message, "\n")
  invisible(object)
}
