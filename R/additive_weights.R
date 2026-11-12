#' Compute additive-model weights, reliability and measurement variance
#'
#' Given the estimated signal variance \eqn{\tau^2} and measurement-error scale
#' \eqn{\kappa} (and optionally exponent \eqn{\beta_1}), computes for every row
#' the quantities a downstream analysis needs.
#'
#' @section In plain terms:
#' Every assay number is "true value + biology + measurement error". The
#' calibration curve tells us the measurement error of each sample
#' (\code{se}); the model estimates how much true between-subject variation
#' (\eqn{\tau^2}) sits on top of it. Because the two \emph{add}, the total
#' variance of a sample is \eqn{\tau^2 + v_i} and its inverse-variance weight
#' \eqn{w_i = 1/(\tau^2 + v_i)} can never exceed \eqn{1/\tau^2}: a perfectly
#' measured sample counts at most as much as the biology allows.
#'
#' The \strong{reliability} \eqn{\lambda_i = \tau^2/(\tau^2 + v_i)} is the share
#' of a sample's variance that is signal rather than measurement error. It runs
#' from 0 (pure noise) to 1 (perfectly measured) and reads like an intraclass
#' correlation.
#'
#' @section Output columns:
#' \describe{
#'   \item{v_meas}{Measurement variance \eqn{v_i = \kappa^2\, se_i^{2\beta_1}}.}
#'   \item{tau2}{Signal variance (the same for every row).}
#'   \item{sigma_total}{\eqn{\sqrt{\tau^2 + v_i}}.}
#'   \item{w}{Raw inverse-variance weight \eqn{1/\sigma_{total}^2}. Use for
#'     estimators where the absolute variance matters (e.g. \code{metafor}'s
#'     \code{vi} is \code{v_meas}, not \code{1/w}).}
#'   \item{reliability}{\eqn{\tau^2 w}, in (0, 1].}
#'   \item{w_norm}{\code{w} divided by its mean over the eligible rows
#'     (quantified rows), so the mean is 1. Use for software that treats
#'     weights as case counts, such as \code{brms}'s \code{weights()}.}
#'   \item{reliability_lo, reliability_hi}{Posterior interval of the
#'     reliability, present when \code{draws} is supplied.}
#'   \item{censor_type, bound}{Carried through when present. Censored rows get
#'     \code{w_norm = NA}: represent them by their censoring bound, not a weight.}
#' }
#'
#' @param df Data frame with the measurement SE column.
#' @param tau2 Signal variance \eqn{\tau^2} (point estimate).
#' @param kappa Measurement-error scale \eqn{\kappa}. Default 1 (the
#'   delta-method value: \code{se} is the true measurement SD).
#' @param beta1 Exponent on \code{se}. Default 1.
#' @param se_col Name of the measurement SE column. Default \code{"se_meas"}.
#' @param censor_col Name of the censoring-type column (optional).
#' @param draws Optional data frame of posterior draws with columns
#'   \code{tau2}, \code{kappa} and \code{beta1}; used only for the reliability
#'   interval.
#' @param level Credible level for the reliability interval.
#'
#' @return \code{df} with the columns above added.
#'
#' @examples
#' d <- data.frame(se_meas = c(0.03, 0.1, 0.3, 1))
#' compute_additive_weights(d, tau2 = 0.16, kappa = 1)
#'
#' @seealso [fit_additive_weight()], [normalize_weights()]
#' @export
compute_additive_weights <- function(df, tau2, kappa = 1, beta1 = 1,
                                     se_col = "se_meas", censor_col = NULL,
                                     draws = NULL, level = 0.95) {
  if (!se_col %in% names(df))
    stop("compute_additive_weights: se_col '", se_col, "' not found.")
  if (!is.finite(tau2) || tau2 <= 0)
    stop("compute_additive_weights: tau2 must be a positive number.")
  se <- df[[se_col]]
  ok <- is.finite(se) & se > 0

  v <- ifelse(ok, kappa^2 * se^(2 * beta1), NA_real_)
  tot <- tau2 + v
  df$v_meas      <- v
  df$tau2        <- tau2
  df$sigma_total <- sqrt(tot)
  df$w           <- 1 / tot
  df$reliability <- tau2 / tot

  elig <- is.finite(df$w)
  if (!is.null(censor_col) && censor_col %in% names(df))
    elig <- elig & (is.na(df[[censor_col]]) | df[[censor_col]] == "none")
  df$w_norm <- NA_real_
  if (any(elig)) df$w_norm[elig] <- df$w[elig] / mean(df$w[elig])

  if (!is.null(draws) && nrow(draws) > 0L && any(ok)) {
    a  <- (1 - level) / 2
    d  <- seq_len(nrow(draws))
    m  <- vapply(d, function(k)
      draws$tau2[k] / (draws$tau2[k] + draws$kappa[k]^2 * se[ok]^(2 * draws$beta1[k])),
      numeric(sum(ok)))
    m  <- matrix(m, nrow = sum(ok))
    q  <- apply(m, 1, stats::quantile, probs = c(a, 1 - a), names = FALSE)
    df$reliability_lo <- NA_real_; df$reliability_hi <- NA_real_
    df$reliability_lo[ok] <- q[1, ]
    df$reliability_hi[ok] <- q[2, ]
  }
  df
}


#' Normalize weights to mean 1 within an analysis subset
#'
#' Rescales a weight column so its mean is 1 over the rows you will actually
#' analyse (optionally within groups).
#'
#' @section Why this matters:
#' Software that treats weights as case counts (\code{brms}'s
#' \code{| weights(w)}, \code{lm(weights = )}) assumes the weights sum to the
#' number of observations. Normalizing over the full fitting frame and then
#' filtering to a subset (one arm, one time point, quantified rows only) breaks
#' that: the subset's mean is no longer 1, and uncertainty is over- or
#' under-stated. Subset first, then call this function.
#'
#' Rescaling never changes estimates from weighted least squares, but it does
#' change standard errors in software that reads weights as counts. Do not use
#' it to make weights fall between 0 and 1; for a bounded quantity use
#' \code{reliability}.
#'
#' @param df Data frame (already subset to the analysis rows).
#' @param w_col Weight column to rescale. Default \code{"w"}.
#' @param by Optional character vector of grouping columns; each group is
#'   normalized separately.
#' @param out_col Name of the new column. Default \code{"w_norm"}.
#'
#' @return \code{df} with \code{out_col} (mean 1 over finite positive weights).
#'
#' @examples
#' d <- data.frame(g = c("a", "a", "b", "b"), w = c(1, 3, 10, 30))
#' normalize_weights(d, by = "g")
#' @export
normalize_weights <- function(df, w_col = "w", by = NULL, out_col = "w_norm") {
  if (!w_col %in% names(df))
    stop("normalize_weights: column '", w_col, "' not found.")
  w <- df[[w_col]]
  g <- if (is.null(by)) rep(1L, nrow(df)) else
    as.integer(factor(do.call(paste, c(df[by], sep = "\r"))))
  out <- rep(NA_real_, nrow(df))
  for (k in unique(g)) {
    i  <- g == k & is.finite(w) & w > 0
    if (any(i)) out[i] <- w[i] / mean(w[i])
  }
  df[[out_col]] <- out
  df
}


#' Measurement variance of each observation
#'
#' Returns the measurement-error variance \eqn{v_i} (not the total variance).
#' This is the right input for \code{metafor::rma.mv(V = )} and
#' \code{rma(vi = )}: those functions add the between-subject variance
#' \eqn{\tau^2} themselves and estimate it by REML, so passing the total
#' variance would count the signal twice.
#'
#' @param x A \code{precision_weights} object (additive model) or a data frame
#'   with a \code{v_meas} column.
#' @param quantified_only Logical: return \code{NA} for censored rows? Default
#'   \code{TRUE}.
#' @return Numeric vector, one value per row of the weight table.
#' @examples
#' d <- data.frame(v_meas = c(0.01, 0.2), censor_type = c("none", "left"))
#' meas_variance(d)
#' @export
meas_variance <- function(x, quantified_only = TRUE) {
  d <- if (inherits(x, "precision_weights")) x$weights else x
  if (!"v_meas" %in% names(d))
    stop("meas_variance: no 'v_meas' column; fit the additive model ",
         "(fit_precision_weights(model = \"additive\")).")
  v <- d$v_meas
  if (quantified_only && "censor_type" %in% names(d))
    v[!(is.na(d$censor_type) | d$censor_type == "none")] <- NA_real_
  v
}


#' Regression calibration: shrink measured values toward their cell mean
#'
#' Replaces each measurement \eqn{y_i} by its best prediction of the true value,
#' \eqn{\mu_{cell} + \lambda_i (y_i - \mu_{cell})}. Use this when the measured
#' concentration is a \emph{predictor} (correlate-of-protection regressions,
#' antigen-versus-antigen correlations).
#'
#' @section In plain terms:
#' Weights do not repair measurement error in an explanatory variable. Noise in
#' X makes a regression slope too small (attenuation, or "regression dilution").
#' A reliability \eqn{\lambda} of 0.8 means roughly 80\% of the variation you
#' see is real and the slope shrinks by that factor. Replacing X by the
#' shrunken value above removes that bias to first order. Also consider
#' \code{brms::me()} for a fully model-based approach.
#'
#' @param y Numeric vector of measurements (log10 scale).
#' @param cell Grouping vector (the design cell, e.g. arm by time point); the
#'   shrinkage target is the reliability-weighted mean of each cell.
#' @param reliability Numeric vector of per-observation reliabilities in (0, 1].
#' @return Numeric vector of calibrated values, same length as \code{y}.
#' @examples
#' set.seed(1)
#' cell <- rep(c("a", "b"), each = 50)
#' rel  <- runif(100, 0.4, 0.95)
#' y    <- rnorm(100, ifelse(cell == "a", 1, 2))
#' head(reliability_correct(y, cell, rel))
#' @export
reliability_correct <- function(y, cell, reliability) {
  stopifnot(length(y) == length(cell), length(y) == length(reliability))
  out <- rep(NA_real_, length(y))
  for (k in unique(cell)) {
    i  <- cell == k & is.finite(y) & is.finite(reliability)
    if (!any(i)) next
    mu <- stats::weighted.mean(y[i], reliability[i])
    out[i] <- mu + reliability[i] * (y[i] - mu)
  }
  out
}


#' Correct a correlation for attenuation by measurement error
#'
#' Divides an observed correlation by the square root of the product of the two
#' variables' mean reliabilities: \eqn{r / \sqrt{\bar\lambda_x \bar\lambda_y}}.
#'
#' @param r Observed correlation.
#' @param rel_x,rel_y Mean reliability of each variable (a number in (0, 1], or
#'   a vector of per-observation reliabilities whose mean is used).
#' @return The disattenuated correlation. Values beyond \eqn{\pm 1} indicate the
#'   reliabilities are over-estimated or the sample is small; they are returned
#'   as computed, with a warning.
#' @examples
#' attenuation_correct(0.42, rel_x = 0.85, rel_y = 0.7)
#' @export
attenuation_correct <- function(r, rel_x, rel_y) {
  lx <- mean(rel_x, na.rm = TRUE); ly <- mean(rel_y, na.rm = TRUE)
  out <- r / sqrt(lx * ly)
  if (any(abs(out) > 1, na.rm = TRUE))
    warning("attenuation_correct: corrected correlation exceeds 1 in absolute ",
            "value; the reliabilities are probably over-estimated.")
  out
}
