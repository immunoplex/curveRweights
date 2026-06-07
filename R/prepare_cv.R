#' Prepare the Precision Index from Calibration Curve Output
#'
#' Computes the location response (\code{yi}), the precision index
#' (\code{cv_i}) used as the predictor in the scale submodel, and its log
#' transform (\code{log_cv}).
#'
#' @section Scale predictor (\code{predictor}):
#' The scale submodel is \eqn{\log(\sigma_i) = \gamma_0 + \gamma_1 \log(cv_i)},
#' i.e. \eqn{\sigma_i = \phi\, cv_i^{\beta_1}}. The quantity placed in
#' \code{cv_i} therefore defines what \eqn{\phi = 1} means.
#' \describe{
#'   \item{\code{"se"}}{Use \code{se_concentration} directly. In the curveR
#'     ecosystem this is the delta-method SD of the back-calculated
#'     concentration on the log10 scale, i.e. the residual SD of \code{yi}.
#'     This is the recommended predictor: it is uncapped, so it preserves the
#'     full precision gradient, and \eqn{\phi = 1} means "se_concentration is a
#'     calibrated residual SD".}
#'   \item{\code{"pcov"}}{Use the (capped) posterior CV. Legacy / foreign-data
#'     behaviour. Lossy because \code{pcov} is censored at \code{cv_x_max};
#'     equivalent to \code{"se"} only up to the constant \eqn{\ln(10)\cdot 100}
#'     (absorbed into \eqn{\phi}) and only where \code{pcov} is not capped.}
#'   \item{\code{"se_over_conc"}}{Use \code{se_col / concentration_col} (a
#'     natural-scale CV). Useful only when no calibration-curve SD/pcov is
#'     available.}
#'   \item{\code{"auto"}}{Backward-compatible default: prefer \code{pcov} when
#'     present and finite, else \code{se_over_conc} (the original behaviour of
#'     this function).}
#' }
#'
#' @section Response scale (\code{response_is_log10}):
#' \code{yi} is the location response and must be on the log10-concentration
#' scale. If \code{concentration_col} already holds log10 concentration (as in
#' a curveR \code{calibration_result} with \code{is_log_independent = TRUE}),
#' set \code{response_is_log10 = TRUE} so it is used as-is. If it holds a
#' natural-scale concentration (the foreign-data convention), leave
#' \code{FALSE} so \code{yi = log10(conc)}.
#'
#' @param df Data frame with observation-level data.
#' @param concentration_col Character: predicted concentration column.
#' @param se_col Character: SE of concentration column.
#' @param pcov_col Character: posterior CV column. \code{NULL} disables the
#'   \code{pcov} source.
#' @param predictor One of \code{"auto"}, \code{"se"}, \code{"pcov"},
#'   \code{"se_over_conc"}. See the *Scale predictor* section. Default
#'   \code{"auto"} reproduces the historical behaviour.
#' @param response_is_log10 Logical: is \code{concentration_col} already on the
#'   log10 scale? Default \code{FALSE}.
#'
#' @return The input data frame with added columns \code{yi}, \code{cv_i},
#'   \code{log_cv}, \code{cv_source}.
#'
#' @examples
#' data(example_assay)
#' dat_sub <- example_assay[example_assay$antigen == "prn" &
#'                          example_assay$feature == "IgG1", ]
#' d <- prepare_cv(dat_sub, pcov_col = "pcov")              # legacy: auto
#' d2 <- prepare_cv(dat_sub, predictor = "se")              # se as predictor
#' head(d[, c("yi", "cv_i", "log_cv", "cv_source")])
#'
#' @export
prepare_cv <- function(df,
                       concentration_col = "predicted_concentration",
                       se_col            = "se_concentration",
                       pcov_col          = "pcov",
                       predictor         = c("auto", "se", "pcov",
                                             "se_over_conc"),
                       response_is_log10 = FALSE) {

  predictor <- match.arg(predictor)

  conc <- df[[concentration_col]]
  n    <- nrow(df)

  # ---- Location response yi (log10-concentration scale) --------------------
  yi <- rep(NA_real_, n)
  if (isTRUE(response_is_log10)) {
    ok_y    <- is.finite(conc)
    yi[ok_y] <- conc[ok_y]
    ok_conc  <- is.finite(conc)          # for se_over_conc validity below
  } else {
    ok_conc       <- is.finite(conc) & conc > 0
    yi[ok_conc]   <- log10(conc[ok_conc])
  }

  se_vec   <- if (se_col   %in% names(df)) df[[se_col]]   else rep(NA_real_, n)
  pcov_vec <- if (!is.null(pcov_col) && pcov_col %in% names(df))
    df[[pcov_col]] else rep(NA_real_, n)

  cv_i      <- rep(NA_real_, n)
  cv_source <- rep(NA_character_, n)

  pick_se <- function() {
    ok <- is.finite(se_vec) & se_vec > 0
    cv_i[ok]      <<- se_vec[ok]
    cv_source[ok] <<- "se"
  }
  pick_pcov <- function() {
    ok <- is.finite(pcov_vec) & pcov_vec > 0
    cv_i[ok]      <<- pcov_vec[ok]
    cv_source[ok] <<- "pcov"
  }
  pick_se_over_conc <- function() {
    ok <- is.finite(se_vec) & se_vec > 0 & ok_conc &
      (if (response_is_log10) TRUE else conc > 0)
    denom <- if (response_is_log10) 10^conc else conc
    cv_i[ok]      <<- se_vec[ok] / denom[ok]
    cv_source[ok] <<- "se_over_conc"
  }

  switch(
    predictor,
    se           = pick_se(),
    pcov         = pick_pcov(),
    se_over_conc = pick_se_over_conc(),
    auto = {
      # historical behaviour: prefer pcov, fall back to se/conc
      use_pcov <- any(is.finite(pcov_vec) & pcov_vec > 0, na.rm = TRUE)
      if (use_pcov) {
        pick_pcov()
        # fill remaining rows from se/conc
        remaining <- is.na(cv_i) & is.finite(se_vec) & se_vec > 0 & ok_conc
        denom <- if (response_is_log10) 10^conc else conc
        cv_i[remaining]      <- se_vec[remaining] / denom[remaining]
        cv_source[remaining] <- "se_over_conc"
      } else {
        pick_se_over_conc()
      }
    }
  )

  log_cv <- rep(NA_real_, n)
  cv_ok  <- is.finite(cv_i) & cv_i > 0
  log_cv[cv_ok] <- log(cv_i[cv_ok])

  df$yi        <- yi
  df$cv_i      <- cv_i
  df$log_cv    <- log_cv
  df$cv_source <- cv_source
  df
}
