#' Prepare the Precision Index from Calibration Curve Output
#'
#' Computes the precision index (cv_i) and its log transform (log_cv) used as
#' the predictor in the scale submodel.  Uses the stored posterior CV (pcov)
#' from the calibration curve when available; falls back to
#' se_concentration / predicted_concentration otherwise.
#'
#' The pcov is preferred because it correctly captures the non-Gaussian
#' posterior near the LLOQ and ULOQ where the standard calibration curve is flat and the
#' delta-method se/conc approximation breaks down.
#'
#' @param df Data frame with observation-level data.
#' @param concentration_col Character: name of the predicted concentration
#'   column.
#' @param se_col Character: name of the SE of concentration column.
#' @param pcov_col Character: name of the posterior CV column.  Set to
#'   `NULL` to force computation from `se_col / concentration_col`.
#'
#' @return The input data frame with additional columns:
#' \describe{
#'   \item{yi}{log10(predicted_concentration).  `NA` when concentration is
#'     non-finite or non-positive.}
#'   \item{cv_i}{Precision index: pcov when available and finite, else
#'     se/conc.  `NA` when neither source is usable.}
#'   \item{log_cv}{log(cv_i).  `NA` when cv_i is non-finite or non-positive.}
#'   \item{cv_source}{Character indicating which source was used for each
#'     row: `"pcov"` or `"se_over_conc"`.}
#' }
#'
#' @examples
#' data(example_assay)
#' dat_sub <- example_assay[example_assay$antigen == "prn" &
#'                          example_assay$feature == "IgG1", ]
#' d <- prepare_cv(dat_sub, pcov_col = "pcov")
#' head(d[, c("yi", "cv_i", "log_cv", "cv_source")])
#'
#' @export
prepare_cv <- function(df,
                       concentration_col = "predicted_concentration",
                       se_col            = "se_concentration",
                       pcov_col          = "pcov") {

  use_pcov <- !is.null(pcov_col) &&
    pcov_col %in% names(df) &&
    any(is.finite(df[[pcov_col]]) & df[[pcov_col]] > 0, na.rm = TRUE)

  # Compute yi safely: clamp non-positive to NA before log10 to avoid NaN
  conc <- df[[concentration_col]]
  yi   <- rep(NA_real_, nrow(df))
  ok_conc <- is.finite(conc) & conc > 0
  yi[ok_conc] <- log10(conc[ok_conc])

  # Compute cv_i: two code paths to avoid .data[[NULL]] when pcov_col is NULL
  if (use_pcov) {
    pcov_vec <- df[[pcov_col]]
    se_vec   <- df[[se_col]]

    pcov_ok  <- is.finite(pcov_vec) & pcov_vec > 0
    se_ok    <- !pcov_ok & is.finite(se_vec) & se_vec > 0 & ok_conc

    cv_i      <- rep(NA_real_, nrow(df))
    cv_source <- rep(NA_character_, nrow(df))

    cv_i[pcov_ok]      <- pcov_vec[pcov_ok]
    cv_source[pcov_ok] <- "pcov"

    cv_i[se_ok]      <- se_vec[se_ok] / conc[se_ok]
    cv_source[se_ok] <- "se_over_conc"
  } else {
    se_vec <- df[[se_col]]
    se_ok  <- is.finite(se_vec) & se_vec > 0 & ok_conc

    cv_i      <- rep(NA_real_, nrow(df))
    cv_source <- rep("se_over_conc", nrow(df))

    cv_i[se_ok] <- se_vec[se_ok] / conc[se_ok]
    cv_source[!se_ok] <- NA_character_
  }

  # Compute log_cv
  log_cv <- rep(NA_real_, nrow(df))
  cv_ok  <- is.finite(cv_i) & cv_i > 0
  log_cv[cv_ok] <- log(cv_i[cv_ok])

  df$yi        <- yi
  df$cv_i      <- cv_i
  df$log_cv    <- log_cv
  df$cv_source <- cv_source

  df
}
