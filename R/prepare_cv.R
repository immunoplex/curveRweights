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
#' @importFrom dplyr mutate if_else n
#' @importFrom rlang .data
#' @export
prepare_cv <- function(df,
                       concentration_col = "predicted_concentration",
                       se_col            = "se_concentration",
                       pcov_col          = "pcov") {

  use_pcov <- !is.null(pcov_col) &&
    pcov_col %in% names(df) &&
    any(is.finite(df[[pcov_col]]) & df[[pcov_col]] > 0, na.rm = TRUE)

  out <- df |>
    dplyr::mutate(
      yi = dplyr::if_else(
        is.finite(.data[[concentration_col]]) & .data[[concentration_col]] > 0,
        log10(.data[[concentration_col]]),
        NA_real_
      ),
      cv_i = dplyr::if_else(
        if (use_pcov) is.finite(.data[[pcov_col]]) & .data[[pcov_col]] > 0
        else          rep(FALSE, dplyr::n()),
        .data[[pcov_col]],
        dplyr::if_else(
          is.finite(.data[[se_col]])            & .data[[se_col]] > 0 &
          is.finite(.data[[concentration_col]]) & .data[[concentration_col]] > 0,
          .data[[se_col]] / .data[[concentration_col]],
          NA_real_
        )
      ),
      log_cv = dplyr::if_else(is.finite(.data$cv_i) & .data$cv_i > 0,
                               log(.data$cv_i), NA_real_),
      cv_source = dplyr::if_else(
        if (use_pcov) is.finite(.data[[pcov_col]]) & .data[[pcov_col]] > 0
        else          rep(FALSE, dplyr::n()),
        "pcov", "se_over_conc"
      )
    )

  out
}
