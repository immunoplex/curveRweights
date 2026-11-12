#' Prepare censored observations for the additive weight model
#'
#' Turns the censoring labels produced by curveRcore
#' (\code{classify_censoring_multiplate()}: \code{censor_type},
#' \code{censor_lower}, \code{censor_upper}) into the columns the additive
#' model's censored likelihood needs.
#'
#' @section In plain terms:
#' A sample outside the quantifiable range of the standard curve is not a
#' noisy number; it is a sample we only know to be \emph{below} (or \emph{above})
#' some limit. Treating its extrapolated back-calculated value as an ordinary
#' measurement is wrong, and so is giving it a weight. Instead the model uses
#' the statement "the true value is below the lower limit" (left-censored),
#' "above the upper limit" (right-censored) or "between the two" (interval).
#' This function builds, for every row, the censoring flag and the limit that
#' statement refers to.
#'
#' For a censored row the model also needs to know how precisely the limit
#' itself is determined. That is the measurement SD at the limit, not the
#' (extrapolated, usually huge) SE of the row's own back-calculated value. By
#' default it is taken from the quantified sample of the same curve whose value
#' is closest to the limit; supply \code{se_bound_col} to use a better source
#' (for example the SE read off the curve's precision grid at the limit).
#'
#' @section Scales:
#' \code{censor_lower}/\code{censor_upper} are on the natural concentration
#' scale and must be on the same scale as the response \code{yi}
#' (see the \code{dilution_col} argument of [fit_additive_weight()]). They are
#' converted to log10 here.
#'
#' @param df Data frame containing \code{yi} (log10 response) and \code{se_col}.
#' @param censor_col Column holding \code{"none"}, \code{"left"},
#'   \code{"interval"}, \code{"right"} (or \code{NA}). \code{NULL} means no
#'   censoring: every row is treated as quantified.
#' @param lower_col,upper_col Natural-scale censoring bounds.
#' @param yi_col Name of the log10 response column.
#' @param se_col Name of the measurement SE column (log10 scale).
#' @param se_bound_col Optional column with the measurement SE at the censoring
#'   limit. \code{NULL} (default) uses the nearest quantified sample of the same
#'   curve.
#' @param curve_col Column identifying the curve (used to find the nearest
#'   quantified sample). Ignored if absent.
#'
#' @return \code{df} with added columns:
#' \describe{
#'   \item{cens_flag}{brms censoring code: 0 quantified, -1 left, 1 right,
#'     2 interval. \code{NA} if the row is censored but its bound is missing
#'     (such rows cannot be used).}
#'   \item{yi}{Overwritten for censored rows by the (log10) censoring limit:
#'     the upper bound for left, the lower bound for right and interval.}
#'   \item{yi2}{Upper bound for interval rows, otherwise equal to \code{yi}.}
#'   \item{se_meas}{Measurement SE used by the model: the row's own SE for
#'     quantified rows, the SE at the limit for censored rows.}
#'   \item{bound}{log10 censoring limit (\code{NA} for quantified rows).}
#' }
#'
#' @examples
#' d <- data.frame(
#'   curve_id = "c1",
#'   yi = c(0.1, 0.8, 1.2, 1.6, NA),
#'   cv_i = c(NA, 0.12, 0.08, 0.15, NA),
#'   censor_type = c("left", "none", "none", "none", "right"),
#'   censor_lower = c(NA, NA, NA, NA, 70),
#'   censor_upper = c(5, NA, NA, NA, NA))
#' prepare_censoring(d, se_col = "cv_i")
#'
#' @seealso [fit_additive_weight()]
#' @export
prepare_censoring <- function(df,
                              censor_col   = "censor_type",
                              lower_col    = "censor_lower",
                              upper_col    = "censor_upper",
                              yi_col       = "yi",
                              se_col       = "cv_i",
                              se_bound_col = NULL,
                              curve_col    = "curve_id") {
  n    <- nrow(df)
  yi   <- df[[yi_col]]
  se   <- df[[se_col]]
  type <- if (!is.null(censor_col) && censor_col %in% names(df))
    as.character(df[[censor_col]]) else rep("none", n)
  # unlabelled rows (NA) are treated as quantified; rows with no response are
  # dropped later because yi is not finite.
  type[is.na(type)] <- "none"
  bad <- setdiff(unique(type), c("none", "left", "right", "interval"))
  if (length(bad))
    stop("prepare_censoring: unknown censor_type value(s): ",
         paste(bad, collapse = ", "),
         ". Expected none, left, right, interval.")

  lg <- function(col) {
    v <- if (!is.null(col) && col %in% names(df))
      suppressWarnings(as.numeric(df[[col]])) else rep(NA_real_, n)
    v[!is.finite(v) | v <= 0] <- NA_real_
    log10(v)
  }
  lo <- lg(lower_col); up <- lg(upper_col)

  flag  <- rep(0L, n)
  bound <- rep(NA_real_, n)
  yi2   <- yi
  i_l <- type == "left";     i_r <- type == "right"; i_i <- type == "interval"
  flag[i_l] <- -1L; flag[i_r] <- 1L; flag[i_i] <- 2L
  bound[i_l] <- up[i_l]
  bound[i_r] <- lo[i_r]
  bound[i_i] <- up[i_i]
  usable <- !(i_l | i_r | i_i) |
    (i_l & is.finite(up)) | (i_r & is.finite(lo)) |
    (i_i & is.finite(lo) & is.finite(up))
  flag[!usable] <- NA_integer_

  cens <- !is.na(flag) & flag != 0L
  yi[i_l & usable] <- up[i_l & usable]
  yi[i_r & usable] <- lo[i_r & usable]
  yi[i_i & usable] <- lo[i_i & usable]
  yi2[i_i & usable] <- up[i_i & usable]
  yi2[!(i_i & usable)] <- yi[!(i_i & usable)]

  se_meas <- se
  if (any(cens)) {
    if (!is.null(se_bound_col) && se_bound_col %in% names(df)) {
      se_meas[cens] <- df[[se_bound_col]][cens]
    } else {
      quant <- !is.na(flag) & flag == 0L & is.finite(df[[yi_col]]) &
        is.finite(se) & se > 0
      qy <- df[[yi_col]]
      curve <- if (!is.null(curve_col) && curve_col %in% names(df))
        as.character(df[[curve_col]]) else rep("all", n)
      nearest_se <- function(target, pool) {
        if (!any(pool)) return(NA_real_)
        se[pool][which.min(abs(qy[pool] - target))]
      }
      for (i in which(cens)) {
        se_meas[i] <- nearest_se(bound[i], quant & curve == curve[i])
        if (is.na(se_meas[i])) se_meas[i] <- nearest_se(bound[i], quant)
      }
    }
    n_drop <- sum(cens & !(is.finite(se_meas) & se_meas > 0))
    if (n_drop > 0)
      warning("prepare_censoring: ", n_drop, " censored row(s) have no usable ",
              "measurement SE at the limit and cannot be used.")
  }
  n_unusable <- sum(is.na(flag))
  if (n_unusable > 0)
    warning("prepare_censoring: ", n_unusable, " censored row(s) lack the ",
            "bound their censor_type requires and cannot be used.")

  df[[yi_col]]   <- yi
  df$yi2         <- yi2
  df$cens_flag   <- flag
  df$se_meas     <- se_meas
  df$bound       <- bound
  df
}
