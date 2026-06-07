#' Build precision-weight input from a curveR calibration result
#'
#' Converts a curveRcore \code{calibration_result} or
#' \code{calibration_result_multiplate} (the object returned by
#' \code{curveRfreq::fit_calibration_freq_multiplate()} or
#' \code{curveRbayes::fit_calibration_bayes()}) into the tidy, standardized
#' data frame consumed by [fit_precision_weights()]. A plain data frame is also
#' accepted (validation/passthrough), which is the escape hatch for the
#' original foreign-data workflow.
#'
#' Extraction from the S3 object is delegated to
#' \code{curveRcore::tidy_samples()} / \code{curveRcore::tidy_grid()} so that
#' the data contract has a single owner; this adapter never reaches into object
#' internals.
#'
#' @param x A \code{calibration_result}, \code{calibration_result_multiplate},
#'   or \code{data.frame}.
#' @param design Character vector of design-group column names that define the
#'   saturated cells (e.g. \code{c("timeperiod", "cohort_arm")}). These must be
#'   present on the input. For a \code{calibration_result(_multiplate)} they are
#'   carried through from the original \code{samples} data frame passed at fit
#'   time --- if they are absent, supply them on the fitting input first, or use
#'   the \code{data.frame} method with a pre-joined frame.
#' @param source \code{"samples"} (default) extracts per-sample predictions;
#'   \code{"grid"} extracts the precision grid (used by [predict_weights()] to
#'   build a continuous weight profile).
#' @param conc_scale Location-response scale: \code{"log10"} (default, uses
#'   \code{predicted_concentration}) or \code{"natural"} (uses
#'   \code{final_concentration} when present, else \code{10^predicted_concentration}).
#' @param include_plate Logical: attach \code{curve_id} as a column and make it
#'   available as a plate grouping. Default \code{TRUE}.
#' @param plate_in_cell Logical: fold \code{curve_id} into the cell factor
#'   alongside \code{design}. Default \code{FALSE} --- plate variation is left
#'   to be absorbed by \code{phi} / the plate random effect (see the plan, §5).
#' @param drop_oor Logical: drop rows with \code{pcov_pass == FALSE}? Default
#'   \code{FALSE} (keep but warn) --- not discarding out-of-range observations
#'   is the point of the package.
#' @param ... Passed to methods.
#'
#' @return A data frame (class \code{weight_data}) with at least: \code{obs_id},
#'   \code{sampleid} (if present), \code{curve_id} (if \code{include_plate}),
#'   \code{concentration} (location response on \code{conc_scale}),
#'   \code{predicted_concentration} (log10, for the estimator), \code{se}
#'   (= \code{se_concentration}, the canonical uncapped scale predictor),
#'   \code{pcov} (reference only), \code{pcov_pass}, the \code{design} columns,
#'   and \code{.cell}. Carries attributes \code{conc_scale},
#'   \code{is_log_independent}, \code{design}, \code{plate_in_cell}.
#'
#' @seealso [fit_precision_weights()], [predict_weights()],
#'   \code{curveRcore::tidy_samples()}
#' @export
as_weight_data <- function(x, design, source = c("samples", "grid"),
                           conc_scale = c("log10", "natural"),
                           include_plate = TRUE, plate_in_cell = FALSE,
                           drop_oor = FALSE, ...)
  UseMethod("as_weight_data")

#' @rdname as_weight_data
#' @export
as_weight_data.calibration_result_multiplate <- function(
    x, design, source = c("samples", "grid"),
    conc_scale = c("log10", "natural"),
    include_plate = TRUE, plate_in_cell = FALSE, drop_oor = FALSE, ...) {

  source     <- match.arg(source)
  conc_scale <- match.arg(conc_scale)

  tidy <- switch(
    source,
    samples = curveRcore::tidy_samples(x),
    grid    = curveRcore::tidy_grid(x)
  )
  is_log_indep <- isTRUE(x$meta$is_log_independent)
  .assemble_weight_data(tidy, design = design, source = source,
                        conc_scale = conc_scale,
                        is_log_independent = is_log_indep,
                        include_plate = include_plate,
                        plate_in_cell = plate_in_cell, drop_oor = drop_oor)
}

#' @rdname as_weight_data
#' @export
as_weight_data.calibration_result <- function(
    x, design, source = c("samples", "grid"),
    conc_scale = c("log10", "natural"),
    include_plate = TRUE, plate_in_cell = FALSE, drop_oor = FALSE, ...) {

  source     <- match.arg(source)
  conc_scale <- match.arg(conc_scale)
  tidy <- switch(source,
                 samples = curveRcore::tidy_samples(x),
                 grid    = curveRcore::tidy_grid(x))
  is_log_indep <- isTRUE(x$meta$is_log_independent)
  .assemble_weight_data(tidy, design = design, source = source,
                        conc_scale = conc_scale,
                        is_log_independent = is_log_indep,
                        include_plate = include_plate,
                        plate_in_cell = plate_in_cell, drop_oor = drop_oor)
}

#' @rdname as_weight_data
#' @export
as_weight_data.data.frame <- function(
    x, design, source = c("samples", "grid"),
    conc_scale = c("log10", "natural"),
    include_plate = TRUE, plate_in_cell = FALSE, drop_oor = FALSE, ...) {

  source     <- match.arg(source)
  conc_scale <- match.arg(conc_scale)
  # A plain frame is assumed already in ecosystem convention if it has
  # se_concentration + predicted_concentration; foreign frames (natural-scale)
  # are detected via the absence of is_log_independent and handled by the
  # caller's chosen conc_scale. We do not know is_log_independent here, so we
  # assume TRUE when a log10 predicted_concentration column is present.
  is_log_indep <- "predicted_concentration" %in% names(x)
  .assemble_weight_data(x, design = design, source = source,
                        conc_scale = conc_scale,
                        is_log_independent = is_log_indep,
                        include_plate = include_plate,
                        plate_in_cell = plate_in_cell, drop_oor = drop_oor)
}


# ---- internal assembly ------------------------------------------------------

.assemble_weight_data <- function(tidy, design, source, conc_scale,
                                  is_log_independent, include_plate,
                                  plate_in_cell, drop_oor) {

  if (!is.data.frame(tidy) || nrow(tidy) == 0L)
    stop("as_weight_data: no rows extracted from the calibration result ",
         "(source = '", source, "'). Were samples provided to the fit?")

  # ---- design columns must be present --------------------------------------
  missing_design <- setdiff(design, names(tidy))
  if (length(missing_design) > 0)
    stop("as_weight_data: design column(s) not found: ",
         paste(missing_design, collapse = ", "), ".\n",
         "  Design metadata (e.g. timeperiod, cohort_arm) is carried through ",
         "from the `samples` data frame passed to the fitting call -- it is ",
         "not invented by the fitter. Add these columns to `samples` before ",
         "fitting, or use as_weight_data() on a pre-joined data.frame.")

  if (!"se_concentration" %in% names(tidy))
    stop("as_weight_data: 'se_concentration' not found in the extracted table.")
  if (!"predicted_concentration" %in% names(tidy))
    stop("as_weight_data: 'predicted_concentration' not found.")

  out <- data.frame(obs_id = seq_len(nrow(tidy)), stringsAsFactors = FALSE)
  if ("sampleid" %in% names(tidy)) out$sampleid <- tidy$sampleid
  if (include_plate && "curve_id" %in% names(tidy))
    out$curve_id <- as.character(tidy$curve_id)

  # location response on the requested scale
  pc <- tidy$predicted_concentration
  if (identical(conc_scale, "natural")) {
    if ("final_concentration" %in% names(tidy)) {
      out$concentration <- tidy$final_concentration
    } else {
      out$concentration <- if (is_log_independent) 10^pc else pc
    }
  } else {
    out$concentration <- if (is_log_independent) pc else log10(pc)
  }

  # estimator inputs
  out$predicted_concentration <- pc          # passed with response_is_log10
  out$se   <- tidy$se_concentration           # canonical uncapped predictor
  out$se_concentration <- tidy$se_concentration
  if ("pcov" %in% names(tidy))      out$pcov      <- tidy$pcov
  if ("pcov_pass" %in% names(tidy)) out$pcov_pass <- tidy$pcov_pass

  # design columns + cell factor
  for (d in design) out[[d]] <- tidy[[d]]
  cell_terms <- design
  if (plate_in_cell && "curve_id" %in% names(out))
    cell_terms <- c(design, "curve_id")
  out$.cell <- interaction(out[cell_terms], drop = TRUE, sep = "|")

  # ---- out-of-range policy --------------------------------------------------
  if ("pcov_pass" %in% names(out)) {
    n_oor <- sum(!isTRUE_vec(out$pcov_pass))
    if (n_oor > 0) {
      if (drop_oor) {
        out <- out[isTRUE_vec(out$pcov_pass), , drop = FALSE]
        message("as_weight_data: dropped ", n_oor,
                " out-of-range rows (pcov_pass == FALSE).")
      } else {
        warning("as_weight_data: ", n_oor, " rows have pcov_pass == FALSE; ",
                "kept (drop_oor = FALSE). They still carry an uncapped se and ",
                "contribute to weighting.")
      }
    }
  }

  # ---- within-cell replication check ---------------------------------------
  tab <- table(out$.cell)
  if (all(tab < 2))
    stop("as_weight_data: every cell is a singleton (max n per cell = ",
         max(tab), "). The saturated cell-means location model needs ",
         "within-cell replication to identify the residual scale. Use coarser ",
         "design grouping.")
  if (nlevels(out$.cell) < 2)
    stop("as_weight_data: only ", nlevels(out$.cell),
         " cell level. Need >= 2.")

  attr(out, "conc_scale")          <- conc_scale
  attr(out, "is_log_independent")  <- is_log_independent
  attr(out, "design")              <- design
  attr(out, "plate_in_cell")       <- plate_in_cell
  class(out) <- c("weight_data", "data.frame")
  out
}

isTRUE_vec <- function(x) !is.na(x) & x
