#' Fit Saturated Weight Models Across Multiple Groups
#'
#' Loops over groups defined by `group_vars` (e.g., antigen x source) and fits
#' one [fit_saturated_weight()] model per group.  Returns a combined data frame
#' with weights and a summary table of scale estimates.
#'
#' Each group gets its own (phi, beta1) because the pcov-to-variance
#' relationship may differ across antigens (different 4PL curve shapes) and
#' standard curve sources (different concentration ranges).
#'
#' @param datg Data frame containing all groups.
#' @param group_vars Character vector: column names defining groups.
#'   Each unique combination gets its own model.
#' @param cell_col Character: name of the saturated cell-means factor
#'   (created externally).
#' @param ... Additional arguments passed to [fit_saturated_weight()],
#'   such as `pcov_col`, `plate_col`, `iter`, `warmup`, `chains`, `cores`.
#'
#' @return Named list:
#' \describe{
#'   \item{data}{Full data frame with `w_saturated` and `w_saturated_norm`
#'     columns added (all groups combined).}
#'   \item{scale_table}{A [tibble][tibble::tibble-package] with one row per
#'     group containing: group key columns, `phi`, `beta1`, credible intervals,
#'     `interpretation`, `n_fit`, `n_eff`, `weight_ratio`.}
#'   \item{fits}{Named list of `brmsfit` objects indexed by group label.}
#'   \item{diagnostics}{Named list of per-group diagnostic lists.}
#' }
#'
#' @examples
#' \donttest{
#' data(example_assay)
#'
#' # Select IgG1 for pertussis antigens
#' dat_igg1 <- example_assay[example_assay$feature == "IgG1" &
#'                           example_assay$antigen %in% c("pt", "fha", "prn"), ]
#' dat_igg1$cell <- interaction(dat_igg1$group_a, dat_igg1$group_b, drop = TRUE)
#'
#' # Fit across antigens (reduced iterations for speed)
#' batch <- fit_saturated_weight_batch(
#'   datg       = dat_igg1,
#'   group_vars = c("antigen"),
#'   cell_col   = "cell",
#'   pcov_col   = "pcov",
#'   plate_col  = "plate",
#'   iter = 1000, warmup = 500, chains = 2, cores = 2
#' )
#'
#' # Scale summary: one row per antigen
#' batch$scale_table
#'
#' # Weighted data for one comparison
#' dat_comparison <- batch$data[batch$data$group_a %in% c("vaccine_a", "vaccine_b") &
#'                              batch$data$group_b == "timepoint_3", ]
#' nrow(dat_comparison)
#' }
#'
#' @seealso [fit_saturated_weight()] for the per-group fitting function,
#'   [apply_saturated_weights()] for applying a saved `scale_table` to new
#'   data without re-fitting.
#'
#' @importFrom dplyr distinct across all_of bind_rows select
#' @importFrom tibble tibble
#' @export
fit_saturated_weight_batch <- function(
    datg,
    group_vars = c("antigen", "source"),
    cell_col   = "cell",
    ...
) {
  for (col in group_vars)
    if (col %in% names(datg))
      datg[[col]] <- as.character(datg[[col]])

  group_keys <- datg |>
    dplyr::distinct(dplyr::across(dplyr::all_of(group_vars)))

  n_groups <- nrow(group_keys)
  message("fit_saturated_weight_batch: ", n_groups, " groups\n")

  all_data   <- vector("list", n_groups)
  all_tables <- vector("list", n_groups)
  all_fits   <- list()
  all_diag   <- list()

  for (i in seq_len(n_groups)) {
    keys      <- group_keys[i, , drop = FALSE]
    key_label <- paste(
      mapply(function(col) paste0(col, "=", keys[[col]]), group_vars),
      collapse = ", "
    )
    message("\n[", i, "/", n_groups, "] ", key_label)

    df_sub <- datg
    for (col in group_vars) {
      val    <- as.character(keys[[col]])
      df_sub <- df_sub[!is.na(df_sub[[col]]) & df_sub[[col]] == val,
                       , drop = FALSE]
    }

    sw <- tryCatch(
      fit_saturated_weight(df = df_sub, cell_col = cell_col, ...),
      error = function(e) {
        warning("Group [", key_label, "] failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(sw)) {
      df_sub$yi               <- NA_real_
      df_sub$cv_i             <- NA_real_
      df_sub$log_cv           <- NA_real_
      df_sub$cv_source        <- NA_character_
      df_sub$sigma_i          <- NA_real_
      df_sub$w_saturated      <- NA_real_
      df_sub$w_saturated_norm <- NA_real_
      all_data[[i]] <- df_sub
      all_tables[[i]] <- dplyr::bind_cols(keys, tibble::tibble(
        phi = NA_real_, beta1 = NA_real_,
        phi_lo = NA_real_, phi_hi = NA_real_,
        beta1_lo = NA_real_, beta1_hi = NA_real_,
        interpretation = "failed",
        n_fit = 0L, n_eff = NA_real_, weight_ratio = NA_real_
      ))
      next
    }

    all_data[[i]]          <- sw$data
    all_fits[[key_label]]  <- sw$fit
    all_diag[[key_label]]  <- sw$diagnostics

    all_tables[[i]] <- dplyr::bind_cols(keys, tibble::tibble(
      phi            = sw$phi,
      beta1          = sw$beta1,
      phi_lo         = sw$phi_CI[["lo"]],
      phi_hi         = sw$phi_CI[["hi"]],
      beta1_lo       = sw$gamma_1_CI[["lo"]],
      beta1_hi       = sw$gamma_1_CI[["hi"]],
      interpretation = sw$interpretation,
      n_fit          = as.integer(sw$diagnostics$n_fit),
      n_eff          = sw$diagnostics$weight$n_eff,
      weight_ratio   = sw$diagnostics$weight$weight_ratio
    ))
  }

  list(
    data        = dplyr::bind_rows(all_data),
    scale_table = dplyr::bind_rows(all_tables),
    fits        = all_fits,
    diagnostics = all_diag
  )
}
