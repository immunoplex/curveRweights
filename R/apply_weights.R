#' Apply Previously Estimated Scale Parameters to New Data
#'
#' When you have scale estimates from a previous fit (via
#' [fit_saturated_weight()] or [fit_saturated_weight_batch()]) and want to
#' apply them to a new or subsetted data frame without re-fitting the brms
#' model.
#'
#' This is the Stage 2 function in the recommended two-stage workflow: fit
#' once on the full dataset (all arms, all timepoints), then apply weights
#' to any comparison subset.
#'
#' @param df Data frame to receive weights.
#' @param scale_table A `scale_table` from [fit_saturated_weight_batch()],
#'   or a one-row tibble/data.frame with at least `phi` and `beta1`
#'   (or `gamma_0` and `gamma_1`).
#' @param group_vars Character vector: column names to match between `df`
#'   and `scale_table` for group-specific scale parameters.
#'   `NULL` = apply a single global (phi, beta1) to all rows.
#' @param concentration_col Character: predicted concentration column name.
#' @param se_col Character: SE of concentration column name.
#' @param pcov_col Character: posterior CV column name.
#'
#' @return The input data frame with added columns: `yi`, `cv_i`, `log_cv`,
#'   `cv_source`, `sigma_i`, `w_saturated`, `w_saturated_norm`.
#'
#' @examples
#' \donttest{
#' data(example_assay)
#'
#' # Fit on one group first
#' dat_prn <- example_assay[example_assay$antigen == "prn" &
#'                          example_assay$feature == "IgG1", ]
#' dat_prn$cell <- interaction(dat_prn$group_a, dat_prn$group_b, drop = TRUE)
#' sw <- fit_saturated_weight(dat_prn, cell_col = "cell", pcov_col = "pcov",
#'                            plate_col = "plate",
#'                            iter = 1000, warmup = 500, chains = 2, cores = 2)
#'
#' # Build a scale_table manually (or use batch$scale_table)
#' st <- data.frame(antigen = "prn", phi = sw$phi, beta1 = sw$beta1)
#'
#' # Apply to new data without re-fitting
#' dat_new <- example_assay[example_assay$antigen == "prn" &
#'                          example_assay$feature == "IgG1" &
#'                          example_assay$group_b == "timepoint_3", ]
#' dat_weighted <- apply_saturated_weights(
#'   df          = dat_new,
#'   scale_table = st,
#'   group_vars  = "antigen",
#'   pcov_col    = "pcov"
#' )
#' summary(dat_weighted$w_saturated_norm)
#' }
#'
#' @seealso [fit_saturated_weight_batch()] for producing the `scale_table`,
#'   [compute_saturated_weights()] for the underlying weight computation.
#'
#' @importFrom dplyr left_join select all_of mutate if_else group_by across
#'   ungroup any_of
#' @importFrom rlang .data
#' @export
apply_saturated_weights <- function(
    df,
    scale_table,
    group_vars        = c("antigen", "source"),
    concentration_col = "predicted_concentration",
    se_col            = "se_concentration",
    pcov_col          = "pcov"
) {
  # Prepare cv columns
  df <- prepare_cv(df, concentration_col = concentration_col,
                   se_col = se_col, pcov_col = pcov_col)

  # Ensure gamma_0 / gamma_1 exist in scale_table
  if (!"gamma_0" %in% names(scale_table) && "phi" %in% names(scale_table))
    scale_table$gamma_0 <- log(scale_table$phi)
  if (!"gamma_1" %in% names(scale_table) && "beta1" %in% names(scale_table))
    scale_table$gamma_1 <- scale_table$beta1

  if (is.null(group_vars) || length(group_vars) == 0) {
    # Single global (phi, beta1)
    g0 <- scale_table$gamma_0[1]
    g1 <- scale_table$gamma_1[1]
    df <- compute_saturated_weights(df, gamma_0 = g0, gamma_1 = g1)
  } else {
    for (col in group_vars) {
      if (col %in% names(df))          df[[col]]          <- as.character(df[[col]])
      if (col %in% names(scale_table)) scale_table[[col]] <- as.character(scale_table[[col]])
    }

    df <- dplyr::left_join(
      df,
      scale_table |>
        dplyr::select(dplyr::all_of(c(group_vars, "gamma_0", "gamma_1"))),
      by = group_vars
    )

    df <- df |>
      dplyr::mutate(
        sigma_i = dplyr::if_else(
          is.finite(.data$log_cv) &
            is.finite(.data$gamma_0) &
            is.finite(.data$gamma_1),
          exp(.data$gamma_0 + .data$gamma_1 * .data$log_cv),
          NA_real_
        ),
        w_saturated = dplyr::if_else(
          is.finite(.data$sigma_i) & .data$sigma_i > 0,
          1 / .data$sigma_i^2,
          NA_real_
        )
      )

    # Normalise within each group
    df <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::mutate(
        w_saturated_norm = dplyr::if_else(
          is.finite(.data$w_saturated) & .data$w_saturated > 0,
          .data$w_saturated / mean(
            .data$w_saturated[is.finite(.data$w_saturated) & .data$w_saturated > 0],
            na.rm = TRUE
          ),
          NA_real_
        )
      ) |>
      dplyr::ungroup()

    df <- df |> dplyr::select(-dplyr::any_of(c("gamma_0", "gamma_1")))
  }

  df
}
