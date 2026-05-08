#' Fit a Bayesian Location-Scale Model with Saturated Location and Shared Scale
#'
#' Fits a brms Gaussian location-scale model where the location uses saturated
#' cell means (one coefficient per level of \code{cell_col}) to absorb ALL
#' systematic variation in the response, and the scale estimates the power-law
#' relationship between the calibration curve's posterior CV (pcov) and residual
#' variance, shared across all cells.
#'
#' The location model uses \code{yi ~ 0 + cell [+ (1|plate)]} (saturated cell
#' means).  The scale model uses \code{log(sigma) = gamma_0 + gamma_1 * log(cv)}.
#'
#' Default prior on gamma_0 is Normal(0, 1), centering phi = exp(gamma_0) at 1.
#' Default prior on gamma_1 is Normal(1, 0.5), centering beta1 at the
#' delta-method prediction of 1, with a 95 percent interval from 0 to 2.
#'
#' @param df Data frame with observation-level data.
#' @param cell_col Character: name of the saturated cell-means factor column.
#'   Create this externally, e.g.:
#'   \code{df$cell <- interaction(df$Arm, df$Timeperiod, drop = TRUE)}.
#' @param concentration_col Character: predicted concentration column name.
#' @param se_col Character: SE of concentration column name.
#' @param pcov_col Character: posterior CV column name.
#'   \code{NULL} = compute from se/conc.
#' @param plate_col Character: plate column name for random intercept.
#'   \code{NULL} = no plate random effect.
#' @param prior_gamma0 A \code{brms::prior} object for gamma_0 (log phi).
#' @param prior_gamma1 A \code{brms::prior} object for gamma_1 (beta1).
#' @param prior_location A \code{brms::prior} object for location (cell mean)
#'   coefficients.
#' @param prior_plate_sd A \code{brms::prior} object for plate RE standard
#'   deviation.
#' @param iter Integer: MCMC iterations per chain.
#' @param warmup Integer: warmup iterations per chain.
#' @param chains Integer: number of MCMC chains.
#' @param cores Integer: parallel cores for chain computation.
#' @param adapt_delta Numeric: target HMC acceptance probability (0-1).
#' @param seed Integer: random seed for reproducibility.
#'
#' @return Named list with elements: \code{fit} (brmsfit object), \code{phi}
#'   (baseline scaling factor), \code{beta1} (precision exponent),
#'   \code{gamma_0} and \code{gamma_1} (log-scale estimates),
#'   \code{gamma_0_CI} and \code{gamma_1_CI} (95 percent credible intervals),
#'   \code{phi_CI} (back-transformed CI for phi),
#'   \code{interpretation} (precision regime label),
#'   \code{effective_se_power} (2 * beta1),
#'   \code{data} (input data with weights added: \code{yi}, \code{cv_i},
#'   \code{log_cv}, \code{sigma_i}, \code{w_saturated}, \code{w_saturated_norm}),
#'   \code{diagnostics} (convergence and weight diagnostics),
#'   \code{cv_diagnostics} (output of \code{diagnose_cv}),
#'   \code{formula} (location and scale formula strings),
#'   \code{priors_used} (priors passed to brms).
#'
#' @examples
#' \donttest{
#' data(example_assay)
#'
#' # Select one antigen/feature group and create the cell factor
#' dat_prn <- example_assay[example_assay$antigen == "prn" &
#'                          example_assay$feature == "IgG1", ]
#' dat_prn$cell <- interaction(dat_prn$group_a, dat_prn$group_b, drop = TRUE)
#'
#' # Fit (reduced iterations for speed; use 4000/1000 for real analysis)
#' sw <- fit_saturated_weight(
#'   df        = dat_prn,
#'   cell_col  = "cell",
#'   pcov_col  = "pcov",
#'   plate_col = "plate",
#'   iter = 1000, warmup = 500, chains = 2, cores = 2
#' )
#'
#' cat("phi =", sw$phi, " beta1 =", sw$beta1, "\n")
#' cat(sw$interpretation, "\n")
#'
#' # Weights are on sw$data
#' summary(sw$data$w_saturated_norm)
#' weight_diagnostics(sw$data$w_saturated)
#' }
#'
#' @seealso \code{\link{fit_saturated_weight_batch}},
#'   \code{\link{apply_saturated_weights}},
#'   \code{\link{prepare_cv}},
#'   \code{\link{diagnose_cv}},
#'   \code{\link{weight_diagnostics}}
#'
#' @importFrom brms brm bf fixef rhat nuts_params VarCorr
#' @importFrom dplyr n_distinct filter
#' @importFrom stats as.formula gaussian
#' @importFrom rlang .data
#' @export
fit_saturated_weight <- function(
    df,
    cell_col          = "cell",
    concentration_col = "predicted_concentration",
    se_col            = "se_concentration",
    pcov_col          = "pcov",
    plate_col         = "plate",
    prior_gamma0      = brms::set_prior("normal(0, 1)",   class = "Intercept",
                                         dpar = "sigma"),
    prior_gamma1      = brms::set_prior("normal(1, 0.5)", class = "b",
                                         dpar = "sigma"),
    prior_location    = brms::set_prior("normal(0, 2)",   class = "b"),
    prior_plate_sd    = brms::set_prior("normal(0, 0.5)", class = "sd"),
    iter              = 4000,
    warmup            = 1000,
    chains            = 4,
    cores             = 4,
    adapt_delta       = 0.95,
    seed              = 42
) {

  # ---- Validate inputs -----------------------------------------------------
  if (!cell_col %in% names(df))
    stop("fit_saturated_weight: cell_col '", cell_col,
         "' not found in data frame. ",
         "Create it externally, e.g.: ",
         "df$cell <- interaction(df$Arm, df$Timeperiod, drop = TRUE)")

  if (!concentration_col %in% names(df))
    stop("fit_saturated_weight: concentration_col '", concentration_col,
         "' not found.")

  if (!se_col %in% names(df))
    stop("fit_saturated_weight: se_col '", se_col, "' not found.")

  # ---- Prepare data --------------------------------------------------------
  df[[cell_col]] <- factor(df[[cell_col]])
  n_cells        <- nlevels(df[[cell_col]])

  if (n_cells < 2)
    stop("fit_saturated_weight: cell_col '", cell_col,
         "' has fewer than 2 levels (", n_cells, "). ",
         "Need at least 2 cells for meaningful residuals.")

  d <- prepare_cv(df,
                  concentration_col = concentration_col,
                  se_col            = se_col,
                  pcov_col          = pcov_col)

  d_fit <- d[is.finite(d$yi) & is.finite(d$cv_i) & d$cv_i > 0 &
               is.finite(d$log_cv), , drop = FALSE]

  n_input   <- nrow(df)
  n_fit     <- nrow(d_fit)
  n_removed <- n_input - n_fit

  message("fit_saturated_weight: ",
          n_fit, " of ", n_input, " observations usable (",
          n_removed, " removed); ",
          n_cells, " cell levels")

  if (n_fit < 20)
    stop("fit_saturated_weight: only ", n_fit,
         " usable observations. Need >= 20.")

  # ---- Diagnose CV variation -----------------------------------------------
  cv_diag      <- diagnose_cv(d_fit)
  use_cv_slope <- cv_diag$use_cv_slope
  message("  cv: ", cv_diag$message)

  # ---- Build formula -------------------------------------------------------
  loc_rhs <- paste0("0 + ", cell_col)

  has_plate <- !is.null(plate_col) &&
    plate_col %in% names(d_fit) &&
    dplyr::n_distinct(d_fit[[plate_col]]) > 1

  if (has_plate)
    loc_rhs <- paste0(loc_rhs, " + (1 | ", plate_col, ")")

  sig_rhs <- if (use_cv_slope) "log_cv" else "1"

  bf_formula <- brms::bf(
    stats::as.formula(paste("yi ~", loc_rhs)),
    stats::as.formula(paste("sigma ~", sig_rhs))
  )

  message("  location: yi ~ ", loc_rhs)
  message("  scale:    sigma ~ ", sig_rhs)

  # ---- Priors --------------------------------------------------------------
  model_priors <- c(
    prior_location,
    prior_gamma0,
    if (use_cv_slope) prior_gamma1,
    if (has_plate)    prior_plate_sd
  )

  # ---- Fit -----------------------------------------------------------------
  message("  fitting brms model (", iter, " iter, ", chains, " chains)...")

  fit <- tryCatch(
    brms::brm(
      formula = bf_formula,
      data    = d_fit,
      family  = stats::gaussian(),
      prior   = model_priors,
      iter    = iter,
      warmup  = warmup,
      chains  = chains,
      cores   = cores,
      seed    = seed,
      control = list(adapt_delta = adapt_delta),
      silent  = 2,
      refresh = 0
    ),
    error = function(e) {
      stop("fit_saturated_weight: brms::brm() failed: ",
           conditionMessage(e))
    }
  )

  # ---- Convergence diagnostics ---------------------------------------------
  rhat_max <- max(brms::rhat(fit), na.rm = TRUE)
  n_div    <- tryCatch(
    sum(brms::nuts_params(fit, pars = "divergent__")$Value),
    error = function(e) NA_integer_
  )

  if (!is.na(rhat_max) && rhat_max > 1.05)
    warning("fit_saturated_weight: Rhat = ", signif(rhat_max, 4),
            " > 1.05 --- chains may not have converged")
  if (!is.na(n_div) && n_div > 0)
    warning("fit_saturated_weight: ", n_div, " divergent transitions")

  # ---- Extract scale parameters (robust to brms naming) --------------------
  sig_fix    <- brms::fixef(fit, dpar = "sigma")
  sig_rows   <- rownames(sig_fix)
  sig_cols   <- colnames(sig_fix)

  # Debug: log what brms returned so subscript errors are diagnosable
  message("  sigma fixef rows: ", paste(sig_rows, collapse = ", "))
  message("  sigma fixef cols: ", paste(sig_cols, collapse = ", "))

  # Find the intercept row: may be "Intercept", "sigma_Intercept", or similar
  int_row <- sig_rows[grepl("Intercept", sig_rows, fixed = TRUE)]
  if (length(int_row) == 0) int_row <- sig_rows[1]  # fallback: first row
  int_row <- int_row[1]

  # Find the estimate column: "Estimate" or possibly different capitalisation
  est_col <- sig_cols[grepl("Estimate", sig_cols, ignore.case = TRUE)]
  if (length(est_col) == 0) est_col <- sig_cols[1]
  est_col <- est_col[1]

  # CI columns
  ci_lo_col <- sig_cols[grepl("2.5", sig_cols, fixed = TRUE)]
  ci_hi_col <- sig_cols[grepl("97.5", sig_cols, fixed = TRUE)]
  ci_lo_col <- if (length(ci_lo_col) > 0) ci_lo_col[1] else NA_character_
  ci_hi_col <- if (length(ci_hi_col) > 0) ci_hi_col[1] else NA_character_

  gamma_0    <- sig_fix[int_row, est_col]
  gamma_0_lo <- if (!is.na(ci_lo_col)) sig_fix[int_row, ci_lo_col] else NA_real_
  gamma_0_hi <- if (!is.na(ci_hi_col)) sig_fix[int_row, ci_hi_col] else NA_real_

  # Find log_cv slope row
  log_cv_row <- sig_rows[grepl("log_cv", sig_rows, fixed = TRUE)]

  if (use_cv_slope && length(log_cv_row) >= 1) {
    log_cv_row <- log_cv_row[1]
    gamma_1    <- sig_fix[log_cv_row, est_col]
    gamma_1_lo <- if (!is.na(ci_lo_col)) sig_fix[log_cv_row, ci_lo_col] else NA_real_
    gamma_1_hi <- if (!is.na(ci_hi_col)) sig_fix[log_cv_row, ci_hi_col] else NA_real_
  } else {
    gamma_1    <- 0
    gamma_1_lo <- NA_real_
    gamma_1_hi <- NA_real_
    if (use_cv_slope)
      warning("fit_saturated_weight: log_cv not found in sigma fixef. ",
              "Available rows: ", paste(sig_rows, collapse = ", "))
  }

  phi      <- exp(gamma_0)
  phi_lo   <- exp(gamma_0_lo)
  phi_hi   <- exp(gamma_0_hi)
  beta1    <- gamma_1
  interp   <- interpret_beta1(beta1, gamma_1_lo, gamma_1_hi)

  message("  phi = ", signif(phi, 4),
          "  [", signif(phi_lo, 3), ", ", signif(phi_hi, 3), "]")
  message("  beta1 = ", signif(beta1, 4),
          if (!is.na(gamma_1_lo))
            paste0("  [", signif(gamma_1_lo, 3), ", ", signif(gamma_1_hi, 3), "]")
          else "")
  message("  interpretation: ", interp)

  # ---- Compute weights on the FULL input data frame ------------------------
  d_out <- compute_saturated_weights(d, gamma_0 = gamma_0, gamma_1 = gamma_1)

  w_diag <- weight_diagnostics(
    d_out$w_saturated[is.finite(d_out$yi)]
  )

  message("  n_eff = ", w_diag$n_eff, " of ", w_diag$n_valid,
          " (ratio = ", w_diag$eff_ratio, ")")
  message("  weight_ratio = ", w_diag$weight_ratio,
          "  gini = ", w_diag$gini)

  # ---- Plate RE variance (if fitted) ---------------------------------------
  plate_sd <- NA_real_
  if (has_plate) {
    vc <- brms::VarCorr(fit)
    if (!is.null(vc[[plate_col]]))
      plate_sd <- vc[[plate_col]]$sd[1, "Estimate"]
  }

  # ---- Return ---------------------------------------------------------------
  list(
    fit                = fit,
    phi                = phi,
    beta1              = beta1,
    gamma_0            = gamma_0,
    gamma_1            = gamma_1,
    gamma_0_CI         = c(lo = gamma_0_lo, hi = gamma_0_hi),
    gamma_1_CI         = c(lo = gamma_1_lo, hi = gamma_1_hi),
    phi_CI             = c(lo = phi_lo, hi = phi_hi),
    interpretation     = interp,
    effective_se_power = 2 * beta1,
    data               = d_out,
    diagnostics        = list(
      n_input     = n_input,
      n_fit       = n_fit,
      n_removed   = n_removed,
      n_cells     = n_cells,
      rhat_max    = signif(rhat_max, 4),
      n_divergent = n_div,
      plate_sd    = if (has_plate) signif(plate_sd, 4) else NULL,
      weight      = w_diag
    ),
    cv_diagnostics     = cv_diag,
    formula            = list(
      location = paste("yi ~", loc_rhs),
      scale    = paste("sigma ~", sig_rhs)
    ),
    priors_used        = model_priors
  )
}
