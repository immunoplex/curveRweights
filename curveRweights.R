# ===============================================================================
# curveRweights.R
#
# curveRweights: Bayesian Precision Weighting from Calibration Curve Uncertainty
#
# A package for estimating observation-level precision weights from the
# posterior coefficient of variation (pcov) of calibration curve predictions,
# using a joint location-scale Bayesian regression model.
#
# CORE IDEA:
#   Immunoassay observations are measured via a calibration curve (e.g., 4PL).
#   Each observation has a posterior CV (pcov) that quantifies how precisely
#   the curve determines that concentration.  This package estimates the
#   power-law relationship between pcov and residual variance:
#
#       sigma_i = phi * pcov_i^beta1
#       w_i     = 1 / sigma_i^2
#
#   phi   = baseline scaling (exp of the log-scale intercept)
#   beta1 = precision exponent (how steeply pcov drives weight differentiation)
#
#   The estimation uses a brms location-scale model with a SATURATED location
#   (cell means for every experimental condition) so that the scale parameters
#   are identified purely from within-cell residual structure and are not
#   contaminated by any misspecified location effects.
#
# WORKFLOW:
#   1. Create the saturated cell-means factor externally:
#        dat$cell <- interaction(dat$Arm, dat$Timeperiod, drop = TRUE)
#
#   2. Fit the saturated weight model:
#        sw <- fit_saturated_weight(dat, cell_col = "cell", ...)
#
#   3. Apply weights to any subset:
#        dat_sub <- sw$data |> filter(Arm %in% c("TT","TdaP"), Timeperiod=="post3rd")
#        # dat_sub now has w_saturated and w_saturated_norm ready for rma() or svyranktest()
#
# DEPENDENCIES:
#   brms (>= 2.19.0), dplyr (>= 1.1.0), tibble, stats
#
# AUTHORS:
#   Michael Scot Zens
#
# LICENSE:
#   AGPL 3.0
#
# REFERENCES:
#   Bürkner P-C (2017). brms: An R Package for Bayesian Multilevel Models
#     Using Stan. Journal of Statistical Software, 80(1), 1-28.
#     DOI: 10.18637/jss.v080.i01
#
#   Bürkner P-C (2018). Advanced Bayesian Multilevel Modeling with the R
#     Package brms. The R Journal, 10(1), 395-411.
#     DOI: 10.32614/RJ-2018-017
#
#   Carpenter B, Gelman A, Hoffman MD, et al. (2017). Stan: A Probabilistic
#     Programming Language. Journal of Statistical Software, 76(1), 1-32.
#     DOI: 10.18637/jss.v076.i01
#
#   DerSimonian R, Laird N (1986). Meta-analysis in clinical trials.
#     Controlled Clinical Trials, 7(3), 177-188. PMID: 3802833
#
#   Higgins JPT, Thompson SG (2002). Quantifying heterogeneity in a
#     meta-analysis. Statistics in Medicine, 21(11), 1539-1558. PMID: 12111919
#
#   Viechtbauer W (2010). Conducting Meta-Analyses in R with the metafor
#     Package. Journal of Statistical Software, 36(3), 1-48.
#     DOI: 10.18637/jss.v036.i03
# ===============================================================================


# ---- 1. PREPARE CV -----------------------------------------------------------

#' Prepare the precision index (cv_i) and log-transformed precision (log_cv)
#'
#' Computes the precision index used as the predictor in the scale submodel.
#' Uses the stored posterior CV (pcov) from the calibration curve when available
#' and finite; falls back to se_concentration / predicted_concentration otherwise.
#'
#' The pcov is preferred because it correctly captures the non-Gaussian posterior
#' near the LLOQ and ULOQ where the 4PL curve is flat and the delta-method
#' se/conc approximation breaks down.
#'
#' @param df Data frame with observation-level data.
#' @param concentration_col Name of the predicted concentration column.
#' @param se_col Name of the SE of concentration column.
#' @param pcov_col Name of the posterior CV column.  Set to NULL to force
#'   computation from se_col / concentration_col.
#'
#' @return Data frame with added columns:
#'   \describe{
#'     \item{yi}{log10(predicted_concentration)}
#'     \item{cv_i}{Precision index: pcov when available, else se/conc}
#'     \item{log_cv}{log(cv_i)}
#'     \item{cv_source}{Character: "pcov" or "se_over_conc"}
#'   }
#'
#' @export
prepare_cv <- function(df,
                       concentration_col = "predicted_concentration",
                       se_col            = "se_concentration",
                       pcov_col          = "pcov") {

  # Determine pcov availability before entering mutate
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
          is.finite(.data[[se_col]])              & .data[[se_col]] > 0 &
          is.finite(.data[[concentration_col]])   & .data[[concentration_col]] > 0,
          .data[[se_col]] / .data[[concentration_col]],
          NA_real_
        )
      ),
      log_cv    = dplyr::if_else(is.finite(cv_i) & cv_i > 0, log(cv_i), NA_real_),
      cv_source = dplyr::if_else(
        if (use_pcov) is.finite(.data[[pcov_col]]) & .data[[pcov_col]] > 0
        else          rep(FALSE, dplyr::n()),
        "pcov", "se_over_conc"
      )
    )

  out
}


# ---- 2. DIAGNOSE CV VARIATION ------------------------------------------------

#' Diagnose whether the precision index has sufficient variation for scale
#' estimation
#'
#' The scale submodel log(sigma) ~ log_cv requires meaningful spread in log_cv
#' across observations.  If all observations have nearly identical pcov (e.g.,
#' all in the well-determined midrange), then beta1 is unidentifiable and the
#' model should fall back to an intercept-only sigma (uniform weights).
#'
#' @param df Data frame with cv_i and log_cv columns (from \code{prepare_cv()}).
#'
#' @return Named list:
#'   \describe{
#'     \item{n_finite}{Number of observations with finite cv_i > 0}
#'     \item{cv_min, cv_median, cv_max}{Summary statistics of cv_i}
#'     \item{log_cv_sd}{Standard deviation of log_cv}
#'     \item{log_cv_range}{Range of log_cv}
#'     \item{use_cv_slope}{Logical: TRUE if sd(log_cv) >= 0.05}
#'     \item{message}{Human-readable summary}
#'   }
#'
#' @export
diagnose_cv <- function(df) {
  ok   <- is.finite(df$cv_i) & df$cv_i > 0 & is.finite(df$log_cv)
  cv   <- df$cv_i[ok]
  lcv  <- df$log_cv[ok]
  n_ok <- sum(ok)

  if (n_ok < 5) {
    return(list(
      n_finite     = n_ok,
      cv_min       = NA_real_, cv_median = NA_real_, cv_max = NA_real_,
      log_cv_sd    = NA_real_, log_cv_range = NA_real_,
      use_cv_slope = FALSE,
      message      = paste("Only", n_ok, "finite cv_i values; need >= 5")
    ))
  }

  lcv_sd  <- sd(lcv)
  lcv_rng <- diff(range(lcv))
  use_sl  <- lcv_sd >= 0.05

  list(
    n_finite     = n_ok,
    cv_min       = signif(min(cv), 4),
    cv_median    = signif(median(cv), 4),
    cv_max       = signif(max(cv), 4),
    log_cv_sd    = signif(lcv_sd, 4),
    log_cv_range = signif(lcv_rng, 4),
    use_cv_slope = use_sl,
    message      = if (use_sl)
      paste0("OK: sd(log_cv) = ", signif(lcv_sd, 3),
             "; beta1 identifiable from ", n_ok, " observations")
    else
      paste0("WARNING: sd(log_cv) = ", signif(lcv_sd, 3),
             " < 0.05; beta1 not identifiable. ",
             "Falling back to intercept-only sigma (uniform weights)")
  )
}


# ---- 3. INTERPRET BETA1 -----------------------------------------------------

#' Interpret the estimated beta1 value
#'
#' Classifies the pcov-to-variance power law into interpretable regimes.
#'
#' @param beta1 Estimated beta1 (gamma_1 from the log-scale model).
#' @param beta1_lo Lower bound of 95% credible interval.
#' @param beta1_hi Upper bound of 95% credible interval.
#'
#' @return Character string describing the precision weighting regime.
#'
#' @details
#' The theoretical prediction from the delta method applied to the 4PL
#' calibration curve is beta1 = 1 (pcov is a direct proxy for residual SD).
#' Departures indicate:
#' \describe{
#'   \item{beta1 < 0.2}{Near-uniform weights; pcov carries little information}
#'   \item{beta1 in [0.2, 0.8]}{Compressed weighting; pcov is attenuated}
#'   \item{beta1 in [0.8, 1.2] with CI spanning this range}{Calibrated}
#'   \item{beta1 > 1.2}{Amplified weighting; pcov effect is stronger than theory}
#' }
#'
#' @export
interpret_beta1 <- function(beta1, beta1_lo = NA_real_, beta1_hi = NA_real_) {
  if (is.na(beta1)) return("not identified (uniform precision)")
  dplyr::case_when(
    !is.na(beta1_lo) & beta1_lo > 0.8 &
      !is.na(beta1_hi) & beta1_hi < 1.2 ~
      "calibrated (pcov ~ residual SD)",
    beta1 < 0.2    ~ "near-uniform weights",
    beta1 < 0.8    ~ "compressed precision weighting",
    beta1 > 1.2    ~ "amplified precision weighting",
    TRUE           ~ "moderate precision weighting"
  )
}


# ---- 4. COMPUTE WEIGHTS FROM SCALE PARAMETERS -------------------------------

#' Compute precision weights from estimated scale parameters
#'
#' Given gamma_0 (log phi) and gamma_1 (beta1), computes the per-observation
#' sigma, raw weight, and mean-normalised weight.
#'
#' @param df Data frame with log_cv column (from \code{prepare_cv()}).
#' @param gamma_0 Estimated intercept of the log-sigma model (log(phi)).
#' @param gamma_1 Estimated slope of the log-sigma model (beta1).
#'   Set to 0 for intercept-only (uniform weights).
#'
#' @return Data frame with added columns:
#'   \describe{
#'     \item{sigma_i}{exp(gamma_0 + gamma_1 * log_cv)}
#'     \item{w_saturated}{1 / sigma_i^2 (raw precision weight)}
#'     \item{w_saturated_norm}{w_saturated / mean(w_saturated) within the data}
#'   }
#'
#' @export
compute_saturated_weights <- function(df, gamma_0, gamma_1 = 0) {

  out <- df |>
    dplyr::mutate(
      sigma_i = dplyr::if_else(
        is.finite(log_cv),
        exp(gamma_0 + gamma_1 * log_cv),
        NA_real_
      ),
      w_saturated = dplyr::if_else(
        is.finite(sigma_i) & sigma_i > 0,
        1 / sigma_i^2,
        NA_real_
      )
    )

  # Mean-normalise: shape is identical to w_saturated, but mean = 1
  # for cross-group comparability
  w_valid <- out$w_saturated[is.finite(out$w_saturated) & out$w_saturated > 0]
  w_mean  <- if (length(w_valid) > 0) mean(w_valid) else NA_real_

  out <- out |>
    dplyr::mutate(
      w_saturated_norm = dplyr::if_else(
        is.finite(w_saturated) & w_saturated > 0 & is.finite(w_mean),
        w_saturated / w_mean,
        NA_real_
      )
    )

  out
}


# ---- 5. WEIGHT DIAGNOSTICS --------------------------------------------------

#' Compute summary diagnostics for a set of precision weights
#'
#' @param w Numeric vector of weights (w_saturated or w_saturated_norm).
#'
#' @return Named list:
#'   \describe{
#'     \item{n_obs}{Total observations (including NA)}
#'     \item{n_valid}{Observations with finite, positive weights}
#'     \item{n_eff}{Effective sample size: [sum(w)]^2 / sum(w^2)}
#'     \item{eff_ratio}{n_eff / n_valid (1.0 = uniform; <1 = weight variation)}
#'     \item{weight_ratio}{max(w) / min(w) among valid weights}
#'     \item{gini}{Gini coefficient of weights (0 = uniform, 1 = concentrated)}
#'   }
#'
#' @export
weight_diagnostics <- function(w) {
  n_obs   <- length(w)
  ok      <- is.finite(w) & w > 0
  w_valid <- w[ok]
  n_valid <- length(w_valid)

  if (n_valid < 2)
    return(list(n_obs = n_obs, n_valid = n_valid,
                n_eff = NA_real_, eff_ratio = NA_real_,
                weight_ratio = NA_real_, gini = NA_real_))

  n_eff <- sum(w_valid)^2 / sum(w_valid^2)
  wr    <- max(w_valid) / min(w_valid)

  # Gini coefficient: mean absolute difference / (2 * mean)
  ws    <- sort(w_valid)
  n     <- length(ws)
  gini  <- sum((2 * seq_len(n) - n - 1) * ws) / (n * sum(ws))

  list(
    n_obs        = n_obs,
    n_valid      = n_valid,
    n_eff        = round(n_eff, 1),
    eff_ratio    = round(n_eff / n_valid, 3),
    weight_ratio = signif(wr, 4),
    gini         = round(gini, 3)
  )
}


# ---- 6. FIT SATURATED WEIGHT MODEL ------------------------------------------

#' Fit a Bayesian location-scale model with saturated location and shared scale
#'
#' The location model uses cell means (one intercept per level of
#' \code{cell_col}) to absorb ALL systematic variation in the response.
#' The scale model estimates the power-law relationship between pcov and
#' residual variance across the entire dataset, shared across all cells.
#'
#' This separation is by design: the scale relationship is a property of
#' the assay (beads, optics, calibration curve), not of the experimental
#' conditions, so it should be estimated from the largest possible dataset
#' while the location model remains fully flexible.
#'
#' @param df Data frame with observation-level data.
#' @param cell_col Name of the saturated cell-means factor column.
#'   Created externally via \code{interaction(Arm, Timeperiod, drop=TRUE)}.
#'   The function uses \code{0 + cell_col} in the brms formula (no global
#'   intercept; one coefficient per cell).
#' @param concentration_col Predicted concentration column.
#' @param se_col SE of concentration column.
#' @param pcov_col Posterior CV column.  NULL = compute from se/conc.
#' @param plate_col Plate column for random intercept.  NULL = no plate RE.
#' @param prior_gamma0 Prior on gamma_0 (log phi).
#'   Default: normal(0, 1) — centered at phi=1, 95pct interval [0.14, 7.4].
#' @param prior_gamma1 Prior on gamma_1 (beta1).
#'   Default: normal(1, 0.5) — centered at theory-predicted value of 1,
#'   95pct interval [0, 2].
#' @param prior_location Prior on location (cell mean) coefficients.
#'   Default: normal(0, 2) — very permissive on log10-concentration scale.
#' @param prior_plate_sd Prior on plate RE standard deviation.
#'   Default: normal(0, 0.5) — half-normal via brms lower-bound convention.
#' @param iter MCMC iterations per chain.
#' @param warmup Warmup iterations per chain.
#' @param chains Number of MCMC chains.
#' @param cores Parallel cores.
#' @param adapt_delta Target HMC acceptance probability.
#' @param seed Random seed for reproducibility.
#'
#' @return Named list:
#'   \describe{
#'     \item{fit}{brms fit object}
#'     \item{phi}{exp(gamma_0): baseline scaling factor}
#'     \item{beta1}{gamma_1: precision exponent}
#'     \item{gamma_0}{Log-scale intercept estimate}
#'     \item{gamma_1}{Log-scale slope estimate}
#'     \item{gamma_0_CI}{95pct credible interval for gamma_0}
#'     \item{gamma_1_CI}{95pct credible interval for gamma_1 (beta1)}
#'     \item{phi_CI}{95pct credible interval for phi (back-transformed)}
#'     \item{interpretation}{Character: precision weighting regime}
#'     \item{data}{Input data frame with w_saturated and w_saturated_norm added}
#'     \item{diagnostics}{List of convergence and weight diagnostics}
#'     \item{cv_diagnostics}{Output of \code{diagnose_cv()}}
#'   }
#'
#' @details
#' MODEL STRUCTURE:
#'
#' Location (saturated cell means):
#'   yi ~ 0 + cell [+ (1 | plate)]
#'
#' Scale (shared across all cells):
#'   log(sigma_i) = gamma_0 + gamma_1 * log(cv_i)
#'
#' PRIORS (on the log scale):
#'
#' gamma_0 ~ Normal(0, 1):
#'   phi = exp(gamma_0).  Center at 0 implies phi = 1 a priori,
#'   meaning pcov is expected to be a reasonable (not exact) proxy
#'   for residual SD on the log10-concentration scale.
#'   The theoretical first-principles value is phi = 1/ln(10) = 0.434,
#'   i.e., gamma_0 = -0.83.  The default prior is wider than this to
#'   accommodate excess biological variance (phi > 1) which is common.
#'
#' gamma_1 ~ Normal(1, 0.5):
#'   beta1 = gamma_1.  Center at 1 is the delta-method prediction:
#'   if pcov is a perfect CV, then sigma proportional to pcov^1 on
#'   the log-concentration scale.  The SD of 0.5 allows the 95pct
#'   interval [0, 2], spanning from "pcov is uninformative" to
#'   "pcov effect is quadratic."
#'
#' For a sensitivity analysis, pass tighter or shifted priors:
#'   prior_gamma1 = brms::prior(normal(1.3, 0.7), class="b", dpar="sigma")
#'   uses the empirical mean from your data as the center.
#'
#' @examples
#' \dontrun{
#' # Create the cell factor externally
#' dat_all$cell <- interaction(dat_all$Arm, dat_all$timeperiod, drop = TRUE)
#'
#' # Fit the saturated weight model
#' sw <- fit_saturated_weight(
#'   df       = dat_all,
#'   cell_col = "cell",
#'   pcov_col = "pcov",
#'   plate_col = "plate"
#' )
#'
#' # Inspect scale estimates
#' cat("phi =", sw$phi, " beta1 =", sw$beta1, "\n")
#' cat(sw$interpretation, "\n")
#' cat("n_eff =", sw$diagnostics$weight$n_eff,
#'     "of", sw$diagnostics$weight$n_valid, "\n")
#'
#' # Use weights for a specific comparison
#' dat_post3rd <- sw$data |>
#'   dplyr::filter(Arm %in% c("TT", "TdaP"), timeperiod == "post3rd")
#'
#' # Weights are ready: dat_post3rd$w_saturated, dat_post3rd$w_saturated_norm
#' }
#'
#' @export
fit_saturated_weight <- function(
    df,
    cell_col          = "cell",
    concentration_col = "predicted_concentration",
    se_col            = "se_concentration",
    pcov_col          = "pcov",
    plate_col         = "plate",
    prior_gamma0      = brms::prior(normal(0, 1),   class = "Intercept",
                                     dpar = "sigma"),
    prior_gamma1      = brms::prior(normal(1, 0.5), class = "b",
                                     dpar = "sigma"),
    prior_location    = brms::prior(normal(0, 2),   class = "b"),
    prior_plate_sd    = brms::prior(normal(0, 0.5), class = "sd"),
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
         "Create it externally, e.g.: df$cell <- interaction(df$Arm, df$Time)")

  if (!concentration_col %in% names(df))
    stop("fit_saturated_weight: concentration_col '", concentration_col,
         "' not found.")

  if (!se_col %in% names(df))
    stop("fit_saturated_weight: se_col '", se_col, "' not found.")

  # ---- Prepare data --------------------------------------------------------
  # Ensure cell_col is a factor
  df[[cell_col]] <- factor(df[[cell_col]])
  n_cells        <- nlevels(df[[cell_col]])

  if (n_cells < 2)
    stop("fit_saturated_weight: cell_col '", cell_col,
         "' has fewer than 2 levels (", n_cells, "). ",
         "Need at least 2 cells for meaningful residuals.")

  # Compute yi, cv_i, log_cv
  d <- prepare_cv(df,
                  concentration_col = concentration_col,
                  se_col            = se_col,
                  pcov_col          = pcov_col)

  # Filter to observations usable for fitting
  d_fit <- d[is.finite(d$yi) & is.finite(d$cv_i) & d$cv_i > 0 &
               is.finite(d$log_cv), , drop = FALSE]

  n_input    <- nrow(df)
  n_fit      <- nrow(d_fit)
  n_removed  <- n_input - n_fit

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
  # Location: saturated cell means (0 + cell removes the global intercept;
  # each cell gets its own coefficient = its own mean)
  loc_rhs <- paste0("0 + ", cell_col)

  # Add plate random intercept if plate_col is present and has > 1 level
  has_plate <- !is.null(plate_col) &&
    plate_col %in% names(d_fit) &&
    dplyr::n_distinct(d_fit[[plate_col]]) > 1

  if (has_plate)
    loc_rhs <- paste0(loc_rhs, " + (1 | ", plate_col, ")")

  # Scale: log(sigma) ~ [1 +] log_cv   or   ~ 1 (intercept-only if flat cv)
  sig_rhs <- if (use_cv_slope) "log_cv" else "1"

  bf_formula <- brms::bf(
    as.formula(paste("yi ~", loc_rhs)),
    as.formula(paste("sigma ~", sig_rhs))
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
      family  = gaussian(),
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
            " > 1.05 — chains may not have converged")
  if (!is.na(n_div) && n_div > 0)
    warning("fit_saturated_weight: ", n_div, " divergent transitions")

  # ---- Extract scale parameters --------------------------------------------
  sig_fix <- brms::fixef(fit, dpar = "sigma")

  gamma_0    <- sig_fix["Intercept", "Estimate"]
  gamma_0_lo <- sig_fix["Intercept", "Q2.5"]
  gamma_0_hi <- sig_fix["Intercept", "Q97.5"]

  # Find log_cv row robustly (brms may prefix it)
  log_cv_row <- rownames(sig_fix)[grepl("log_cv", rownames(sig_fix),
                                         fixed = TRUE)]

  if (use_cv_slope && length(log_cv_row) == 1) {
    gamma_1    <- sig_fix[log_cv_row, "Estimate"]
    gamma_1_lo <- sig_fix[log_cv_row, "Q2.5"]
    gamma_1_hi <- sig_fix[log_cv_row, "Q97.5"]
  } else {
    gamma_1    <- 0
    gamma_1_lo <- NA_real_
    gamma_1_hi <- NA_real_
    if (use_cv_slope)
      warning("fit_saturated_weight: log_cv not found in sigma fixef. ",
              "Available rows: ",
              paste(rownames(sig_fix), collapse = ", "))
  }

  phi        <- exp(gamma_0)
  phi_lo     <- exp(gamma_0_lo)
  phi_hi     <- exp(gamma_0_hi)
  beta1      <- gamma_1
  beta1_lo   <- gamma_1_lo
  beta1_hi   <- gamma_1_hi
  interp     <- interpret_beta1(beta1, beta1_lo, beta1_hi)

  message("  phi = ", signif(phi, 4),
          "  [", signif(phi_lo, 3), ", ", signif(phi_hi, 3), "]")
  message("  beta1 = ", signif(beta1, 4),
          if (!is.na(beta1_lo))
            paste0("  [", signif(beta1_lo, 3), ", ", signif(beta1_hi, 3), "]")
          else "")
  message("  interpretation: ", interp)

  # ---- Compute weights on the FULL input data frame ------------------------
  # Weights are computed for ALL rows in the original data frame, not just
  # the rows that passed the fitting filter.  Rows with non-finite cv_i
  # get NA weights (transparent to downstream na.rm logic).
  d_out <- compute_saturated_weights(d, gamma_0 = gamma_0, gamma_1 = gamma_1)

  # Weight diagnostics on the fitted subset
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
    fit             = fit,
    phi             = phi,
    beta1           = beta1,
    gamma_0         = gamma_0,
    gamma_1         = gamma_1,
    gamma_0_CI      = c(lo = gamma_0_lo, hi = gamma_0_hi),
    gamma_1_CI      = c(lo = gamma_1_lo, hi = gamma_1_hi),
    phi_CI          = c(lo = phi_lo, hi = phi_hi),
    interpretation  = interp,
    effective_se_power = 2 * beta1,
    data            = d_out,
    diagnostics     = list(
      n_input   = n_input,
      n_fit     = n_fit,
      n_removed = n_removed,
      n_cells   = n_cells,
      rhat_max  = signif(rhat_max, 4),
      n_divergent = n_div,
      plate_sd  = if (has_plate) signif(plate_sd, 4) else NULL,
      weight    = w_diag
    ),
    cv_diagnostics  = cv_diag,
    formula         = list(
      location = paste("yi ~", loc_rhs),
      scale    = paste("sigma ~", sig_rhs)
    ),
    priors_used     = model_priors
  )
}


# ---- 7. BATCH FITTING ACROSS GROUPS -----------------------------------------

#' Fit saturated weight models across multiple groups
#'
#' Loops over groups defined by \code{group_vars} (e.g., antigen x source)
#' and fits one \code{fit_saturated_weight()} model per group.  Returns a
#' combined data frame with weights and a summary table of scale estimates.
#'
#' @param datg Full data frame.
#' @param group_vars Grouping column names (each unique combination gets
#'   its own phi and beta1).
#' @param cell_col Name of the saturated cell-means factor column.
#' @param ... Additional arguments passed to \code{fit_saturated_weight()}.
#'
#' @return Named list:
#'   \describe{
#'     \item{data}{Full data frame with w_saturated and w_saturated_norm}
#'     \item{scale_table}{Tibble: one row per group with phi, beta1, CIs,
#'       interpretation, n_fit, n_eff, weight_ratio}
#'     \item{fits}{Named list of brms fit objects (indexed by group label)}
#'     \item{diagnostics}{Named list of per-group diagnostic lists}
#'   }
#'
#' @examples
#' \dontrun{
#' dat_all$cell <- interaction(dat_all$Arm, dat_all$timeperiod, drop = TRUE)
#'
#' batch <- fit_saturated_weight_batch(
#'   datg       = dat_all,
#'   group_vars = c("antigen", "source"),
#'   cell_col   = "cell",
#'   pcov_col   = "pcov",
#'   plate_col  = "plate",
#'   iter = 4000, warmup = 1000, chains = 4, cores = 4
#' )
#'
#' # Scale summary
#' batch$scale_table
#'
#' # Weighted data for one comparison
#' batch$data |>
#'   dplyr::filter(Arm %in% c("TT", "TdaP"), timeperiod == "post3rd")
#' }
#'
#' @export
fit_saturated_weight_batch <- function(
    datg,
    group_vars = c("antigen", "source"),
    cell_col   = "cell",
    ...
) {
  # Coerce group columns to character for safe matching
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
    key_label <- paste(mapply(function(col) paste0(col, "=", keys[[col]]),
                              group_vars), collapse = ", ")

    message("\n[", i, "/", n_groups, "] ", key_label)

    # Subset
    df_sub <- datg
    for (col in group_vars) {
      val    <- as.character(keys[[col]])
      df_sub <- df_sub[!is.na(df_sub[[col]]) & df_sub[[col]] == val,
                       , drop = FALSE]
    }

    # Fit
    sw <- tryCatch(
      fit_saturated_weight(df = df_sub, cell_col = cell_col, ...),
      error = function(e) {
        warning("Group [", key_label, "] failed: ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(sw)) {
      # Return data with NA weights for failed groups
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
      phi_lo         = sw$phi_CI["lo"],
      phi_hi         = sw$phi_CI["hi"],
      beta1_lo       = sw$gamma_1_CI["lo"],
      beta1_hi       = sw$gamma_1_CI["hi"],
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


# ---- 8. APPLY WEIGHTS FROM A PREVIOUS FIT -----------------------------------

#' Apply previously estimated scale parameters to new data
#'
#' When you have scale estimates from a full dataset (via
#' \code{fit_saturated_weight()} or \code{fit_saturated_weight_batch()})
#' and want to apply them to a new or subsetted data frame without re-fitting.
#'
#' @param df New data frame.
#' @param scale_table A scale_table from \code{fit_saturated_weight_batch()},
#'   or a one-row tibble with at least gamma_0 and gamma_1 (or phi and beta1).
#' @param group_vars Column names to match between df and scale_table.
#'   If NULL, a single set of (phi, beta1) is applied to all rows.
#' @param concentration_col Predicted concentration column.
#' @param se_col SE of concentration column.
#' @param pcov_col Posterior CV column.
#'
#' @return Data frame with yi, cv_i, log_cv, sigma_i, w_saturated,
#'   w_saturated_norm columns added.
#'
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
    # Join group-specific scale parameters, then compute per-row
    for (col in group_vars) {
      if (col %in% names(df))          df[[col]]          <- as.character(df[[col]])
      if (col %in% names(scale_table)) scale_table[[col]] <- as.character(scale_table[[col]])
    }

    df <- dplyr::left_join(
      df,
      scale_table |> dplyr::select(dplyr::all_of(c(group_vars, "gamma_0", "gamma_1"))),
      by = group_vars
    )

    df <- df |>
      dplyr::mutate(
        sigma_i = dplyr::if_else(
          is.finite(log_cv) & is.finite(gamma_0) & is.finite(gamma_1),
          exp(gamma_0 + gamma_1 * log_cv),
          NA_real_
        ),
        w_saturated = dplyr::if_else(
          is.finite(sigma_i) & sigma_i > 0,
          1 / sigma_i^2,
          NA_real_
        )
      )

    # Normalise within each group
    df <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::mutate(
        w_saturated_norm = dplyr::if_else(
          is.finite(w_saturated) & w_saturated > 0,
          w_saturated / mean(w_saturated[is.finite(w_saturated) & w_saturated > 0],
                             na.rm = TRUE),
          NA_real_
        )
      ) |>
      dplyr::ungroup()

    # Clean up joined columns
    df <- df |> dplyr::select(-dplyr::any_of(c("gamma_0", "gamma_1")))
  }

  df
}
