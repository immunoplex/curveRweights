#' Fit the additive measurement-error model (signal + measurement variance)
#'
#' Fits a Bayesian Gaussian model in which each observation's variance is the
#' \emph{sum} of true between-subject variance and measurement variance,
#' \deqn{\sigma_i^2 = \tau^2 + \kappa^2\, se_i^{2\beta_1},}
#' with a censored likelihood for samples outside the quantifiable range. It is
#' the default estimator of [fit_precision_weights()] from curveRweights 0.3.0.
#'
#' @section In plain terms:
#' An assay value varies between subjects for two different reasons. Subjects
#' truly differ (biology; this is the signal, variance \eqn{\tau^2}), and the
#' calibration curve measures each sample with some error (variance
#' \eqn{v_i}, known to us through \code{se_i}). The two add. The older
#' multiplicative model of curveRweights 0.2.x instead made the \emph{whole}
#' spread proportional to \code{se}, which implies that a precisely measured
#' sample has little biological variation, so it could be given a weight many
#' times larger than its neighbours. In the additive model the weight is at most
#' \eqn{1/\tau^2}; no sample can count for more than a perfectly measured one.
#'
#' \eqn{\kappa} is a calibration factor: \eqn{\kappa = 1} means the reported
#' \code{se} is exactly the measurement SD (the delta-method value). \eqn{\beta_1}
#' is fixed at 1 by default; \eqn{\kappa} and \eqn{\tau} are separated only
#' through the spread of \code{se} across samples, so freeing \eqn{\beta_1} as
#' well makes the model much harder to identify. Use \code{fix_beta1 = FALSE}
#' only with plenty of data over a wide range of \code{se}.
#'
#' @section Model:
#' Location: \code{yi ~ 0 + cell + (1|patient) + (1|plate) [+ (1|sample)]}, where
#' \code{cell} is the saturated design factor (e.g. arm by time point). The
#' patient term absorbs correlation between repeated samples of the same
#' subject (otherwise counted as noise); the optional sample term does the same
#' for several dilutions or replicates of one physical sample.
#'
#' Scale: \eqn{\sigma_i = \sqrt{\tau^2 + \kappa^2 se_i^{2\beta_1}}} via a
#' non-linear \code{sigma} formula. Rows with \code{censor_type} \code{left},
#' \code{right} or \code{interval} contribute a cumulative-probability term at
#' their censoring limit (see [prepare_censoring()]); quantified rows contribute
#' a normal density.
#'
#' Default priors: location \code{normal(0, 2)}; \eqn{\log\tau \sim
#' N(\log 0.3, 1)}; \eqn{\log\kappa \sim N(0, 0.5)} (centred on the delta-method
#' value); \eqn{\beta_1 \sim N(1, 0.25)} when free; group SDs
#' \code{normal(0, 0.5)}.
#'
#' @param df Data frame with observation-level data.
#' @param cell_col Name of the saturated cell factor column.
#' @param concentration_col Predicted concentration column.
#' @param se_col Measurement SE column (log10 scale; the curveR
#'   \code{se_concentration}).
#' @param plate_col Plate (or curve) grouping column; \code{NULL} for none.
#' @param patient_col Subject grouping column; \code{NULL} for none. Strongly
#'   recommended when subjects are measured more than once.
#' @param sample_col Optional column identifying a physical sample measured
#'   several times (dilutions/replicates).
#' @param censor_col Name of the censoring-type column (\code{"none"},
#'   \code{"left"}, \code{"right"}, \code{"interval"}). \code{NULL} (default)
#'   treats every row as quantified.
#' @param lower_col,upper_col Natural-scale censoring bounds.
#' @param se_bound_col Optional measurement SE at the censoring limit; see
#'   [prepare_censoring()].
#' @param dilution_col Optional dilution-factor column. If supplied, the response
#'   is converted to the final-concentration scale (\code{yi + log10(dilution)}),
#'   which is the scale curveRcore's censoring bounds are on. Supply it whenever
#'   you use curveRcore censoring bounds and dilutions differ.
#' @param response_is_log10 Is \code{concentration_col} already log10?
#' @param fix_beta1 Fix \eqn{\beta_1 = 1} (default) or estimate it.
#' @param prior_location,prior_ltau,prior_lkappa,prior_beta1,prior_group_sd
#'   \code{brms::prior} objects.
#' @param iter,warmup,chains,cores,adapt_delta,seed MCMC controls.
#' @param backend Passed to \code{brms::brm()}; \code{NULL} uses the brms default.
#' @param n_draws_weights Number of posterior draws used for the reliability
#'   interval (thinned from the fit).
#'
#' @return A list with \code{fit}, \code{tau}, \code{tau2}, \code{kappa},
#'   \code{beta1}, their credible intervals (\code{tau_CI}, \code{kappa_CI},
#'   \code{beta1_CI}), overall \code{reliability} and \code{reliability_CI}
#'   (\eqn{\tau^2/(\tau^2 + \bar v)} over quantified rows), \code{data} (the full
#'   input with the weight columns of [compute_additive_weights()]),
#'   \code{diagnostics}, \code{formula} and \code{priors_used}.
#'
#' @seealso [fit_precision_weights()], [compute_additive_weights()],
#'   [prepare_censoring()]
#'
#' @importFrom brms brm bf nlf fixef rhat nuts_params VarCorr
#' @export
fit_additive_weight <- function(
    df,
    cell_col          = "cell",
    concentration_col = "predicted_concentration",
    se_col            = "se_concentration",
    plate_col         = "plate",
    patient_col       = NULL,
    sample_col        = NULL,
    censor_col        = NULL,
    lower_col         = "censor_lower",
    upper_col         = "censor_upper",
    se_bound_col      = NULL,
    dilution_col      = NULL,
    response_is_log10 = FALSE,
    fix_beta1         = TRUE,
    prior_location    = brms::set_prior("normal(0, 2)", class = "b"),
    prior_ltau        = brms::set_prior("normal(-1.2, 1)", class = "b", nlpar = "ltau"),
    prior_lkappa      = brms::set_prior("normal(0, 0.5)",  class = "b", nlpar = "lkappa"),
    prior_beta1       = brms::set_prior("normal(1, 0.25)", class = "b", nlpar = "beta1"),
    prior_group_sd    = brms::set_prior("normal(0, 0.5)", class = "sd"),
    iter              = 2000,
    warmup            = 1000,
    chains            = 4,
    cores             = 4,
    adapt_delta       = 0.95,
    seed              = 42,
    backend           = NULL,
    n_draws_weights   = 1000
) {

  # ---- Validate ------------------------------------------------------------
  need <- c(cell_col, concentration_col, se_col)
  miss <- setdiff(need, names(df))
  if (length(miss))
    stop("fit_additive_weight: column(s) not found: ", paste(miss, collapse = ", "))
  for (g in c(plate_col, patient_col, sample_col))
    if (!is.null(g) && !g %in% names(df))
      stop("fit_additive_weight: grouping column '", g, "' not found.")
  if (!is.null(dilution_col) && !dilution_col %in% names(df))
    stop("fit_additive_weight: dilution_col '", dilution_col, "' not found.")

  df[[cell_col]] <- factor(df[[cell_col]])
  n_cells <- nlevels(df[[cell_col]])
  if (n_cells < 2)
    stop("fit_additive_weight: cell_col '", cell_col, "' has fewer than 2 levels.")

  # ---- Response, precision, censoring --------------------------------------
  d <- prepare_cv(df, concentration_col = concentration_col, se_col = se_col,
                  pcov_col = NULL, predictor = "se",
                  response_is_log10 = response_is_log10)
  if (!is.null(dilution_col)) {
    ds <- .dilution_shift(d[[dilution_col]])
    d$yi <- d$yi + ds$shift
    if (ds$n_unknown > 0)
      warning("fit_additive_weight: ", ds$n_unknown, " row(s) have no numeric ",
              "dilution in '", dilution_col, "' and were treated as dilution 1.")
  }
  d <- prepare_censoring(d, censor_col = censor_col, lower_col = lower_col,
                         upper_col = upper_col, se_col = "cv_i",
                         se_bound_col = se_bound_col)

  use <- is.finite(d$yi) & is.finite(d$se_meas) & d$se_meas > 0 &
    !is.na(d$cens_flag) & !is.na(d[[cell_col]])
  d_fit <- d[use, , drop = FALSE]
  n_input <- nrow(df); n_fit <- nrow(d_fit)
  n_cens  <- table(factor(d_fit$cens_flag, levels = c(0, -1, 1, 2),
                          labels = c("quantified", "left", "right", "interval")))
  message("fit_additive_weight: ", n_fit, " of ", n_input, " observations usable (",
          n_input - n_fit, " removed); ", n_cells, " cells; censored: ",
          sum(n_cens[-1]))
  if (n_fit < 20)
    stop("fit_additive_weight: only ", n_fit, " usable observations. Need >= 20.")
  if (sum(n_cens["quantified"]) < 10)
    stop("fit_additive_weight: fewer than 10 quantified observations; the ",
         "signal and measurement variance cannot be estimated.")

  # ---- Formula --------------------------------------------------------------
  has_grp <- function(g) !is.null(g) && g %in% names(d_fit) &&
    dplyr::n_distinct(d_fit[[g]]) > 1
  loc_rhs <- paste0("0 + ", cell_col)
  groups  <- c(patient = patient_col, plate = plate_col, sample = sample_col)
  groups  <- groups[vapply(groups, has_grp, logical(1))]
  for (g in groups) loc_rhs <- paste0(loc_rhs, " + (1 | ", g, ")")
  any_cens <- any(d_fit$cens_flag != 0L)
  lhs <- if (any_cens) "yi | cens(cens_flag, yi2)" else "yi"

  sig <- if (fix_beta1)
    "sigma ~ log(sqrt(exp(2 * ltau) + exp(2 * lkappa) * se_meas^2))"
  else
    "sigma ~ log(sqrt(exp(2 * ltau) + exp(2 * lkappa) * se_meas^(2 * beta1)))"
  bf_formula <- if (fix_beta1)
    brms::bf(stats::as.formula(paste(lhs, "~", loc_rhs)),
             brms::nlf(stats::as.formula(sig)),
             ltau ~ 1, lkappa ~ 1)
  else
    brms::bf(stats::as.formula(paste(lhs, "~", loc_rhs)),
             brms::nlf(stats::as.formula(sig)),
             ltau ~ 1, lkappa ~ 1, beta1 ~ 1)
  message("  location: ", lhs, " ~ ", loc_rhs)
  message("  scale:    ", sig)

  priors <- c(prior_location, prior_ltau, prior_lkappa,
              if (!fix_beta1) prior_beta1,
              if (length(groups)) prior_group_sd)

  # ---- Fit ------------------------------------------------------------------
  message("  fitting brms model (", iter, " iter, ", chains, " chains)...")
  args <- list(formula = bf_formula, data = d_fit, family = stats::gaussian(),
               prior = priors, iter = iter, warmup = warmup, chains = chains,
               cores = cores, seed = seed,
               control = list(adapt_delta = adapt_delta), silent = 2, refresh = 0)
  if (!is.null(backend)) args$backend <- backend
  fit <- tryCatch(do.call(brms::brm, args),
                  error = function(e)
                    stop("fit_additive_weight: brms::brm() failed: ",
                         conditionMessage(e)))

  rhat_max <- max(brms::rhat(fit), na.rm = TRUE)
  n_div <- tryCatch(sum(brms::nuts_params(fit, pars = "divergent__")$Value),
                    error = function(e) NA_integer_)
  if (!is.na(rhat_max) && rhat_max > 1.05)
    warning("fit_additive_weight: Rhat = ", signif(rhat_max, 4),
            " > 1.05 --- chains may not have converged")
  if (!is.na(n_div) && n_div > 0)
    warning("fit_additive_weight: ", n_div, " divergent transitions")

  # ---- Parameters -----------------------------------------------------------
  vars <- c("b_ltau_Intercept", "b_lkappa_Intercept",
            if (!fix_beta1) "b_beta1_Intercept")
  dr <- as.data.frame(brms::as_draws_df(fit, variable = vars))
  ltau <- dr$b_ltau_Intercept; lkap <- dr$b_lkappa_Intercept
  beta_d <- if (fix_beta1) rep(1, nrow(dr)) else dr$b_beta1_Intercept
  tau_d  <- exp(ltau); kappa_d <- exp(lkap)
  ci <- function(x) stats::setNames(
    stats::quantile(x, c(0.025, 0.975), names = FALSE), c("lo", "hi"))
  tau   <- exp(mean(ltau));  kappa <- exp(mean(lkap))
  beta1 <- mean(beta_d)
  tau2  <- tau^2

  keep <- unique(round(seq(1, nrow(dr), length.out = min(n_draws_weights, nrow(dr)))))
  draws_w <- data.frame(tau2 = tau_d[keep]^2, kappa = kappa_d[keep],
                        beta1 = beta_d[keep])

  # identifiability: if the data barely move log(kappa) away from its prior the
  # signal/measurement split is being decided by the prior.
  sd_prior <- 0.5
  kappa_ratio <- stats::sd(lkap) / sd_prior
  ident <- if (kappa_ratio > 0.8)
    paste0("kappa is weakly identified (posterior SD of log kappa is ",
           round(100 * kappa_ratio), "% of its prior SD): the split between ",
           "signal and measurement variance is largely set by the prior. ",
           "Consider more samples, or replicate measurements of the same sample.")
  else "kappa identified from the data"

  message("  tau = ", signif(tau, 3), " [", signif(ci(tau_d)[1], 3), ", ",
          signif(ci(tau_d)[2], 3), "]   kappa = ", signif(kappa, 3), " [",
          signif(ci(kappa_d)[1], 3), ", ", signif(ci(kappa_d)[2], 3), "]")
  message("  ", ident)

  # ---- Weights on the full input frame --------------------------------------
  d_out <- compute_additive_weights(d, tau2 = tau2, kappa = kappa, beta1 = beta1,
                                    se_col = "se_meas",
                                    censor_col = if (!is.null(censor_col) &&
                                                     censor_col %in% names(d)) censor_col,
                                    draws = draws_w)
  q_rows <- d_out$cens_flag %in% 0L & is.finite(d_out$v_meas)
  vbar_d <- vapply(seq_len(nrow(draws_w)), function(k)
    mean(draws_w$kappa[k]^2 * d_out$se_meas[q_rows]^(2 * draws_w$beta1[k])),
    numeric(1))
  rel_d <- draws_w$tau2 / (draws_w$tau2 + vbar_d)
  rel   <- tau2 / (tau2 + mean(d_out$v_meas[q_rows]))
  w_diag <- weight_diagnostics(d_out$w_norm)
  message("  overall reliability = ", signif(rel, 3), "   weight_ratio = ",
          w_diag$weight_ratio, "   n_eff = ", w_diag$n_eff)

  vc <- brms::VarCorr(fit)
  gsd <- vapply(groups, function(g)
    if (!is.null(vc[[g]])) unname(vc[[g]]$sd[1, "Estimate"]) else NA_real_,
    numeric(1))

  list(
    fit = fit, tau = tau, tau2 = tau2, kappa = kappa, beta1 = beta1,
    tau_CI = ci(tau_d), kappa_CI = ci(kappa_d),
    beta1_CI = if (fix_beta1) c(lo = NA_real_, hi = NA_real_) else ci(beta_d),
    reliability = rel, reliability_CI = ci(rel_d),
    data = d_out,
    diagnostics = list(
      n_input = n_input, n_fit = n_fit, n_removed = n_input - n_fit,
      n_cells = n_cells, n_by_censoring = as.list(n_cens),
      rhat_max = signif(rhat_max, 4), n_divergent = n_div,
      group_sd = if (length(gsd)) signif(gsd, 4) else NULL,
      fix_beta1 = fix_beta1, identifiability = ident,
      kappa_sd_ratio = signif(kappa_ratio, 3), weight = w_diag),
    formula = list(location = paste(lhs, "~", loc_rhs), scale = sig),
    priors_used = priors
  )
}


# log10 dilution shift that puts a pre-dilution log10 response on the
# final-concentration scale. Text such as "2000" is accepted; a missing or
# non-numeric dilution (e.g. the '__none__' placeholder stored by calib_samples)
# counts as dilution 1 rather than silently dropping the row.
.dilution_shift <- function(x) {
  d   <- suppressWarnings(as.numeric(as.character(x)))
  bad <- !is.finite(d) | d <= 0
  shift <- numeric(length(d))
  shift[!bad] <- log10(d[!bad])
  list(shift = shift, n_unknown = sum(bad))
}
