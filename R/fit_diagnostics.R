# Diagnostics shared by the multiplicative and the additive fit, so the two
# models are judged on the same footing: did the sampler converge, did the data
# (rather than the prior) set the key scale parameter, and does the fitted
# variance model describe the residuals across the range of the measurement se.

# SD of a normal() prior, read from a brms prior object; NA if it is not a
# normal prior with a numeric SD.
.prior_sd <- function(prior) {
  p <- tryCatch(as.character(prior$prior[1]), error = function(e) NA_character_)
  if (is.na(p)) return(NA_real_)
  m <- regmatches(p, regexec("^\\s*normal\\(\\s*[^,]+,\\s*([0-9.eE+-]+)\\s*\\)\\s*$", p))[[1]]
  if (length(m) < 2) return(NA_real_)
  suppressWarnings(as.numeric(m[2]))
}

# Convergence: worst R-hat, number of divergent transitions, and the lowest
# effective sample size (bulk or tail) over all sampled parameters.
.sampling_diagnostics <- function(fit) {
  rhat_max <- suppressWarnings(max(brms::rhat(fit), na.rm = TRUE))
  if (!is.finite(rhat_max)) rhat_max <- NA_real_
  n_div <- tryCatch(as.integer(sum(brms::nuts_params(fit, pars = "divergent__")$Value)),
                    error = function(e) NA_integer_)
  ess_min <- tryCatch({
    dr <- posterior::as_draws_array(fit)
    dr <- posterior::subset_draws(
      dr, variable = setdiff(posterior::variables(dr), c("lp__", "lprior")))
    # posterior warns when it caps an ESS (very short chains); the cap is fine here
    s <- suppressWarnings(posterior::summarise_draws(dr, "ess_bulk", "ess_tail"))
    v <- c(s$ess_bulk, s$ess_tail)
    if (any(is.finite(v))) min(v, na.rm = TRUE) else NA_real_
  }, error = function(e) NA_real_)
  list(rhat_max = rhat_max, n_divergent = n_div, ess_min = ess_min)
}

# Posterior mean of the location (mu) for the rows the model was fitted to,
# from up to `ndraws` evenly spaced draws (deterministic, no RNG).
.posterior_mu <- function(fit, ndraws = 200L) {
  tryCatch({
    n <- posterior::ndraws(posterior::as_draws_array(fit, variable = "lp__"))
    ids <- unique(round(seq(1, n, length.out = min(ndraws, n))))
    colMeans(brms::posterior_linpred(fit, dpar = "mu", draw_ids = ids))
  }, error = function(e) NULL)
}

# Variance-model check. z = (y - mu) / sigma should have SD 1 everywhere if the
# variance model is right. resid_sd is a robust SD of z (1.4826 * MAD) over all
# rows; resid_sd_ratio is the robust SD in the lowest-se fifth over that in the
# highest-se fifth. A ratio above 1 means the low-se samples scatter more than
# the model allows, i.e. they are given too much weight; below 1, the high-se
# samples are. NA when there are fewer than 10 rows in either fifth.
.resid_check <- function(y, mu, sigma, se) {
  ok <- is.finite(y) & is.finite(mu) & is.finite(sigma) & sigma > 0 & is.finite(se)
  out <- list(resid_sd = NA_real_, resid_sd_ratio = NA_real_)
  if (sum(ok) < 20) return(out)
  z <- (y[ok] - mu[ok]) / sigma[ok]; s <- se[ok]
  rsd <- function(x) 1.4826 * stats::median(abs(x - stats::median(x)))
  out$resid_sd <- rsd(z)
  q <- stats::quantile(s, c(0.2, 0.8), names = FALSE)
  lo <- z[s <= q[1]]; hi <- z[s >= q[2]]
  if (length(lo) >= 10 && length(hi) >= 10 && rsd(hi) > 0)
    out$resid_sd_ratio <- rsd(lo) / rsd(hi)
  out
}
