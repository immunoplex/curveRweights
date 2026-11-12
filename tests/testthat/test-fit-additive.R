# Additive-model fitter; requires brms + a Stan backend (skipped otherwise).

skip_if_no_stan <- function() {
  testthat::skip_if_not_installed("brms")
  backend_ok <- requireNamespace("rstan", quietly = TRUE) ||
    requireNamespace("cmdstanr", quietly = TRUE)
  testthat::skip_if_not(backend_ok, "No Stan backend (rstan/cmdstanr)")
}

make_add_sim <- function(n_cell = 6, n_per = 50, tau = 0.4, kappa = 1,
                         cens_q = 0.08, seed = 11) {
  set.seed(seed)
  cells <- paste0("c", seq_len(n_cell))
  df <- expand.grid(rep = seq_len(n_per), cell = cells,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df$patientid <- rep(paste0("p", seq_len(n_per)), n_cell)
  mu  <- stats::setNames(seq(0.3, 2.3, length.out = n_cell), cells)
  pat <- stats::setNames(stats::rnorm(n_per, 0, 0.1), paste0("p", seq_len(n_per)))
  df$se_concentration <- exp(stats::rnorm(nrow(df), log(0.12), 0.7))
  y <- mu[df$cell] + pat[df$patientid] +
    stats::rnorm(nrow(df), 0, sqrt(tau^2 + (kappa * df$se_concentration)^2))
  lo <- stats::quantile(y, cens_q)
  df$censor_type  <- ifelse(y < lo, "left", "none")
  df$censor_upper <- ifelse(y < lo, 10^lo, NA_real_)
  df$censor_lower <- NA_real_
  df$predicted_concentration <- ifelse(y < lo, y - 0.5, y)  # extrapolated, ignored
  df$pcov <- NA_real_
  df$timeperiod <- df$cell
  df$truth_y <- y
  df
}

test_that("fit_additive_weight recovers tau and handles left-censoring", {
  skip_if_no_stan()
  df <- make_add_sim(tau = 0.4, kappa = 1)
  sw <- fit_additive_weight(df, cell_col = "cell", plate_col = NULL,
                            patient_col = "patientid", censor_col = "censor_type",
                            response_is_log10 = TRUE,
                            iter = 1000, warmup = 500, chains = 2, cores = 2,
                            seed = 5)
  expect_equal(sw$tau, 0.4, tolerance = 0.25)
  expect_true(sw$tau_CI[["lo"]] < 0.4 && 0.4 < sw$tau_CI[["hi"]])
  expect_true(sw$reliability > 0 && sw$reliability <= 1)
  expect_equal(sw$beta1, 1)                              # fixed by default
  expect_equal(sw$diagnostics$n_by_censoring$left, sum(df$censor_type == "left"))
  d <- sw$data
  expect_true(all(d$w <= 1 / sw$tau2 + 1e-9, na.rm = TRUE))   # bounded weights
  expect_true(all(is.na(d$w_norm[d$censor_type == "left"])))
  expect_equal(mean(d$w_norm, na.rm = TRUE), 1)
  expect_true(all(c("v_meas", "reliability", "reliability_lo", "reliability_hi",
                    "sigma_total", "bound") %in% names(d)))
})

test_that("censored likelihood reduces low-cell bias vs treating the bound as a point", {
  skip_if_no_stan()
  df <- make_add_sim(cens_q = 0.15, seed = 21)
  f <- function(col) fit_additive_weight(
    df, cell_col = "cell", plate_col = NULL, patient_col = "patientid",
    censor_col = col, response_is_log10 = TRUE,
    iter = 800, warmup = 400, chains = 2, cores = 2, seed = 3)
  with_c <- f("censor_type"); no_c <- f(NULL)
  truth <- mean(df$truth_y[df$cell == "c1"])
  est <- function(sw) unname(brms::fixef(sw$fit)["cellc1", "Estimate"])
  # the point-treated fit uses the (back-calculated, pulled-down) value
  expect_lt(abs(est(with_c) - truth), abs(est(no_c) - truth) + 0.05)
})

test_that("fit_precision_weights(model = \"additive\") returns the additive result", {
  skip_if_no_stan()
  df <- make_add_sim(n_per = 40, seed = 4)
  df$predicted_concentration <- df$truth_y
  df$censor_type <- "none"
  wd <- as_weight_data(df, design = "timeperiod", include_plate = FALSE)
  expect_true("patientid" %in% names(wd))                 # carried through
  pw <- fit_precision_weights(wd, model = "additive",
                              iter = 600, warmup = 300, chains = 2, cores = 2)
  expect_s3_class(pw, "precision_weights")
  expect_equal(pw$estimates$model, "additive")
  expect_true(all(c("v_meas", "tau2", "sigma", "w", "reliability", "w_norm",
                    "censor_type", "bound") %in% names(pw$weights)))
  expect_output(print(pw), "additive")
  expect_output(summary(pw), "tau\\^2")
  # meas_variance feeds metafor: v_meas, never the total variance
  expect_equal(meas_variance(pw), pw$weights$v_meas)
  # predict_weights reproduces the table for the fitted data
  pr <- predict_weights(pw, newdata = wd)
  expect_equal(pr$reliability, pw$weights$reliability, tolerance = 1e-8)
  # join_weights carries the new columns
  df$sampleid <- wd$sampleid <- paste0("S", seq_len(nrow(df)))
  pw2 <- pw; pw2$weights$sampleid <- wd$sampleid
  j <- join_weights(df, pw2, by = "sampleid")
  expect_true(all(c("w", "reliability", "v_meas") %in% names(j)))
})

test_that("multiplicative model remains available and flagged", {
  skip_if_no_stan()
  df <- make_add_sim(n_per = 30, seed = 8)
  df$predicted_concentration <- df$truth_y
  wd <- as_weight_data(df, design = "timeperiod", include_plate = FALSE)
  pw <- fit_precision_weights(wd, model = "multiplicative",
                              iter = 500, warmup = 250, chains = 2, cores = 2)
  expect_equal(pw$estimates$model, "multiplicative")
  expect_true(is.finite(pw$estimates$phi))
  expect_output(print(pw), "multiplicative")
})
