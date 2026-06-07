# These tests require brms + a Stan backend; skipped otherwise.

skip_if_no_stan <- function() {
  testthat::skip_if_not_installed("brms")
  backend_ok <- requireNamespace("rstan", quietly = TRUE) ||
    requireNamespace("cmdstanr", quietly = TRUE)
  testthat::skip_if_not(backend_ok, "No Stan backend (rstan/cmdstanr)")
}

make_sim <- function(n_cell = 6, n_per = 12, phi = 1.3, beta1 = 1.0,
                     seed = 1) {
  set.seed(seed)
  cells <- paste0("cell", seq_len(n_cell))
  df <- expand.grid(rep = seq_len(n_per), cell = cells,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # per-obs se on the log10 scale, spread enough that sd(log se) >> 0.05
  df$se_concentration <- exp(stats::rnorm(nrow(df), mean = log(0.15), sd = 0.6))
  mu <- stats::setNames(stats::rnorm(n_cell, 0, 1), cells)
  sigma <- phi * df$se_concentration^beta1
  df$predicted_concentration <- mu[df$cell] + stats::rnorm(nrow(df), 0, sigma)
  df$pcov      <- df$se_concentration * log(10) * 100
  df$pcov_pass <- TRUE
  df$timeperiod <- df$cell
  df
}

test_that("fit_precision_weights recovers phi and beta1 (se predictor)", {
  skip_if_no_stan()
  df <- make_sim(phi = 1.3, beta1 = 1.0)
  wd <- as_weight_data(df, design = "timeperiod", include_plate = FALSE)
  pw <- fit_precision_weights(wd, scale_predictor = "se",
                              iter = 1000, warmup = 500, chains = 2, cores = 2,
                              seed = 7)
  expect_s3_class(pw, "precision_weights")
  expect_true(is.finite(pw$estimates$phi))
  expect_true(is.finite(pw$estimates$beta1))
  # generous tolerances for short chains
  expect_equal(pw$estimates$phi,   1.3, tolerance = 0.5)
  expect_equal(pw$estimates$beta1, 1.0, tolerance = 0.5)
  expect_true(all(c("obs_id", "se", "sigma", "w") %in% names(pw$weights)))
})

test_that("se vs uncapped-pcov give the same beta1, phi off by the constant", {
  skip_if_no_stan()
  df  <- make_sim(phi = 1.2, beta1 = 1.1, seed = 3)
  wd  <- as_weight_data(df, design = "timeperiod", include_plate = FALSE)
  pw_se   <- fit_precision_weights(wd, scale_predictor = "se",
                                   iter = 800, warmup = 400, chains = 2, seed = 9)
  pw_pcov <- fit_precision_weights(wd, scale_predictor = "pcov",
                                   iter = 800, warmup = 400, chains = 2, seed = 9)
  expect_equal(pw_se$estimates$beta1, pw_pcov$estimates$beta1, tolerance = 0.2)
  # phi_pcov = phi_se * (ln10*100)^(-beta1) ; check log relationship loosely
  k <- log(10) * 100
  implied <- pw_se$estimates$gamma_0 - pw_pcov$estimates$gamma_0
  expect_equal(implied, pw_pcov$estimates$beta1 * log(k), tolerance = 0.5)
})

test_that("join_weights and predict_weights round-trip keys", {
  skip_if_no_stan()
  df <- make_sim()
  df$sampleid <- paste0("S", seq_len(nrow(df)))
  wd <- as_weight_data(df, design = "timeperiod", include_plate = FALSE)
  pw <- fit_precision_weights(wd, iter = 600, warmup = 300, chains = 2)
  joined <- join_weights(df, pw, by = "sampleid")
  expect_true(all(c("w", "sigma") %in% names(joined)))
  prof <- predict_weights(pw, newdata = wd)
  expect_true(all(c("sigma", "w") %in% names(prof)))
})
