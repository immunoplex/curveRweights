# i-spi-compute (the worker that I-SPI submits jobs to) calls these functions
# with their default arguments. These tests pin the 0.2.x defaults so that a
# new curveRweights release cannot change what the worker computes or writes.

test_that("0.2.x defaults are unchanged", {
  expect_equal(eval(formals(fit_precision_weights)$model)[1], "multiplicative")
  expect_equal(eval(formals(prepare_cv)$predictor)[1], "auto")
  expect_equal(eval(formals(fit_saturated_weight)$predictor)[1], "auto")
  expect_equal(eval(formals(apply_saturated_weights)$predictor)[1], "auto")
  # the worker's call shape: named MCMC args, no `model`
  expect_true(all(c("design", "scale_predictor", "plate_col", "...") %in%
                    names(formals(fit_precision_weights))))
})

test_that("the multiplicative result keeps the shape the worker reads", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not(requireNamespace("rstan", quietly = TRUE) ||
                          requireNamespace("cmdstanr", quietly = TRUE))
  set.seed(2)
  df <- expand.grid(rep = 1:20, cell = paste0("c", 1:4), KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
  df$se_concentration <- exp(stats::rnorm(nrow(df), log(0.15), 0.6))
  df$predicted_concentration <- stats::rnorm(nrow(df), as.integer(factor(df$cell)),
                                             df$se_concentration)
  df$timeperiod <- df$cell
  wd <- as_weight_data(df, design = "timeperiod", include_plate = FALSE)
  pw <- fit_precision_weights(wd, iter = 400, warmup = 200, chains = 2, cores = 2)
  e <- pw$estimates
  # exactly the fields worker_weights.R reads when it writes calib_weights_fit
  expect_true(all(c("phi", "phi_CI", "beta1", "beta1_CI", "interpretation",
                    "diagnostics") %in% names(e)))
  expect_true(all(c("lo", "hi") %in% names(e$phi_CI)))
  expect_true(all(c("n_fit", "weight") %in% names(e$diagnostics)))
  expect_true(all(c("n_eff", "weight_ratio") %in% names(e$diagnostics$weight)))
  # exactly the columns it writes to calib_weights
  expect_true(all(c("obs_id", "se", "pcov", "sigma", "w", "w_norm") %in%
                    names(pw$weights)))
  expect_equal(e$model, "multiplicative")
})
