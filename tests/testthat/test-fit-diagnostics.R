# Diagnostics shared by both models (convergence, prior dominance, residual check).

test_that(".prior_sd reads the SD of a normal() prior", {
  testthat::skip_if_not_installed("brms")
  expect_equal(.prior_sd(brms::set_prior("normal(1, 0.5)", class = "b")), 0.5)
  expect_equal(.prior_sd(brms::set_prior("normal(0,2)", class = "b")), 2)
  expect_true(is.na(.prior_sd(brms::set_prior("student_t(3, 0, 1)", class = "b"))))
})

test_that(".resid_check finds a variance model that ignores real variation", {
  set.seed(1)
  n  <- 2000
  se <- exp(stats::rnorm(n, log(0.1), 0.7))
  tau <- 0.3
  y  <- stats::rnorm(n, 0, sqrt(tau^2 + se^2))
  # right variance model: scaled residuals have SD ~1 at every se
  ok <- .resid_check(y, rep(0, n), sqrt(tau^2 + se^2), se)
  expect_equal(ok$resid_sd, 1, tolerance = 0.1)
  expect_equal(ok$resid_sd_ratio, 1, tolerance = 0.2)
  # multiplicative-style model (sigma proportional to se, no tau): the low-se
  # samples scatter far more than the model allows, i.e. are over-weighted
  bad <- .resid_check(y, rep(0, n), 1.5 * se, se)
  expect_gt(bad$resid_sd_ratio, 2)
  # too few rows: NA, no error
  na <- .resid_check(y[1:10], rep(0, 10), se[1:10], se[1:10])
  expect_true(is.na(na$resid_sd) && is.na(na$resid_sd_ratio))
})

test_that("both models report the same diagnostic set", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not(requireNamespace("rstan", quietly = TRUE) ||
                          requireNamespace("cmdstanr", quietly = TRUE))
  set.seed(3)
  df <- expand.grid(rep = 1:40, cell = paste0("c", 1:4), KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
  df$patientid <- rep(paste0("p", 1:40), 4)
  df$se_concentration <- exp(stats::rnorm(nrow(df), log(0.12), 0.7))
  df$predicted_concentration <- stats::rnorm(
    nrow(df), as.integer(factor(df$cell)), sqrt(0.2^2 + df$se_concentration^2))
  df$censor_type <- "none"
  df$timeperiod <- df$cell
  wd <- as_weight_data(df, design = "timeperiod", include_plate = FALSE)
  keys <- c("rhat_max", "n_divergent", "ess_min", "ident_sd_ratio",
            "resid_sd", "resid_sd_ratio")
  for (m in c("multiplicative", "additive")) {
    pw <- quiet_ess(fit_precision_weights(wd, model = m, iter = 500, warmup = 250,
                                          chains = 2, cores = 2, seed = 4))
    dg <- pw$estimates$diagnostics
    expect_true(all(keys %in% names(dg)), info = m)
    expect_true(is.finite(dg$rhat_max) && dg$rhat_max >= 0.99, info = m)
    expect_true(is.finite(dg$ess_min) && dg$ess_min > 0, info = m)
    expect_true(is.finite(dg$n_divergent), info = m)
    expect_true(is.finite(dg$ident_sd_ratio) && dg$ident_sd_ratio > 0, info = m)
    expect_true(is.finite(dg$resid_sd) && dg$resid_sd > 0.5 && dg$resid_sd < 2, info = m)
    expect_true(is.finite(dg$resid_sd_ratio), info = m)
  }
})
