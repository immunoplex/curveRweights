# Additive-model helpers: no Stan needed.

test_that("compute_additive_weights: bounded weights and reliability in (0, 1]", {
  d <- data.frame(se_meas = c(1e-6, 0.03, 0.1, 0.3, 1, 5))
  o <- compute_additive_weights(d, tau2 = 0.16, kappa = 1)
  expect_true(all(o$reliability > 0 & o$reliability <= 1))
  expect_true(all(o$w <= 1 / 0.16 + 1e-9))          # bounded above by 1/tau2
  expect_equal(o$reliability, 0.16 * o$w)           # reliability = tau2 * w
  expect_equal(o$w[1], 1 / 0.16, tolerance = 1e-6)  # perfectly measured -> cap
  expect_true(all(diff(o$w) < 0))                   # noisier -> smaller weight
  expect_equal(mean(o$w_norm), 1)
})

test_that("additive weight ratio is bounded where 1/v is not", {
  se <- c(0.01, 0.3)
  w_add  <- 1 / (0.16 + se^2)
  w_mult <- 1 / se^2
  expect_lt(max(w_add) / min(w_add), 2)       # <= (tau2 + v_max) / tau2
  expect_gt(max(w_mult) / min(w_mult), 100)
})

test_that("weights are invariant to rescaling in the sense of ratios", {
  d <- data.frame(se_meas = c(0.05, 0.2, 0.4))
  a <- compute_additive_weights(d, tau2 = 0.1)
  expect_equal(a$w / sum(a$w), a$w_norm / sum(a$w_norm))
})

test_that("censored rows get no w_norm and the mean over the rest is 1", {
  d <- data.frame(se_meas = c(0.05, 0.1, 0.2, 0.3),
                  censor_type = c("none", "none", "left", NA))
  o <- compute_additive_weights(d, tau2 = 0.1, censor_col = "censor_type")
  expect_true(is.na(o$w_norm[3]))
  expect_equal(mean(o$w_norm[c(1, 2, 4)]), 1)
  expect_true(is.finite(o$w[3]))
})

test_that("reliability interval is returned when draws are supplied", {
  d <- data.frame(se_meas = c(0.05, 0.2, 0.4))
  dr <- data.frame(tau2 = c(0.1, 0.15, 0.2), kappa = c(0.9, 1, 1.1), beta1 = 1)
  o <- compute_additive_weights(d, tau2 = 0.15, draws = dr)
  expect_true(all(o$reliability_lo <= o$reliability_hi))
  expect_true(all(o$reliability_lo > 0 & o$reliability_hi <= 1))
})

test_that("compute_additive_weights validates input", {
  expect_error(compute_additive_weights(data.frame(x = 1), tau2 = 0.1), "se_col")
  expect_error(compute_additive_weights(data.frame(se_meas = 1), tau2 = 0), "positive")
})

test_that("normalize_weights gives mean 1 within the subset and groups", {
  d <- data.frame(g = c("a", "a", "b", "b"), w = c(1, 3, 10, 30))
  o <- normalize_weights(d, by = "g")
  expect_equal(as.vector(tapply(o$w_norm, o$g, mean)), c(1, 1))
  sub <- d[d$w > 1, ]
  expect_equal(mean(normalize_weights(sub)$w_norm), 1)
  expect_error(normalize_weights(d, w_col = "nope"), "not found")
})

test_that("meas_variance returns v_meas and blanks censored rows", {
  d <- data.frame(v_meas = c(0.01, 0.2, 0.3),
                  censor_type = c("none", "left", NA))
  expect_equal(meas_variance(d), c(0.01, NA, 0.3))
  expect_equal(meas_variance(d, quantified_only = FALSE), d$v_meas)
  expect_error(meas_variance(data.frame(w = 1)), "v_meas")
})

test_that("reliability_correct shrinks toward the cell mean", {
  y <- c(1, 3, 10, 14); cell <- c("a", "a", "b", "b")
  full <- reliability_correct(y, cell, rep(1, 4))
  expect_equal(full, y)                                   # lambda = 1: unchanged
  none <- reliability_correct(y, cell, rep(1e-9, 4))
  expect_equal(none, c(2, 2, 12, 12), tolerance = 1e-6)   # lambda -> 0: cell mean
})

test_that("attenuation_correct divides by sqrt(lambda_x lambda_y)", {
  expect_equal(attenuation_correct(0.42, 0.85, 0.7), 0.42 / sqrt(0.85 * 0.7))
  expect_warning(attenuation_correct(0.9, 0.5, 0.5), "exceeds 1")
})

test_that("prepare_censoring builds flags, bounds and bound SEs", {
  d <- data.frame(
    curve_id = "c1",
    yi   = c(0.1, 0.8, 1.2, 1.6, NA, 0.3),
    cv_i = c(5, 0.12, 0.08, 0.15, NA, 4),
    censor_type  = c("left", "none", "none", "none", "right", "interval"),
    censor_lower = c(NA, NA, NA, NA, 70, 2),
    censor_upper = c(5, NA, NA, NA, NA, 6))
  o <- prepare_censoring(d, se_col = "cv_i")
  expect_equal(o$cens_flag, c(-1L, 0L, 0L, 0L, 1L, 2L))
  expect_equal(o$yi[1], log10(5))            # left: upper bound
  expect_equal(o$yi[5], log10(70))           # right: lower bound
  expect_equal(o$yi[6], log10(2)); expect_equal(o$yi2[6], log10(6))
  # SE at the limit comes from the nearest quantified sample, not the row's own
  expect_equal(o$se_meas[1], 0.12)           # nearest quantified yi to log10(5)=0.70 is 0.8
  expect_equal(o$se_meas[5], 0.15)           # nearest to log10(70)=1.85 is 1.6
  expect_equal(o$se_meas[2:4], d$cv_i[2:4])
})

test_that("prepare_censoring: no censor column means all quantified", {
  d <- data.frame(yi = c(1, 2), cv_i = c(0.1, 0.2))
  o <- prepare_censoring(d, censor_col = NULL, se_col = "cv_i")
  expect_equal(o$cens_flag, c(0L, 0L))
  expect_equal(o$se_meas, d$cv_i)
})

test_that("prepare_censoring rejects unknown types and flags missing bounds", {
  d <- data.frame(yi = 1, cv_i = 0.1, censor_type = "weird")
  expect_error(prepare_censoring(d, se_col = "cv_i"), "unknown censor_type")
  d2 <- data.frame(yi = c(1, 2), cv_i = c(0.1, 0.1),
                   censor_type = c("none", "left"),
                   censor_lower = NA_real_, censor_upper = NA_real_)
  expect_warning(o <- prepare_censoring(d2, se_col = "cv_i"), "lack the bound")
  expect_true(is.na(o$cens_flag[2]))
})

test_that(".dilution_shift handles text, placeholders and non-positive values", {
  x <- c("2000", "100", "__none__", NA, "0", "-5", "1")
  o <- .dilution_shift(x)
  expect_equal(o$shift, c(log10(2000), 2, 0, 0, 0, 0, 0))
  expect_equal(o$n_unknown, 4L)
  expect_equal(.dilution_shift(c(10, 1000))$n_unknown, 0L)
})
