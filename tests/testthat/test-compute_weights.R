test_that("compute_saturated_weights computes sigma correctly", {
  df <- data.frame(log_cv = log(c(0.05, 0.10, 0.20)))

  gamma_0 <- log(1.5)  # phi = 1.5
  gamma_1 <- 1.0       # beta1 = 1

  result <- compute_saturated_weights(df, gamma_0 = gamma_0, gamma_1 = gamma_1)

  # sigma_i = exp(gamma_0 + gamma_1 * log_cv) = phi * cv_i^beta1
  expected_sigma <- 1.5 * c(0.05, 0.10, 0.20)^1.0
  expect_equal(result$sigma_i, expected_sigma, tolerance = 1e-10)
})

test_that("compute_saturated_weights computes w = 1/sigma^2", {
  df <- data.frame(log_cv = log(c(0.05, 0.10, 0.20)))

  gamma_0 <- 0    # phi = 1

  gamma_1 <- 1.0  # beta1 = 1

  result <- compute_saturated_weights(df, gamma_0 = gamma_0, gamma_1 = gamma_1)

  expected_w <- 1 / (c(0.05, 0.10, 0.20)^2)
  expect_equal(result$w_saturated, expected_w, tolerance = 1e-8)
})

test_that("compute_saturated_weights normalises to mean = 1", {
  df <- data.frame(log_cv = log(c(0.05, 0.10, 0.15, 0.20, 0.30)))

  result <- compute_saturated_weights(df, gamma_0 = 0.3, gamma_1 = 1.2)

  w_norm <- result$w_saturated_norm
  expect_equal(mean(w_norm), 1.0, tolerance = 1e-10)
})

test_that("compute_saturated_weights handles gamma_1 = 0 (uniform)", {
  df <- data.frame(log_cv = log(c(0.05, 0.10, 0.20, 0.40)))

  result <- compute_saturated_weights(df, gamma_0 = 0.5, gamma_1 = 0)

  # All sigma_i should be exp(0.5) = constant
  expect_equal(length(unique(result$sigma_i)), 1)
  # All weights should be equal
  expect_equal(length(unique(result$w_saturated)), 1)
  # All normalised weights should be 1
  expect_equal(result$w_saturated_norm, rep(1, 4), tolerance = 1e-10)
})

test_that("compute_saturated_weights handles NA log_cv", {
  df <- data.frame(log_cv = c(log(0.1), NA, log(0.2)))

  result <- compute_saturated_weights(df, gamma_0 = 0, gamma_1 = 1)

  expect_true(is.finite(result$sigma_i[1]))
  expect_true(is.na(result$sigma_i[2]))
  expect_true(is.finite(result$sigma_i[3]))
  expect_true(is.na(result$w_saturated[2]))
  expect_true(is.na(result$w_saturated_norm[2]))
})

test_that("higher pcov gives lower weight", {
  df <- data.frame(log_cv = log(c(0.03, 0.10, 0.30)))

  result <- compute_saturated_weights(df, gamma_0 = 0, gamma_1 = 1)

  expect_true(result$w_saturated[1] > result$w_saturated[2])
  expect_true(result$w_saturated[2] > result$w_saturated[3])
})

test_that("higher beta1 amplifies weight differentiation", {
  df <- data.frame(log_cv = log(c(0.05, 0.20)))

  # beta1 = 0.5 (compressed)
  r_low  <- compute_saturated_weights(df, gamma_0 = 0, gamma_1 = 0.5)
  ratio_low <- r_low$w_saturated[1] / r_low$w_saturated[2]

  # beta1 = 2.0 (amplified)
  r_high <- compute_saturated_weights(df, gamma_0 = 0, gamma_1 = 2.0)
  ratio_high <- r_high$w_saturated[1] / r_high$w_saturated[2]

  expect_true(ratio_high > ratio_low)
})
