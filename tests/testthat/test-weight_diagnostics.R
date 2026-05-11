test_that("weight_diagnostics returns n_eff = n for uniform weights", {
  w <- rep(1, 100)
  result <- weight_diagnostics(w)

  expect_equal(result$n_obs, 100)
  expect_equal(result$n_valid, 100)
  expect_equal(result$n_eff, 100)
  expect_equal(result$eff_ratio, 1.0)
  expect_equal(result$weight_ratio, 1.0)
  expect_equal(result$gini, 0)
})

test_that("weight_diagnostics n_eff < n for non-uniform weights", {
  w <- c(rep(10, 10), rep(1, 90))
  result <- weight_diagnostics(w)

  expect_true(result$n_eff < 100)
  expect_true(result$eff_ratio < 1.0)
  expect_true(result$weight_ratio > 1)
  expect_true(result$gini > 0)
})

test_that("weight_diagnostics handles NA weights", {
  w <- c(1, 2, 3, NA, NA)
  result <- weight_diagnostics(w)

  expect_equal(result$n_obs, 5)
  expect_equal(result$n_valid, 3)
})

test_that("weight_diagnostics handles all-NA weights", {
  w <- rep(NA_real_, 5)
  result <- weight_diagnostics(w)

  expect_equal(result$n_obs, 5)
  expect_equal(result$n_valid, 0)
  expect_true(is.na(result$n_eff))
})

test_that("weight_diagnostics handles single weight", {
  w <- c(5)
  result <- weight_diagnostics(w)

  expect_equal(result$n_valid, 1)
  expect_true(is.na(result$n_eff))
})

test_that("weight_diagnostics weight_ratio is max/min", {
  w <- c(2, 4, 8)
  result <- weight_diagnostics(w)

  expect_equal(result$weight_ratio, signif(8 / 2, 4))
})

test_that("weight_diagnostics n_eff formula is correct", {
  # Manual check: w = c(3, 1, 1, 1)
  # sum(w) = 6, sum(w^2) = 9+1+1+1 = 12
  # n_eff = 36/12 = 3.0
  w <- c(3, 1, 1, 1)
  result <- weight_diagnostics(w)

  expect_equal(result$n_eff, 3.0)
})
