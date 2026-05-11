test_that("diagnose_cv returns use_cv_slope = TRUE for spread data", {
  df <- data.frame(
    cv_i   = c(0.03, 0.05, 0.10, 0.20, 0.40, 0.60),
    log_cv = log(c(0.03, 0.05, 0.10, 0.20, 0.40, 0.60))
  )
  result <- diagnose_cv(df)

  expect_true(result$use_cv_slope)
  expect_equal(result$n_finite, 6)
  expect_true(result$log_cv_sd > 0.05)
  expect_true(grepl("OK", result$message))
})

test_that("diagnose_cv returns use_cv_slope = FALSE for flat data", {
  # All pcov values nearly identical
  df <- data.frame(
    cv_i   = c(0.100, 0.101, 0.099, 0.100, 0.100),
    log_cv = log(c(0.100, 0.101, 0.099, 0.100, 0.100))
  )
  result <- diagnose_cv(df)

  expect_false(result$use_cv_slope)
  expect_true(result$log_cv_sd < 0.05)
  expect_true(grepl("WARNING", result$message))
})

test_that("diagnose_cv handles too few observations", {
  df <- data.frame(
    cv_i   = c(0.1, 0.2),
    log_cv = log(c(0.1, 0.2))
  )
  result <- diagnose_cv(df)

  expect_false(result$use_cv_slope)
  expect_equal(result$n_finite, 2)
  expect_true(is.na(result$log_cv_sd))
})

test_that("diagnose_cv handles NA values", {
  df <- data.frame(
    cv_i   = c(0.05, NA, 0.10, NA, 0.20, 0.30, 0.40),
    log_cv = log(c(0.05, NA, 0.10, NA, 0.20, 0.30, 0.40))
  )
  result <- diagnose_cv(df)

  expect_equal(result$n_finite, 5)  # 2 NAs excluded
  expect_true(result$use_cv_slope)
})

test_that("diagnose_cv summary statistics are correct", {
  vals <- c(0.05, 0.10, 0.20, 0.30, 0.50)
  df <- data.frame(cv_i = vals, log_cv = log(vals))
  result <- diagnose_cv(df)

  expect_equal(result$cv_min, signif(min(vals), 4))
  expect_equal(result$cv_median, signif(median(vals), 4))
  expect_equal(result$cv_max, signif(max(vals), 4))
})
