test_that("prepare_cv computes yi correctly", {
  df <- data.frame(
    predicted_concentration = c(10, 100, 1000, 0, -5, NA),
    se_concentration        = c(1, 5, 50, 1, 1, 1),
    pcov                    = c(0.1, 0.05, 0.05, 0.1, 0.1, 0.1)
  )
  result <- prepare_cv(df, pcov_col = "pcov")

  expect_equal(result$yi[1], log10(10))
  expect_equal(result$yi[2], log10(100))
  expect_equal(result$yi[3], log10(1000))
  expect_true(is.na(result$yi[4]))   # conc = 0
  expect_true(is.na(result$yi[5]))   # conc < 0
  expect_true(is.na(result$yi[6]))   # conc = NA
})

test_that("prepare_cv prefers pcov over se/conc", {
  df <- data.frame(
    predicted_concentration = c(100, 100),
    se_concentration        = c(10, 10),
    pcov                    = c(0.15, 0.15)
  )
  result <- prepare_cv(df, pcov_col = "pcov")

  # Should use pcov (0.15), not se/conc (10/100 = 0.10)
  expect_equal(result$cv_i[1], 0.15)
  expect_equal(result$cv_source[1], "pcov")
})

test_that("prepare_cv falls back to se/conc when pcov is NULL", {
  df <- data.frame(
    predicted_concentration = c(100, 200),
    se_concentration        = c(10, 20)
  )
  result <- prepare_cv(df, pcov_col = NULL)

  expect_equal(result$cv_i[1], 10 / 100)
  expect_equal(result$cv_i[2], 20 / 200)
  expect_equal(result$cv_source[1], "se_over_conc")
})

test_that("prepare_cv falls back to se/conc when pcov is NA", {
  df <- data.frame(
    predicted_concentration = c(100),
    se_concentration        = c(10),
    pcov                    = c(NA_real_)
  )
  result <- prepare_cv(df, pcov_col = "pcov")

  expect_equal(result$cv_i[1], 0.10)
  expect_equal(result$cv_source[1], "se_over_conc")
})

test_that("prepare_cv computes log_cv correctly", {
  df <- data.frame(
    predicted_concentration = c(100),
    se_concentration        = c(10),
    pcov                    = c(0.2)
  )
  result <- prepare_cv(df, pcov_col = "pcov")

  expect_equal(result$log_cv[1], log(0.2))
})

test_that("prepare_cv handles all-NA pcov gracefully", {
  df <- data.frame(
    predicted_concentration = c(100, 200),
    se_concentration        = c(10, 20),
    pcov                    = c(NA_real_, NA_real_)
  )
  result <- prepare_cv(df, pcov_col = "pcov")

  expect_equal(result$cv_i[1], 0.10)
  expect_equal(result$cv_i[2], 0.10)
  expect_true(all(result$cv_source == "se_over_conc"))
})
