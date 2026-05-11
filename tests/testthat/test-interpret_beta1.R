test_that("interpret_beta1 classifies calibrated correctly", {
  result <- interpret_beta1(1.0, beta1_lo = 0.85, beta1_hi = 1.15)
  expect_equal(result, "calibrated (pcov ~ residual SD)")
})

test_that("interpret_beta1 classifies amplified correctly", {
  result <- interpret_beta1(2.0, beta1_lo = 1.5, beta1_hi = 2.5)
  expect_equal(result, "amplified precision weighting")
})

test_that("interpret_beta1 classifies compressed correctly", {
  result <- interpret_beta1(0.5, beta1_lo = 0.3, beta1_hi = 0.7)
  expect_equal(result, "compressed precision weighting")
})

test_that("interpret_beta1 classifies near-uniform correctly", {
  result <- interpret_beta1(0.1)
  expect_equal(result, "near-uniform weights")
})

test_that("interpret_beta1 handles NA", {
  result <- interpret_beta1(NA_real_)
  expect_equal(result, "not identified (uniform precision)")
})

test_that("interpret_beta1 classifies moderate when CI straddles 1", {
  # beta1 = 1.0 but CI is wide [0.5, 1.5] — not fully calibrated
  result <- interpret_beta1(1.0, beta1_lo = 0.5, beta1_hi = 1.5)
  expect_equal(result, "moderate precision weighting")
})

test_that("interpret_beta1 requires both CI bounds for calibrated", {
  # Only one CI bound provided — can't confirm calibrated
  result <- interpret_beta1(1.0, beta1_lo = 0.9, beta1_hi = NA_real_)
  expect_equal(result, "moderate precision weighting")
})
