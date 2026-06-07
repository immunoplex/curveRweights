test_that("prepare_cv: response_is_log10 avoids the double-log", {
  df <- data.frame(predicted_concentration = c(-1, 0, 1, 2),  # log10 scale
                   se_concentration        = c(0.3, 0.2, 0.1, 0.4))
  d <- prepare_cv(df, predictor = "se", response_is_log10 = TRUE)
  # yi should equal the supplied log10 values, NOT log10(of them)
  expect_equal(d$yi, df$predicted_concentration)
  # cv_i should be se_concentration verbatim under predictor = "se"
  expect_equal(d$cv_i, df$se_concentration)
  expect_true(all(d$cv_source[is.finite(d$cv_i)] == "se"))
})

test_that("prepare_cv: natural-scale default still logs the concentration", {
  df <- data.frame(predicted_concentration = c(0.1, 1, 10, 100),
                   se_concentration        = c(0.05, 0.1, 0.2, 0.3),
                   pcov                     = c(50, 20, 10, 8))
  d <- prepare_cv(df)  # auto + response_is_log10 = FALSE (legacy)
  expect_equal(d$yi, log10(df$predicted_concentration))
  # auto prefers pcov
  expect_equal(d$cv_i, df$pcov)
})

test_that("prepare_cv: se vs uncapped-pcov differ only by ln(10)*100 in log_cv", {
  df <- data.frame(predicted_concentration = c(-1, 0, 1),
                   se_concentration        = c(0.3, 0.2, 0.1))
  df$pcov <- df$se_concentration * log(10) * 100   # uncapped pcov
  d_se   <- prepare_cv(df, predictor = "se",   response_is_log10 = TRUE)
  d_pcov <- prepare_cv(df, predictor = "pcov", response_is_log10 = TRUE)
  # log_cv differs by an additive constant -> identical slope (beta1) identified
  diffs <- d_pcov$log_cv - d_se$log_cv
  expect_equal(diff(diffs), c(0, 0), tolerance = 1e-12)
  expect_equal(unique(round(diffs, 10)), round(log(log(10) * 100), 10))
})

test_that("as_weight_data.data.frame builds .cell, carries se, validates design", {
  df <- data.frame(
    sampleid = paste0("S", 1:6),
    predicted_concentration = c(-1, 0, 1, -1, 0, 1),  # log10
    se_concentration        = c(0.3, 0.2, 0.1, 0.25, 0.18, 0.12),
    pcov                    = c(70, 46, 23, 58, 41, 28),
    pcov_pass               = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    timeperiod              = c("t1", "t1", "t1", "t2", "t2", "t2"),
    cohort_arm              = c("a", "a", "a", "a", "a", "a")
  )
  wd <- suppressWarnings(
    as_weight_data(df, design = c("timeperiod", "cohort_arm"),
                   include_plate = FALSE)
  )
  expect_s3_class(wd, "weight_data")
  expect_true(all(c("obs_id", "concentration", "se", "pcov", ".cell") %in% names(wd)))
  expect_equal(wd$se, df$se_concentration)
  expect_identical(attr(wd, "design"), c("timeperiod", "cohort_arm"))
  expect_equal(nlevels(wd$.cell), 2L)  # t1|a and t2|a
})

test_that("as_weight_data errors on missing design columns", {
  df <- data.frame(predicted_concentration = 1:4, se_concentration = rep(.1, 4))
  expect_error(as_weight_data(df, design = "timeperiod"),
               "design column")
})

test_that("as_weight_data errors when every cell is a singleton", {
  df <- data.frame(
    predicted_concentration = c(-1, 0, 1, 2),
    se_concentration        = c(.3, .2, .1, .15),
    grp                     = c("a", "b", "c", "d")   # all unique -> singletons
  )
  expect_error(as_weight_data(df, design = "grp", include_plate = FALSE),
               "singleton")
})

test_that("as_weight_data warns but keeps out-of-range rows by default", {
  df <- data.frame(
    predicted_concentration = c(-1, 0, 1, -1, 0, 1),
    se_concentration        = c(.3, .2, .1, .25, .18, .12),
    pcov_pass               = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    g                       = c("a", "a", "a", "b", "b", "b")
  )
  expect_warning(
    wd <- as_weight_data(df, design = "g", include_plate = FALSE),
    "pcov_pass == FALSE"
  )
  expect_equal(nrow(wd), 6L)            # nothing dropped
  wd2 <- suppressMessages(
    as_weight_data(df, design = "g", include_plate = FALSE, drop_oor = TRUE)
  )
  expect_equal(nrow(wd2), 5L)           # one dropped
})
