# Tests for calculate_auc() and the me_result S3 class
# (utils_auc.R + me_result.R)

# ---- calculate_auc() ---------------------------------------------------------

test_that("calculate_auc returns 0 for a single point", {
  # Single point → no interval → AUC is 0 or NA depending on implementation
  result <- calculate_auc(1, 100)
  expect_true(result == 0 || is.na(result))
})

test_that("calculate_auc computes correct trapezoidal area", {
  # Rectangle: 2 points, constant volume of 10, from time 0 to 5 → AUC = 50
  expect_equal(calculate_auc(c(0, 5), c(10, 10)), 50)
  # Triangle: from 0→10, time 0→2 → AUC = (0+10)/2 * 2 = 10

  expect_equal(calculate_auc(c(0, 2), c(0, 10)), 10)
})

test_that("calculate_auc handles unsorted time values", {
  # Should sort internally, same result either way
  auc_sorted   <- calculate_auc(c(0, 1, 2), c(0, 5, 10))
  auc_unsorted <- calculate_auc(c(2, 0, 1), c(10, 0, 5))
  expect_equal(auc_sorted, auc_unsorted)
})

test_that("calculate_auc removes NAs and still computes", {
  # NA in positions 2 and 4 → only uses (0,0), (2,10), (3,15)
  auc <- calculate_auc(c(0, NA, 2, 3, NA), c(0, NA, 10, 15, NA))
  # (0+10)/2*2 + (10+15)/2*1 = 10 + 12.5 = 22.5
  expect_equal(auc, 22.5)
})

test_that("calculate_auc returns NA for empty / all-NA input", {
  expect_true(is.na(calculate_auc(numeric(0), numeric(0))))
  expect_true(is.na(calculate_auc(c(NA, NA), c(NA, NA))))
})

test_that("calculate_auc validates mismatched lengths", {
  expect_error(calculate_auc(1:3, 1:2), "same length")
})

# ---- new_me_result() ---------------------------------------------------------

test_that("new_me_result creates a valid me_result object", {
  res <- new_me_result(
    analysis_type = "test",
    data   = data.frame(x = 1:3),
    results = list(a = 1),
    plots  = list(),
    summary = data.frame(stat = "mean", value = 2)
  )
  expect_s3_class(res, "me_result")
  expect_equal(res$analysis_type, "test")
  expect_true(!is.null(res$timestamp))
})

test_that("print.me_result runs without error", {
  res <- new_me_result("test", data.frame(), list(), list(), data.frame())
  expect_output(print(res), "Analysis Result")
})

test_that("summary.me_result returns the summary slot", {
  s <- data.frame(stat = "mean", value = 42)
  res <- new_me_result("test", data.frame(), list(), list(), s)
  expect_equal(summary(res), s)
})

test_that("plot.me_result returns plots", {
  p1 <- ggplot2::ggplot()
  res <- new_me_result("test", data.frame(), list(), list(growth = p1), data.frame())
  expect_identical(plot(res, "growth"), p1)
  expect_true(is.list(plot(res)))
})

# ---- export_diagnostics() ----------------------------------------------------

test_that("export_diagnostics works with an me_result wrapping lm", {
  m <- lm(mpg ~ cyl + wt, data = mtcars)
  res <- new_me_result("test", mtcars, list(model = m), list(), data.frame())
  df <- export_diagnostics(res)
  expect_true(is.data.frame(df))
  expect_true("AIC" %in% df$metric)
  expect_true("R_squared" %in% df$metric)
})

test_that("export_diagnostics writes CSV when file is supplied", {
  m  <- lm(mpg ~ hp, data = mtcars)
  res <- new_me_result("test", mtcars, list(model = m), list(), data.frame())
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf))
  export_diagnostics(res, file = tf)
  expect_true(file.exists(tf))
  csv <- read.csv(tf)
  expect_true("BIC" %in% csv$metric)
})

# ---- tumor_doubling_time() ---------------------------------------------------

test_that("tumor_doubling_time returns correct structure", {
  # Exponential growth: V0 * exp(beta * t) → doubling time = log(2) / beta
  set.seed(42)
  df <- data.frame(
    ID        = rep(c("M1", "M2"), each = 6),
    Treatment = rep(c("Control", "Drug"), each = 6),
    Day       = rep(0:5, 2),
    Volume    = c(100 * exp(0.1 * 0:5),   # growth rate ≈ 0.1 → doubling ≈ 6.93
                  100 * exp(0.05 * 0:5))   # growth rate ≈ 0.05 → doubling ≈ 13.86
  )
  res <- tumor_doubling_time(df, time_column = "Day", volume_column = "Volume",
                             id_column = "ID", treatment_column = "Treatment")
  expect_s3_class(res, "data.frame")
  expect_true("doubling_time" %in% names(res))
  expect_equal(nrow(res), 2)
  # Doubling time for M1 should be close to log(2)/0.1 ≈ 6.93

  m1 <- res[res$ID == "M1", ]
  expect_equal(m1$doubling_time, log(2) / 0.1, tolerance = 0.05)
})
