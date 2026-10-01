# Tests for calculate_auc() (utils_auc.R) and tumor_doubling_time()
# (tumor_doubling_time.R).

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
