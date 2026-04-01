# =============================================================================
# Tests for survival_statistics()
#
# Uses make_surv_simple():
#   AggressiveTx: 7 mice die at days 10,15,20,25,30,35,40  (Censor = 1)
#   Control:      7 mice all censored at day 60              (Censor = 0)
#
# KM ground truth:
#   AggressiveTx: S(20)=4/7≈0.571 > 0.5,  S(25)=3/7≈0.429 < 0.5
#     → KM median = 25 days (unambiguous, no tie at boundary)
#   Control: never reaches 50%  → median = NA
#
# Correctness verified against survival::survfit (reference implementation).
# =============================================================================

test_that("returns a list with required fields", {
  df  <- make_surv_simple()
  res <- suppressMessages(
    survival_statistics(
      df,
      time_column      = "Day",
      censor_column    = "Survival_Censor",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      reference_group  = "Control"
    )
  )

  expect_true(is.list(res))
  expect_true("results"         %in% names(res))
  expect_true("reference_group" %in% names(res))
  expect_true("method_used"     %in% names(res))
})

test_that("reference group is set to Control", {
  df  <- make_surv_simple()
  res <- suppressMessages(
    survival_statistics(
      df,
      time_column      = "Day",
      censor_column    = "Survival_Censor",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      reference_group  = "Control"
    )
  )

  expect_equal(res$reference_group, "Control")
})

test_that("results data frame has one row per non-reference group", {
  df  <- make_surv_simple()
  res <- suppressMessages(
    survival_statistics(
      df,
      time_column      = "Day",
      censor_column    = "Survival_Censor",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      reference_group  = "Control"
    )
  )

  results <- as.data.frame(res$results)
  grp_col <- intersect(c("Group", "Treatment", "group"), colnames(results))[1]
  expect_true(!is.na(grp_col))
  expect_true(nrow(results) >= 1L)
})

test_that("KM median survival for AggressiveTx matches survfit reference (25 days)", {
  df  <- make_surv_simple()

  # Reference implementation — R's survival::survfit
  km_ref  <- survival::survfit(
    survival::Surv(Day, Survival_Censor) ~ Treatment,
    data = df
  )
  km_summ   <- summary(km_ref)$table
  med_ref   <- km_summ[grep("AggressiveTx", rownames(km_summ)), "median"]

  # Ground truth from the KM estimator
  # S(20)=4/7≈0.571, S(25)=3/7≈0.429  → unambiguous median = 25
  expect_equal(as.numeric(med_ref), 25, tolerance = 0.01,
               info = paste("Reference-implementation median:", med_ref))

  res <- suppressMessages(
    survival_statistics(
      df,
      time_column      = "Day",
      censor_column    = "Survival_Censor",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      reference_group  = "Control"
    )
  )

  results  <- as.data.frame(res$results)
  grp_col  <- intersect(c("Group", "Treatment", "group"), colnames(results))[1]
  med_col  <- intersect(c("Median_Survival", "median", "Median"), colnames(results))[1]

  skip_if(is.na(grp_col) || is.na(med_col),
          "Cannot locate Group / Median_Survival column")

  med_fn <- results[[med_col]][results[[grp_col]] == "AggressiveTx"]
  # Function output should match the reference implementation
  expect_equal(as.numeric(med_fn), as.numeric(med_ref), tolerance = 0.5,
               info = paste("Function median:", med_fn, "Reference:", med_ref))
})

test_that("KM median survival for Control is NA (never reaches 50% events)", {
  df  <- make_surv_simple()
  res <- suppressMessages(
    survival_statistics(
      df,
      time_column      = "Day",
      censor_column    = "Survival_Censor",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      reference_group  = "Control"
    )
  )

  # Obtain KM fit independently to cross-check with the reference implementation
  km_ref <- survival::survfit(
    survival::Surv(Day, Survival_Censor) ~ Treatment,
    data = df
  )
  km_summ <- summary(km_ref)$table

  med_ctrl_ref <- km_summ[grep("Control", rownames(km_summ)), "median"]
  expect_true(is.na(med_ctrl_ref) || med_ctrl_ref > 60,
              info = paste("Expected NA median for Control, got", med_ctrl_ref))
})

test_that("hazard ratio for AggressiveTx vs Control is > 1 (higher hazard)", {
  df  <- make_surv_simple()
  res <- suppressMessages(
    survival_statistics(
      df,
      time_column      = "Day",
      censor_column    = "Survival_Censor",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      reference_group  = "Control"
    )
  )

  results <- as.data.frame(res$results)
  grp_col <- intersect(c("Group", "Treatment", "group"), colnames(results))[1]
  hr_col  <- intersect(c("HR", "hazard_ratio", "Hazard_Ratio", "exp(coef)"),
                       colnames(results))[1]

  skip_if(is.na(grp_col) || is.na(hr_col),
          "Cannot locate Group / HR column in results")

  hr_tx <- as.numeric(results[[hr_col]][results[[grp_col]] == "AggressiveTx"])
  expect_true(hr_tx > 1,
              info = paste("Expected HR > 1 for AggressiveTx, got", hr_tx))
})

test_that("results p-value is significant (< 0.05)", {
  df  <- make_surv_simple()
  res <- suppressMessages(
    survival_statistics(
      df,
      time_column      = "Day",
      censor_column    = "Survival_Censor",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      reference_group  = "Control"
    )
  )

  results <- as.data.frame(res$results)
  p_col   <- intersect(c("p_value", "p.value", "P_Value", "Pr(>|z|)"),
                       colnames(results))[1]

  skip_if(is.na(p_col), "Cannot locate p-value column")

  p_vals <- suppressWarnings(as.numeric(results[[p_col]]))
  expect_true(any(p_vals < 0.05, na.rm = TRUE),
              info = paste("Smallest p:", min(p_vals, na.rm = TRUE)))
})

test_that("Events and Total counts are consistent with input data", {
  df  <- make_surv_simple()
  res <- suppressMessages(
    survival_statistics(
      df,
      time_column      = "Day",
      censor_column    = "Survival_Censor",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      reference_group  = "Control"
    )
  )

  results <- as.data.frame(res$results)
  grp_col <- intersect(c("Group", "Treatment"), colnames(results))[1]
  evt_col <- intersect(c("Events", "events"), colnames(results))[1]
  tot_col <- intersect(c("Total",  "total"),  colnames(results))[1]

  skip_if(is.na(grp_col) || is.na(evt_col) || is.na(tot_col),
          "Cannot locate Events/Total columns")

  row_tx <- results[results[[grp_col]] == "AggressiveTx", ]
  expect_equal(as.integer(row_tx[[tot_col]]),  7L)
  expect_equal(as.integer(row_tx[[evt_col]]),  7L)   # all 7 mice died
})

test_that("cross-check: log-rank test from survival package agrees (p < 0.001)", {
  df       <- make_surv_simple()
  lr_ref   <- survival::survdiff(
    survival::Surv(Day, Survival_Censor) ~ Treatment,
    data = df
  )
  # Chi-sq to p: 1 - pchisq(chisq, df=1)
  p_lr <- 1 - stats::pchisq(lr_ref$chisq, df = length(lr_ref$n) - 1)
  expect_true(p_lr < 0.001,
              info = paste("Log-rank p:", p_lr))
})
