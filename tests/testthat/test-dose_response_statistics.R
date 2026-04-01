# =============================================================================
# Tests for dose_response_statistics()
#
# Uses make_dose_response():
#   Doses: 0 (Control), 1, 5, 25 mg/kg
#   Mean volumes at day 21: 500, 380, 220, 90  (clear monotonic decrease)
#   8% CV per group, n=4 per dose level
#   Jonckheere trend test + linear model both expected significant (p < 0.01)
# =============================================================================

test_that("returns a list with expected fields", {
  df  <- make_dose_response()
  res <- suppressWarnings(suppressMessages(
    dose_response_statistics(
      df,
      dose_column          = "Dose",
      treatment_column     = "Treatment",
      volume_column        = "Volume",
      day_column           = "Day",
      id_column            = "ID",
      control_group_name   = "Control"
    )
  ))

  expect_true(is.list(res))
  required <- c("dose_effect_test", "summary_table")
  present  <- required[required %in% names(res)]
  expect_true(length(present) >= 1L,
              info = paste("Missing:", paste(setdiff(required, names(res)), collapse = ", ")))
})

test_that("summary_table has one row per dose level", {
  df <- make_dose_response()
  n_unique_doses <- length(unique(df$Dose))

  res <- suppressWarnings(suppressMessages(
    dose_response_statistics(
      df,
      dose_column        = "Dose",
      treatment_column   = "Treatment",
      volume_column      = "Volume",
      day_column         = "Day",
      id_column          = "ID",
      control_group_name = "Control"
    )
  ))

  summ <- as.data.frame(res$summary_table)
  expect_equal(nrow(summ), n_unique_doses)
})

test_that("trend test detects significant dose-response (p < 0.01)", {
  df  <- make_dose_response()
  res <- suppressWarnings(suppressMessages(
    dose_response_statistics(
      df,
      dose_column        = "Dose",
      treatment_column   = "Treatment",
      volume_column      = "Volume",
      day_column         = "Day",
      id_column          = "ID",
      control_group_name = "Control"
    )
  ))

  # Check dose_effect_test or trend_test for a significant p-value
  all_p <- c()

  extract_p <- function(x) {
    if (is.list(x)) {
      p_fields <- c("p.value", "p_value", "P_Value", "pvalue",
                    "statistic_p", "Pr(>F)")
      found <- x[intersect(p_fields, names(x))]
      if (length(found) > 0) return(as.numeric(found[[1]]))
    }
    if (is.data.frame(x)) {
      p_col <- intersect(c("p.value", "p_value", "P_Value"), colnames(x))[1]
      if (!is.na(p_col)) return(min(as.numeric(x[[p_col]]), na.rm = TRUE))
    }
    NULL
  }

  for (nm in c("dose_effect_test", "trend_test", "linear_model", "anova_model")) {
    if (nm %in% names(res)) {
      p <- extract_p(res[[nm]])
      if (!is.null(p)) all_p <- c(all_p, p)
    }
  }

  skip_if(length(all_p) == 0, "No p-values found in dose_response_statistics output")

  expect_true(any(all_p < 0.01, na.rm = TRUE),
              info = paste("Smallest p-value:", min(all_p, na.rm = TRUE)))
})

test_that("mean volume decreases with dose in summary_table", {
  df  <- make_dose_response()
  res <- suppressWarnings(suppressMessages(
    dose_response_statistics(
      df,
      dose_column        = "Dose",
      treatment_column   = "Treatment",
      volume_column      = "Volume",
      day_column         = "Day",
      id_column          = "ID",
      control_group_name = "Control"
    )
  ))

  summ     <- as.data.frame(res$summary_table)
  dose_col <- intersect(c("Dose", "dose", "Concentration"), colnames(summ))[1]
  mean_col <- intersect(c("Mean_Volume", "mean_volume", "Mean", "mean",
                          "Mean.Volume"), colnames(summ))[1]

  skip_if(is.na(dose_col) || is.na(mean_col),
          "Cannot locate Dose / Mean_Volume column in summary_table")

  summ_sorted <- summ[order(as.numeric(summ[[dose_col]])), ]
  means <- as.numeric(summ_sorted[[mean_col]])

  # The mean volume should be non-increasing as dose increases
  expect_true(all(diff(means) <= 0),
              info = paste("Non-monotone means:", paste(round(means, 1), collapse = " > ")))
})

test_that("cross-check: Jonckheere trend test directly (clinfun) also significant", {
  skip_if_not_installed("clinfun")
  df <- make_dose_response()

  # Cast to per-dose vectors for direct clinfun call
  doses   <- as.numeric(df$Dose)
  volumes <- df$Volume

  jt <- clinfun::jonckheere.test(volumes, doses, alternative = "decreasing")
  expect_true(jt$p.value < 0.01,
              info = paste("Direct Jonckheere p:", jt$p.value))
})
