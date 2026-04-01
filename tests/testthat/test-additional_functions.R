# =============================================================================
# Tests for additional exported functions:
#   - analyze_drug_synergy_over_time()
#   - generate_summary_statistics()
#   - prepare_dose_data()
#   - repeated_measures_anova()
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: multi-timepoint synergy data for analyze_drug_synergy_over_time
#
# Extends make_combo_four_group() concept to use the same group names as the
# synergy fixtures (Control, DrugA, DrugB, DrugA+DrugB) with multiple time
# points (days 0, 7, 14, 21).
# ---------------------------------------------------------------------------
make_synergy_multi_timepoint <- function() {
  set.seed(42)
  groups <- c("Control", "DrugA", "DrugB", "DrugA+DrugB")
  # Control grows fast; combo grows slowest => increasing synergy over time

  slopes <- c(0.12, 0.07, 0.08, 0.02)
  n_per  <- 6
  days   <- c(0, 7, 14, 21)

  do.call(rbind, lapply(seq_along(groups), function(g) {
    do.call(rbind, lapply(seq_len(n_per), function(m) {
      prefix <- substr(gsub("[^A-Za-z0-9]", "", groups[g]), 1, 3)
      id     <- paste0(prefix, sprintf("%02d", m))
      cage   <- paste0(prefix, ceiling(m / 2))
      noise  <- rnorm(length(days), 0, 0.05)
      vol    <- exp(log(200) + slopes[g] * days + noise)
      data.frame(
        ID        = id,
        Cage      = cage,
        Treatment = groups[g],
        Day       = days,
        Volume    = round(vol, 2),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

# ===========================================================================
# 1. analyze_drug_synergy_over_time
# ===========================================================================

call_synergy_ot <- function(df) {
  suppressWarnings(suppressMessages(
    analyze_drug_synergy_over_time(
      df,
      treatment_column = "Treatment",
      volume_column    = "Volume",
      time_column      = "Day",
      drug_a_name      = "DrugA",
      drug_b_name      = "DrugB",
      combo_name       = "DrugA+DrugB",
      control_name     = "Control",
      verbose          = FALSE
    )
  ))
}

test_that("analyze_drug_synergy_over_time returns expected list structure", {
  df  <- make_synergy_multi_timepoint()
  res <- call_synergy_ot(df)

  expect_true(is.list(res))
  required <- c("timepoint_results", "synergy_summary",
                "peak_ci_synergy", "peak_bliss_synergy",
                "drug_a_name", "drug_b_name", "combo_name")
  expect_true(all(required %in% names(res)),
              info = paste("Missing:", paste(setdiff(required, names(res)), collapse = ", ")))
})

test_that("analyze_drug_synergy_over_time synergy_summary is a data frame with expected columns", {
  df  <- make_synergy_multi_timepoint()
  res <- call_synergy_ot(df)

  ss <- res$synergy_summary
  expect_s3_class(ss, "data.frame")
  # At least some of the key columns should be present
  expected_cols <- c("Time_Point", "TGI_Combo", "Bliss_Expected_TGI",
                     "Combination_Index", "Synergy_Assessment")
  present <- expected_cols[expected_cols %in% colnames(ss)]
  expect_true(length(present) >= 3,
              info = paste("Missing columns:", paste(setdiff(expected_cols, colnames(ss)), collapse = ", ")))
})

test_that("analyze_drug_synergy_over_time analyses multiple time points", {
  df  <- make_synergy_multi_timepoint()
  res <- call_synergy_ot(df)

  # Day 0 may be skipped (all volumes ~ equal), but we should have at least 2 time points
  expect_true(nrow(res$synergy_summary) >= 2)
  expect_true(length(res$timepoint_results) >= 2)
})

test_that("analyze_drug_synergy_over_time peak synergy rows are single-row data frames", {
  df  <- make_synergy_multi_timepoint()
  res <- call_synergy_ot(df)

  expect_equal(nrow(res$peak_ci_synergy), 1)
  expect_equal(nrow(res$peak_bliss_synergy), 1)
})

test_that("analyze_drug_synergy_over_time respects min/max time point filters", {
  df  <- make_synergy_multi_timepoint()
  res <- suppressWarnings(suppressMessages(
    analyze_drug_synergy_over_time(
      df,
      drug_a_name  = "DrugA",
      drug_b_name  = "DrugB",
      combo_name   = "DrugA+DrugB",
      min_time_point = 7,
      max_time_point = 14,
      verbose      = FALSE
    )
  ))

  # All analysed time points should fall within [7, 14]
  expect_true(all(res$synergy_summary$Time_Point >= 7))
  expect_true(all(res$synergy_summary$Time_Point <= 14))
})

test_that("analyze_drug_synergy_over_time stores group names", {
  df  <- make_synergy_multi_timepoint()
  res <- call_synergy_ot(df)

  expect_equal(res$drug_a_name, "DrugA")
  expect_equal(res$drug_b_name, "DrugB")
  expect_equal(res$combo_name, "DrugA+DrugB")
})

# ===========================================================================
# 2. generate_summary_statistics
# ===========================================================================

test_that("generate_summary_statistics returns summary with expected columns", {
  df <- make_dose_response()
  ss <- suppressMessages(generate_summary_statistics(df, dose_column = "Dose", volume_column = "Volume"))

  expect_s3_class(ss, "data.frame")
  expected_cols <- c("mean_volume", "median_volume", "sd_volume", "n", "sem_volume")
  expect_true(all(expected_cols %in% colnames(ss)),
              info = paste("Missing:", paste(setdiff(expected_cols, colnames(ss)), collapse = ", ")))
})

test_that("generate_summary_statistics groups by dose", {
  df <- make_dose_response()
  ss <- suppressMessages(generate_summary_statistics(df, dose_column = "Dose", volume_column = "Volume"))

  # Should have one row per unique dose level (0, 1, 5, 25)
  expect_equal(nrow(ss), length(unique(df$Dose)))
})

test_that("generate_summary_statistics computes correct n per group", {
  df <- make_dose_response()
  ss <- suppressMessages(generate_summary_statistics(df, dose_column = "Dose", volume_column = "Volume"))

  # Each dose level has 4 mice
  expect_true(all(ss$n == 4))
})

test_that("generate_summary_statistics mean decreases with dose", {
  df <- make_dose_response()
  ss <- suppressMessages(generate_summary_statistics(df, dose_column = "Dose", volume_column = "Volume"))

  ss <- ss[order(ss$Dose), ]
  # Higher doses -> lower mean volume for this fixture
  expect_true(ss$mean_volume[1] > ss$mean_volume[nrow(ss)])
})

# ===========================================================================
# 3. prepare_dose_data
# ===========================================================================

test_that("prepare_dose_data returns a data frame", {
  df <- make_dose_response()
  ad <- prepare_dose_data(df, dose_column = "Dose", treatment_column = "Treatment",
                          volume_column = "Volume", day_column = "Day", id_column = "ID")
  expect_s3_class(ad, "data.frame")
})

test_that("prepare_dose_data ensures Dose column is numeric", {
  df <- make_dose_response()
  df$Dose <- as.character(df$Dose)  # force character
  ad <- prepare_dose_data(df, dose_column = "Dose", treatment_column = "Treatment",
                          volume_column = "Volume", day_column = "Day", id_column = "ID")
  expect_true(is.numeric(ad$Dose))
})

test_that("prepare_dose_data filters to specified time point", {
  # Create multi-timepoint dose data
  df <- make_dose_response()
  df2 <- df
  df2$Day <- 14L
  df2$Volume <- df2$Volume * 0.8
  multi <- rbind(df, df2)

  ad <- prepare_dose_data(multi, dose_column = "Dose", treatment_column = "Treatment",
                          volume_column = "Volume", day_column = "Day", id_column = "ID",
                          time_point = 21)
  expect_true(all(ad$Day == 21))
})

test_that("prepare_dose_data uses last time point when time_point is NULL", {
  # Create multi-timepoint dose data
  df <- make_dose_response()
  df2 <- df
  df2$Day <- 14L
  df2$Volume <- df2$Volume * 0.8
  multi <- rbind(df2, df)

  ad <- prepare_dose_data(multi, dose_column = "Dose", treatment_column = "Treatment",
                          volume_column = "Volume", day_column = "Day", id_column = "ID",
                          time_point = NULL)
  # Each mouse should have one row at day 21 (the latest day)
  expect_true(all(ad$Day == 21))
})

test_that("prepare_dose_data errors on non-existent time point", {
  df <- make_dose_response()
  expect_error(
    prepare_dose_data(df, dose_column = "Dose", treatment_column = "Treatment",
                      volume_column = "Volume", day_column = "Day", id_column = "ID",
                      time_point = 999),
    "No data found"
  )
})

# ===========================================================================
# 4. repeated_measures_anova
# ===========================================================================

test_that("repeated_measures_anova returns an me_result", {
  skip_if_not_installed("lmerTest")

  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    repeated_measures_anova(df, time_column = "Day", volume_column = "Volume",
                            treatment_column = "Treatment", id_column = "ID",
                            transform = "log")
  ))

  expect_s3_class(res, "me_result")
})

test_that("repeated_measures_anova results contain anova_table", {
  skip_if_not_installed("lmerTest")

  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    repeated_measures_anova(df, time_column = "Day", volume_column = "Volume",
                            treatment_column = "Treatment", id_column = "ID",
                            transform = "log")
  ))

  expect_true("anova_table" %in% names(res$results))
  expect_s3_class(res$results$anova_table, "data.frame")
})

test_that("repeated_measures_anova detects significant interaction (make_tg_simple)", {
  skip_if_not_installed("lmerTest")

  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    repeated_measures_anova(df, time_column = "Day", volume_column = "Volume",
                            treatment_column = "Treatment", id_column = "ID",
                            transform = "log")
  ))

  anova_tbl <- res$results$anova_table
  # The interaction row (Treatment:Day) should have a very small p-value
  p_col <- grep("Pr|p.value|p\\.value", colnames(anova_tbl), value = TRUE, ignore.case = TRUE)
  if (length(p_col) > 0) {
    interaction_row <- grep(":", rownames(anova_tbl))
    if (length(interaction_row) > 0) {
      expect_true(anova_tbl[interaction_row[1], p_col[1]] < 0.05)
    }
  }
})

test_that("repeated_measures_anova respects transform argument", {
  skip_if_not_installed("lmerTest")

  df <- make_tg_simple()
  res_log <- suppressWarnings(suppressMessages(
    repeated_measures_anova(df, transform = "log")
  ))
  res_sqrt <- suppressWarnings(suppressMessages(
    repeated_measures_anova(df, transform = "sqrt")
  ))

  expect_equal(res_log$results$transform, "log")
  expect_equal(res_sqrt$results$transform, "sqrt")
})

test_that("repeated_measures_anova errors on missing columns", {
  skip_if_not_installed("lmerTest")

  df <- make_tg_simple()
  expect_error(
    repeated_measures_anova(df, volume_column = "NonExistent"),
    "Missing columns"
  )
})
