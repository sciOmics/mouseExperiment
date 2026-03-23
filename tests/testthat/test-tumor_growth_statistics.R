# =============================================================================
# Tests for tumor_growth_statistics()
#
# Uses make_tg_simple() — two-group dataset with a very high SNR (≈42) so
# every assertion about significance is stable across platforms.
# =============================================================================

# ---------------------------------------------------------------------------
# LME4 path
# ---------------------------------------------------------------------------
test_that("lme4: returns a list with all expected top-level fields", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      model_type       = "lme4",
      transform        = "log",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))

  required <- c("model", "anova", "pairwise_comparisons",
                "treatment_effects", "growth_rates", "summary")
  expect_true(all(required %in% names(res)),
              info = paste("Missing:", paste(setdiff(required, names(res)), collapse = ", ")))
})

test_that("lme4: ANOVA is a data frame with at least one significant row (p < 0.001)", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      model_type       = "lme4",
      transform        = "log",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))

  anova_df <- as.data.frame(res$anova)
  expect_true(is.data.frame(anova_df))
  expect_true(nrow(anova_df) >= 1L)

  # At least one p-value must be significant — the Day:Treatment interaction
  # or the Treatment main effect (slope difference is massive at SNR≈42).
  p_cols <- grep("Pr|p.value|p_value", colnames(anova_df),
                 value = TRUE, ignore.case = TRUE)
  if (length(p_cols) > 0) {
    p_vals <- unlist(anova_df[, p_cols, drop = FALSE])
    p_vals <- suppressWarnings(as.numeric(p_vals))
    expect_true(any(p_vals < 0.001, na.rm = TRUE),
                info = paste("Smallest p-value:", min(p_vals, na.rm = TRUE)))
  }
})

test_that("lme4: pairwise comparison for TreatmentA > Control has correct direction", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      model_type       = "lme4",
      transform        = "log",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))

  pw <- as.data.frame(res$pairwise_comparisons)
  expect_true(nrow(pw) >= 1L)

  # The contrast name should mention TreatmentA
  contrast_col <- intersect(c("contrast", "Contrast", "comparison"), colnames(pw))[1]
  est_col      <- intersect(c("estimate", "Estimate"), colnames(pw))[1]
  p_col        <- intersect(c("p.value", "p_value", "adj.p.value"), colnames(pw))[1]

  if (!is.na(contrast_col) && !is.na(est_col)) {
    row_a <- grepl("TreatmentA", ignore.case = TRUE, pw[[contrast_col]])
    if (any(row_a)) {
      est <- as.numeric(pw[[est_col]][which(row_a)[1]])
      # TreatmentA grows faster → positive pairwise estimate vs Control
      expect_true(est > 0,
                  info = paste("Expected positive estimate for TreatmentA, got", est))
    }
  }

  if (!is.na(p_col)) {
    expect_true(any(as.numeric(pw[[p_col]]) < 0.05, na.rm = TRUE),
                info = "No significant pairwise comparison found")
  }
})

test_that("lme4: treatment_effects has one row per group and numeric estimates", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      model_type       = "lme4",
      transform        = "log",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))

  eff <- as.data.frame(res$treatment_effects)
  expect_true(is.data.frame(eff))
  expect_equal(nrow(eff), 2L, info = "Expected 2 rows (Control & TreatmentA)")

  num_cols <- sapply(eff, is.numeric)
  expect_true(any(num_cols), info = "No numeric columns in treatment_effects")
})

test_that("lme4: growth_rates contains one row per subject", {
  df  <- make_tg_simple()
  n_mice <- length(unique(df$ID))
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      model_type       = "lme4",
      transform        = "log",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))

  gr <- as.data.frame(res$growth_rates)
  expect_true(is.data.frame(gr))
  expect_equal(nrow(gr), n_mice)
})

test_that("lme4: mean growth rate of TreatmentA > mean growth rate of Control", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      model_type       = "lme4",
      transform        = "log",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))

  gr <- as.data.frame(res$growth_rates)
  tx_col   <- intersect(c("Treatment", "treatment", "Group"), colnames(gr))[1]
  rate_col <- intersect(c("Growth_Rate", "rate", "slope", "Slope"), colnames(gr))[1]

  if (!is.na(tx_col) && !is.na(rate_col)) {
    mean_ctrl <- mean(gr[[rate_col]][gr[[tx_col]] == "Control"],    na.rm = TRUE)
    mean_tx   <- mean(gr[[rate_col]][gr[[tx_col]] == "TreatmentA"], na.rm = TRUE)
    expect_true(mean_tx > mean_ctrl,
                info = paste("Control mean rate:", mean_ctrl,
                             "TreatmentA mean rate:", mean_tx))
  }
})

# ---------------------------------------------------------------------------
# AUC model path
# ---------------------------------------------------------------------------
test_that("auc model: returns auc_analysis and pairwise_comparisons fields", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      model_type       = "auc",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))

  expect_true(!isTRUE(res$error),
              info = paste("AUC model errored:", res$message))

  # Either auc_analysis or treatment_effects must be present
  has_auc  <- "auc_analysis"      %in% names(res)
  has_eff  <- "treatment_effects" %in% names(res)
  expect_true(has_auc || has_eff,
              info = "Neither auc_analysis nor treatment_effects returned by AUC model")
})

test_that("auc model: mean AUC of TreatmentA > Control", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      cage_column      = "Cage",
      id_column        = "ID",
      model_type       = "auc",
      reference_group  = "Control",
      plots            = FALSE,
      verbose          = FALSE
    )
  ))

  skip_if(isTRUE(res$error), "AUC model returned an error")

  # auc_analysis$individual has per-mouse AUCs
  if (!is.null(res$auc_analysis) && !is.null(res$auc_analysis$individual)) {
    ind <- as.data.frame(res$auc_analysis$individual)
    grp_col <- intersect(c("Group", "Treatment", "treatment"), colnames(ind))[1]
    auc_col <- intersect(c("AUC", "auc"), colnames(ind))[1]
    if (!is.na(grp_col) && !is.na(auc_col)) {
      mu_ctrl <- mean(ind[[auc_col]][ind[[grp_col]] == "Control"],    na.rm = TRUE)
      mu_tx   <- mean(ind[[auc_col]][ind[[grp_col]] == "TreatmentA"], na.rm = TRUE)
      expect_true(mu_tx > mu_ctrl,
                  info = paste("Control mean AUC:", mu_ctrl,
                               "TreatmentA mean AUC:", mu_tx))
    }
  }
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
test_that("errors when required columns are missing", {
  df <- data.frame(ID = "M1", Treatment = "Control", Day = 0, Volume = 200)
  expect_error(
    tumor_growth_statistics(df, cage_column = "Cage", plots = FALSE),
    regexp = "Cage|column"
  )
})
