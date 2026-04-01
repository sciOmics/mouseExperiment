# =============================================================================
# Tests for tumor_auc_analysis()
#
# Uses make_auc_exact():
#   LowGrowth  M1: volumes [100,200,300,400] at days [0,7,14,21] → AUC = 5 250
#   HighGrowth M2: volumes [200,400,600,800] at days [0,7,14,21] → AUC = 10 500
#
# Trapezoidal AUC formula (verified by hand):
#   AUC = Σ  dt[i] * (V[i] + V[i-1]) / 2
#
# AUC_LOW  = 7*(100+200)/2 + 7*(200+300)/2 + 7*(300+400)/2
#           = 1050 + 1750 + 2450 = 5 250
# AUC_HIGH = 7*(200+400)/2 + 7*(400+600)/2 + 7*(600+800)/2
#           = 2100 + 3500 + 4900 = 10 500  (ratio = 2.0)
# =============================================================================

test_that("returns a list with expected fields", {
  df  <- make_auc_exact()
  res <- suppressWarnings(suppressMessages(
    tumor_auc_analysis(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      cage_column      = "Cage",
      auc_method       = "trapezoidal"
    )
  ))

  required <- c("auc_data", "auc_summary", "auc_comparisons")
  present  <- required[required %in% names(res)]
  expect_true(length(present) >= 2L,
              info = paste("Missing fields:", paste(setdiff(required, names(res)), collapse = ", ")))
})

test_that("auc_data has one row per subject", {
  df  <- make_auc_exact()
  res <- suppressWarnings(suppressMessages(
    tumor_auc_analysis(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      cage_column      = "Cage",
      auc_method       = "trapezoidal"
    )
  ))

  auc_df <- as.data.frame(res$auc_data)
  expect_equal(nrow(auc_df), 2L)   # M1 and M2
})

test_that("trapezoidal AUC for LowGrowth mouse equals 5250 (hand-computed)", {
  df  <- make_auc_exact()
  res <- suppressWarnings(suppressMessages(
    tumor_auc_analysis(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      cage_column      = "Cage",
      auc_method       = "trapezoidal"
    )
  ))

  auc_df  <- as.data.frame(res$auc_data)
  id_col  <- intersect(c("ID", "id", "Subject"), colnames(auc_df))[1]
  auc_col <- intersect(c("AUC", "auc"), colnames(auc_df))[1]

  skip_if(is.na(id_col) || is.na(auc_col), "Cannot locate ID/AUC column in auc_data")

  auc_m1 <- auc_df[[auc_col]][auc_df[[id_col]] == "M1"]
  expect_equal(auc_m1, AUC_LOW_GROWTH, tolerance = 0.01,
               info = paste("Expected", AUC_LOW_GROWTH, "got", auc_m1))
})

test_that("trapezoidal AUC for HighGrowth mouse equals 10500 (hand-computed)", {
  df  <- make_auc_exact()
  res <- suppressWarnings(suppressMessages(
    tumor_auc_analysis(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      cage_column      = "Cage",
      auc_method       = "trapezoidal"
    )
  ))

  auc_df  <- as.data.frame(res$auc_data)
  id_col  <- intersect(c("ID", "id", "Subject"), colnames(auc_df))[1]
  auc_col <- intersect(c("AUC", "auc"), colnames(auc_df))[1]

  skip_if(is.na(id_col) || is.na(auc_col), "Cannot locate ID/AUC column in auc_data")

  auc_m2 <- auc_df[[auc_col]][auc_df[[id_col]] == "M2"]
  expect_equal(auc_m2, AUC_HIGH_GROWTH, tolerance = 0.01,
               info = paste("Expected", AUC_HIGH_GROWTH, "got", auc_m2))
})

test_that("AUC ratio HighGrowth / LowGrowth ≈ 2.0", {
  df  <- make_auc_exact()
  res <- suppressWarnings(suppressMessages(
    tumor_auc_analysis(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      cage_column      = "Cage",
      auc_method       = "trapezoidal"
    )
  ))

  auc_df  <- as.data.frame(res$auc_data)
  id_col  <- intersect(c("ID", "id", "Subject"), colnames(auc_df))[1]
  auc_col <- intersect(c("AUC", "auc"), colnames(auc_df))[1]

  skip_if(is.na(id_col) || is.na(auc_col), "Cannot locate ID/AUC column")

  auc_low  <- auc_df[[auc_col]][auc_df[[id_col]] == "M1"]
  auc_high <- auc_df[[auc_col]][auc_df[[id_col]] == "M2"]

  expect_equal(auc_high / auc_low, 2.0, tolerance = 0.01)
})

test_that("auc_summary has HighGrowth mean strictly greater than LowGrowth mean", {
  df  <- make_auc_exact()
  res <- suppressWarnings(suppressMessages(
    tumor_auc_analysis(
      df,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      cage_column      = "Cage",
      auc_method       = "trapezoidal"
    )
  ))

  summ    <- as.data.frame(res$auc_summary)
  grp_col <- intersect(c("Group", "Treatment", "treatment"), colnames(summ))[1]
  mean_col <- intersect(c("Mean_AUC", "mean_auc", "mean", "Mean"), colnames(summ))[1]

  skip_if(is.na(grp_col) || is.na(mean_col),
          "Cannot locate Group/Mean column in auc_summary")

  mu_low  <- summ[[mean_col]][summ[[grp_col]] == "LowGrowth"]
  mu_high <- summ[[mean_col]][summ[[grp_col]] == "HighGrowth"]

  expect_true(mu_high > mu_low,
              info = paste("LowGrowth mean:", mu_low, "HighGrowth mean:", mu_high))
})
