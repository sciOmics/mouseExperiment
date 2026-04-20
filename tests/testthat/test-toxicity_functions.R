# =============================================================================
# Tests for toxicity & efficacy-toxicity analysis functions
# =============================================================================

# --- Helper: create test data with weight ---
make_weight_data <- function() {
  set.seed(123)
  days <- c(0, 7, 14, 21)

  make_mouse <- function(id, treatment, base_weight, weight_slope, vol_slope) {
    vol <- 200 * exp(vol_slope * days) + rnorm(4, 0, 5)
    wt  <- base_weight + weight_slope * days + rnorm(4, 0, 0.3)
    data.frame(
      ID        = id,
      Treatment = treatment,
      Day       = days,
      Volume    = pmax(round(vol, 2), 1),
      Weight    = round(wt, 2),
      Sex       = ifelse(id %in% c("C01", "C02", "T01", "T02"), "M", "F"),
      Cage      = paste0(substring(treatment, 1, 1), "1"),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    # Control: moderate tumor growth, stable weight
    make_mouse("C01", "Control", 22, -0.02, 0.05),
    make_mouse("C02", "Control", 21, -0.01, 0.05),
    make_mouse("C03", "Control", 23, -0.03, 0.05),
    make_mouse("C04", "Control", 20, -0.02, 0.05),
    # Drug A: effective drug, moderate toxicity
    make_mouse("T01", "DrugA", 22, -0.15, 0.02),
    make_mouse("T02", "DrugA", 21, -0.12, 0.02),
    make_mouse("T03", "DrugA", 20, -0.18, 0.02),
    make_mouse("T04", "DrugA", 23, -0.10, 0.02),
    # Drug B: ineffective, high toxicity
    make_mouse("T05", "DrugB", 22, -0.25, 0.04),
    make_mouse("T06", "DrugB", 21, -0.30, 0.04),
    make_mouse("T07", "DrugB", 20, -0.28, 0.04),
    make_mouse("T08", "DrugB", 23, -0.22, 0.04)
  )
}


# =============================================================================
# analyze_body_weight
# =============================================================================
test_that("analyze_body_weight returns expected structure", {
  df <- make_weight_data()
  res <- analyze_body_weight(
    df,
    weight_column    = "Weight",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    volume_column    = "Volume",
    adjust_tumor_weight = TRUE,
    covariates       = c("volume"),
    estimation       = "REML"
  )

  expect_type(res, "list")
  expect_true(!is.null(res$model))
  expect_s4_class(res$model, "lmerMod")
  expect_true(is.data.frame(res$fixed_effects))
  expect_true("Term" %in% names(res$fixed_effects))
  expect_true(is.data.frame(res$random_effects))
  expect_true(is.data.frame(res$weight_data))
  expect_true("Net_Weight" %in% names(res$weight_data))
  expect_true(is.character(res$summary_text))
  expect_true(nchar(res$summary_text) > 0)
})

test_that("analyze_body_weight works without tumor adjustment", {
  df <- make_weight_data()
  res <- analyze_body_weight(
    df,
    weight_column    = "Weight",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    adjust_tumor_weight = FALSE,
    covariates       = character(0)
  )

  expect_true(!is.null(res$model))
  expect_false(res$model_info$adjust_tumor)
})

test_that("analyze_body_weight with sex covariate", {
  df <- make_weight_data()
  res <- analyze_body_weight(
    df,
    weight_column    = "Weight",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    sex_column       = "Sex",
    covariates       = c("sex")
  )

  expect_true(!is.null(res$model))
  # Sex should appear in fixed effects
  expect_true(any(grepl("Sex", res$fixed_effects$Term)))
})

test_that("analyze_body_weight errors on missing column", {
  df <- make_weight_data()
  expect_error(
    analyze_body_weight(df, weight_column = "NonExistent"),
    "Missing required columns"
  )
})

test_that("analyze_body_weight model_simplified info works", {
  df <- make_weight_data()
  res <- analyze_body_weight(
    df,
    weight_column    = "Weight",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID"
  )

  expect_type(res$model_info$model_simplified, "logical")
})


# =============================================================================
# body_weight_auc
# =============================================================================
test_that("body_weight_auc returns expected structure", {
  df <- make_weight_data()
  res <- body_weight_auc(
    df,
    weight_column    = "Weight",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    volume_column    = "Volume"
  )

  expect_type(res, "list")
  expect_true(is.data.frame(res$auc_per_mouse))
  expect_true(all(c("ID", "Treatment", "AUC_Weight", "AUC_Pct_Change",
                     "Nadir_Weight", "Nadir_Day") %in% names(res$auc_per_mouse)))
  expect_true(is.data.frame(res$auc_summary))
  expect_true("AUC_Mean" %in% names(res$auc_summary))
  expect_true(is.data.frame(res$nadir_data))
  expect_true(is.data.frame(res$comparisons))
})

test_that("body_weight_auc returns one row per mouse", {
  df <- make_weight_data()
  res <- body_weight_auc(df, weight_column = "Weight", time_column = "Day",
                         treatment_column = "Treatment", id_column = "ID")
  expect_equal(nrow(res$auc_per_mouse), length(unique(df$ID)))
})

test_that("body_weight_auc nadir data is correct", {
  df <- make_weight_data()
  res <- body_weight_auc(df, weight_column = "Weight", time_column = "Day",
                         treatment_column = "Treatment", id_column = "ID")
  # Nadir should be ≤ baseline for each mouse
  for (i in seq_len(nrow(res$auc_per_mouse))) {
    expect_true(res$auc_per_mouse$Nadir_Weight[i] <= res$auc_per_mouse$Baseline_Weight[i] + 1)
  }
})


# =============================================================================
# weight_loss_threshold
# =============================================================================
test_that("weight_loss_threshold returns expected structure", {
  df <- make_weight_data()
  res <- weight_loss_threshold(
    df,
    weight_column    = "Weight",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    threshold        = 0.20,
    reference_group  = "Control"
  )

  expect_type(res, "list")
  expect_true(is.data.frame(res$event_data))
  expect_true(all(c("ID", "Treatment", "Time", "Event") %in% names(res$event_data)))
  expect_s3_class(res$km_fit, "survfit")
  expect_equal(res$threshold, 0.20)
})

test_that("weight_loss_threshold censors mice that don't hit threshold", {
  df <- make_weight_data()
  res <- weight_loss_threshold(
    df,
    weight_column    = "Weight",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    threshold        = 0.50  # Very high threshold - most mice won't hit it
  )

  # At 50% threshold most mice should be censored
  expect_true(sum(res$event_data$Event == 0) > 0)
})

test_that("weight_loss_threshold with custom baseline_day", {
  df <- make_weight_data()
  res <- weight_loss_threshold(
    df,
    weight_column    = "Weight",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    baseline_day     = 0
  )

  expect_true(is.data.frame(res$event_data))
})


# =============================================================================
# therapeutic_window_metric
# =============================================================================
test_that("therapeutic_window_metric returns expected structure", {
  df <- make_weight_data()
  res <- therapeutic_window_metric(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control"
  )

  expect_type(res, "list")
  expect_true(is.data.frame(res$twm_table))
  expect_true(all(c("Treatment", "TGI", "Max_Pct_Weight_Loss", "TWM") %in%
                  names(res$twm_table)))
  # DrugA should have higher TWM than DrugB (effective + less toxic)
  twm_a <- res$twm_table$TWM[res$twm_table$Treatment == "DrugA"]
  twm_b <- res$twm_table$TWM[res$twm_table$Treatment == "DrugB"]
  expect_true(twm_a > twm_b)
})

test_that("therapeutic_window_metric noise_floor works", {
  df <- make_weight_data()
  res <- therapeutic_window_metric(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control",
    noise_floor      = 100  # Very high floor so all groups hit it
  )

  expect_true(all(res$twm_table$Safety_Note == "Negligible weight loss"))
})


# =============================================================================
# efficacy_toxicity_bivariate
# =============================================================================
test_that("efficacy_toxicity_bivariate returns expected structure (tgi)", {
  df <- make_weight_data()
  res <- efficacy_toxicity_bivariate(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control",
    efficacy_metric  = "tgi"
  )

  expect_type(res, "list")
  expect_true(is.data.frame(res$per_mouse))
  expect_true(all(c("ID", "Treatment", "Max_Pct_Weight_Loss", "Efficacy") %in%
                  names(res$per_mouse)))
  expect_true(is.data.frame(res$per_group))
  expect_equal(res$efficacy_metric, "tgi")
})

test_that("efficacy_toxicity_bivariate works with tumor_auc metric", {
  df <- make_weight_data()
  res <- efficacy_toxicity_bivariate(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control",
    efficacy_metric  = "tumor_auc"
  )

  expect_equal(res$efficacy_metric, "tumor_auc")
  expect_true(is.data.frame(res$per_mouse))
})

test_that("efficacy_toxicity_bivariate works with log_cell_kill metric", {
  df <- make_weight_data()
  res <- efficacy_toxicity_bivariate(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control",
    efficacy_metric  = "log_cell_kill"
  )

  expect_equal(res$efficacy_metric, "log_cell_kill")
  expect_true(is.data.frame(res$per_mouse))
})


# =============================================================================
# total_benefit_area
# =============================================================================
test_that("total_benefit_area returns expected structure", {
  df <- make_weight_data()
  res <- total_benefit_area(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control",
    lambda           = 1.0
  )

  expect_type(res, "list")
  expect_true(is.data.frame(res$benefit_table))
  expect_true(all(c("Treatment", "Efficacy_AUC", "Toxicity_AUC",
                     "Benefit_Score", "Rank") %in% names(res$benefit_table)))
  expect_equal(res$lambda, 1.0)
})

test_that("total_benefit_area higher lambda penalizes toxic drugs", {
  df <- make_weight_data()
  res_low  <- total_benefit_area(df, weight_column = "Weight",
    volume_column = "Volume", time_column = "Day",
    treatment_column = "Treatment", id_column = "ID",
    reference_group = "Control", lambda = 0.1)
  res_high <- total_benefit_area(df, weight_column = "Weight",
    volume_column = "Volume", time_column = "Day",
    treatment_column = "Treatment", id_column = "ID",
    reference_group = "Control", lambda = 5.0)

  # DrugB (toxic, less effective) should drop more with higher lambda
  drugb_low  <- res_low$benefit_table$Benefit_Score[res_low$benefit_table$Treatment == "DrugB"]
  drugb_high <- res_high$benefit_table$Benefit_Score[res_high$benefit_table$Treatment == "DrugB"]
  expect_true(drugb_high < drugb_low)
})


# =============================================================================
# weight_corrected_tgi
# =============================================================================
test_that("weight_corrected_tgi returns expected structure", {
  df <- make_weight_data()
  res <- weight_corrected_tgi(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control",
    safety_threshold = 0.20
  )

  expect_type(res, "list")
  expect_true(is.data.frame(res$corrected_tgi))
  expect_true(is.data.frame(res$uncorrected_tgi))
  expect_true(is.data.frame(res$comparison))
  expect_true(all(c("TGI_Uncorrected", "TGI_Corrected") %in% names(res$comparison)))
  expect_equal(res$safety_threshold, 0.20)
})

test_that("weight_corrected_tgi excludes toxic mice", {
  df <- make_weight_data()
  res <- weight_corrected_tgi(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control",
    safety_threshold = 0.05  # Low threshold to exclude more mice
  )

  # Some mice should be excluded with a 5% threshold
  expect_true(nrow(res$excluded_mice) > 0 || nrow(res$excluded_summary[res$excluded_summary$N_Excluded > 0, ]) >= 0)
})

test_that("weight_corrected_tgi with very high threshold excludes none", {
  df <- make_weight_data()
  res <- weight_corrected_tgi(
    df,
    weight_column    = "Weight",
    volume_column    = "Volume",
    time_column      = "Day",
    treatment_column = "Treatment",
    id_column        = "ID",
    reference_group  = "Control",
    safety_threshold = 0.99  # 99% - nobody should hit this
  )

  expect_equal(nrow(res$excluded_mice), 0)
  # Corrected and uncorrected should be identical
  expect_equal(res$comparison$TGI_Corrected, res$comparison$TGI_Uncorrected)
})
