# =============================================================================
# Tests for toxicity & efficacy-toxicity analysis functions
# =============================================================================

# --- Helper: create test data with weight ---
# make_weight_data() moved to helper-fixtures.R (CODE_REVIEW.md K.8) so the
# toxicity fixture is shared and maintained in one place. It is deliberately
# NOT merged with make_bw_simple(): this one carries Volume and uses four
# timepoints, so they are different fixtures, not duplicates.


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
  expect_true(all(c("Treatment", "TGI", "Mean_Pct_Weight_Loss", "TWM") %in%
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


