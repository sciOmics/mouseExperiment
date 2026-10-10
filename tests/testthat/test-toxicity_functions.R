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
    volume_units     = "mm3",   # required with adjustment since v0.25.0 (R20.83)
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
    reference_group  = "Control",
    volume_units     = "mm3",
    boot_seed        = 1
  )

  expect_type(res, "list")
  win <- res$window_table
  expect_true(is.data.frame(win))
  expect_true(all(c("Treatment", "N_Animals", "TGI", "TGI_Lower", "TGI_Upper",
                    "Worst_Loss", "Worst_Loss_Lower", "Worst_Loss_Upper",
                    "N_Over_Threshold", "Tolerability") %in% names(win)))
  expect_identical(win$Treatment[1], "Control")          # reference first
  # DrugA is effective and less toxic; DrugB is neither.
  a <- win[win$Treatment == "DrugA", ]; b <- win[win$Treatment == "DrugB", ]
  expect_gt(a$TGI, b$TGI)
  expect_lt(a$Worst_Loss, b$Worst_Loss)
  expect_identical(a$Tolerability, "Tolerated")
  expect_identical(b$Tolerability, "Not tolerated")
  expect_equal(b$N_Over_Threshold, 4L)
})

test_that("therapeutic_window_metric tolerability threshold is in percent", {
  df <- make_weight_data()
  tw <- function(...) therapeutic_window_metric(
    df, reference_group = "Control", volume_units = "mm3", boot_seed = 1, ...)
  expect_error(tw(tolerability_threshold = 0.2), "percent")
  # At 10 %, DrugA's worst loss (about 15 %) is no longer tolerated.
  win <- tw(tolerability_threshold = 10)$window_table
  expect_identical(win$Tolerability[win$Treatment == "DrugA"], "Not tolerated")
  expect_identical(win$Tolerability[win$Treatment == "Control"], "Tolerated")
})


