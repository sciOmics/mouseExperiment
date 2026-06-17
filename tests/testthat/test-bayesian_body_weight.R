# =============================================================================
# Tests for bayesian_body_weight()
#
# Uses make_bw_simple() (2 groups × 5 mice × 5 time-points).
# All tests are skipped when brms is not installed.
# 2 chains × 500 iterations keeps each run under 90 s on CI while still
# exercising every code path through the return list.
# =============================================================================

skip_bayes_bw <- function() {
  skip_if_not_installed("brms")
}

local({
  .cached_result <- NULL
  get_bw_result <- function() {
    if (is.null(.cached_result)) {
      df <- make_bw_simple()
      .cached_result <<- suppressWarnings(suppressMessages(
        bayesian_body_weight(
          df,
          weight_column                = "Weight",
          time_column                  = "Day",
          treatment_column             = "Treatment",
          id_column                    = "ID",
          cage_column                  = "Cage",
          transform                    = "none",
          random_effects_specification = "intercept_only",
          reference_group              = "Control",
          prior_strength               = "weakly_informative",
          n_chains                     = 2L,
          n_iter                       = 500L,
          seed                         = 42L,
          return_model                 = TRUE,
          plots                        = FALSE,
          verbose                      = FALSE
        )
      ))
    }
    .cached_result
  }

  # ── Return structure ────────────────────────────────────────────────────────
  test_that("bayesian_body_weight: returns list with all expected fields", {
    skip_bayes_bw()
    res <- get_bw_result()
    required <- c(
      "model", "model_type_used", "transform_used", "summary",
      "posterior_summary", "treatment_effects", "pairwise_comparisons",
      "mcmc_diagnostics", "weight_loss_summary", "weight_data"
    )
    expect_true(
      all(required %in% names(res)),
      info = paste("Missing:", paste(setdiff(required, names(res)), collapse = ", "))
    )
  })

  test_that("bayesian_body_weight: model_type_used is 'bayes_bw'", {
    skip_bayes_bw()
    expect_equal(get_bw_result()$model_type_used, "bayes_bw")
  })

  test_that("bayesian_body_weight: transform_used matches requested transform", {
    skip_bayes_bw()
    expect_equal(get_bw_result()$transform_used, "none")
  })

  test_that("bayesian_body_weight: model is a brmsfit when return_model = TRUE", {
    skip_bayes_bw()
    expect_s3_class(get_bw_result()$model, "brmsfit")
  })

  # ── treatment_effects ───────────────────────────────────────────────────────
  test_that("bayesian_body_weight: treatment_effects has required columns", {
    skip_bayes_bw()
    te <- get_bw_result()$treatment_effects
    required_cols <- c("Group", "Adjusted_Mean", "Lower_CrI", "Upper_CrI", "Note")
    expect_s3_class(te, "data.frame")
    expect_true(
      all(required_cols %in% colnames(te)),
      info = paste("Missing:", paste(setdiff(required_cols, colnames(te)), collapse = ", "))
    )
  })

  test_that("bayesian_body_weight: reference group is labelled correctly", {
    skip_bayes_bw()
    te  <- get_bw_result()$treatment_effects
    ref <- te[te$Note == "Reference group", ]
    expect_equal(nrow(ref), 1L)
    expect_equal(ref$Group, "Control")
  })

  # ── mcmc_diagnostics ───────────────────────────────────────────────────────
  test_that("bayesian_body_weight: convergence per current Stan recommendations", {
    skip_bayes_bw()
    diag <- get_bw_result()$mcmc_diagnostics
    expect_s3_class(diag, "data.frame")
    expect_true("Rhat" %in% colnames(diag))
    expect_true("Converged" %in% colnames(diag))
    expect_type(diag$Converged, "logical")
    expect_true(all(diag$Rhat <= 1.01, na.rm = TRUE),
                info = "Rhat above 1.01 (Vehtari 2021 threshold)")
    expect_true(all(diag$Bulk_ESS >= 400, na.rm = TRUE),
                info = "Bulk_ESS below 400 (Stan recommendation)")
    expect_true(all(diag$Tail_ESS >= 400, na.rm = TRUE),
                info = "Tail_ESS below 400 (Stan recommendation)")
  })

  # ── posterior_summary ──────────────────────────────────────────────────────
  test_that("bayesian_body_weight: posterior_summary is a non-empty data frame with Rhat", {
    skip_bayes_bw()
    ps <- get_bw_result()$posterior_summary
    expect_s3_class(ps, "data.frame")
    expect_gt(nrow(ps), 0L)
    expect_true("Rhat" %in% colnames(ps))
  })

  # ── weight_loss_summary ────────────────────────────────────────────────────
  test_that("bayesian_body_weight: weight_loss_summary has required columns", {
    skip_bayes_bw()
    wls <- get_bw_result()$weight_loss_summary
    expect_s3_class(wls, "data.frame")
    required_cols <- c("Group", "Pct_Change", "Lower_CrI", "Upper_CrI",
                       "Day_Baseline", "Day_End")
    expect_true(
      all(required_cols %in% colnames(wls)),
      info = paste("Missing:", paste(setdiff(required_cols, colnames(wls)), collapse = ", "))
    )
  })

  test_that("bayesian_body_weight: weight_loss_summary has one row per group", {
    skip_bayes_bw()
    wls <- get_bw_result()$weight_loss_summary
    df  <- make_bw_simple()
    expect_equal(nrow(wls), length(unique(df$Treatment)))
  })

  test_that("bayesian_body_weight: TreatmentA has negative Pct_Change (weight loss)", {
    skip_bayes_bw()
    wls   <- get_bw_result()$weight_loss_summary
    tx_row <- wls[wls$Group == "TreatmentA", ]
    expect_equal(nrow(tx_row), 1L)
    expect_lt(tx_row$Pct_Change, 0)
  })

  # ── weight_data passthrough ────────────────────────────────────────────────
  test_that("bayesian_body_weight: weight_data mirrors input row count", {
    skip_bayes_bw()
    wd <- get_bw_result()$weight_data
    df <- make_bw_simple()
    expect_s3_class(wd, "data.frame")
    expect_equal(nrow(wd), nrow(df))
  })

  # ── summary metadata ───────────────────────────────────────────────────────
  test_that("bayesian_body_weight: summary contains expected metadata", {
    skip_bayes_bw()
    s <- get_bw_result()$summary
    expect_equal(s$data_description$reference_group, "Control")
    expect_equal(s$model_specification$transform, "none")
    expect_equal(s$model_specification$prior_strength, "weakly_informative")
  })

  # ── Input validation ───────────────────────────────────────────────────────
  test_that("bayesian_body_weight: errors on missing required column", {
    skip_bayes_bw()
    df <- make_bw_simple()
    df$Weight <- NULL
    expect_error(
      bayesian_body_weight(df, weight_column = "Weight"),
      "Missing required columns"
    )
  })

  test_that("bayesian_body_weight: errors on invalid reference group", {
    skip_bayes_bw()
    expect_error(
      suppressMessages(bayesian_body_weight(
        make_bw_simple(),
        reference_group = "NoSuchGroup",
        n_chains = 2L, n_iter = 200L, plots = FALSE, verbose = FALSE
      )),
      "not found"
    )
  })
})

# ── return_model = FALSE ────────────────────────────────────────────────────

test_that("bayesian_body_weight: model is NULL when return_model = FALSE", {
  skip_bayes_bw()
  res <- suppressWarnings(suppressMessages(
    bayesian_body_weight(
      make_bw_simple(),
      reference_group = "Control",
      n_chains = 2L, n_iter = 500L,
      return_model = FALSE, plots = FALSE, verbose = FALSE, seed = 1L
    )
  ))
  expect_null(res$model)
})

# ── plots = TRUE smoke test ─────────────────────────────────────────────────

test_that("bayesian_body_weight: weight_trajectory_plot is a ggplot when plots = TRUE", {
  skip_bayes_bw()
  skip_if_not_installed("bayesplot")
  res <- suppressWarnings(suppressMessages(
    bayesian_body_weight(
      make_bw_simple(),
      reference_group = "Control",
      n_chains = 2L, n_iter = 500L,
      return_model = FALSE, plots = TRUE, verbose = FALSE, seed = 2L
    )
  ))
  expect_s3_class(res$weight_trajectory_plot, "gg")
})
