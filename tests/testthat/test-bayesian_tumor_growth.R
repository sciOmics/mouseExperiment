# =============================================================================
# Tests for bayesian_tumor_growth()
#
# Uses make_tg_simple() (2 groups × 6 mice × 4 time-points).
# All tests are skipped when brms / bayesplot are not installed.
# We use 2 chains × 500 iterations with adapt_delta defaults so each run
# completes in < 60 s on CI and still exercises the full code path.
# =============================================================================

skip_bayes <- function() {
  skip_if_not_installed("brms")
  skip_if_not_installed("bayesplot")
}

# Fit once and cache in local env to avoid redundant compilation per test_that
local({
  .cached_result <- NULL
  get_bayes_result <- function() {
    if (is.null(.cached_result)) {
      df <- make_tg_simple()
      .cached_result <<- suppressWarnings(suppressMessages(
        bayesian_tumor_growth(
          df,
          time_column      = "Day",
          volume_column    = "Volume",
          treatment_column = "Treatment",
          id_column        = "ID",
          transform        = "log",
          random_effects_specification = "intercept_only",
          reference_group  = "Control",
          n_chains         = 2L,
          n_iter           = 500L,
          seed             = 42L,
          return_model     = TRUE,
          plots            = FALSE,
          verbose          = FALSE
        )
      ))
    }
    .cached_result
  }

  test_that("bayesian_tumor_growth: returns list with all expected top-level fields", {
    skip_bayes()
    res <- get_bayes_result()
    required <- c(
      "model", "model_type_used", "transform_used", "summary",
      "posterior_summary", "treatment_effects", "pairwise_comparisons",
      "mcmc_diagnostics", "growth_rates", "data_summary"
    )
    expect_true(
      all(required %in% names(res)),
      info = paste("Missing:", paste(setdiff(required, names(res)), collapse = ", "))
    )
  })

  test_that("bayesian_tumor_growth: model_type_used is 'bayes_tg'", {
    skip_bayes()
    res <- get_bayes_result()
    expect_equal(res$model_type_used, "bayes_tg")
  })

  test_that("bayesian_tumor_growth: treatment_effects uses CrI column naming", {
    skip_bayes()
    res <- get_bayes_result()
    required_cols <- c("Group", "Adjusted_Mean", "Lower_CrI", "Upper_CrI")
    te <- res$treatment_effects
    expect_s3_class(te, "data.frame")
    expect_true(
      all(required_cols %in% colnames(te)),
      info = paste("Missing:", paste(setdiff(required_cols, colnames(te)), collapse = ", "))
    )
  })

  test_that("bayesian_tumor_growth: pairwise_comparisons is a non-empty data frame", {
    skip_bayes()
    res <- get_bayes_result()
    expect_s3_class(res$pairwise_comparisons, "data.frame")
    expect_gt(nrow(res$pairwise_comparisons), 0)
  })

  test_that("bayesian_tumor_growth: convergence per current Stan recommendations", {
    skip_bayes()
    res <- get_bayes_result()
    diag <- res$mcmc_diagnostics
    expect_s3_class(diag, "data.frame")
    expect_true("Rhat" %in% colnames(diag),
                info = "mcmc_diagnostics must contain an Rhat column")
    # Vehtari et al. 2021 / current Stan recommendation: Rhat <= 1.01
    expect_true(all(diag$Rhat <= 1.01, na.rm = TRUE),
                info = paste("Rhat above 1.01:",
                             paste(round(diag$Rhat[diag$Rhat > 1.01], 4),
                                   collapse = ", ")))
    # Bulk / Tail ESS minimum thresholds (Stan default suggestion ~400 each)
    expect_true("Bulk_ESS" %in% colnames(diag))
    expect_true("Tail_ESS" %in% colnames(diag))
    expect_true(all(diag$Bulk_ESS >= 400, na.rm = TRUE),
                info = paste("Bulk_ESS below 400:",
                             paste(diag$Bulk_ESS[diag$Bulk_ESS < 400],
                                   collapse = ", ")))
    expect_true(all(diag$Tail_ESS >= 400, na.rm = TRUE),
                info = paste("Tail_ESS below 400:",
                             paste(diag$Tail_ESS[diag$Tail_ESS < 400],
                                   collapse = ", ")))
  })

  test_that("bayesian_tumor_growth: mcmc_diagnostics has Converged column", {
    skip_bayes()
    res <- get_bayes_result()
    expect_true("Converged" %in% colnames(res$mcmc_diagnostics))
  })

  test_that("bayesian_tumor_growth: posterior_summary is a data frame with at least one row", {
    skip_bayes()
    res <- get_bayes_result()
    expect_s3_class(res$posterior_summary, "data.frame")
    expect_gt(nrow(res$posterior_summary), 0)
  })

  test_that("bayesian_tumor_growth: growth_rates has same structure as lme4 path", {
    skip_bayes()
    res <- get_bayes_result()
    expect_s3_class(res$growth_rates, "data.frame")
    expect_true("growth_rate" %in% colnames(res$growth_rates))
  })

  test_that("bayesian_tumor_growth: model is a brmsfit object when return_model=TRUE", {
    skip_bayes()
    res <- get_bayes_result()
    expect_s3_class(res$model, "brmsfit")
  })

  test_that("bayesian_tumor_growth: model is NULL when return_model=FALSE", {
    skip_bayes()
    df <- make_tg_simple()
    res <- suppressWarnings(suppressMessages(
      bayesian_tumor_growth(
        df,
        time_column      = "Day",
        volume_column    = "Volume",
        treatment_column = "Treatment",
        id_column        = "ID",
        n_chains         = 2L,
        n_iter           = 500L,
        seed             = 42L,
        return_model     = FALSE,
        plots            = FALSE
      )
    ))
    expect_null(res$model)
  })

  test_that("bayesian_tumor_growth: stops with informative error when brms is unavailable", {
    skip_if(requireNamespace("brms", quietly = TRUE),
            "brms is installed — cannot test missing-package path")
    df <- make_tg_simple()
    expect_error(
      bayesian_tumor_growth(df),
      "brms"
    )
  })
})
