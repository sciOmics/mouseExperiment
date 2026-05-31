# =============================================================================
# Tests for bayesian_dose_response()
#
# Uses make_dose_response() fixture: Control (dose=0, mean vol≈500) +
# DrugX at doses 1, 5, 25 (mean vols ≈380, 220, 90), 4 mice/dose.
# Ground-truth: EC50 ≈ 4–5, Emax ≈ 0.88, Hill ≈ 1.5–2.
# All tests are skipped when brms is not installed.
# 2 chains × 600 iterations keeps each run under 90 s on CI.
# =============================================================================

skip_bayes_dr <- function() {
  skip_if_not_installed("brms")
}

local({
  .cached_result <- NULL
  get_dr_result <- function() {
    if (is.null(.cached_result)) {
      df <- make_dose_response()
      .cached_result <<- suppressWarnings(suppressMessages(
        bayesian_dose_response(
          df,
          dose_column      = "Dose",
          volume_column    = "Volume",
          treatment_column = "Treatment",
          id_column        = "ID",
          day_column       = "Day",
          endpoint_day     = 21L,
          reference_group  = "Control",
          prior_strength   = "skeptical",
          n_chains         = 2L,
          n_iter           = 600L,
          seed             = 42L,
          return_model     = TRUE,
          plots            = FALSE,
          verbose          = FALSE
        )
      ))
    }
    .cached_result
  }

  # ── Return structure ────────────────────────────────────────────────────────
  test_that("bayesian_dose_response: returns list with all expected fields", {
    skip_bayes_dr()
    res <- get_dr_result()
    required <- c(
      "model", "model_type_used", "summary", "posterior_summary",
      "dr_parameters", "mcmc_diagnostics", "dose_response_summary",
      "tgi_data", "control_mean_volume"
    )
    expect_true(
      all(required %in% names(res)),
      info = paste("Missing:", paste(setdiff(required, names(res)), collapse = ", "))
    )
  })

  test_that("bayesian_dose_response: model_type_used is 'bayes_dr'", {
    skip_bayes_dr()
    expect_equal(get_dr_result()$model_type_used, "bayes_dr")
  })

  test_that("bayesian_dose_response: model is brmsfit when return_model = TRUE", {
    skip_bayes_dr()
    expect_s3_class(get_dr_result()$model, "brmsfit")
  })

  # ── dr_parameters ──────────────────────────────────────────────────────────
  test_that("bayesian_dose_response: dr_parameters has EC50, Emax, Hill rows", {
    skip_bayes_dr()
    dp <- get_dr_result()$dr_parameters
    expect_s3_class(dp, "data.frame")
    expect_true(all(c("EC50", "Emax", "Hill") %in% dp$Parameter))
  })

  test_that("bayesian_dose_response: Emax posterior median is in (0, 1)", {
    skip_bayes_dr()
    dp   <- get_dr_result()$dr_parameters
    emax <- dp$Median[dp$Parameter == "Emax"]
    expect_gt(emax, 0)
    expect_lt(emax, 1)
  })

  test_that("bayesian_dose_response: EC50 is positive", {
    skip_bayes_dr()
    dp   <- get_dr_result()$dr_parameters
    ec50 <- dp$Median[dp$Parameter == "EC50"]
    expect_gt(ec50, 0)
  })

  test_that("bayesian_dose_response: Hill is positive", {
    skip_bayes_dr()
    dp   <- get_dr_result()$dr_parameters
    hill <- dp$Median[dp$Parameter == "Hill"]
    expect_gt(hill, 0)
  })

  test_that("bayesian_dose_response: EC50 credible interval brackets median", {
    skip_bayes_dr()
    dp   <- get_dr_result()$dr_parameters
    ec50 <- dp[dp$Parameter == "EC50", ]
    expect_lt(ec50$Lower_CrI, ec50$Median)
    expect_gt(ec50$Upper_CrI, ec50$Median)
  })

  # ── dose_response_summary ──────────────────────────────────────────────────
  test_that("bayesian_dose_response: dose_response_summary has required columns", {
    skip_bayes_dr()
    drs <- get_dr_result()$dose_response_summary
    required_cols <- c("Dose", "N", "Observed_TGI_Mean",
                       "Predicted_Median", "Predicted_Lower", "Predicted_Upper")
    expect_s3_class(drs, "data.frame")
    expect_true(all(required_cols %in% colnames(drs)))
  })

  test_that("bayesian_dose_response: dose_response_summary has one row per treated dose", {
    skip_bayes_dr()
    df  <- make_dose_response()
    drs <- get_dr_result()$dose_response_summary
    n_treated_doses <- length(unique(df$Dose[df$Dose > 0]))
    expect_equal(nrow(drs), n_treated_doses)
  })

  test_that("bayesian_dose_response: observed TGI increases with dose", {
    skip_bayes_dr()
    drs <- get_dr_result()$dose_response_summary
    drs <- drs[order(drs$Dose), ]
    expect_true(all(diff(drs$Observed_TGI_Mean) > 0))
  })

  # ── tgi_data and control_mean_volume ───────────────────────────────────────
  test_that("bayesian_dose_response: tgi_data excludes control and has TGI column", {
    skip_bayes_dr()
    td <- get_dr_result()$tgi_data
    expect_s3_class(td, "data.frame")
    expect_true("TGI" %in% colnames(td))
    expect_true(all(td$Dose > 0))
  })

  test_that("bayesian_dose_response: control_mean_volume is positive and finite", {
    skip_bayes_dr()
    cmv <- get_dr_result()$control_mean_volume
    expect_true(is.finite(cmv) && cmv > 0)
  })

  # ── mcmc_diagnostics ───────────────────────────────────────────────────────
  test_that("bayesian_dose_response: convergence per current Stan recommendations", {
    skip_bayes_dr()
    diag <- get_dr_result()$mcmc_diagnostics
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

  # ── summary metadata ───────────────────────────────────────────────────────
  test_that("bayesian_dose_response: summary has expected metadata", {
    skip_bayes_dr()
    s <- get_dr_result()$summary
    expect_equal(s$data_description$reference_group, "Control")
    expect_equal(s$data_description$endpoint_day, 21)
    expect_equal(s$model_specification$prior_strength, "skeptical")
  })

  # ── Input validation ───────────────────────────────────────────────────────
  test_that("bayesian_dose_response: errors on missing required column", {
    skip_bayes_dr()
    df <- make_dose_response()
    df$Dose <- NULL
    expect_error(
      bayesian_dose_response(df, dose_column = "Dose"),
      "Missing required columns"
    )
  })

  test_that("bayesian_dose_response: errors on invalid reference group", {
    skip_bayes_dr()
    expect_error(
      bayesian_dose_response(
        make_dose_response(),
        reference_group = "NoSuchGroup",
        n_chains = 2L, n_iter = 200L, plots = FALSE, verbose = FALSE
      ),
      "not found"
    )
  })
})

# ── return_model = FALSE ────────────────────────────────────────────────────

test_that("bayesian_dose_response: model is NULL when return_model = FALSE", {
  skip_bayes_dr()
  res <- suppressWarnings(suppressMessages(
    bayesian_dose_response(
      make_dose_response(),
      endpoint_day = 21L, reference_group = "Control",
      n_chains = 2L, n_iter = 600L, seed = 1L,
      return_model = FALSE, plots = FALSE, verbose = FALSE
    )
  ))
  expect_null(res$model)
})

# ── plots = TRUE smoke test ─────────────────────────────────────────────────

test_that("bayesian_dose_response: dose_response_curve_plot is a ggplot when plots = TRUE", {
  skip_bayes_dr()
  res <- suppressWarnings(suppressMessages(
    bayesian_dose_response(
      make_dose_response(),
      endpoint_day = 21L, reference_group = "Control",
      n_chains = 2L, n_iter = 600L, seed = 2L,
      return_model = FALSE, plots = TRUE, verbose = FALSE
    )
  ))
  expect_s3_class(res$dose_response_curve_plot, "gg")
})
