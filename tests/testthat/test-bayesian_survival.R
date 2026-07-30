# =============================================================================
# Tests for bayesian_survival()
#
# Uses make_surv_simple() (2 groups: AggressiveTx with 7 events,
# Control with 7 censored at day 60).
# All tests are skipped when brms / bayesplot are not installed.
# 2 chains × 500 iterations keeps each run under 90 s on CI while still
# exercising every code path through the return list.
# =============================================================================

# brms, bayesplot, gamm4 and mgcv are hard Imports as of v0.10.0, so this
# helper can no longer skip. Retained as a no-op because the call sites are
# numerous, and because a skip here is exactly what let bayesian_synergy()
# stay broken for five releases (CODE_REVIEW.md R3-L).
skip_bayes_surv <- function() invisible(NULL)

local({
  .cached_result <- NULL
  get_surv_result <- function() {
    if (is.null(.cached_result)) {
      df <- make_surv_simple()
      .cached_result <<- suppressWarnings(suppressMessages(
        bayesian_survival(
          df,
          time_column      = "Day",
          event_column     = "Survival_Censor",
          treatment_column = "Treatment",
          id_column        = "ID",
          cage_column      = "Cage",
          family           = "weibull",
          include_cage_effect  = FALSE,   # only 1 cage per group in fixture
          reference_group  = "Control",
          prior_strength   = "weakly_informative",
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

  # ── Return structure ────────────────────────────────────────────────────────
  test_that("bayesian_survival: returns list with all expected top-level fields", {
    skip_bayes_surv()
    res <- get_surv_result()
    required <- c(
      "model", "model_type_used", "family_used", "frailty_used",
      "summary", "treatment_effects", "posterior_summary",
      "mcmc_diagnostics", "survival_data"
    )
    expect_true(
      all(required %in% names(res)),
      info = paste("Missing:", paste(setdiff(required, names(res)), collapse = ", "))
    )
  })

  test_that("bayesian_survival: model_type_used is 'bayes_survival'", {
    skip_bayes_surv()
    expect_equal(get_surv_result()$model_type_used, "bayes_survival")
  })

  test_that("bayesian_survival: family_used matches requested family", {
    skip_bayes_surv()
    expect_equal(get_surv_result()$family_used, "weibull")
  })

  # ── treatment_effects table ─────────────────────────────────────────────────
  test_that("bayesian_survival: treatment_effects has required columns", {
    skip_bayes_surv()
    te <- get_surv_result()$treatment_effects
    required_cols <- c(
      "Group", "Time_Ratio", "Lower_CrI", "Upper_CrI",
      "HR", "Median_Survival", "Events", "Total", "Event_Rate", "Note"
    )
    expect_s3_class(te, "data.frame")
    expect_true(
      all(required_cols %in% colnames(te)),
      info = paste("Missing:", paste(setdiff(required_cols, colnames(te)), collapse = ", "))
    )
  })

  test_that("bayesian_survival: reference group has Time_Ratio = 1 and HR = 1", {
    skip_bayes_surv()
    te  <- get_surv_result()$treatment_effects
    ref <- te[te$Note == "Reference group", ]
    expect_equal(nrow(ref), 1L)
    expect_equal(ref$Time_Ratio, 1)
    expect_equal(ref$HR, 1)
  })

  test_that("bayesian_survival: non-reference group TR credible interval excludes 1", {
    skip_bayes_surv()
    te      <- get_surv_result()$treatment_effects
    non_ref <- te[te$Note != "Reference group", ]
    expect_equal(nrow(non_ref), 1L)
    # AggressiveTx has all events; TR should be < 1 (shorter survival)
    expect_lt(non_ref$Time_Ratio, 1)
    expect_lt(non_ref$Upper_CrI,   1)
  })

  test_that("bayesian_survival: event counts match input data", {
    skip_bayes_surv()
    df <- make_surv_simple()
    te <- get_surv_result()$treatment_effects
    expect_equal(te$Events[te$Group == "AggressiveTx"], 7L)
    expect_equal(te$Events[te$Group == "Control"],      0L)
  })

  test_that("bayesian_survival: Median_Survival is positive and finite", {
    skip_bayes_surv()
    te <- get_surv_result()$treatment_effects
    expect_true(all(is.finite(te$Median_Survival) & te$Median_Survival > 0))
  })

  # ── posterior_summary ──────────────────────────────────────────────────────
  test_that("bayesian_survival: posterior_summary is a data frame with Rhat column", {
    skip_bayes_surv()
    ps <- get_surv_result()$posterior_summary
    expect_s3_class(ps, "data.frame")
    expect_true("Rhat" %in% colnames(ps))
  })

  # ── mcmc_diagnostics ───────────────────────────────────────────────────────
  test_that("bayesian_survival: mcmc_diagnostics has Converged column", {
    skip_bayes_surv()
    diag <- get_surv_result()$mcmc_diagnostics
    expect_true("Converged" %in% colnames(diag))
    expect_type(diag$Converged, "logical")
  })

  # ── survival_data passthrough ───────────────────────────────────────────────
  test_that("bayesian_survival: survival_data mirrors input structure", {
    skip_bayes_surv()
    sd <- get_surv_result()$survival_data
    expect_s3_class(sd, "data.frame")
    expect_true(all(c("Time", "Event", "Treatment") %in% colnames(sd)))
    expect_equal(nrow(sd), nrow(make_surv_simple()))
  })

  # ── Input validation ───────────────────────────────────────────────────────
  test_that("bayesian_survival: errors on missing required column", {
    skip_bayes_surv()
    df <- make_surv_simple()
    df$Day <- NULL
    expect_error(
      bayesian_survival(df, time_column = "Day"),
      "Missing required columns"
    )
  })

  test_that("bayesian_survival: errors on invalid reference group", {
    skip_bayes_surv()
    expect_error(
      suppressMessages(bayesian_survival(
        make_surv_simple(),
        reference_group = "NonExistentGroup",
        n_chains = 2L, n_iter = 200L, plots = FALSE, verbose = FALSE
      )),
      "not found"
    )
  })

  # ── summary metadata ───────────────────────────────────────────────────────
  test_that("bayesian_survival: summary contains expected metadata fields", {
    skip_bayes_surv()
    s <- get_surv_result()$summary
    expect_equal(s$data_description$total_events, 7L)
    expect_equal(s$data_description$reference_group, "Control")
    expect_equal(s$model_specification$family, "weibull")
  })
})

# ── Family-specific: lognormal ──────────────────────────────────────────────

test_that("bayesian_survival: lognormal family runs and returns HR = NA for non-reference", {
  skip_bayes_surv()
  df <- make_surv_simple()
  res <- suppressWarnings(suppressMessages(
    bayesian_survival(
      df, family = "lognormal", reference_group = "Control",
      include_cage_effect = FALSE, n_chains = 2L, n_iter = 500L,
      plots = FALSE, verbose = FALSE, seed = 1L
    )
  ))
  expect_equal(res$family_used, "lognormal")
  non_ref <- res$treatment_effects[res$treatment_effects$Note != "Reference group", ]
  expect_true(is.na(non_ref$HR))
})

# ── Frailty ────────────────────────────────────────────────────────────────

test_that("bayesian_survival: frailty_used is FALSE when no cage_column supplied", {
  skip_bayes_surv()
  df <- make_surv_simple()
  res <- suppressWarnings(suppressMessages(
    bayesian_survival(
      df, cage_column = NULL, reference_group = "Control",
      n_chains = 2L, n_iter = 500L, plots = FALSE, verbose = FALSE, seed = 2L
    )
  ))
  expect_false(res$frailty_used)
})
