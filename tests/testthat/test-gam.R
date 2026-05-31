# =============================================================================
# Tests for the GAM (gamm4) model option added to
#   - tumor_growth_statistics(model_type = "gam")
#   - analyze_body_weight(model_type = "gam")
#   - bayesian_tumor_growth(model_type = "gam")
#
# All tests skip cleanly when gamm4 (or brms, for the Bayesian path) is not
# installed so CI without these heavy deps stays green.
# =============================================================================

skip_gam <- function() {
  skip_if_not_installed("gamm4")
  skip_if_not_installed("mgcv")
  skip_if_not_installed("emmeans")
}

# Build longitudinal data with a non-linear "treated" trajectory so GAM is the
# right tool. Control grows exponentially; Treatment grows exponentially then
# plateaus from day 14 onward — a pattern a linear LMM can't represent.
make_nonlinear_tg <- function() {
  set.seed(202605)
  days <- c(0, 4, 8, 12, 16, 20, 24)

  ctrl_rate <- 0.08
  make_mouse <- function(id, cage, treatment) {
    noise <- rnorm(length(days), 0, 0.06)
    log_vol <- if (treatment == "Control") {
      log(200) + ctrl_rate * days + noise
    } else {
      # exponential then plateau at day 14
      pre  <- 0.05 * pmin(days, 14)
      post <- 0.005 * pmax(days - 14, 0)
      log(200) + pre + post + noise
    }
    data.frame(
      ID        = id,
      Cage      = cage,
      Treatment = treatment,
      Day       = days,
      Volume    = round(exp(log_vol), 2),
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_mouse("C01", "C1", "Control"),
    make_mouse("C02", "C1", "Control"),
    make_mouse("C03", "C2", "Control"),
    make_mouse("C04", "C2", "Control"),
    make_mouse("C05", "C3", "Control"),
    make_mouse("C06", "C3", "Control"),
    make_mouse("T01", "T1", "TreatmentA"),
    make_mouse("T02", "T1", "TreatmentA"),
    make_mouse("T03", "T2", "TreatmentA"),
    make_mouse("T04", "T2", "TreatmentA"),
    make_mouse("T05", "T3", "TreatmentA"),
    make_mouse("T06", "T3", "TreatmentA")
  )
}

# ---------------------------------------------------------------------------
# tumor_growth_statistics(model_type = "gam")
# ---------------------------------------------------------------------------

test_that("tumor_growth_statistics(gam) returns the same top-level field set as lme4", {
  skip_gam()
  df  <- make_nonlinear_tg()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df, model_type = "gam", transform = "log",
      reference_group = "Control", plots = FALSE, verbose = FALSE
    )
  ))

  required <- c("model", "model_type_used", "anova", "summary",
                "pairwise_comparisons", "treatment_effects",
                "treatment_effects_over_time", "growth_rates",
                "cage_analysis", "model_selection", "diagnostics",
                "data_summary")
  expect_true(all(required %in% names(res)),
              info = paste("Missing:", paste(setdiff(required, names(res)),
                                             collapse = ", ")))
  expect_equal(res$model_type_used, "gam")
})

test_that("tumor_growth_statistics(gam) treatment_effects has one row per group", {
  skip_gam()
  df  <- make_nonlinear_tg()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df, model_type = "gam", transform = "log",
      reference_group = "Control", plots = FALSE
    )
  ))

  te <- as.data.frame(res$treatment_effects)
  expect_equal(nrow(te), 2L)
  expect_true(all(c("Group", "Adjusted_Mean", "Lower_CL", "Upper_CL")
                  %in% colnames(te)))
})

test_that("tumor_growth_statistics(gam) pairwise table has one row per quantile day", {
  skip_gam()
  df  <- make_nonlinear_tg()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df, model_type = "gam", transform = "log",
      reference_group = "Control", plots = FALSE
    )
  ))

  pw <- as.data.frame(res$pairwise_comparisons)
  # Five study-day quantiles, one Treatment − Control contrast → 5 rows
  expect_equal(nrow(pw), 5L)
  expect_true("Day" %in% colnames(pw))
  expect_true("estimate" %in% colnames(pw))
})

test_that("tumor_growth_statistics(gam) diagnostics include k_check and deviance_explained", {
  skip_gam()
  df  <- make_nonlinear_tg()
  res <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(
      df, model_type = "gam", include_diagnostics = TRUE,
      plots = FALSE, reference_group = "Control"
    )
  ))
  expect_true(!is.null(res$diagnostics$k_check))
  expect_true(is.numeric(res$diagnostics$deviance_explained))
})

test_that("polynomial_degree argument was removed from tumor_growth_statistics", {
  skip_gam()
  df <- make_nonlinear_tg()
  expect_error(
    suppressWarnings(suppressMessages(
      tumor_growth_statistics(
        df, model_type = "lme4", polynomial_degree = 2,
        reference_group = "Control", plots = FALSE
      )
    )),
    regexp = "unused argument"
  )
})

# ---------------------------------------------------------------------------
# analyze_body_weight(model_type = "gam")
# ---------------------------------------------------------------------------

test_that("analyze_body_weight(gam) returns expected structure", {
  skip_gam()
  df <- make_bw_simple()  # from helper-fixtures.R
  res <- suppressWarnings(suppressMessages(
    analyze_body_weight(
      df, weight_column = "Weight", time_column = "Day",
      treatment_column = "Treatment", id_column = "ID",
      cage_column = "Cage", model_type = "gam",
      adjust_tumor_weight = FALSE, covariates = character(0)
    )
  ))

  expect_equal(res$model_info$model_type, "gam")
  expect_true(!is.null(res$model$gam))
  expect_true(!is.null(res$model$mer))
  expect_true(is.numeric(res$model_info$smoother_k))
  expect_gte(res$model_info$smoother_k, 3L)
})

# ---------------------------------------------------------------------------
# bayesian_tumor_growth(model_type = "gam") — heavy, gated by brms
# ---------------------------------------------------------------------------

test_that("bayesian_tumor_growth(gam) returns model_type_used 'bayes_tg_gam'", {
  skip_gam()
  skip_if_not_installed("brms")

  df <- make_nonlinear_tg()
  res <- suppressWarnings(suppressMessages(
    bayesian_tumor_growth(
      df, model_type = "gam",
      reference_group = "Control",
      n_chains = 2L, n_iter = 250L, seed = 42L,
      return_model = FALSE, plots = FALSE, verbose = FALSE
    )
  ))

  expect_equal(res$model_type_used, "bayes_tg_gam")
  expect_s3_class(res$mcmc_diagnostics, "data.frame")
})
