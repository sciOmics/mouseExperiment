# =============================================================================
# apriori_power_simulation() — CODE_REVIEW.md K.13
#
# Surfaced by the first working coverage baseline (K.10): this file measured
# **0.00%**. The function is exported, is wired into the dashboard's Power module,
# and had no test of any kind.
#
# That matters beyond the number. Round 1 §2.2 rewrote its data-generating model
# (linear scale -> log scale, because simulating linear and fitting log gave
# anti-conservative power), §3.7 inverted its p-value extraction (the LRT branch
# was unreachable), and §J.4 turned `baseline_sd` from a ghost parameter into a
# functional one. Three substantive fixes to a function nothing exercised — the
# same shape as R3.35 and R3.37, and only luck that it still works.
#
# Simulations are stochastic, so these are properties that must hold rather than
# fixed numbers: monotonicity in N and in effect size, correct behaviour under the
# null, and every argument having an observable effect (the K.4 bug class).
# n_simulations is kept low; each replicate fits two mixed models.
# =============================================================================

test_that("K.13: returns the documented structure", {
  res <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = c(6L, 10L), n_simulations = 20L, seed = 1L)))

  expect_type(res, "list")
  expect_true("power_table" %in% names(res))
  expect_s3_class(res$power_table, "data.frame")
  expect_true(all(c("N_Per_Group", "Power") %in% names(res$power_table)))
  expect_equal(nrow(res$power_table), 2L)
  # Power is a proportion.
  expect_true(all(res$power_table$Power >= 0 & res$power_table$Power <= 1))
})

test_that("K.13: power increases with sample size", {
  # The defining property. If this fails the simulation is not measuring power.
  res <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = c(4L, 20L), n_simulations = 40L,
    treatment_effect = 0.06, seed = 7L)))

  p_small <- res$power_table$Power[res$power_table$N_Per_Group == 4L]
  p_large <- res$power_table$Power[res$power_table$N_Per_Group == 20L]
  expect_gt(p_large, p_small)
})

test_that("K.13: power increases with effect size", {
  run <- function(eff) suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 40L, treatment_effect = eff,
    seed = 11L)))$power_table$Power

  expect_gt(run(0.10), run(0.02))
})

test_that("K.13: a null effect gives power near the nominal alpha", {
  # With treatment_effect = 0 the "power" is the type-I error rate, which should
  # sit near alpha rather than near 0 or 1. This is the check that would catch a
  # sign error or a broken test-extraction path -- which §3.7 found once already,
  # when the LRT branch was unreachable.
  res <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 10L, n_simulations = 60L, treatment_effect = 0,
    alpha = 0.05, seed = 3L)))

  # Generous bounds: 60 replicates give a wide interval around 0.05.
  expect_lt(res$power_table$Power, 0.25)
  expect_gte(res$power_table$Power, 0)
})

test_that("K.13: alpha changes the answer", {
  # K.4 bug class -- a documented argument that does not reach the computation.
  res <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 40L, treatment_effect = 0.05,
    alpha = c(0.01, 0.20), seed = 5L)))

  expect_equal(nrow(res$power_table), 2L)
  p01 <- res$power_table$Power[res$power_table$Alpha == 0.01]
  p20 <- res$power_table$Power[res$power_table$Alpha == 0.20]
  # A laxer threshold cannot yield less power.
  expect_gte(p20, p01)
})

test_that("K.14: baseline_sd cannot move the reported power, by construction", {
  # §J.4 found baseline_sd stored in params but never read, and wired it into the
  # data-generating process with the comment that this "makes the baseline_sd
  # argument materially affect the simulation". It does not, and no wiring could:
  # the test is on Treatment:Day, a contrast of slopes, while baseline_sd injects
  # only per-mouse intercept variation, which (Day | ID) absorbs exactly.
  #
  # This test pins the invariance rather than asserting the fix worked, because
  # the invariance is the correct behaviour. If it ever breaks, either the fitted
  # model lost its random intercept or the estimand changed -- both worth knowing.
  run <- function(bsd) suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 30L, baseline_sd = bsd,
    treatment_effect = 0.05, seed = 2L)))$power_table$Power

  expect_equal(run(1), run(200),
               info = "slope contrast must stay orthogonal to baseline spread")
})

test_that("K.14: baseline_sd and random_intercept_sd are aliased", {
  # Both are drawn independently and summed into the same per-mouse offset, so
  # only their root-sum-square is identifiable. Documented so nobody spends
  # calibration effort trying to separate them.
  a <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 20L, baseline_sd = 0,
    random_intercept_sd = 0.30, treatment_effect = 0.05, seed = 4L)))
  b <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 20L, baseline_sd = 30,
    random_intercept_sd = 0.30, treatment_effect = 0.05, seed = 4L)))
  expect_equal(a$power_table$Power, b$power_table$Power)
})

# ---- K.15: dropout ---------------------------------------------------------

test_that("K.15: dropout_limit removes observations and is reported", {
  res <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 20L, control_growth_rate = 0.13,
    treatment_effect = 0.02, timepoints = seq(0, 28, 4),
    dropout_limit = 700, seed = 21L)))

  expect_false(is.null(res$attrition))
  expect_true(all(c("Pct_Obs_Lost", "Pct_Mice_Truncated", "Dropout_Limit") %in%
                    names(res$attrition)))
  expect_gt(res$attrition$Pct_Obs_Lost, 0)
  expect_identical(res$params$dropout_limit, 700)
})

test_that("K.15: the default is complete data, so attrition is exactly zero", {
  # Inf must preserve the historical behaviour -- this argument is additive.
  res <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 20L, timepoints = seq(0, 28, 4),
    seed = 21L)))
  expect_equal(res$attrition$Pct_Obs_Lost, 0)
  expect_equal(res$attrition$Pct_Mice_Truncated, 0)
})

test_that("K.15: dropout lowers power relative to complete data", {
  # The documented direction: MAR dropout costs information, not validity.
  run <- function(lim) suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 120L, control_growth_rate = 0.13,
    treatment_effect = 0.018, random_slope_sd = 0.012, residual_sd = 0.10,
    random_intercept_sd = 0.20, timepoints = seq(0, 28, 4),
    dropout_limit = lim, seed = 31L)))$power_table$Power

  expect_lt(run(500), run(Inf))
})

test_that("K.15: an invalid dropout_limit is rejected", {
  expect_error(apriori_power_simulation(dropout_limit = -1),
               regexp = "dropout_limit")
  expect_error(apriori_power_simulation(dropout_limit = c(100, 200)),
               regexp = "dropout_limit")
})

test_that("K.13: results are reproducible given a seed", {
  run <- function() suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_simulations = 20L, seed = 99L)))$power_table$Power
  expect_equal(run(), run())
})

test_that("K.13: invalid arguments are rejected", {
  expect_error(apriori_power_simulation(n_groups = 1L), regexp = "n_groups")
  expect_error(apriori_power_simulation(n_simulations = 2L),
               regexp = "n_simulations")
  expect_error(apriori_power_simulation(timepoints = 0), regexp = "timepoints")
})

test_that("K.13: multi-group designs are supported (J.5 sibling)", {
  res <- suppressWarnings(suppressMessages(apriori_power_simulation(
    n_per_group = 8L, n_groups = 3L, n_simulations = 20L,
    treatment_effect = 0.08, seed = 13L)))
  expect_s3_class(res$power_table, "data.frame")
  expect_true(all(is.finite(res$power_table$Power)))
})
