# Regression tests for CODE_REVIEW.md Round 20 fixes, labelled by finding ID.
# Each asserts observable behaviour, not the presence of code: two earlier
# "fixes" in this package existed in the source and never executed.

# ---- R20.7 ------------------------------------------------------------------

test_that("R20.7: k >= 3 powers each treated-vs-control comparison", {
  n_for <- function(k, d) {
    apriori_power_analysis(effect_size = d, n_groups = k, alpha = 0.05,
                           target_power = 0.8)$scenario_table$Required_N
  }
  # power.t.test at alpha / (k - 1), the reference column in R20.7. The old
  # omnibus path returned 27, 8, 10 and 7.
  expect_equal(n_for(3, 0.5), 78)
  expect_equal(n_for(3, 1.0), 21)
  expect_equal(n_for(4, 0.8), 35)
  expect_equal(n_for(4, 1.0), 23)

  # At the old recommendation (n = 8 for k = 3, d = 1) each comparison has about
  # a third of the intended power; a 20,000-study simulation of Bonferroni
  # t-tests against control gave 0.338.
  fp <- apriori_power_analysis(effect_size = 1, n_groups = 3, alpha = 0.05,
                               n_per_group = 8, mode = "find_power")
  expect_equal(fp$scenario_table$Achieved_Power, 0.338, tolerance = 0.005)
  expect_equal(fp$scenario_table$Alpha_Per_Comparison, 0.025)
  expect_match(fp$method_note, "treated-vs-control comparison")
})

test_that("R20.7: the default correction changes nothing with two groups", {
  a <- apriori_power_analysis(effect_size = 0.8, n_groups = 2, alpha = 0.05,
                              target_power = 0.8)
  b <- apriori_power_analysis(effect_size = 0.8, n_groups = 2, alpha = 0.05,
                              target_power = 0.8, p_adjust_method = "none")
  expect_identical(a$scenario_table$Required_N, b$scenario_table$Required_N)
  expect_null(a$method_note)
})

test_that("R20.7: the SD sensitivity table uses the per-comparison alpha", {
  r <- apriori_power_analysis(delta = 0.5, pooled_sd = 0.5, n_groups = 4,
                              alpha = 0.05, target_power = 0.8)
  mid <- r$sensitivity_table[r$sensitivity_table$SD_Change == "0%", ]
  # d = 1 at the assumed SD, so the middle row must match the scenario table.
  expect_equal(mid$Required_N, r$scenario_table$Required_N[1])
  expect_equal(mid$Alpha_Per_Comparison, 0.05 / 3)
})

# ---- R20.8 ------------------------------------------------------------------

# Four arms with log-linear growth; `rates` are per-day growth rates.
r20_synergy_df <- function(rates, n = 8, days = seq(0, 28, by = 3.5), seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(names(rates), function(arm) {
    do.call(rbind, lapply(seq_len(n), function(i) {
      b0 <- stats::rnorm(1, log(100), 0.1)
      r  <- rates[[arm]] + stats::rnorm(1, 0, 0.005)
      data.frame(ID = paste(arm, i), Treatment = arm, Day = days,
                 Volume = exp(b0 + r * days + stats::rnorm(length(days), 0, 0.08)))
    }))
  }))
}

test_that("R20.8: an agent that accelerates growth yields no Bliss quantities", {
  # The R20.8 scenario: DrugA accelerates growth, DrugB is inert, the
  # combination equals control.
  df <- r20_synergy_df(c(Control = 0.10, DrugA = 0.13, DrugB = 0.10, Combo = 0.10))
  res <- suppressWarnings(analyze_drug_synergy(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE, n_boot = 200, boot_seed = 1))

  expect_false(res$evaluable)
  expect_true(is.na(res$bliss_independence$synergy))      # was TRUE
  expect_true(is.na(res$bliss_independence$difference))
  expect_true(is.na(res$bliss_independence$expected_effect))
  excess <- res$synergy_ci[res$synergy_ci$Metric == "Bliss_Excess_FE", ]
  expect_true(is.na(excess$CI_Lower) && is.na(excess$CI_Upper))  # was [0.80, 3.79]
  expect_match(res$overall_assessment, "Not evaluable")

  ot <- suppressWarnings(suppressMessages(analyze_drug_synergy_over_time(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE)))
  expect_false(any(ot$synergy_summary$Evaluable))
  expect_equal(nrow(ot$peak_bliss_synergy), 0L)            # was day 28, 173.8 %

  p <- suppressWarnings(plot_synergy_trend(ot))
  built <- suppressWarnings(ggplot2::ggplot_build(p))
  ribbon_rows <- vapply(seq_along(p$layers), function(i) {
    if (inherits(p$layers[[i]]$geom, "GeomRibbon")) nrow(built$data[[i]]) else 0L
  }, integer(1))
  expect_equal(sum(ribbon_rows), 0L)                       # was drawn on all 9 days
})

test_that("R20.8: inhibitory agents are unaffected", {
  df <- r20_synergy_df(c(Control = 0.10, DrugA = 0.07, DrugB = 0.075, Combo = 0.03))
  res <- suppressWarnings(analyze_drug_synergy(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE, n_boot = 200, boot_seed = 1))
  expect_true(res$evaluable)
  expect_true(is.logical(res$bliss_independence$synergy) &&
                !is.na(res$bliss_independence$synergy))
  excess <- res$synergy_ci[res$synergy_ci$Metric == "Bliss_Excess_FE", ]
  expect_true(is.finite(excess$CI_Lower) && is.finite(excess$CI_Upper))

  ot <- suppressWarnings(suppressMessages(analyze_drug_synergy_over_time(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE)))
  expect_equal(nrow(ot$peak_bliss_synergy), 1L)
  expect_true(ot$synergy_summary$Evaluable[
    ot$synergy_summary$Time_Point == ot$peak_bliss_synergy$Time_Point])
})
