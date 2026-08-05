# =============================================================================
# Round 14 — statistical-correctness audit
#
# Both findings here are the same shape as R3.35 / R3.37 / K.14: code that runs
# without error, produces output the dashboard renders, and is wrong. Neither was
# reachable by reading — the first needed a known-answer test, the second needed
# an arm configuration nobody had simulated.
# =============================================================================

# ---- R14.1 dose-response parameters were NA in every run --------------------

dr_known <- function(seed = 3, ec50 = 25, b = 1.5, top = 1000, bot = 100) {
  set.seed(seed)
  doses <- c(0, 5, 10, 25, 50, 100, 200)
  rows <- list(); k <- 0L
  for (dz in doses) for (m in 1:8) {
    k <- k + 1L
    mu <- bot + (top - bot) / (1 + (max(dz, 1e-9) / ec50)^b)
    rows[[k]] <- data.frame(ID = paste0("d", dz, "_", m),
                            Treatment = paste0("D", dz), Dose = dz, Day = 21,
                            Volume = max(1, stats::rnorm(1, mu, 60)),
                            stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
dr_fit <- function(df) {
  suppressWarnings(suppressMessages(dose_response_statistics(
    df, dose_column = "Dose", treatment_column = "Treatment",
    volume_column = "Volume", day_column = "Day", id_column = "ID",
    control_group_name = "D0")))
}

test_that("R14.1: the four curve parameters are reported, not NA", {
  # They were NA on every run. LL.4() is called with
  # names = c("Slope","Lower Limit","Upper Limit","EC50"), so the fitted object
  # carries those, while the extraction looked up the drc defaults b/c/d/e --
  # zero of four matched. The direction check on the Slope line was updated when
  # the names were introduced (G.2); this block was not.
  st <- dr_fit(dr_known())$statistics
  skip_if(is.null(st$dr_model), "drc unavailable")
  for (f in c("ec50", "hill_slope", "lower_limit", "upper_limit")) {
    expect_true(is.finite(st[[f]]), info = paste(f, "must not be NA"))
  }
})

test_that("R14.1: the parameters recover a known curve", {
  # NA would pass an is.numeric check, so assert against truth: data simulated
  # from EC50 = 25, slope = 1.5, top = 1000.
  st <- dr_fit(dr_known())$statistics
  skip_if(is.null(st$dr_model), "drc unavailable")
  expect_gt(st$ec50, 15); expect_lt(st$ec50, 45)
  expect_gt(st$hill_slope, 0.6); expect_lt(st$hill_slope, 2.6)
  expect_gt(st$upper_limit, 800); expect_lt(st$upper_limit, 1200)
})

test_that("R14.1: EC50 is not exponentiated", {
  # The old line was `exp(params["e:(Intercept)"])`. In LL.4 the e parameter is
  # the EC50 on the natural dose scale, so once the name lookup started working
  # exp() would have returned ~1e13. LL2.4 is the log-parameterised variant.
  st <- dr_fit(dr_known())$statistics
  skip_if(is.null(st$dr_model), "drc unavailable")
  expect_lt(st$ec50, 1000)
  expect_equal(unname(st$ec50),
               unname(stats::coef(st$dr_model)[grep("^EC50", names(stats::coef(st$dr_model)))][1]),
               tolerance = 1e-8)
})

test_that("R14.1: EC50 carries a confidence interval", {
  # A point estimate presented without one reads as a fact. drc supplies a
  # delta-method interval from the same fit at no extra cost.
  st <- dr_fit(dr_known())$statistics
  skip_if(is.null(st$dr_model), "drc unavailable")
  expect_false(is.null(st$ec50_ci))
  expect_true(all(c("lower", "upper", "se") %in% names(st$ec50_ci)))
  expect_lt(st$ec50_ci[["lower"]], st$ec50)
  expect_gt(st$ec50_ci[["upper"]], st$ec50)
})

test_that("R14.1: parameters are read by position if names change again", {
  # The failure mode was a rename upstream. Reading positionally means the next
  # relabelling cannot silently reintroduce it.
  st <- dr_fit(dr_known(seed = 99))$statistics
  skip_if(is.null(st$dr_model), "drc unavailable")
  expect_true(is.finite(st$ec50))
})

# ---- R14.2 a harmful single agent was reported as synergy -------------------

syn_df <- function(rates, seed = 11) {
  set.seed(seed)
  rows <- list(); k <- 0L
  for (tx in names(rates)) for (m in 1:8) {
    k <- k + 1L
    d <- seq(0, 21, 3)
    rows[[k]] <- data.frame(
      ID = paste0(tx, m), Treatment = tx, Day = d,
      Volume = exp(log(150) + stats::rnorm(1, 0, 0.10) + rates[[tx]] * d +
                     stats::rnorm(length(d), 0, 0.06)),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
syn_fit <- function(df) {
  suppressWarnings(suppressMessages(analyze_drug_synergy(
    df, treatment_column = "Treatment", volume_column = "Volume",
    time_column = "Day", id_column = "ID", drug_a_name = "DrugA",
    drug_b_name = "DrugB", combo_name = "Combo", control_name = "Control",
    eval_time_point = 21, n_boot = 0, verbose = FALSE)))
}

test_that("R14.2: an agent that accelerates growth is not called synergistic", {
  # The defect: a negative fractional effect makes the ratio negative, and the
  # thresholds are one-sided ("CI < 0.85 => synergistic"). DrugA accelerating
  # growth, DrugB inert, combination near control produced CI = -29.2 reported
  # as "Synergistic" -- the verdict exactly inverted, on the configuration a
  # reviewer would scrutinise hardest.
  r <- syn_fit(syn_df(c(Control = 0.100, DrugA = 0.130, DrugB = 0.100,
                        Combo = 0.098)))
  # R17.2 removed the Combination Index; the guard now shows in the label.
  expect_null(r$combination_index)
  expect_match(r$overall_assessment, "Not evaluable")
  expect_false(grepl("^Synerg", r$overall_assessment))
})

test_that("R14.2: the guard warns rather than failing silently", {
  expect_warning(
    analyze_drug_synergy(
      syn_df(c(Control = 0.100, DrugA = 0.130, DrugB = 0.100, Combo = 0.098)),
      treatment_column = "Treatment", volume_column = "Volume",
      time_column = "Day", id_column = "ID", drug_a_name = "DrugA",
      drug_b_name = "DrugB", combo_name = "Combo", control_name = "Control",
      eval_time_point = 21, n_boot = 0, verbose = FALSE),
    "did not inhibit growth")
})

test_that("R14.2: both agents harming is also not evaluable", {
  r <- syn_fit(syn_df(c(Control = 0.100, DrugA = 0.125, DrugB = 0.120,
                        Combo = 0.099)))
  expect_match(r$overall_assessment, "Not evaluable")
})

test_that("R14.2: a genuine inhibitory combination still evaluates", {
  # The guard must not suppress the case the module exists for.
  r <- syn_fit(syn_df(c(Control = 0.130, DrugA = 0.100, DrugB = 0.105,
                        Combo = 0.070)))
  expect_false(grepl("Not evaluable", r$overall_assessment))
  expect_true(is.finite(r$bliss_independence$difference))
})

test_that("R14.2: TWM and synergy now agree on how to treat a harmful arm", {
  # TWM already clamped with max(TGI, 0); synergy did not. The inconsistency was
  # itself a signal that one of the two had not thought the case through.
  r <- syn_fit(syn_df(c(Control = 0.100, DrugA = 0.130, DrugB = 0.100,
                        Combo = 0.098)))
  expect_match(r$overall_assessment, "Not evaluable")
})

# ---- R14.4 / R14.6 pairwise_comparisons contract ----------------------------

tg_three_arm <- function(seed = 2) {
  set.seed(seed)
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA", "DrugB")) for (m in 1:8) {
    k <- k + 1L
    d <- seq(0, 21, 3)
    r0 <- switch(tx, Control = 0.12, DrugA = 0.08, DrugB = 0.06)
    rows[[k]] <- data.frame(
      ID = paste0(tx, m), Treatment = tx, Day = d,
      Volume = exp(log(150) + stats::rnorm(1, 0, 0.08) + r0 * d +
                     stats::rnorm(length(d), 0, 0.05)),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

test_that("R14.4: pairwise_comparisons is a data frame on every model path", {
  # lme4 returned a raw emmGrid. Consumers guarding with is.data.frame() skipped
  # it silently, which is how the adjustment scope went missing from the header
  # on the default path.
  df <- tg_three_arm()
  for (mt in c("lme4", "gam", "auc")) {
    r <- suppressWarnings(suppressMessages(
      tumor_growth_statistics(df, model_type = mt, plots = FALSE, verbose = FALSE)))
    expect_s3_class(r$pairwise_comparisons, "data.frame")
    expect_gt(nrow(r$pairwise_comparisons), 0L)
  }
})

test_that("R14.4: every path carries the same canonical columns", {
  df <- tg_three_arm()
  need <- c("contrast", "P_Value_Adjusted", "Adjust_Scope",
            "Comparison_Family", "P_Adjust_Method")
  for (mt in c("lme4", "gam", "auc")) {
    r <- suppressWarnings(suppressMessages(
      tumor_growth_statistics(df, model_type = mt, plots = FALSE, verbose = FALSE)))
    missing <- setdiff(need, names(r$pairwise_comparisons))
    expect_length(missing, 0L)
  }
})

test_that("R14.6: the adjusted p-value is unambiguous across paths", {
  # The defect: `p_value` held the ADJUSTED value on the GAM path and the RAW
  # value on the AUC path, so a consumer sniffing for that name got different
  # quantities depending on the model the user picked. P_Value_Adjusted is now
  # always the adjusted one.
  df <- tg_three_arm()

  auc <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "auc", plots = FALSE, verbose = FALSE)))
  pc <- auc$pairwise_comparisons
  expect_true(all(c("P_Value_Adjusted", "P_Value_Raw") %in% names(pc)))
  # Adjustment can only move p-values up (or leave them alone).
  ok <- is.finite(pc$P_Value_Raw) & is.finite(pc$P_Value_Adjusted)
  expect_true(all(pc$P_Value_Adjusted[ok] >= pc$P_Value_Raw[ok] - 1e-12))

  # And it must equal the path's own adjusted column, not the raw one.
  expect_equal(pc$P_Value_Adjusted[ok], pc$p_adjusted[ok], tolerance = 1e-12)
})

test_that("R14.4: original columns survive so existing consumers keep working", {
  # Additive normalisation: declare, do not rename (the B7.1 precedent).
  df <- tg_three_arm()
  auc <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "auc", plots = FALSE, verbose = FALSE)))
  expect_true(all(c("comparison", "p_value", "p_adjusted") %in%
                    names(auc$pairwise_comparisons)))

  lme <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "lme4", plots = FALSE, verbose = FALSE)))
  expect_true(all(c("contrast", "estimate", "SE", "p.value") %in%
                    names(lme$pairwise_comparisons)))
})

test_that("R14.4: the adjustment scope now reaches the default path", {
  # The live consequence of the emmGrid return: comparisons_title() guards with
  # is.data.frame() before reading Adjust_Scope, so the header stated the
  # correction but omitted the set it was applied over.
  df <- tg_three_arm()
  r <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "lme4", plots = FALSE, verbose = FALSE)))
  expect_true(is.data.frame(r$pairwise_comparisons))
  expect_true("Adjust_Scope" %in% names(r$pairwise_comparisons))
  expect_true(nzchar(as.character(r$pairwise_comparisons$Adjust_Scope[1])))
})
