# =============================================================================
# Backend contract tests — CODE_REVIEW.md R3.D9 / I.2 / R4.D2
#
# The dashboard reads specific named fields out of backend result objects. Round 3
# added several (attrition, synergy_ci, twm_ci, censoring_summary, variance_test,
# permutation_test) and the render code added in v26.07.004 depends on their exact
# names. Nothing verified those names, which is the same gap that let the suite
# pass against a backend six releases stale (R4.D2).
#
# These tests call the backend the way the dashboard does and assert the contract.
# They are deliberately about *field presence and shape*, not statistical values —
# the statistics are the backend's own tests' job. A rename upstream should break
# here rather than silently blank a panel.
# =============================================================================

# A study with the real dropout mechanism: per-animal growth rates, animals removed
# when volume crosses a limit, so attrition is non-zero and the estimand matters.
contract_df <- function(seed = 42, limit = 2000) {
  set.seed(seed)
  rows <- list(); cg <- 0L
  arms <- list(Control = 0.115, DrugA = 0.075, DrugB = 0.060, Combo = 0.030)
  for (tx in names(arms)) {
    for (c_i in 1:2) {
      cg <- cg + 1L
      for (m in 1:4) {
        rate <- stats::rnorm(1, arms[[tx]], 0.02)
        days <- seq(0, 28, 4)
        v <- exp(log(150) + stats::rnorm(1, 0, 0.15) + rate * days +
                   stats::rnorm(length(days), 0, 0.08))
        hit  <- which(v >= limit)
        keep <- if (length(hit)) seq_len(hit[1]) else seq_along(days)
        rows[[length(rows) + 1L]] <- data.frame(
          ID = paste0("C", cg, "_m", m), Cage = paste0("C", cg), Treatment = tx,
          Day = days[keep], Volume = v[keep],
          Weight = 25 - 0.02 * days[keep] + stats::rnorm(length(keep), 0, 0.15),
          stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

test_that("contract: synergy result carries the fields the Synergy tab renders", {
  res <- suppressWarnings(suppressMessages(mouseExperiment::analyze_drug_synergy(
    contract_df(), treatment_column = "Treatment", volume_column = "Volume",
    time_column = "Day", id_column = "ID",
    drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", eval_time_point = 28,
    endpoint_method = "model", n_boot = 200, boot_seed = 1, verbose = FALSE)))

  # syn_ci_table / syn_ci_verdict_note
  expect_false(is.null(res$synergy_ci))
  expect_true(all(c("Metric", "Estimate", "CI_Lower", "CI_Upper") %in%
                    names(res$synergy_ci)))
  expect_true("Bliss_Excess_FE" %in% res$synergy_ci$Metric,
              info = "syn_ci_verdict_note keys on this exact metric name")

  # syn_attrition_table / syn_attrition_heading
  expect_false(is.null(res$attrition))
  expect_true(all(c("Treatment", "N_Enrolled", "N_At_Endpoint", "Pct_Lost") %in%
                    names(res$attrition)))
  expect_false(is.null(res$endpoint_method))

  # The existing metrics table still keys on Metric / Value.
  expect_true(all(c("Metric", "Value") %in% names(res$synergy_metrics)))
})

test_that("contract: TWM result carries the ranking bounds and attrition", {
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::therapeutic_window_metric(
      contract_df(), volume_column = "Volume", weight_column = "Weight",
      time_column = "Day", treatment_column = "Treatment", id_column = "ID",
      reference_group = "Control", adjust_tumor_weight = FALSE,
      endpoint_method = "model", n_boot = 200, boot_seed = 1)))

  # tox_twm_ranking_note reads these off twm_table directly.
  expect_true(all(c("Treatment", "TGI", "TWM", "TWM_Lower", "TWM_Upper") %in%
                    names(res$twm_table)))
  # tox_twm_attrition_ui / _table
  expect_false(is.null(res$attrition))
  expect_true(all(c("N_Enrolled", "N_At_Endpoint", "Pct_Lost") %in%
                    names(res$attrition)))
  expect_false(is.null(res$endpoint_day))
  expect_false(is.null(res$endpoint_method))
})

test_that("contract: weight-loss result carries the competing-risks fields", {
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::weight_loss_threshold(
      contract_df(), weight_column = "Weight", time_column = "Day",
      treatment_column = "Treatment", id_column = "ID", cage_column = "Cage",
      adjust_tumor_weight = FALSE, threshold = 0.05)))

  # tox_km_competing_note / tox_km_censoring_table
  expect_false(is.null(res$n_competing_risk))
  expect_true(is.numeric(res$n_competing_risk))
  expect_false(is.null(res$censoring_summary))
  expect_true(all(c("Treatment", "Censor_Type", "Freq") %in%
                    names(res$censoring_summary)))
  # The KM renderer still re-fits from event_data on Time / Event / Treatment.
  expect_true(all(c("Time", "Event", "Treatment") %in% names(res$event_data)))
})

test_that("contract: AUC result carries variance_test for the ANOVA tab", {
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(
      contract_df(), model_type = "auc", plots = FALSE, verbose = FALSE)))

  # tg_variance_test_ui reads p_value; the table renders the whole frame.
  expect_false(is.null(res$variance_test))
  expect_true("p_value" %in% names(res$variance_test))
  # The Welch omnibus replaced aov(); the ANOVA tab renders this frame.
  expect_true("Method" %in% names(res$anova))
  expect_match(res$anova$Method[1], "Welch")
})

test_that("contract: permutation result carries every field the panel reads", {
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(
      contract_df(), plots = FALSE, verbose = FALSE,
      permutation_test = mouseExperiment::perm_spec(
        unit = "mouse", n_perm = 49, seed = 1))))

  pt <- res$permutation_test
  expect_false(is.null(pt))
  # tg_perm_test_ui reads each of these by name.
  for (f in c("p_value", "unit", "n_units", "n_perm_converged", "exhaustive",
              "lrt_p_asymptotic", "min_attainable_p", "n_distinct_statistics",
              "method")) {
    expect_true(f %in% names(pt), info = paste("missing field:", f))
  }
  expect_true(is.character(pt$method))
  expect_true(pt$unit %in% c("mouse", "cage"))

  # Absent unless requested, which is what the "Not run" branch relies on.
  plain <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(contract_df(), plots = FALSE,
                                             verbose = FALSE)))
  expect_null(plain$permutation_test)
})

test_that("contract: comparison provenance fields the summary panel prints exist", {
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(
      contract_df(), plots = FALSE, verbose = FALSE,
      comparison_family = "all_pairs", p_adjust_method = "tukey")))

  # comparisons_title() and the COMPARISONS summary block read these.
  expect_identical(res$comparison_family, "all_pairs")
  expect_identical(res$p_adjust_method_used, "tukey")
  # The CAGE STRUCTURE block reads these.
  expect_false(is.null(res$summary$model_specification$cage_structure))
  expect_false(is.null(res$summary$model_specification$cage_placement))
  expect_false(is.null(res$cage_analysis$structure))
})

# ---- R4.D5: interval columns resolved from meta, not sniffed -----------------

test_that("R4.D5: the old sniff list missed the column the backend actually uses", {
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(contract_df(), plots = FALSE,
                                             verbose = FALSE)))
  te <- res$treatment_effects

  # The historical sniff list. This is the bug: treatment_effects carries
  # "Lower_CL", which matches none of these, so has_ci was always FALSE and the
  # Treatment Effects plot drew SE-derived bars instead of the model CI.
  legacy <- intersect(c("lower.CL", "lower_CL", "Lower", "lower"), colnames(te))[1]
  expect_true(is.na(legacy),
              info = "if this ever matches, the historical bug is gone upstream")
  expect_true("Lower_CL" %in% colnames(te))

  # The replacement resolves it from the backend's own declaration.
  ic <- resolve_interval_cols(te, res)
  expect_identical(ic$lower, "Lower_CL")
  expect_identical(ic$upper, "Upper_CL")
})

test_that("R4.D5: resolution falls back cleanly for a meta-less result", {
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(contract_df(), plots = FALSE,
                                             verbose = FALSE)))
  te <- res$treatment_effects

  # An older backend returns no meta; the widened fallback list must still find
  # the columns rather than silently yielding NA.
  ic <- resolve_interval_cols(te, result = NULL)
  expect_identical(ic$lower, "Lower_CL")
  expect_identical(ic$upper, "Upper_CL")

  # A Bayesian-style table uses CrI names; both routes must handle it.
  bayes_like <- data.frame(Group = "A", Adjusted_Mean = 1,
                           Lower_CrI = 0.5, Upper_CrI = 1.5)
  ic2 <- resolve_interval_cols(bayes_like, result = NULL)
  expect_identical(ic2$lower, "Lower_CrI")
})

test_that("R4.D5: meta declares which kind of interval it is", {
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(contract_df(), plots = FALSE,
                                             verbose = FALSE)))
  # A consumer can now label the axis correctly without inspecting names:
  # "confidence" vs "credible" is a real distinction, not a formatting choice.
  expect_identical(res$meta$inference, "frequentist")
  expect_identical(res$meta$interval_type, "confidence")
  expect_identical(res$meta$estimate_scale, "log volume")
})

# ---- R3.D6: scale metadata must record the outcome, not the request ---------

test_that("R3.D6: the AUC path reports the transform it actually applied", {
  # The divergence that makes reading input$tg_transform wrong: AUC ignores the
  # transform argument and computes on raw volumes, so a run submitted with the
  # default transform = "log" must not be labelled "log".
  auc <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(
      contract_df(), model_type = "auc", transform = "log",
      plots = FALSE, verbose = FALSE)))

  expect_identical(auc$meta$transform_used, "none",
                   info = "AUC computes on raw volumes regardless of `transform`")
  expect_identical(auc$meta$model_type_used, "auc")

  # And the lme4 path does honour it, so the field is not simply always "none".
  lme <- suppressWarnings(suppressMessages(
    mouseExperiment::tumor_growth_statistics(
      contract_df(), model_type = "lme4", transform = "log",
      plots = FALSE, verbose = FALSE)))
  expect_identical(lme$meta$transform_used, "log")
})

test_that("R3.D6: the module prefers meta over the input when recording scale", {
  # Guards the resolution order in mod_tumor_growth.R. Reading the input would
  # tag AUC output "log" and mis-scale any consumer not gated on model_type.
  src <- paste(readLines(testthat::test_path("..", "..", "R",
                                             "mod_tumor_growth.R"),
                         warn = FALSE), collapse = "\n")
  expect_match(src, "result\\$transform_used\\s*<-\\s*result\\$meta\\$transform_used")
  expect_match(src, "result\\$model_type_used\\s*<-\\s*result\\$meta\\$model_type_used")
})
