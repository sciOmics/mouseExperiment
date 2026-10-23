# The R script on the other analysis tabs (CODE_REVIEW.md R20-AA): each test
# runs the tab through testServer(), downloads only the script (which writes
# the data it reads), runs it in a clean environment, and requires the
# script's tables to equal the dashboard's.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

sample_rv <- function(file, ...) {
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", file, package = "mouseExperiment"),
                                    stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL, settings = list(...))
}
master_rv <- function() sample_rv("all_features_synthetic_data.csv",
  col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
  volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
  col_weight = "Weight", col_dose = "Dose", col_survival_event = "Survival_Censor",
  control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse")
combo_rv <- function() sample_rv("combo_treatment_synthetic_data.csv",
  col_id = "ID", col_treatment = "Treatment", col_cage = "Cage", col_dose = "",
  time_method = "from_date", col_date_for_day = "Date", reference_date = "2025-03-24",
  date_order = "auto", volume_method = "calculate", col_length = "Length", col_width = "Width",
  volume_formula = "modified_ellipsoid", control_group = "Control", volume_units = "mm3",
  randomisation_unit = "mouse")

# Tick only the R script on tab `prefix`, download, unzip and run it.
run_tab_script <- function(output, session, prefix) {
  groups <- c("tables", "plots", "diagnostics", "data", "report")
  ins <- stats::setNames(rep(list(character(0)), length(groups)),
                         paste0(prefix, "_dl_group_", groups))
  ins[[paste0(prefix, "_dl_group_code")]] <- "script"
  do.call(session$setInputs, ins)
  z <- output[[paste0(prefix, "_dl_zip")]]
  dir <- tempfile("script_")
  utils::unzip(z, exdir = dir)
  env <- new.env(parent = globalenv())
  old <- setwd(dir)
  on.exit(setwd(old))
  # R22.19: R's default generator, as in a fresh session. The dashboard's
  # worker draws with L'Ecuyer-CMRG; a seeded result must not depend on it.
  old_kind <- RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]), add = TRUE)
  suppressMessages(suppressWarnings(utils::capture.output(
    sys.source("analysis.R", envir = env, keep.source = FALSE))))
  list(env = env, files = utils::unzip(z, list = TRUE)$Name,
       script = readLines(file.path(dir, "analysis.R")))
}
same_table <- function(a, b, label) {
  expect_false(is.null(a), info = label)
  expect_equal(as.data.frame(a), as.data.frame(b), tolerance = 1e-8, ignore_attr = TRUE,
               info = label)
}

test_that("R20-AA: Survival's script reproduces the hazard ratios", {
  rv <- master_rv(); shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = rv), {
    session$setInputs(run_surv_analysis = 1)
    res <- rv$surv_results
    expect_false(isTRUE(res$error))
    run <- run_tab_script(output, session, "surv")
    expect_true(all(c("analysis.R", "data/analysed_data.rds") %in% run$files))
    expect_true(any(grepl("result <- survival_statistics(", run$script, fixed = TRUE)))
    same_table(run$env$result$results, res$results, "results")
  }))
})

test_that("R20-AA: Survival's script uses the events the dashboard derived from the volume limit", {
  rv <- master_rv(); shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = rv), {
    session$setInputs(surv_event_source = "volume", surv_vol_threshold = 1500)
    session$setInputs(run_surv_analysis = 1)
    res <- rv$surv_results
    skip_if(!identical(res$event_src_used, "volume"), "the volume-limit events were not used")
    run <- run_tab_script(output, session, "surv")
    expect_true(any(grepl("volume limit of 1500", run$script, fixed = TRUE)))
    same_table(run$env$result$results, res$results, "results")
  }))
})

test_that("R20-AA: Toxicity's body-weight script repeats the model and the threshold", {
  rv <- master_rv(); shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE, tox_tumor_density = 1,
                      tox_covariates = "volume", tox_estimation = "REML", tox_threshold = 20,
                      tox_baseline_day = "first", run_tox_analysis = 1)
    res <- rv$tox_results
    skip_if(is.null(res$mixed_model), "the body-weight model did not finish here")
    run <- run_tab_script(output, session, "tox")
    expect_true(any(grepl("bw <- analyze_body_weight(", run$script, fixed = TRUE)))
    expect_true(any(grepl("wl <- weight_loss_threshold(", run$script, fixed = TRUE)))
    # Both calls read the one data file.
    expect_equal(sum(grepl("read_analysed(", run$script, fixed = TRUE)), 1L)
    same_table(run$env$bw$fixed_effects, res$mixed_model$fixed_effects, "fixed effects")
    same_table(run$env$bw$pairwise_comparisons, res$mixed_model$pairwise_comparisons,
               "rate comparisons")
    same_table(run$env$wl$event_data, res$weight_loss$event_data, "weight-loss events")
  }))
})

test_that("R20-AA: Toxicity's therapeutic-window script reproduces the window", {
  rv <- master_rv(); shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "et", tox_adjust_tumor = TRUE, tox_tumor_density = 1,
                      tox_endpoint_method = "model", tox_n_boot = 200, tox_tol_threshold = 20,
                      run_tox_analysis = 1)
    res <- rv$tox_results
    expect_false(isTRUE(res$twm$error))
    run <- run_tab_script(output, session, "tox")
    expect_true(any(grepl("boot_seed = 20260730L", run$script, fixed = TRUE)))
    same_table(run$env$twm$window_table, res$twm$window_table, "window")
  }))
})

test_that("R20-AA: Drug Synergy's script reproduces one day's metrics, with the bootstrap", {
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = combo_rv()), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "28", syn_endpoint_method = "model", syn_n_boot = 50,
                      run_syn_analysis = 1)
    res <- rv$syn_results
    expect_false(isTRUE(res$error))
    run <- run_tab_script(output, session, "syn")
    # Days and volumes came from dates and calipers; the data hold them.
    same_table(run$env$result$synergy_metrics, res$synergy_metrics, "metrics")
    same_table(run$env$result$synergy_ci, res$synergy_ci, "intervals")
  }))
})

test_that("R20-AA: Drug Synergy's script reproduces every evaluable day", {
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = combo_rv()), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "all", syn_endpoint_method = "model", syn_n_boot = 0,
                      run_syn_analysis = 1)
    res <- rv$syn_results
    expect_false(isTRUE(res$error))
    run <- run_tab_script(output, session, "syn")
    expect_true(any(grepl("analyze_drug_synergy_over_time(", run$script, fixed = TRUE)))
    same_table(run$env$result$synergy_summary, res$synergy_summary, "summary")
  }))
})

test_that("R20-AA: Dose Response's script reproduces the curve and the TGI", {
  rv <- master_rv(); shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    session$setInputs(dr_treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
                      dr_time_point = "default", run_dr_analysis = 1)
    res <- rv$dr_results
    expect_false(isTRUE(res$error))
    run <- run_tab_script(output, session, "dr")
    expect_true(any(grepl('treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High")',
                          run$script, fixed = TRUE)))
    same_table(run$env$result$tgi_table, res$tgi_table, "TGI")
    same_table(g("dr_stats_table")(run$env$result), g("dr_stats_table")(res), "statistics")
  }))
})

test_that("R20-AA: Power's analytic script stacks one call per effect size", {
  shiny::testServer(g("power_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(pw_mode = "apriori_analytic", pw_ap_es_mode = "cohens_d",
                      pw_ap_cohens_d = "0.8, 1.2", pw_ap_n_groups = "3", pw_ap_submode = "find_n",
                      pw_alpha = "0.05", pw_power = "0.80", pw_ap_p_adjust = "bonferroni",
                      pw_ap_n_comparisons = 2, pw_ap_dropout = 10, run_pw_analysis = 1)
    res <- rv$ap_results
    expect_length(res$call_records, 2L)
    run <- run_tab_script(output, session, "pw")
    # No data: the script has none to read.
    expect_identical(run$files, c("README.txt", "analysis.R"))
    expect_true(any(grepl("It needs no data", run$script, fixed = TRUE)))
    same_table(run$env$scenario_table, res$scenario_table, "scenarios")
  })
})

test_that("R20-AA: Power's delta-and-SD script includes the SD sensitivity", {
  shiny::testServer(g("power_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(pw_mode = "apriori_analytic", pw_ap_es_mode = "delta_sd",
                      pw_ap_delta = 0.5, pw_ap_pooled_sd = 0.6, pw_ap_n_groups = "2",
                      pw_ap_submode = "find_n", pw_alpha = "0.05", pw_power = "0.80",
                      pw_ap_p_adjust = "bonferroni", pw_ap_n_comparisons = 1, pw_ap_dropout = 0,
                      run_pw_analysis = 1)
    res <- rv$ap_results
    run <- run_tab_script(output, session, "pw")
    same_table(run$env$power$scenario_table, res$scenario_table, "scenarios")
    same_table(run$env$power$sensitivity_table, res$sensitivity_table, "sensitivity")
  })
})

test_that("R20-AA: Power's LMM simulation is seeded, and its script gives the same power", {
  shiny::testServer(g("power_server"), args = list(rv = shiny::reactiveValues()), {
    # A small effect, so the power is not near 0 or 1 under any stream and a
    # different random stream shows (R22.19).
    session$setInputs(pw_mode = "apriori_lmm", pw_lmm_n_per_group = "5", pw_lmm_n_groups = 2,
                      pw_lmm_n_sims = 40, pw_lmm_n_timepoints = 5, pw_lmm_day_interval = 3,
                      pw_lmm_tx_effect = 0.04, pw_alpha = "0.05", run_pw_analysis = 1)
    res <- rv$ap_results
    skip_if(!identical(res$mode, "apriori_lmm"), "the simulation did not finish here")
    expect_true(all(res$scenario_table$Power > 0.1 & res$scenario_table$Power < 0.9))
    run <- run_tab_script(output, session, "pw")
    expect_true(any(grepl("seed = 20260730L", run$script, fixed = TRUE)))
    same_table(run$env$sim$power_table, res$scenario_table, "simulated power")
  })
})
