# Regression tests for step 4 of the Round 20/21 implementation (CODE_REVIEW.md
# R20-D.L): the evaluable-day rule relayed to users, one source for TGI, the
# model-based AUC and random slopes by default. Each test fails on v26.09.003
# run against backend 0.26.0.

skip_if_not_installed("shiny", minimum_version = "1.5.0")
skip_if_not_installed("mouseExperiment", minimum_version = "0.26.0")

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}

demo_csv <- function(f) {
  utils::read.csv(system.file("sample_data", f, package = "mouseExperimentDashboard"),
                  stringsAsFactors = FALSE)
}

# The Combo demo: the control arm thins out after day 28 (3 of 8 on day 30).
combo_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = demo_csv("combo_treatment_synthetic_data.csv"), data_valid = TRUE,
    session_id = NULL,
    settings = list(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", time_method = "from_date",
                    col_date_for_day = "Date", reference_date = "2025-03-24",
                    date_order = "auto", volume_method = "calculate",
                    col_length = "Length", col_width = "Width",
                    volume_formula = "modified_ellipsoid",
                    control_group = "Control", volume_units = "mm3",
                    randomisation_unit = "mouse"))
}

master_rv <- function(...) {
  shiny::reactiveValues(
    uploaded_data = demo_csv("master_synthetic_data.csv"), data_valid = TRUE,
    session_id = NULL,
    settings = utils::modifyList(list(
      col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
      volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
      col_weight = "Weight", col_dose = "Dose", col_survival_event = "Survival_Censor",
      control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse"),
      list(...)))
}

# ---- R20-K: the rule, relayed ---------------------------------------------------

test_that("R20-K: the dashboard states the backend's evaluable-day rule, word for word", {
  d <- data.frame(ID = rep(1:3, each = 2), Treatment = "A", Day = rep(c(0, 7), 3))
  expect_identical(g("evaluable_rule_text")(), mouseExperiment::evaluable_days(d)$rule)
})

test_that("R20-K: the Study Design panel states the rule", {
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(demo_dataset = "master", load_demo_data = 1)
    expect_match(html_of(output$study_design_ui), "Evaluable days")
  })
})

test_that("R20-K: Synergy offers only evaluable days, and names the day and the exclusions", {
  rv <- combo_rv()
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = rv), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1")
    ui <- html_of(output$syn_time_point_ui)
    expect_match(ui, "All evaluable days")
    expect_match(ui, 'value="28"')
    expect_false(grepl('value="30"', ui))             # 3 of 8 controls on study
    expect_false(grepl('value="32"', ui))
    session$setInputs(syn_time_point = "28", syn_endpoint_method = "model",
                      syn_n_boot = 200, run_syn_analysis = 1)
    expect_equal(rv$syn_results$eval_time_point, 28)
    note <- html_of(output$syn_evaluability_note)
    expect_match(note, "Evaluated on day 28, the last day on which all four arms are evaluable")
    expect_match(note, "day 30 \\(Control 3/8\\)")
  }))
})

test_that("R20.D37: the Synergy TGI plot draws the backend's TGIs, not a re-derivation", {
  rv <- combo_rv()
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = rv), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "all", syn_endpoint_method = "model",
                      syn_n_boot = 200, run_syn_analysis = 1)
    ss <- rv$syn_results$synergy_summary
    expect_equal(max(ss$Time_Point), 28)              # evaluable days only
    long <- g("syn_tgi_long")(ss, rv$syn_plot_inputs)
    a <- long[long$Treatment == "aPD1", ]
    expect_equal(a$mean_tgi, ss$TGI_Drug_A)
    expect_equal(a$lo, ss$TGI_Drug_A_Lower)
    expect_match(html_of(output$syn_evaluability_note), "Analysed on every evaluable day")
  }))
})

test_that("R20-K: the therapeutic window names its day and says the ratio has no interval", {
  rv <- master_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "et", tox_adjust_tumor = TRUE,
                      tox_tumor_density = 1, tox_endpoint_method = "model",
                      tox_n_boot = 200, run_tox_analysis = 1)
    twm <- rv$tox_results$twm
    expect_false(isTRUE(twm$error))
    expect_match(html_of(output$tox_twm_evaluability_note),
                 paste0("Evaluated on day ", twm$endpoint_day))
    expect_match(html_of(output$tox_twm_ranking_note), "no interval")
  }))
})

test_that("R20.17 / R20.D48: Dose Response analyses the last evaluable day and shows the backend's TGI", {
  rv <- master_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    ui <- html_of(output$dr_time_point_ui)
    expect_match(ui, "Last evaluable day")
    expect_false(grepl("All dates", ui, fixed = TRUE))
    session$setInputs(dr_time_point = "default", run_dr_analysis = 1)
    res <- rv$dr_results
    expect_false(isTRUE(res$error))
    ev <- mouseExperiment::evaluable_days(
      transform(demo_csv("master_synthetic_data.csv"), DoseGroup = as.character(Dose)),
      treatment_column = "DoseGroup", cage_column = "Cage", volume_column = "Volume")
    expect_equal(res$endpoint_day, ev$last_day)
    expect_false(is.null(res$tgi_table))
    expect_match(html_of(output$dr_evaluability_note),
                 paste0("Evaluated on day ", ev$last_day))
    expect_true(nzchar(output$dr_plot_tgi$src))
  }))
})

# ---- R20.5: random slopes by default ------------------------------------------

test_that("R20.5: Tumour Growth defaults to random slopes and reports the structure used", {
  ui <- html_of(g("tumorGrowthUI")("tg"))
  expect_match(ui, '<option value="slope" selected>')
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = master_rv()), {
    session$setInputs(tg_model_type = "lme4", tg_random_effects = "slope",
                      tg_p_adjust = "dunnett", run_tg_analysis = 1)
    re <- rv$tg_results$random_effects
    expect_identical(re$requested, "slope")
    expect_true("F value" %in% names(rv$tg_results$anova))
    expect_match(paste(output$tg_summary, collapse = " "), "Satterthwaite")
  }))
})

# ---- R20.4: the model-based AUC -----------------------------------------------

test_that("R20.4: the AUC comparisons and summary are the backend's model-based ones", {
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = master_rv()), {
    session$setInputs(tg_model_type = "auc", tg_p_adjust = "holm", run_tg_analysis = 1)
    res <- rv$tg_results
    expect_false(isTRUE(res$error))
    pw <- res$all_pairwise_comparisons
    expect_true(all(c("ratio", "ratio.lower", "ratio.upper", "p.value") %in% names(pw)))
    expect_true(all(pw$ratio.lower <= pw$ratio & pw$ratio <= pw$ratio.upper))
    txt <- paste(output$tg_summary, collapse = " ")
    expect_match(txt, "Window: day 0 to day 16")
  }))
})
