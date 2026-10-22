# "Download Results" replaces "Download All" (2026-10-22, CODE_REVIEW.md
# R20-D.X): the user chooses which tables, plots, diagnostics, data and report
# go into one ZIP, with a README. The tests fail on 26.10.008.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}
all_features_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                                package = "mouseExperimentDashboard"),
                                    stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "Cage", col_weight = "Weight",
                    control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
}
zip_names <- function(path) utils::unzip(path, list = TRUE)$Name
read_zipped <- function(path, name) {
  d <- tempfile(); on.exit(unlink(d, recursive = TRUE))
  utils::unzip(path, files = name, exdir = d)
  readLines(file.path(d, name), warn = FALSE)
}

test_that("write_results_zip() puts each kind in its folder, in the chosen plot format", {
  item <- g("dl_item")
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  items <- list(
    item("t", "tables", "A table", "a_table", "table", function() data.frame(a = 1:2)),
    item("p", "plots", "A plot", "a_plot", "plot", function() p),
    item("d", "diagnostics", "Unavailable", "nothing", "plot", function() NULL),
    item("r", "report", "Report", "report.html", "file",
         function(path) writeLines("<html></html>", path)))
  # SVG uses Cairo, which the server image has; a Mac without XQuartz does not.
  svg_ok <- tryCatch(suppressWarnings({
    grDevices::svg(tempfile(fileext = ".svg")); grDevices::dev.off(); TRUE
  }), error = function(e) FALSE)
  for (fmt in c("png", "pdf", if (svg_ok) "svg")) {
    z <- tempfile(fileext = ".zip")
    g("write_results_zip")(z, items, fmt = fmt, header = "Header line")
    expect_setequal(zip_names(z), c("README.txt", "tables/a_table.csv",
                                    paste0("plots/a_plot.", fmt), "report.html"))
    readme <- read_zipped(z, "README.txt")
    expect_identical(readme[1], "Header line")
    expect_true(any(grepl("tables/a_table.csv", readme, fixed = TRUE)))
    # A chosen item that could not be made is named, not silently dropped.
    expect_true(any(grepl("Unavailable", readme, fixed = TRUE)))
  }
})

test_that("Tumour Growth lists what the run produced, and exports only what is ticked", {
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    # Before a run there is nothing to choose.
    expect_match(html_of(output$tg_dl_choices), "Run an analysis first", fixed = TRUE)
    session$setInputs(tg_model_type = "lme4", tg_transform = "log",
                      tg_comparison_type = "reference", tg_p_adjust = "dunnett")
    session$setInputs(run_tg_analysis = 1)
    ui <- html_of(output$tg_dl_choices)
    for (lbl in c("Tables (CSV)", "Plots", "Diagnostics", "Data (CSV)", "Report",
                  "ANOVA", "Comparisons", "Growth curves", "Residual Q-Q",
                  "Analysed data", "HTML report")) {
      expect_match(ui, lbl, fixed = TRUE, info = lbl)
    }
    # LME4: no AUC tables, no Bayesian diagnostics.
    expect_false(grepl("Per-animal AUC", ui, fixed = TRUE))
    expect_false(grepl("MCMC", ui, fixed = TRUE))

    session$setInputs(tg_dl_group_tables = c("anova", "comparisons"),
                      tg_dl_group_plots = "plot_growth",
                      tg_dl_group_diagnostics = "diag_qq",
                      tg_dl_group_data = character(0),
                      tg_dl_group_report = "report",
                      tg_dl_plot_format = "pdf")
    z <- output$tg_dl_zip
    expect_setequal(zip_names(z), c("README.txt", "tables/anova.csv", "tables/comparisons.csv",
                                    "plots/growth_curves.pdf", "diagnostics/residual_qq.pdf",
                                    "report.html"))
    readme <- read_zipped(z, "README.txt")
    expect_true(any(grepl("model_type: lme4", readme, fixed = TRUE)))
    expect_true(any(grepl("mouseExperiment ", readme, fixed = TRUE)))
  }))
})

test_that("an AUC run offers the AUC tables under their own names", {
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "auc", tg_comparison_type = "reference",
                      tg_p_adjust = "holm")
    session$setInputs(run_tg_analysis = 1)
    ui <- html_of(output$tg_dl_choices)
    expect_match(ui, "AUC by arm", fixed = TRUE)
    expect_match(ui, "Per-animal AUC", fixed = TRUE)
  }))
})

test_that("the tab is called Download Results and offers PNG, PDF and SVG", {
  ui <- html_of(g("tumorGrowthUI")("tg"))
  expect_match(ui, "Download Results", fixed = TRUE)
  expect_false(grepl("Download All", ui, fixed = TRUE))
  for (f in c('value="png"', 'value="pdf"', 'value="svg"')) expect_match(ui, f, fixed = TRUE)
})

# ---- The other analysis tabs (26.10.009, CODE_REVIEW.md R20-D.Y) ----------------

master_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                                package = "mouseExperimentDashboard"),
                                    stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "Cage", col_weight = "Weight",
                    col_dose = "Dose", col_survival_event = "Survival_Censor",
                    control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
}
combo_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "combo_treatment_synthetic_data.csv",
                                                package = "mouseExperimentDashboard"),
                                    stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", time_method = "from_date", col_date_for_day = "Date",
                    reference_date = "2025-03-24", date_order = "auto",
                    volume_method = "calculate", col_length = "Length", col_width = "Width",
                    volume_formula = "modified_ellipsoid", control_group = "Control",
                    volume_units = "mm3", randomisation_unit = "mouse"))
}
expect_listed <- function(ui, labels) for (l in labels) expect_match(ui, l, fixed = TRUE, info = l)

test_that("Survival: the Cox run's tables, plots and diagnostics go into the ZIP", {
  rv <- master_rv()
  # Without the two arms that have no events, the Cox model (not Firth's) is
  # used, so the proportional-hazards test exists to export.
  d <- shiny::isolate(rv$uploaded_data)
  shiny::isolate(rv$uploaded_data <- d[!(d$Treatment == "Drug_A High" |
                                         d$Treatment == "Drug_A Mid + Drug_B"), ])
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = rv), {
    session$setInputs(run_surv_analysis = 1)
    expect_listed(html_of(output$surv_dl_choices),
                  c("Hazard ratios", "Risk table", "Kaplan-Meier curves", "Forest plot",
                    "Proportional-hazards test (CSV)", "Schoenfeld residuals"))
    session$setInputs(surv_dl_group_tables = c("results", "risk"),
                      surv_dl_group_plots = c("plot_km", "plot_forest"),
                      surv_dl_group_diagnostics = c("ph_table", "ph_plot"))
    expect_setequal(zip_names(output$surv_dl_zip),
                    c("README.txt", "tables/survival_results.csv", "tables/risk_table.csv",
                      "plots/kaplan_meier.png", "plots/forest_plot.png",
                      "diagnostics/ph_test.csv", "diagnostics/schoenfeld_residuals.png"))
  }))
})

test_that("Toxicity: Download Results replaces Download All, with the plots and diagnostics", {
  ui <- html_of(g("toxicityUI")("tox"))
  expect_match(ui, "Download Results", fixed = TRUE)
  expect_false(grepl("Download All", ui, fixed = TRUE))
  rv <- master_rv(); shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE, tox_tumor_density = 1,
                      tox_covariates = "volume", tox_estimation = "REML", tox_threshold = 20,
                      tox_baseline_day = "first", run_tox_analysis = 1)
    expect_listed(html_of(output$tox_dl_choices),
                  c("Weight change per day", "Rate comparisons", "Weight-loss events",
                    "Body-weight trajectories", "Time to weight loss", "Residual Q-Q",
                    "Body-weight data"))
    session$setInputs(tox_dl_group_tables = "rates",
                      tox_dl_group_plots = c("plot_traj", "plot_km"),
                      tox_dl_group_diagnostics = "diag_qq", tox_dl_group_data = "data")
    expect_setequal(zip_names(output$tox_dl_zip),
                    c("README.txt", "tables/weight_change_rates.csv",
                      "plots/body_weight_trajectories.png", "plots/weight_loss_km.png",
                      "diagnostics/residual_qq.png", "data/body_weight_data.csv"))
  }))
})

test_that("Synergy: the metrics and the plot go into the ZIP", {
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = combo_rv()), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "28", syn_endpoint_method = "model", syn_n_boot = 0,
                      run_syn_analysis = 1)
    expect_listed(html_of(output$syn_dl_choices), c("Synergy metrics", "Synergy plot"))
    invisible(output$syn_plot_main)
    session$setInputs(syn_dl_group_tables = "metrics", syn_dl_group_plots = "plot_main")
    expect_setequal(zip_names(output$syn_dl_zip),
                    c("README.txt", "tables/synergy_metrics.csv", "plots/synergy_plot.png"))
  }))
})

test_that("Dose Response: the statistics table is the one the single download gave", {
  rv <- master_rv(); shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    session$setInputs(dr_treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
                      dr_time_point = "default", run_dr_analysis = 1)
    expect_listed(html_of(output$dr_dl_choices),
                  c("Summary by dose", "Dose-response statistics", "TGI by dose"))
    session$setInputs(dr_dl_group_tables = "stats")
    z <- output$dr_dl_zip
    expect_true("tables/dose_response_statistics.csv" %in% zip_names(z))
    single <- utils::read.csv(output$download_dr_stats)
    d <- tempfile(); utils::unzip(z, "tables/dose_response_statistics.csv", exdir = d)
    expect_equal(utils::read.csv(file.path(d, "tables/dose_response_statistics.csv")), single)
  }))
})

test_that("Power: the analytic run's tables go into the ZIP", {
  shiny::testServer(g("power_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(pw_mode = "apriori_analytic", pw_ap_es_mode = "cohens_d",
                      pw_ap_cohens_d = "1", pw_ap_n_groups = "2", pw_ap_submode = "find_n",
                      pw_alpha = "0.05", pw_power = "0.80", pw_ap_p_adjust = "bonferroni",
                      pw_ap_n_comparisons = 1, pw_ap_dropout = 0, run_pw_analysis = 1)
    expect_listed(html_of(output$pw_dl_choices), c("Scenarios", "Power curves"))
    session$setInputs(pw_dl_group_tables = "scenarios")
    expect_setequal(zip_names(output$pw_dl_zip), c("README.txt", "tables/scenarios.csv"))
  })
})
