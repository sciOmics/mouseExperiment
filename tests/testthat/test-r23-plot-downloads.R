# Download Results draws every plot with the current settings, also a plot
# that was never shown (CODE_REVIEW.md Round 23, R23.14, R23-H). Drug Synergy,
# Dose Response and Power exported the last plot drawn on screen, so a plot
# type that had not been opened was missing from the ZIP. The tests below
# never read a plot output before downloading.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
zip_names <- function(path) utils::unzip(path, list = TRUE)$Name
sample_rv <- function(file, ...) {
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", file, package = "mouseExperiment"),
                                    stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL, settings = list(...))
}
no_groups <- function(prefix) {
  groups <- c("tables", "plots", "diagnostics", "data", "report", "code")
  stats::setNames(rep(list(character(0)), length(groups)), paste0(prefix, "_dl_group_", groups))
}

test_that("R23.14: Drug Synergy's plots go into the ZIP, styled, without being shown first", {
  rv <- sample_rv("combo_treatment_synthetic_data.csv",
    col_id = "ID", col_treatment = "Treatment", col_cage = "Cage", col_dose = "",
    time_method = "from_date", col_date_for_day = "Date", reference_date = "2025-03-24",
    date_order = "auto", volume_method = "calculate", col_length = "Length", col_width = "Width",
    volume_formula = "modified_ellipsoid", control_group = "Control", volume_units = "mm3",
    randomisation_unit = "mouse")
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = rv), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "28", syn_endpoint_method = "model", syn_n_boot = 0,
                      syn_plot_title = "Bliss at day 28", run_syn_analysis = 1)
    ins <- no_groups("syn")
    ins$syn_dl_group_plots <- "plot_main"
    do.call(session$setInputs, ins)
    expect_true("plots/synergy_plot.png" %in% zip_names(output$syn_dl_zip))
    expect_identical(syn_main_plot_r()$labels$title, "Bliss at day 28")
    # A single day has no TGI-over-time plot, and none is offered.
    expect_false(grepl("tgi_by_arm", paste(as.character(output$syn_dl_choices), collapse = ""), fixed = TRUE))
    # Over every evaluable day: the trend and TGI plots, never shown either.
    session$setInputs(syn_time_point = "all", run_syn_analysis = 2)
    ins$syn_dl_group_plots <- c("plot_main", "plot_trend", "plot_tgi")
    do.call(session$setInputs, ins)
    expect_true(all(c("plots/synergy_trend.png", "plots/tgi_by_arm.png") %in% zip_names(output$syn_dl_zip)))
    expect_identical(syn_tgi_plot_r()$labels$title, "Bliss at day 28")
  }))
})

test_that("R23.14: Dose Response's two plots go into the ZIP without being shown first", {
  rv <- sample_rv("all_features_synthetic_data.csv",
    col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
    col_weight = "Weight", col_dose = "Dose", col_survival_event = "Survival_Censor",
    control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse")
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    session$setInputs(dr_treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
                      dr_time_point = "default", dr_plot_title = "Dose series", run_dr_analysis = 1)
    ins <- no_groups("dr")
    ins$dr_dl_group_plots <- c("plot_scatter", "plot_tgi")
    do.call(session$setInputs, ins)
    expect_true(all(c("plots/dose_response.png", "plots/tgi_by_dose.png") %in% zip_names(output$dr_dl_zip)))
    expect_identical(dr_scatter_plot_r()$labels$title, "Dose series")
  }))
})

test_that("R23.14: Power's curves go into the ZIP without being shown first", {
  shiny::testServer(g("power_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(pw_mode = "apriori_analytic", pw_ap_es_mode = "cohens_d",
                      pw_ap_cohens_d = "1", pw_ap_n_groups = "2", pw_ap_submode = "find_n",
                      pw_alpha = "0.05", pw_power = "0.80", pw_ap_p_adjust = "bonferroni",
                      pw_ap_n_comparisons = 1, pw_ap_dropout = 0, run_pw_analysis = 1)
    ins <- no_groups("pw")
    ins$pw_dl_group_plots <- "plot_curves"
    do.call(session$setInputs, ins)
    expect_true("plots/power_curves.png" %in% zip_names(output$pw_dl_zip))
  })
})

test_that("R23.14: a plot that only explains why there is nothing to draw is not exported", {
  m <- g("message_plot")("Nothing to draw")
  expect_s3_class(m, "ggplot")
  expect_null(g("shown_plot")(m))
  p <- ggplot2::ggplot()
  expect_identical(g("shown_plot")(p), p)
  # Help no longer tells users to open a plot before downloading it.
  help <- paste(as.character(g("helpUI")()), collapse = "")
  expect_false(grepl("exported as last shown", help, fixed = TRUE))
})
