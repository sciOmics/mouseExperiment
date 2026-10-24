# Round 22 proposal F9 (CODE_REVIEW.md R22-J, 0.48.0): the R script redraws the
# plots, with the dashboard's styling, through the exported style_plot().

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

test_that("F9: style_plot() applies each option, and leaves the rest as they were", {
  p <- plot_tumor_growth(master_synthetic_data, cage_column = "Cage")
  s <- style_plot(p, theme = "classic", title = "Growth", title_align = "left", text_size = 15,
                  x_label = "Study day", legend_position = "inplot_tr", y_scale = "log10",
                  x_limits = c(0, 30))
  expect_identical(s$labels$title, "Growth")
  expect_identical(s$labels$x, "Study day")
  expect_identical(s$labels$y, p$labels$y)                       # left as it was
  expect_identical(s$theme$plot.title$hjust, 0)
  expect_identical(s$theme$text$size, 15)
  expect_identical(s$theme$legend.position, "inside")
  expect_identical(s$coordinates$limits$x, c(0, 30))
  expect_true(any(vapply(s$scales$scales, function(sc) "y" %in% sc$aesthetics &&
                                             identical(sc$trans$name, "log-10"), NA)))
  # A later call does not undo an earlier one (the legend, then nothing else).
  s2 <- style_plot(s, legend_position = "bottom")
  expect_identical(s2$theme$plot.title$hjust, 0)
  expect_identical(s2$theme$legend.position, "bottom")
  expect_error(style_plot(p, theme = "dark"), "should be one of")
})

test_that("F9: Tumor Growth's script redraws the growth curves and the waterfall as on the dashboard", {
  rv <- shiny::reactiveValues(
    uploaded_data = master_synthetic_data, data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference", tg_p_adjust = "holm",
                      tg_resp_min_day = 10, tg_delay_volume = 1000, tg_delay_n_boot = 0,
                      tg_waterfall_fill = "treatment", tg_waterfall_cap = 150)
    session$setInputs(run_tg_analysis = 1)
    groups <- rv$tg_plot_groups
    cols <- stats::setNames(grDevices::hcl.colors(length(groups)), groups)
    do.call(session$setInputs, c(
      stats::setNames(as.list(unname(cols)), paste0("tg_color_", make.names(groups))),
      stats::setNames(as.list(as.character(seq_along(groups) + 14L)), paste0("tg_shape_", make.names(groups))),
      list(tg_plot_title = "Tumour growth", tg_opt_text_size = 13, tg_opt_legend_pos = "bottom",
           tg_opt_yscale = "log10")))
    session$elapse(1000)   # the labels are debounced
    ins <- stats::setNames(rep(list(character(0)), 5),
                           paste0("tg_dl_group_", c("tables", "plots", "diagnostics", "data", "report")))
    do.call(session$setInputs, c(ins, list(tg_dl_group_code = "script")))
    dir <- tempfile("plots_")
    utils::unzip(output$tg_dl_zip, exdir = dir)
    script <- readLines(file.path(dir, "analysis.R"))
    expect_true(any(grepl("p_growth <- plot_tumor_growth(", script, fixed = TRUE)))
    expect_true(any(grepl("p_growth <- style_plot(", script, fixed = TRUE)))
    expect_true(any(grepl("p_waterfall <- plot_response_waterfall(", script, fixed = TRUE)))
    expect_true(file.exists(file.path(dir, "data", "plot_data.rds")))
    expect_true(any(grepl("use the dashboard's options", script, fixed = TRUE)))
    env <- new.env(parent = globalenv())
    old <- setwd(dir); on.exit(setwd(old))
    grDevices::pdf(NULL)
    suppressMessages(suppressWarnings(utils::capture.output(
      sys.source("analysis.R", envir = env, keep.source = FALSE))))
    grDevices::dev.off()
    p <- env$p_growth
    expect_s3_class(p, "ggplot")
    expect_identical(p$labels$title, "Tumour growth")
    expect_identical(p$theme$text$size, 13)
    expect_identical(p$theme$legend.position, "bottom")
    colour <- Filter(function(sc) "colour" %in% sc$aesthetics, p$scales$scales)[[1]]
    expect_identical(colour$palette(length(cols)), cols)
    expect_s3_class(env$p_waterfall, "ggplot")
    expect_match(env$p_waterfall$labels$caption, "above 150%", fixed = TRUE)
  }))
})

test_that("F9: Drug Synergy's script redraws the Bliss plot with the dashboard's colours and title", {
  path <- system.file("sample_data", "combo_treatment_synthetic_data.csv", package = "mouseExperiment")
  rv <- shiny::reactiveValues(
    uploaded_data = read_study_csv(path), data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage", col_dose = "",
                    time_method = "from_date", col_date_for_day = "Date", reference_date = "2025-03-24",
                    date_order = "auto", volume_method = "calculate", col_length = "Length",
                    col_width = "Width", volume_formula = "modified_ellipsoid", control_group = "Control",
                    volume_units = "mm3", randomisation_unit = "mouse"))
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = rv), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "28", syn_endpoint_method = "model", syn_n_boot = 0,
                      syn_plot_title = "Bliss at day 28", syn_opt_text_size = 14,
                      run_syn_analysis = 1)
    expect_false(isTRUE(rv$syn_results$error))
    invisible(output$syn_plot_main)            # drawn, as when the tab is shown
    ins <- stats::setNames(rep(list(character(0)), 5),
                           paste0("syn_dl_group_", c("tables", "plots", "diagnostics", "data", "report")))
    do.call(session$setInputs, c(ins, list(syn_dl_group_code = "script")))
    dir <- tempfile("syn_")
    utils::unzip(output$syn_dl_zip, exdir = dir)
    script <- readLines(file.path(dir, "analysis.R"))
    expect_true(any(grepl("p_synergy <- plot_drug_synergy(", script, fixed = TRUE)))
    env <- new.env(parent = globalenv())
    old <- setwd(dir); on.exit(setwd(old))
    grDevices::pdf(NULL)
    suppressMessages(suppressWarnings(utils::capture.output(
      sys.source("analysis.R", envir = env, keep.source = FALSE))))
    grDevices::dev.off()
    expect_identical(env$p_synergy$labels$title, "Bliss at day 28")
    expect_identical(env$p_synergy$theme$text$size, 14)
    fill <- Filter(function(sc) "fill" %in% sc$aesthetics, env$p_synergy$scales$scales)[[1]]
    expect_identical(names(fill$palette(2)), c("Observed", "Expected"))
  }))
})
