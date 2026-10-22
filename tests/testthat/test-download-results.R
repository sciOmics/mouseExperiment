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
