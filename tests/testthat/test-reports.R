# Round 22 proposals F4 and F10 (CODE_REVIEW.md R22-H, 0.46.0): a report for
# every tab and for the whole study, and every package version in the ZIP.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}
all_features_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                                package = "mouseExperiment"),
                                    stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "Cage", col_weight = "Weight",
                    col_survival_event = "Survival_Censor", control_group = "Vehicle",
                    volume_units = "mm3", randomisation_unit = "mouse"))
}
download_only <- function(session, prefix, group, ids) {
  groups <- c("tables", "plots", "diagnostics", "data", "code", "report")
  ins <- stats::setNames(rep(list(character(0)), length(groups)), paste0(prefix, "_dl_group_", groups))
  ins[[paste0(prefix, "_dl_group_", group)]] <- ids
  do.call(session$setInputs, ins)
}
unzip_to <- function(z) { d <- tempfile("zip_"); utils::unzip(z, exdir = d); d }

test_that("F4: a report holds the tab's settings, tables and plots, the plots embedded", {
  item <- g("dl_item")
  items <- list(
    item("t", "tables", "Effects", "effects", "table", function() data.frame(Arm = c("A", "B"), TGI = c(10, 20)),
         "Each arm's TGI"),
    item("p", "plots", "A plot", "a_plot", "plot", function() ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                                                                             ggplot2::aes(x, y)) + ggplot2::geom_point()),
    item("d", "data", "Data", "data", "table", function() data.frame(x = 1)))
  sec <- g("report_section")("Toxicity", items, header = c("Settings of the run", "  threshold: 15"),
                              size = list(width = 400, height = 300, dpi = 72))
  expect_match(sec, "<h1>Toxicity</h1>", fixed = TRUE)
  expect_match(sec, "threshold: 15", fixed = TRUE)
  expect_match(sec, "<th>Arm</th><th>TGI</th>", fixed = TRUE)
  expect_match(sec, "data:image/png;base64,", fixed = TRUE)
  expect_false(grepl(">Data</h3>", sec, fixed = TRUE))   # the data are not in the report
  f <- tempfile(fileext = ".html")
  g("report_document")(f, "Toxicity report", sec)
  expect_match(paste(readLines(f), collapse = ""), "<title>Toxicity report</title>", fixed = TRUE)
})

test_that("F4: Survival, which had no report, offers one built from its own results", {
  rv <- all_features_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = rv), {
    session$setInputs(surv_event_source = "column")
    session$setInputs(run_surv_analysis = 1)
    expect_false(isTRUE(rv$surv_results$error))
    expect_match(html_of(output$surv_dl_choices), "HTML report", fixed = TRUE)
    download_only(session, "surv", "report", "report")
    html <- paste(readLines(file.path(unzip_to(output$surv_dl_zip), "report.html")), collapse = "")
    expect_match(html, "<h1>Survival</h1>", fixed = TRUE)
    expect_match(html, "data:image/png;base64,", fixed = TRUE)   # the Kaplan-Meier plot
    expect_match(html, "<table>", fixed = TRUE)
    # The tab is registered for the study report.
    expect_true("survival" %in% names(session$userData$dl_registry))
  }))
})

test_that("F4: the study report has the data check and every registered tab with a result", {
  item <- g("dl_item")
  d <- utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                   package = "mouseExperiment"))
  d$Volume[3] <- -d$Volume[3]
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    rv$uploaded_data <- d
    rv$data_valid <- TRUE
    session$setInputs(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage",
                      col_day = "Day", time_method = "direct_day", volume_method = "direct",
                      col_volume = "Volume")
    read <- function() paste(readLines(output$download_study_report), collapse = "")
    expect_match(read(), "No analysis has been run", fixed = TRUE)
    session$userData$dl_registry <- list(
      toxicity = list(title = function() "Toxicity",
                      items = function() list(item("t", "tables", "Window", "window", "table",
                                                   function() data.frame(Arm = "A", TGI = 40))),
                      header = function() "Settings of the run", size = function() list()),
      synergy = list(title = function() "Drug Synergy", items = function() list(),
                     header = function() character(0), size = function() list()))
    html <- read()
    expect_match(html, "<h1>Data</h1>", fixed = TRUE)
    expect_match(html, "Negative volume", fixed = TRUE)
    expect_match(html, "<h1>Toxicity</h1>", fixed = TRUE)
    expect_false(grepl("Drug Synergy", html, fixed = TRUE))   # no result, no section
  })
})

test_that("F10: the R script's ZIP holds an renv lockfile of every package version", {
  lock <- g("package_lock")()
  expect_identical(lock$Packages$brms$Version, as.character(utils::packageVersion("brms")))
  expect_identical(lock$Packages$mouseExperiment$Version,
                   utils::packageDescription("mouseExperiment")$Version)
  expect_identical(lock$R$Version, paste(R.version$major, R.version$minor, sep = "."))
  expect_false("stats" %in% names(lock$Packages))          # base packages come with R
  expect_true(all(c("lme4", "ggplot2", "survival") %in% names(lock$Packages)))
  rv <- shiny::reactiveValues(uploaded_data = NULL, data_valid = FALSE, session_id = NULL)
  suppressWarnings(shiny::testServer(g("power_server"), args = list(rv = rv), {
    session$setInputs(pw_mode = "apriori_analytic", pw_ap_es_mode = "cohens_d", pw_ap_cohens_d = "1",
                      pw_ap_n_groups = "2", pw_alpha = "0.05", pw_power = "0.80",
                      pw_ap_p_adjust = "none", pw_ap_n_comparisons = 1, pw_ap_dropout = 0,
                      pw_ap_submode = "find_n")
    session$setInputs(run_pw_analysis = 1)
    download_only(session, "pw", "code", "script")
    dir <- unzip_to(output$pw_dl_zip)
    expect_true(file.exists(file.path(dir, "packages.lock")))
    zipped <- jsonlite::read_json(file.path(dir, "packages.lock"))
    expect_identical(zipped$Packages$brms$Source, lock$Packages$brms$Source)
    expect_true(any(grepl("packages.lock", readLines(file.path(dir, "analysis.R")), fixed = TRUE)))
  }))
})
