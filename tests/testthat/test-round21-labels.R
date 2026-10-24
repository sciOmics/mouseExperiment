# Regression tests for the Round 21 labelling Majors (CODE_REVIEW.md R20-D.Q):
# Q-Q plots (R20.D42), forest-plot intervals (R20.D43), the Comparisons tab
# (R20.D44), the Data Summary scale (R20.D45), the report's settings (R20.D46),
# the Power seed SD (R20.D51) and analytics payloads (R20.D12). Each test fails
# on v26.10.002.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}

all_features_path <- function() {
  system.file("sample_data", c("all_features_synthetic_data.csv",
                               "master_synthetic_data.csv"),
              package = "mouseExperiment")[1]
}
all_features_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(all_features_path(), stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "Cage",
                    col_weight = "Weight", control_group = "Vehicle",
                    volume_units = "mm3", randomisation_unit = "mouse"))
}
run_tg <- function(session, ...) {
  inputs <- utils::modifyList(list(tg_model_type = "lme4", tg_transform = "log",
                                   tg_comparison_type = "reference",
                                   tg_p_adjust = "dunnett"), list(...))
  do.call(session$setInputs, inputs)
  session$setInputs(run_tg_analysis = 1)
}

# ---- R20.D42 ---------------------------------------------------------------------

test_that("R20.D42: the residual Q-Q is the backend's standardised plot", {
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    run_tg(session)
    qd <- ggplot2::ggplot_build(rv$tg_results$diag_qq_plot)$data[[1]]
    # The dashboard's copy plotted raw residuals (SD about 0.09) against y = x.
    expect_equal(stats::sd(qd$y), 1, tolerance = 1e-6)
  }))
})

# ---- R20.D43 / R20.D44 -------------------------------------------------------------

test_that("R20.D43: the forest plot uses the adjusted intervals and says so", {
  pc <- data.frame(contrast = "A - Vehicle", estimate = -0.5, SE = 0.1, df = 20,
                   lower.CL = -0.7, upper.CL = -0.3,
                   CI_Adjusted_Lower = -0.75, CI_Adjusted_Upper = -0.25,
                   CI_Adjustment = "mvt method for 5 estimates")
  expect_match(g("forest_ci_note")(pc, "CI_Adjusted_Lower"), "adjusted .*mvt.*like the p-values")
  expect_match(g("forest_ci_note")(transform(pc, CI_Adjustment = "none (95% per comparison)"),
                                   "CI_Adjusted_Lower"), "no correction applied")
  expect_match(g("forest_ci_note")(pc[, 1:6], NA_character_), "not adjusted")

  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    run_tg(session)
    rc <- rv$tg_results$reference_comparisons_dunnett
    expect_true(all(c("CI_Adjusted_Lower", "CI_Adjusted_Upper") %in% names(rc)))
    # The displayed forest plot (R20.D61 removed an unused second copy).
    p <- tg_forest_plot()
    expect_s3_class(p, "ggplot")
    expect_match(p$labels$caption, "adjusted")
  }))
})

test_that("R20.D44: the Comparisons tab shows the run's family, titled, with the run's correction", {
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    run_tg(session, tg_p_adjust = "none")
    help <- html_of(output$tg_comparisons_tab_info)
    expect_match(help, "not adjusted for multiple comparisons")
    expect_false(grepl("Dunnett correction for", help, fixed = TRUE))
    tbl <- html_of(output$tg_comparisons_table)
    expect_match(tbl, "Comparisons: vs Reference Group \\(unadjusted\\)")
    # Switching the dropdown without re-running keeps this run's rows and says so.
    session$setInputs(tg_comparison_type = "tukey")
    tbl <- html_of(output$tg_comparisons_table)
    expect_match(tbl, "vs Reference Group")
    expect_match(tbl, "run the analysis again")
    expect_false(grepl("not available for tukey", tbl, fixed = TRUE))
  }))
})

# ---- R20.D45 -------------------------------------------------------------------------

test_that("R20.D45: the Data Summary is on the volume scale, with units in the headers", {
  ds <- data.frame(Treatment = "Vehicle", Day = 20, mean_volume = 3991.6,
                   sd_volume = NA, n = 1L, geometric_mean_volume = 3991.6, sem_volume = NA)
  shown <- g("data_summary_display")(ds, "mm\u00b3")
  expect_identical(names(shown), c("Treatment", "Day", "N", "Mean (mm\u00b3)", "SD", "SEM",
                                   "Geometric mean (mm\u00b3)"))
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    run_tg(session)
    d <- rv$tg_results$data_summary
    raw <- rv$uploaded_data
    v0 <- raw$Volume[raw$Treatment == "Vehicle" & raw$Day == 0]
    # v26.10.002 showed log-scale means here (4.3 for 75.9 mm3).
    expect_equal(d$mean_volume[d$Treatment == "Vehicle" & d$Day == 0], mean(v0, na.rm = TRUE))
    expect_match(html_of(output$tg_data_summary_table), "Mean \\(mm\u00b3\\)")
  }))
})

# ---- R20.D46 -------------------------------------------------------------------------

test_that("R20.D46: the report records the run's settings, not the sidebar's", {
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    run_tg(session)
    # Change the sidebar without running again.
    session$setInputs(tg_model_type = "auc", tg_transform = "sqrt")
    html <- paste(readLines(output$download_tg_report, warn = FALSE), collapse = "\n")
    # 0.46.0 (F4): the settings are the README header's, in the report's
    # opening block.
    expect_match(html, "  model_type: lme4\n", fixed = TRUE)
    expect_match(html, "  transform: log\n", fixed = TRUE)
    expect_false(grepl("transform: sqrt", html, fixed = TRUE))
    expect_match(html, "  reference_group: Vehicle\n", fixed = TRUE)
  }))
})

# ---- R20.D51 -------------------------------------------------------------------------

test_that("R20.D51: the seed SD is taken at one endpoint day, not pooled over every day", {
  d <- utils::read.csv(all_features_path(), stringsAsFactors = FALSE)
  s <- g("seed_sd_at_day")(d, "Treatment", "Day", "Volume", "ID", "Cage", day = 14)
  expect_equal(s$day, 14)
  # Every observation pooled (the old seed) was 617.8 on this demo.
  old <- local({
    sds <- tapply(d$Volume, d$Treatment, stats::sd, na.rm = TRUE)
    ns  <- tapply(d$Volume, d$Treatment, function(x) sum(!is.na(x)))
    sqrt(sum((ns - 1) * sds^2) / sum(ns - 1))
  })
  expect_lt(s$sd, old / 2)
  # By hand: each animal's last volume on or before day 14.
  dd <- d[d$Day <= 14 & is.finite(d$Volume), ]
  dd <- dd[order(dd$Treatment, dd$ID, dd$Cage, dd$Day), ]
  last <- dd[!duplicated(paste(dd$Treatment, dd$ID, dd$Cage), fromLast = TRUE), ]
  sds <- tapply(last$Volume, last$Treatment, stats::sd)
  ns  <- tapply(last$Volume, last$Treatment, length)
  expect_equal(s$sd, sqrt(sum((ns - 1) * sds^2) / sum(ns - 1)))
  # No day: the last day every arm still has 3 animals measured.
  auto <- g("seed_sd_at_day")(d, "Treatment", "Day", "Volume", "ID", "Cage")
  n3 <- tapply(d$ID[d$Day == auto$day & is.finite(d$Volume)],
               d$Treatment[d$Day == auto$day & is.finite(d$Volume)], length)
  expect_true(all(n3 >= 3))
})

# ---- R20.D12 -------------------------------------------------------------------------

test_that("R20.D12: no study string reaches the analytics events", {
  js <- paste(readLines(system.file("www", "js", "umami-events.js",
                                    package = "mouseExperiment")), collapse = "\n")
  for (field in c("reference_group:", "error_message:", "treatment_col:",
                  "time_col:", "volume_col:")) {
    expect_false(grepl(field, js, fixed = TRUE), info = field)
  }
  expect_match(js, "properties.error = true;", fixed = TRUE)
  src <- dash_all_sources()
  expect_false(grepl("trackError\\([^)]*%s", src))
  expect_false(grepl("trackColumnMapping\\('%s'", src))
  # Staging traffic is not recorded under the production website id.
  expect_match(src, "`data-domains` = \"mouseexperiment.com,www.mouseexperiment.com\"",
               fixed = TRUE)
})
