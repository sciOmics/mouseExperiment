# Sidebar explanations moved into help_tip() tooltips (26.10.007, CODE_REVIEW.md
# R20-D.V), and the necrotic-handling selector that appeared only after the
# first run. The tests fail on 26.10.006.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}
all_features_rv <- function(...) {
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                                package = "mouseExperiment"),
                                    stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = utils::modifyList(list(
      col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
      volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
      col_weight = "Weight", control_group = "Vehicle", volume_units = "mm3",
      randomisation_unit = "mouse"), list(...)))
}
run_tg <- function(session, n = 1) {
  session$setInputs(tg_model_type = "lme4", tg_transform = "log",
                    tg_comparison_type = "reference", tg_p_adjust = "dunnett")
  session$setInputs(run_tg_analysis = n)
}

test_that("26.10.007: the Tumour Growth sidebar explains its options in tooltips", {
  html <- paste(as.character(g("tumorGrowthUI")("tg")), collapse = "\n")
  side <- substr(html, 1L, regexpr("run_tg_analysis", html, fixed = TRUE))  # up to Run
  for (lbl in c("Statistical Model", "Post-hoc Comparisons", "Run randomisation test",
                "Multiple-comparison correction", "Data Transformation",
                "Cage Effect Handling", "Random Effects", "Prior Strength")) {
    expect_match(side, paste0(lbl, "[^<]*<span class=\"help-tip"), info = lbl)
  }
  # R23.9: Extrapolation points is a plot control now, with its tooltip.
  expect_false(grepl("Extrapolation Points", side, fixed = TRUE))
  expect_match(html, "Extrapolation points[^<]*<span class=\"help-tip")
  # No helpText() paragraphs are left between the controls.
  expect_false(grepl("help-block", side, fixed = TRUE))
  # Internal review references are not user-facing text.
  expect_false(grepl("CODE_REVIEW", html, fixed = TRUE))
})

test_that("R20-D.V: the necrotic selector shows once mapped, and the first run uses it", {
  rv <- all_features_rv(col_necrotic = "Necrotic")
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    # Before any run. 26.10.006 drew it only after a run had set rv$necrotic_col.
    ui <- html_of(output$necrotic_handling_ui)
    expect_match(ui, "tg_necrotic_handling", fixed = TRUE)
    expect_match(ui, '<option value="exclude" selected>', fixed = TRUE)

    # A first run, with the selector's input never set. 26.10.006 sent "none".
    run_tg(session)
    first <- rv$tg_results
    expect_identical(first$run_settings$necrotic_handling, "exclude")
    session$setInputs(tg_necrotic_handling = "exclude")
    run_tg(session, n = 2)
    expect_identical(rv$tg_results$run_settings$necrotic_handling, "exclude")
    expect_identical(stats::nobs(first$model), stats::nobs(rv$tg_results$model))
  }))
})

test_that("R20-D.V: without a necrotic column there is no selector and nothing to handle", {
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    expect_identical(html_of(output$necrotic_handling_ui), "")
    run_tg(session)
    expect_identical(rv$tg_results$run_settings$necrotic_handling, "none")
  }))
})

test_that("26.10.009: the other analysis tabs explain their options in tooltips too", {
  side_of <- function(ui_fn, id) {
    html <- paste(as.character(g(ui_fn)(id)), collapse = "\n")
    run  <- regexpr('id="[a-z]+-run_[a-z]+_analysis"', html)
    substr(html, 1L, run)                              # the sidebar, up to Run
  }
  expect_tip <- function(side, lbl)
    expect_match(side, paste0(lbl, "[^<]*(</strong>)?\\s*<span class=\"help-tip"), info = lbl)
  sv <- side_of("survivalUI", "surv")
  for (l in c("Statistical Model", "Endpoint")) expect_tip(sv, l)
  tx <- side_of("toxicityUI", "tox")
  # R23.19: no Estimation method choice.
  for (l in c("Statistical Model", "Endpoint estimand",
              "Interval draws", "Tolerability threshold")) expect_tip(tx, l)
  dr <- side_of("doseResponseUI", "dr")
  expect_tip(dr, "Pairwise Comparison Type")
  pw <- side_of("powerAnalysisUI", "pw")
  for (l in c("Analysis Mode", "Multiplicity correction", "Expected attrition",
              "Baseline tumor volume \\(SD\\)", "Control log-growth rate",
              "Euthanasia volume limit", "Number of simulations")) expect_tip(pw, l)
  # No helpText() paragraphs or small-print notes left between the controls.
  for (side in list(sv, tx, dr, pw)) expect_false(grepl("help-block", side, fixed = TRUE))
  expect_false(grepl("Benchmarks: small", pw, fixed = TRUE) &&
               !grepl("help-tip[^>]*title=\"Benchmarks", pw))
})
