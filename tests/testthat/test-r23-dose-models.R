# Which dose-response curve: linear, Emax, 3PL and 4PL compared by AIC
# (CODE_REVIEW.md Round 23, R23.26, R23-J).

skip_if_not_installed("drc")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

sim_dose <- function(doses, f, sd = 8, n = 8, seed = 1) {
  set.seed(seed)
  d <- data.frame(Dose = rep(doses, each = n))
  d$Volume <- f(d$Dose) + stats::rnorm(nrow(d), 0, sd)
  d
}
hill <- function(x) 20 + (100 - 20) / (1 + (x / 10)^1.5)

test_that("R23.26: lm and drc log-likelihoods are on the same scale, so their AIC compares", {
  d <- sim_dose(c(0.01, 1, 3, 10, 30), hill)
  gauss <- function(r) { n <- length(r); -n / 2 * (log(2 * pi * sum(r^2) / n) + 1) }
  m4 <- drc::drm(Volume ~ Dose, data = d, fct = drc::LL.4())
  lin <- stats::lm(Volume ~ Dose, data = d)
  expect_equal(as.numeric(stats::logLik(m4)), gauss(stats::residuals(m4)))
  expect_equal(as.numeric(stats::logLik(lin)), gauss(stats::residuals(lin)))
  expect_identical(attr(stats::logLik(m4), "df"), 5)     # 4 curve parameters + the variance
  expect_identical(attr(stats::logLik(lin), "df"), 3)
})

test_that("R23.26: the candidates, the 4PL only with 5 doses, and the simplest within 2 AIC", {
  sel <- g("dr_model_selection")
  six <- sim_dose(c(0.01, 1, 3, 10, 30, 100), hill)
  s6 <- sel(six, "Dose", "Volume", n_dose_levels = 6L)
  expect_identical(s6$table$Model, c("Linear", "Emax", "3PL", "4PL"))
  expect_identical(s6$table$Parameters, c(2L, 3L, 3L, 4L))
  expect_equal(min(s6$table$Delta_AIC), 0)
  # The curve's ED50 is near the truth (10); the linear model has none.
  expect_true(is.na(s6$table$ED50[1]))
  expect_lt(abs(log(s6$table$ED50[s6$table$Model == "4PL"] / 10)), log(2))
  # The rule: the fewest parameters among the models within 2 AIC of the best.
  near <- s6$table[s6$table$Delta_AIC < 2, ]
  expect_identical(s6$selected, near$Model[order(near$Parameters, near$AIC)][1])
  expect_identical(sum(s6$table$Selected), 1L)
  expect_identical(nrow(s6$curve), 200L)
  expect_true(all(is.finite(s6$curve$Fitted)))
  # A sigmoid truth: a curve beats the straight line by far.
  expect_gt(s6$table$Delta_AIC[1], 10)
  # Four dose levels: no 4PL.
  four <- sim_dose(c(0.01, 3, 10, 30), hill)
  expect_false("4PL" %in% sel(four, "Dose", "Volume", n_dose_levels = 4L)$table$Model)
  # A straight-line truth: the simplest model, the linear one, is selected.
  lin <- sim_dose(c(0.01, 1, 3, 10, 30, 100), function(x) 100 - 0.5 * x, sd = 5)
  expect_identical(sel(lin, "Dose", "Volume", n_dose_levels = 6L)$selected, "Linear")
})

test_that("R23.26: dose_response_statistics() reports the comparison and the selected curve", {
  dr <- suppressWarnings(suppressMessages(dose_response_statistics(
    master_synthetic_data, treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
    reference_group = "Vehicle", cage_column = "Cage", verbose = FALSE)))
  st <- dr$statistics
  expect_true(is.data.frame(st$model_selection))
  expect_true(st$selected_model %in% st$model_selection$Model)
  expect_false("4PL" %in% st$model_selection$Model)      # 4 dose levels on this demo
  expect_identical(names(st$selected_curve), c("Dose", "Fitted"))
  expect_match(st$aic_comparison_note, "can be compared", fixed = TRUE)
})

test_that("R23.26: the dashboard shows the comparison, draws the selected curve, and exports both", {
  skip_if_not_installed("shiny", minimum_version = "1.5.0")
  rv <- shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                                package = "mouseExperiment"), stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_weight = "Weight", col_dose = "Dose", col_survival_event = "Survival_Censor",
                    control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse"))
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    session$setInputs(dr_treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
                      dr_time_point = "default", run_dr_analysis = 1)
    st <- rv$dr_results$statistics
    expect_true(is.data.frame(st$model_selection))
    expect_match(paste(as.character(output$dr_model_selection), collapse = ""), "<th>Delta_AIC", fixed = TRUE)
    expect_match(output$dr_summary, "Which dose-response curve", fixed = TRUE)
    # The default trend line is the selected curve.
    p <- dr_scatter_plot_r()
    lines <- Filter(function(l) inherits(l$geom, "GeomLine"), p$layers)
    expect_length(lines, 1L)
    expect_identical(lines[[1]]$data, st$selected_curve)
    session$setInputs(dr_trend_method = "linear")
    expect_length(Filter(function(l) inherits(l$geom, "GeomSmooth"), dr_scatter_plot_r()$layers), 1L)
    groups <- c("tables", "plots", "diagnostics", "data", "report")
    ins <- stats::setNames(rep(list(character(0)), length(groups)), paste0("dr_dl_group_", groups))
    ins$dr_dl_group_tables <- "model_selection"
    ins$dr_dl_group_code <- "script"
    do.call(session$setInputs, ins)
    dir <- tempfile("dr_")
    utils::unzip(output$dr_dl_zip, exdir = dir)
    expect_true(file.exists(file.path(dir, "tables", "dose_response_models.csv")))
    expect_true(any(grepl("print(result$statistics$model_selection)",
                          readLines(file.path(dir, "analysis.R")), fixed = TRUE)))
  }))
})
