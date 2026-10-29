# The maintainer's manual review of the dashboard, the small items
# (CODE_REVIEW.md Round 23, R23-B): tab order, the Bayesian tabs, the Stan
# backend menu, the extrapolation control, the drill-down, wording, Synergy's
# legend, Toxicity's covariates and estimation, and interval labels.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}
ui_html <- function(f) paste(as.character(g(f)("m")), collapse = "\n")
# The values of the first tabset's tabs, in order.
tab_values <- function(html) {
  ul <- regmatches(html, regexpr('(?s)<ul class="nav nav-tabs[^>]*>.*?</ul>', html, perl = TRUE))
  sub('data-value="', "", regmatches(ul, gregexpr('data-value="[^"]*', ul))[[1]], fixed = TRUE)
}
tg_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = master_synthetic_data, data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
}

# ---- R23.13: one tab order -------------------------------------------------------------------

test_that("R23.13: every analysis tab runs Info, Summary, results, Plots, Diagnostics, Download", {
  for (f in c("tumorGrowthUI", "survivalUI", "toxicityUI", "synergyUI", "doseResponseUI")) {
    tabs <- tab_values(ui_html(f))
    expect_identical(tabs[1:2], c("Info", "Summary"), info = f)
    expect_identical(utils::tail(tabs, 1), "__download_results__", info = f)
    expect_identical(utils::tail(tabs, 2)[1], "Diagnostics", info = f)
    # Plots comes after every results tab: Diagnostics, and Survival's
    # Prior vs Posterior (a plot shown for the Bayesian model), may follow.
    after <- tabs[seq(match("Plots", tabs) + 1L, length(tabs))]
    expect_true(all(after %in% c("Prior vs Posterior", "Diagnostics", "__download_results__")),
                info = paste(f, paste(tabs, collapse = ", ")))
  }
  pw <- tab_values(ui_html("powerAnalysisUI"))
  expect_identical(pw[1], "Info")
  expect_identical(utils::tail(pw, 1), "__download_results__")
})

# ---- R23.4, R23.5: the Bayesian tabs ------------------------------------------------------

test_that("R23.4/R23.5: a Bayesian fit swaps ANOVA for Posterior Summary; no tab after Download", {
  tg <- ui_html("tumorGrowthUI")
  tabs <- tab_values(tg)
  expect_identical(tabs[3:4], c("ANOVA", "Posterior Summary"))
  expect_false(any(grepl("Posterior Check", c(tabs, tab_values(ui_html("survivalUI"))))))
  # The posterior predictive check is in the Bayesian diagnostics.
  diag <- substr(tg, regexpr('data-value="Diagnostics"', tg), nchar(tg))
  expect_match(diag, "m-tg_pp_check_plot", fixed = TRUE)
  expect_match(diag, "Posterior predictive check", fixed = TRUE)
  sv <- ui_html("survivalUI")
  sdiag <- substr(sv, regexpr('data-value="Diagnostics"', sv), nchar(sv))
  expect_match(sdiag, "m-surv_pp_check_plot", fixed = TRUE)
  # The server hides and shows; it no longer appends tabs.
  src <- dash_source("R", "mod_tumor_growth.R")
  expect_match(src, 'hideTab(inputId = "tg_results_tabs", target = "ANOVA")', fixed = TRUE)
  expect_false(grepl("appendTab", src, fixed = TRUE))
  expect_false(grepl("appendTab", dash_source("R", "mod_survival.R"), fixed = TRUE))
  # Each column of the posterior summary is explained beside it (R23.7).
  help <- html_of(g("posterior_summary_help")("log tumour volume"))
  for (col in c("Estimate", "Est.Error", "Lower_95_CrI", "Rhat", "Bulk_ESS", "Tail_ESS"))
    expect_match(help, col, fixed = TRUE)
})

# ---- R23.3, R23.10, R23.20: no backend menu, no compile notes ---------------------------------

test_that("R23.3/R23.10/R23.20: no Stan backend menu; compiling is explained in Help only", {
  for (f in c("tumorGrowthUI", "survivalUI")) {
    h <- ui_html(f)
    expect_false(grepl("bayes_backend", h, fixed = TRUE), info = f)
    expect_false(grepl("cmdstanr", h, fixed = TRUE), info = f)
    expect_false(grepl("The first fit of a model compiles it", h, fixed = TRUE), info = f)
  }
  expect_match(dash_source("R", "mod_tumor_growth.R"), 'backend           = "rstan"', fixed = TRUE)
  expect_match(dash_source("R", "mod_survival.R"), 'backend        = "rstan"', fixed = TRUE)
  for (f in c("tg_resume.R", "surv_resume.R"))
    expect_false(grepl("compiles it", dash_source("R", f), fixed = TRUE), info = f)
  expect_match(paste(as.character(g("helpUI")()), collapse = ""), "first fit of a model compiles it",
               fixed = TRUE)
})

# ---- R23.9: the extrapolation control is a plot control ----------------------------------------

test_that("R23.9: Extrapolation points sits under 'Show extrapolated curves/points'", {
  h <- ui_html("tumorGrowthUI")
  plots_at <- regexpr('data-value="Plots"', h)
  ep_at <- regexpr('id="m-tg_extrapolation_points"', h)
  expect_gt(ep_at, plots_at)
  expect_gt(ep_at, regexpr('id="m-tg_opt_show_extrap"', h))
  expect_false(grepl("Extrapolation Points (plots only)", h, fixed = TRUE))
  # It still drives the plot and its note.
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = tg_rv()), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_p_adjust = "holm",
                      tg_extrapolation_points = 2, run_tg_analysis = 1)
    expect_match(html_of(output$tg_extrap_note), "No Event/Censor column is mapped")
  }))
})

# ---- R23.8: several mice in the drill-down -------------------------------------------------------

test_that("R23.8: the drill-down draws every selected animal, by arm, each with its rate", {
  r <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    master_synthetic_data, reference_group = "Vehicle", cage_column = "Cage", plots = FALSE)))
  rates <- as.data.frame(r$growth_rates)
  one_arm <- which(rates$Treatment == rates$Treatment[1])[1:3]
  two_arms <- c(one_arm[1], which(rates$Treatment != rates$Treatment[1])[1])
  draw <- function(rows, ...) g("drilldown_plot")(master_synthetic_data, rates, rows, "ID", "Treatment",
                                                   "Day", "Volume", "Cage", ...)
  p <- draw(one_arm)
  animals <- p$layers[[2]]$data
  expect_identical(nlevels(animals$Animal), 3L)
  expect_true(all(grepl("/day$", levels(animals$Animal))))
  expect_false(inherits(p$facet, "FacetWrap"))
  # Cage-mates that share an ear tag stay apart (R20.D57): each has its own rows.
  per <- table(animals$Animal)
  for (k in seq_along(one_arm))
    expect_identical(as.integer(per[[k]]),
                     nrow(g("drilldown_rows")(master_synthetic_data, rates[one_arm[k], ], "ID",
                                              "Treatment", "Cage")))
  expect_true(inherits(draw(two_arms)$facet, "FacetWrap"))
  single <- draw(one_arm[1])
  expect_match(single$labels$subtitle, "Growth rate", fixed = TRUE)
  capped <- draw(seq_len(nrow(rates)), max_animals = 4L)
  expect_identical(nlevels(capped$layers[[2]]$data$Animal), 4L)
  expect_match(capped$labels$subtitle, "more not drawn", fixed = TRUE)
  expect_null(draw(integer(0)))
  # The table lets several rows be selected.
  expect_match(dash_source("R", "mod_tumor_growth.R"),
               'selection = list(mode = "multiple", target = "row")', fixed = TRUE)
})

# ---- R23.16, R23.17: Synergy ------------------------------------------------------------------

test_that("R23.16/R23.17: no simulation figures in the tooltip; the shading has a legend", {
  src <- dash_source("R", "mod_synergy.R")
  expect_false(grepl("23 points", src, fixed = TRUE))
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
                      syn_time_point = "all", syn_endpoint_method = "model", syn_n_boot = 0,
                      syn_opt_show_shading = TRUE, run_syn_analysis = 1)
    expect_true(isTRUE(rv$syn_over_time))
    invisible(output$syn_plot_trend)
    p <- current_syn_trend_plot()
    fill <- Filter(function(sc) "fill" %in% sc$aesthetics, p$scales$scales)
    expect_length(fill, 1L)
    expect_identical(fill[[1]]$breaks,
                     c("Synergy (combination above Bliss)", "Antagonism (combination below Bliss)"))
    expect_identical(ggplot2::get_guide_data(p, "fill")$.label,
                     c("Synergy (combination above Bliss)", "Antagonism (combination below Bliss)"))
  }))
})

# ---- R23.18, R23.19: Toxicity ------------------------------------------------------------------

test_that("R23.18/R23.19: no covariate by default, no Tumor Volume while mass is subtracted, no ML", {
  h <- ui_html("toxicityUI")
  expect_false(grepl('id="m-tox_estimation"', h, fixed = TRUE))
  covs <- substr(h, regexpr('id="m-tox_covariates"', h), regexpr('id="m-tox_covariate_note"', h))
  expect_gt(nchar(covs), 100)
  expect_false(grepl("checked", covs, fixed = TRUE))
  rv <- shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                                package = "mouseExperiment"), stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_weight = "Weight", control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
  shiny::isolate(g("ensure_preprocessed")(rv))
  shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE)
    expect_match(html_of(output$tox_covariate_note), "not offered as a covariate")
    session$setInputs(tox_adjust_tumor = FALSE)
    expect_identical(html_of(output$tox_covariate_note), "")
  })
  src <- dash_source("R", "mod_toxicity.R")
  expect_match(src, 'if (isTRUE(input$tox_adjust_tumor)) available <- setdiff(available, "volume")',
               fixed = TRUE)
  expect_match(src, 'est_method   <- "REML"', fixed = TRUE)
})

# ---- R23.23: intervals say their level ----------------------------------------------------------

test_that("R23.23: interval columns are shown with their level", {
  lab <- g("interval_col_labels")
  expect_identical(lab(c("Contrast", "lower.CL", "upper.CL", "CI_Lower", "ci_upper", "asymp.LCL")),
                   c("Contrast", "Lower 95% CI", "Upper 95% CI", "Lower 95% CI", "Upper 95% CI",
                     "Lower 95% CI"))
  expect_identical(lab(c("CI_Adjusted_Lower", "Lower_CrI", "lower.HPD")),
                   c("Lower 95% CI (adjusted)", "Lower 95% CrI", "Lower 95% HPD"))
  expect_identical(lab(c("lower.CL", "Lower_95_CrI", "GR 95% CI Lower Bound"), 0.9),
                   c("Lower 90% CI", "Lower_95_CrI", "GR 95% CI Lower Bound"))
  d <- data.frame(contrast = c("a", "b"), estimate = 1:2, lower.CL = 0:1, upper.CL = 2:3)
  for (rn in c(FALSE, TRUE)) {
    head <- g("dt_table")(d, rownames = rn)$x$container
    expect_match(head, "<th>Lower 95% CI</th>", fixed = TRUE)
    expect_false(grepl("lower.CL", head, fixed = TRUE))
  }
  attr(d, "conf_level") <- 0.9
  expect_match(g("dt_table")(d)$x$container, "<th>Upper 90% CI</th>", fixed = TRUE)
  # A table given its own column labels keeps them.
  expect_match(g("dt_table")(d, colnames = c("A", "B", "C", "D"))$x$container, "<th>C</th>", fixed = TRUE)
  # Every dashboard table goes through dt_table(): the one DT::datatable()
  # call left is dt_table()'s own (and its documentation's).
  n_calls <- function(x) lengths(regmatches(x, gregexpr("DT::datatable(", x, fixed = TRUE)))
  expect_identical(n_calls(dash_all_sources()), n_calls(dash_source("R", "helpers_format.R")))
})
