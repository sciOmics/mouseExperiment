# Backend 0.30.0 (package CODE_REVIEW.md R20-T) changed what four results hold:
# a fixed zero floor that needs the volume unit (R20.32), an ANOVA led by the
# growth-rate test (R20.31), body-weight comparisons of the rate of weight change
# (R20.44), and per-arm simulated power (R20.33). These tests check that the
# dashboard passes what the backend needs and shows the results as they are.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}

# The all-features demo, under either of its file names.
all_features_rv <- function(units = "mm3") {
  path <- system.file("sample_data", c("all_features_synthetic_data.csv",
                                       "master_synthetic_data.csv"),
                      package = "mouseExperiment")[1]
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(path, stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "Cage",
                    col_weight = "Weight", control_group = "Vehicle",
                    volume_units = units, randomisation_unit = "mouse"))
}

# ---- R20.32: the unit reaches the zero floor ------------------------------------

test_that("R20.32: the declared unit is the backend's volume_units; undeclared is mm3", {
  unit_of <- function(u) g("backend_volume_units")(list(settings = list(volume_units = u)))
  expect_identical(unit_of("cm3"), "cm3")
  expect_identical(unit_of("mm3"), "mm3")
  expect_identical(unit_of(""), "mm3")
  expect_identical(unit_of(NULL), "mm3")
})

test_that("R20.32: the summary states the floor on the log scale, and only there", {
  line <- g("zero_floor_line")
  zf <- list(value = 2, units = "mm3", n_floored = 3L)
  expect_identical(line(list(zero_floor = zf, transform_used = "log")),
                   "  Zero volumes: set to 2 mm³ before the log (3 values)")
  expect_identical(line(list(zero_floor = zf, transform_used = "sqrt")), character(0))
  expect_identical(line(list(transform_used = "log")), character(0))
})

test_that("R20.32: tumour growth runs with the declared cm3 and reports its floor", {
  rv <- all_features_rv("cm3")
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", run_tg_analysis = 1)
    expect_false(isTRUE(rv$tg_results$error))
    expect_identical(rv$tg_results$zero_floor$units, "cm3")
    expect_equal(rv$tg_results$zero_floor$value, 0.002)
    expect_match(paste(output$tg_summary, collapse = " "),
                 "Zero volumes: set to 0.002 cm³ before the log")
  }))
})

# ---- R20.31: the ANOVA help matches the table -------------------------------------

test_that("R20.31: the ANOVA help describes the growth-rate test first, not a Type III table", {
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", run_tg_analysis = 1)
    a <- rv$tg_results$anova
    expect_true("Tests" %in% names(a))
    expect_match(rownames(a)[1], ":")                    # treatment x time leads
    help <- html_of(output$tg_anova_tab_info)
    expect_match(help, "growth rates differ")
    expect_false(grepl("Type III", help, fixed = TRUE))
  }))
})

# ---- R20.44: body weight shows rates and their comparisons ------------------------

test_that("R20.44: body weight shows each arm's rate (g/day) and the rate comparisons", {
  rv <- all_features_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE,
                      tox_tumor_density = 1, tox_covariates = "volume",
                      tox_estimation = "REML", tox_threshold = 20,
                      tox_baseline_day = "first", run_tox_analysis = 1)
    mm <- rv$tox_results$mixed_model
    expect_false(isTRUE(mm$error))
    expect_true("Weight_Change_Per_Day" %in% names(mm$emmeans_table))
    # DT renders server-side: the output carries the header, the CSV the rows.
    expect_match(html_of(output$tox_rate_comparisons), "contrast")
    expect_match(html_of(output$tox_bw_mixed_model_tab_info), "g/day")
    csv <- utils::read.csv(output$download_tox_rate_comparisons_csv)
    expect_equal(nrow(csv), length(unique(rv$uploaded_data$Treatment)) - 1L)
    expect_true(all(grepl("Vehicle", csv$contrast)))
  }))
})

# ---- R20.33: simulated power is per arm, with Wilson limits -----------------------

test_that("R20.33: the Power summary says what the simulated power refers to", {
  st <- data.frame(N_Per_Group = c(5L, 8L), N_Groups = 3L, Alpha = 0.05,
                   Alpha_Per_Comparison = 0.025, Power = c(0.6, 0.95),
                   Power_Lower = c(0.45, 0.83), Power_Upper = c(0.74, 0.99),
                   SE_Power = c(0.08, 0.03), Power_Any = c(0.8, 1),
                   Power_All = c(0.4, 0.9), Power_Omnibus = c(0.8, 1),
                   N_Significant = c(24L, 38L), N_Simulations = 20L,
                   N_Failed = c(0L, 1L))
  txt <- html_of(g("lmm_power_note")(st))
  expect_match(txt, "a given treated arm is significant")
  expect_match(txt, "0.025 for 2 treated arms")
  expect_match(txt, "1 simulated study could not be fitted")
  # A table from an older backend has no per-comparison alpha: no note.
  expect_null(g("lmm_power_note")(st[, c("N_Per_Group", "Power", "Alpha")]))

  res <- list(scenario_table = st, power_curve_data = st, sensitivity_table = NULL,
              mode = "apriori_lmm", ap_submode = "find_power",
              lmm_params = list(dropout_limit = Inf), lmm_attrition = NULL)
  rv <- shiny::reactiveValues(uploaded_data = NULL, data_valid = TRUE,
                              settings = list(), ap_results = res)
  shiny::testServer(g("power_server"), args = list(rv = rv), {
    expect_match(html_of(output$pw_ap_summary_ui), "Power_All")
    expect_true(nzchar(output$pw_ap_plot_curves$src))   # the Wilson band draws
  })
})
