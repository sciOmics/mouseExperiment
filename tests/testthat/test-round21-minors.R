# Regression tests for the Round 21 Minors and R20.D53 (CODE_REVIEW.md R20-D.S).
# Each test fails on v26.10.004.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
has <- function(nm) exists(nm, envir = .pkg, inherits = FALSE)

html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}
src_dir <- testthat::test_path("..", "..")
read_src <- function(f) paste(readLines(file.path(src_dir, "R", f), warn = FALSE), collapse = "\n")
demo <- function(f) utils::read.csv(system.file("sample_data", f, package = "mouseExperimentDashboard"),
                                    stringsAsFactors = FALSE)
all_features_rv <- function(...) {
  shiny::reactiveValues(
    uploaded_data = demo("all_features_synthetic_data.csv"), data_valid = TRUE, session_id = NULL,
    settings = utils::modifyList(list(
      col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
      volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
      col_weight = "Weight", control_group = "Vehicle", volume_units = "mm3",
      randomisation_unit = "mouse"), list(...)))
}
combo_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = demo("combo_treatment_synthetic_data.csv"), data_valid = TRUE,
    session_id = NULL,
    settings = list(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", time_method = "from_date", col_date_for_day = "Date",
                    reference_date = "2025-03-24", date_order = "auto",
                    volume_method = "calculate", col_length = "Length", col_width = "Width",
                    volume_formula = "modified_ellipsoid", control_group = "Control",
                    volume_units = "mm3", randomisation_unit = "mouse"))
}

# ---- R20.D20 ------------------------------------------------------------------------

test_that("R20.D20: the synergy bars are coloured, not grey", {
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = combo_rv()), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "28", syn_endpoint_method = "model",
                      syn_n_boot = 0, run_syn_analysis = 1)
    invisible(output$syn_plot_main)
    p <- current_syn_main_plot()
    fills <- unique(ggplot2::ggplot_build(p)$data[[1]]$fill)
    expect_length(fills, 2L)                       # Observed and Expected
    expect_false(any(fills %in% c("grey50", "#7F7F7F")))
  }))
})

# ---- R20.D22 / R20.D24 / R20.D25 / R20.D26 / R20.D61: code hygiene ----------------------

test_that("R20.D22/D24/D25/D26/D61: values reach JS escaped; dead code and files are gone", {
  all_src <- paste(vapply(list.files(file.path(src_dir, "R")), read_src, ""), collapse = "\n")
  expect_false(grepl("runjs\\(sprintf\\([^)]*'%s'", all_src))
  expect_false(file.exists(file.path(src_dir, "global.R")))
  expect_false(file.exists(file.path(src_dir, "app_standalone.R")))
  for (fn in c("appender_buffer", "get_log_entries", "safe_plot", "plot_export_ui",
               "extract_summary_table", "homeUI", "log_download")) {
    expect_false(has(fn), info = fn)
  }
  desc <- read.dcf(file.path(src_dir, "DESCRIPTION"))
  expect_match(desc[, "Depends"], "R \\(>= 4.1.0\\)")
  expect_false("Date" %in% colnames(desc))
  ign <- readLines(file.path(src_dir, ".Rbuildignore"))
  expect_true(all(c("^\\.waylog$", "^CODE_REVIEW\\.md$", "^docs$") %in% ign))
  # R20.D61: "NA" as a word, not any "na".
  expect_identical(g("translate_error")("could not find the column name 'x'"),
                   "could not find the column name 'x'")
  expect_match(g("translate_error")("NA/NaN/Inf in 'y'"), "missing or invalid values")
})

# ---- R20.D23: uploads ------------------------------------------------------------------

test_that("R20.D23: semicolon, tab, Latin-1 and UTF-16 files load; a shifted row is named", {
  rd <- g("read_study_csv")
  tmp <- function(lines, enc = "UTF-8") {
    f <- tempfile(fileext = ".csv")
    con <- file(f, open = "w", encoding = enc); writeLines(lines, con); close(con)
    f
  }
  semi <- rd(tmp(c("ID;Day;Volume;Treatment", "1;0;12,5;Vehicle", "2;3;20,25;Drug")))
  expect_equal(semi$Volume, c(12.5, 20.25))         # decimal comma
  tab <- rd(tmp(c("ID\tDay\tVolume", "1\t0\t12.5", "2\t3\t20")))
  expect_identical(ncol(tab), 3L)                   # was one column
  lat <- rd(tmp(c("ID,Treatment,Volume", "1,Vehicle,10", "2,Médicament,12"), enc = "latin1"))
  expect_identical(lat$Treatment[2], "Médicament")
  u16 <- tempfile(fileext = ".csv")
  bytes <- c(as.raw(c(0xFF, 0xFE)),
             iconv("ID,Volume\n1,10\n2,12\n", "UTF-8", "UTF-16LE", toRaw = TRUE)[[1]])
  writeBin(bytes, u16)
  expect_equal(rd(u16)$Volume, c(10, 12))
  expect_error(rd(tmp(c("ID,Treatment,Volume", "1,Vehicle,10", "2,Drug A, 5 mg,12"))),
               "Line 3 has 4 fields but the header has 3")
  # The shipped demos read as before.
  f <- system.file("sample_data", "all_features_synthetic_data.csv",
                   package = "mouseExperimentDashboard")
  expect_equal(rd(f), utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("NA", "")))
  expect_match(read_src("app.R"), "shiny.maxRequestSize = 30 \\* 1024\\^2")
})

# ---- R20.D29 / R20.D54: survival labels -------------------------------------------------

test_that("R20.D29/D54: the KM p-value says what it is; the survival help is correct", {
  d <- data.frame(Time = c(5, 8, 12, 14, 20, 20, 9, 20), Event = c(1, 1, 1, 0, 0, 0, 1, 0),
                  Treatment = rep(c("A", "B"), each = 4))
  expect_match(g("km_pval_label")(d), "^Log-rank, all arms \\(unadjusted\\): p ")
  s <- read_src("mod_survival.R")
  expect_false(grepl("does not include 1.0, the difference is statistically significant", s, fixed = TRUE))
  expect_false(grepl("automatically treated as events", s, fixed = TRUE))
  expect_match(s, "judge significance by P_Value", fixed = TRUE)
})

# ---- R20.D30: explainer text --------------------------------------------------------------

test_that("R20.D30: cage, random-effects and treatment-effects help match what the backend does", {
  u <- read_src("ui_modules.R")
  expect_match(u, "Cage is separable from treatment and enters as a fixed effect", fixed = TRUE)
  expect_false(grepl("'Intercept only' (default)", u, fixed = TRUE))
  expect_false(grepl("five quantile days", u, fixed = TRUE))
  expect_false(grepl("tg_variance_test", u, fixed = TRUE))
})

# ---- R20.D32 / R20.D59: Bayesian diagnostics -----------------------------------------------

test_that("R20.D59/D32: unavailable sampler diagnostics are not shown as clean", {
  panel <- function(nuts) html_of(g("bayes_diagnostics_panel")(list(nuts_diagnostics = nuts)))
  na_nuts <- data.frame(n_divergent = NA_integer_, n_max_treedepth = NA_integer_,
                        ebfmi_min = NA_real_, clean = NA)
  txt <- panel(na_nuts)
  expect_match(txt, "Sampler diagnostics unavailable")
  expect_false(grepl("diagnostics clean", txt, fixed = TRUE))
  bad <- data.frame(n_divergent = 3L, n_max_treedepth = 0L, ebfmi_min = 0.5, clean = FALSE)
  txt <- panel(bad)
  expect_match(txt, "3 divergent transition")
  expect_false(grepl("adapt_delta", txt, fixed = TRUE))
  ok <- data.frame(n_divergent = 0L, n_max_treedepth = 0L, ebfmi_min = 0.9, clean = TRUE)
  expect_match(panel(ok), "diagnostics clean")
  expect_identical(g("forest_ci_note")(data.frame(Lower_CrI = 1), "Lower_CrI"),
                   "Intervals: 95% posterior credible intervals")
})

# ---- R20.D53 / R20.D55: the growth-curve plot ------------------------------------------------

test_that("R20.D53: typing a title does not rebuild the growth plot", {
  calls <- 0L
  real <- mouseExperiment::plot_tumor_growth
  testthat::local_mocked_bindings(
    plot_tumor_growth = function(...) { calls <<- calls + 1L; real(...) },
    .package = "mouseExperiment")
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = all_features_rv()), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", run_tg_analysis = 1)
    invisible(output$tg_plot_growth)
    before <- calls
    for (t in c("M", "My", "My s", "My study")) {
      session$setInputs(tg_plot_title = t)
      invisible(output$tg_plot_growth)
    }
    session$elapse(1000)
    invisible(output$tg_plot_growth)
    expect_identical(calls, before)                 # one build per keystroke before
    expect_identical(tg_growth_plot()$labels$title, "My study")
  }))
})

test_that("R20.D55: the displayed growth plot shows necrotic markers; extrapolation says why not", {
  suppressWarnings(shiny::testServer(g("tumor_growth_server"),
                                     args = list(rv = all_features_rv(col_necrotic = "Necrotic")), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", run_tg_analysis = 1)
    p <- tg_growth_plot()
    shapes <- unlist(lapply(p$layers, function(l) l$aes_params$shape))
    expect_true(4 %in% shapes)                     # the x markers
    expect_match(paste(p$labels$caption, collapse = " "), "necrotic")
    session$setInputs(tg_extrapolation_points = 3)
    expect_match(html_of(output$tg_extrap_note), "No Event/Censor column is mapped")
  }))
})

# ---- R20.D57 / R20.D58 / R20.D60 -------------------------------------------------------------

test_that("R20.D57: the growth-rate 'geometric mean' bar is the mean of the slopes", {
  p <- read_src("mod_tumor_growth_plots.R")
  expect_match(p, "geometric_mean = mean(x), mean(x))", fixed = TRUE)
})

test_that("R20.D58: the dropped Tumor Volume covariate is explained", {
  rv <- all_features_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE, tox_covariates = "volume")
    expect_match(html_of(output$tox_covariate_note), "Tumor Volume covariate is not used")
    session$setInputs(tox_adjust_tumor = FALSE)
    expect_identical(html_of(output$tox_covariate_note), "")
  })
})

test_that("R20.D60: the Power summary uses the run's target power", {
  st <- data.frame(Effect_Size = 1, N_Groups = 2, Alpha = 0.05, Alpha_Per_Comparison = 0.05,
                   N_Comparisons = 1, Target_Power = 0.95, Required_N = 27, Enroll_N = 27,
                   Dropout_Rate = 0)
  res <- list(scenario_table = st, power_curve_data = NULL, sensitivity_table = NULL,
              mode = "apriori_analytic", ap_submode = "find_n", target_power = 0.95,
              p_adjust_method = "bonferroni", n_comparisons = 1, dropout_rate = 0)
  rv <- shiny::reactiveValues(uploaded_data = NULL, data_valid = TRUE, settings = list(),
                              ap_results = res)
  shiny::testServer(g("power_server"), args = list(rv = rv), {
    txt <- html_of(output$pw_ap_summary_ui)
    expect_match(txt, "target power = 95%")         # always said 80 %
    expect_match(txt, "find N for target power")
  })
})
