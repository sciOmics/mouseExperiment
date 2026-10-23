# The R script in Download Results repeats the run with the package
# (CODE_REVIEW.md R20-Z). Each test runs Tumor Growth through the tab, downloads
# the script with the data it reads, runs the script in a clean environment,
# and requires the script's tables to equal the dashboard's. That round trip is
# what keeps the script honest as the module changes.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

tg_rv <- function(...) {
  settings <- utils::modifyList(
    list(col_id = "ID", col_day = "Day", time_method = "day",
         col_volume = "Volume", volume_method = "precalculated",
         col_treatment = "Treatment", col_cage = "Cage", col_weight = "Weight",
         control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse"),
    list(...))
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                                package = "mouseExperiment"),
                                    stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL, settings = settings)
}

# Download only the script (the data must come with it), unzip, and run it.
run_downloaded_script <- function(output, session) {
  session$setInputs(tg_dl_group_tables = character(0), tg_dl_group_plots = character(0),
                    tg_dl_group_diagnostics = character(0), tg_dl_group_data = character(0),
                    tg_dl_group_report = character(0), tg_dl_group_code = "script")
  z <- output$tg_dl_zip
  dir <- tempfile("script_")
  utils::unzip(z, exdir = dir)
  files <- utils::unzip(z, list = TRUE)$Name
  env <- new.env(parent = globalenv())
  old <- setwd(dir)
  on.exit(setwd(old))
  suppressMessages(suppressWarnings(utils::capture.output(
    sys.source("analysis.R", envir = env, keep.source = FALSE))))
  list(files = files, env = env, script = readLines(file.path(dir, "analysis.R")))
}

same_table <- function(a, b, label) {
  expect_false(is.null(a), info = label)
  expect_equal(as.data.frame(a), as.data.frame(b), tolerance = 1e-8,
               ignore_attr = TRUE, info = label)
}

test_that("R20-Z: an LME4 run's script reproduces its tables", {
  rv <- tg_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference",
                      tg_p_adjust = "holm")
    session$setInputs(run_tg_analysis = 1)
    res <- rv$tg_results
    expect_false(isTRUE(res$error))
    expect_match(paste(as.character(output$tg_dl_choices), collapse = " "), "R script", fixed = TRUE)

    run <- run_downloaded_script(output, session)
    # The script comes with the data it reads, although only the script was ticked.
    expect_setequal(run$files, c("README.txt", "analysis.R", "data/analysed_data.csv",
                                 "data/analysed_data.rds"))
    # Every argument the dashboard passed is written, NULLs included: a NULL
    # can differ from the default (cage_column defaults to "Cage").
    expect_true(any(grepl('p_adjust_method = "holm"', run$script, fixed = TRUE)))
    expect_true(any(grepl("dose_column = NULL", run$script, fixed = TRUE)))
    for (f in c("anova", "pairwise_comparisons", "treatment_effects"))
      same_table(run$env$result[[f]], res[[f]], f)
  }))
})

test_that("R20-Z: an AUC run on calculated volumes, with necrosis handling, round-trips", {
  rv <- tg_rv(volume_method = "calculate", col_length = "Length", col_width = "Width",
              col_height = "Height", volume_formula = "ellipsoid", col_necrotic = "Necrotic")
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "auc", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference",
                      tg_p_adjust = "holm", tg_necrotic_handling = "exclude")
    session$setInputs(run_tg_analysis = 1)
    res <- rv$tg_results
    expect_false(isTRUE(res$error))
    run <- run_downloaded_script(output, session)
    # The calculated volume is in the data the script reads.
    expect_true(any(grepl('volume_column = "Volume"', run$script, fixed = TRUE)))
    expect_true(any(grepl('model_type = "auc"', run$script, fixed = TRUE)))
    for (f in c("pairwise_comparisons", "treatment_effects"))
      same_table(run$env$result[[f]], res[[f]], f)
  }))
})

test_that("R20-Z: the randomisation test is written as the call that ran", {
  rv <- tg_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference",
                      tg_p_adjust = "holm", tg_run_perm = TRUE, tg_perm_n = 19L)
    session$setInputs(run_tg_analysis = 1)
    res <- rv$tg_results
    # Outside a server the background worker runs in the session, so the
    # test's result is already in.
    skip_if(is.null(res$permutation_test), "the randomisation test did not finish here")
    run <- run_downloaded_script(output, session)
    # perm_spec() is written as the call the dashboard evaluated, with its unit,
    # number of permutations and seed.
    expect_true(any(grepl('perm_spec(unit = "mouse", n_perm = 19L, seed = 20260730L)',
                          run$script, fixed = TRUE)))
    expect_equal(run$env$perm$p_value, res$permutation_test$p_value)
    expect_equal(run$env$perm$observed_statistic, res$permutation_test$observed_statistic,
                 tolerance = 1e-8)

    # Before the test's result arrives, the script says so instead of running it.
    res$permutation_test <- NULL
    s <- g("tg_reproduce_script")(res)
    expect_true(any(grepl("had not finished", s, fixed = TRUE)))
    expect_false(any(grepl("perm <-", s, fixed = TRUE)))
  }))
})

test_that("R20-Z: a Bayesian run's script names the fit the worker made", {
  # No fit here (the Bayesian tests cover the model); the script is written
  # from the worker's own argument list.
  data <- utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                      package = "mouseExperiment"), stringsAsFactors = FALSE)
  params <- list(data = data,
                 col_names = list(time_column = "Day", volume_column = "Volume",
                                  id_column = "ID", treatment_column = "Treatment",
                                  cage_column = NULL),
                 reference_group = "Vehicle", random_effects = "slope",
                 prior_strength = "skeptical", manual_priors = NULL,
                 n_chains = 2L, n_warmup = 300L, n_iter = 200L, backend = "rstan",
                 necrotic_handling = "none", volume_units = "mm3")
  task <- g(".tg_bayes_task")(params)
  res <- list(model_type_used = "bayes",
              call_record = g("record_call")("bayesian_tumor_growth", task$args),
              posterior_summary = data.frame(x = 1))
  s <- g("tg_reproduce_script")(res)
  expect_silent(parse(text = s))
  expect_true(any(grepl("result <- bayesian_tumor_growth(", s, fixed = TRUE)))
  expect_true(any(grepl("n_chains = 2L", s, fixed = TRUE)))
  expect_true(any(grepl("seed is the package default", s, fixed = TRUE)))
  # The script's call has exactly the worker's arguments.
  call <- parse(text = s)[[grep("bayesian_tumor_growth", vapply(parse(text = s), function(e)
    paste(deparse(e), collapse = ""), ""))[1]]]
  expect_setequal(setdiff(names(as.list(call[[3]]))[-1], "df"), names(task$args)[-1])
})

test_that("R20-Z: the endpoint model keeps its structure when the data move by 1e-13", {
  # Found by the round trip above. The calculated volumes, written with 17
  # significant digits and read back, differ from the originals in 98 of 604
  # values by about 1e-13 (R's parser is not correctly rounded on every
  # platform). On those data the correlated model failed lme4's gradient
  # check where the original passed, so the AUC fell back to uncorrelated
  # random effects and its intervals moved by up to 5%. The endpoint model now
  # restarts from the estimates first, as the main model does (R20-T), and the
  # results agree to about 1e-6.
  d0 <- utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                    package = "mouseExperiment"), stringsAsFactors = FALSE)
  d <- suppressMessages(mouseExperiment::calculate_volume(
    d0, length_column = "Length", width_column = "Width", height_column = "Height",
    formula = "ellipsoid"))
  d2 <- d
  d2$Volume <- as.numeric(sprintf("%.17g", d$Volume))
  skip_if(all(d2$Volume == d$Volume), "this platform reads 17 digits back exactly")
  expect_lt(max(abs(d$Volume - d2$Volume)), 1e-9)
  fit <- function(x) suppressWarnings(suppressMessages(mouseExperiment::tumor_growth_statistics(
    x, model_type = "auc", reference_group = "Vehicle", p_adjust_method = "holm",
    necrotic_column = "Necrotic", necrotic_handling = "exclude", plots = FALSE)))
  num <- c("AUC", "Lower_CL", "Upper_CL")
  a <- as.matrix(fit(d)$treatment_effects[, num])
  b <- as.matrix(fit(d2)$treatment_effects[, num])
  expect_lt(max(abs(a / b - 1)), 1e-4)
})
