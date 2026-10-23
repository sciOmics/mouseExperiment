# Regression tests for R20.D6 (CODE_REVIEW.md R20-D.R): the power simulation,
# the randomisation test and the body-weight model go through the background
# worker's queue instead of running in the Shiny process, where they froze every
# session on the server. With the worker held busy (mocked), each run returns at
# once with its job queued; freed, the result arrives. Each test fails on
# v26.10.003, where the work ran inline.

skip_if_not_installed("shiny", minimum_version = "1.5.0")
skip_if_not_installed("future", minimum_version = "1.49.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

# The scheduler is driven by hand: no `later` timer is left running.
local_registry <- function(env = parent.frame()) {
  old_opt  <- options(mE.job_pump = FALSE)
  old_plan <- future::plan()
  g(".reset_job_registry")()
  withr::defer({
    g(".reset_job_registry")()
    future::plan(old_plan)
    options(old_opt)
  }, envir = env)
}

sid <- strrep("d", 64L)
jobs_of <- function() unname(g("list_jobs")(sid))

all_features_rv <- function() {
  path <- system.file("sample_data", c("all_features_synthetic_data.csv",
                                       "master_synthetic_data.csv"),
                      package = "mouseExperiment")[1]
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(path, stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = sid,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "Cage",
                    col_weight = "Weight", control_group = "Vehicle",
                    volume_units = "mm3", randomisation_unit = "mouse"))
}

test_that("R20.D6: a worker task brings its warnings back with the value", {
  f <- function() { warning("check the units"); list(a = 1) }
  environment(f) <- baseenv()
  t <- g("background_task")(f, list())
  v <- do.call(t$fun, t$args)
  expect_equal(v$a, 1)
  expect_identical(v$analysis_warnings, "check the units")
  # A NULL result keeps its reason.
  h <- function() { warning("skipped: too few cages"); NULL }
  environment(h) <- baseenv()
  t <- g("background_task")(h, list())
  expect_identical(do.call(t$fun, t$args)$analysis_warnings, "skipped: too few cages")
})

test_that("R20.D6: the power simulation is queued on the worker, then delivered", {
  local_registry()
  future::plan(future::sequential)
  free <- 0L
  testthat::local_mocked_bindings(.free_workers = function() free,
                                  .package = "mouseExperiment")
  rv <- shiny::reactiveValues(uploaded_data = NULL, data_valid = FALSE,
                              settings = list(), session_id = sid)
  suppressWarnings(shiny::testServer(g("power_server"), args = list(rv = rv), {
    session$setInputs(pw_mode = "apriori_lmm", pw_alpha = "0.05",
                      pw_lmm_n_per_group = "4", pw_lmm_n_sims = 10,
                      pw_lmm_n_groups = 2, pw_lmm_n_timepoints = 6)
    session$setInputs(run_pw_analysis = 1)
    j <- jobs_of()
    expect_length(j, 1L)
    expect_identical(j[[1]]$kind, "power_sim")
    expect_identical(j[[1]]$status, "queued")       # v26.10.003 ran it inline
    expect_null(rv$ap_results)
    free <<- 1L
    for (i in 1:3) session$elapse(2100)
    expect_identical(rv$ap_results$mode, "apriori_lmm")
    expect_true(all(c("Power", "Power_Lower") %in% names(rv$ap_results$scenario_table)))
    expect_length(jobs_of(), 0L)                    # retired once delivered
  }))
})

test_that("R20.D6: the randomisation test runs on the worker after the fit is shown", {
  local_registry()
  future::plan(future::sequential)
  free <- 0L
  testthat::local_mocked_bindings(.free_workers = function() free,
                                  .package = "mouseExperiment")
  rv <- all_features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log",
                      tg_run_perm = TRUE, tg_perm_n = 99)
    session$setInputs(run_tg_analysis = 1)
    # The fit is on screen at once; the test waits in the queue.
    expect_false(isTRUE(rv$tg_results$error))
    expect_null(rv$tg_results$permutation_test)
    expect_identical(rv$tg_perm_status, "running")
    expect_identical(jobs_of()[[1]]$kind, "tg_perm")
    expect_match(paste(output$tg_perm_test_ui$html, collapse = " "), "background worker")
    free <<- 1L
    for (i in 1:3) session$elapse(2100)
    pt <- rv$tg_results$permutation_test
    expect_true(is.numeric(pt$p_value))
    expect_null(rv$tg_perm_status)
    expect_match(paste(output$tg_perm_test_ui$html, collapse = " "), "Permutation p =")
  }))
})

test_that("R20.D6: a randomisation test that arrives after a newer run is discarded", {
  rv <- shiny::reactiveValues(tg_results = list(run_token = "new"), tg_perm_status = "running")
  shiny::isolate(g(".finalize_perm_success")(list(p_value = 0.01),
                                             list(run_token = "old"), NULL, NULL, NULL, rv))
  expect_null(shiny::isolate(rv$tg_results$permutation_test))
})

test_that("R20.D6: the body-weight model is queued on the worker, then delivered", {
  local_registry()
  future::plan(future::sequential)
  free <- 0L
  testthat::local_mocked_bindings(.free_workers = function() free,
                                  .package = "mouseExperiment")
  rv <- all_features_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE,
                      tox_tumor_density = 1, tox_covariates = "volume",
                      tox_estimation = "REML", tox_threshold = 20,
                      tox_baseline_day = "first")
    session$setInputs(run_tox_analysis = 1)
    expect_identical(jobs_of()[[1]]$kind, "bw_model")
    expect_identical(jobs_of()[[1]]$status, "queued")
    expect_null(rv$tox_results)
    free <<- 1L
    for (i in 1:3) session$elapse(2100)
    mm <- rv$tox_results$mixed_model
    expect_false(isTRUE(mm$error))
    expect_true("Weight_Change_Per_Day" %in% names(mm$emmeans_table))
    expect_false(isTRUE(rv$tox_results$weight_loss$error))  # computed after the model
    expect_gt(nrow(rv$tox_results$trajectory_data), 0L)
  }))
})
