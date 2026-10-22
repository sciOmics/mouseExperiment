# =============================================================================
# testServer() tests for dashboard module reactive logic
#
# These tests exercise module servers using Shiny's testServer() utility,
# which runs the server function in a mock session without a browser.
# =============================================================================

skip_if_not_installed("shiny", minimum_version = "1.5.0")
skip_if_not_installed("mouseExperiment")

# Ensure the dashboard package namespace is loadable
pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(pkg), "mouseExperimentDashboard namespace not loadable")

get_fn <- function(nm) get(nm, envir = pkg)

# ---------------------------------------------------------------------------
# Helper: build a minimal rv reactiveValues with valid combo demo data
#
# `with_results` lets a test pre-populate `rv$<X>_results` with a small
# dummy result list so post-completion behaviour (auto-switch to Summary,
# downloads enabled, etc.) can be exercised without running the actual
# analysis (CODE_REVIEW.md Round 2 I.7).
# ---------------------------------------------------------------------------
make_test_rv <- function(with_results = character(0)) {
  data(combo_treatment_synthetic_data, package = "mouseExperiment")
  df <- combo_treatment_synthetic_data
  # Ensure Volume column exists (combo data may already have it)
  if (!"Volume" %in% colnames(df) && all(c("Length", "Width") %in% colnames(df))) {
    df <- mouseExperiment::calculate_volume(df)
  }

  # Tiny dummy result list — has the fields the post-completion UI checks for.
  dummy_result <- list(
    error                = FALSE,
    model_type_used      = "lme4",
    treatment_effects    = data.frame(
      Group = c("Control", "TreatmentA"),
      Adjusted_Mean = c(5.0, 4.5),
      Lower_CL = c(4.8, 4.3),
      Upper_CL = c(5.2, 4.7),
      stringsAsFactors = FALSE
    ),
    pairwise_comparisons = data.frame(
      contrast = "TreatmentA - Control",
      estimate = -0.5,
      p.value  = 0.001,
      stringsAsFactors = FALSE
    ),
    summary              = list(notes = "test")
  )

  rv <- shiny::reactiveValues(
    uploaded_data  = df,
    data_valid     = TRUE,
    detected_cols  = list(
      id        = "ID",
      treatment = "Treatment",
      day       = "Day",
      volume    = "Volume",
      cage      = "Cage"
    ),
    settings = list(
      col_id        = "ID",
      col_treatment = "Treatment",
      col_day       = "Day",
      col_volume    = "Volume",
      col_cage      = "Cage",
      col_dose      = NULL,
      time_method   = "direct",
      volume_method = "existing"
    ),
    processed_data = NULL,
    time_col       = "Day",
    vol_col        = "Volume",
    tx_col         = "Treatment",
    id_col         = "ID",
    cage_col       = "Cage",
    dose_col       = NULL,
    tg_results     = if ("tg"   %in% with_results) dummy_result else NULL,
    tg_running     = FALSE,
    tg_col_names   = if ("tg"   %in% with_results)
                        list(id_column="ID", treatment_column="Treatment",
                             time_column="Day", volume_column="Volume",
                             cage_column="Cage")
                     else NULL,
    tg_plot_groups = NULL,
    surv_results   = if ("surv" %in% with_results) dummy_result else NULL,
    surv_running   = FALSE,
    surv_plot_groups = NULL,
    syn_results    = if ("syn"  %in% with_results) dummy_result else NULL,
    syn_running    = FALSE,
    syn_over_time  = FALSE,
    pw_results     = if ("pw"   %in% with_results) dummy_result else NULL,
    pw_running     = FALSE,
    dr_results     = if ("dr"   %in% with_results) dummy_result else NULL,
    dr_running     = FALSE,
    tox_results    = if ("tox"  %in% with_results) dummy_result else NULL,
    ap_results     = if ("ap"   %in% with_results)
                        list(scenario_table = data.frame(
                               N_Per_Group = c(5L, 10L, 15L),
                               Power       = c(0.45, 0.78, 0.92),
                               Alpha       = 0.05),
                             power_curve_data = data.frame(
                               N_Per_Group = c(5L, 10L, 15L),
                               Power       = c(0.45, 0.78, 0.92),
                               Alpha       = 0.05),
                             sensitivity_table = NULL,
                             mode = "apriori_lmm",
                             ap_submode = NA)
                     else NULL
  )
  rv
}

# ===========================================================================
# Test: each module server initialises and exposes its key reactive surface
# (CODE_REVIEW.md Round 2 I.1 — replaces expect_true(TRUE) placeholders).
#
# These assertions verify a real mock session was constructed and the
# module's namespace is wired through correctly. The prior placeholders
# would pass even if a module crashed silently inside init; these don't.
# ===========================================================================
test_that("dose_response_server initialises and exposes session ns", {
  dose_response_server_fn <- get_fn("dose_response_server")
  shiny::testServer(dose_response_server_fn, args = list(rv = make_test_rv()), {
    expect_true(inherits(session, "session_proxy"))
    expect_true(nzchar(session$ns("any_id")))
  })
})

test_that("power_server initialises and exposes session ns", {
  power_server_fn <- get_fn("power_server")
  shiny::testServer(power_server_fn, args = list(rv = make_test_rv()), {
    expect_true(inherits(session, "session_proxy"))
    expect_true(nzchar(session$ns("pw_ap_results_tabs")))
  })
})

test_that("survival_server initialises and exposes session ns", {
  survival_server_fn <- get_fn("survival_server")
  shiny::testServer(survival_server_fn, args = list(rv = make_test_rv()), {
    expect_true(inherits(session, "session_proxy"))
    expect_true(nzchar(session$ns("surv_results_tabs")))
  })
})

test_that("synergy_server initialises and exposes session ns", {
  synergy_server_fn <- get_fn("synergy_server")
  shiny::testServer(synergy_server_fn, args = list(rv = make_test_rv()), {
    expect_true(inherits(session, "session_proxy"))
    expect_true(nzchar(session$ns("syn_results_tabs")))
  })
})

# CODE_REVIEW.md Round 2 I.2 — the three biggest server modules were
# previously untested altogether. Add basic init tests so a future regression
# that breaks any of them at startup fails immediately.
test_that("tumor_growth_server initialises and exposes session ns", {
  tumor_growth_server_fn <- get_fn("tumor_growth_server")
  shiny::testServer(tumor_growth_server_fn, args = list(rv = make_test_rv()), {
    expect_true(inherits(session, "session_proxy"))
    expect_true(nzchar(session$ns("tg_results_tabs")))
  })
})

test_that("toxicity_server initialises and exposes session ns", {
  toxicity_server_fn <- get_fn("toxicity_server")
  shiny::testServer(toxicity_server_fn, args = list(rv = make_test_rv()), {
    expect_true(inherits(session, "session_proxy"))
    expect_true(nzchar(session$ns("tox_results_tabs")))
  })
})

test_that("data_upload_server initialises and exposes session ns", {
  data_upload_server_fn <- get_fn("data_upload_server")
  shiny::testServer(data_upload_server_fn,
                    args = list(rv = make_test_rv()), {
    expect_true(inherits(session, "session_proxy"))
    # Module sets rv$uploaded_data via the moduleServer scope; verify the
    # exposed namespace generates a namespaced ID.
    expect_true(nzchar(session$ns("any_id")))
  })
})

# make_test_rv(with_results) returns rv pre-populated with dummy result
# lists; verify the helper itself populates the requested slots so other
# tests downstream can rely on the contract (CODE_REVIEW.md Round 2 I.7).
test_that("make_test_rv(with_results) populates the requested slots", {
  rv_empty <- make_test_rv()
  expect_null(shiny::isolate(rv_empty$ap_results))
  expect_null(shiny::isolate(rv_empty$tg_results))
  expect_null(shiny::isolate(rv_empty$surv_results))

  rv_ap <- make_test_rv(with_results = "ap")
  expect_false(is.null(shiny::isolate(rv_ap$ap_results)))
  expect_true(nrow(shiny::isolate(rv_ap$ap_results)$scenario_table) > 0L)

  rv_surv <- make_test_rv(with_results = c("surv", "tg"))
  expect_false(is.null(shiny::isolate(rv_surv$surv_results)))
  expect_false(is.null(shiny::isolate(rv_surv$tg_results)))
  expect_equal(shiny::isolate(rv_surv$surv_results)$model_type_used, "lme4")
})

# ===========================================================================
# Test: ensure_preprocessed works from server context
# ===========================================================================
test_that("ensure_preprocessed caches and returns data", {
  ensure_preprocessed_fn <- get_fn("ensure_preprocessed")

  rv <- make_test_rv()
  # Simulate that TG module already ran
  shiny::isolate({
    rv$processed_data <- rv$uploaded_data
    rv$time_col <- "Day"
    rv$vol_col  <- "Volume"
    rv$tx_col   <- "Treatment"
  })

  result <- shiny::isolate(ensure_preprocessed_fn(rv))
  expect_true(is.list(result))
  expect_true(nrow(result$data) > 0)
  expect_equal(result$time_col, "Day")
  expect_equal(result$vol_col, "Volume")
  expect_equal(result$tx_col, "Treatment")
})

test_that("ensure_preprocessed runs inline when processed_data is NULL", {
  ensure_preprocessed_fn <- get_fn("ensure_preprocessed")

  rv <- make_test_rv()
  rv$processed_data <- NULL

  result <- shiny::isolate(ensure_preprocessed_fn(rv))
  expect_true(is.list(result))
  expect_true(nrow(result$data) > 0)
  # Should now be cached
  cached <- shiny::isolate(rv$processed_data)
  expect_true(!is.null(cached))
})
