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
# ---------------------------------------------------------------------------
make_test_rv <- function() {
  data(combo_treatment_synthetic_data, package = "mouseExperiment")
  df <- combo_treatment_synthetic_data
  # Ensure Volume column exists (combo data may already have it)
  if (!"Volume" %in% colnames(df) && all(c("Length", "Width") %in% colnames(df))) {
    df <- mouseExperiment::calculate_volume(df)
  }

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
    tg_results     = NULL,
    tg_running     = FALSE,
    tg_col_names   = NULL,
    tg_plot_groups = NULL,
    surv_results   = NULL,
    surv_running   = FALSE,
    surv_plot_groups = NULL,
    syn_results    = NULL,
    syn_running    = FALSE,
    syn_over_time  = FALSE,
    pw_results     = NULL,
    pw_running     = FALSE,
    dr_results     = NULL,
    dr_running     = FALSE
  )
  rv
}

# ===========================================================================
# Test: dose_response_server reactives
# ===========================================================================
test_that("dose_response_server initialises without error", {
  dose_response_server_fn <- get_fn("dose_response_server")

  shiny::testServer(dose_response_server_fn, args = list(rv = make_test_rv()), {
    # Module should initialise without errors
    expect_true(TRUE)
  })
})

# ===========================================================================
# Test: power_server initialises without error
# ===========================================================================
test_that("power_server initialises without error", {
  power_server_fn <- get_fn("power_server")

  shiny::testServer(power_server_fn, args = list(rv = make_test_rv()), {
    expect_true(TRUE)
  })
})

# ===========================================================================
# Test: survival_server initialises without error
# ===========================================================================
test_that("survival_server initialises without error", {
  survival_server_fn <- get_fn("survival_server")

  shiny::testServer(survival_server_fn, args = list(rv = make_test_rv()), {
    expect_true(TRUE)
  })
})

# ===========================================================================
# Test: synergy_server initialises without error
# ===========================================================================
test_that("synergy_server initialises without error", {
  synergy_server_fn <- get_fn("synergy_server")

  shiny::testServer(synergy_server_fn, args = list(rv = make_test_rv()), {
    expect_true(TRUE)
  })
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
