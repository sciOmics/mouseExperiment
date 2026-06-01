# =============================================================================
# Shinytest2: integration tests for the dashboard server modules
#
# These tests start a real (headless) Shiny session and interact with it
# programmatically.  They require:
#   - shinytest2 >= 0.2.0
#   - chromote / Chrome accessible on PATH  (skipped on CI if absent)
#   - mouseExperiment installed
#   - mouseExperimentDashboard package namespace loadable
#
# If any requirement is missing the entire file skips gracefully.
# =============================================================================

skip_if_not_installed("shinytest2")
skip_if_not_installed("shiny")
skip_if_not_installed("chromote",
                       message = "chromote not installed — skipping Shiny app tests")

# Check Chrome / Chromium is available
skip_if(
  tryCatch({ chromote::find_chrome(); FALSE }, error = function(e) TRUE),
  message = "Chrome/Chromium not found — skipping Shiny app tests"
)

# Build the minimal test app from helper-app.R
app_obj <- build_test_app()
skip_if(is.null(app_obj), "Could not build test app — package namespace unavailable")

# ===========================================================================
# Helpers
# ===========================================================================

start_driver <- function(timeout = 15000) {
  shinytest2::AppDriver$new(
    app_obj,
    name              = "dashboard-data-tab",
    timeout           = timeout,
    wait              = TRUE
  )
}

# ===========================================================================
# Test: app starts and renders the Data tab
# ===========================================================================
test_that("app starts and renders the Data tab without error", {
  app <- start_driver()
  on.exit(app$stop(), add = TRUE)

  expect_no_error(app$wait_for_idle(timeout = 10000))
  html <- app$get_html("#shiny-disconnected-overlay", timeout = 3000) %||% ""
  expect_false(grepl("disconnected", html, ignore.case = TRUE),
               info = "App disconnected immediately after start")
})

# ===========================================================================
# Test: demo data combo loads successfully
# ===========================================================================
test_that("Load Demo Data (combo) sets data_valid to TRUE", {
  app <- start_driver()
  on.exit(app$stop(), add = TRUE)

  # Select the combo dataset
  app$set_inputs(demo_dataset = "combo", wait_ = FALSE)

  # Click the Load Demo Data button and wait for reactivity to settle
  app$click("load_demo_data")
  app$wait_for_idle(timeout = 15000)

  valid_text <- app$get_text("#test_data_valid")
  expect_equal(trimws(valid_text), "TRUE",
               info = paste("Expected data_valid=TRUE after loading combo demo, got:", valid_text))
})

test_that("Load Demo Data (combo) loads non-zero rows", {
  app <- start_driver()
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(demo_dataset = "combo", wait_ = FALSE)
  app$click("load_demo_data")
  app$wait_for_idle(timeout = 15000)

  nrow_text <- app$get_text("#test_nrow")
  nrow_val  <- suppressWarnings(as.integer(trimws(nrow_text)))
  expect_true(!is.na(nrow_val) && nrow_val > 0,
              info = paste("Expected nrow > 0, got:", nrow_text))
})

# ===========================================================================
# Test: demo data dose_levels loads successfully
# ===========================================================================
test_that("Load Demo Data (dose_levels) sets data_valid to TRUE", {
  app <- start_driver()
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(demo_dataset = "dose_levels", wait_ = FALSE)
  app$click("load_demo_data")
  app$wait_for_idle(timeout = 15000)

  valid_text <- app$get_text("#test_data_valid")
  expect_equal(trimws(valid_text), "TRUE",
               info = paste("Expected data_valid=TRUE after loading dose_levels demo, got:", valid_text))
})

# ===========================================================================
# Test: success message visible after loading combo
# ===========================================================================
test_that("success alert appears after loading combo demo data", {
  app <- start_driver()
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(demo_dataset = "combo", wait_ = FALSE)
  app$click("load_demo_data")
  app$wait_for_idle(timeout = 15000)

  html <- app$get_html("body")
  expect_true(
    grepl("alert-success", html, fixed = TRUE) ||
    grepl("Loaded|loaded|success", html, ignore.case = TRUE),
    info = "Expected success message in page HTML after loading demo data"
  )
})

# ===========================================================================
# Test: column mapping UI appears after loading data
# ===========================================================================
test_that("column mapping UI renders after loading demo data", {
  app <- start_driver()
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(demo_dataset = "combo", wait_ = FALSE)
  app$click("load_demo_data")
  app$wait_for_idle(timeout = 15000)

  html <- app$get_html("body")
  # Column mapping selectors should be present
  has_col_mapping <- grepl("col_id|col_treatment|col_day|Column Mapping",
                           html, ignore.case = TRUE)
  expect_true(has_col_mapping,
              info = "Expected column mapping UI to appear after loading demo data")
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ===========================================================================
# Test: end-to-end upload → run TG LME4 → verify results populate
# (CODE_REVIEW.md Round 2 I.5)
# ===========================================================================
#
# Extends shinytest2 coverage beyond the Data tab. Drives a complete
# user workflow: load combo demo, switch to the TG tab, accept defaults,
# run the LME4 fit, wait for completion, assert tg_results is non-NULL.
#
# LME4 on the combo fixture completes in < 10s. Bayesian / GAMM paths
# are NOT exercised here because they require brms + Stan and would
# blow CI budgets. They're covered by mouseExperiment's unit tests.

test_that("end-to-end: load combo demo + run TG LME4 populates results", {
  skip_if_not_installed("mouseExperiment")

  workflow_app <- build_workflow_test_app()
  skip_if(is.null(workflow_app),
          "Could not build workflow test app — package namespace incomplete")

  app <- shinytest2::AppDriver$new(
    workflow_app,
    name    = "dashboard-workflow",
    timeout = 30000,
    wait    = TRUE
  )
  on.exit(app$stop(), add = TRUE)

  # Step 1 — load combo demo data
  app$set_inputs(demo_dataset = "combo", wait_ = FALSE)
  app$click("load_demo_data")
  app$wait_for_idle(timeout = 15000)
  expect_equal(trimws(app$get_text("#test_data_valid")), "TRUE",
               info = "Data should load successfully")

  # Step 2 — switch to TG tab
  app$set_inputs(main_tabs = "TG")
  app$wait_for_idle(timeout = 5000)

  # Step 3 — click Run; LME4 is the default model_type
  app$click("tg-run_tg_analysis")

  # Step 4 — poll up to 30s for results to populate
  deadline <- Sys.time() + 30
  populated <- FALSE
  while (Sys.time() < deadline) {
    app$wait_for_idle(timeout = 2000)
    tg_state <- trimws(app$get_text("#test_tg_has_results") %||% "FALSE")
    if (identical(tg_state, "TRUE")) {
      populated <- TRUE
      break
    }
    Sys.sleep(1)
  }
  expect_true(populated,
              info = "TG LME4 analysis should populate rv$tg_results within 30s")
})
