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
