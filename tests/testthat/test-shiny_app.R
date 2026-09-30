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
# R20.D11: skip_if_not_installed() has no `message` argument, so this line
# errored at file level whenever shinytest2 was installed and none of the tests
# below ever ran.
skip_if_not_installed("chromote")

# Check Chrome / Chromium is available
skip_if(
  tryCatch({ chromote::find_chrome(); FALSE }, error = function(e) TRUE),
  message = "Chrome/Chromium not found — skipping Shiny app tests"
)

# The apps live in tests/testthat/apps/<name>/app.R and load the package
# themselves: AppDriver runs them in a separate R process, which cannot see a
# package loaded from source in this one (R20.D11).
upload_app <- test_path("apps", "data-upload")

# ===========================================================================
# Helpers
# ===========================================================================

start_driver <- function(timeout = 15000) {
  shinytest2::AppDriver$new(
    upload_app,
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
  html <- app$get_html("#shiny-disconnected-overlay") %||% ""
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
  app$set_inputs(`data-demo_dataset` = "combo", wait_ = FALSE)

  # Click the Load Demo Data button and wait for reactivity to settle
  app$click("data-load_demo_data")
  app$wait_for_idle(timeout = 15000)

  valid_text <- app$get_text("#test_data_valid")
  expect_equal(trimws(valid_text), "TRUE",
               info = paste("Expected data_valid=TRUE after loading combo demo, got:", valid_text))
})

test_that("Load Demo Data (combo) loads non-zero rows", {
  app <- start_driver()
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(`data-demo_dataset` = "combo", wait_ = FALSE)
  app$click("data-load_demo_data")
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

  app$set_inputs(`data-demo_dataset` = "dose_levels", wait_ = FALSE)
  app$click("data-load_demo_data")
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

  app$set_inputs(`data-demo_dataset` = "combo", wait_ = FALSE)
  app$click("data-load_demo_data")
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

  app$set_inputs(`data-demo_dataset` = "combo", wait_ = FALSE)
  app$click("data-load_demo_data")
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
# Test: the full dashboard starts a working session (R20.D3)
# ===========================================================================
# The first version of the R20.D3 fix read a reactive value in the top-level
# server body. That closed every session at start-up, while the page still
# loaded with HTTP 200 and testServer() -- which runs the body inside a reactive
# context -- passed.
test_that("the full dashboard starts a working session", {
  app <- shinytest2::AppDriver$new(
    test_path("apps", "dashboard"),
    name         = "dashboard-full",
    timeout      = 30000,
    load_timeout = 90000,
    wait         = TRUE
  )
  on.exit(app$stop(), add = TRUE)

  expect_true(isTRUE(app$get_js("Shiny.shinyapp.isConnected()")))
  # The server sends the disclaimer modal as the session starts.
  expect_true(isTRUE(app$get_js("$('.modal-dialog').length > 0")))

  # And the modules respond: loading a demo maps the weight column.
  app$run_js("$('.modal-footer button').click();")
  app$set_inputs(`data-demo_dataset` = "master", wait_ = FALSE)
  app$click("data-load_demo_data")
  app$wait_for_js("$('#data-col_weight').val() === 'Weight'", timeout = 30000)
  expect_true(isTRUE(app$get_js("Shiny.shinyapp.isConnected()")))
})

# ===========================================================================
# Test: end-to-end upload → run TG LME4 → verify results populate
# (CODE_REVIEW.md Round 2 I.5)
# ===========================================================================
#
# Extends shinytest2 coverage beyond the Data tab. Drives a complete
# user workflow: load combo demo, switch to the TG tab, accept defaults,
# run the LME4 fit, wait for completion, assert tg_results is non-NULL.
#
# LME4 on the combo fixture completes in < 10s. Bayesian paths
# are NOT exercised here because they require brms + Stan and would
# blow CI budgets. They're covered by mouseExperiment's unit tests.

test_that("end-to-end: load combo demo + run TG LME4 populates results", {
  skip_if_not_installed("mouseExperiment")

  app <- shinytest2::AppDriver$new(
    test_path("apps", "workflow"),
    name    = "dashboard-workflow",
    timeout = 30000,
    wait    = TRUE
  )
  on.exit(app$stop(), add = TRUE)

  # Step 1 — load combo demo data
  app$set_inputs(`data-demo_dataset` = "combo", wait_ = FALSE)
  app$click("data-load_demo_data")
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
