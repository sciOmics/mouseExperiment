# ============================================================================
# Tests for show_analysis_modal(), the auto-switch-to-Summary behaviour added
# in v26.05.022, and the cache-busting ?v=<version> query string added in
# v26.05.021. CODE_REVIEW.md Round 2 I.3.
# ============================================================================

skip_if_not_installed("shiny", minimum_version = "1.5.0")
skip_if_not_installed("mouseExperiment")

pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(pkg), "mouseExperimentDashboard namespace not loadable")
get_fn <- function(nm) get(nm, envir = pkg)

# ---------------------------------------------------------------------------
# show_analysis_modal()
# ---------------------------------------------------------------------------

test_that("show_analysis_modal fails outside an active session (expected)", {
  # showModal() requires a current session; outside one it raises. The
  # exact error message has varied across Shiny versions; assert only
  # that *some* error is raised so the function isn't silently no-op'd.
  show_fn <- get_fn("show_analysis_modal")
  expect_error(show_fn(message = "Hi", detail = "test"))
})

test_that("show_analysis_modal default arguments are present and stringy", {
  show_fn <- get_fn("show_analysis_modal")
  fmls <- formals(show_fn)
  expect_true("message" %in% names(fmls))
  expect_true("detail"  %in% names(fmls))
  # Defaults should be non-empty strings so calls with no args still
  # display *something*.
  default_message <- eval(fmls$message)
  default_detail  <- eval(fmls$detail)
  expect_true(is.character(default_message) && nzchar(default_message))
  expect_true(is.character(default_detail)  && nzchar(default_detail))
})

# ---------------------------------------------------------------------------
# auto-switch to "Summary" tab on analysis completion (v26.05.022)
# ---------------------------------------------------------------------------

test_that("updateTabsetPanel(... selected = 'Summary') is wired in every module", {
  # The auto-switch behaviour is implemented as
  #   updateTabsetPanel(session, "<X>_results_tabs", selected = "Summary")
  # called once per completion path. A regression that drops one of these
  # would leave the user staring at the wrong tab. Easiest assertion: each
  # module file contains a corresponding updateTabsetPanel call.
  module_files <- c(
    "R/mod_tumor_growth.R",
    "R/mod_survival.R",
    "R/mod_synergy.R",
    "R/mod_dose_response.R",
    "R/mod_toxicity.R",
    "R/mod_power.R"
  )
  for (f in module_files) {
    txt <- strsplit(dash_source(f), "\n", fixed = TRUE)[[1]]   # skips outside the source tree
    # R20.D6: a run that finishes on the background worker switches the tab
    # from its finalizer, which lives beside the module.
    helper <- c("R/mod_toxicity.R" = "R/tox_background.R")[f]
    if (!is.na(helper)) {
      txt <- c(txt, strsplit(dash_source(helper), "\n", fixed = TRUE)[[1]])
    }
    has_switch <- any(grepl(
      "updateTabsetPanel\\(session,.*selected\\s*=\\s*[\"']Summary[\"']", txt
    ))
    expect_true(has_switch,
                info = paste0(f, " is missing updateTabsetPanel(..., selected = 'Summary')"))
  }
})

# ---------------------------------------------------------------------------
# Cache-busting ?v=<version> on static asset includes (v26.05.021)
# ---------------------------------------------------------------------------

test_that("static asset includes carry a ?v=<version> query string", {
  # The cache-bust mechanism appends a version-tagged query to every
  # <link>/<script> in the head so browsers fetch fresh assets after
  # every release. We grep the source for the pattern; the actual head
  # construction is too verbose to mock here.
  app_path <- file.path(dash_source_dir(), "R", "app.R")
  expect_true(file.exists(app_path), info = "R/app.R must exist")
  txt <- readLines(app_path, warn = FALSE)

  # Look for the "v=" cache-bust string OR a paste-style call that
  # appends APP_VERSION onto an asset URL.
  has_cache_bust <- any(grepl("\\?v=", txt)) ||
    any(grepl("paste0?\\(.*\\.css.*[Vv]ersion", txt)) ||
    any(grepl("paste0?\\(.*\\.js.*[Vv]ersion", txt))
  expect_true(has_cache_bust,
              info = "R/app.R should append ?v=<version> to static asset URLs (Round 2 A.4)")
})
