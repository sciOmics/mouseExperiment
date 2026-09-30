# Tests for the resume-banner module (CODE_REVIEW.md F.7).

skip_if_not_installed("shiny")
skip_if_not_installed("promises")

pkg <- tryCatch(getNamespace("mouseExperimentDashboard"),
                error = function(e) NULL)
skip_if(is.null(pkg), "mouseExperimentDashboard namespace not loadable")

g <- function(nm) get(nm, envir = pkg)

resolved_promise <- function(value) {
  promises::promise(function(resolve, reject) resolve(value))
}
flush_promises <- function() later::run_now(timeoutSecs = 1)

# ---------------------------------------------------------------------------
# Banner renders nothing when sid has no jobs
# ---------------------------------------------------------------------------
test_that("banner is empty for a sid with no jobs", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)

  shiny::testServer(g("resume_banner_server"),
                    args = list(
                      session_id = shiny::reactive(strrep("a", 64L)),
                      on_resume  = function(job) NULL),
                    {
    # No jobs registered → banner output is NULL.
    expect_null(output$banner)
  })
})

test_that("banner is empty when sid is SID_NONE", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)

  shiny::testServer(g("resume_banner_server"),
                    args = list(
                      session_id = shiny::reactive(g("SID_NONE")),
                      on_resume  = function(job) NULL),
                    {
    expect_null(output$banner)
  })
})

# ---------------------------------------------------------------------------
# Banner renders content for a sid with jobs
# ---------------------------------------------------------------------------
test_that("banner renders job entries for the active sid", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)

  sid <- strrep("b", 64L)
  g("launch_job")(sid, "bayes_tg", resolved_promise(list(value = 1)))
  flush_promises()

  shiny::testServer(g("resume_banner_server"),
                    args = list(
                      session_id = shiny::reactive(sid),
                      on_resume  = function(job) NULL),
                    {
    html <- paste(as.character(output$banner), collapse = "")
    expect_match(html, "Bayesian Tumor Growth", fixed = TRUE)
    expect_match(html, "completed", fixed = TRUE)
    expect_match(html, "Resume",    fixed = TRUE)
  })
})

# ---------------------------------------------------------------------------
# Sid leak audit — the banner HTML never contains the sid string
# ---------------------------------------------------------------------------
test_that("banner HTML never contains the sid value", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)

  sid <- paste0(strrep("c", 32L), strrep("d", 32L))  # 64 hex chars
  g("launch_job")(sid, "bayes_tg", resolved_promise(list(v = 1)))
  flush_promises()

  shiny::testServer(g("resume_banner_server"),
                    args = list(
                      session_id = shiny::reactive(sid),
                      on_resume  = function(job) NULL),
                    {
    html <- paste(as.character(output$banner), collapse = "")
    expect_false(grepl(sid, html, fixed = TRUE))
  })
})

# ---------------------------------------------------------------------------
# Resume dispatch
# ---------------------------------------------------------------------------
test_that("register_resume_handler + dispatch_resume route by kind", {
  g(".reset_resume_handlers")()
  on.exit(g(".reset_resume_handlers")(), add = TRUE)

  received <- NULL
  g("register_resume_handler")("bayes_tg",
    function(job, session, rv) received <<- job)

  fake_job <- list(kind = "bayes_tg", id = "job_test", result = "ok")
  g("dispatch_resume")(fake_job, session = NULL, rv = NULL)
  expect_equal(received$id, "job_test")
})

test_that("dispatch_resume on unknown kind no-ops gracefully", {
  g(".reset_resume_handlers")()
  on.exit(g(".reset_resume_handlers")(), add = TRUE)
  # No handler registered for "bayes_xyz" — should not error.
  fake_job <- list(kind = "bayes_xyz", id = "x")
  # showNotification needs a Shiny session normally; in unit-test mode
  # it falls back silently. Just assert no error.
  expect_silent(
    suppressWarnings(g("dispatch_resume")(fake_job, NULL, NULL))
  )
})

# ---------------------------------------------------------------------------
# TTL sweeper idempotency
# ---------------------------------------------------------------------------
test_that("start_ttl_sweeper is idempotent", {
  g(".reset_ttl_sweeper")()
  on.exit(g(".reset_ttl_sweeper")(), add = TRUE)
  # First call sets the started flag.
  g("start_ttl_sweeper")(interval_seconds = 9999L)  # huge so it doesn't fire
  expect_true(g(".sweeper_state")$started)
  # Second call is a no-op.
  g("start_ttl_sweeper")(interval_seconds = 9999L)
  expect_true(g(".sweeper_state")$started)
})

test_that(".reset_ttl_sweeper removes the pending callback", {
  # A callback left on the global later loop made testServer() hang in every
  # later test file that read an output: MockShinySession waits for an empty
  # loop (shiny:::wait_for_it) on each output read.
  later::run_now(0)
  g(".reset_ttl_sweeper")()
  g("start_ttl_sweeper")(interval_seconds = 9999L)
  expect_false(later::loop_empty())
  g(".reset_ttl_sweeper")()
  expect_true(later::loop_empty())
})
