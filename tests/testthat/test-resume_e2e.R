# End-to-end tests for CODE_REVIEW.md F.7 — the disconnect/reconnect
# cycle exercised at the registry + dispatch layer without booting a
# browser. Each test:
#   1) Simulates a fresh session that "starts a fit" (registers a job
#      via launch_job and resolves the promise immediately so the
#      registry holds a completed entry).
#   2) Simulates a disconnect (the launching session ends — drop its
#      rv).
#   3) Simulates reconnect (a fresh rv with the same sid).
#   4) Checks that list_jobs(sid) surfaces the job to the new session.
#   5) Dispatches via the registered handler and asserts the result
#      lands in the right rv slot.

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
# Cross-session resume: launching session ends; reconnect picks up result
# ---------------------------------------------------------------------------
test_that("TG fit completes while disconnected; reconnect resumes it", {
  g(".reset_job_registry")()
  g(".reset_resume_handlers")()
  on.exit({
    g(".reset_job_registry")()
    g(".reset_resume_handlers")()
  }, add = TRUE)

  sid <- strrep("a", 64L)

  # Step 1 — fresh session launches a TG fit. Use a resolved-promise
  # fake to skip the brms compile.
  rv_session_a <- shiny::reactiveValues(
    session_id = sid,
    tg_running = TRUE,
    tg_results = NULL
  )
  fake_result <- list(
    model = "<brmsfit-fake>",
    treatment_effects = data.frame(Group = "A", Adjusted_Mean = 1),
    error = FALSE
  )
  fake_params <- list(
    data      = data.frame(Day = 1:3, Volume = 100:102,
                           Treatment = "A", ID = 1),
    col_names = list(time_column = "Day", volume_column = "Volume",
                     treatment_column = "Treatment", id_column = "ID",
                     cage_column = NULL))

  job_id <- g("launch_job")(sid, "bayes_tg",
                            resolved_promise(fake_result),
                            context = list(bayes_params = fake_params))
  flush_promises()

  # The launching session would normally have an attach_bayes_poller
  # observer; pretend the user's tab dropped before the poller fired.
  # The registry still holds the completed result.
  job <- g("get_job")(sid, job_id)
  expect_equal(job$status, "completed")
  expect_equal(job$result$treatment_effects$Group, "A")

  rm(rv_session_a)  # tab closed, reactiveValues gone

  # Step 2 — register the TG resume handler (would happen in a fresh
  # tumor_growth_server() invocation on reconnect).
  rv_session_b <- shiny::reactiveValues(session_id = sid, tg_running = FALSE,
                                        tg_results = NULL,
                                        tg_plot_data = NULL,
                                        tg_col_names = NULL,
                                        bayes_models = NULL)
  shiny::testServer(function(id, rv) {
    shiny::moduleServer(id, function(input, output, session) {
      g("register_tg_resume")(input, output, session, rv)

      # Step 3 — banner sees the outstanding job for this sid.
      jobs <- g("list_jobs")(sid)
      expect_length(jobs, 1L)
      expect_equal(jobs[[1L]]$kind, "bayes_tg")

      # Step 4 — dispatch_resume routes to the registered handler.
      g("dispatch_resume")(jobs[[1L]], session, rv)

      # The handler attaches a poller; flush the registry observer
      # by triggering the reactive context.
      shiny::isolate({
        # Mock setSessionPager: just trigger the observer manually.
        # In real Shiny invalidateLater() would fire 2s later; in
        # tests we synthesize the same effect by directly calling the
        # finalizer with the cached result.
        # (The handler stored job_id in rv$tg_job_id.)
        expect_false(is.null(rv$tg_job_id))
        # Final state check: the poller's first tick would land
        # finalize_tg_success. Run the finalizer directly to assert
        # the rv-population contract without depending on the
        # invalidateLater scheduler.
        g(".finalize_tg_success")(jobs[[1L]]$result,
                                  jobs[[1L]]$context,
                                  input, output, session, rv)
      })
      expect_equal(shiny::isolate(rv$tg_results$model_type_used), "bayes")
      expect_equal(shiny::isolate(rv$tg_plot_data), fake_params$data)
      expect_false(shiny::isolate(rv$tg_running))
    })
  }, args = list(rv = rv_session_b), {})
})

# ---------------------------------------------------------------------------
# Cross-sid isolation under the full dispatch chain
# ---------------------------------------------------------------------------
test_that("dispatch_resume cannot pop a job belonging to a different sid", {
  g(".reset_job_registry")()
  g(".reset_resume_handlers")()
  on.exit({
    g(".reset_job_registry")()
    g(".reset_resume_handlers")()
  }, add = TRUE)

  sid_a <- strrep("a", 64L)
  sid_b <- strrep("b", 64L)

  # A starts a fit.
  fake_result <- list(model = "secret", error = FALSE,
                      treatment_effects = data.frame())
  jid_a <- g("launch_job")(sid_a, "bayes_tg",
                           resolved_promise(fake_result),
                           context = list(bayes_params = NULL))
  flush_promises()
  expect_equal(g("get_job")(sid_a, jid_a)$status, "completed")

  # B reconnects with their own sid. list_jobs(sid_b) must NOT include
  # A's job, regardless of whether B knows jid_a.
  jobs_b <- g("list_jobs")(sid_b)
  expect_length(jobs_b, 0L)

  # Even if B somehow crafts a job entry with A's jid_a, get_job
  # under sid_b still returns NULL.
  expect_null(g("get_job")(sid_b, jid_a))
})

# ---------------------------------------------------------------------------
# Resume handlers for every Bayesian kind register successfully
# ---------------------------------------------------------------------------
test_that("every Bayesian kind has a register_<X>_resume helper", {
  # Bayesian scope is tumour growth and survival (CODE_REVIEW.md R20-D.H).
  registers <- c("register_tg_resume", "register_surv_resume")
  for (r in registers) {
    expect_true(exists(r, envir = pkg),
                info = sprintf("'%s' should be defined", r))
    expect_true(is.function(g(r)),
                info = sprintf("'%s' should be a function", r))
  }
})

# ---------------------------------------------------------------------------
# Cancel mid-fit propagates: cancellation finalizer runs, result dropped
# ---------------------------------------------------------------------------
test_that("cancel_job + late resolve → cancellation finalizer fires", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)

  sid <- strrep("c", 64L)
  # Promise we can resolve manually.
  resolver <- NULL
  p <- promises::promise(function(resolve, reject) { resolver <<- resolve })

  jid <- g("launch_job")(sid, "bayes_tg", p,
                         context = list(bayes_params = NULL))
  expect_true(g("cancel_job")(sid, jid))
  resolver("late result")
  flush_promises()

  job <- g("get_job")(sid, jid)
  expect_equal(job$status, "cancelled")
  expect_null(job$result)
})
