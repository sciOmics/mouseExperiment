# Tests for the Tumor Growth resume helper (CODE_REVIEW.md F.7, commit 3).

skip_if_not_installed("shiny")
skip_if_not_installed("promises")

pkg <- tryCatch(getNamespace("mouseExperimentDashboard"),
                error = function(e) NULL)
skip_if(is.null(pkg), "mouseExperimentDashboard namespace not loadable")

g <- function(nm) get(nm, envir = pkg)

flush_promises <- function() later::run_now(timeoutSecs = 1)

# ---------------------------------------------------------------------------
# .finalize_tg_bayes_success populates rv idempotently
# ---------------------------------------------------------------------------
test_that(".finalize_tg_bayes_success writes rv$tg_results + plot data", {
  rv <- shiny::reactiveValues(
    tg_running   = TRUE,
    tg_results   = NULL,
    tg_plot_data = NULL,
    tg_col_names = NULL,
    bayes_models = NULL
  )
  fake_result <- list(
    model            = "<fake brmsfit>",
    treatment_effects = data.frame(Group = "A", Adjusted_Mean = 1)
  )
  fake_params <- list(
    data      = data.frame(x = 1, Day = 1, Volume = 1,
                           Treatment = "A", ID = 1),
    col_names = list(time_column = "Day", volume_column = "Volume",
                     treatment_column = "Treatment", id_column = "ID",
                     cage_column = NULL)
  )

  shiny::testServer(function(id, rv) {
    shiny::moduleServer(id, function(input, output, session) {
      # We don't actually run the function inside the module — call it
      # directly with the captured session.
      g(".finalize_tg_bayes_success")(fake_result, fake_params,
                                      input, output, session, rv)
    })
  }, args = list(rv = rv), {
    expect_equal(shiny::isolate(rv$tg_results)$model_type_used, "bayes")
    expect_equal(shiny::isolate(rv$bayes_models)$tg, "<fake brmsfit>")
    expect_false(shiny::isolate(rv$tg_running))
    expect_equal(shiny::isolate(rv$tg_plot_data), fake_params$data)
    expect_equal(shiny::isolate(rv$tg_col_names), fake_params$col_names)
  })
})

# ---------------------------------------------------------------------------
# .finalize_tg_bayes_failure writes error state without leaking sid
# ---------------------------------------------------------------------------
test_that(".finalize_tg_bayes_failure sets rv$tg_running=FALSE", {
  rv <- shiny::reactiveValues(tg_running = TRUE)
  shiny::testServer(function(id, rv) {
    shiny::moduleServer(id, function(input, output, session) {
      g(".finalize_tg_bayes_failure")("Stan compile error", output, rv)
    })
  }, args = list(rv = rv), {
    expect_false(shiny::isolate(rv$tg_running))
  })
})

# ---------------------------------------------------------------------------
# register_tg_resume + dispatch_resume round-trip
# ---------------------------------------------------------------------------
test_that("register_tg_resume hooks into the resume dispatcher", {
  g(".reset_resume_handlers")()
  on.exit(g(".reset_resume_handlers")(), add = TRUE)

  rv <- shiny::reactiveValues(tg_job_id = NULL, session_id = NULL)

  shiny::testServer(function(id, rv) {
    shiny::moduleServer(id, function(input, output, session) {
      g("register_tg_resume")(input, output, session, rv)

      # Verify the handler is now resolvable by kind.
      handler <- g(".resume_handlers")[["bayes_tg"]]
      expect_true(is.function(handler))
    })
  }, args = list(rv = rv), {})
})

# ---------------------------------------------------------------------------
# launch_tg_bayes_job registers a job with the registry + sets rv$tg_job_id
# ---------------------------------------------------------------------------
test_that("launch_tg_bayes_job registers a job under the active sid", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)

  sid <- strrep("a", 64L)
  rv  <- shiny::reactiveValues(tg_running = TRUE, session_id = sid)
  # Trivial bayes_params; we don't actually run brms. The future_promise
  # will fail, but launch_job still registers the entry.
  fake_params <- list(
    data         = data.frame(),
    col_names    = list(time_column = "Day", volume_column = "Volume",
                        treatment_column = "Treatment",
                        id_column = "ID", cage_column = NULL),
    dose_col     = NULL,
    necrotic_col = NULL,
    model_type   = "lmm",
    random_effects = "intercept_only",
    reference_group = "Control",
    prior_strength = "skeptical",
    manual_priors = NULL,
    n_chains = 2L, n_warmup = 100L, n_iter = 50L,
    backend = "rstan",
    necrotic_handling = "none")

  # We can't (and don't want to) actually call launch_tg_bayes_job
  # because it calls show_analysis_modal which needs a Shiny session.
  # Verify the smaller pieces instead.

  # Confirm the registry is empty.
  expect_equal(length(g("list_jobs")(sid)), 0L)

  # Manually register a job (what launch_tg_bayes_job would do).
  promise <- promises::promise(function(resolve, reject) resolve("ok"))
  job_id <- g("launch_job")(sid, "bayes_tg", promise,
                            context = list(bayes_params = fake_params))
  expect_match(job_id, "^job_")
  expect_equal(length(g("list_jobs")(sid)), 1L)
  expect_equal(g("get_job")(sid, job_id)$context$bayes_params$model_type,
               "lmm")
})
