# Regression tests for step 6 of the Round 20/21 implementation (CODE_REVIEW.md
# R20-D.N): Bayesian robustness. The fits run on a real background worker
# (R20.D62), wait in a queue with one fit per browser and can be cancelled
# (R20.D7), Dismiss and a lost job no longer strand the modal (R20.D5), Resume
# delivers its result (R20.D63) to the right browser (R20.D64), and the
# diagnostics no longer need the fitted model (R20.D18). Each test fails on
# v26.09.005.

skip_if_not_installed("shiny")
skip_if_not_installed("promises")
skip_if_not_installed("future", minimum_version = "1.49.0")

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

# Every test drives the scheduler by hand: no `later` timer is left running.
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

task_of <- function(fun, ...) {
  environment(fun) <- globalenv()
  list(fun = fun, args = list(...))
}

sid_a <- strrep("a", 64L)
sid_b <- strrep("b", 64L)

# ---- R20.D62: a real background worker --------------------------------------------

test_that("R20.D62: Bayesian fits run in a background process, not the Shiny one", {
  skip_on_cran()
  local_registry()
  # What runDashboard() used: future runs a one-worker multisession plan in
  # the calling process, so every fit froze every session.
  future::plan(future::multisession, workers = 1L)
  expect_identical(future::value(future::future(Sys.getpid())), Sys.getpid())
  g("use_background_worker")()
  expect_false(identical(future::value(future::future(Sys.getpid())), Sys.getpid()))
})

# ---- R20.D7: the queue, one fit per browser, a cancel that stops the fit ----------

test_that("R20.D7: fits wait for the worker in arrival order, then run", {
  local_registry()
  future::plan(future::sequential)
  free <- 0L
  testthat::local_mocked_bindings(.free_workers = function() free,
                                  .package = "mouseExperimentDashboard")
  j1 <- g("launch_job")(sid_a, "bayes_tg", task = task_of(function(x) x + 1, x = 1))
  j2 <- g("launch_job")(sid_b, "bayes_surv", task = task_of(function(x) x * 10, x = 4))
  expect_identical(g("get_job")(sid_a, j1)$status, "queued")
  expect_identical(g("jobs_ahead")(sid_a, j1), 0L)
  expect_identical(g("jobs_ahead")(sid_b, j2), 1L)       # behind sid_a's fit
  expect_identical(g("active_job")(sid_a)$id, j1)
  expect_null(g("active_job")(strrep("c", 64L)))
  expect_match(g("job_status_text")(g("get_job")(sid_b, j2), 1L),
               "1 fit ahead.*one at a time")
  free <- 1L                                               # the worker frees up
  g(".pump_jobs")(schedule = FALSE)
  g(".pump_jobs")(schedule = FALSE)
  expect_identical(g("get_job")(sid_a, j1)$status, "completed")
  expect_equal(g("get_job")(sid_a, j1)$result, 2)
  expect_equal(g("get_job")(sid_b, j2)$result, 40)
})

test_that("R20.D7: a cancelled queued fit never starts", {
  local_registry()
  future::plan(future::sequential)
  free <- 0L
  testthat::local_mocked_bindings(.free_workers = function() free,
                                  .package = "mouseExperimentDashboard")
  j <- g("launch_job")(sid_a, "bayes_tg", task = task_of(function() "ran"))
  expect_true(g("cancel_job")(sid_a, j))
  free <- 1L
  g(".pump_jobs")(schedule = FALSE)
  expect_identical(g("get_job")(sid_a, j)$status, "cancelled")
  expect_null(g("get_job")(sid_a, j)$result)
})

test_that("R20.D7: cancelling a running fit interrupts it and frees the worker", {
  skip_on_cran()
  local_registry()
  g("use_background_worker")()
  j <- g("launch_job")(sid_a, "bayes_tg", task = task_of(function() { Sys.sleep(60); "late" }))
  expect_identical(g("get_job")(sid_a, j)$status, "running")
  t0 <- Sys.time()
  expect_true(g("cancel_job")(sid_a, j))
  while (future::nbrOfFreeWorkers() < 1L &&
         difftime(Sys.time(), t0, units = "secs") < 20) Sys.sleep(0.2)
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 20)  # not 60
  j2 <- g("launch_job")(sid_a, "bayes_tg", task = task_of(function() 42))
  while (!identical(g("get_job")(sid_a, j2)$status, "completed") &&
         difftime(Sys.time(), t0, units = "secs") < 40) {
    Sys.sleep(0.2)
    g(".pump_jobs")(schedule = FALSE)
  }
  expect_equal(g("get_job")(sid_a, j2)$result, 42)
})

test_that("R20.D7: a second Run from the same browser does not queue a second fit", {
  local_registry()
  testthat::local_mocked_bindings(.free_workers = function() 0L,
                                  .package = "mouseExperimentDashboard")
  rv <- shiny::reactiveValues(session_id = sid_a, tg_running = TRUE, tg_job_id = NULL)
  params <- list(data = data.frame(), col_names = list(), n_chains = 2L,
                 n_warmup = 100L, n_iter = 100L, backend = "rstan")
  shiny::testServer(function(id, rv) {
    shiny::moduleServer(id, function(input, output, session) {
      g("launch_tg_bayes_job")(sid_a, params, input, output, session, rv)
      shiny::isolate(rv$tg_running <- TRUE)
      g("launch_tg_bayes_job")(sid_a, params, input, output, session, rv)
    })
  }, args = list(rv = rv), {
    expect_length(g("list_jobs")(sid_a), 1L)               # was two full fits
    expect_false(shiny::isolate(rv$tg_running))
  })
  # The task asks for no model back (R20.D18).
  expect_false(g(".tg_bayes_task")(params)$args$return_model)
  expect_false(g(".surv_bayes_task")(list())$args$return_model)
})

# ---- R20.D5: Dismiss and a lost job --------------------------------------------------

test_that("R20.D5: Dismiss clears finished jobs only, and a lost job closes its modal", {
  local_registry()
  done <- g("launch_job")(sid_a, "bayes_tg",
                          promises::promise(function(resolve, reject) resolve("ok")))
  later::run_now(timeoutSecs = 1)
  live <- g("launch_job")(sid_a, "bayes_surv", promises::promise(function(resolve, reject) NULL))
  shiny::testServer(g("resume_banner_server"),
                    args = list(session_id = shiny::reactive(sid_a),
                                on_resume = function(job) NULL), {
    session$setInputs(dismiss_clicked = list(n = 1))
    expect_null(g("get_job")(sid_a, done))
    expect_identical(g("get_job")(sid_a, live)$status, "running")  # was retired
  })

  rv <- shiny::reactiveValues(session_id = sid_a, tg_job_id = "job_gone",
                              tg_running = TRUE)
  shiny::testServer(function(id, rv) {
    shiny::moduleServer(id, function(input, output, session) {
      g("attach_bayes_poller")(input, output, session, rv, "tg_job_id",
                               g(".finalize_tg_success"), g(".finalize_tg_failure"),
                               g(".finalize_tg_cancellation"), g(".finalize_tg_lost"))
    })
  }, args = list(rv = rv), {
    session$flushReact()
    expect_null(shiny::isolate(rv$tg_job_id))
    expect_false(shiny::isolate(rv$tg_running))           # the modal was left open
    expect_match(paste(as.character(output$tg_status_msg$html), collapse = " "),
                 "dismissed or expired")
  })
})

# ---- R20.D63 / R20.D64: Resume delivers to the right browser -------------------------

test_that("R20.D63: Resume delivers a finished fit's result", {
  local_registry()
  fake <- list(treatment_effects = data.frame(Group = "A", Adjusted_Mean = 1),
               error = FALSE)
  jid <- g("launch_job")(sid_a, "bayes_tg",
                         promises::promise(function(resolve, reject) resolve(fake)),
                         context = list(bayes_params = list(
                           data = data.frame(Day = 1), col_names = list())))
  later::run_now(timeoutSecs = 1)
  app <- function(input, output, session) {
    rv <- shiny::reactiveValues(session_id = sid_a, tg_running = FALSE,
                                tg_results = NULL, tg_job_id = NULL)
    session$userData$rv <- rv
    # The handler is registered in the module, the click dispatched from the
    # app's session (their userData is shared).
    shiny::moduleServer("tg", function(input, output, session) {
      g("register_tg_resume")(input, output, session, rv)
    })
    g("resume_banner_server")("rb", session_id = shiny::reactive(sid_a),
                              on_resume = function(job) g("dispatch_resume")(job, session, rv))
  }
  shiny::testServer(app, {
    session$setInputs(`rb-resume_clicked` = list(jid = jid, n = 1))
    session$flushReact()
    rv <- session$userData$rv
    expect_equal(shiny::isolate(rv$tg_results)$treatment_effects$Group, "A")  # was lost
    expect_null(g("get_job")(sid_a, jid))                  # retired after delivery
  })
})

test_that("R20.D64: a resume handler belongs to its own browser session", {
  a <- list(userData = new.env())
  b <- list(userData = new.env())
  got <- character(0)
  g("register_resume_handler")("bayes_tg", function(job, s, rv) got <<- c(got, "A"),
                               session = a)
  g("register_resume_handler")("bayes_tg", function(job, s, rv) got <<- c(got, "B"),
                               session = b)
  g("dispatch_resume")(list(kind = "bayes_tg", id = "j"), a, NULL)
  expect_identical(got, "A")                               # the later session's handler ran
})

# ---- R20.D18: diagnostics without the fitted model ----------------------------------

test_that("R20.D18: the Bayesian diagnostics draw from the result, not a kept model", {
  skip_if_not_installed("posterior")
  res <- list(model_type_used = "bayes", error = FALSE,
              pp_check_plot   = ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                                                ggplot2::aes(x, y)) + ggplot2::geom_point(),
              posterior_draws = posterior::example_draws(),
              loo_diagnostics = data.frame(k_threshold = 0.64, n_high_k = 1L,
                                           pareto_k = I(list(c(0.66, 0.2)))))
  rv <- shiny::reactiveValues(tg_results = res, data_valid = TRUE, settings = list())
  shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    expect_true(nzchar(output$tg_pp_check_plot$src))
    expect_true(nzchar(output$tg_mcmc_trace_plot$src))
    expect_true(nzchar(output$tg_mcmc_rank_plot$src))
  })
  # The backend's draw-dependent Pareto-k threshold, not a fixed 0.7.
  expect_equal(nrow(g("bayes_influential_obs")(res)), 1L)
})

test_that("R20.D7: the progress window's Cancel stops this tab's fit", {
  local_registry()
  testthat::local_mocked_bindings(.free_workers = function() 0L,
                                  .package = "mouseExperimentDashboard")
  jid <- g("launch_job")(sid_a, "bayes_tg", task = task_of(function() "ran"))
  rv <- shiny::reactiveValues(session_id = sid_a, tg_job_id = jid, tg_running = TRUE)
  shiny::testServer(function(id, rv) {
    shiny::moduleServer(id, function(input, output, session) {
      g("register_tg_resume")(input, output, session, rv)
    })
  }, args = list(rv = rv), {
    session$setInputs(tg_cancel_bayes = 1)     # the button in the progress window
    expect_identical(g("get_job")(sid_a, jid)$status, "cancelled")
  })
})
