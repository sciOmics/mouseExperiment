# Regression tests for step 7 of the Round 20/21 implementation: the VPS
# Round 3 items that live in the dashboard. In the container the proxy-timeout
# reminder fires and no browser is launched (R20.D28 / VPS V3.22), and the
# Bayesian worker, not the Shiny process, is the kernel's preferred OOM victim
# (VPS V3.12). The first three tests fail on v26.09.006; the last checks that
# a job still returns its value through the change.

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

test_that("R20.D28 / V3.22: in the container, runDashboard() reminds about the proxy and opens no browser", {
  # What the production container runs: R -e "runDashboard(host = '0.0.0.0')",
  # with neither SHINY_PORT nor SHINY_SERVER_VERSION set.
  withr::local_envvar(SHINY_PORT = NA, SHINY_SERVER_VERSION = NA)
  got <- NULL
  testthat::local_mocked_bindings(
    runApp = function(appDir, ...) { got <<- list(...); invisible() },
    .package = "shiny")
  testthat::local_mocked_bindings(
    dashboard_app = function() structure(list(), class = "shiny.appobj"),
    use_background_worker = function() invisible(),
    .package = "mouseExperimentDashboard")
  msgs <- testthat::capture_messages(g("runDashboard")(host = "0.0.0.0", port = 3838))
  expect_true(any(grepl("proxy_read_timeout >= 3600s", msgs, fixed = TRUE)))  # was silent
  expect_false(got$launch.browser)                                              # was TRUE
  expect_identical(got$host, "0.0.0.0")
})

test_that("R20.D28: a local run stays quiet about the proxy", {
  withr::local_envvar(SHINY_PORT = NA, SHINY_SERVER_VERSION = NA)
  expect_silent(g("startup_proxy_reminder")("127.0.0.1"))
  withr::local_envvar(SHINY_PORT = "3838")                    # under shiny-server
  expect_message(g("startup_proxy_reminder")("127.0.0.1"), "reverse proxy")
})

test_that("V3.12: a background worker raises its own OOM score; the Shiny process never does", {
  f <- withr::local_tempfile(lines = "0")
  # A worker: another process than the one that launched the job.
  expect_true(g(".prefer_oom_kill")(main_pid = -1L, path = f))
  expect_identical(readLines(f), "500")
  # A sequential plan runs the task in the Shiny process: left alone.
  writeLines("0", f)
  expect_false(g(".prefer_oom_kill")(main_pid = Sys.getpid(), path = f))
  expect_identical(readLines(f), "0")
  # No procfs (macOS, Windows): nothing to do, and no error.
  expect_false(g(".prefer_oom_kill")(main_pid = -1L, path = file.path(tempdir(), "absent")))
})

test_that("V3.12: a job still runs, and returns its value, through the OOM preference", {
  old_opt <- options(mE.job_pump = FALSE)
  old_plan <- future::plan(future::sequential)
  g(".reset_job_registry")()
  withr::defer({ g(".reset_job_registry")(); future::plan(old_plan); options(old_opt) })
  testthat::local_mocked_bindings(.free_workers = function() 1L,
                                  .package = "mouseExperimentDashboard")
  fun <- function(x) x * 2
  environment(fun) <- globalenv()
  sid <- strrep("a", 64L)
  j <- g("launch_job")(sid, "bayes_tg", task = list(fun = fun, args = list(x = 21)))
  g(".pump_jobs")(schedule = FALSE)
  expect_identical(g("get_job")(sid, j)$status, "completed")
  expect_equal(g("get_job")(sid, j)$result, 42)
})
