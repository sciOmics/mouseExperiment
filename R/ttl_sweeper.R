#' Periodic registry TTL sweeper (CODE_REVIEW.md F.7).
#'
#' Runs `evict_stale_jobs()` on a recurring schedule via `later::later`
#' so the package-level `.job_registry` doesn't grow without bound. Each
#' tick:
#'   * drops completed / failed / cancelled jobs older than
#'     `max_done_age` (default 30 minutes),
#'   * drops still-running jobs older than `max_run_age` (default
#'     4 hours — comfortably exceeds the longest expected Bayesian fit).
#'
#' The sweeper is process-scoped: one instance per R process is enough,
#' shared across every Shiny session. `start_ttl_sweeper()` is
#' idempotent — calling it from `runDashboard()` twice doesn't double up.

# Track whether the sweeper has been started so repeat
# `start_ttl_sweeper()` calls are no-ops.
.sweeper_state <- new.env(parent = emptyenv())
.sweeper_state$started <- FALSE

#' Start the registry TTL sweeper.
#'
#' Idempotent. Should be called once per R process when the dashboard
#' app boots (from `runDashboard()`).
#'
#' @param interval_seconds Seconds between sweeps. Default 60.
#' @param max_done_age,max_run_age Passed through to
#'   [evict_stale_jobs()].
#' @return `invisible(NULL)`.
#' @keywords internal
start_ttl_sweeper <- function(interval_seconds = 60L,
                              max_done_age     = 1800L,
                              max_run_age      = 14400L) {
  if (isTRUE(.sweeper_state$started)) return(invisible(NULL))
  if (!requireNamespace("later", quietly = TRUE)) {
    warning("Package 'later' not available; TTL sweeper not started. ",
            "The job registry will grow over time.", call. = FALSE)
    return(invisible(NULL))
  }
  .sweeper_state$started <- TRUE
  do_sweep <- function() {
    tryCatch(
      evict_stale_jobs(max_done_age = max_done_age,
                       max_run_age  = max_run_age),
      error = function(e) {
        # Don't kill the recurring task on a transient error; just
        # log and continue. The next tick tries again.
        tryCatch(log_warn("ttl_sweep_error | %s", conditionMessage(e)),
                 error = function(e2) NULL)
        0L
      })
    later::later(do_sweep, delay = interval_seconds)
  }
  later::later(do_sweep, delay = interval_seconds)
  invisible(NULL)
}

#' Reset the sweeper state. Test-only.
#' @keywords internal
.reset_ttl_sweeper <- function() {
  .sweeper_state$started <- FALSE
  invisible(NULL)
}
