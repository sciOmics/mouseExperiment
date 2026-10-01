#' Process-level job registry for resumable Bayesian fits (CODE_REVIEW.md F.7).
#'
#' Holds the long-running fits *outside* the per-session reactive scope so a
#' fit survives tab disconnect / NPM proxy timeout / browser sleep. The
#' registry is a package-level environment keyed by a session-identity
#' cookie sid; lookups are always sid-scoped so one user can't enumerate or
#' read another's jobs.
#'
#' API:
#'   launch_job(sid, kind, task = , context = )  queue a fit; returns job_id
#'   launch_job(sid, kind, promise, context)     register a promise (tests)
#'   get_job(sid, job_id)                        one job or NULL
#'   list_jobs(sid)                              named list of jobs for that sid
#'   active_job(sid)                             the sid's queued or running job
#'   jobs_ahead(sid, job_id)                     fits that will run before it
#'   cancel_job(sid, job_id)                     drop a queued job; interrupt a running one
#'   evict_stale_jobs(max_done_age = 1800, max_run_age = 14400)
#'
#' Scheduling (CODE_REVIEW.md R20.D7, R20.D62). A fit is a `task`, a function
#' and its arguments. It waits as "queued" until the future backend has a free
#' worker, then runs as a future ("running") and ends "completed" or "failed".
#' `.pump_jobs()` starts queued jobs in arrival order and collects finished
#' ones; it runs on a `later` timer while any job is pending, and on every
#' poller tick. Because the registry holds the Future object, cancelling a
#' running fit interrupts its worker (`future::cancel()`); the dashboard used
#' `promises::future_promise()`, which hides the future, so a cancelled fit
#' kept the worker busy until it finished.
#'
#' Isolation guarantees (CODE_REVIEW.md F.7, layer 1):
#'   * Every read API takes sid as the first argument. There is no
#'     list_all_jobs / get_job(job_id_only).
#'   * Job IDs are never access tokens — they're indices within an sid's
#'     namespace. Two sids could hold the same job_id and still see
#'     different jobs.
#'   * The registry env is package-private (the dot-prefix). No exported
#'     accessor returns the env itself. `jobs_ahead()` reveals only how many
#'     other fits are waiting, not whose.

# Package-private environment. Initialized empty; populated as users
# launch jobs. Survives across Shiny sessions because it's in the R
# process, not in any session's reactiveValues.
.job_registry <- new.env(parent = emptyenv())
.job_registry$sessions <- list()

# ---- Helpers (package-private) --------------------------------------------

.new_job_id <- function() {
  # 16 bytes of randomness; collision probability across the registry's
  # lifetime is vanishing. Hex-encoded for grep-friendly logs.
  paste0("job_", paste0(format(as.hexmode(sample.int(255L, 16L, replace = TRUE)),
                                width = 2L),
                        collapse = ""))
}

.registry_put <- function(sid, job_id, job) {
  if (is.null(.job_registry$sessions[[sid]])) {
    .job_registry$sessions[[sid]] <- list()
  }
  .job_registry$sessions[[sid]][[job_id]] <- job
  invisible(NULL)
}

.registry_set <- function(sid, job_id, field, value) {
  if (is.null(.job_registry$sessions[[sid]])) return(invisible(NULL))
  if (is.null(.job_registry$sessions[[sid]][[job_id]])) return(invisible(NULL))
  .job_registry$sessions[[sid]][[job_id]][[field]] <- value
  invisible(NULL)
}

# ---- Public-within-package API --------------------------------------------

#' Launch a job: queue a task, or register a promise.
#'
#' A `task` is `list(fun = <function>, args = <list>)`; it is queued and run
#' as a future when a worker is free (see `.pump_jobs()`). A `promise` is
#' registered as already running, with then() handlers writing its outcome
#' back; the tests use it with resolved promises.
#'
#' Race note: if `cancel_job()` runs before a promise resolves, the
#' onFulfilled handler detects status == "cancelled" and discards the
#' result instead of overwriting the cancellation.
#'
#' @param sid Session-identity cookie value. Must be non-NULL non-empty.
#' @param kind String tag for the job's analysis kind ("bayes_tg" or
#'   "bayes_surv"). Used by the banner to label and route resumes.
#' @param promise A `promises::promise` object, or NULL when `task` is given.
#' @param context Optional named list of metadata the resume handler
#'   will need (e.g. `bayes_params` for re-populating plot data). Not
#'   inspected by the registry; passed through to the handler verbatim.
#' @param task `list(fun, args)` to run as a future, or NULL.
#' @return The new job_id (character).
#' @keywords internal
launch_job <- function(sid, kind, promise = NULL, context = NULL, task = NULL) {
  stopifnot(is.character(sid), length(sid) == 1L, nzchar(sid))
  stopifnot(is.character(kind), length(kind) == 1L, nzchar(kind))
  if (is.null(task) == is.null(promise)) {
    stop("Supply exactly one of `task` or `promise`.", call. = FALSE)
  }
  if (!is.null(promise) && !inherits(promise, "promise")) {
    stop("`promise` must be a promises::promise object.", call. = FALSE)
  }
  if (!is.null(task) && (!is.function(task$fun) || !is.list(task$args))) {
    stop("`task` must be list(fun = <function>, args = <list>).", call. = FALSE)
  }
  job_id <- .new_job_id()
  now <- Sys.time()
  job <- list(
    id           = job_id,
    kind         = kind,
    sid          = sid,
    queued_at    = now,
    started_at   = if (is.null(task)) now else NA,
    completed_at = NA,
    status       = if (is.null(task)) "running" else "queued",
    result       = NULL,
    error        = NULL,
    context      = context,
    task         = task,
    future       = NULL
  )
  .registry_put(sid, job_id, job)

  if (!is.null(promise)) {
    # Use promises::then() rather than %...>% so we don't leave a dangling
    # reference in the caller's scope.
    promises::then(promise,
      onFulfilled = function(value) {
        cur <- get_job(sid, job_id)
        if (is.null(cur) || cur$status == "cancelled") return(invisible(NULL))
        .registry_set(sid, job_id, "status",       "completed")
        .registry_set(sid, job_id, "result",       value)
        .registry_set(sid, job_id, "completed_at", Sys.time())
      },
      onRejected = function(err) {
        cur <- get_job(sid, job_id)
        if (is.null(cur) || cur$status == "cancelled") return(invisible(NULL))
        msg <- tryCatch(conditionMessage(err), error = function(e) as.character(err))
        .registry_set(sid, job_id, "status",       "failed")
        .registry_set(sid, job_id, "error",        msg)
        .registry_set(sid, job_id, "completed_at", Sys.time())
      })
  } else {
    .pump_jobs()
  }
  job_id
}

#' The sid's queued or running job, or NULL
#'
#' The dashboard runs one Bayesian fit per browser at a time (R20.D7): a
#' second Run click used to queue a second full fit behind the first.
#' @keywords internal
active_job <- function(sid) {
  jobs <- list_jobs(sid)
  for (job in jobs) if (job$status %in% c("queued", "running")) return(job)
  NULL
}

#' How many fits will run before this one
#'
#' The running fits plus the jobs queued earlier, across every sid; the count
#' says how long the queue is without saying whose fits are in it.
#' @keywords internal
jobs_ahead <- function(sid, job_id) {
  me <- get_job(sid, job_id)
  if (is.null(me) || !identical(me$status, "queued")) return(0L)
  n <- 0L
  for (bucket in .job_registry$sessions) for (job in bucket) {
    if (identical(job$status, "running")) n <- n + 1L
    else if (identical(job$status, "queued") &&
             as.numeric(job$queued_at) < as.numeric(me$queued_at)) n <- n + 1L
  }
  n
}

# ---- Scheduler (package-private; CODE_REVIEW.md R20.D7) --------------------

# Interrupt a running future and release its worker. future counts a
# cancelled future's worker as busy until the future is queried, so it is
# queried here; otherwise the next fit waited for the cancelled one anyway.
.cancel_future <- function(f) {
  if (!inherits(f, "Future")) return(invisible(NULL))
  tryCatch(future::cancel(f, interrupt = TRUE), error = function(e) NULL)
  tryCatch(future::resolved(f), error = function(e) NULL)
  tryCatch(future::value(f), error = function(e) NULL)
  invisible(NULL)
}

.free_workers <- function() {
  n <- tryCatch(future::nbrOfFreeWorkers(), error = function(e) 0L)
  if (length(n) != 1L || is.na(n)) 0L else n
}

# Jobs the pump has work for: queued ones, and running futures. A
# promise-based job resolves through its own then() handlers.
.jobs_pending <- function() {
  for (bucket in .job_registry$sessions) for (job in bucket) {
    if (identical(job$status, "queued") ||
        (identical(job$status, "running") && inherits(job$future, "Future"))) {
      return(TRUE)
    }
  }
  FALSE
}

.next_queued_job <- function() {
  best <- NULL
  for (bucket in .job_registry$sessions) for (job in bucket) {
    if (identical(job$status, "queued") &&
        (is.null(best) || as.numeric(job$queued_at) < as.numeric(best$queued_at))) {
      best <- job
    }
  }
  best
}

#' Make this process the kernel's preferred OOM victim
#'
#' VPS CODE_REVIEW.md V3.12: a Bayesian fit's model compile can take 2.5 GB.
#' When the container reaches its memory limit, the kernel should kill the
#' background worker (or the compiler it started, which inherits the setting),
#' and the fit fails with an error, rather than the Shiny process that serves
#' every session. Raising one's own `oom_score_adj` needs no privilege. Never
#' applied to the process that launched the job: with a sequential plan the
#' task runs in the Shiny process itself.
#'
#' @param main_pid PID of the process that launched the job.
#' @param path The procfs file (a parameter for tests).
#' @return `TRUE` (invisibly) when the score was raised.
#' @keywords internal
.prefer_oom_kill <- function(main_pid, path = "/proc/self/oom_score_adj") {
  if (Sys.getpid() == main_pid || !file.exists(path)) return(invisible(FALSE))
  invisible(!inherits(try(writeLines("500", path), silent = TRUE), "try-error"))
}

.start_job <- function(job) {
  task <- job$task
  # Base functions only, so the worker needs no namespace to run it.
  prefer_oom_kill <- .prefer_oom_kill
  environment(prefer_oom_kill) <- baseenv()
  f <- tryCatch(
    future::future(quote({ prefer_oom_kill(main_pid); do.call(fun, args) }),
                   substitute = FALSE,
                   globals  = list(fun = task$fun, args = task$args,
                                   prefer_oom_kill = prefer_oom_kill,
                                   main_pid = Sys.getpid()),
                   packages = task$packages,
                   seed     = TRUE),
    error = function(e) e)
  if (inherits(f, "error")) {
    .registry_set(job$sid, job$id, "status",       "failed")
    .registry_set(job$sid, job$id, "error",        conditionMessage(f))
    .registry_set(job$sid, job$id, "completed_at", Sys.time())
  } else {
    .registry_set(job$sid, job$id, "status",     "running")
    .registry_set(job$sid, job$id, "started_at", Sys.time())
    .registry_set(job$sid, job$id, "future",     f)
  }
  .registry_set(job$sid, job$id, "task", NULL)
  invisible(NULL)
}

.collect_job <- function(job) {
  f <- job$future
  if (!inherits(f, "Future")) return(invisible(NULL))
  if (!isTRUE(tryCatch(future::resolved(f), error = function(e) TRUE))) {
    return(invisible(NULL))
  }
  out <- tryCatch(
    list(ok = TRUE, value = suppressWarnings(suppressMessages(future::value(f)))),
    error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  cur <- get_job(job$sid, job$id)
  if (is.null(cur) || !identical(cur$status, "running")) return(invisible(NULL))
  if (isTRUE(out$ok)) {
    .registry_set(job$sid, job$id, "status", "completed")
    .registry_set(job$sid, job$id, "result", out$value)
  } else {
    .registry_set(job$sid, job$id, "status", "failed")
    .registry_set(job$sid, job$id, "error",  out$msg)
  }
  .registry_set(job$sid, job$id, "completed_at", Sys.time())
  .registry_set(job$sid, job$id, "future", NULL)
  invisible(NULL)
}

#' Start queued jobs and collect finished ones
#'
#' Idempotent and cheap; safe to call from every poller tick. With
#' `schedule = TRUE` it re-arms itself on a `later` timer while any job is
#' queued or running (option `mE.job_pump`, default TRUE; the tests turn it
#' off and call it directly).
#' @keywords internal
.pump_jobs <- function(schedule = getOption("mE.job_pump", TRUE)) {
  for (bucket in .job_registry$sessions) for (job in bucket) {
    if (identical(job$status, "running")) .collect_job(job)
  }
  repeat {
    nxt <- .next_queued_job()
    if (is.null(nxt) || .free_workers() < 1L) break
    .start_job(nxt)
  }
  if (isTRUE(schedule) && .jobs_pending() && is.null(.job_registry$pump) &&
      requireNamespace("later", quietly = TRUE)) {
    .job_registry$pump <- later::later(function() {
      .job_registry$pump <- NULL
      .pump_jobs()
    }, delay = 1)
  }
  invisible(NULL)
}

#' Retrieve a single job by sid + job_id.
#'
#' Returns `NULL` if the (sid, job_id) pair doesn't exist — including
#' when a malicious caller supplies a job_id that belongs to a different
#' sid. This is the isolation primitive.
#'
#' @param sid    Session-identity cookie value.
#' @param job_id Job ID returned by [launch_job()].
#' @return One job entry (list) or `NULL`.
#' @keywords internal
get_job <- function(sid, job_id) {
  if (!is.character(sid) || length(sid) != 1L || !nzchar(sid)) return(NULL)
  if (!is.character(job_id) || length(job_id) != 1L || !nzchar(job_id)) return(NULL)
  bucket <- .job_registry$sessions[[sid]]
  if (is.null(bucket)) return(NULL)
  bucket[[job_id]]
}

#' List all jobs belonging to a single sid.
#'
#' Returns an empty list when the sid is unknown or has no jobs. Never
#' returns jobs from any other sid.
#'
#' @param sid Session-identity cookie value.
#' @return Named list of jobs (possibly empty).
#' @keywords internal
list_jobs <- function(sid) {
  if (!is.character(sid) || length(sid) != 1L || !nzchar(sid)) return(list())
  bucket <- .job_registry$sessions[[sid]]
  if (is.null(bucket)) list() else bucket
}

#' Cancel a job.
#'
#' A queued job is dropped before it starts. A running job's future is
#' cancelled and its worker interrupted (future >= 1.49), so the worker is
#' free for the next fit; a promise-based job's late result is discarded by
#' its onFulfilled handler. Idempotent: cancelling a completed, failed or
#' cancelled job is a no-op.
#'
#' @param sid    Session-identity cookie value.
#' @param job_id Job ID returned by [launch_job()].
#' @return `TRUE` when the job's status was changed; `FALSE` when the
#'   job didn't exist or was already terminal.
#' @keywords internal
cancel_job <- function(sid, job_id) {
  cur <- get_job(sid, job_id)
  if (is.null(cur)) return(FALSE)
  if (cur$status %in% c("completed", "failed", "cancelled")) return(FALSE)
  if (inherits(cur$future, "Future")) {
    .cancel_future(cur$future)
  }
  .registry_set(sid, job_id, "status",       "cancelled")
  .registry_set(sid, job_id, "completed_at", Sys.time())
  .registry_set(sid, job_id, "future",       NULL)
  .registry_set(sid, job_id, "task",         NULL)
  .pump_jobs()
  TRUE
}

#' Remove a job from the registry.
#'
#' Called by a poller once it has delivered a job's outcome, and by the
#' banner's "Dismiss finished". A job that is still queued or running is
#' cancelled first, so a removed job cannot keep the worker busy.
#' @keywords internal
retire_job <- function(sid, job_id) {
  cur <- get_job(sid, job_id)
  if (is.null(cur)) return(FALSE)
  if (cur$status %in% c("queued", "running")) cancel_job(sid, job_id)
  bucket <- .job_registry$sessions[[sid]]
  bucket[[job_id]] <- NULL
  .job_registry$sessions[[sid]] <- bucket
  if (length(bucket) == 0L) {
    .job_registry$sessions[[sid]] <- NULL
  }
  TRUE
}

#' Evict stale jobs from the registry.
#'
#' Completed / failed / cancelled jobs older than `max_done_age` seconds
#' are dropped. Running jobs older than `max_run_age` are also dropped
#' (assumed abandoned — the user never came back). Called periodically
#' by the TTL sweeper.
#'
#' @param max_done_age Seconds to retain completed / failed / cancelled
#'   jobs after they finished. Default 1800 (30 min).
#' @param max_run_age  Seconds to retain still-running jobs. Default
#'   14400 (4 hours). Should comfortably exceed the longest expected
#'   Bayesian fit time.
#' @return Number of jobs evicted.
#' @keywords internal
evict_stale_jobs <- function(max_done_age = 1800, max_run_age = 14400) {
  now <- Sys.time()
  n_evicted <- 0L
  for (sid in names(.job_registry$sessions)) {
    bucket <- .job_registry$sessions[[sid]]
    keep <- vapply(bucket, function(job) {
      if (job$status %in% c("completed", "failed", "cancelled")) {
        # Use queued_at if completed_at is NA (shouldn't happen, but be safe)
        end_t <- if (is.na(job$completed_at)) job$queued_at else job$completed_at
        as.numeric(difftime(now, end_t, units = "secs")) < max_done_age
      } else {
        # queued or running — judge by when it was queued
        as.numeric(difftime(now, job$queued_at, units = "secs")) < max_run_age
      }
    }, logical(1L))
    # An evicted fit that is still running would keep the worker busy.
    for (job in bucket[!keep]) {
      if (inherits(job$future, "Future")) {
        .cancel_future(job$future)
      }
    }
    bucket <- bucket[keep]
    n_evicted <- n_evicted + sum(!keep)
    if (length(bucket) == 0L) {
      .job_registry$sessions[[sid]] <- NULL
    } else {
      .job_registry$sessions[[sid]] <- bucket
    }
  }
  n_evicted
}

#' Reset the registry. Test-only.
#' @keywords internal
.reset_job_registry <- function() {
  for (bucket in .job_registry$sessions) for (job in bucket) {
    if (inherits(job$future, "Future")) {
      .cancel_future(job$future)
    }
  }
  if (is.function(.job_registry$pump)) .job_registry$pump()
  .job_registry$pump     <- NULL
  .job_registry$sessions <- list()
  invisible(NULL)
}
