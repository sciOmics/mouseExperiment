#' Process-level job registry for resumable Bayesian fits (CODE_REVIEW.md F.7).
#'
#' Holds references to long-running futures *outside* the per-session
#' reactive scope so a fit survives tab disconnect / NPM proxy timeout /
#' browser sleep. The registry is a package-level environment keyed by
#' a session-identity cookie sid; lookups are always sid-scoped so one
#' user can't enumerate or read another's jobs.
#'
#' API:
#'   launch_job(sid, kind, promise)       returns job_id; wires then() handlers
#'   get_job(sid, job_id)                  one job or NULL
#'   list_jobs(sid)                        named list of jobs for that sid
#'   cancel_job(sid, job_id)               marks status = "cancelled"
#'   evict_stale_jobs(max_done_age = 1800, max_run_age = 14400)
#'
#' Lifecycle: every job starts with status = "running"; the promise's
#' fulfilled / rejected callbacks rewrite the status to "completed" or
#' "failed" and stash the result / error in the registry. A user-driven
#' cancel sets status = "cancelled"; the original future continues
#' running in its worker (R futures don't support clean cancellation) but
#' its result is dropped on resolution rather than written back.
#'
#' Isolation guarantees (CODE_REVIEW.md F.7, layer 1):
#'   * Every read API takes sid as the first argument. There is no
#'     list_all_jobs / get_job(job_id_only).
#'   * Job IDs are never access tokens — they're indices within an sid's
#'     namespace. Two sids could hold the same job_id and still see
#'     different jobs.
#'   * The registry env is package-private (the dot-prefix). No exported
#'     accessor returns the env itself.

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

#' Launch a job: register it and wire then() handlers to the registry.
#'
#' The caller constructs the `promise` itself (typically with
#' `promises::future_promise({...}, globals = list(...))`). This helper
#' attaches the registry-write handlers and stores the entry. Returns the
#' new job_id so the caller can store it in their reactiveValues for
#' polling.
#'
#' Race note: if `cancel_job()` runs before the future resolves, the
#' onFulfilled handler detects status == "cancelled" and discards the
#' result instead of overwriting the cancellation.
#'
#' @param sid Session-identity cookie value. Must be non-NULL non-empty.
#' @param kind String tag for the job's analysis kind ("bayes_tg",
#'   "bayes_bw", etc.). Used by the banner to label and route resumes.
#' @param promise A `promises::promise` object — already constructed by
#'   the caller.
#' @return The new job_id (character).
#' @keywords internal
launch_job <- function(sid, kind, promise) {
  stopifnot(is.character(sid), length(sid) == 1L, nzchar(sid))
  stopifnot(is.character(kind), length(kind) == 1L, nzchar(kind))
  if (!inherits(promise, "promise")) {
    stop("`promise` must be a promises::promise object.", call. = FALSE)
  }
  job_id <- .new_job_id()
  job <- list(
    id           = job_id,
    kind         = kind,
    sid          = sid,
    started_at   = Sys.time(),
    completed_at = NA,
    status       = "running",
    result       = NULL,
    error        = NULL
  )
  .registry_put(sid, job_id, job)

  # Wire callbacks. Use promises::then() rather than %...>% so we don't
  # leave a dangling reference in the caller's scope.
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

  job_id
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
#' Marks status = "cancelled". The underlying future continues to
#' completion in its worker (R futures don't support clean
#' cancellation), but its result is discarded by the onFulfilled
#' callback rather than written into the registry. Idempotent — calling
#' cancel on a completed/failed/cancelled job is a no-op.
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
  .registry_set(sid, job_id, "status",       "cancelled")
  .registry_set(sid, job_id, "completed_at", Sys.time())
  TRUE
}

#' Mark a job as retrieved (so the banner stops listing it).
#'
#' Called by the resume banner after a successful Resume + tab-switch.
#' Removes the entry from the registry to free memory.
#' @keywords internal
retire_job <- function(sid, job_id) {
  if (is.null(get_job(sid, job_id))) return(FALSE)
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
        # Use started_at if completed_at is NA (shouldn't happen, but be safe)
        end_t <- if (is.na(job$completed_at)) job$started_at else job$completed_at
        as.numeric(difftime(now, end_t, units = "secs")) < max_done_age
      } else {
        # status == "running" — judge by started_at
        as.numeric(difftime(now, job$started_at, units = "secs")) < max_run_age
      }
    }, logical(1L))
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
  .job_registry$sessions <- list()
  invisible(NULL)
}
