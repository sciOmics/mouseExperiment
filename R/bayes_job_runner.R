# Shared helpers for routing Bayesian fits through the job registry
# (CODE_REVIEW.md F.7).
#
# Per-module helpers (`tg_resume.R`, `surv_resume.R`) compose
# these two primitives — the polling observer and the resume
# registration — so each module's resume file stays small and the
# shared "poll → terminal status → call finalizer → retire → stop
# polling" logic lives in one place.

#' Attach the polling observer that consumes a single job.
#'
#' Polls `get_job(rv$session_id, rv[[job_id_field]])` every 2 s; on
#' terminal status, dispatches to the matching finalizer, retires the
#' registry entry, and clears the job-id field so polling stops. Each
#' module's `launch_<X>_bayes_job` and `register_<X>_resume` call this
#' with the module's finalizers + job_id_field name.
#'
#' Idempotent: re-attaching the observer for the same job_id is safe
#' (each observer sees the registry transition; only the first one
#' actually retires the entry).
#'
#' @param input,output,session,rv  Module-scoped Shiny objects.
#' @param job_id_field             String name of the rv slot holding
#'                                 the active job_id (e.g. "tg_job_id").
#' @param finalize_success         `function(result, context, input,
#'                                 output, session, rv)`.
#' @param finalize_failure         `function(err_msg, input, output,
#'                                 session, rv)`.
#' @param finalize_cancellation    `function(input, output, session,
#'                                 rv)`.
#' @param finalize_lost            `function(input, output, session, rv)`,
#'                                 called when the job has left the registry
#'                                 before its outcome was delivered
#'                                 (dismissed or evicted).
#' @keywords internal
attach_bayes_poller <- function(input, output, session, rv,
                                job_id_field,
                                finalize_success,
                                finalize_failure,
                                finalize_cancellation,
                                finalize_lost = NULL) {
  shiny::observe({
    job_id <- rv[[job_id_field]]
    sid    <- rv$session_id
    if (is.null(job_id) || is.null(sid) || identical(sid, SID_NONE)) return()
    shiny::invalidateLater(2000L, session)
    .pump_jobs()
    job <- get_job(sid, job_id)
    if (is.null(job)) {
      # R20.D5: the job was dismissed or evicted while this tab waited. It
      # used to return here and leave the modal open for good.
      rv[[job_id_field]] <- NULL
      if (is.function(finalize_lost)) finalize_lost(input, output, session, rv)
      return()
    }
    if (job$status %in% c("queued", "running")) {
      # R20.D7: the queue position or elapsed time, in the modal.
      try(shinyjs::html(id = "mE-job-status", asis = TRUE,
                        html = job_status_text(job, jobs_ahead(sid, job_id))),
          silent = TRUE)
      return()
    }
    if (job$status == "completed") {
      result <- job$result
      if (is.list(result) &&
          !is.null(job$started_at) && !is.null(job$completed_at)) {
        result$elapsed_seconds <- as.numeric(
          difftime(job$completed_at, job$started_at, units = "secs"))
      }
      finalize_success(result, job$context, input, output, session, rv)
    } else if (job$status == "failed") {
      finalize_failure(job$error, input, output, session, rv)
    } else if (job$status == "cancelled") {
      finalize_cancellation(input, output, session, rv)
    }
    retire_job(sid, job_id)
    rv[[job_id_field]] <- NULL
  })
}

#' Cancel this module's fit from the progress window's button
#'
#' @param input,rv Module scope.
#' @param cancel_input The button's (unnamespaced) input id.
#' @param job_id_field Name of the rv slot holding the active job id.
#' @keywords internal
observe_bayes_cancel <- function(input, rv, cancel_input, job_id_field) {
  shiny::observeEvent(input[[cancel_input]], {
    sid <- rv$session_id
    jid <- rv[[job_id_field]]
    if (!is.null(sid) && !is.null(jid)) cancel_job(sid, jid)
  })   # an action button's initial 0 counts as NULL, so it does not fire
}

#' Register a per-kind resume handler.
#'
#' Wraps `register_resume_handler()` with the module-scoped poller
#' attach and the tab switch. Captures the module's input / output /
#' session / rv via closure so the handler can re-attach the poller
#' in module scope on Resume.
#'
#' @param kind                String tag ("bayes_tg", …) used to look
#'                            up this handler from the banner.
#' @param tab_name            navbarPage tab to switch to on Resume.
#' @param modal_title,modal_detail  Strings shown by
#'                            `show_analysis_modal()` when the job is
#'                            still running. Detail is prepended with
#'                            "Reconnected to a fit started in a
#'                            previous session." automatically.
#' @param input,output,session,rv  Module scope, captured via closure.
#' @param job_id_field        Name of the rv slot for this module's
#'                            active job_id.
#' @param finalize_success,finalize_failure,finalize_cancellation
#'                            See [attach_bayes_poller()].
#' @keywords internal
register_bayes_resume <- function(kind,
                                  tab_name,
                                  modal_title,
                                  modal_detail,
                                  input, output, session, rv,
                                  job_id_field,
                                  finalize_success,
                                  finalize_failure,
                                  finalize_cancellation,
                                  finalize_lost = NULL,
                                  cancel_input = NULL) {
  if (!is.null(cancel_input)) observe_bayes_cancel(input, rv, cancel_input, job_id_field)
  # R20.D64: registered on this session, so a Resume clicked in another
  # browser cannot reach this one's handler.
  register_resume_handler(kind, session = session, function(job, app_session, app_rv) {
    rv[[job_id_field]] <- job$id
    attach_bayes_poller(input, output, session, rv,
                        job_id_field,
                        finalize_success,
                        finalize_failure,
                        finalize_cancellation,
                        finalize_lost)
    shiny::updateNavbarPage(
      session  = app_session,
      inputId  = "main_nav",
      selected = tab_name)
    # No window for a job that runs while its results are on screen (R20.D6).
    if (job$status %in% c("queued", "running") && !is.null(modal_title)) {
      show_analysis_modal(
        modal_title,
        paste("Reconnected to an analysis started in a previous session.",
              modal_detail),
        cancel_id = if (!is.null(cancel_input)) session$ns(cancel_input))
    }
  })
}

# ---- Long frequentist analyses on the worker (CODE_REVIEW.md R20.D6) --------
#
# The Bayesian fits already ran on the background worker (R20.D62). The power
# simulation (121 s with its defaults, and reachable before any upload), the
# randomisation test (22 s at 999 permutations, 145 s at 9,999) and the
# body-weight model (15 s) ran in the Shiny process, and every session on the
# server froze until they finished. They now use the same queue: one job per
# browser, a progress window with Cancel, and the module's finalizers.

#' Run a call on the worker, bringing its warnings back with the value
#'
#' The registry collects a future's value with its warnings suppressed, and
#' the dashboard shows backend warnings (R20.D34), so they travel in the
#' result as `analysis_warnings`, the field `safe_analysis()` uses.
#' @keywords internal
.run_capturing <- function(fun, args) {
  w <- character()
  value <- withCallingHandlers(do.call(fun, args),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  if (length(w) && (is.list(value) || is.null(value))) {
    if (is.null(value)) value <- list()   # a NULL result keeps its reason
    value$analysis_warnings <- unique(w)
  }
  value
}

#' A registry task that runs `fun(args)` on the worker
#'
#' @param fun A function the worker can call: a package function such as
#'   `mouseExperiment::apriori_power_simulation`, or one whose environment is
#'   the base environment.
#' @param args Named list of arguments.
#' @param packages Packages the worker attaches.
#' @return `list(fun, args, packages)` for [launch_job()].
#' @keywords internal
background_task <- function(fun, args, packages = "mouseExperiment") {
  run <- .run_capturing
  environment(run) <- baseenv()
  list(fun = run, args = list(fun = fun, args = args), packages = packages)
}

#' Launch a long analysis on the worker
#'
#' [launch_tg_bayes_job()], made general. One job per browser (R20.D7): with
#' one already queued or running, it says so and returns NULL.
#'
#' @param kind Registry kind, labelled by `.kind_to_label()`.
#' @param task From [background_task()].
#' @param context Passed to the finalizers.
#' @param modal_title,modal_detail The progress window, or `modal_title = NULL`
#'   for none (the randomisation test runs while its results are on screen).
#' @param cancel_input Unnamespaced id of the window's Cancel button.
#' @inheritParams attach_bayes_poller
#' @return The job id, or NULL when this browser already has a job.
#' @keywords internal
launch_background_job <- function(kind, task, context = NULL,
                                  input, output, session, rv, job_id_field,
                                  modal_title = NULL, modal_detail = NULL,
                                  cancel_input = NULL,
                                  finalize_success, finalize_failure,
                                  finalize_cancellation, finalize_lost = NULL) {
  sid <- rv$session_id
  if (is.null(sid) || identical(sid, SID_NONE)) {
    sid <- make_local_sid()
    rv$session_id <- sid
  }
  if (!is.null(busy <- active_job(sid))) {
    shiny::showNotification(sprintf(paste0(
      "%s from this browser is already %s. Wait for it, or cancel it from the ",
      "banner at the top of the page."),
      .kind_to_label(busy$kind),
      if (identical(busy$status, "queued")) "waiting" else "running"),
      type = "warning", duration = 8)
    return(invisible(NULL))
  }
  if (!is.null(modal_title)) {
    show_analysis_modal(modal_title, modal_detail,
                        cancel_id = if (!is.null(cancel_input)) session$ns(cancel_input))
  }
  job_id <- launch_job(sid, kind, task = task, context = context)
  rv[[job_id_field]] <- job_id
  attach_bayes_poller(input, output, session, rv, job_id_field,
                      finalize_success, finalize_failure,
                      finalize_cancellation, finalize_lost)
  invisible(job_id)
}

#' A failed background job as `safe_analysis()` reports a failure
#' @keywords internal
background_failure <- function(err_msg, module) {
  log_error(sprintf("analysis_failed | module=%s | error=%s", module, err_msg))
  list(error = TRUE, message = translate_error(err_msg %||% "unknown error"))
}
