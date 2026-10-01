#' Resume-banner module (CODE_REVIEW.md F.7).
#'
#' UI module that lists outstanding / completed / failed jobs for the
#' current sid and offers per-job Resume / Cancel actions plus a
#' Dismiss-all. Renders nothing when the sid has no jobs — so the
#' banner is invisible on fresh sessions.
#'
#' Wiring contract:
#'   * `session_id` is a reactive that returns the current sid (set by
#'     `app.R` from `get_sid(session)` once the cookie arrives).
#'   * `on_resume` is a callback `function(job)` provided by the parent
#'     scope. It's responsible for populating `rv$<kind>_results` and
#'     switching the navbar to the right tab (the banner is content-
#'     unaware; it just hands the job back).
#'
#' Click delegation: per-row buttons use a single delegated JS handler
#' that sets one of three module-scoped inputs (`resume_clicked`,
#' `cancel_clicked`, `dismiss_clicked`) each with an event-priority
#' message so repeat clicks fire even when the value is the same. This
#' avoids registering an observer per dynamic button ID.

#' Resume-banner UI
#' @param id moduleServer namespace.
#' @keywords internal
resume_banner_ui <- function(id = "resume_banner") {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("banner")),
    # One delegated click handler per action. Buttons inside the banner
    # carry a `data-jid` attribute that the handler reads. `priority:
    # 'event'` ensures repeat clicks on the same button still fire.
    shiny::tags$script(shiny::HTML(sprintf(paste0(
      "$(document).on('click', '.mE-resume-btn', function() {",
      "  Shiny.setInputValue('%s', ",
      "    {jid: $(this).data('jid'), n: Math.random()}, ",
      "    {priority: 'event'});",
      "});",
      "$(document).on('click', '.mE-cancel-btn', function() {",
      "  Shiny.setInputValue('%s', ",
      "    {jid: $(this).data('jid'), n: Math.random()}, ",
      "    {priority: 'event'});",
      "});",
      "$(document).on('click', '.mE-dismiss-btn', function() {",
      "  Shiny.setInputValue('%s', ",
      "    {n: Math.random()}, ",
      "    {priority: 'event'});",
      "});"
    ),
    ns("resume_clicked"),
    ns("cancel_clicked"),
    ns("dismiss_clicked"))))
  )
}

#' Resume-banner server
#'
#' @param id           moduleServer namespace (must match `resume_banner_ui`).
#' @param session_id   Reactive returning the current sid (character) or
#'                     `SID_NONE` while the cookie hasn't arrived yet.
#' @param on_resume    `function(job)` callback fired when the user
#'                     clicks Resume on `job`. Responsible for moving
#'                     `job$result` into the correct `rv$<kind>_results`
#'                     slot and switching the navbar tab.
#' @keywords internal
resume_banner_server <- function(id = "resume_banner",
                                 session_id,
                                 on_resume) {
  shiny::moduleServer(id, function(input, output, session) {

    # Bump to invalidate the banner UI after any state mutation
    # (resume / cancel / dismiss / periodic poll).
    refresh <- shiny::reactiveVal(0L)

    output$banner <- shiny::renderUI({
      refresh()                                    # depend
      sid <- session_id()
      if (is.null(sid) || identical(sid, SID_NONE)) return(NULL)
      jobs <- list_jobs(sid)
      if (length(jobs) == 0L) return(NULL)
      # The registry is not reactive: re-read it while there is anything to
      # show, so a finished or delivered fit's line updates on its own.
      shiny::invalidateLater(5000L, session)
      .build_resume_banner_ui(jobs, sid)
    })

    # ── Action handlers ──────────────────────────────────────────────
    shiny::observeEvent(input$resume_clicked, {
      sid <- session_id()
      jid <- input$resume_clicked$jid
      if (is.null(sid) || identical(sid, SID_NONE) || is.null(jid)) return()
      job <- get_job(sid, jid)
      if (is.null(job)) return()
      # R20.D63: the job is not retired here. The resume handler attaches a
      # poller that delivers the result and then retires the job; retiring
      # first meant the poller found nothing and the result was lost.
      tryCatch(on_resume(job), error = function(e) {
        shiny::showNotification(
          sprintf("Failed to resume job: %s", conditionMessage(e)),
          type = "error")
      })
      refresh(refresh() + 1L)
    })

    shiny::observeEvent(input$cancel_clicked, {
      sid <- session_id()
      jid <- input$cancel_clicked$jid
      if (is.null(sid) || identical(sid, SID_NONE) || is.null(jid)) return()
      cancel_job(sid, jid)
      refresh(refresh() + 1L)
    })

    shiny::observeEvent(input$dismiss_clicked, {
      sid <- session_id()
      if (is.null(sid) || identical(sid, SID_NONE)) return()
      # R20.D5: only finished jobs. Dismissing a running job stranded the
      # modal of the tab waiting for it and threw its result away; a running
      # or queued fit is stopped with its own Cancel button.
      for (job in list_jobs(sid)) {
        if (job$status %in% c("completed", "failed", "cancelled")) retire_job(sid, job$id)
      }
      refresh(refresh() + 1L)
    })
  })
}

# ---- UI builders (package-private) ----------------------------------------

.build_resume_banner_ui <- function(jobs, sid = NULL) {
  # Sort: running, then queued, then completed (newest first), then failed,
  # then cancelled. Within each status, newest first.
  status_rank <- match(vapply(jobs, function(j) j$status, character(1L)),
                       c("running", "queued", "completed", "failed", "cancelled"))
  queued_t    <- vapply(jobs, function(j) as.numeric(j$queued_at %||% j$started_at),
                        numeric(1L))
  ord <- order(status_rank, -queued_t)
  jobs <- jobs[ord]

  shiny::div(
    class = "alert alert-info mE-resume-banner",
    role  = "alert",
    style = "margin: 12px; padding: 12px 16px;",
    shiny::tags$strong(
      shiny::icon("circle-info"),
      sprintf(" You have %d %s from a previous session:",
              length(jobs),
              if (length(jobs) > 1L) "analyses" else "analysis")
    ),
    shiny::tags$ul(
      class = "list-unstyled",
      style = "margin: 10px 0 6px 24px;",
      lapply(jobs, .render_one_resume_job, sid = sid)
    ),
    if (any(vapply(jobs, function(j) j$status %in% c("completed", "failed", "cancelled"),
                   logical(1L))))
      shiny::tags$button(
        class            = "btn btn-default btn-sm mE-dismiss-btn",
        `data-dismiss`   = "finished",
        type             = "button",
        "Dismiss finished"
      )
  )
}

.render_one_resume_job <- function(job, sid = NULL) {
  since <- if (!is.null(job$started_at) && !is.na(job$started_at)) job$started_at
           else job$queued_at %||% Sys.time()
  elapsed_sec <- as.numeric(difftime(
    if (is.na(job$completed_at)) Sys.time() else job$completed_at,
    since, units = "secs"))
  elapsed_str <- if (elapsed_sec >= 60) {
    sprintf("%d min %d s", as.integer(elapsed_sec) %/% 60L,
                            as.integer(elapsed_sec) %% 60L)
  } else {
    sprintf("%d s", as.integer(elapsed_sec))
  }
  status_lbl <- switch(job$status,
    queued    = {
      ahead <- if (!is.null(sid)) jobs_ahead(sid, job$id) else 0L
      sprintf("waiting for the fitting worker (%d %s ahead)", ahead,
              if (ahead == 1L) "fit" else "fits")
    },
    running   = sprintf("running (%s elapsed)", elapsed_str),
    completed = sprintf("completed %s ago", elapsed_str),
    failed    = "failed",
    cancelled = "cancelled",
    job$status)
  kind_lbl <- .kind_to_label(job$kind)
  shiny::tags$li(
    style = "margin-bottom: 6px;",
    shiny::tags$strong(kind_lbl), " — ", status_lbl,
    shiny::span(style = "margin-left: 12px;",
      if (job$status %in% c("queued", "running", "completed")) {
        shiny::tags$button(
          class            = "btn btn-primary btn-xs mE-resume-btn",
          `data-jid`       = job$id,
          type             = "button",
          style            = "margin-right: 4px;",
          "Resume"
        )
      },
      if (job$status %in% c("queued", "running")) {
        shiny::tags$button(
          class      = "btn btn-default btn-xs mE-cancel-btn",
          `data-jid` = job$id,
          type       = "button",
          "Cancel"
        )
      }
    )
  )
}

# Map registry `kind` tags to user-facing labels. New kinds added in
# future commits register here.
.kind_to_label <- function(kind) {
  switch(kind,
    bayes_tg    = "Bayesian Tumor Growth",
    bayes_surv  = "Bayesian Survival",
    kind)
}

# ---- Resume dispatch (package-private) ------------------------------------
#
# Each module registers an on-resume handler for its `kind` at server
# startup via `register_resume_handler()`. The banner's on_resume
# callback then dispatches by kind. Decouples the banner from any
# specific module — the banner doesn't import the module servers.

#' Register a resume handler for a specific job kind.
#'
#' Called by each Bayesian module's server at startup. The handler is
#' invoked when the user clicks Resume on a job of this kind; it's
#' responsible for moving `job$result` into the right `rv$<kind>_results`
#' slot and switching the navbar tab.
#'
#' CODE_REVIEW.md R20.D64: the handlers lived in one environment for the
#' whole R process, keyed by kind, so the browser session that started last
#' owned every handler. A Resume clicked in an earlier session ran the newer
#' session's handler, which looked for the job under its own sid, found
#' nothing, and did nothing. They are now kept in the session's `userData`,
#' which a module's session shares with the app's.
#'
#' @param kind    String tag (e.g. "bayes_tg") matching the `kind`
#'                passed to [launch_job()].
#' @param handler `function(job, session, rv)` — receives the registry
#'                job entry plus the active session and rv.
#' @param session The session the handler belongs to.
#' @return `invisible(NULL)`.
#' @keywords internal
register_resume_handler <- function(kind, handler,
                                    session = shiny::getDefaultReactiveDomain()) {
  stopifnot(is.character(kind), length(kind) == 1L, nzchar(kind),
            is.function(handler), !is.null(session))
  if (!is.environment(session$userData$mE_resume_handlers)) {
    session$userData$mE_resume_handlers <- new.env(parent = emptyenv())
  }
  session$userData$mE_resume_handlers[[kind]] <- handler
  invisible(NULL)
}

#' Dispatch a resume click to the registered handler for the job's kind.
#'
#' Called by the banner module's `on_resume` callback. When no handler
#' is registered for the kind (commit 2 / placeholder state), shows a
#' notification rather than silently failing.
#'
#' @keywords internal
dispatch_resume <- function(job, session, rv) {
  handlers <- session$userData$mE_resume_handlers
  handler  <- if (is.environment(handlers)) handlers[[job$kind]]
  if (is.null(handler)) {
    # showNotification needs an active Shiny session. Skip silently when
    # called from unit tests where no session exists.
    if (!is.null(shiny::getDefaultReactiveDomain())) {
      shiny::showNotification(
        sprintf("No resume handler is registered for '%s' yet.", job$kind),
        type = "warning", duration = 5)
    }
    return(invisible(NULL))
  }
  handler(job, session, rv)
}

#' Reset a session's registered resume handlers. Test-only.
#' @keywords internal
.reset_resume_handlers <- function(session = shiny::getDefaultReactiveDomain()) {
  if (!is.null(session)) session$userData$mE_resume_handlers <- NULL
  invisible(NULL)
}
