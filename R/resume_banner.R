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
      .build_resume_banner_ui(jobs)
    })

    # ── Action handlers ──────────────────────────────────────────────
    shiny::observeEvent(input$resume_clicked, {
      sid <- session_id()
      jid <- input$resume_clicked$jid
      if (is.null(sid) || identical(sid, SID_NONE) || is.null(jid)) return()
      job <- get_job(sid, jid)
      if (is.null(job)) return()
      tryCatch(on_resume(job), error = function(e) {
        shiny::showNotification(
          sprintf("Failed to resume job: %s", conditionMessage(e)),
          type = "error")
      })
      retire_job(sid, jid)
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
      for (jid in names(list_jobs(sid))) retire_job(sid, jid)
      refresh(refresh() + 1L)
    })
  })
}

# ---- UI builders (package-private) ----------------------------------------

.build_resume_banner_ui <- function(jobs) {
  # Sort: running first, then completed (newest first), then failed,
  # then cancelled. Within each status, newest first.
  status_rank <- match(vapply(jobs, function(j) j$status, character(1L)),
                       c("running", "completed", "failed", "cancelled"))
  started_t   <- vapply(jobs, function(j) as.numeric(j$started_at),
                        numeric(1L))
  ord <- order(status_rank, -started_t)
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
      lapply(jobs, .render_one_resume_job)
    ),
    shiny::tags$button(
      class            = "btn btn-default btn-sm mE-dismiss-btn",
      `data-dismiss`   = "all",
      type             = "button",
      "Dismiss all"
    )
  )
}

.render_one_resume_job <- function(job) {
  elapsed_sec <- as.numeric(difftime(
    if (is.na(job$completed_at)) Sys.time() else job$completed_at,
    job$started_at, units = "secs"))
  elapsed_str <- if (elapsed_sec >= 60) {
    sprintf("%d min %d s", as.integer(elapsed_sec) %/% 60L,
                            as.integer(elapsed_sec) %% 60L)
  } else {
    sprintf("%d s", as.integer(elapsed_sec))
  }
  status_lbl <- switch(job$status,
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
      if (job$status %in% c("running", "completed")) {
        shiny::tags$button(
          class            = "btn btn-primary btn-xs mE-resume-btn",
          `data-jid`       = job$id,
          type             = "button",
          style            = "margin-right: 4px;",
          "Resume"
        )
      },
      if (job$status == "running") {
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
    bayes_bw    = "Bayesian Body Weight",
    bayes_surv  = "Bayesian Survival",
    bayes_syn   = "Bayesian Synergy",
    bayes_twm   = "Bayesian Therapeutic Window",
    bayes_dr    = "Bayesian Dose-Response",
    kind)
}

# ---- Resume dispatch (package-private) ------------------------------------
#
# Each module registers an on-resume handler for its `kind` at server
# startup via `register_resume_handler()`. The banner's on_resume
# callback then dispatches by kind. Decouples the banner from any
# specific module — the banner doesn't import the module servers.

.resume_handlers <- new.env(parent = emptyenv())

#' Register a resume handler for a specific job kind.
#'
#' Called by each Bayesian module's server at startup. The handler is
#' invoked when the user clicks Resume on a job of this kind; it's
#' responsible for moving `job$result` into the right `rv$<kind>_results`
#' slot and switching the navbar tab.
#'
#' @param kind    String tag (e.g. "bayes_tg") matching the `kind`
#'                passed to [launch_job()].
#' @param handler `function(job, session, rv)` — receives the registry
#'                job entry plus the active session and rv.
#' @return `invisible(NULL)`.
#' @keywords internal
register_resume_handler <- function(kind, handler) {
  stopifnot(is.character(kind), length(kind) == 1L, nzchar(kind),
            is.function(handler))
  .resume_handlers[[kind]] <- handler
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
  handler <- .resume_handlers[[job$kind]]
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

#' Reset registered resume handlers. Test-only.
#' @keywords internal
.reset_resume_handlers <- function() {
  rm(list = ls(.resume_handlers, all.names = TRUE), envir = .resume_handlers)
  invisible(NULL)
}
