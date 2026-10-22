# UI chrome: modals, alerts, error translation and JS escaping
#
# Split out of the former `helpers.R` under CODE_REVIEW.md D.1.
#
# Presentation and error handling. `safe_analysis()` and `translate_error()` are
# the ones that matter: they decide what a user is told when an analysis fails,
# and a translation that swallows the cause turns a diagnosable error into a
# shrug. `js_string()` escapes values interpolated into `runjs()` calls.

# ---- Bayesian progress modal ------------------------------------------------

#' Show a modal dialog while a Bayesian fit is running
#'
#' Uses shiny::showModal() which flushes to the browser immediately, before
#' the future_promise() computation starts. Dismiss with shiny::removeModal().
#'
#' @param message Title line shown in the dialog.
#' @param detail  Subtitle / detail line.
#' @keywords internal
show_analysis_modal <- function(
    message = "Analysis running\u2026",
    detail  = "Results will appear automatically when complete.",
    cancel_id = NULL) {
  shiny::showModal(shiny::modalDialog(
    title     = NULL,
    # R20.D7: a Bayesian fit can be cancelled from here; the banner only
    # lists fits when the page loads.
    footer    = if (!is.null(cancel_id))
                  shiny::actionButton(cancel_id, "Cancel", class = "btn-default btn-sm"),
    easyClose = FALSE,
    size      = "s",
    shiny::tags$div(
      style = "text-align:center; padding:24px 8px 16px;",
      shiny::tags$i(
        class = "fa fa-spinner fa-spin fa-2x",
        style = "color:var(--primary-color); margin-bottom:16px; display:block;"
      ),
      shiny::tags$h4(message, style = "font-weight:600; margin-bottom:16px;"),
      shiny::tags$div(
        class = "progress",
        style = "height:6px; border-radius:3px; margin-bottom:16px;",
        shiny::tags$div(
          class = "progress-bar progress-bar-striped active",
          role  = "progressbar",
          style = "width:100%;"
        )
      ),
      shiny::tags$p(detail, style = "color:#888; font-size:0.9em; margin-bottom:0;"),
      # Filled in by the Bayesian job poller: queue position or elapsed time.
      shiny::tags$p(id = "mE-job-status",
                    style = "color:#666; font-size:0.85em; margin:8px 0 0;")
    )
  ))
}

#' The analysis modal's status line for a queued or running job
#'
#' CODE_REVIEW.md R20.D7: fits run one at a time on the server's background
#' worker, so a fit can wait behind another user's. Say so, with the queue
#' position, instead of a spinner that looks like the fit itself.
#' @param job A registry job.
#' @param ahead Fits that will run before it ([jobs_ahead()]).
#' @keywords internal
job_status_text <- function(job, ahead = 0L) {
  if (identical(job$status, "queued")) {
    return(sprintf(paste0("Waiting for the fitting worker: %d %s ahead. ",
                          "Bayesian fits run one at a time on this server."),
                   ahead, if (ahead == 1L) "fit" else "fits"))
  }
  secs <- as.integer(difftime(Sys.time(), job$started_at, units = "secs"))
  if (is.na(secs)) secs <- 0L
  sprintf("Fitting: %d:%02d elapsed. Bayesian fits run one at a time on this server.",
          secs %/% 60L, secs %% 60L)
}

# ---- Alert helpers (consolidated validate/need pattern) ----------------------

#' Create a Bootstrap alert div for module status messages
#'
#' @param msg  Character message to display.
#' @param type One of "danger", "warning", "info", "success".
#' @return A \code{shiny::div} tag.
#' @keywords internal
status_alert <- function(msg, type = "danger") {
  icon_name <- switch(type,
    danger  = "exclamation-circle",
    warning = "exclamation-triangle",
    info    = "spinner",
    success = "check",
    "info-circle"
  )
  # R20.D21: "spinner fa-spin" as the icon *name* logged "does not correspond
  # to a known icon" on every info alert; the spin is a class.
  shiny::div(class = paste0("alert alert-", type),
             shiny::icon(icon_name, class = if (type == "info") "fa-spin"), msg)
}

#' Multiple-comparison corrections offered for a tumour-growth model
#'
#' CODE_REVIEW.md R20.D39 -- the AUC path compares arms with independent Welch
#' t-tests, which have no common error term or joint covariance, so Dunnett and
#' Tukey are undefined there and the backend rejects them. The dropdown's
#' default was Dunnett, so choosing AUC and clicking Run failed every time.
#'
#' @param model_type `"lme4"`, `"auc"` or `"bayes"`.
#' @return Named character vector for `selectInput(choices =)`.
#' @keywords internal
tg_p_adjust_choices <- function(model_type = "lme4") {
  joint <- c("Dunnett (exact, vs reference)" = "dunnett",
             "Tukey HSD (all pairs)"         = "tukey")
  family <- c("Bonferroni"             = "bonferroni",
              "Holm"                   = "holm",
              "Benjamini-Hochberg FDR" = "fdr",
              "None (unadjusted)"      = "none")
  if (identical(model_type, "auc")) family else c(joint, family)
}

#' Translate an R error message into one a user can act on
#'
#' @param msg The error message.
#' @return A character string: the plain-language message, or `msg` itself
#'   when no translation applies.
#' @keywords internal
translate_error <- function(msg) {
  if (grepl("not found in|not found$|object .* not found", msg, perl = TRUE))
    return("A required column was not found. Check your column mappings in Data Upload.")
  if (grepl("singular|singular fit", msg, ignore.case = TRUE))
    return(paste0(
      "Model fitting failed (singular fit). ",
      "Try increasing sample size or simplifying the random effects structure."
    ))
  if (grepl("subscript out of bounds", msg, fixed = TRUE))
    return(paste0(
      "An internal indexing error occurred. ",
      "Check that your data is complete and all required columns are mapped."
    ))
  if (grepl("non-conformable", msg, fixed = TRUE))
    return("Matrix dimension mismatch \u2014 check that all treatment groups have sufficient observations.")
  # R20.D61: whole words, case-sensitive; ignore.case matched any "na", as in
  # "name" or "analysis".
  if (grepl("\\bNaN\\b|\\bNA\\b|missing values", msg, perl = TRUE))
    return(paste0(
      "The model encountered missing or invalid values. ",
      "Check your data for NA values or zero-variance groups."
    ))
  msg
}

toggle_downloads <- function(session, ids, enabled) {
  fn <- if (isTRUE(enabled)) shinyjs::enable else shinyjs::disable
  invisible(lapply(ids, function(id) fn(session$ns(id))))
}

#' Safely execute analysis function with error handling
#'
#' @param expr Expression to execute
#' @param module Optional module name for structured logging.
#'
#' @return Result or error message
#' @keywords internal
safe_analysis <- function(expr, module = "unknown") {
  t0 <- proc.time()[["elapsed"]]
  # CODE_REVIEW.md R20.D34 -- there was no warning handler, so every backend
  # warning was dropped: a 1000x unit error ran to "Analysis completed
  # successfully!" while the package warned about it. Collect them, log them,
  # and attach them to the result for completion_status() to show.
  warns <- character(0)
  tryCatch(
    {
      result <- withCallingHandlers(expr, warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      elapsed <- proc.time()[["elapsed"]] - t0
      log_info(sprintf("analysis_complete | module=%s | elapsed_s=%.2f", module, round(elapsed, 2)))
      warns <- unique(warns)
      for (w in warns) log_warn(sprintf("analysis_warning | module=%s | %s", module, w))
      if (length(warns) && is.list(result)) result$analysis_warnings <- warns
      result
    },
    error = function(e) {
      elapsed <- proc.time()[["elapsed"]] - t0
      raw_msg <- conditionMessage(e)
      log_error(sprintf("analysis_failed | module=%s | elapsed_s=%.2f | error=%s",
                        module, round(elapsed, 2), raw_msg))
      list(error = TRUE, message = translate_error(raw_msg))
    }
  )
}

js_string <- function(x) {
  jsonlite::toJSON(as.character(x %||% ""), auto_unbox = TRUE)
}

#' Completion message that carries the backend's warnings
#'
#' CODE_REVIEW.md R20.D34. A run that produced warnings is reported as
#' "completed with warnings", listing them (the first five), rather than as a
#' plain success.
#'
#' @param results A list of result lists (one element for most tabs; Toxicity
#'   runs several analyses).
#' @param success_msg Message for a clean run.
#' @return A status alert.
#' @keywords internal
completion_status <- function(results, success_msg = "Analysis completed successfully!") {
  warns <- unique(unlist(lapply(results, function(r) {
    if (is.list(r)) r$analysis_warnings else NULL
  })))
  if (!length(warns)) return(status_alert(success_msg, "success"))
  shown <- utils::head(warns, 5L)
  shiny::div(class = "alert alert-warning",
    shiny::icon("exclamation-triangle"),
    " ", shiny::strong(sprintf("%s with %d warning(s):",
                               sub("\\s*successfully!?$", "", success_msg),
                               length(warns))),
    shiny::tags$ul(lapply(shown, shiny::tags$li)),
    if (length(warns) > 5L) shiny::tags$small(length(warns) - 5L, " more in the server log."))
}

#' The study's control group, in one analysis's labels
#'
#' The study-design panel on Data Upload records one control group, a value of
#' the treatment column, and every tab uses it (it replaced the per-tab
#' reference selectors). Tumour Growth and Survival label arms "Treatment Dose"
#' when a dose column is mapped, so the control maps to its matching label.
#'
#' @param rv Shared reactiveValues.
#' @param labels Group labels the analysis uses.
#' @return One label, or NULL when there are no labels.
#' @keywords internal
study_reference <- function(rv, labels) {
  labels <- unique(as.character(stats::na.omit(labels)))
  if (!length(labels)) return(NULL)
  ctrl <- rv$settings$control_group %||% ""
  if (nzchar(ctrl)) {
    if (ctrl %in% labels) return(ctrl)
    with_dose <- labels[startsWith(labels, paste0(ctrl, " "))]
    if (length(with_dose)) return(sort(with_dose)[1L])
  }
  pick_reference_group(labels)
}

#' Read-only note naming a tab's reference group
#'
#' The reference is no longer chosen per tab (R20-D.H). The note says which
#' group is used and, when the study's control group is not in the analysis
#' data (for example, excluded by the subset filter), that another group
#' stands in for it.
#'
#' @param rv Shared reactiveValues.
#' @param ref The reference label in use, from [study_reference()].
#' @keywords internal
reference_note <- function(rv, ref) {
  where <- "the control group under Study Design on Data Upload"
  if (is.null(ref) || !nzchar(ref)) {
    return(shiny::tags$p(class = "text-muted", style = "font-size:12px;",
                         "Reference group: set ", where, "."))
  }
  ctrl <- rv$settings$control_group %||% ""
  is_ctrl <- nzchar(ctrl) && (identical(ref, ctrl) || startsWith(ref, paste0(ctrl, " ")))
  shiny::tags$p(style = "font-size:12px;", "Reference group: ", shiny::strong(ref),
    if (is_ctrl) paste0(" (", where, ").")
    else if (nzchar(ctrl)) shiny::tags$span(class = "text-warning", paste0(
      ". The control group '", ctrl, "' is not in this analysis's data ",
      "(is it excluded by the subset filter?), so this group stands in for it."))
    else paste0(". Set ", where, "."))
}

#' The declared volume unit for labels
#'
#' @param rv Shared reactiveValues.
#' @return "mm³", "cm³", or "units not declared".
#' @keywords internal
volume_unit_label <- function(rv) {
  switch(rv$settings$volume_units %||% "",
         mm3 = "mm\u00b3", cm3 = "cm\u00b3", "units not declared")
}

#' The declared volume unit, as the backend's `volume_units` argument
#'
#' Backend 0.30.0 (R20.32) floors zero volumes at 2 mm³ or 0.002 cm³ before
#' the log, so tumour growth needs the unit. Undeclared units are treated as
#' mm³, the backend's default; it warns when the volumes look like cm³.
#'
#' @param rv Shared reactiveValues.
#' @return "mm3" or "cm3".
#' @keywords internal
backend_volume_units <- function(rv) {
  if (identical(rv$settings$volume_units, "cm3")) "cm3" else "mm3"
}

#' The evaluable-day rule, as the backend states it
#'
#' Maintainer decision R20-K (confirmed 2026-09-30): synergy, the therapeutic
#' window and dose-response report only days on which every arm compared is
#' evaluable. The wording matches `mouseExperiment::evaluable_days()$rule`; a
#' test keeps them in step.
#' @keywords internal
evaluable_rule_text <- function() {
  paste0("An arm is evaluable on a day when at least 50% of its enrolled ",
         "animals, and at least 3 animals, are still on study (measured on ",
         "that day or later).")
}

#' Which day an endpoint result used, and which days the rule excluded
#'
#' Every result that uses an evaluable day names the day and lists the days
#' and arms the rule excluded (R20-K: the rule must be relayed to users).
#'
#' @param ev The `evaluability` element of a backend result.
#' @param day The day used; `NULL` for an over-time result (every evaluable
#'   day was analysed).
#' @param arms_label The arms compared, as a plural subject for "are
#'   evaluable", e.g. "all four arms".
#' @return A note, or `NULL` when `ev` is missing.
#' @keywords internal
evaluability_note <- function(ev, day = NULL, arms_label = "all arms") {
  if (is.null(ev) || is.null(ev$table)) return(NULL)
  used <- if (!is.null(day)) {
    last <- isTRUE(all.equal(as.numeric(day), as.numeric(ev$last_day)))
    paste0("Evaluated on day ", format(day),
           if (last) paste0(", the last day on which ", arms_label, " are evaluable.")
           else paste0(", an evaluable day (the last is day ", format(ev$last_day), ")."))
  } else if (length(ev$days)) {
    paste0("Analysed on every evaluable day: days ", format(min(ev$days)), "\u2013",
           format(max(ev$days)), ".")
  } else "No day is evaluable."
  excl <- ev$excluded
  shiny::div(class = "alert alert-info", style = "font-size:0.9em; padding:8px 12px;",
    shiny::icon("calendar-check"), " ", used, " ",
    shiny::tags$small(ev$rule %||% evaluable_rule_text()),
    if (!is.null(excl) && nrow(excl)) shiny::tags$div(style = "margin-top:4px;",
      shiny::strong("Not evaluated: "),
      paste(sprintf("day %s (%s)", format(excl$Day), excl$Arms), collapse = "; "),
      "."))
}
