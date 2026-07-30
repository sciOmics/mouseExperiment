#' UI chrome: modals, alerts, error translation and JS escaping
#'
#' Split out of the former `helpers.R` under CODE_REVIEW.md D.1.
#'
#' Presentation and error handling. `safe_analysis()` and `translate_error()` are
#' the ones that matter: they decide what a user is told when an analysis fails,
#' and a translation that swallows the cause turns a diagnosable error into a
#' shrug. `js_string()` escapes values interpolated into `runjs()` calls.

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
    message = "Analysis running…",
    detail  = "Results will appear automatically when complete.") {
  shiny::showModal(shiny::modalDialog(
    title     = NULL,
    footer    = NULL,
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
      shiny::tags$p(detail, style = "color:#888; font-size:0.9em; margin-bottom:0;")
    )
  ))
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
  icon_cls <- if (type == "info") paste(icon_name, "fa-spin") else icon_name
  shiny::div(class = paste0("alert alert-", type),
             shiny::icon(icon_cls), msg)
}

#' Safely execute analysis function with error handling
#'
#' @param expr Expression to execute
#' @param module Optional module name for structured logging.
#'
#' @return Result or error message
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
    return("Matrix dimension mismatch — check that all treatment groups have sufficient observations.")
  if (grepl("NaN|NA|missing values", msg, ignore.case = TRUE))
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

safe_analysis <- function(expr, module = "unknown") {
  t0 <- proc.time()[["elapsed"]]
  tryCatch(
    {
      result <- expr
      elapsed <- proc.time()[["elapsed"]] - t0
      log_info(sprintf("analysis_complete | module=%s | elapsed_s=%.2f", module, round(elapsed, 2)))
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
