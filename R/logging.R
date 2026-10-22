#' Structured logging for the mouseExperiment dashboard
#'
#' Provides a lightweight logging layer using the \pkg{logger} package.
#' In a Docker container, logs are emitted to \code{stdout} as JSON lines
#' and captured automatically by the container's log driver.
#'
#' @name logging
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

#' Initialise the dashboard logger
#'
#' Called once from the server function to configure
#' the logger namespace, threshold, layout, and appender.
#'
#' @param level Log threshold (default \code{logger::INFO}).
#'   Use \code{logger::DEBUG} for verbose development output.
#' @keywords internal
setup_logging <- function(level = logger::INFO) {
  # Use a dedicated namespace to avoid interfering with other packages
  logger::log_threshold(level, namespace = "dashboard")
  logger::log_layout(logger::layout_json(), namespace = "dashboard")
  logger::log_appender(logger::appender_stdout, namespace = "dashboard")
}

# ---------------------------------------------------------------------------
# Convenience wrappers — all write to the "dashboard" namespace
# ---------------------------------------------------------------------------

#' @keywords internal
log_info <- function(msg, ...) {
  logger::log_info(logger::skip_formatter(msg), namespace = "dashboard")
}

#' @keywords internal
log_warn <- function(msg, ...) {
  logger::log_warn(logger::skip_formatter(msg), namespace = "dashboard")
}

#' @keywords internal
log_error <- function(msg, ...) {
  logger::log_error(logger::skip_formatter(msg), namespace = "dashboard")
}

#' @keywords internal
log_debug <- function(msg, ...) {
  logger::log_debug(logger::skip_formatter(msg), namespace = "dashboard")
}

# ---------------------------------------------------------------------------
# Structured event helpers
# ---------------------------------------------------------------------------

#' Log a session lifecycle event (start / end)
#'
#' @param session Shiny session object.
#' @param event   One of \code{"start"} or \code{"end"}.
#' @keywords internal
log_session <- function(session, event = "start") {
  ua <- tryCatch(session$request$HTTP_USER_AGENT, error = function(e) "unknown")
  log_info(sprintf("session_%s | session=%s | user_agent=%s",
                   event, session$token, ua %||% "unknown"))
}

#' Log a data-upload event
#'
#' @param session  Shiny session object.
#' @param filename Uploaded filename.
#' @param nrow     Number of rows.
#' @param ncol     Number of columns.
#' @keywords internal
log_upload <- function(session, filename, nrow, ncol) {
  log_info(sprintf("data_upload | session=%s | file=%s | rows=%d | cols=%d",
                   session$token, filename, nrow, ncol))
}



#' Log an error caught by a tryCatch handler
#'
#' @param module  Module name or context label.
#' @param context Brief description of what was attempted.
#' @param error   Condition object or character message.
#' @keywords internal
log_caught_error <- function(module, context, error) {
  err <- if (inherits(error, "condition")) conditionMessage(error) else as.character(error)
  log_error(sprintf("caught_error | module=%s | context=%s | error=%s",
                    module, context, err))
}



#' Log a tab navigation event
#'
#' @param session Shiny session object.
#' @param tab     Tab name.
#' @keywords internal
log_navigation <- function(session, tab) {
  log_debug(sprintf("navigation | session=%s | tab=%s",
                    session$token, tab))
}

# R20.D25: every log line was also kept in a process-wide 500-line buffer for
# the Analysis Log panel, which was removed (I.4). Nothing read it; logging
# goes to stdout only (setup_logging()).

