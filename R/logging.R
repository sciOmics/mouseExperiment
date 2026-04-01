#' Structured Logging for mouseExperimentDashboard
#'
#' Provides a lightweight logging layer using the \pkg{logger} package.
#' In a Docker container, logs are emitted to \code{stdout} as JSON lines
#' and captured automatically by the container's log driver.
#'
#' @name logging
#' @keywords internal

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

#' Initialise the dashboard logger
#'
#' Called once from the server function (or \code{.onLoad()}) to configure
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

#' Log the start of an analysis run
#'
#' @param session Shiny session object.
#' @param module  Module name (e.g. \code{"tumor_growth"}).
#' @param params  Named list of key analysis parameters.
#' @keywords internal
log_analysis_start <- function(session, module, params = list()) {
  p <- paste(names(params), "=", params, collapse = "; ")
  log_info(sprintf("analysis_start | session=%s | module=%s | params=%s",
                   session$token, module, p))
}

#' Log the completion of an analysis run
#'
#' @param session     Shiny session object.
#' @param module      Module name.
#' @param success     Logical — did the analysis succeed?
#' @param elapsed_sec Wall-clock seconds.
#' @param message     Optional error / status message.
#' @keywords internal
log_analysis_end <- function(session, module, success, elapsed_sec, message = NULL) {
  secs <- round(elapsed_sec, 2)
  if (isTRUE(success)) {
    log_info(sprintf("analysis_complete | session=%s | module=%s | elapsed_s=%.2f",
                     session$token, module, secs))
  } else {
    log_error(sprintf("analysis_failed | session=%s | module=%s | elapsed_s=%.2f | error=%s",
                       session$token, module, secs, message %||% "unknown"))
  }
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

#' Log a plot render event
#'
#' @param session   Shiny session object.
#' @param module    Module name.
#' @param plot_type Which plot was rendered.
#' @keywords internal
log_plot_render <- function(session, module, plot_type) {
  log_debug(sprintf("plot_render | session=%s | module=%s | plot=%s",
                    session$token, module, plot_type))
}

#' Log a download event
#'
#' @param session  Shiny session object.
#' @param module   Module name.
#' @param what     Description of downloaded artifact.
#' @keywords internal
log_download <- function(session, module, what) {
  log_info(sprintf("download | session=%s | module=%s | artifact=%s",
                   session$token, module, what))
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

# ---------------------------------------------------------------------------
# In-process log buffer for the Analysis Log UI panel
# ---------------------------------------------------------------------------

# Private environment holding the circular buffer.
# Shared across all Shiny sessions within a single R process
# (suitable for single-user / Docker deployments).
.log_buf <- new.env(parent = emptyenv())
.log_buf$entries  <- character(0L)
.log_buf$max_size <- 500L

#' Custom logger appender that also writes to the in-process buffer
#'
#' Wraps \code{logger::appender_stdout} and additionally stores every
#' formatted log line in \code{.log_buf$entries}.
#' @param lines Character vector of log lines.
#' @keywords internal
appender_buffer <- function(lines) {
  # Write to stdout (Docker log driver)
  logger::appender_stdout(lines)
  # Also store in ring buffer
  .log_buf$entries <- c(.log_buf$entries, lines)
  if (length(.log_buf$entries) > .log_buf$max_size) {
    .log_buf$entries <- tail(.log_buf$entries, .log_buf$max_size)
  }
}

#' Set up logging with the buffer appender
#'
#' Like \code{setup_logging()} but uses the tee-to-buffer appender so log
#' messages are available for the Analysis Log UI panel.
#'
#' @param level Log threshold (default \code{logger::INFO}).
#' @keywords internal
setup_logging_with_buffer <- function(level = logger::INFO) {
  logger::log_threshold(level, namespace = "dashboard")
  logger::log_layout(logger::layout_json(), namespace = "dashboard")
  logger::log_appender(appender_buffer, namespace = "dashboard")
}

#' Retrieve buffered log entries
#'
#' Returns the current contents of the in-process log buffer as a character
#' vector (newest last).  Used by the Analysis Log panel.
#'
#' @return Character vector of log lines (may be empty).
#' @keywords internal
get_log_entries <- function() {
  .log_buf$entries
}

#' Clear the in-process log buffer
#' @keywords internal
clear_log_buffer <- function() {
  .log_buf$entries <- character(0L)
}
