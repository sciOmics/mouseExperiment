#' Session Management Module
#'
#' Server-side logic for the Session Save / Restore controls in the Data
#' Upload tab. (The Analysis Log viewer was removed in v26.05.020; the
#' \code{logger} package continues to write structured log lines to
#' stdout — visible via \code{journalctl} / \code{docker logs} — so the
#' operational signal is preserved without a user-facing tab.)
#'
#' @name log_session_modules
#' @keywords internal

# ---------------------------------------------------------------------------
# Session Save / Restore
# ---------------------------------------------------------------------------

#' Session Save/Restore UI
#'
#' Download / upload buttons for saving and restoring session state.
#' Intended to be placed inside the Data Upload tab.
#'
#' @param ns Namespace function from the enclosing module.
#' @return A \code{tagList}.
#' @keywords internal
session_save_restore_ui <- function(ns) {
  wellPanel(
    h5(strong(icon("floppy-disk"), " Session Save / Restore")),
    p(style = "font-size:12px; color:#6c757d;",
      "Save your current analysis state (settings, results, plots) ",
      "and reload it later to continue where you left off."),
    fluidRow(
      column(6,
        downloadButton(ns("save_session"), "Save Session",
                       class = "btn-sm btn-success", style = "width:100%;")
      ),
      column(6,
        fileInput(ns("restore_session"), NULL,
                  accept = ".rds",
                  buttonLabel = "Restore Session",
                  placeholder = "No file")
      )
    )
  )
}

#' Session save handler — call from inside the data upload module server
#'
#' @param input   Module input.
#' @param output  Module output.
#' @param session Module session.
#' @param rv      Shared reactiveValues.
#' @keywords internal
setup_session_save_restore <- function(input, output, session, rv) {

  output$save_session <- downloadHandler(
    filename = function() {
      paste0("mouseExperiment_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      # Snapshot serialisable state only (skip reactive internals, ggplot envs)
      snapshot <- list(
        uploaded_data  = isolate(rv$uploaded_data),
        data_valid     = isolate(rv$data_valid),
        detected_cols  = isolate(rv$detected_cols),
        settings       = isolate(rv$settings),
        processed_data = isolate(rv$processed_data),
        time_col       = isolate(rv$time_col),
        vol_col        = isolate(rv$vol_col),
        tx_col         = isolate(rv$tx_col),
        id_col         = isolate(rv$id_col),
        cage_col       = isolate(rv$cage_col),
        dose_col       = isolate(rv$dose_col),
        tg_col_names   = isolate(rv$tg_col_names),
        tg_plot_groups = isolate(rv$tg_plot_groups),
        surv_plot_groups = isolate(rv$surv_plot_groups),
        saved_at       = Sys.time()
      )
      saveRDS(snapshot, file)
      log_info("session_save | file=mouseExperiment_session.rds")
      shinyjs::runjs("window.trackEvent('session_saved', {format: 'rds'});")
    }
  )

  observeEvent(input$restore_session, {
    file_info <- input$restore_session
    req(file_info)
    tryCatch({
      snapshot <- readRDS(file_info$datapath)
      if (!is.list(snapshot) || is.null(snapshot$uploaded_data)) {
        showNotification("Invalid session file.", type = "error")
        shinyjs::runjs("window.trackError('session', 'Invalid session file', 'session_restore');")
        return()
      }
      # Restore keys
      for (key in names(snapshot)) {
        if (key != "saved_at") rv[[key]] <- snapshot[[key]]
      }
      shinyjs::runjs("window.trackEvent('session_restored', {format: 'rds'});")
      showNotification(
        paste("Session restored from", format(snapshot$saved_at, "%Y-%m-%d %H:%M")),
        type = "message", duration = 5
      )
      log_info(sprintf("session_restore | saved_at=%s", snapshot$saved_at))
    }, error = function(e) {
      showNotification(paste("Failed to restore session:", e$message), type = "error")
      shinyjs::runjs(sprintf("window.trackError('session', %s, 'session_restore');",
                             js_string(e$message)))
      log_caught_error("session", "restore", e)
    })
  })
}
