#' Analysis Log and Session Management Modules
#'
#' Server-side logic for the Analysis Log viewer and Session Save/Restore.
#'
#' @name log_session_modules
#' @keywords internal

# ---------------------------------------------------------------------------
# Analysis Log Module
# ---------------------------------------------------------------------------

#' Analysis Log Server
#'
#' Provides a live view of the in-process log buffer (see \code{logging.R}).
#' Users can refresh, clear, and download the log.
#'
#' @param id Module namespace id.
#' @param rv Shared reactiveValues.
#' @keywords internal
analysis_log_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

    # Auto-refresh every 5 seconds, or on manual refresh
    log_trigger <- reactiveTimer(5000)

    log_text <- reactive({
      log_trigger()
      input$refresh_log   # also invalidate on button click
      entries <- get_log_entries()
      if (length(entries) == 0L) return("No log entries yet.")
      formatted <- vapply(rev(entries), function(line) {
        tryCatch({
          p <- jsonlite::fromJSON(line)
          sprintf("[%s] %-5s  %s", p$time, p$level, p$msg)
        }, error = function(e) line)
      }, character(1L))
      paste(formatted, collapse = "\n")
    })

    output$log_output <- renderText(log_text())

    observeEvent(input$clear_log, {
      clear_log_buffer()
    })

    output$download_log <- downloadHandler(
      filename = function() {
        paste0("analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
      },
      content = function(file) {
        writeLines(get_log_entries(), file)
      }
    )
  })
}

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
    }
  )

  observeEvent(input$restore_session, {
    file_info <- input$restore_session
    req(file_info)
    tryCatch({
      snapshot <- readRDS(file_info$datapath)
      if (!is.list(snapshot) || is.null(snapshot$uploaded_data)) {
        showNotification("Invalid session file.", type = "error")
        return()
      }
      # Restore keys
      for (key in names(snapshot)) {
        if (key != "saved_at") rv[[key]] <- snapshot[[key]]
      }
      showNotification(
        paste("Session restored from", format(snapshot$saved_at, "%Y-%m-%d %H:%M")),
        type = "message", duration = 5
      )
      log_info(sprintf("session_restore | saved_at=%s", snapshot$saved_at))
    }, error = function(e) {
      showNotification(paste("Failed to restore session:", e$message), type = "error")
      log_caught_error("session", "restore", e)
    })
  })
}
