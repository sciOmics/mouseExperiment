#' Session Management Module
#'
#' Server-side logic for the Session Save / Restore controls in the Data
#' Upload tab. (The Analysis Log viewer was removed in v26.05.020; the
#' \code{logger} package continues to write structured log lines to
#' stdout — visible via \code{journalctl} / \code{docker logs} — so the
#' operational signal is preserved without a user-facing tab.)
#'
#' @name mod_log_session
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
    ),
    # CODE_REVIEW.md Round 2 F.2 — settings-only preset save/load.
    # Smaller JSON file portable across datasets; doesn't carry uploaded
    # data or results.
    hr(style = "margin: 12px 0;"),
    h6(strong(icon("sliders-h"), " Settings preset (column mapping + methods)")),
    p(style = "font-size:11px; color:#6c757d; margin: 4px 0 8px 0;",
      "Save just your column mappings and time / volume / formula settings ",
      "to a small JSON file. Apply on a future dataset with the same column ",
      "names to skip the mapping step."),
    fluidRow(
      column(6,
        downloadButton(ns("save_preset"), "Save Preset",
                       class = "btn-sm btn-info", style = "width:100%;")
      ),
      column(6,
        fileInput(ns("load_preset"), NULL,
                  accept = ".json",
                  buttonLabel = "Load Preset",
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

  # ── Settings preset save / load (CODE_REVIEW.md Round 2 F.2) ─────────────
  # The preset contains only the column mappings + time/volume methods —
  # small enough to be portable across datasets with the same headers.

  output$save_preset <- downloadHandler(
    filename = function() {
      paste0("mouseExperiment_preset_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".json")
    },
    content = function(file) {
      preset <- list(
        format_version = "1",
        saved_at       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        settings       = isolate(rv$settings) %||% list()
      )
      jsonlite::write_json(preset, file, pretty = TRUE, auto_unbox = TRUE)
      log_info("preset_save | file=mouseExperiment_preset.json")
      shinyjs::runjs("window.trackEvent('preset_saved', {format: 'json'});")
    }
  )

  observeEvent(input$load_preset, {
    file_info <- input$load_preset
    req(file_info)
    tryCatch({
      preset <- jsonlite::read_json(file_info$datapath, simplifyVector = TRUE)
      if (!is.list(preset) || is.null(preset$settings)) {
        showNotification("Invalid preset file (missing settings).", type = "error")
        return()
      }

      # Apply each column-mapping / methods setting via update*Input. The
      # observe() in data_upload that mirrors inputs into rv$settings will
      # pick up the changes automatically.
      st <- preset$settings
      updateSelectInput(session, "col_treatment",     selected = st$col_treatment     %||% "")
      updateSelectInput(session, "col_id",            selected = st$col_id            %||% "")
      updateSelectInput(session, "col_day",           selected = st$col_day           %||% "")
      updateSelectInput(session, "col_date_for_day",  selected = st$col_date_for_day  %||% "")
      updateSelectInput(session, "col_volume",        selected = st$col_volume        %||% "")
      updateSelectInput(session, "col_length",        selected = st$col_length        %||% "")
      updateSelectInput(session, "col_width",         selected = st$col_width         %||% "")
      updateSelectInput(session, "col_height",        selected = st$col_height        %||% "")
      updateSelectInput(session, "col_cage",          selected = st$col_cage          %||% "")
      updateSelectInput(session, "col_dose",          selected = st$col_dose          %||% "")
      updateSelectInput(session, "col_survival_event",selected = st$col_survival_event%||% "")
      updateSelectInput(session, "col_weight",        selected = st$col_weight        %||% "")
      updateSelectInput(session, "col_sex",           selected = st$col_sex           %||% "")
      updateSelectInput(session, "col_necrotic",      selected = st$col_necrotic      %||% "")
      if (!is.null(st$time_method))    updateRadioButtons(session, "time_method",    selected = st$time_method)
      if (!is.null(st$volume_method))  updateRadioButtons(session, "volume_method",  selected = st$volume_method)
      if (!is.null(st$volume_formula)) updateSelectInput(session,  "volume_formula", selected = st$volume_formula)
      if (!is.null(st$reference_date)) updateTextInput(session,    "reference_date", value    = st$reference_date)

      saved_at <- preset$saved_at %||% "(unknown)"
      showNotification(
        paste("Preset loaded — saved at", saved_at),
        type = "message", duration = 5
      )
      shinyjs::runjs("window.trackEvent('preset_loaded', {format: 'json'});")
      log_info(sprintf("preset_load | saved_at=%s", saved_at))
    }, error = function(e) {
      showNotification(paste("Failed to load preset:", e$message), type = "error")
      log_caught_error("preset", "load", e)
    })
  })
}


# ---------------------------------------------------------------------------
# Subset Filter UI (CODE_REVIEW.md Round 2 F.5)
# ---------------------------------------------------------------------------

#' Subset Filter UI — exclude mice / cages / treatments from every analysis
#'
#' Renders a collapsed wellPanel below the data preview on the Data Upload
#' tab. Lets users exclude specific IDs, treatments, or cages without
#' re-uploading the CSV. The exclusions are written to
#' \code{rv$exclude_ids / exclude_treatments / exclude_cages} and applied
#' transparently by \code{apply_subset_filter()} at the top of every
#' analysis module.
#'
#' @param ns Namespace function from the enclosing module.
#' @return A \code{tagList}.
#' @keywords internal
subset_filter_ui <- function(ns) {
  wellPanel(
    h5(strong(icon("filter"), " Subset Filter"),
       tags$small(style = "color:#6c757d; margin-left:8px;",
                  "Exclude mice, treatments, or cages from all downstream analyses.")),
    fluidRow(
      column(4,
        selectizeInput(ns("exclude_ids"), "Exclude IDs",
                       choices = NULL, multiple = TRUE,
                       options = list(placeholder = "(none excluded)"))
      ),
      column(4,
        selectizeInput(ns("exclude_treatments"), "Exclude Treatments",
                       choices = NULL, multiple = TRUE,
                       options = list(placeholder = "(none excluded)"))
      ),
      column(4,
        selectizeInput(ns("exclude_cages"), "Exclude Cages",
                       choices = NULL, multiple = TRUE,
                       options = list(placeholder = "(none excluded)"))
      )
    ),
    div(style = "margin-top: 8px;",
      actionButton(ns("apply_filter"), "Apply Filtering",
                   icon = icon("filter"), class = "btn-primary"),
      actionButton(ns("clear_filter"), "Clear",
                   icon = icon("times"), class = "btn-default",
                   style = "margin-left: 6px;"),
      uiOutput(ns("subset_filter_status"), inline = TRUE)
    ),
    p(style = "font-size:11px; color:#6c757d; margin:8px 0 0 0;",
      "Edit the lists above, then click ", strong("Apply Filtering"),
      ". The Dataset Summary and Data Preview below update to reflect ",
      "what's currently applied. Filter state survives session save / restore.")
  )
}

#' Server-side handler for the Subset Filter — populates choices, syncs to rv
#'
#' Called from inside the data_upload module server. Populates the three
#' selectizeInput choices once data is uploaded and column mappings are
#' resolved; mirrors the user's selections into \code{rv$exclude_*} so
#' downstream modules can read them via \code{apply_subset_filter()}.
#' @keywords internal
setup_subset_filter <- function(input, output, session, rv) {

  # Populate choices whenever uploaded_data + column mappings change.
  observe({
    df <- rv$uploaded_data
    req(df)
    id_col <- rv$settings$col_id
    tx_col <- rv$settings$col_treatment
    cg_col <- rv$settings$col_cage

    ids <- if (!is.null(id_col) && nzchar(id_col) && id_col %in% colnames(df))
      sort(unique(as.character(df[[id_col]]))) else character(0)
    txs <- if (!is.null(tx_col) && nzchar(tx_col) && tx_col %in% colnames(df))
      sort(unique(as.character(df[[tx_col]]))) else character(0)
    cgs <- if (!is.null(cg_col) && nzchar(cg_col) && cg_col %in% colnames(df))
      sort(unique(as.character(df[[cg_col]]))) else character(0)

    updateSelectizeInput(session, "exclude_ids",
                        choices = ids,
                        selected = intersect(isolate(rv$exclude_ids), ids),
                        server = TRUE)
    updateSelectizeInput(session, "exclude_treatments",
                        choices = txs,
                        selected = intersect(isolate(rv$exclude_treatments), txs),
                        server = TRUE)
    updateSelectizeInput(session, "exclude_cages",
                        choices = cgs,
                        selected = intersect(isolate(rv$exclude_cages), cgs),
                        server = TRUE)
  })

  # Explicit Apply / Clear behaviour. The selectizeInputs are staging
  # the user's pending choice; only on Apply do we commit those choices
  # to rv$exclude_* and invalidate the preprocessing cache so analyses
  # and the Data Preview pick up the new exclusion state.
  observeEvent(input$apply_filter, {
    rv$exclude_ids        <- input$exclude_ids        %||% character(0)
    rv$exclude_treatments <- input$exclude_treatments %||% character(0)
    rv$exclude_cages      <- input$exclude_cages      %||% character(0)
    rv$processed_data     <- NULL
    rv$filter_applied_at  <- Sys.time()
    n <- length(rv$exclude_ids) + length(rv$exclude_treatments) +
         length(rv$exclude_cages)
    shiny::showNotification(
      sprintf("Filter applied: %d exclusion(s) active.", n),
      type = "message", duration = 3)
  })

  observeEvent(input$clear_filter, {
    updateSelectizeInput(session, "exclude_ids",
                        selected = character(0))
    updateSelectizeInput(session, "exclude_treatments",
                        selected = character(0))
    updateSelectizeInput(session, "exclude_cages",
                        selected = character(0))
    rv$exclude_ids        <- character(0)
    rv$exclude_treatments <- character(0)
    rv$exclude_cages      <- character(0)
    rv$processed_data     <- NULL
    rv$filter_applied_at  <- Sys.time()
    shiny::showNotification("Filter cleared.", type = "message",
                            duration = 2)
  })

  # Status line — shows what's currently applied. Recomputes whenever
  # the user clicks Apply or Clear (via the filter_applied_at tick) or
  # when rv$uploaded_data changes (e.g., new file or session restore).
  output$subset_filter_status <- renderUI({
    rv$filter_applied_at  # tick dependency
    n_ids <- length(rv$exclude_ids        %||% character(0))
    n_tx  <- length(rv$exclude_treatments %||% character(0))
    n_cg  <- length(rv$exclude_cages      %||% character(0))
    if (n_ids + n_tx + n_cg == 0L) {
      return(tags$span(
        style = "color:#6c757d; margin-left: 14px; font-size: 12px;",
        "No filter applied — full dataset is shown below."))
    }
    if (!is.null(rv$uploaded_data) && nrow(rv$uploaded_data) > 0L) {
      n_full <- nrow(rv$uploaded_data)
      n_filt <- nrow(apply_subset_filter(rv$uploaded_data, rv))
      tags$span(
        style = "color: #155724; margin-left: 14px; font-size: 12px;",
        icon("check-circle"),
        sprintf(" Active filter — %d of %d rows shown (%d excluded).",
                n_filt, n_full, n_full - n_filt))
    } else {
      tags$span(
        style = "color: #155724; margin-left: 14px; font-size: 12px;",
        icon("check-circle"),
        sprintf(" Active filter — %d ID(s), %d treatment(s), %d cage(s) excluded.",
                n_ids, n_tx, n_cg))
    }
  })
}
