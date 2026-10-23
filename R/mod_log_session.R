#' Settings Preset and Subset Filter Module
#'
#' Server-side logic for the settings preset and subset filter in the Data
#' Upload tab. The Session Save / Restore controls were removed in
#' v26.09.001 (CODE_REVIEW.md R20.D10, R20-D.H): the saved file carried only
#' the uploaded data, column mappings and plot-group settings -- no results,
#' plots or filters -- and restoring meant deserialising an uploaded R object.
#' Re-uploading the CSV and loading a settings preset gives the same result.
#'
#' @name mod_log_session
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# Settings preset (CODE_REVIEW.md Round 2 F.2)
# ---------------------------------------------------------------------------

#' Settings preset UI
#'
#' Save / load buttons for the column-mapping and method settings preset.
#' Intended to be placed inside the Data Upload tab.
#'
#' @param ns Namespace function from the enclosing module.
#' @return A \code{tagList}.
#' @keywords internal
settings_preset_ui <- function(ns) {
  wellPanel(
    h6(strong(icon("sliders-h"), " Settings preset (column mapping + methods)")),
    p(style = "font-size:11px; color:#6c757d; margin: 4px 0 8px 0;",
      "Save your column mappings, time / volume / formula settings, study ",
      "design and subset filter to a small JSON file. Load it after ",
      "uploading a dataset with the same column names to skip these steps."),
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

#' Settings preset handlers -- call from inside the data upload module server
#'
#' @param input   Module input.
#' @param output  Module output.
#' @param session Module session.
#' @param rv      Shared reactiveValues.
#' @param study_defaults The Data Upload module's reactiveVal of initial Study
#'   Design selections; a loaded preset writes to it so the control group
#'   survives the treatment-column update that follows.
#' @keywords internal
setup_settings_preset <- function(input, output, session, rv,
                                  study_defaults = NULL) {

  # ── Settings preset save / load (CODE_REVIEW.md Round 2 F.2) ─────────────
  # The preset contains only the column mappings + time/volume methods —
  # small enough to be portable across datasets with the same headers.

  output$save_preset <- downloadHandler(
    filename = function() {
      paste0("mouseExperiment_preset_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".json")
    },
    content = function(file) {
      # Version 2 (R20-D.H) adds the study design, which is part of
      # rv$settings, and the subset filter.
      preset <- list(
        format_version = "2",
        saved_at       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        settings       = isolate(rv$settings) %||% list(),
        filters        = list(
          exclude_ids        = I(isolate(rv$exclude_ids)        %||% character(0)),
          exclude_treatments = I(isolate(rv$exclude_treatments) %||% character(0)),
          exclude_cages      = I(isolate(rv$exclude_cages)      %||% character(0)))
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
      # The column selectors and the Study Design panel exist only once a
      # dataset is loaded; updates sent before that were silently lost.
      if (!isTRUE(rv$data_valid)) {
        showNotification(paste("Load a dataset first: a preset applies to",
                               "the loaded data's columns."), type = "warning")
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
      updateSelectInput(session, "col_removal_reason", selected = st$col_removal_reason %||% "")
      if (!is.null(st$time_method))    updateRadioButtons(session, "time_method",    selected = st$time_method)
      if (!is.null(st$volume_method))  updateRadioButtons(session, "volume_method",  selected = st$volume_method)
      if (!is.null(st$volume_formula)) updateSelectInput(session,  "volume_formula", selected = st$volume_formula)
      if (!is.null(st$reference_date)) updateTextInput(session,    "reference_date", value    = st$reference_date)
      one_of <- function(x, allowed) {
        x <- as.character(unlist(x))
        if (length(x) == 1L && x %in% allowed) x else NULL
      }
      if (!is.null(v <- one_of(st$date_order, c("auto", "dmy", "mdy"))))
        updateSelectInput(session, "date_order", selected = v)

      # Study design (R20-D.H)
      design <- list(volume_units = one_of(st$volume_units, c("mm3", "cm3")),
                     randomisation_unit = one_of(st$randomisation_unit, c("mouse", "cage")),
                     cage_labels = one_of(st$cage_labels, c("study", "arm")),
                     control_group = {
                       cg <- as.character(unlist(st$control_group))
                       if (length(cg) == 1L && nzchar(cg)) cg else NULL
                     })
      design <- design[!vapply(design, is.null, logical(1L))]
      # Removal-reason categories: kept for when the mapping panel is drawn,
      # and applied to the selectors already on screen.
      rmap <- st$removal_reason_map
      if (is.list(rmap) && length(rmap)) {
        rmap <- rmap[vapply(rmap, function(x) length(x) == 1L &&
                              as.character(x) %in% removal_reason_categories,
                            logical(1L))]
        rmap <- lapply(rmap, as.character)
        design$removal_reason_map <- rmap
        rcol <- as.character(unlist(st$col_removal_reason))
        vals <- if (length(rcol) == 1L) removal_reason_values(rv$uploaded_data, rcol)
                else character(0)
        for (i in seq_along(vals)) if (!is.null(rmap[[vals[i]]]))
          updateSelectInput(session, paste0("reason_cat_", i), selected = rmap[[vals[i]]])
      }
      if (is.function(study_defaults)) {
        study_defaults(utils::modifyList(isolate(study_defaults()), design))
      }
      if (!is.null(design$volume_units))
        updateRadioButtons(session, "volume_units", selected = design$volume_units)
      if (!is.null(design$randomisation_unit))
        updateRadioButtons(session, "randomisation_unit", selected = design$randomisation_unit)
      if (!is.null(design$cage_labels))
        updateRadioButtons(session, "cage_labels", selected = design$cage_labels)
      if (!is.null(design$control_group))
        updateSelectInput(session, "control_group", selected = design$control_group)

      # Subset filter (format 2)
      if (is.list(preset$filters)) {
        chr <- function(x) as.character(unlist(x))
        rv$exclude_ids        <- chr(preset$filters$exclude_ids)
        rv$exclude_treatments <- chr(preset$filters$exclude_treatments)
        rv$exclude_cages      <- chr(preset$filters$exclude_cages)
        rv$processed_data     <- NULL
        rv$filter_applied_at  <- Sys.time()
      }

      saved_at <- preset$saved_at %||% "(unknown)"
      showNotification(
        paste("Preset loaded \u2014 saved at", saved_at),
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
        selectizeInput(ns("exclude_ids"), "Exclude animals",
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
      "what's currently applied. Filters apply to every analysis until you clear them.")
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

  # Populate choices whenever uploaded_data + column mappings change, and
  # after Apply / Clear / a preset load (filter_applied_at), which set the
  # selections from rv$exclude_*.
  observe({
    rv$filter_applied_at
    df <- rv$uploaded_data
    req(df)
    id_col <- rv$settings$col_id
    tx_col <- rv$settings$col_treatment
    cg_col <- rv$settings$col_cage

    # One entry per animal (treatment + ID + cage), labelled "3 (Vehicle,
    # cage C1)": a reused ear tag is a different animal in each arm.
    ids <- if (has_col(df, id_col) && has_col(df, tx_col)) {
      cg  <- if (has_col(df, cg_col)) df[[cg_col]]
      arm <- arm_labels(df, tx_col, rv$settings$col_dose)
      key <- animal_key(arm, df[[id_col]], cg)
      lab <- animal_label(arm, df[[id_col]], cg)
      u   <- !duplicated(key)
      idv <- as.character(df[[id_col]])[u]
      o   <- order(arm[u], suppressWarnings(as.numeric(idv)), idv)
      stats::setNames(key[u][o], lab[u][o])
    } else character(0)
    txs <- if (!is.null(tx_col) && nzchar(tx_col) && tx_col %in% colnames(df))
      sort(unique(as.character(df[[tx_col]]))) else character(0)
    # R22.3: with cages numbered within each arm, each arm's cage.
    cgs <- cage_filter_choices(df, rv$settings, isolate(rv$exclude_cages) %||% character(0))

    updateSelectizeInput(session, "exclude_ids",
                        choices = ids,
                        selected = intersect(isolate(rv$exclude_ids), unname(ids)),
                        server = TRUE)
    updateSelectizeInput(session, "exclude_treatments",
                        choices = txs,
                        selected = intersect(isolate(rv$exclude_treatments), txs),
                        server = TRUE)
    updateSelectizeInput(session, "exclude_cages",
                        choices = cgs$choices,
                        selected = cgs$selected,
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
  # when rv$uploaded_data changes (e.g., a new file).
  output$subset_filter_status <- renderUI({
    rv$filter_applied_at  # tick dependency
    n_ids <- length(rv$exclude_ids        %||% character(0))
    n_tx  <- length(rv$exclude_treatments %||% character(0))
    n_cg  <- length(rv$exclude_cages      %||% character(0))
    if (n_ids + n_tx + n_cg == 0L) {
      return(tags$span(
        style = "color:#6c757d; margin-left: 14px; font-size: 12px;",
        "No filter applied \u2014 full dataset is shown below."))
    }
    if (!is.null(rv$uploaded_data) && nrow(rv$uploaded_data) > 0L) {
      n_full <- nrow(rv$uploaded_data)
      n_filt <- nrow(apply_subset_filter(rv$uploaded_data, rv))
      tags$span(
        style = "color: #155724; margin-left: 14px; font-size: 12px;",
        icon("check-circle"),
        sprintf(" Active filter \u2014 %d of %d rows shown (%d excluded).",
                n_filt, n_full, n_full - n_filt))
    } else {
      tags$span(
        style = "color: #155724; margin-left: 14px; font-size: 12px;",
        icon("check-circle"),
        sprintf(" Active filter \u2014 %d animal(s), %d treatment(s), %d cage(s) excluded.",
                n_ids, n_tx, n_cg))
    }
  })
}
