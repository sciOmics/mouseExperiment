#' Server Modules for mouseExperimentDashboard

#' Data Upload Server Logic
#' @keywords internal
data_upload_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
  
  # Reactive value to store uploaded data
  rv$uploaded_data <- NULL
  rv$data_valid <- FALSE
  rv$detected_cols <- NULL
  
  # Demo data loader
  observeEvent(input$load_demo_data, {
    req(nzchar(input$demo_dataset))

    demo_files <- list(
      combo = list(
        data     = "combo_treatment_synthetic_data.csv",
        ref_date = "2025-03-24",
        label    = "Combo Treatment demo"
      ),
      dose_levels = list(
        data     = "dose_levels_synthetic_data.csv",
        ref_date = "2025-03-24",
        label    = "Dose Levels demo"
      ),
      necrotic = list(
        data     = "necrotic_synthetic_data.csv",
        ref_date = "",
        label    = "Necrotic Data demo"
      ),
      weight = list(
        data     = "weight_synthetic_data.csv",
        ref_date = "2025-03-01",
        label    = "Weight Data demo"
      )
    )

    chosen   <- demo_files[[input$demo_dataset]]
    pkg_name <- "mouseExperimentDashboard"

    # Locate files: installed package first, then development source tree
    find_demo <- function(filename) {
      p <- system.file("sample_data", filename, package = pkg_name)
      if (nzchar(p)) return(p)
      # fallback for devtools::load_all() / running from source
      local_p <- file.path("inst", "sample_data", filename)
      if (file.exists(local_p)) return(local_p)
      NULL
    }

    data_path <- find_demo(chosen$data)
    if (is.null(data_path)) {
      output$demo_load_msg <- renderUI({
        status_alert("Demo file not found. Please reinstall the package.", "danger")
      })
      return()
    }

    # Load main data file
    result <- validate_uploaded_file(data_path, "csv")
    if (!result$valid) {
      output$demo_load_msg <- renderUI({
        status_alert(result$message, "danger")
      })
      return()
    }

    rv$uploaded_data  <- result$data
    rv$data_valid     <- TRUE
    rv$detected_cols  <- detect_columns(colnames(result$data))

    log_upload(session, chosen$label, nrow(result$data), ncol(result$data))

    # Set reference date for this dataset
    updateTextInput(session, "reference_date", value = chosen$ref_date)

    output$file_validation_msg <- renderUI({
      status_alert(paste(chosen$label, "loaded:",
                nrow(result$data), "rows,",
                ncol(result$data), "cols."), "success")
    })

    output$demo_load_msg <- renderUI({
      div(class="alert alert-success", style="margin:0;",
          icon("check"), " Loaded!")
    })

  })

  # File upload handler
  observeEvent(input$data_file, {
    req(input$data_file)
    file_info <- input$data_file
    file_path <- file_info$datapath
    # Extract file extension manually
    file_parts <- strsplit(file_info$name, "\\.")[[1]]
    file_ext <- tolower(file_parts[length(file_parts)])
    
    # Validate and load file
    result <- validate_uploaded_file(file_path, file_ext)
    
    if (result$valid) {
      rv$uploaded_data <- result$data
      rv$data_valid <- TRUE
      rv$detected_cols <- detect_columns(colnames(result$data))

      log_upload(session, file_info$name, nrow(result$data), ncol(result$data))

      # Display validation message
      output$file_validation_msg <- renderUI({
        status_alert(result$message, "success")
      })
      
    } else {
      rv$data_valid <- FALSE
      rv$uploaded_data <- NULL
      output$file_validation_msg <- renderUI({
        status_alert(result$message, "danger")
      })
    }
  })
  
  # File info output
  output$file_info <- renderText({
    if (is.null(rv$uploaded_data)) {
      "No file loaded yet"
    } else {
      paste("Rows:", nrow(rv$uploaded_data), "\n",
            "Columns:", ncol(rv$uploaded_data), "\n",
            "Column names:", paste(colnames(rv$uploaded_data), collapse = ", "))
    }
  })
  
  # Column mapping UI
  output$column_mapping_ui <- renderUI({
    if (!rv$data_valid) return(NULL)
    column_mapping_ui(colnames(rv$uploaded_data), rv$detected_cols, ns = session$ns)
  })
  
  # ---------------------------------------------------------------------------
  # Reactive: compute preview data with calculated columns applied
  # Mirrors the same transformations used in tumor_growth_server so the preview
  # shows Volume and Day exactly as analysis will see them.
  # ---------------------------------------------------------------------------
  preview_data <- reactive({
    req(rv$uploaded_data, rv$data_valid)
    data <- rv$uploaded_data
    
    # Apply date -> Day calculation if selected
    if (!is.null(input$time_method) && input$time_method == "from_date" &&
        !is.null(input$col_date_for_day) && input$col_date_for_day != "") {
      tryCatch({
        date_info <- detect_date_format(data, input$col_date_for_day,
                                        input$reference_date %||% "")
        data <- mouseExperiment::calculate_dates(
          date_info$data,
          start_date  = date_info$start_date,
          date_column = input$col_date_for_day,
          date_format = date_info$date_format
        )
      }, error = function(e) {
        message("[TG Preview] Date calculation failed: ", conditionMessage(e))
      })
    }
    
    # Apply volume calculation if selected
    if (!is.null(input$volume_method) && input$volume_method == "calculate" &&
        !is.null(input$col_length) && input$col_length != "" &&
        !is.null(input$col_width)  && input$col_width  != "") {
      tryCatch({
        data <- mouseExperiment::calculate_volume(
          data,
          length_column = input$col_length,
          width_column  = input$col_width,
          formula       = input$volume_formula
        )
      }, error = function(e) {
        message("[TG Preview] Volume calculation failed: ", conditionMessage(e))
      })
    }
    
    data
  })

  # Detect rows where the user has width > length before auto-correction.
  # Returns a named integer vector: names are the 1-based row indices, values
  # are those same indices (for easy display and DT highlighting).
  transposed_rows <- reactive({
    req(rv$uploaded_data, rv$data_valid)
    if (!isTRUE(input$volume_method == "calculate")) return(integer(0))
    len_col <- input$col_length %||% ""
    wid_col <- input$col_width  %||% ""
    if (!nzchar(len_col) || !nzchar(wid_col)) return(integer(0))
    d <- rv$uploaded_data
    if (!len_col %in% colnames(d) || !wid_col %in% colnames(d)) return(integer(0))
    len_vals <- suppressWarnings(as.numeric(d[[len_col]]))
    wid_vals <- suppressWarnings(as.numeric(d[[wid_col]]))
    which(!is.na(len_vals) & !is.na(wid_vals) & wid_vals > len_vals)
  })

  # ---------------------------------------------------------------------------
  # Dataset Summary
  # ---------------------------------------------------------------------------
  output$data_summary_ui <- renderUI({
    if (!rv$data_valid) return(NULL)
    
    data    <- preview_data()
    tx_col  <- if (!is.null(input$col_treatment) && input$col_treatment != "") input$col_treatment else NULL
    dose_col <- if (!is.null(input$col_dose) && input$col_dose != "" && input$col_dose != "None") input$col_dose else NULL
    cage_col <- if (!is.null(input$col_cage) && input$col_cage != "") input$col_cage else NULL
    id_col  <- if (!is.null(input$col_id) && input$col_id != "") input$col_id else NULL
    
    # Determine the time column used in the preview
    time_col <- if (!is.null(input$time_method) && input$time_method == "from_date") "Day"
                else if (!is.null(input$col_day) && input$col_day != "") input$col_day
                else NULL
    
    vol_col <- if (!is.null(input$volume_method) && input$volume_method == "calculate") "Volume"
               else if (!is.null(input$col_volume) && input$col_volume != "") input$col_volume
               else NULL
    
    # --- Treatment groups ---
    tx_box <- if (has_col(data, tx_col)) {
      if (has_col(data, dose_col)) {
        # Show treatment × dose breakdown as a table
        td <- unique(data[, c(tx_col, dose_col)])
        td[[tx_col]]   <- as.character(td[[tx_col]])
        td[[dose_col]] <- as.character(td[[dose_col]])
        td <- td[order(td[[tx_col]], td[[dose_col]]), ]
        # Collapse doses per treatment into a single row
        dose_per_tx <- tapply(td[[dose_col]], td[[tx_col]], function(d)
          paste(sort(unique(d)), collapse = ", "))
        tx_rows <- lapply(sort(names(dose_per_tx)), function(g)
          tags$tr(tags$td(g), tags$td(dose_per_tx[[g]]))
        )
        div(class = "panel panel-default",
          div(class = "panel-heading", tags$b("Treatment Groups")),
          div(class = "panel-body",
            tags$table(
              class = "table table-condensed table-bordered",
              style = "width:auto;",
              tags$thead(tags$tr(tags$th("Treatment"), tags$th("Dose(s)"))),
              tags$tbody(tx_rows)
            )
          )
        )
      } else {
        groups <- sort(unique(na.omit(as.character(data[[tx_col]]))))
        div(class = "panel panel-default",
          div(class = "panel-heading", tags$b("Treatment Groups")),
          div(class = "panel-body",
            tags$ul(lapply(groups, tags$li))
          )
        )
      }
    } else div(class = "text-muted", "(Treatment column not configured)")
    
    # --- Cage IDs removed --- #

    # --- Day / date range ---
    time_box <- if (has_col(data, time_col)) {
      vals <- suppressWarnings(as.numeric(data[[time_col]]))
      vals <- vals[!is.na(vals)]
      if (length(vals) > 0) {
        div(class = "panel panel-default",
          div(class = "panel-heading", tags$b("Day Range")),
          div(class = "panel-body",
            p(paste0("Min: ", min(vals), " \u2014 Max: ", max(vals),
                    " (", length(unique(vals)), " unique timepoints)"))
          )
        )
      } else NULL
    } else NULL
    
    # --- Mice per group ---
    mice_box <- if (has_col(data, tx_col) && has_col(data, id_col)) {

      has_dose <- has_col(data, dose_col)

      inner <- if (has_col(data, cage_col)) {
        # Cage × treatment (× dose) table
        key_cols <- c(id_col, cage_col, tx_col, if (has_dose) dose_col)
        ct <- unique(data[, key_cols])
        group_cols <- c("Cage", "Treatment", if (has_dose) "Dose")
        ct_renamed <- data.frame(
          Cage      = as.character(ct[[cage_col]]),
          Treatment = as.character(ct[[tx_col]]),
          stringsAsFactors = FALSE
        )
        if (has_dose) ct_renamed$Dose <- as.character(ct[[dose_col]])
        ct_counts <- as.data.frame(table(ct_renamed[, group_cols, drop = FALSE]))
        ct_counts <- ct_counts[ct_counts$Freq > 0, ]
        ct_counts <- ct_counts[do.call(order, ct_counts[, group_cols, drop = FALSE]), ]
        ct_rows <- lapply(seq_len(nrow(ct_counts)), function(i) {
          cells <- lapply(group_cols, function(col) tags$td(as.character(ct_counts[[col]][i])))
          do.call(tags$tr, c(cells, list(tags$td(ct_counts$Freq[i]))))
        })
        header_cells <- lapply(c(group_cols, "Mice"), tags$th)
        tags$table(
          class = "table table-condensed table-bordered",
          style = "width:auto;",
          tags$thead(do.call(tags$tr, header_cells)),
          tags$tbody(ct_rows)
        )
      } else {
        # Per-treatment (+ dose) when no cage column configured
        key_cols <- c(id_col, tx_col, if (has_dose) dose_col)
        ct <- unique(data[, key_cols])
        if (has_dose) {
          group_df <- data.frame(
            Treatment = as.character(ct[[tx_col]]),
            Dose      = as.character(ct[[dose_col]]),
            stringsAsFactors = FALSE
          )
          counts <- as.data.frame(table(group_df))
          counts <- counts[counts$Freq > 0, ]
          counts <- counts[order(counts$Treatment, counts$Dose), ]
          tx_rows <- lapply(seq_len(nrow(counts)), function(i)
            tags$tr(tags$td(counts$Treatment[i]),
                    tags$td(counts$Dose[i]),
                    tags$td(counts$Freq[i]))
          )
          tags$table(
            class = "table table-condensed table-bordered",
            style = "width:auto;",
            tags$thead(tags$tr(tags$th("Treatment"), tags$th("Dose"), tags$th("Mice"))),
            tags$tbody(tx_rows)
          )
        } else {
          tx_counts <- tapply(
            as.character(data[[id_col]]),
            as.character(data[[tx_col]]),
            function(x) length(unique(x))
          )
          tx_rows <- lapply(sort(names(tx_counts)), function(g)
            tags$tr(tags$td(g), tags$td(tx_counts[[g]]))
          )
          tags$table(
            class = "table table-condensed table-bordered",
            style = "width:auto;",
            tags$thead(tags$tr(tags$th("Treatment"), tags$th("Unique Mice"))),
            tags$tbody(tx_rows)
          )
        }
      }

      div(class = "panel panel-default",
        div(class = "panel-heading", tags$b("Cage Info")),
        div(class = "panel-body", inner)
      )
    } else NULL
    
    # --- Transposed measurement warning ---
    tr           <- transposed_rows()
    raw_data     <- rv$uploaded_data
    len_col_raw  <- if (!is.null(input$col_length) && input$col_length != "") input$col_length else NULL
    wid_col_raw  <- if (!is.null(input$col_width)  && input$col_width  != "") input$col_width  else NULL
    id_col_disp  <- if (!is.null(input$col_id)     && input$col_id     != "") input$col_id     else NULL
    time_col_disp <- if (has_col(raw_data, time_col)) time_col else NULL

    transposed_box <- if (length(tr) > 0 && !is.null(len_col_raw) && !is.null(wid_col_raw)) {
      tbl_rows <- lapply(tr, function(i) {
        len_v <- raw_data[[len_col_raw]][i]
        wid_v <- raw_data[[wid_col_raw]][i]
        id_v  <- if (has_col(raw_data, id_col_disp))   raw_data[[id_col_disp]][i]   else ""
        day_v <- if (has_col(raw_data, time_col_disp)) raw_data[[time_col_disp]][i] else ""
        tags$tr(
          tags$td(i),
          tags$td(if (nzchar(trimws(as.character(id_v))))  as.character(id_v)  else tags$em("\u2014")),
          tags$td(if (nzchar(trimws(as.character(day_v)))) as.character(day_v) else tags$em("\u2014")),
          tags$td(style = "color:#c0392b; font-weight:bold;", round(as.numeric(wid_v), 3)),
          tags$td(round(as.numeric(len_v), 3)),
          tags$td(style = "color:#27ae60;",
                  paste(round(as.numeric(len_v), 3), "/", round(as.numeric(wid_v), 3)))
        )
      })
      div(class = "alert alert-warning", style = "margin-top: 10px;",
        icon("exclamation-triangle"), " ",
        strong(paste(length(tr), "measurement row(s) auto-corrected: width value exceeded length")),
        tags$br(),
        tags$small(
          "The dimensions were automatically swapped so length \u2265 width before calculating volume. ",
          "These rows are highlighted in yellow in the Data Preview below."
        ),
        tags$br(), tags$br(),
        div(style = "overflow-x:auto;",
          tags$table(
            class = "table table-condensed table-bordered",
            style = "width:auto; background:#fff;",
            tags$thead(
              tags$tr(
                tags$th("Row #"),
                tags$th(if (!is.null(id_col_disp))   id_col_disp   else "ID"),
                tags$th(if (!is.null(time_col_disp)) time_col_disp else "Day"),
                tags$th(style = "color:#c0392b;",
                        paste0(wid_col_raw, " (was larger — swapped to width)")),
                tags$th(paste0(len_col_raw, " (was smaller — swapped to length)")),
                tags$th("Corrected: Length / Width")
              )
            ),
            tags$tbody(tbl_rows)
          )
        )
      )
    } else NULL

    tagList(
      wellPanel(
        h4("Dataset Summary"),
        p(strong(paste(nrow(data), "rows \u00d7", ncol(data), "columns"))),
        transposed_box,
        fluidRow(
          column(4, tx_box),
          column(4, time_box),
          column(4, mice_box)
        )
      )
    )
  })
  
  # Data preview
  output$data_preview_ui <- renderUI({
    if (!rv$data_valid) return(NULL)
    data_preview_ui(ns = session$ns)
  })
  
  output$preview_summary <- renderText({
    if (!rv$data_valid) return("")
    data <- preview_data()
    paste("Dataset:", nrow(data), "rows ×", ncol(data), "columns")
  })
  
  output$data_table <- renderDT({
    req(rv$data_valid)
    tr   <- transposed_rows()
    opts <- list(scrollX = TRUE, scrollY = "600px", paging = FALSE, dom = "ti")
    if (length(tr) > 0) {
      js_idx <- paste(tr - 1L, collapse = ",")
      opts$rowCallback <- DT::JS(sprintf(
        "function(row, data, displayNum, displayIndex, dataIndex) {
           if ([%s].indexOf(dataIndex) >= 0) {
             $(row).css({'background-color': '#fff3cd', 'font-weight': 'bold'});
             $(row).attr('title',
               'Auto-corrected: width exceeded length \\u2014 dimensions were swapped before volume calculation');
           }
         }", js_idx))
    }
    datatable(preview_data(), options = opts)
  })
  
  output$download_preview_data <- downloadHandler(
    filename = function() paste0("processed_data_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content  = function(file) write.csv(preview_data(), file, row.names = FALSE)
  )
  
  # Reset button
  observeEvent(input$reset_data, {
    # Clear reactive values
    rv$uploaded_data <- NULL
    rv$data_valid <- FALSE
    rv$detected_cols <- NULL
    
    # Clear all outputs by rendering empty content
    output$file_validation_msg <- renderUI({
      div(class = "alert alert-info",
          icon("info-circle"), "Ready to upload a new file")
    })
    
    output$file_info <- renderText({
      "No file loaded yet"
    })
    
    output$column_mapping_ui <- renderUI({
      NULL
    })
    
    output$data_preview_ui <- renderUI({
      NULL
    })
    
    output$data_summary_ui <- renderUI({
      NULL
    })
    
    output$preview_summary <- renderText({
      ""
    })
    
    output$data_table <- renderDT({
      datatable(data.frame(Message = "No data loaded"), 
                options = list(scrollX = TRUE))
    })
    
    # Reset the file input widget
    shinyjs::reset(session$ns("data_file"))
  })
  

  # Proceed button is always clickable - validation enforced at analysis time
  # No need to disable/enable the button

  # ---------------------------------------------------------------------------
  # Live validation feedback on column mapping inputs (shinyFeedback)
  # ---------------------------------------------------------------------------

  # Helper: show success/warning on a selectInput depending on whether a
  # non-empty value is selected.
  validate_col_input <- function(input_id) {
    val <- input[[input_id]]
    if (!is.null(val) && nzchar(val)) {
      shinyFeedback::hideFeedback(input_id)
      shinyFeedback::showFeedbackSuccess(input_id, text = "")
    } else {
      shinyFeedback::hideFeedback(input_id)
      shinyFeedback::showFeedbackWarning(input_id, text = "Required")
    }
  }

  # Watch required column selectors (only after data is uploaded)
  observe({
    req(rv$data_valid)
    validate_col_input("col_id")
    validate_col_input("col_treatment")
  })

  observe({
    req(rv$data_valid)
    if ((input$time_method %||% "") == "direct_day") {
      validate_col_input("col_day")
    } else {
      validate_col_input("col_date_for_day")
    }
  })

  observe({
    req(rv$data_valid)
    if ((input$volume_method %||% "") == "direct") {
      validate_col_input("col_volume")
    } else {
      validate_col_input("col_length")
      validate_col_input("col_width")
    }
  })

  # ---------------------------------------------------------------------------
  # Sync column mapping inputs into rv$settings for cross-module access.
  # Other modules (TG, survival, synergy, etc.) read rv$settings$* instead of
  # input$col_* so they work correctly under moduleServer() namespacing.
  # ---------------------------------------------------------------------------
  observe({
    rv$settings <- list(
      col_treatment     = input$col_treatment %||% "",
      col_id            = input$col_id %||% "",
      col_day           = input$col_day %||% "",
      col_date_for_day  = input$col_date_for_day %||% "",
      col_volume        = input$col_volume %||% "",
      col_length        = input$col_length %||% "",
      col_width         = input$col_width %||% "",
      col_cage          = input$col_cage %||% "",
      col_dose          = input$col_dose %||% "",
      col_survival_event = input$col_survival_event %||% "",
      time_method       = input$time_method %||% "direct_day",
      volume_method     = input$volume_method %||% "direct",
      volume_formula    = input$volume_formula %||% "ellipsoid",
      reference_date    = input$reference_date %||% ""
    )
  })

  # Session save / restore handlers
  setup_session_save_restore(input, output, session, rv)

  })
}
