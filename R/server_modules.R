#' Server Modules for mouseExperimentDashboard

#' Data Upload Server Logic
#' @keywords internal
data_upload_server <- function(input, output, session, rv) {
  
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
        div(class="alert alert-danger", icon("exclamation-circle"),
            "Demo file not found. Please reinstall the package.")
      })
      return()
    }

    # Load main data file
    result <- validate_uploaded_file(data_path, "csv")
    if (!result$valid) {
      output$demo_load_msg <- renderUI({
        div(class="alert alert-danger", icon("exclamation-circle"), result$message)
      })
      return()
    }

    rv$uploaded_data  <- result$data
    rv$data_valid     <- TRUE
    rv$detected_cols  <- detect_columns(colnames(result$data))

    # Set reference date for this dataset
    updateTextInput(session, "reference_date", value = chosen$ref_date)

    output$file_validation_msg <- renderUI({
      div(class="alert alert-success",
          icon("check"),
          paste(chosen$label, "loaded:",
                nrow(result$data), "rows,",
                ncol(result$data), "cols."))
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
      
      # Display validation message
      output$file_validation_msg <- renderUI({
        div(class = "alert alert-success",
            icon("check"), result$message)
      })
      
    } else {
      cat("DEBUG: Validation failed - ", result$message, "\n")
      rv$data_valid <- FALSE
      rv$uploaded_data <- NULL
      output$file_validation_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), result$message)
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
    column_mapping_ui(colnames(rv$uploaded_data), rv$detected_cols)
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
        first_val    <- as.character(data[1, input$col_date_for_day])
        ref_date_str <- input$reference_date %||% ""

        # Detect format and normalise both data dates and start_date
        if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_val)) {
          fmt <- "%m/%d/%Y"
        } else if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", first_val)) {
          fmt <- "%Y-%m-%d"
        } else if (grepl("^\\d{1,2}-[A-Za-z]{3}$", first_val)) {
          fmt <- "%d-%b"
        } else if (grepl("^[A-Za-z]+\\s+\\d{1,2}$", first_val)) {
          # "Month Day" (e.g. "March 1") — inject year from reference date
          ref_parsed <- tryCatch(
            as.Date(ref_date_str,
                    tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d-%b-%Y",
                                   "%B %d %Y", "%B %d, %Y")),
            error = function(e) Sys.Date()
          )
          if (is.na(ref_parsed)) ref_parsed <- Sys.Date()
          ref_year <- as.integer(format(ref_parsed, "%Y"))
          data[[input$col_date_for_day]] <-
            paste(data[[input$col_date_for_day]], ref_year)
          ref_date_str <- format(ref_parsed, "%B %d %Y")
          fmt <- "%B %d %Y"
        } else {
          fmt <- NULL
        }

        data <- mouseExperiment::calculate_dates(
          data,
          start_date  = ref_date_str,
          date_column = input$col_date_for_day,
          date_format = fmt
        )
      }, error = function(e) NULL)
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
      }, error = function(e) NULL)
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
    tx_box <- if (!is.null(tx_col) && tx_col %in% colnames(data)) {
      if (!is.null(dose_col) && dose_col %in% colnames(data)) {
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
    time_box <- if (!is.null(time_col) && time_col %in% colnames(data)) {
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
    mice_box <- if (!is.null(tx_col) && tx_col %in% colnames(data) &&
                    !is.null(id_col)  && id_col  %in% colnames(data)) {

      has_dose <- !is.null(dose_col) && dose_col %in% colnames(data)

      inner <- if (!is.null(cage_col) && cage_col %in% colnames(data)) {
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
    time_col_disp <- if (!is.null(time_col) && !is.null(time_col) && time_col %in% colnames(raw_data)) time_col else NULL

    transposed_box <- if (length(tr) > 0 && !is.null(len_col_raw) && !is.null(wid_col_raw)) {
      tbl_rows <- lapply(tr, function(i) {
        len_v <- raw_data[[len_col_raw]][i]
        wid_v <- raw_data[[wid_col_raw]][i]
        id_v  <- if (!is.null(id_col_disp)   && id_col_disp   %in% colnames(raw_data)) raw_data[[id_col_disp]][i]   else ""
        day_v <- if (!is.null(time_col_disp) && time_col_disp %in% colnames(raw_data)) raw_data[[time_col_disp]][i] else ""
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
    data_preview_ui()
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
    shinyjs::reset("data_file")
  })
  

  
  # Proceed button is always clickable - validation enforced at analysis time
  # No need to disable/enable the button
}

#' Tumor Growth Analysis Server Logic
#' @keywords internal
tumor_growth_server <- function(input, output, session, rv) {
  
  rv$tg_results <- NULL
  rv$tg_running <- FALSE
  
  # Update reference group choices when data is available or column mapping changes.
  # When a dose column is also mapped, group choices are shown as "Treatment Dose"
  # combinations so each dose level is independently selectable as the reference.
  observe({
    if (!is.null(rv$uploaded_data) && rv$data_valid) {
      # Prefer the user-selected treatment column; fall back to auto-detected one
      treatment_col <- if (!is.null(input$col_treatment) && nzchar(input$col_treatment))
                         input$col_treatment
                       else
                         rv$detected_cols$treatment
      dose_col <- if (!is.null(input$col_dose) && nzchar(input$col_dose))
                    input$col_dose
                  else NULL

      cat("DEBUG: Updating reference group choices. treatment_col=", treatment_col,
          " dose_col=", if (is.null(dose_col)) "NULL" else dose_col, "\n")

      if (!is.null(treatment_col) && nzchar(treatment_col)) {
        df <- rv$uploaded_data

        if (!is.null(dose_col) && dose_col %in% colnames(df)) {
          # Build unique treatment+dose combination labels
          treatment_values <- sort(unique(na.omit(
            paste(as.character(df[[treatment_col]]),
                  as.character(df[[dose_col]]), sep = " ")
          )))
        } else {
          treatment_values <- sort(unique(na.omit(as.character(df[[treatment_col]]))))
        }

        if (length(treatment_values) > 0) {
          new_choices <- c("-- Select a reference group --" = "",
                          setNames(treatment_values, treatment_values))
          default_ref <- pick_reference_group(treatment_values)
          updateSelectInput(session, "tg_reference", choices = new_choices,
                            selected = default_ref)
          cat("  Reference group choices updated. Default:", default_ref, "\n")
        }
      }
    }
  })
  
  # Run analysis button handler
  observeEvent(input$run_tg_analysis, {
    cat("DEBUG: Run Analysis button clicked\n")
    req(rv$uploaded_data, rv$data_valid)
    cat("DEBUG: Data requirements met, starting analysis\n")
    
    # Validate reference group selection
    if (input$tg_reference == "" || is.null(input$tg_reference)) {
      output$tg_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error: Please select a reference group before running analysis")
      })
      return()
    }
    
    rv$tg_running <- TRUE
    output$tg_status_msg <- renderUI({
      div(class = "alert alert-info",
          icon("spinner", class = "fa-spin"), "Running analysis...")
    })
    
    # Prepare data
    data <- rv$uploaded_data
    
    # Handle time column - either use direct day or calculate from date
    if (input$time_method == "from_date") {
      # Calculate days from start date
      start_date <- input$reference_date
      
      # Detect date format from the actual data in the date column
      first_date <- as.character(data[1, input$col_date_for_day])
      
      # Determine date format by checking the pattern
      date_format <- NULL
      if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_date)) {
        date_format <- "%m/%d/%Y"  # MM/DD/YYYY format
      } else if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", first_date)) {
        date_format <- "%Y-%m-%d"  # YYYY-MM-DD format
      } else if (grepl("^\\d{1,2}-[A-Za-z]{3}$", first_date)) {
        date_format <- "%d-%b"     # DD-Mon format
      } else if (grepl("^[A-Za-z]+\\s+\\d{1,2}$", first_date)) {
        # "Month Day" format (e.g. "March 1", "January 15") – no year in data.
        # Extract year from the reference date so we can build complete timestamps.
        ref_date_parsed <- tryCatch(
          as.Date(start_date,
                  tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d-%b-%Y",
                                 "%B %d %Y", "%B %d, %Y")),
          error = function(e) Sys.Date()
        )
        if (is.na(ref_date_parsed)) ref_date_parsed <- Sys.Date()
        ref_year <- as.integer(format(ref_date_parsed, "%Y"))
        # Append year to every data date:  "March 1" -> "March 1 2025"
        data[[input$col_date_for_day]] <-
          paste(data[[input$col_date_for_day]], ref_year)
        # Reformat start_date to the same pattern so calculate_dates
        # parses both sides identically
        start_date  <- format(ref_date_parsed, "%B %d %Y")
        date_format <- "%B %d %Y"
      }
      
      data <- mouseExperiment::calculate_dates(
        data,
        start_date = start_date,
        date_column = input$col_date_for_day,
        date_format = date_format
      )
      time_col <- "Day"  # calculate_dates produces a "Day" column
    } else {
      time_col <- input$col_day
    }
    
    # Apply volume calculation if needed
    if (input$volume_method == "calculate") {
      len_col <- input$col_length %||% ""
      wid_col <- input$col_width  %||% ""
      if (!nzchar(len_col) || !nzchar(wid_col)) {
        output$tg_status_msg <- renderUI({
          div(class = "alert alert-danger",
              icon("exclamation-circle"),
              "Please select both Length and Width columns before running analysis.")
        })
        rv$tg_running <- FALSE
        return()
      }
      if (!len_col %in% colnames(data) || !wid_col %in% colnames(data)) {
        output$tg_status_msg <- renderUI({
          div(class = "alert alert-danger",
              icon("exclamation-circle"),
              paste0("Column not found in data: '",
                     if (!len_col %in% colnames(data)) len_col else wid_col, "'"))
        })
        rv$tg_running <- FALSE
        return()
      }
      data <- mouseExperiment::calculate_volume(
        data,
        length_column = len_col,
        width_column  = wid_col,
        formula       = input$volume_formula
      )
      vol_col <- "Volume"
    } else {
      vol_col <- input$col_volume
    }
    
    # Store preprocessed data and resolved column names for reuse by other modules
    rv$processed_data <- data
    rv$time_col       <- time_col
    rv$vol_col        <- vol_col
    rv$tx_col         <- input$col_treatment
    rv$id_col         <- input$col_id
    rv$cage_col       <- if (!is.null(input$col_cage) && input$col_cage != "") input$col_cage else NULL
    rv$dose_col       <- if (!is.null(input$col_dose) && input$col_dose != "") input$col_dose else NULL

    # When a dose column is mapped, combine treatment + dose into a single
    # group label so each dose level is treated as its own group by all analyses.
    if (!is.null(rv$dose_col)) {
      data[["TreatmentGroup"]] <- paste(data[[input$col_treatment]],
                                        data[[rv$dose_col]], sep = " ")
      rv$processed_data <- data
      rv$tx_col         <- "TreatmentGroup"
    }

    # Prepare column names
    col_names <- list(
      time_column = time_col,
      volume_column = vol_col,
      id_column = input$col_id,
      treatment_column = rv$tx_col,
      cage_column = if (input$col_cage != "") input$col_cage else NULL
    )
    
    # Run analysis
    result <- safe_analysis({
      mouseExperiment::tumor_growth_statistics(
        df = data,
        time_column = col_names$time_column,
        volume_column = col_names$volume_column,
        id_column = col_names$id_column,
        treatment_column = col_names$treatment_column,
        dose_column = rv$dose_col,
        cage_column = col_names$cage_column,
        model_type = input$tg_model_type,
        transform = input$tg_transform,
        polynomial_degree = 1L,
        random_effects_specification = input$tg_random_effects,
        handle_cage_effects = input$tg_cage_effect,
        reference_group = input$tg_reference,  # Already validated to not be empty
        extrapolation_points = if (!is.null(input$tg_extrapolation_points) &&
                                    !is.na(input$tg_extrapolation_points))
                                 as.integer(input$tg_extrapolation_points) else 0L,
        return_model = TRUE,
        include_diagnostics = TRUE,
        plots = FALSE
      )
    })
    
    # Generate plots from the results
    if (!isTRUE(result$error)) {
      tryCatch({
        # Plot 1: Tumor growth curves from original data.
        # Enable extrapolation when extrapolation_points > 0.
        # If a survival/censor column is mapped, use it directly.
        # Otherwise synthesize one by flagging mice that disappear before
        # the study's final timepoint.
        if (!is.null(result$growth_rates)) {
          surv_col_mapped <- !is.null(input$col_survival_event) &&
                             input$col_survival_event != "" &&
                             input$col_survival_event != "None" &&
                             input$col_survival_event %in% colnames(data)
          extra_pts <- if (!is.null(input$tg_extrapolation_points) &&
                           !is.na(input$tg_extrapolation_points))
                         as.integer(input$tg_extrapolation_points) else 0L
          do_extrapolate <- extra_pts > 0L

          # Resolve which data frame + censor column to use
          plot_data  <- data
          surv_col   <- NULL
          if (do_extrapolate) {
            if (surv_col_mapped) {
              surv_col <- input$col_survival_event
            } else {
              syn <- synthesize_censor_column(
                data     = data,
                id_col   = col_names$id_column,
                time_col = col_names$time_column
              )
              if (syn$n_flagged > 0L) {
                plot_data <- syn$data
                surv_col  <- syn$col
                cat(sprintf("INFO: Auto-detected %d early-dropout mice; synthesized censor column '%s'\n",
                            syn$n_flagged, syn$col))
              } else {
                # All mice survived to study end — nothing to extrapolate
                do_extrapolate <- FALSE
              }
            }
          }

          result$tumor_growth_plot <- mouseExperiment::plot_tumor_growth(
            df = plot_data,
            volume_column = col_names$volume_column,
            day_column = col_names$time_column,
            treatment_column = if (!is.null(rv$dose_col)) input$col_treatment else col_names$treatment_column,
            dose_column = rv$dose_col,
            ID_column = col_names$id_column,
            cage_column = col_names$cage_column,
            survival_column = surv_col,
            extrapolate_volumes = do_extrapolate,
            extrapolation_points = if (extra_pts > 0L) extra_pts else "all",
            group_summary_line = TRUE
          ) + ggplot2::theme_classic()
        }
        
        # Plot 2: Growth rate distribution
        if (!is.null(result$growth_rates) && nrow(result$growth_rates) > 0) {
          result$growth_rate_plot <- mouseExperiment::plot_growth_rate(
            growth_data = result$growth_rates,
            show_mean = TRUE,
            error_bar_type = "SEM"
          ) + ggplot2::theme_classic()
        }
        
        # Plot 3: AUC / treatment effects plot
        # When model_type="auc" the result contains auc_analysis$individual with
        # per-subject AUC values — use the package's plot_auc() to match the demo.
        # For lme4 models fall back to a bar chart of adjusted means.
        if (!is.null(result$auc_analysis) && !is.null(result$auc_analysis$individual)) {
          tryCatch({
            extra_pts2 <- if (!is.null(input$tg_extrapolation_points) &&
                              !is.na(input$tg_extrapolation_points))
                            as.integer(input$tg_extrapolation_points) else 0L
            auc_ind <- result$auc_analysis$individual
            # Safety: if the Group column values don't match known treatment groups
            # (can happen when IDs with underscores confuse the package's string-split),
            # re-derive Group by joining on the original ID column.
            known_groups <- unique(as.character(data[[col_names$treatment_column]]))
            if (!all(auc_ind$Group %in% known_groups)) {
              id_tx <- unique(data[, c(col_names$id_column, col_names$treatment_column),
                                   drop = FALSE])
              merged <- merge(auc_ind, id_tx,
                              by.x = "ID", by.y = col_names$id_column,
                              all.x = TRUE)
              new_grp <- merged[[col_names$treatment_column]]
              if (!is.null(new_grp) && !all(is.na(new_grp))) {
                auc_ind$Group    <- as.character(new_grp)
                auc_ind$Treatment <- as.character(new_grp)
                result$auc_analysis$individual <- auc_ind
              }
            }
            result$auc_plot <- mouseExperiment::plot_auc(
              result$auc_analysis$individual,
              show_mean = TRUE,
              error_bar_type = "CI",
              extrapolated_column = if (extra_pts2 > 0L &&
                                       "Extrapolated" %in% colnames(result$auc_analysis$individual))
                                      "Extrapolated" else NULL
            ) + ggplot2::theme_classic()
          }, error = function(e) {
            cat("WARNING: plot_auc() failed:", e$message, "\n")
          })
        }
        if (is.null(result$auc_plot) && !is.null(result$treatment_effects)) {
          effects_df <- as.data.frame(result$treatment_effects)
          if (nrow(effects_df) > 0) {
            # treatment_effects uses "Group" (lme4) or "Treatment" (AUC) as the label column;
            # do NOT fall back to rownames() since rbind() resets them to "1","2","3"...
            label_col <- intersect(c("Group", "Treatment"), colnames(effects_df))[1]
            # Choose a meaningful numeric column: prefer adjusted mean / mean AUC over
            # count/DF/Reference columns
            preferred <- c("Adjusted_Mean", "Mean_AUC", "emmean", "estimate", "mean")
            effect_col <- intersect(preferred, colnames(effects_df))[1]
            if (is.na(effect_col)) {
              # fallback: first numeric that isn't N / DF / Reference
              skip_cols <- c("N", "DF", "Reference")
              num_cols <- colnames(effects_df)[sapply(effects_df, is.numeric)]
              num_cols <- setdiff(num_cols, c(skip_cols, label_col))
              effect_col <- if (length(num_cols) > 0) num_cols[1] else NA
            }
            if (!is.na(label_col) && !is.na(effect_col)) {
              # Detect CI or SE columns for error bars
              lower_col <- intersect(c("lower.CL", "lower_CL", "Lower", "lower"), colnames(effects_df))[1]
              upper_col <- intersect(c("upper.CL", "upper_CL", "Upper", "upper"), colnames(effects_df))[1]
              se_col    <- intersect(c("SE", "se", "SEM", "sem"), colnames(effects_df))[1]
              has_ci    <- !is.na(lower_col) && !is.na(upper_col)
              has_se    <- !is.na(se_col)

              plot_df <- data.frame(
                Treatment = as.character(effects_df[[label_col]]),
                Effect    = as.numeric(effects_df[[effect_col]]),
                stringsAsFactors = FALSE
              )
              if (has_ci) {
                plot_df$Lower <- as.numeric(effects_df[[lower_col]])
                plot_df$Upper <- as.numeric(effects_df[[upper_col]])
              } else if (has_se) {
                se_vals <- as.numeric(effects_df[[se_col]])
                plot_df$Lower <- plot_df$Effect - 1.96 * se_vals
                plot_df$Upper <- plot_df$Effect + 1.96 * se_vals
              }

              # Build beeswarm using per-subject growth rates when available,
              # overlaid with the LME adjusted mean ± CI from treatment_effects.
              gr_df <- tryCatch(as.data.frame(result$growth_rates), error = function(e) NULL)
              grp_col_gr <- if (!is.null(gr_df)) intersect(c("Group","Treatment","TreatmentGroup"), colnames(gr_df))[1] else NA
              val_col_gr <- if (!is.null(gr_df)) intersect(c("growth_rate","GrowthRate"), colnames(gr_df))[1] else NA

              if (!is.null(gr_df) && !is.na(grp_col_gr) && !is.na(val_col_gr)) {
                gr_df$plot_group <- as.character(gr_df[[grp_col_gr]])
                gr_df$plot_val   <- as.numeric(gr_df[[val_col_gr]])
                grp_levels <- as.character(plot_df$Treatment)
                gr_df$plot_group <- factor(gr_df$plot_group, levels = grp_levels)

                ci_half_fn <- function(x) {
                  n <- sum(!is.na(x))
                  qt(0.975, df = max(n - 1L, 1L)) * sd(x, na.rm = TRUE) / sqrt(n)
                }

                p_fallback <- ggplot2::ggplot(gr_df,
                    ggplot2::aes(x = plot_group, y = plot_val, colour = plot_group)) +
                  # Layer 1: 95% CI error bars with horizontal end caps
                  ggplot2::stat_summary(
                    fun.min = function(x) mean(x, na.rm = TRUE) - ci_half_fn(x),
                    fun.max = function(x) mean(x, na.rm = TRUE) + ci_half_fn(x),
                    geom = "errorbar", colour = "black", width = 0.25, linewidth = 0.9
                  ) +
                  # Layer 2: Mean marker as horizontal bar (degenerate crossbar)
                  ggplot2::stat_summary(
                    fun = function(x) mean(x, na.rm = TRUE),
                    fun.min = function(x) mean(x, na.rm = TRUE),
                    fun.max = function(x) mean(x, na.rm = TRUE),
                    geom = "crossbar", colour = "black", width = 0.25, linewidth = 0.9
                  ) +
                  # Layer 3: Individual points with horizontal jitter only
                  ggplot2::geom_point(
                    position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 42),
                    size = 2.5, alpha = 0.85
                  ) +
                  ggplot2::labs(title = "Growth Rate Comparison",
                               x = "Treatment", y = "Growth Rate (log-volume / day)") +
                  ggplot2::theme_classic() +
                  ggplot2::theme(legend.position = "none")
              } else {
                # No individual data — point + errorbar (still not a bar chart)
                p_fallback <- ggplot2::ggplot(plot_df,
                                              ggplot2::aes(x = reorder(Treatment, Effect),
                                                           y = Effect, colour = Treatment)) +
                  ggplot2::geom_point(size = 4) +
                  ggplot2::labs(title = "Treatment Effects (Adjusted Means)",
                               x = "Treatment", y = effect_col) +
                  ggplot2::theme_classic() +
                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                                legend.position = "none")

                if ("Lower" %in% colnames(plot_df)) {
                  p_fallback <- p_fallback +
                    ggplot2::geom_errorbar(
                      ggplot2::aes(ymin = Lower, ymax = Upper),
                      width = 0.2, colour = "black", linewidth = 0.8
                    )
                }
              }
              result$auc_plot <- p_fallback
            }
          }
        }
        
        # Plot 4: QQ-plot + residuals-vs-fitted for diagnostics tab
        if (!is.null(result$diagnostics) && !is.null(result$diagnostics$residuals)) {
          resid_vec <- result$diagnostics$residuals$residuals
          # QQ plot
          qq_data <- qqnorm(resid_vec, plot.it = FALSE)
          qq_df <- data.frame(theoretical = qq_data$x, sample = qq_data$y)
          qq_p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
            ggplot2::geom_point(alpha = 0.6) +
            ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
            ggplot2::theme_classic() +
            ggplot2::labs(title = "Q-Q Plot of Residuals",
                         x = "Theoretical Quantiles", y = "Sample Quantiles")
          result$caterpillar_plot <- qq_p  # keep for Plots tab
          result$diag_qq_plot     <- qq_p
        } else if (!is.null(result$model)) {
          tryCatch({
            resid_vec <- stats::residuals(result$model)
            qq_data   <- qqnorm(resid_vec, plot.it = FALSE)
            qq_df     <- data.frame(theoretical = qq_data$x, sample = qq_data$y)
            qq_p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = sample)) +
              ggplot2::geom_point(alpha = 0.6) +
              ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red") +
              ggplot2::theme_classic() +
              ggplot2::labs(title = "Q-Q Plot of Residuals",
                           x = "Theoretical Quantiles", y = "Sample Quantiles")
            result$caterpillar_plot <- qq_p
            result$diag_qq_plot     <- qq_p
          }, error = function(e) cat("WARNING: QQ plot from model failed:", e$message, "\n"))
        }
        # Residuals vs Fitted
        if (!is.null(result$model)) {
          tryCatch({
            fit_vals  <- stats::fitted(result$model)
            resid_vals <- stats::residuals(result$model)
            rf_df <- data.frame(Fitted = fit_vals, Residuals = resid_vals)
            result$diag_resid_fitted_plot <- ggplot2::ggplot(
              rf_df, ggplot2::aes(x = Fitted, y = Residuals)
            ) +
              ggplot2::geom_point(alpha = 0.5) +
              ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
              ggplot2::geom_smooth(method = "loess", se = FALSE,
                                  colour = "steelblue", linewidth = 0.8) +
              ggplot2::theme_classic() +
              ggplot2::labs(title = "Residuals vs Fitted",
                           x = "Fitted Values", y = "Residuals")
          }, error = function(e) cat("WARNING: resid-vs-fitted plot failed:", e$message, "\n"))
        }
        
      }, error = function(e) {
        cat("WARNING: Error generating plots:", e$message, "\n")
      })
    }
    
    # Generate all comparison types for comparisons tab
    if (!isTRUE(result$error) && !is.null(result$model)) {
      tryCatch({
        # Extract fixed effects terms and treatment comparisons from emmeans
        treatment_col_name <- col_names$treatment_column
        
        # Get emmeans object for treatment
        emm <- emmeans::emmeans(result$model, as.formula(paste("~", treatment_col_name)))
        emm_levels <- levels(emm)[[1]]   # character vector of factor levels
        
        # Type 1: All pairwise comparisons with Tukey adjustment
        # infer=c(TRUE,TRUE) returns both CIs and test statistics (t.ratio, p.value)
        result$all_pairwise_comparisons <- as.data.frame(
          summary(pairs(emm, adjust = "tukey"), infer = c(TRUE, TRUE))
        )
        
        # Type 2: vs Reference Group with Dunnett correction
        ref_group <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference))
          input$tg_reference else emm_levels[1]
        ref_idx <- which(emm_levels == ref_group)
        if (length(ref_idx) == 1) {
          result$reference_comparisons_dunnett <- as.data.frame(
            summary(emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_idx, adjust = "dunnett"),
                    infer = c(TRUE, TRUE))
          )
        }

        # Plot 5: Effect size forest plot — built here so reference_comparisons_dunnett exists
        pc <- if (!is.null(result$reference_comparisons_dunnett)) result$reference_comparisons_dunnett
              else if (!is.null(result$all_pairwise_comparisons)) result$all_pairwise_comparisons
              else NULL
        if (!is.null(pc)) {
          tryCatch({
            pc_df <- as.data.frame(pc)
            label_col   <- intersect(c("contrast", "Contrast", "comparison", "pair"), colnames(pc_df))[1]
            est_col     <- intersect(c("estimate", "Estimate"), colnames(pc_df))[1]
            se_col      <- intersect(c("SE"), colnames(pc_df))[1]
            lower_col   <- intersect(c("lower.CL", "Lower", "lower_CL", "lower"), colnames(pc_df))[1]
            upper_col   <- intersect(c("upper.CL", "Upper", "upper_CL", "upper"), colnames(pc_df))[1]
            p_col       <- intersect(c("p.value", "p_value", "p.adj", "adj.p.value", "P_Value"), colnames(pc_df))[1]
            if (!is.na(label_col) && !is.na(est_col)) {
              if (is.na(lower_col) && !is.na(se_col)) {
                pc_df$lower.CL <- pc_df[[est_col]] - 1.96 * pc_df[[se_col]]
                pc_df$upper.CL <- pc_df[[est_col]] + 1.96 * pc_df[[se_col]]
                lower_col <- "lower.CL"; upper_col <- "upper.CL"
              }
              p_label <- if (!is.na(p_col)) vapply(pc_df[[p_col]], format_pval, character(1)) else NULL
              using_dunnett <- !is.null(result$reference_comparisons_dunnett) &&
                               identical(pc, result$reference_comparisons_dunnett)
              ref_grp_ef <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference))
                              input$tg_reference else NULL
              contrast_labels <- pc_df[[label_col]]
              if (isTRUE(using_dunnett) && !is.null(ref_grp_ef)) {
                suffix_pat <- paste0("\\s*[-\u2013]\\s*", gsub("([+()])", "\\\\\\1", ref_grp_ef), "$")
                clean_labels <- sub(suffix_pat, "", contrast_labels)
                if (all(clean_labels != contrast_labels)) contrast_labels <- clean_labels
              }
              forest_df <- data.frame(
                Contrast = contrast_labels,
                Estimate = pc_df[[est_col]],
                Lower    = if (!is.na(lower_col)) pc_df[[lower_col]] else pc_df[[est_col]],
                Upper    = if (!is.na(upper_col)) pc_df[[upper_col]] else pc_df[[est_col]],
                P_Label  = if (is.null(p_label)) "" else paste0("p=", p_label),
                stringsAsFactors = FALSE
              )
              if (!is.null(ref_grp_ef)) {
                forest_df <- rbind(forest_df, data.frame(
                  Contrast = paste0(ref_grp_ef, "\n(Reference)"),
                  Estimate = 0, Lower = NA_real_, Upper = NA_real_,
                  P_Label  = "", stringsAsFactors = FALSE
                ))
              }
              forest_df$is_ref <- seq_len(nrow(forest_df)) == nrow(forest_df) & !is.null(ref_grp_ef)
              all_vals <- c(forest_df$Estimate[!forest_df$is_ref],
                            forest_df$Lower[!forest_df$is_ref],
                            forest_df$Upper[!forest_df$is_ref])
              max_abs  <- ceiling(max(abs(all_vals), na.rm = TRUE))
              x_lim    <- c(-max_abs, max_abs * 1.20)
              x_breaks <- pretty(c(-max_abs, max_abs), n = 6)
              x_breaks <- x_breaks[x_breaks >= -max_abs & x_breaks <= max_abs]
              result$effect_forest_plot <- ggplot2::ggplot(
                forest_df,
                ggplot2::aes(x = Estimate, y = stats::reorder(Contrast, -Estimate))
              ) +
                ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.6) +
                ggplot2::geom_errorbarh(
                  data = ~ subset(., !is_ref),
                  ggplot2::aes(xmin = Lower, xmax = Upper),
                  height = 0.25, colour = "steelblue"
                ) +
                ggplot2::geom_point(data = ~ subset(., !is_ref), size = 3, colour = "steelblue") +
                ggplot2::geom_point(data = ~ subset(., is_ref),  size = 4, colour = "grey50", shape = 18) +
                ggplot2::geom_text(
                  data = ~ subset(., !is_ref & nzchar(P_Label) & !is.na(Upper)),
                  ggplot2::aes(label = P_Label, x = Upper),
                  hjust = -0.2, size = 3.2
                ) +
                ggplot2::scale_x_continuous(limits = x_lim, breaks = x_breaks) +
                ggplot2::labs(
                  title = "Effect Size Forest Plot",
                  x = if (isTRUE(using_dunnett)) "Mean Difference (treatment \u2212 reference)" else "Mean Difference",
                  y = NULL
                ) +
                ggplot2::theme_classic() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
            }
          }, error = function(e) cat("WARNING: effect forest plot failed:", e$message, "\n"))
        }

      }, error = function(e) {
        cat("WARNING: Error generating comparison types:", e$message, "\n")
      })
    }

    # Persist the settings used so renderDTs can condition on them correctly
    # (the backend result$methods list is not reliably populated for all paths)
    result$transform_used  <- input$tg_transform   %||% "none"
    result$model_type_used <- input$tg_model_type  %||% "lme4"

    rv$tg_results <- result
    rv$tg_running <- FALSE

    # Store the data and column names needed for reactive plot re-rendering
    rv$tg_plot_data   <- data
    rv$tg_col_names   <- col_names
    rv$tg_plot_groups <- if (!is.null(rv$dose_col) && rv$dose_col %in% colnames(data) &&
                              !is.null(input$col_treatment) && input$col_treatment %in% colnames(data)) {
      # plot_tumor_growth builds Group as paste(treatment, dose, sep = " - Dose: ")
      # so colour-picker names must match that format exactly
      unique(paste(data[[input$col_treatment]], data[[rv$dose_col]], sep = " - Dose: "))
    } else if (!is.null(data) && !is.null(col_names$treatment_column) &&
                col_names$treatment_column %in% colnames(data)) {
      unique(as.character(data[[col_names$treatment_column]]))
    } else NULL
    
    cat("DEBUG: Analysis complete. Plots and comparisons generated.\n")
    
    # Update status
    if (isTRUE(result$error)) {
      output$tg_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error:", result$message)
      })
    } else {
      output$tg_status_msg <- renderUI({
        div(class = "alert alert-success",
            icon("check"), "Analysis completed successfully!")
      })
    }
  })
  
  # Summary output
  output$tg_summary <- renderText({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) {
      return(paste("Error:", rv$tg_results$message))
    }
    
    # Build a summary text from the summary object
    summary_list <- rv$tg_results$summary
    
    # Format summary information
    summary_lines <- c(
      "=== TUMOR GROWTH ANALYSIS SUMMARY ===",
      "",
      "MODEL SPECIFICATION:",
      paste("  Model Type:", summary_list$analysis_type),
      paste("  Transformation:", paste(summary_list$methods$volume_transformation)),
      "",
      "DATA DESCRIPTION:",
      paste("  Subjects:", summary_list$data_description$subjects),
      paste("  Treatment Groups:", summary_list$data_description$treatment_groups),
      paste("  Time Points:", summary_list$data_description$time_points),
      paste("  Reference Group:", summary_list$data_description$reference_group),
      "",
      "MODEL DETAILS:",
      paste("  Fixed Effects:", summary_list$model_specification$fixed_effects),
      paste("  Random Effects:", summary_list$model_specification$random_effects),
      paste("  Polynomial Degree:", summary_list$model_specification$polynomial_degree %||% 1L),
      "",
      "NOTES:",
      paste("  ", paste(summary_list$notes, collapse="\n  "))
    )
    
    paste(summary_lines, collapse = "\n")
  })

  # Update plot selector choices based on model type (AUC plot only in AUC mode)
  observeEvent(input$tg_model_type, {
    base_choices <- c(
      "Tumor Growth Curves"       = "growth_curves",
      "Growth Rates Plot"         = "growth_rate",
      "Waterfall Plot"            = "waterfall",
      "Effect Size Forest Plot"   = "effect_forest"
    )
    if (isTRUE(input$tg_model_type == "auc")) {
      base_choices <- append(
        base_choices,
        c("Area Under the Curve Plot" = "auc"),
        after = 2L  # insert after "Growth Rates Plot"
      )
    }
    current <- isolate(input$tg_plot_selector) %||% "growth_curves"
    selected <- if (current %in% base_choices) current else "growth_curves"
    updateSelectInput(session, "tg_plot_selector",
      choices  = base_choices,
      selected = selected
    )
  }, ignoreNULL = FALSE)

  # Conditionally render AUC tabs only in AUC mode
  output$tg_auc_tabs <- renderUI({
    if (isTRUE(input$tg_model_type == "auc")) {
      list(
        tabPanel("AUC Individual Data",
          tab_info(
            strong("Per-animal area under the growth curve (AUC) "), "computed for each subject across their observed measurement timepoints. Only displayed in AUC analysis mode.",
            tags$ul(
              tags$li(strong("ID"), " — individual animal identifier."),
              tags$li(strong("Treatment"), " — treatment group."),
              tags$li(strong("Cage"), " — cage identifier."),
              tags$li(strong("Group"), " — group assignment."),
              tags$li(strong("AUC"), " — area under the volume-time curve (mm³ × days)."),
              tags$li(strong("Last_Day"), " — final measurement timepoint for this animal."),
              tags$li(strong("First_Day"), " — initial measurement timepoint for this animal."),
              tags$li(strong("Extrapolated"), " — logical flag indicating if AUC was extended beyond observed data."),
              tags$li(strong("NumPoints"), " — number of observed volume measurements used to compute AUC.")
            )
          ),
          shinycssloaders::withSpinner(
            DTOutput("tg_auc_individual_table")
          ),
          downloadButton("download_tg_auc_individual", "Download AUC Individual Data")
        ),
        tabPanel("AUC Summary",
          tab_info(
            strong("Group-level AUC summary statistics"), " showing mean area under the growth curve, standard deviation, sample size, and standard error for each treatment group. Only displayed in AUC analysis mode.",
            tags$ul(
              tags$li(strong("Treatment"), " — treatment group name."),
              tags$li(strong("Mean_AUC"), " — group mean area under curve (mm³ × days)."),
              tags$li(strong("SD_AUC"), " — standard deviation of AUC within the group."),
              tags$li(strong("N"), " — number of animals in the group with valid AUC data."),
              tags$li(strong("SEM_AUC"), " — standard error of the mean AUC (SD / √N).")
            ),
            em("Use these summaries for descriptive group comparisons and to verify sample sizes and variability before statistical testing.")
          ),
          shinycssloaders::withSpinner(
            DTOutput("tg_auc_summary_table")
          ),
          downloadButton("download_tg_auc_summary", "Download AUC Summary")
        )
      )
    }
  })

  # ANOVA table
  output$tg_download_plot_btn <- renderUI({
    btn_id <- switch(input$tg_plot_selector %||% "growth_curves",
      growth_curves  = "download_tg_plot_growth",
      growth_rate    = "download_tg_plot_rates",
      auc            = "download_tg_plot_auc",
      waterfall      = "download_tg_plot_waterfall",
      effect_forest  = "download_tg_plot_effect_forest",
      "download_tg_plot_growth"
    )
    downloadButton(btn_id, "Download")
  })

  output$tg_anova_tab_info <- renderUI({
    model_type <- input$tg_model_type %||% "lme4"
    if (model_type == "auc") {
      tab_info(
        strong("One-way ANOVA "), "on per-animal AUC (area under the growth curve) values. Tests whether mean AUC differs significantly across treatment groups.",
        tags$ul(
          tags$li(strong("Df"), " — degrees of freedom (Treatment row = number of groups minus 1; Residuals row = total animals minus number of groups)."),
          tags$li(strong("Sum Sq"), " — sum of squares: variance attributable to that source."),
          tags$li(strong("Mean Sq"), " — Sum Sq ÷ Df."),
          tags$li(strong("F value"), " — ratio of the Treatment Mean Sq to the Residuals Mean Sq. Larger values indicate more separation between group means relative to within-group variability."),
          tags$li(strong("Pr(>F)"), " — p-value. Values below 0.05 indicate a statistically significant difference in mean AUC across treatment groups.")
        ),
        "This test assumes independent observations and equal variances across groups. For pairwise differences see the Comparisons tab."
      )
    } else {
      tab_info(
        strong("Type III ANOVA "), "from the linear mixed-effects model. Tests whether each fixed effect (treatment group, time, and their interaction) significantly explains variance in tumor volume after accounting for all other terms.",
        tags$ul(
          tags$li(strong("Chisq"), " — likelihood-ratio chi-square statistic for the effect. Larger values indicate stronger evidence against the null hypothesis."),
          tags$li(strong("Df"), " — degrees of freedom for the test (number of parameters being tested)."),
          tags$li(strong("Pr(>Chisq)"), " — p-value. Values below 0.05 indicate a statistically significant effect at the 5% level.")
        ),
        "The ", em("Treatment:Time"), " interaction row is the most biologically informative: a significant interaction means treatment groups diverge in their growth trajectories over time, not just in overall volume."
      )
    }
  })

  output$tg_anova_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$anova)) {
      return(datatable(data.frame(Message = "ANOVA table not available")))
    }
    
    anova_df <- as.data.frame(rv$tg_results$anova)
    anova_df <- format_df_for_display(anova_df)
    
    datatable(anova_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Comparisons title (reactive based on selection)
  output$tg_comparisons_title <- renderUI({
    comp_type <- input$tg_comparison_type
    
    title_text <- switch(comp_type,
      "reference" = h4("Pairwise Comparisons vs Reference Group"),
      "pairwise" = h4("All Pairwise Comparisons (Unadjusted)"),
      "tukey" = h4("All Pairwise Comparisons (Tukey HSD Adjusted)"),
      h4("Pairwise Comparisons")
    )
    
    tagList(
      title_text,
      p(class = "text-muted", style = "margin-bottom: 15px;",
        em("Change comparison type in the 'Post-hoc Comparisons' dropdown in Analysis Parameters"))
    )
  })
  
  # Comparisons table
  output$tg_comparisons_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) {
      return(datatable(data.frame(Message = "Comparisons not available")))
    }
    
    # Select comparison type based on user input
    comp_type <- input$tg_comparison_type
    
    comp_df <- NULL
    table_title <- NULL
    
    if (comp_type == "reference") {
      # vs Reference with Dunnett correction (preferred), fallback to unadjusted
      if (!is.null(rv$tg_results$reference_comparisons_dunnett)) {
        comp_df <- as.data.frame(rv$tg_results$reference_comparisons_dunnett)
        table_title <- "Comparisons: vs Reference Group (Dunnett-Adjusted)"
      } else if (!is.null(rv$tg_results$pairwise_comparisons)) {
        comp_df <- as.data.frame(rv$tg_results$pairwise_comparisons)
        table_title <- "Comparisons: vs Reference Group"
      }
    } else if (comp_type == "tukey") {
      # All pairwise with Tukey HSD adjustment
      if (!is.null(rv$tg_results$all_pairwise_comparisons)) {
        comp_df <- as.data.frame(rv$tg_results$all_pairwise_comparisons)
        table_title <- "Comparisons: All Pairwise (Tukey HSD Adjusted)"
      }
    }
    
    if (is.null(comp_df)) {
      return(datatable(data.frame(Message = paste("Comparisons not available for", comp_type))))
    }

    # Add fold-change back-transform for LME4 log-scale models
    transform_used   <- rv$tg_results$transform_used  %||% "none"
    model_type_used  <- rv$tg_results$model_type_used %||% "lme4"
    if (model_type_used == "lme4" && transform_used == "log") {
      est_col <- intersect(c("estimate", "Estimate"), colnames(comp_df))[1]
      lo_col  <- intersect(c("lower.CL", "lower_CL", "Lower"), colnames(comp_df))[1]
      hi_col  <- intersect(c("upper.CL", "upper_CL", "Upper"), colnames(comp_df))[1]
      if (!is.na(est_col)) {
        comp_df$Fold_Change <- round(exp(as.numeric(comp_df[[est_col]])), 3)
        if (!is.na(lo_col)) comp_df$Fold_Change_Lower <- round(exp(as.numeric(comp_df[[lo_col]])), 3)
        if (!is.na(hi_col)) comp_df$Fold_Change_Upper <- round(exp(as.numeric(comp_df[[hi_col]])), 3)
      }
    }

    # Enforce column order: contrast, df, estimate, lower.CL, upper.CL, SE,
    # t.ratio, p.value, Fold_Change, Fold_Change_Lower, Fold_Change_Upper, then any remaining cols
    preferred_order <- c("contrast", "Contrast",
                         "df",
                         "estimate", "Estimate",
                         "lower.CL", "lower_CL",
                         "upper.CL", "upper_CL",
                         "SE",
                         "t.ratio", "z.ratio",
                         "p.value",
                         "Fold_Change", "Fold_Change_Lower", "Fold_Change_Upper")
    ordered_cols <- intersect(preferred_order, colnames(comp_df))
    rest_cols    <- setdiff(colnames(comp_df), ordered_cols)
    comp_df      <- comp_df[, c(ordered_cols, rest_cols), drop = FALSE]

    comp_df <- format_df_for_display(comp_df)
    datatable(comp_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Treatment effects table
  output$tg_effects_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$treatment_effects)) {
      return(datatable(data.frame(Message = "Treatment effects not available")))
    }

    effects_df <- as.data.frame(rv$tg_results$treatment_effects)
    transform_used   <- rv$tg_results$transform_used  %||% "none"
    model_type_used  <- rv$tg_results$model_type_used %||% "lme4"

    if (model_type_used == "lme4") {
      # Keep only the columns the user needs; drop any extras / _orig columns
      keep <- c("Group", "DF", "Adjusted_Mean", "Lower_CL", "Upper_CL", "SE", "Note")
      effects_df <- effects_df[, intersect(keep, colnames(effects_df)), drop = FALSE]

      # Enforce display order
      effects_df <- effects_df[, intersect(keep, colnames(effects_df)), drop = FALSE]

      effects_df <- format_df_for_display(effects_df)

      scale_note <- switch(transform_used, "log" = " (log scale)", "sqrt" = " (sqrt scale)", "")
      colnames(effects_df) <- gsub("^Adjusted_Mean$", paste0("Adj Mean", scale_note), colnames(effects_df))
      colnames(effects_df) <- gsub("^Lower_CL$",      paste0("Lower CI", scale_note), colnames(effects_df))
      colnames(effects_df) <- gsub("^Upper_CL$",      paste0("Upper CI", scale_note), colnames(effects_df))
    } else {
      # Non-lme4 models: show all available columns, formatted
      effects_df <- format_df_for_display(effects_df)
    }

    datatable(effects_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Growth rates table
  output$tg_growth_rates_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$growth_rates)) {
      return(datatable(data.frame(Message = "Growth rates not available")))
    }
    
    rates_df <- as.data.frame(rv$tg_results$growth_rates)
    if ("growth_rate" %in% colnames(rates_df)) {
      # Back-transformed point estimates
      rates_df$daily_fold_change  <- exp(rates_df$growth_rate)
      rates_df$pct_change_per_day <- (exp(rates_df$growth_rate) - 1) * 100

      # Back-transformed 95% CIs (from per-animal lm confint)
      if ("growth_rate_lower" %in% colnames(rates_df)) {
        rates_df$fold_change_lower    <- exp(rates_df$growth_rate_lower)
        rates_df$fold_change_upper    <- exp(rates_df$growth_rate_upper)
        rates_df$pct_change_lower     <- (exp(rates_df$growth_rate_lower) - 1) * 100
        rates_df$pct_change_upper     <- (exp(rates_df$growth_rate_upper) - 1) * 100
      }

      # Enforce column order:
      # growth_rate, growth_rate_lower, growth_rate_upper,
      # daily_fold_change, fold_change_lower, fold_change_upper,
      # pct_change_per_day, pct_change_lower, pct_change_upper,
      # R_squared, then ID/Treatment/Cage
      preferred <- c("growth_rate", "growth_rate_lower", "growth_rate_upper",
                     "daily_fold_change", "fold_change_lower", "fold_change_upper",
                     "pct_change_per_day", "pct_change_lower", "pct_change_upper",
                     "R_squared")
      ordered_cols <- intersect(preferred, colnames(rates_df))
      rest_cols    <- setdiff(colnames(rates_df), ordered_cols)
      rates_df     <- rates_df[, c(rest_cols, ordered_cols), drop = FALSE]
    }
    rates_df <- format_df_for_display(rates_df)
    # Human-readable column headers
    rename_map <- c(
      growth_rate        = "Growth Rate",
      growth_rate_lower  = "GR 95% CI Lower Bound",
      growth_rate_upper  = "GR 95% CI Upper Bound",
      daily_fold_change  = "Fold Change",
      fold_change_lower  = "FC 95% CI Lower Bound",
      fold_change_upper  = "FC 95% CI Upper Bound",
      pct_change_per_day = "Percent Change",
      pct_change_lower   = "PC 95% CI Lower Bound",
      pct_change_upper   = "PC 95% CI Upper Bound",
      R_squared          = "R²"
    )
    new_names <- colnames(rates_df)
    for (old in names(rename_map)) new_names[new_names == old] <- rename_map[[old]]
    datatable(rates_df,
      colnames = new_names,
      options  = list(scrollX = TRUE, paging = FALSE, searching = FALSE, scrollY = "500px", scrollCollapse = TRUE)
    )
  })
  
  # AUC individual data table
  output$tg_auc_individual_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$auc_analysis$individual)) {
      return(datatable(data.frame(Message = "AUC individual data not available")))
    }
    
    auc_ind_df <- as.data.frame(rv$tg_results$auc_analysis$individual)
    auc_ind_df <- format_df_for_display(auc_ind_df)
    datatable(auc_ind_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # AUC summary table
  output$tg_auc_summary_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$auc_analysis$summary)) {
      return(datatable(data.frame(Message = "AUC summary not available")))
    }
    
    auc_sum_df <- as.data.frame(rv$tg_results$auc_analysis$summary)
    auc_sum_df <- format_df_for_display(auc_sum_df)
    datatable(auc_sum_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Data summary table
  output$tg_data_summary_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$data_summary)) {
      return(datatable(data.frame(Message = "Data summary not available")))
    }
    
    data_sum_df <- as.data.frame(rv$tg_results$data_summary)
    data_sum_df <- format_df_for_display(data_sum_df)
    datatable(data_sum_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Group order and color combined control for Plots tab
  output$tg_group_order_color_ui <- renderUI({
    req(rv$tg_plot_groups)
    group_order_and_color_ui("tg", rv$tg_plot_groups)
  })

  # Growth curves plot  — reactive: re-renders on colour/option changes
  tg_growth_plot <- reactive({
    tryCatch({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    req(rv$tg_plot_data, rv$tg_col_names)
    data      <- rv$tg_plot_data
    col_names <- rv$tg_col_names
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    colors    <- build_color_vector(input, "tg_color", ordered_grps)
    # For dose data, apply factor level order on the original treatment column;
    # plot_tumor_growth will build Group = 'Treatment - Dose: X' internally.
    # For non-dose data, apply order on the treatment column directly.
    if (isTRUE(length(ordered_grps) > 0)) {
      if (!is.null(rv$dose_col) && rv$dose_col %in% colnames(data) &&
          !is.null(input$col_treatment) && input$col_treatment %in% colnames(data)) {
        # Derive treatment-level order from the combined group order
        tx_order <- unique(sub(" - Dose: .*$", "", ordered_grps))
        data[[input$col_treatment]] <- factor(data[[input$col_treatment]], levels = tx_order)
      } else if (isTRUE(col_names$treatment_column %in% colnames(data))) {
        data[[col_names$treatment_column]] <- factor(data[[col_names$treatment_column]],
                                                     levels = ordered_grps)
      }
    }

    surv_col_mapped <- local({
      col <- input$col_survival_event
      is.character(col) && length(col) == 1L &&
        nzchar(col) && col != "None" && isTRUE(col %in% colnames(data))
    })
    extra_pts <- local({
      ep <- input$tg_extrapolation_points
      if (is.numeric(ep) && length(ep) == 1L && !is.na(ep)) as.integer(ep) else 0L
    })
    do_extrapolate <- isTRUE(surv_col_mapped) && extra_pts > 0L && isTRUE(input$tg_opt_show_extrap)

    p <- tryCatch({
      mouseExperiment::plot_tumor_growth(
        df                  = data,
        volume_column       = col_names$volume_column,
        day_column          = col_names$time_column,
        # When a dose column is present, pass the original treatment column so
        # plot_tumor_growth can build 'Group = Treatment - Dose: X' itself,
        # matching the names used in rv$tg_plot_groups and the colour pickers.
        treatment_column    = if (!is.null(rv$dose_col)) input$col_treatment else col_names$treatment_column,
        dose_column         = rv$dose_col,
        ID_column           = col_names$id_column,
        cage_column         = col_names$cage_column,
        survival_column     = if (isTRUE(surv_col_mapped)) input$col_survival_event else NULL,
        extrapolate_volumes = do_extrapolate,
        extrapolation_points = if (extra_pts > 0L) extra_pts else "all",
        group_summary_line  = isTRUE(input$tg_opt_group_line),
        point_size          = input$tg_opt_growth_pt_size %||% 2
      ) + ggplot2::theme_classic()
    }, error = function(e) {
      if (!is.null(rv$tg_results$tumor_growth_plot)) return(rv$tg_results$tumor_growth_plot)
      NULL
    })
    req(!is.null(p))
    shapes <- build_shape_vector(input, "tg_shape", ordered_grps)
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_color_manual(values = colors)
      # scale_fill_manual keeps dual-fill shapes (e.g. pch 25, down-triangle)
      # rendered with the correct group colour interior
      p <- p + ggplot2::scale_fill_manual(values = colors)
    }
    if (!is.null(shapes)) {
      p <- p + ggplot2::scale_shape_manual(values = shapes)
    }
    # Apply text size, title/subtitle, and legend position
    text_size  <- input$tg_opt_text_size %||% 16
    legend_pos <- input$tg_opt_legend_pos %||% "right"
    p <- apply_plot_labels(p, input$tg_plot_title %||% "", input$tg_plot_subtitle %||% "",
                           input$tg_plot_title_align %||% "center", text_size,
                           input$tg_plot_xaxis_label %||% "", input$tg_plot_yaxis_label %||% "",
                           input$tg_plot_legend_label %||% "")
    p <- apply_legend_pos(p, legend_pos)
    # Apply y-axis scale
    yscale <- input$tg_opt_yscale %||% "linear"
    if (yscale == "log10") p <- p + ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.05)))
    # Apply x-axis start and ensure 5% right padding.
    # coord_cartesian is always applied so the explicit x_upper value
    # (data max + 5%) is honoured regardless of scale expand settings.
    x_start <- local({
      v <- input$tg_opt_xstart
      if (is.numeric(v) && length(v) == 1L && !is.na(v) && v > 0) v else 0
    })
    x_upper <- tryCatch({
      time_col <- rv$tg_col_names$time_column
      x_max <- max(rv$tg_plot_data[[time_col]], na.rm = TRUE)
      x_max + 0.05 * (x_max - x_start)
    }, error = function(e) NA)
    if (!is.na(x_upper))
      p <- suppressMessages(p + ggplot2::coord_cartesian(xlim = c(x_start, x_upper)))
    # Annotate hollow circles = extrapolated when extrapolation is active
    if (isTRUE(do_extrapolate))
      p <- p + ggplot2::labs(caption = "Hollow circles / dashed lines = extrapolated data points after mouse death")
    p
    }, error = function(e) {
      msg <- conditionMessage(e)
      calls <- paste(sapply(sys.calls(), deparse), collapse = " -> ")
      cat("ERROR in tg_plot_growth:", msg, "\nCALLS:", calls, "\n")
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
          label = paste0("Plot error:\n", msg),
          hjust = 0.5, vjust = 0.5, size = 4, color = "red") +
        ggplot2::theme_void()
    })
  })

  output$tg_plot_growth <- renderPlot({
    tg_growth_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)
  
  tg_rates_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    req(rv$tg_results$growth_rates)
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    colors     <- build_color_vector(input, "tg_color", ordered_grps)
    shapes     <- build_shape_vector(input, "tg_shape", ordered_grps)
    center_bar <- input$tg_opt_center_bar %||% "mean"
    show_mean  <- center_bar == "mean"
    err_type   <- input$tg_opt_error_bar  %||% "SEM"
    pt_size    <- input$tg_opt_point_size %||% 2.5
    jit_width  <- input$tg_opt_jitter_width %||% 0.2
    text_size  <- input$tg_opt_text_size %||% 12
    legend_pos <- input$tg_opt_legend_pos %||% "auto"
    rates      <- rv$tg_results$growth_rates
    grp_col_r  <- intersect(c("Group", "Treatment"), colnames(rates))[1]
    # Derive treatment-level order for plot_growth_rate (strip dose suffix if present)
    gr_order <- if (length(ordered_grps) > 0) {
      raw <- unique(sub(" - Dose: .*$", "", ordered_grps))
      raw[raw %in% as.character(unique(rates[["Treatment"]]))]
    } else NULL
    # Also factor the rates data for any stat_summary overlays added after the main plot
    if (!is.null(gr_order) && length(gr_order) > 0 && "Treatment" %in% colnames(rates))
      rates$Treatment <- factor(rates$Treatment, levels = gr_order)
    p <- tryCatch({
      mouseExperiment::plot_growth_rate(
        growth_data    = rates,
        colors         = colors,
        group_order    = if (length(gr_order) > 0) gr_order else NULL,
        show_mean      = show_mean,
        error_bar_type = err_type,
        point_size     = pt_size,
        jitter_width   = jit_width
      ) + ggplot2::theme_classic()
    }, error = function(e) rv$tg_results$growth_rate_plot)
    req(!is.null(p))
    # Extra center bar overlays for median / geometric mean
    if (center_bar %in% c("median", "geometric_mean") && !is.na(grp_col_r)) {
      cfun <- if (center_bar == "median") median else
                function(x) exp(mean(log(pmax(x, 1e-9, na.rm = TRUE))))
      p <- p + ggplot2::stat_summary(
        data = rates,
        ggplot2::aes(x = .data[[grp_col_r]], y = growth_rate,
                     color = .data[[grp_col_r]], group = .data[[grp_col_r]]),
        fun.min = cfun, fun = cfun, fun.max = cfun,
        geom = "crossbar", width = 0.4, linewidth = 0.8
      )
    }
    if (!is.null(shapes)) p <- p + ggplot2::scale_shape_manual(values = shapes)
    p <- apply_plot_labels(p, input$tg_plot_title %||% "", input$tg_plot_subtitle %||% "",
                           input$tg_plot_title_align %||% "center", text_size,
                           input$tg_plot_xaxis_label %||% "", input$tg_plot_yaxis_label %||% "",
                           input$tg_plot_legend_label %||% "")
    p <- apply_legend_pos(p, legend_pos)
    p
  })

  output$tg_plot_rates <- renderPlot({
    tg_rates_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)
  
  # AUC/Treatment effects plot — reactive
  tg_auc_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    # Shared controls
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    colors     <- build_color_vector(input, "tg_color", ordered_grps)
    shapes     <- build_shape_vector(input, "tg_shape", ordered_grps)
    center_bar <- input$tg_opt_center_bar %||% "mean"
    err_type   <- input$tg_opt_error_bar  %||% "SEM"
    pt_size    <- input$tg_opt_point_size %||% 2.5
    jit_width  <- input$tg_opt_jitter_width %||% 0.1
    text_size  <- input$tg_opt_text_size %||% 16
    legend_pos <- input$tg_opt_legend_pos %||% "right"

    # Helper: build a summary df (per group: center, lower, upper)
    make_summary <- function(df, grp_col, val_col) {
      groups <- levels(factor(df[[grp_col]]))
      do.call(rbind, lapply(groups, function(g) {
        x <- df[[val_col]][as.character(df[[grp_col]]) == g]
        x <- x[!is.na(x)]
        n <- length(x)
        ctr <- switch(center_bar,
          mean           = mean(x),
          median         = median(x),
          geometric_mean = exp(mean(log(pmax(x, 1e-9)))),
                           mean(x))
        half <- switch(err_type,
          SEM = sd(x) / sqrt(n),
          SD  = sd(x),
          CI  = qt(0.975, df = max(n - 1L, 1L)) * sd(x) / sqrt(n),
                0)
        data.frame(grp = g, ctr = ctr, lo = ctr - half, hi = ctr + half,
                   stringsAsFactors = FALSE)
      }))
    }

    # Helper: y-axis limits — floor at 0 (unless data dips below), 10% top buffer
    y_limits <- function(point_vals, lo_vals = NULL, hi_vals = NULL) {
      all_vals <- c(point_vals, lo_vals, hi_vals)
      all_vals <- all_vals[!is.na(all_vals)]
      ymin <- min(0, min(all_vals))
      ymax <- max(all_vals) * 1.10
      c(ymin, ymax)
    }

    # Helper: apply colours + labels + legend
    finish_plot <- function(p) {
      if (!is.null(colors)) p <- p + ggplot2::scale_colour_manual(values = colors)
      if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
      if (!is.null(shapes)) p <- p + ggplot2::scale_shape_manual(values = shapes)
      p <- apply_plot_labels(p, input$tg_plot_title %||% "", input$tg_plot_subtitle %||% "",
                             input$tg_plot_title_align %||% "center", text_size = text_size,
                             xlab = input$tg_plot_xaxis_label %||% "",
                             ylab = input$tg_plot_yaxis_label %||% "",
                             legend_lab = input$tg_plot_legend_label %||% "")
      apply_legend_pos(p, legend_pos)
    }

    # ── Path A: AUC model ── per-animal AUC values available ──────────────
    if (!is.null(rv$tg_results$auc_analysis) &&
        !is.null(rv$tg_results$auc_analysis$individual)) {

      auc_ind   <- rv$tg_results$auc_analysis$individual
      grp_col_a <- intersect(c("Treatment", "Group"), colnames(auc_ind))[1]
      req(!is.na(grp_col_a))

      # Factor ordering
      if (length(ordered_grps) > 0) {
        auc_levels <- unique(sub(" - Dose: .*$", "", ordered_grps))
        auc_levels <- auc_levels[auc_levels %in% as.character(unique(auc_ind[[grp_col_a]]))]
        if (length(auc_levels) > 0)
          auc_ind[[grp_col_a]] <- factor(auc_ind[[grp_col_a]], levels = auc_levels)
      }

      summ <- make_summary(auc_ind, grp_col_a, "AUC")
      summ$grp <- factor(summ$grp, levels = levels(factor(auc_ind[[grp_col_a]])))

      p <- ggplot2::ggplot(auc_ind,
             ggplot2::aes(x = .data[[grp_col_a]], y = .data[["AUC"]],
                          colour = .data[[grp_col_a]], shape = .data[[grp_col_a]])) +
        ggplot2::theme_classic() +
        ggplot2::labs(x = "Treatment", y = "AUC", title = "Area Under the Curve") +
        ggplot2::theme(legend.title = ggplot2::element_blank())

      if (err_type != "none")
        p <- p + ggplot2::geom_errorbar(
          data = summ,
          ggplot2::aes(x = grp, ymin = lo, ymax = hi),
          inherit.aes = FALSE,
          colour = "black", width = 0.25, linewidth = 0.9)

      if (center_bar != "none")
        p <- p + ggplot2::geom_crossbar(
          data = summ,
          ggplot2::aes(x = grp, y = ctr, ymin = ctr, ymax = ctr),
          inherit.aes = FALSE,
          colour = "black", width = 0.25, linewidth = 0.9)

      p <- p + ggplot2::geom_point(
        position = ggplot2::position_jitter(width = jit_width, height = 0, seed = 42),
        size = pt_size, alpha = 0.85)

      ylims <- y_limits(auc_ind[["AUC"]], summ$lo, summ$hi)
      p <- p + ggplot2::coord_cartesian(ylim = ylims)
      return(finish_plot(p))
    }

    # ── Path B: LME4 model ── same visual style as Path A, using growth rates ──
    rates <- rv$tg_results$growth_rates
    if (!is.null(rates) && nrow(rates) > 0) {
      grp_col_b <- intersect(c("Treatment", "Group"), colnames(rates))[1]
      val_col_b <- intersect(c("growth_rate", "GrowthRate"), colnames(rates))[1]
      if (!is.na(grp_col_b) && !is.na(val_col_b)) {
        gr_order <- if (length(ordered_grps) > 0) {
          raw <- unique(sub(" - Dose: .*$", "", ordered_grps))
          raw[raw %in% as.character(unique(rates[[grp_col_b]]))]
        } else NULL
        if (!is.null(gr_order) && length(gr_order) > 0)
          rates[[grp_col_b]] <- factor(rates[[grp_col_b]], levels = gr_order)

        summ_b <- make_summary(rates, grp_col_b, val_col_b)
        summ_b$grp <- factor(summ_b$grp, levels = levels(factor(rates[[grp_col_b]])))

        p <- ggplot2::ggplot(rates,
               ggplot2::aes(x = .data[[grp_col_b]], y = .data[[val_col_b]],
                            colour = .data[[grp_col_b]], shape = .data[[grp_col_b]])) +
          ggplot2::theme_classic() +
          ggplot2::labs(x = "Treatment", y = "Growth Rate (log-volume / day)",
                        title = "Growth Rate Comparison") +
          ggplot2::theme(legend.title = ggplot2::element_blank())

        if (err_type != "none")
          p <- p + ggplot2::geom_errorbar(
            data = summ_b,
            ggplot2::aes(x = grp, ymin = lo, ymax = hi),
            inherit.aes = FALSE,
            colour = "black", width = 0.25, linewidth = 0.9)

        if (center_bar != "none")
          p <- p + ggplot2::geom_crossbar(
            data = summ_b,
            ggplot2::aes(x = grp, y = ctr, ymin = ctr, ymax = ctr),
            inherit.aes = FALSE,
            colour = "black", width = 0.25, linewidth = 0.9)

        p <- p + ggplot2::geom_point(
          position = ggplot2::position_jitter(width = jit_width, height = 0, seed = 42),
          size = pt_size, alpha = 0.85)

        ylims <- y_limits(rates[[val_col_b]], summ_b$lo, summ_b$hi)
        p <- p + ggplot2::coord_cartesian(ylim = ylims)
        return(finish_plot(p))
      }
    }

    # ── Path C: last resort — stored plot ─────────────────────────────────
    p <- rv$tg_results$auc_plot
    if (is.null(p)) return(NULL)
    finish_plot(p)
  })

  output$tg_plot_auc <- renderPlot({
    tg_auc_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)
  
  # ── Waterfall Plot ────────────────────────────────────────────────────────
  tg_wf_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    rates <- rv$tg_results$growth_rates
    req(!is.null(rates), nrow(rates) > 0)

    tryCatch({
      wf_data <- rates |>
        dplyr::filter(!is.na(.data$growth_rate)) |>
        dplyr::mutate(
          animal_label = if ("Cage" %in% colnames(rates))
                           paste0(.data$Cage, "/", .data$ID)
                         else
                           as.character(.data$ID),
          pct_per_day = (exp(.data$growth_rate) - 1) * 100
        ) |>
        dplyr::arrange(dplyr::desc(.data$pct_per_day)) |>
        dplyr::mutate(rank = dplyr::row_number())

      req(nrow(wf_data) > 0)

      trt_col      <- "Treatment"
      ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
      colors       <- build_color_vector(input, "tg_color", ordered_grps)
      text_size    <- input$tg_opt_text_size %||% 16

      p <- ggplot2::ggplot(
             wf_data,
             ggplot2::aes(x = .data$rank,
                          y = .data$pct_per_day,
                          fill = .data[[trt_col]])) +
        ggplot2::geom_col(width = 0.8, colour = "grey30", linewidth = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.6) +
        ggplot2::scale_fill_manual(values = colors, name = trt_col) +
        ggplot2::scale_x_continuous(
          breaks = wf_data$rank,
          labels = wf_data$animal_label,
          expand = c(0.01, 0)) +
        ggplot2::labs(
          x = "Animal",
          y = "Tumor growth rate (% change per day)") +
        ggplot2::theme_classic(base_size = text_size) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

      p <- apply_legend_pos(p, input$tg_opt_legend_pos %||% "right")
      apply_plot_labels(p,
        input$tg_plot_title    %||% "",
        input$tg_plot_subtitle %||% "",
        input$tg_plot_title_align %||% "center",
        text_size,
        input$tg_plot_xaxis_label %||% "",
        input$tg_plot_yaxis_label %||% "",
        input$tg_plot_legend_label %||% "")

    }, error = function(e) {
      cat("ERROR in tg_wf_plot:", conditionMessage(e), "\n")
      NULL
    })
  })

  output$tg_plot_waterfall <- renderPlot({
    p <- tg_wf_plot()
    req(!is.null(p))
    p
  }, width  = function() input$tg_plot_width  %||% 800,
     height = function() input$tg_plot_height %||% 550)

  # Diagnostics/QQ plot
  output$tg_plot_caterpillar <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    
    if (!is.null(rv$tg_results$caterpillar_plot)) {
      return(rv$tg_results$caterpillar_plot)
    } else {
      plot(1, 1, type="n", main="QQ plot generating...", xlab="", ylab="")
    }
  })

  # Effect size forest plot — rebuilt from raw comparison data each render
  # so group order changes take effect without re-running analysis
  tg_forest_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    # Pick best comparison source (Dunnett preferred)
    using_dunnett <- !is.null(rv$tg_results$reference_comparisons_dunnett)
    pc <- if (using_dunnett) rv$tg_results$reference_comparisons_dunnett
          else rv$tg_results$all_pairwise_comparisons
    if (is.null(pc)) return(NULL)

    text_size    <- input$tg_opt_text_size %||% 12
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    grps_clean   <- unique(sub(" - Dose: .*$", "", ordered_grps))
    colors       <- build_color_vector(input, "tg_color", ordered_grps)
    shapes       <- build_shape_vector(input, "tg_shape", ordered_grps)

    pc_df     <- as.data.frame(pc)
    label_col <- intersect(c("contrast", "Contrast", "comparison", "pair"), colnames(pc_df))[1]
    est_col   <- intersect(c("estimate", "Estimate"), colnames(pc_df))[1]
    se_col    <- intersect(c("SE"), colnames(pc_df))[1]
    lower_col <- intersect(c("lower.CL", "Lower", "lower_CL", "lower"), colnames(pc_df))[1]
    upper_col <- intersect(c("upper.CL", "Upper", "upper_CL", "upper"), colnames(pc_df))[1]
    p_col     <- intersect(c("p.value", "p_value", "p.adj", "adj.p.value", "P_Value"), colnames(pc_df))[1]
    if (is.na(label_col) || is.na(est_col)) return(NULL)

    if (is.na(lower_col) && !is.na(se_col)) {
      pc_df$lower.CL <- pc_df[[est_col]] - 1.96 * pc_df[[se_col]]
      pc_df$upper.CL <- pc_df[[est_col]] + 1.96 * pc_df[[se_col]]
      lower_col <- "lower.CL"; upper_col <- "upper.CL"
    }

    p_label <- if (!is.na(p_col)) vapply(pc_df[[p_col]], format_pval, character(1)) else NULL
    ref_grp <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference)) input$tg_reference else NULL

    contrast_labels <- as.character(pc_df[[label_col]])
    if (isTRUE(using_dunnett) && !is.null(ref_grp)) {
      suffix_pat <- paste0("\\s*[-\u2013]\\s*", gsub("([+()])", "\\\\\\1", ref_grp), "$")
      clean_tmp  <- sub(suffix_pat, "", contrast_labels)
      if (all(clean_tmp != contrast_labels)) contrast_labels <- clean_tmp
    }

    forest_df <- data.frame(
      Contrast = contrast_labels,
      Estimate = as.numeric(pc_df[[est_col]]),
      Lower    = if (!is.na(lower_col)) as.numeric(pc_df[[lower_col]]) else as.numeric(pc_df[[est_col]]),
      Upper    = if (!is.na(upper_col)) as.numeric(pc_df[[upper_col]]) else as.numeric(pc_df[[est_col]]),
      P_Label  = if (is.null(p_label)) "" else paste0("p=", p_label),
      is_ref   = FALSE,
      stringsAsFactors = FALSE
    )
    if (!is.null(ref_grp)) {
      forest_df <- rbind(forest_df, data.frame(
        Contrast = paste0(ref_grp, "\n(Reference)"),
        Estimate = 0, Lower = NA_real_, Upper = NA_real_,
        P_Label = "", is_ref = TRUE,
        stringsAsFactors = FALSE
      ))
    }

    # Sort non-reference rows by user group order; reference always at bottom of plot
    non_ref <- forest_df[!forest_df$is_ref, , drop = FALSE]
    ref_row <- forest_df[ forest_df$is_ref, , drop = FALSE]
    if (length(grps_clean) > 0) {
      rank_by_grp <- match(non_ref$Contrast, grps_clean)
      rank_by_grp[is.na(rank_by_grp)] <- .Machine$integer.max
      non_ref <- non_ref[order(rank_by_grp), , drop = FALSE]
    } else {
      non_ref <- non_ref[order(-non_ref$Estimate), , drop = FALSE]
    }
    # Factor levels bottom-to-top: reference at bottom, first group-order item at top
    forest_df <- rbind(ref_row, non_ref)
    forest_df$Contrast <- factor(
      forest_df$Contrast,
      levels = c(as.character(ref_row$Contrast), rev(as.character(non_ref$Contrast)))
    )

    show_pval   <- isTRUE(input$tg_opt_forest_show_pval   %||% TRUE)
    show_effect <- isTRUE(input$tg_opt_forest_show_effect %||% FALSE)
    show_ci     <- isTRUE(input$tg_opt_forest_show_ci     %||% FALSE)
    any_ann     <- show_pval || show_effect || show_ci

    # Build per-row annotation label from whichever elements are enabled
    forest_df$Ann_Label <- vapply(seq_len(nrow(forest_df)), function(j) {
      if (forest_df$is_ref[j]) return("")
      parts <- character(0)
      if (show_effect)
        parts <- c(parts, paste0("est: ", signif(forest_df$Estimate[j], 3)))
      if (show_ci && !is.na(forest_df$Lower[j]))
        parts <- c(parts,
          paste0("[", signif(forest_df$Lower[j], 3), ", ", signif(forest_df$Upper[j], 3), "]"))
      if (show_pval && nzchar(forest_df$P_Label[j]))
        parts <- c(parts, forest_df$P_Label[j])
      paste(parts, collapse = "\n")
    }, character(1L))

    # Build per-contrast colour/shape maps (contrast labels correspond to group names)
    contrast_lvls <- setdiff(as.character(levels(forest_df$Contrast)),
                             as.character(ref_row$Contrast))
    forest_colours <- if (!is.null(colors)) {
      vals <- ifelse(contrast_lvls %in% names(colors), unname(colors[contrast_lvls]), "steelblue")
      setNames(vals, contrast_lvls)
    } else setNames(rep("steelblue", length(contrast_lvls)), contrast_lvls)
    forest_shapes <- if (!is.null(shapes)) {
      vals <- ifelse(contrast_lvls %in% names(shapes), unname(shapes[contrast_lvls]), 16L)
      setNames(as.integer(vals), contrast_lvls)
    } else setNames(rep(16L, length(contrast_lvls)), contrast_lvls)

    all_vals <- c(non_ref$Estimate, non_ref$Lower, non_ref$Upper)
    max_abs  <- ceiling(max(abs(all_vals), na.rm = TRUE))
    # Widen right margin when annotations are shown so text is not clipped
    x_lim    <- c(-max_abs, max_abs * if (any_ann) 1.55 else 1.20)
    x_breaks <- pretty(c(-max_abs, max_abs), n = 6)
    x_breaks <- x_breaks[x_breaks >= -max_abs & x_breaks <= max_abs]

    x_lab <- if (isTRUE(using_dunnett)) "Mean Difference (treatment \u2212 reference)" else "Mean Difference"
    p <- ggplot2::ggplot(forest_df, ggplot2::aes(x = Estimate, y = Contrast)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.6) +
      ggplot2::geom_errorbarh(
        data = ~ subset(., !is_ref),
        ggplot2::aes(xmin = Lower, xmax = Upper, colour = Contrast),
        height = 0.25) +
      ggplot2::geom_point(data = ~ subset(., !is_ref),
        ggplot2::aes(colour = Contrast, shape = Contrast), size = 3) +
      ggplot2::geom_point(data = ~ subset(., is_ref),  size = 4, colour = "grey50", shape = 18) +
      ggplot2::geom_text(
        data = ~ subset(., !is_ref & nzchar(Ann_Label) & !is.na(Upper)),
        ggplot2::aes(label = Ann_Label, x = Upper),
        hjust = -0.15, vjust = 0.5, size = text_size * 0.22, lineheight = 0.9) +
      ggplot2::scale_colour_manual(values = forest_colours, guide = "none") +
      ggplot2::scale_shape_manual(values = forest_shapes, guide = "none") +
      ggplot2::scale_x_continuous(limits = x_lim, breaks = x_breaks) +
      ggplot2::labs(x = x_lab, y = NULL, title = "Effect Size Forest Plot") +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
    apply_plot_labels(p,
      input$tg_plot_title    %||% "",
      input$tg_plot_subtitle %||% "",
      input$tg_plot_title_align %||% "center",
      text_size,
      input$tg_plot_xaxis_label  %||% "",
      input$tg_plot_yaxis_label  %||% "",
      input$tg_plot_legend_label %||% "")
  })

  output$tg_plot_effect_forest <- renderPlot({
    p <- tg_forest_plot()
    if (!is.null(p)) {
      p
    } else {
      plot(1, 1, type="n", main="Effect size forest plot not available\n(requires pairwise comparison data)", xlab="", ylab="")
    }
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)
  
  # Diagnostics
  # Diagnostics tab: Q-Q plot
  output$tg_diag_qq <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    if (!is.null(rv$tg_results$diag_qq_plot)) {
      return(rv$tg_results$diag_qq_plot)
    }
    plot(1, 1, type = "n", main = "Q-Q plot not available", xlab = "", ylab = "")
  })

  # Diagnostics tab: Residuals vs Fitted
  output$tg_diag_resid_fitted <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    if (!is.null(rv$tg_results$diag_resid_fitted_plot)) {
      return(rv$tg_results$diag_resid_fitted_plot)
    }
    plot(1, 1, type = "n", main = "Residuals vs Fitted not available", xlab = "", ylab = "")
  })
  
  # Download handlers (placeholder implementations)
  output$download_tg_anova <- downloadHandler(
    filename = "anova_results.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$anova)) {
        write.csv(as.data.frame(rv$tg_results$anova), file, row.names = TRUE)
      }
    }
  )
  
  output$download_tg_comparisons <- downloadHandler(
    filename = "pairwise_comparisons.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$pairwise_comparisons)) {
        write.csv(as.data.frame(rv$tg_results$pairwise_comparisons), file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_effects <- downloadHandler(
    filename = "treatment_effects.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$treatment_effects)) {
        write.csv(as.data.frame(rv$tg_results$treatment_effects), file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_rates <- downloadHandler(
    filename = "growth_rates.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$growth_rates)) {
        df <- as.data.frame(rv$tg_results$growth_rates)
        if ("growth_rate" %in% colnames(df)) {
          df$daily_fold_change  <- exp(df$growth_rate)
          df$pct_change_per_day <- (exp(df$growth_rate) - 1) * 100
        }
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_auc_individual <- downloadHandler(
    filename = "auc_individual_data.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$auc_analysis$individual)) {
        write.csv(as.data.frame(rv$tg_results$auc_analysis$individual), file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_auc_summary <- downloadHandler(
    filename = "auc_summary.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$auc_analysis$summary)) {
        write.csv(as.data.frame(rv$tg_results$auc_analysis$summary), file, row.names = FALSE)
      }
    }
  )
  
  output$download_tg_data_summary <- downloadHandler(
    filename = "data_summary.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$data_summary)) {
        write.csv(as.data.frame(rv$tg_results$data_summary), file, row.names = FALSE)
      }
    }
  )
  
  # Plot download handlers - save as PNG
  output$download_tg_plot_growth <- downloadHandler(
    filename = "tumor_growth_curves.png",
    content = function(file) {
      p <- tg_growth_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_rates <- downloadHandler(
    filename = "growth_rate_distribution.png",
    content = function(file) {
      p <- tg_rates_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_auc <- downloadHandler(
    filename = "treatment_effects.png",
    content = function(file) {
      p <- tg_auc_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_caterpillar <- downloadHandler(
    filename = "qq_plot.png",
    content = function(file) {
      if (!is.null(rv$tg_results$caterpillar_plot)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(rv$tg_results$caterpillar_plot)
        dev.off()
      }
    }
  )

  output$download_tg_plot_effect_forest <- downloadHandler(
    filename = "effect_size_forest_plot.png",
    content = function(file) {
      p <- tg_forest_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 900,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )

  output$download_tg_plot_waterfall <- downloadHandler(
    filename = "waterfall_plot.png",
    content = function(file) {
      p <- tg_wf_plot()
      png(file, width  = input$tg_plot_width  %||% 800,
                height = input$tg_plot_height %||% 550,
                units  = "px", res = input$tg_plot_dpi %||% 100)
      tryCatch(print(p), finally = dev.off())
    }
  )

}
