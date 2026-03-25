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
        first_val  <- as.character(data[1, input$col_date_for_day])
        fmt <- if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_val)) "%m/%d/%Y"
               else if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", first_val)) "%Y-%m-%d"
               else if (grepl("^\\d{1,2}-[A-Za-z]{3}$", first_val)) "%d-%b"
               else NULL
        data <- mouseExperiment::calculate_dates(
          data,
          start_date  = input$reference_date,
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
    
    tagList(
      wellPanel(
        h4("Dataset Summary"),
        p(strong(paste(nrow(data), "rows ×", ncol(data), "columns"))),
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
    datatable(preview_data(),
             options = list(scrollX = TRUE, scrollY = "600px", paging = FALSE, dom = "ti"))
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
      data <- mouseExperiment::calculate_volume(
        data,
        length_column = input$col_length,
        width_column = input$col_width,
        formula = input$volume_formula
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
                # Individual data beeswarm + mean-CI overlay from LME
                gr_df$plot_group <- as.character(gr_df[[grp_col_gr]])
                gr_df$plot_val   <- as.numeric(gr_df[[val_col_gr]])
                # Keep group ordering consistent with treatment_effects
                grp_levels <- as.character(plot_df$Treatment)
                gr_df$plot_group <- factor(gr_df$plot_group, levels = grp_levels)
                plot_df$Treatment_f <- factor(plot_df$Treatment, levels = grp_levels)

                p_fallback <- ggplot2::ggplot(gr_df,
                                              ggplot2::aes(x = plot_group, y = plot_val, colour = plot_group)) +
                  ggbeeswarm::geom_beeswarm(size = 2.5, alpha = 0.8, cex = 2.5) +
                  ggplot2::labs(title = "Individual Growth Rates by Treatment Group",
                               x = "Treatment", y = "Growth Rate") +
                  ggplot2::theme_classic() +
                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
                                legend.position = "none")

                if ("Lower" %in% colnames(plot_df)) {
                  p_fallback <- p_fallback +
                    ggplot2::geom_pointrange(
                      data = plot_df,
                      ggplot2::aes(x = Treatment_f, y = Effect, ymin = Lower, ymax = Upper),
                      inherit.aes = FALSE,
                      colour = "black", size = 0.8, linewidth = 1.1
                    )
                } else {
                  p_fallback <- p_fallback +
                    ggplot2::stat_summary(fun = mean, geom = "crossbar",
                                         width = 0.45, colour = "black", linewidth = 0.7)
                }
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
        
        # Plot 5: Effect size forest plot.
        # Use reference_comparisons_dunnett (treatment - control, correct sign)
        # in preference to all_pairwise_comparisons (alphabetical pairs whose
        # sign depends on factor level order and is often backwards).
        pc <- if (!is.null(result$reference_comparisons_dunnett)) result$reference_comparisons_dunnett
              else if (!is.null(result$all_pairwise_comparisons)) result$all_pairwise_comparisons
              else if (!is.null(result$pairwise_comparisons)) result$pairwise_comparisons
              else NULL
        if (!is.null(pc)) {
          tryCatch({
            pc_df <- as.data.frame(pc)
            # Standardise column names from emmeans / model output
            label_col   <- intersect(c("contrast", "Contrast", "comparison", "pair"), colnames(pc_df))[1]
            est_col     <- intersect(c("estimate", "Estimate"), colnames(pc_df))[1]
            se_col      <- intersect(c("SE"), colnames(pc_df))[1]
            lower_col   <- intersect(c("lower.CL", "Lower", "lower_CL", "lower"), colnames(pc_df))[1]
            upper_col   <- intersect(c("upper.CL", "Upper", "upper_CL", "upper"), colnames(pc_df))[1]
            p_col       <- intersect(c("p.value", "p_value", "p.adj", "adj.p.value", "P_Value"), colnames(pc_df))[1]

            if (!is.na(label_col) && !is.na(est_col)) {
              # Compute CIs from SE if explicit CI columns are absent
              if (is.na(lower_col) && !is.na(se_col)) {
                pc_df$lower.CL <- pc_df[[est_col]] - 1.96 * pc_df[[se_col]]
                pc_df$upper.CL <- pc_df[[est_col]] + 1.96 * pc_df[[se_col]]
                lower_col <- "lower.CL"; upper_col <- "upper.CL"
              }
              p_label <- if (!is.na(p_col))
                           vapply(pc_df[[p_col]], format_pval, character(1))
                         else NULL

              # For dunnett contrasts the label is already "TreatA - Reference";
              # strip the reference name from the label so it reads more cleanly
              # in the forest plot (e.g. "HDACi - Control" → "HDACi").
              using_dunnett <- !is.null(result$reference_comparisons_dunnett) &&
                               identical(pc, result$reference_comparisons_dunnett)
              ref_grp_ef <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference))
                              input$tg_reference else NULL
              contrast_labels <- pc_df[[label_col]]
              if (isTRUE(using_dunnett) && !is.null(ref_grp_ef)) {
                # emmeans trt.vs.ctrl labels: "HDACi - Control" → "HDACi"
                suffix_pat <- paste0("\\s*[-–]\\s*", gsub("([+()])", "\\\\\\1", ref_grp_ef), "$")
                clean_labels <- sub(suffix_pat, "", contrast_labels)
                # Only apply the strip if every label had the suffix (safety check)
                if (all(clean_labels != contrast_labels)) contrast_labels <- clean_labels
              }

              forest_df <- data.frame(
                Contrast  = contrast_labels,
                Estimate  = pc_df[[est_col]],
                Lower     = if (!is.na(lower_col)) pc_df[[lower_col]] else pc_df[[est_col]],
                Upper     = if (!is.na(upper_col)) pc_df[[upper_col]] else pc_df[[est_col]],
                P_Label   = if (is.null(p_label)) "" else paste0("p=", p_label),
                stringsAsFactors = FALSE
              )
              # Add reference group row at Estimate = 0
              if (!is.null(ref_grp_ef)) {
                forest_df <- rbind(forest_df, data.frame(
                  Contrast = paste0(ref_grp_ef, "\n(Reference)"),
                  Estimate = 0, Lower = NA_real_, Upper = NA_real_,
                  P_Label  = "", stringsAsFactors = FALSE
                ))
              }
              forest_df$is_ref <- seq_len(nrow(forest_df)) == nrow(forest_df) & !is.null(ref_grp_ef)

              # Balanced symmetric axis: same extent either side of 0,
              # rounded up to the nearest whole unit (mirrors HR forest plot logic)
              all_vals <- c(forest_df$Estimate[!forest_df$is_ref],
                            forest_df$Lower[!forest_df$is_ref],
                            forest_df$Upper[!forest_df$is_ref])
              max_abs  <- ceiling(max(abs(all_vals), na.rm = TRUE))
              # Add 20% padding on the right so p-value labels aren't clipped
              x_lim <- c(-max_abs, max_abs * 1.20)
              x_breaks <- pretty(c(-max_abs, max_abs), n = 6)
              x_breaks <- x_breaks[x_breaks >= -max_abs & x_breaks <= max_abs]

              result$effect_forest_plot <- ggplot2::ggplot(
                forest_df,
                # reorder by -Estimate: most negative (most effective) at top
                ggplot2::aes(x = Estimate, y = stats::reorder(Contrast, -Estimate))
              ) +
                ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.6) +
                ggplot2::geom_errorbarh(
                  data = ~ subset(., !is_ref),
                  ggplot2::aes(xmin = Lower, xmax = Upper),
                  height = 0.25, colour = "steelblue"
                ) +
                ggplot2::geom_point(
                  data = ~ subset(., !is_ref), size = 3, colour = "steelblue"
                ) +
                ggplot2::geom_point(
                  data = ~ subset(., is_ref), size = 4, colour = "grey50", shape = 18
                ) +
                ggplot2::geom_text(
                  data = ~ subset(., !is_ref & nzchar(P_Label) & !is.na(Upper)),
                  ggplot2::aes(label = P_Label, x = Upper),
                  hjust = -0.2, size = 3.2
                ) +
                ggplot2::scale_x_continuous(
                  limits = x_lim,
                  breaks = x_breaks
                ) +
                ggplot2::labs(
                  title = "Effect Size Forest Plot",
                  x = if (isTRUE(using_dunnett))
                        "Mean Difference (treatment \u2212 reference)"
                      else "Mean Difference",
                  y = NULL
                ) +
                ggplot2::theme_classic() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
            }
          }, error = function(e) {
            cat("WARNING: effect forest plot failed:", e$message, "\n")
          })
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
        
        # Type 1: All pairwise comparisons with Tukey adjustment
        result$all_pairwise_comparisons <- as.data.frame(
          pairs(emm, adjust = "tukey")
        )
        
        # Type 2: vs Reference Group with Dunnett correction
        ref_group <- if (!is.null(input$tg_reference) && nzchar(input$tg_reference))
          input$tg_reference else levels(emm)[[1]]
        ref_idx <- which(levels(emm) == ref_group)
        if (length(ref_idx) == 1) {
          result$reference_comparisons_dunnett <- as.data.frame(
            emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_idx, adjust = "dunnett")
          )
        }
        
      }, error = function(e) {
        cat("WARNING: Error generating comparison types:", e$message, "\n")
      })
    }
    
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
  
  # ANOVA table
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
    
    datatable(anova_df, options = list(scrollX = TRUE, pageLength = 10))
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
    transform_used   <- rv$tg_results$methods$volume_transformation %||% "none"
    model_type_used  <- input$tg_model_type %||% "lme4"
    if (model_type_used == "lme4" && transform_used == "log") {
      est_col <- intersect(c("estimate", "Estimate"), colnames(comp_df))[1]
      lo_col  <- intersect(c("lower.CL", "lower_CL", "Lower"), colnames(comp_df))[1]
      hi_col  <- intersect(c("upper.CL", "upper_CL", "Upper"), colnames(comp_df))[1]
      if (!is.na(est_col)) {
        comp_df$Fold_Change <- round(exp(as.numeric(comp_df[[est_col]])), 3)
        if (!is.na(lo_col)) comp_df$Fold_Change_Lower <- round(exp(as.numeric(comp_df[[lo_col]])), 3)
        if (!is.na(hi_col)) comp_df$Fold_Change_Upper <- round(exp(as.numeric(comp_df[[hi_col]])), 3)
        # Insert fold-change columns right after the estimate column
        est_idx  <- which(colnames(comp_df) == est_col)
        fc_idxs  <- match(c("Fold_Change", "Fold_Change_Lower", "Fold_Change_Upper"),
                          colnames(comp_df))
        fc_idxs  <- fc_idxs[!is.na(fc_idxs)]
        rest     <- setdiff(seq_len(ncol(comp_df)), c(est_idx, fc_idxs))
        comp_df  <- comp_df[, c(rest[rest <= est_idx], est_idx, fc_idxs, rest[rest > est_idx])]
      }
    }

    comp_df <- format_df_for_display(comp_df)
    datatable(comp_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Treatment effects table
  output$tg_effects_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$treatment_effects)) {
      return(datatable(data.frame(Message = "Treatment effects not available")))
    }
    
    effects_df <- as.data.frame(rv$tg_results$treatment_effects)

    # Add back-transformed columns when LME4 model used a transformation
    transform_used <- rv$tg_results$methods$volume_transformation %||% "none"
    bt_cols <- intersect(c("Adjusted_Mean", "Lower_CL", "Upper_CL"), colnames(effects_df))
    if (transform_used %in% c("log", "sqrt") && length(bt_cols) > 0) {
      for (col in bt_cols) {
        effects_df[[paste0(col, "_orig")]] <- round(back_transform(effects_df[[col]], transform_used), 3)
      }
      # Interleave _orig columns immediately after their source
      col_order <- character(0)
      for (col in colnames(effects_df)) {
        col_order <- c(col_order, col)
        orig <- paste0(col, "_orig")
        if (orig %in% colnames(effects_df)) col_order <- c(col_order, orig)
      }
      effects_df <- effects_df[, col_order]
    }

    effects_df <- format_df_for_display(effects_df)
    # Improve column headers
    scale_note <- switch(transform_used, "log" = " (log scale)", "sqrt" = " (sqrt scale)", "")
    orig_note  <- " (original units)"
    colnames(effects_df) <- gsub("Adjusted_Mean$",       paste0("Adj_Mean", scale_note), colnames(effects_df))
    colnames(effects_df) <- gsub("Adjusted_Mean_orig$",  paste0("Adj_Mean", orig_note),  colnames(effects_df))
    colnames(effects_df) <- gsub("Lower_CL$",            paste0("Lower_CL", scale_note), colnames(effects_df))
    colnames(effects_df) <- gsub("Lower_CL_orig$",       paste0("Lower_CL", orig_note),  colnames(effects_df))
    colnames(effects_df) <- gsub("Upper_CL$",            paste0("Upper_CL", scale_note), colnames(effects_df))
    colnames(effects_df) <- gsub("Upper_CL_orig$",       paste0("Upper_CL", orig_note),  colnames(effects_df))
    datatable(effects_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Growth rates table
  output$tg_growth_rates_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$growth_rates)) {
      return(datatable(data.frame(Message = "Growth rates not available")))
    }
    
    rates_df <- as.data.frame(rv$tg_results$growth_rates)
    if ("growth_rate" %in% colnames(rates_df)) {
      rates_df$daily_fold_change  <- exp(rates_df$growth_rate)
      rates_df$pct_change_per_day <- (exp(rates_df$growth_rate) - 1) * 100
      # Reorder: put both new columns right after growth_rate
      gr_idx    <- which(colnames(rates_df) == "growth_rate")
      new_idxs  <- match(c("daily_fold_change", "pct_change_per_day"), colnames(rates_df))
      other     <- setdiff(seq_len(ncol(rates_df)), c(gr_idx, new_idxs))
      rates_df  <- rates_df[, c(other[other < gr_idx], gr_idx, new_idxs, other[other > gr_idx])]
    }
    rates_df <- format_df_for_display(rates_df)
    datatable(rates_df,
      options  = list(scrollX = TRUE, paging = FALSE, scrollY = "500px", scrollCollapse = TRUE),
      colnames = setNames(colnames(rates_df),
        ifelse(colnames(rates_df) == "growth_rate",          "growth_rate (log-scale/day)",
        ifelse(colnames(rates_df) == "daily_fold_change",    "daily_fold_change (×/day)",
        ifelse(colnames(rates_df) == "pct_change_per_day",   "% change/day",
        colnames(rates_df))))))
  })
  
  # ---- Group order drag-and-drop for Plots tab ----
  output$tg_group_order_ui <- renderUI({
    req(rv$tg_plot_groups)
    group_order_ui("tg", rv$tg_plot_groups)
  })

  # ---- Colour picker UI for Plots tab ----
  output$tg_color_pickers_ui <- renderUI({
    req(rv$tg_plot_groups)
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    group_color_pickers_ui("tg_color", ordered_grps)
  })

  # Growth curves plot  — reactive: re-renders on colour/option changes
  output$tg_plot_growth <- renderPlot({
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
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_color_manual(values = colors)
    }
    # Apply text size and legend position
    text_size  <- input$tg_opt_text_size %||% 12
    legend_pos <- input$tg_opt_legend_pos %||% "auto"
    p <- p + ggplot2::theme(text = ggplot2::element_text(size = text_size))
    p <- apply_legend_pos(p, legend_pos)
    # Apply y-axis scale
    yscale <- input$tg_opt_yscale %||% "linear"
    if (yscale == "log10") p <- p + ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.05)))
    # Apply x-axis start
    x_start <- local({
      v <- input$tg_opt_xstart
      if (is.numeric(v) && length(v) == 1L && !is.na(v) && v > 0) v else NULL
    })
    if (!is.null(x_start)) p <- p + ggplot2::coord_cartesian(xlim = c(x_start, NA))
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
  
  # Growth rate distribution plot — reactive
  output$tg_plot_rates <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    req(rv$tg_results$growth_rates)
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    colors     <- build_color_vector(input, "tg_color", ordered_grps)
    center_bar <- input$tg_opt_center_bar %||% "mean"
    show_mean  <- center_bar == "mean"
    err_type   <- input$tg_opt_error_bar  %||% "SEM"
    pt_size    <- input$tg_opt_point_size %||% 2.5
    jit_width  <- input$tg_opt_jitter_width %||% 0.2
    text_size  <- input$tg_opt_text_size %||% 12
    legend_pos <- input$tg_opt_legend_pos %||% "auto"
    rates      <- rv$tg_results$growth_rates
    # Apply user group order as factor levels
    grp_col_r  <- intersect(c("Group", "Treatment"), colnames(rates))[1]
    if (isTRUE(!is.na(grp_col_r)) && isTRUE(length(ordered_grps) > 0))
      rates[[grp_col_r]] <- factor(rates[[grp_col_r]], levels = ordered_grps)
    p <- tryCatch({
      mouseExperiment::plot_growth_rate(
        growth_data    = rates,
        colors         = colors,
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
    p <- p + ggplot2::theme(text = ggplot2::element_text(size = text_size))
    p <- apply_legend_pos(p, legend_pos)
    p
  })
  
  # AUC/Treatment effects plot — reactive (beeswarm)
  output$tg_plot_auc <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    # Fall back to stored plot if AUC individual data not available
    if (is.null(rv$tg_results$auc_analysis) ||
        is.null(rv$tg_results$auc_analysis$individual)) {
      p <- rv$tg_results$auc_plot
      req(!is.null(p))
      return(p)
    }
    ordered_grps <- get_ordered_groups(input, "tg", rv$tg_plot_groups)
    colors     <- build_color_vector(input, "tg_color", ordered_grps)
    center_bar <- input$tg_opt_center_bar %||% "mean"
    err_type   <- input$tg_opt_error_bar  %||% "SEM"
    pt_size    <- input$tg_opt_point_size %||% 2.5
    text_size  <- input$tg_opt_text_size %||% 12
    legend_pos <- input$tg_opt_legend_pos %||% "auto"
    auc_ind    <- rv$tg_results$auc_analysis$individual
    grp_col_a  <- intersect(c("Treatment", "Group"), colnames(auc_ind))[1]
    if (isTRUE(!is.na(grp_col_a)) && isTRUE(length(ordered_grps) > 0)) {
      # The AUC individual data uses plain treatment names (no " - Dose: X" suffix).
      # When dose data is loaded, ordered_grps contains composite names — strip
      # the dose portion so the factor levels match the actual Treatment values.
      auc_levels <- unique(sub(" - Dose: .*$", "", ordered_grps))
      auc_levels <- auc_levels[auc_levels %in% as.character(unique(auc_ind[[grp_col_a]]))]
      if (length(auc_levels) > 0)
        auc_ind[[grp_col_a]] <- factor(auc_ind[[grp_col_a]], levels = auc_levels)
    }
    p <- tryCatch({
      p_b <- ggplot2::ggplot(auc_ind,
               ggplot2::aes(x = .data[[grp_col_a]], y = .data[["AUC"]],
                            color = .data[[grp_col_a]])) +
        ggplot2::theme_classic() +
        ggplot2::labs(x = NULL, y = "AUC", title = "AUC Comparison") +
        ggplot2::theme(legend.title = ggplot2::element_blank())
      # Error bars — drawn first so individual points render on top
      if (err_type != "none") {
        sem_fn  <- function(x) mean(x, na.rm=TRUE) + c(-1,1) * sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))
        sd_fn   <- function(x) mean(x, na.rm=TRUE) + c(-1,1) * sd(x, na.rm=TRUE)
        ci_fn   <- function(x) mean(x, na.rm=TRUE) + c(-1,1) * qt(0.975, df=max(sum(!is.na(x))-1,1)) * sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))
        eb_fn   <- switch(err_type, SEM=sem_fn, SD=sd_fn, CI=ci_fn, sem_fn)
        p_b <- p_b + ggplot2::stat_summary(
          ggplot2::aes(group = .data[[grp_col_a]]),
          fun.min = function(x) eb_fn(x)[1],
          fun     = function(x) mean(x, na.rm=TRUE),
          fun.max = function(x) eb_fn(x)[2],
          geom = "errorbar", width = 0.4, linewidth = 0.8, color = "black"
        )
      }
      # Center bar — degenerate crossbar (ymin = y = ymax) renders as a
      # horizontal bar at the centre value; drawn before points so it stays
      # under the individual data cloud.
      if (center_bar != "none") {
        cfun <- switch(center_bar,
          mean           = function(x) mean(x, na.rm=TRUE),
          median         = function(x) median(x, na.rm=TRUE),
          geometric_mean = function(x) exp(mean(log(pmax(x, 1e-9, na.rm=TRUE)))),
          function(x) mean(x, na.rm=TRUE)
        )
        p_b <- p_b + ggplot2::stat_summary(
          fun.min = cfun, fun = cfun, fun.max = cfun,
          geom = "crossbar", width = 0.5, linewidth = 0.8,
          ggplot2::aes(group = .data[[grp_col_a]], color = .data[[grp_col_a]])
        )
      }
      # Individual points — use geom_jitter (matches plot_auc() backend) so
      # the y aesthetic is never affected by beeswarm position computation.
      # height = 0 keeps all points at their true AUC values.
      p_b <- p_b + ggplot2::geom_jitter(
        position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 123),
        size = pt_size, alpha = 0.85
      )
      if (!is.null(colors)) p_b <- p_b + ggplot2::scale_color_manual(values = colors)
      p_b
    }, error = function(e) rv$tg_results$auc_plot)
    req(!is.null(p))
    p <- p + ggplot2::theme(text = ggplot2::element_text(size = text_size))
    p <- apply_legend_pos(p, legend_pos)
    p
  })
  
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

  # Effect size forest plot
  output$tg_plot_effect_forest <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    if (!is.null(rv$tg_results$effect_forest_plot)) {
      return(rv$tg_results$effect_forest_plot)
    } else {
      plot(1, 1, type="n", main="Effect size forest plot not available\n(requires pairwise comparison data)", xlab="", ylab="")
    }
  })
  
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
  
  # Plot download handlers - save as PNG
  output$download_tg_plot_growth <- downloadHandler(
    filename = "tumor_growth_curves.png",
    content = function(file) {
      if (!is.null(rv$tg_results$tumor_growth_plot)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(rv$tg_results$tumor_growth_plot)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_rates <- downloadHandler(
    filename = "growth_rate_distribution.png",
    content = function(file) {
      if (!is.null(rv$tg_results$growth_rate_plot)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(rv$tg_results$growth_rate_plot)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_auc <- downloadHandler(
    filename = "treatment_effects.png",
    content = function(file) {
      if (!is.null(rv$tg_results$auc_plot)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(rv$tg_results$auc_plot)
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
      if (!is.null(rv$tg_results$effect_forest_plot)) {
        png(file, width = input$tg_plot_width %||% 900,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(rv$tg_results$effect_forest_plot)
        dev.off()
      }
    }
  )
  
}
