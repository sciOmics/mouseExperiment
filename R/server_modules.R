#' Server Modules for mouseExperimentDashboard

#' Data Upload Server Logic
#' @keywords internal
data_upload_server <- function(input, output, session, rv) {
  
  # Reactive value to store uploaded data
  rv$uploaded_data <- NULL
  rv$data_valid <- FALSE
  rv$detected_cols <- NULL
  
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
      groups <- sort(unique(na.omit(as.character(data[[tx_col]]))))
      div(class = "panel panel-default",
        div(class = "panel-heading", tags$b("Treatment Groups")),
        div(class = "panel-body",
          tags$ul(lapply(groups, tags$li))
        )
      )
    } else div(class = "text-muted", "(Treatment column not configured)")
    
    # --- Cage IDs ---
    cage_box <- if (!is.null(cage_col) && cage_col %in% colnames(data)) {
      cages <- sort(unique(na.omit(as.character(data[[cage_col]]))))
      div(class = "panel panel-default",
        div(class = "panel-heading", tags$b("Cage IDs")),
        div(class = "panel-body",
          p(paste(cages, collapse = ", "))
        )
      )
    } else div(class = "text-muted", "(Cage column not configured)")
    
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
      
      inner <- if (!is.null(cage_col) && cage_col %in% colnames(data)) {
        # Cage × treatment table (includes both dimensions)
        ct <- unique(data[, c(id_col, cage_col, tx_col)])
        ct_counts <- as.data.frame(table(
          Cage = as.character(ct[[cage_col]]),
          Treatment = as.character(ct[[tx_col]])
        ))
        ct_counts <- ct_counts[ct_counts$Freq > 0, ]
        ct_rows <- lapply(seq_len(nrow(ct_counts)), function(i)
          tags$tr(tags$td(ct_counts$Cage[i]),
                  tags$td(ct_counts$Treatment[i]),
                  tags$td(ct_counts$Freq[i]))
        )
        tags$table(
          class = "table table-condensed table-bordered",
          style = "width:auto;",
          tags$thead(tags$tr(tags$th("Cage"), tags$th("Treatment"), tags$th("Mice"))),
          tags$tbody(ct_rows)
        )
      } else {
        # Fallback: per-treatment only when no cage column configured
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
      
      div(class = "panel panel-default",
        div(class = "panel-heading", tags$b("Mice per Group")),
        div(class = "panel-body", inner)
      )
    } else NULL
    
    tagList(
      hr(),
      wellPanel(
        h4("Dataset Summary"),
        p(strong(paste(nrow(data), "rows ×", ncol(data), "columns"))),
        fluidRow(
          column(3, tx_box),
          column(3, cage_box),
          column(3, time_box),
          column(3, mice_box)
        )
      ),
      hr()
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
             options = list(scrollX = TRUE, scrollY = "400px", pageLength = 25, paging = TRUE))
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
    rv$schedule_data <- NULL
    output$schedule_validation_msg <- renderUI(NULL)
    output$schedule_preview_ui <- renderUI(NULL)
  })
  
  # Treatment schedule file upload handler (optional)
  observeEvent(input$schedule_file, {
    req(input$schedule_file)
    sch_info <- input$schedule_file
    sch_path <- sch_info$datapath
    sch_parts <- strsplit(sch_info$name, "\\.")[[1]]
    sch_ext   <- tolower(sch_parts[length(sch_parts)])
    
    sch_result <- validate_uploaded_file(sch_path, sch_ext)
    
    if (sch_result$valid) {
      rv$schedule_data <- sch_result$data
      output$schedule_validation_msg <- renderUI({
        div(class = "alert alert-success",
            icon("check"),
            paste("Schedule loaded:", nrow(sch_result$data), "rows,",
                  ncol(sch_result$data), "cols."))
      })
    } else {
      rv$schedule_data <- NULL
      output$schedule_validation_msg <- renderUI(NULL)
      output$schedule_preview_ui <- renderUI(NULL)
      output$schedule_validation_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), sch_result$message)
      })
    }
  })
  
  # Schedule data preview — with Day column calculated if applicable
  preview_schedule <- reactive({
    req(rv$schedule_data)
    sch <- rv$schedule_data

    # Apply date -> Day only when:
    #  1. main data is using "from_date" mode
    #  2. the schedule doesn't already have a "Day" column
    #  3. the schedule has a column matching the date column chosen for the main data
    if (!is.null(input$time_method) && input$time_method == "from_date" &&
        !"Day" %in% colnames(sch) &&
        !is.null(input$col_date_for_day) && input$col_date_for_day != "" &&
        input$col_date_for_day %in% colnames(sch)) {
      tryCatch({
        first_val <- as.character(sch[1, input$col_date_for_day])
        fmt <- if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_val))  "%m/%d/%Y"
               else if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", first_val)) "%Y-%m-%d"
               else if (grepl("^\\d{1,2}-[A-Za-z]{3}$", first_val))  "%d-%b"
               else NULL
        sch <- mouseExperiment::calculate_dates(
          sch,
          start_date  = input$reference_date,
          date_column = input$col_date_for_day,
          date_format = fmt
        )
      }, error = function(e) NULL)
    }
    sch
  })

  # Schedule data preview
  output$schedule_preview_ui <- renderUI({
    if (is.null(rv$schedule_data)) return(NULL)
    tagList(
      hr(),
      wellPanel(
        h4("Treatment Schedule Preview"),
        fluidRow(
          column(12,
            p(strong(textOutput("schedule_preview_summary"))),
            DTOutput("schedule_table")
          )
        )
      )
    )
  })

  output$schedule_preview_summary <- renderText({
    if (is.null(rv$schedule_data)) return("")
    sch <- preview_schedule()
    paste("Schedule:", nrow(sch), "rows \u00d7",
          ncol(sch), "columns \u2014",
          paste(colnames(sch), collapse = ", "))
  })

  output$schedule_table <- renderDT({
    req(rv$schedule_data)
    datatable(preview_schedule(),
             options = list(scrollX = TRUE, scrollY = "300px",
                           pageLength = 25, paging = TRUE))
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
        cage_column = col_names$cage_column,
        model_type = input$tg_model_type,
        transform = input$tg_transform,
        polynomial_degree = input$tg_poly_degree,
        random_effects_specification = input$tg_random_effects,
        handle_cage_effects = input$tg_cage_effect,
        reference_group = input$tg_reference,  # Already validated to not be empty
        return_model = TRUE,
        include_diagnostics = TRUE,
        plots = FALSE
      )
    })
    
    # Generate plots from the results
    if (!isTRUE(result$error)) {
      tryCatch({
        # Plot 1: Tumor growth curves from original data
        if (!is.null(result$growth_rates)) {
          result$tumor_growth_plot <- mouseExperiment::plot_tumor_growth(
            df = data,
            volume_column = col_names$volume_column,
            day_column = col_names$time_column,
            treatment_column = col_names$treatment_column,
            ID_column = col_names$id_column,
            cage_column = col_names$cage_column,
            survival_column = if (!is.null(col_names$survival_column)) col_names$survival_column else NULL,
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
            result$auc_plot <- mouseExperiment::plot_auc(
              result$auc_analysis$individual,
              show_mean = TRUE,
              error_bar_type = "CI"
            ) + ggplot2::theme_classic()
          }, error = function(e) {
            cat("WARNING: plot_auc() failed:", e$message, "\n")
          })
        }
        if (is.null(result$auc_plot) && !is.null(result$treatment_effects)) {
          effects_df <- as.data.frame(result$treatment_effects)
          if (nrow(effects_df) > 0) {
            numeric_cols <- which(sapply(effects_df, is.numeric))
            if (length(numeric_cols) > 0) {
              effect_col <- numeric_cols[1]
              plot_df <- data.frame(
                Treatment = rownames(effects_df),
                Effect = as.numeric(effects_df[[effect_col]])
              )
              result$auc_plot <- ggplot2::ggplot(plot_df,
                                                  ggplot2::aes(x = reorder(Treatment, Effect),
                                                               y = Effect, fill = Treatment)) +
                ggplot2::geom_col(width = 0.6) +
                ggplot2::labs(title = "Treatment Effects (Adjusted Means)",
                             x = "Treatment", y = "Effect Size") +
                ggplot2::theme_classic() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                              legend.position = "none")
            }
          }
        }
        
        # Plot 4: QQ-plot for residual diagnostics (caterpillar equivalent)
        if (!is.null(result$diagnostics) && !is.null(result$diagnostics$residuals)) {
          residuals <- result$diagnostics$residuals$residuals
          # Create QQ plot data
          qq_data <- qqnorm(residuals, plot.it = FALSE)
          qq_df <- data.frame(
            theoretical = qq_data$x,
            sample = qq_data$y
          )
          
          result$caterpillar_plot <- ggplot2::ggplot(qq_df,
                                                      ggplot2::aes(x = theoretical,
                                                                   y = sample)) +
            ggplot2::geom_point(alpha = 0.6) +
            ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
            ggplot2::theme_classic() +
            ggplot2::labs(title = "Q-Q Plot of Residuals",
                         x = "Theoretical Quantiles",
                         y = "Sample Quantiles")
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
        
        # Type 2: Simple (non-adjusted) pairwise comparisons
        result$simple_pairwise_comparisons <- as.data.frame(
          pairs(emm, adjust = "none")
        )
        
      }, error = function(e) {
        cat("WARNING: Error generating comparison types:", e$message, "\n")
      })
    }
    
    rv$tg_results <- result
    rv$tg_running <- FALSE
    
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
      paste("  Polynomial Degree:", summary_list$model_specification$polynomial_degree),
      "",
      "NOTES:",
      paste("  ", paste(summary_list$notes, collapse="\n  "))
    )
    
    paste(summary_lines, collapse = "\n")
  })
  
  # ANOVA table
  output$tg_anova_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$anova)) {
      return(datatable(data.frame(Message = "ANOVA table not available")))
    }
    
    anova_df <- as.data.frame(rv$tg_results$anova)
    # Round numeric columns for better display
    numeric_cols <- sapply(anova_df, is.numeric)
    anova_df[numeric_cols] <- round(anova_df[numeric_cols], 4)
    
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
      # Reference group comparisons (default from mouseExperiment)
      if (!is.null(rv$tg_results$pairwise_comparisons)) {
        comp_df <- as.data.frame(rv$tg_results$pairwise_comparisons)
        table_title <- "Comparisons: vs Reference Group"
      }
    } else if (comp_type == "pairwise") {
      # All pairwise comparisons (unadjusted)
      if (!is.null(rv$tg_results$simple_pairwise_comparisons)) {
        comp_df <- as.data.frame(rv$tg_results$simple_pairwise_comparisons)
        table_title <- "Comparisons: All Pairwise (Unadjusted)"
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
    
    # Round numeric columns for better display
    numeric_cols <- sapply(comp_df, is.numeric)
    comp_df[numeric_cols] <- round(comp_df[numeric_cols], 4)
    
    # Create and return datatable
    datatable(comp_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Treatment effects table
  output$tg_effects_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$treatment_effects)) {
      return(datatable(data.frame(Message = "Treatment effects not available")))
    }
    
    effects_df <- as.data.frame(rv$tg_results$treatment_effects)
    # Round numeric columns
    numeric_cols <- sapply(effects_df, is.numeric)
    effects_df[numeric_cols] <- round(effects_df[numeric_cols], 4)
    
    datatable(effects_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Growth rates table
  output$tg_growth_rates_table <- renderDT({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error) || is.null(rv$tg_results$growth_rates)) {
      return(datatable(data.frame(Message = "Growth rates not available")))
    }
    
    rates_df <- as.data.frame(rv$tg_results$growth_rates)
    # Round numeric columns
    numeric_cols <- sapply(rates_df, is.numeric)
    rates_df[numeric_cols] <- round(rates_df[numeric_cols], 4)
    
    datatable(rates_df, options = list(scrollX = TRUE, pageLength = 15))
  })
  
  # Growth curves plot
  output$tg_plot_growth <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    
    if (!is.null(rv$tg_results$tumor_growth_plot)) {
      return(rv$tg_results$tumor_growth_plot)
    } else {
      plot(1, 1, type="n", main="Growth curves plot - generating...", xlab="", ylab="")
    }
  })
  
  # Growth rate distribution plot
  output$tg_plot_rates <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    
    if (!is.null(rv$tg_results$growth_rate_plot)) {
      return(rv$tg_results$growth_rate_plot)
    } else {
      plot(1, 1, type="n", main="Growth rate distribution not available", xlab="", ylab="")
    }
  })
  
  # AUC/Treatment effects plot
  output$tg_plot_auc <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    
    if (!is.null(rv$tg_results$auc_plot)) {
      return(rv$tg_results$auc_plot)
    } else {
      plot(1, 1, type="n", main="Treatment effects plot generating...", xlab="", ylab="")
    }
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
  
  # Diagnostics
  output$tg_diagnostics <- renderText({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return("Diagnostics not available")
    
    if (!is.null(rv$tg_results$diagnostics)) {
      paste(capture.output(print(rv$tg_results$diagnostics)), collapse = "\n")
    } else {
      "No diagnostic information available"
    }
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
        write.csv(as.data.frame(rv$tg_results$growth_rates), file, row.names = FALSE)
      }
    }
  )
  
  # Plot download handlers - save as PNG
  output$download_tg_plot_growth <- downloadHandler(
    filename = "tumor_growth_curves.png",
    content = function(file) {
      if (!is.null(rv$tg_results$tumor_growth_plot)) {
        png(file, width = 800, height = 600)
        print(rv$tg_results$tumor_growth_plot)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_rates <- downloadHandler(
    filename = "growth_rate_distribution.png",
    content = function(file) {
      if (!is.null(rv$tg_results$growth_rate_plot)) {
        png(file, width = 800, height = 600)
        print(rv$tg_results$growth_rate_plot)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_auc <- downloadHandler(
    filename = "treatment_effects.png",
    content = function(file) {
      if (!is.null(rv$tg_results$auc_plot)) {
        png(file, width = 800, height = 600)
        print(rv$tg_results$auc_plot)
        dev.off()
      }
    }
  )
  
  output$download_tg_plot_caterpillar <- downloadHandler(
    filename = "qq_plot.png",
    content = function(file) {
      if (!is.null(rv$tg_results$caterpillar_plot)) {
        png(file, width = 800, height = 600)
        print(rv$tg_results$caterpillar_plot)
        dev.off()
      }
    }
  )
  
  # Treatment schedule plot
  output$tg_plot_treatment_schedule <- renderPlot({
    if (is.null(rv$schedule_data)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "No treatment schedule uploaded (see Step 1b in Data Upload)")
      return(invisible(NULL))
    }
    req(rv$data_valid)
    
    # Use the preprocessed data (with Day + Volume cols) if available;
    # otherwise fall back to raw uploaded data
    tg_data  <- if (!is.null(rv$processed_data)) rv$processed_data else rv$uploaded_data
    time_col <- if (!is.null(rv$time_col)) rv$time_col else
                  if (!is.null(input$col_day) && input$col_day != "") input$col_day else "Day"
    tx_col   <- if (!is.null(rv$tx_col)) rv$tx_col else
                  if (!is.null(input$col_treatment) && input$col_treatment != "") input$col_treatment else "Treatment"
    
    # Determine which Day/Treatment columns exist in the schedule file
    sch <- rv$schedule_data
    sch_time_col <- if (time_col %in% colnames(sch)) time_col else
                     if ("Day" %in% colnames(sch)) "Day" else colnames(sch)[1]
    sch_tx_col   <- if (tx_col %in% colnames(sch)) tx_col else
                     if ("Treatment" %in% colnames(sch)) "Treatment" else colnames(sch)[2]
    
    # Coerce Day column to numeric (CSV uploads read all columns as character)
    sch[[sch_time_col]] <- suppressWarnings(as.numeric(sch[[sch_time_col]]))
    if (time_col %in% colnames(tg_data)) {
      tg_data[[time_col]] <- suppressWarnings(as.numeric(tg_data[[time_col]]))
    }
    
    tryCatch(
      mouseExperiment::plot_treatments(
        treatment_data    = sch,
        tumor_growth_data = tg_data,
        day_column        = sch_time_col,
        treatment_column  = sch_tx_col
      ),
      error = function(e) {
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = paste("Treatment schedule plot error:", e$message))
      }
    )
  })
  
  output$download_tg_plot_schedule <- downloadHandler(
    filename = "treatment_schedule.png",
    content = function(file) {
      req(rv$schedule_data, rv$data_valid)
      tg_data  <- if (!is.null(rv$processed_data)) rv$processed_data else rv$uploaded_data
      time_col <- if (!is.null(rv$time_col)) rv$time_col else "Day"
      tx_col   <- if (!is.null(rv$tx_col)) rv$tx_col else "Treatment"
      sch <- rv$schedule_data
      sch_time_col <- if (time_col %in% colnames(sch)) time_col else
                       if ("Day" %in% colnames(sch)) "Day" else colnames(sch)[1]
      sch_tx_col   <- if (tx_col %in% colnames(sch)) tx_col else
                       if ("Treatment" %in% colnames(sch)) "Treatment" else colnames(sch)[2]
      sch[[sch_time_col]] <- suppressWarnings(as.numeric(sch[[sch_time_col]]))
      if (time_col %in% colnames(tg_data)) {
        tg_data[[time_col]] <- suppressWarnings(as.numeric(tg_data[[time_col]]))
      }
      png(file, width = 1000, height = 400)
      tryCatch(
        print(mouseExperiment::plot_treatments(
          treatment_data    = sch,
          tumor_growth_data = tg_data,
          day_column        = sch_time_col,
          treatment_column  = sch_tx_col
        )),
        error = function(e) cat("Error generating schedule plot:", e$message, "\n")
      )
      dev.off()
    }
  )
}
