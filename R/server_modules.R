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
  
  # Data preview
  output$data_preview_ui <- renderUI({
    if (!rv$data_valid) return(NULL)
    data_preview_ui()
  })
  
  output$preview_summary <- renderText({
    if (!rv$data_valid) return("")
    paste("Showing first 20 rows of", nrow(rv$uploaded_data), 
          "total rows")
  })
  
  output$data_table <- renderDT({
    req(rv$uploaded_data)
    datatable(head(rv$uploaded_data, 20),
             options = list(scrollX = TRUE, pageLength = 5))
  })
  
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
  
  # Update reference group choices when data is available
  observe({
    cat("DEBUG: Checking if reference group should be updated...\n")
    cat("  Data valid?", rv$data_valid, "\n")
    cat("  Data null?", is.null(rv$uploaded_data), "\n")
    
    if (!is.null(rv$uploaded_data) && rv$data_valid) {
      treatment_col <- rv$detected_cols$treatment
      cat("  Treatment column:", treatment_col, "\n")
      cat("  All detected columns:", paste(names(rv$detected_cols), "=", 
                                          unlist(rv$detected_cols), collapse="; "), "\n")
      
      if (!is.null(treatment_col) && treatment_col != "") {
        # Get the actual column from the data
        col_data <- rv$uploaded_data[[treatment_col]]
        cat("  Column data class:", class(col_data), "\n")
        cat("  Column data length:", length(col_data), "\n")
        cat("  First 5 values:", paste(head(col_data, 5), collapse=", "), "\n")
        cat("  Any NA?", any(is.na(col_data)), "\n")
        
        treatment_values <- unique(na.omit(col_data))
        cat("  Unique treatment values:", paste(treatment_values, collapse=", "), "\n")
        cat("  Number of unique values:", length(treatment_values), "\n")
        
        if (length(treatment_values) > 0) {
          # Update with proper choice format
          new_choices <- c("-- Select a reference group --" = "", 
                          setNames(as.character(treatment_values), 
                                  as.character(treatment_values)))
          updateSelectInput(session, "tg_reference", choices = new_choices)
          cat("  Reference group updated with choices\n")
        } else {
          cat("  WARNING: No treatment values found after unique/na.omit\n")
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
      data <- mouseExperiment::calculate_dates(
        data,
        start_date = start_date,
        date_column = input$col_date_for_day,
        date_format = "%m/%d/%Y"  # Assuming MM/DD/YYYY format based on sample data
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
    
    # Prepare column names
    col_names <- list(
      time_column = time_col,
      volume_column = vol_col,
      id_column = input$col_id,
      treatment_column = input$col_treatment,
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
        plots = TRUE
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
          )
        }
        
        # Plot 2: Growth rate distribution
        if (!is.null(result$growth_rates) && nrow(result$growth_rates) > 0) {
          result$growth_rate_plot <- mouseExperiment::plot_growth_rate(
            growth_data = result$growth_rates,
            show_mean = TRUE,
            error_bar_type = "SEM"
          )
        }
        
        # Plot 3: Treatment effects plot (AUC equivalent for LME4 analysis)
        if (!is.null(result$treatment_effects)) {
          effects_df <- as.data.frame(result$treatment_effects)
          if (nrow(effects_df) > 0) {
            # Try to convert to a plottable format
            if (ncol(effects_df) >= 1) {
              # Get the first numeric column as the effect size
              numeric_cols <- which(sapply(effects_df, is.numeric))
              if (length(numeric_cols) > 0) {
                effect_col <- numeric_cols[1]
                plot_df <- data.frame(
                  Treatment = rownames(effects_df),
                  Effect = as.numeric(effects_df[[effect_col]])
                )
                
                result$auc_plot <- ggplot2::ggplot(plot_df, 
                                                    ggplot2::aes(x = Treatment, 
                                                                 y = Effect, 
                                                                 fill = Treatment)) +
                  ggplot2::geom_col() +
                  ggplot2::theme_minimal() +
                  ggplot2::labs(title = "Treatment Effects (Adjusted Means)",
                              x = "Treatment",
                              y = "Effect Size",
                              fill = "Treatment") +
                  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                               legend.position = "none")
              }
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
            ggplot2::theme_minimal() +
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
          emmeans::pairs(emm, adjust = "tukey")
        )
        
        # Type 2: Simple (non-adjusted) pairwise comparisons
        result$simple_pairwise_comparisons <- as.data.frame(
          emmeans::pairs(emm, adjust = "none")
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
}
