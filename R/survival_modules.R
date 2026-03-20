#' Survival Analysis UI
#' @keywords internal
survivalUI <- function() {
  fluidPage(
    h2("Survival Analysis"),
    
    p("Analyze time-to-event data using Kaplan-Meier survival curves and proportional hazards models. ",
      "Ensure data is uploaded and configured in the Data Upload tab."),
    
    fluidRow(
      column(3,
        wellPanel(
          h4("Analysis Parameters"),
          
          checkboxInput("surv_include_risk_table",
                       "Show Risk Table",
                       value = TRUE),
          
          selectInput("surv_palette",
                     "Color Palette",
                     choices = c("Default" = "hue",
                                "Colorblind" = "colorblind",
                                "Viridis" = "viridis",
                                "Set1" = "Set1",
                                "Set2" = "Set2"),
                     selected = "hue"),
          
          uiOutput("surv_reference_group_ui"),
          
          br(),
          
          actionButton("run_surv_analysis", 
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),
          
          uiOutput("surv_status_msg")
        )
      ),
      
      column(9,
        tabsetPanel(
          id = "surv_results_tabs",
          
          tabPanel("Summary",
            shinycssloaders::withSpinner(
              verbatimTextOutput("surv_summary")
            )
          ),
          
          tabPanel("Cox Model",
            shinycssloaders::withSpinner(
              verbatimTextOutput("surv_cox_summary")
            ),
            shinycssloaders::withSpinner(
              DTOutput("surv_cox_table")
            ),
            downloadButton("download_surv_cox", "Download Results")
          ),
          
          tabPanel("Kaplan-Meier Curve",
            shinycssloaders::withSpinner(
              plotOutput("surv_plot_km", height = "500px")
            ),
            downloadButton("download_surv_plot_km", "Download Plot")
          ),
          
          tabPanel("Risk Table",
            shinycssloaders::withSpinner(
              DTOutput("surv_risk_table")
            ),
            downloadButton("download_surv_risk_table", "Download Table")
          )
        )
      )
    )
  )
}

#' Survival Analysis Server Logic
#' @keywords internal
survival_server <- function(input, output, session, rv) {
  
  rv$surv_results <- NULL
  rv$surv_running <- FALSE
  
  # Dynamically generate reference group dropdown based on treatment groups
  output$surv_reference_group_ui <- renderUI({
    req(rv$uploaded_data, rv$data_valid)
    
    # Get unique treatment values
    if (!is.null(input$col_treatment) && input$col_treatment != "") {
      treatment_col <- rv$uploaded_data[[input$col_treatment]]
      if (!is.null(treatment_col)) {
        treatment_groups <- sort(unique(as.character(na.omit(treatment_col))))
        
        if (length(treatment_groups) > 0) {
          return(
            selectInput("surv_reference_group",
                       "Reference Group",
                       choices = treatment_groups,
                       selected = treatment_groups[1])
          )
        }
      }
    }
    
    # Fallback if no data
    selectInput("surv_reference_group",
               "Reference Group",
               choices = "",
               selected = "")
  })
  
  # Run survival analysis
  observeEvent(input$run_surv_analysis, {
    req(rv$uploaded_data, rv$data_valid)
    
    # Validate that required columns are configured from data upload
    if (input$col_treatment == "" || input$col_treatment == "-- Select a reference group --") {
      output$surv_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error: Please select a Treatment Group on the Data Upload tab")
      })
      return()
    }
    

    if (is.null(input$col_survival_event) || input$col_survival_event == "" || input$col_survival_event == "None") {
      output$surv_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error: Please select an Event/Censor Status Column on the Data Upload tab")
      })
      return()
    }
    
    rv$surv_running <- TRUE
    output$surv_status_msg <- renderUI({
      div(class = "alert alert-info",
          icon("spinner", class = "fa-spin"), "Running survival analysis...")
    })
    
    # Prepare data - apply same time processing as tumor growth analysis
    data <- rv$uploaded_data
    
    # Handle time column - either use direct day or calculate from date
    # This mirrors the tumor growth analysis logic
    if (input$time_method == "from_date") {
      # Calculate days from reference date
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
    
    # For survival analysis, reduce data to one row per animal
    # Get unique animals (ID/Cage combination)
    if ("ID" %in% colnames(data) && "Cage" %in% colnames(data)) {
      # Aggregate by ID and Cage - take the last (most recent) measurement
      unique_animals <- unique(data[, c("ID", "Cage"), drop = FALSE])
      survival_data <- data.frame()
      
      for (i in 1:nrow(unique_animals)) {
        animal_rows <- data[
          data$ID == unique_animals[i, "ID"] & 
          data$Cage == unique_animals[i, "Cage"],
        ]
        # Take the last row for this animal
        if (nrow(animal_rows) > 0) {
          survival_data <- rbind(survival_data, animal_rows[nrow(animal_rows), ])
        }
      }
      
      # Validate aggregation
      if (nrow(survival_data) == 0) {
        result <- list(error = TRUE, message = "Data aggregation resulted in empty dataset. Check ID and Cage columns.")
        rv$surv_results <- result
        rv$surv_running <- FALSE
        output$surv_status_msg <- renderUI({
          div(class = "alert alert-danger",
              icon("exclamation-circle"), 
              "Error:", result$message)
        })
        return()
      }
      
      data <- survival_data
    } else {
      # If no ID/Cage columns, try to operate on data as-is (single measurement per animal)
      if (nrow(data) == 0) {
        result <- list(error = TRUE, message = "No data available for analysis.")
        rv$surv_results <- result
        rv$surv_running <- FALSE
        output$surv_status_msg <- renderUI({
          div(class = "alert alert-danger",
              icon("exclamation-circle"), 
              "Error: No data available for analysis")
        })
        return()
      }
    }
    
    # Run analysis with processed data and calculated/selected time column
    # Add diagnostic output
    cat("\n=== SURVIVAL ANALYSIS DIAGNOSTICS ===\n")
    cat("Data rows:", nrow(data), "\n")
    cat("Data columns:", paste(colnames(data), collapse = ", "), "\n")
    cat("Time column specified:", time_col, "\n")
    cat("Time column exists:", time_col %in% colnames(data), "\n")
    cat("Censor column specified:", input$col_survival_event, "\n")
    cat("Censor column exists:", input$col_survival_event %in% colnames(data), "\n")
    cat("Treatment column specified:", input$col_treatment, "\n")
    cat("Treatment column exists:", input$col_treatment %in% colnames(data), "\n")
    
    if (time_col %in% colnames(data)) {
      cat("Time values sample:", paste(head(data[[time_col]], 5), collapse = ", "), "\n")
      cat("Time NAs:", sum(is.na(data[[time_col]])), "\n")
    }
    
    if (input$col_survival_event %in% colnames(data)) {
      cat("Censor values sample:", paste(head(data[[input$col_survival_event]], 5), collapse = ", "), "\n")
      cat("Censor NAs:", sum(is.na(data[[input$col_survival_event]])), "\n")
    }
    cat("====================================\n\n")
    
    result <- safe_analysis({
      mouseExperiment::survival_statistics(
        df = data,
        time_column = time_col,
        censor_column = input$col_survival_event,
        treatment_column = input$col_treatment,
        reference_group = input$surv_reference_group,
        firth_correction = TRUE
      )
    })
    
    # Store time calculation info for later display
    result$time_method <- input$time_method
    result$time_column_used <- time_col
    result$reference_date <- input$reference_date
    
    # Generate Kaplan-Meier plot if analysis succeeded
    if (!isTRUE(result$error) && !is.null(result$survival_data)) {
      tryCatch({
        # Use the survival_data returned by survival_statistics (has Time, Event, Treatment columns)
        surv_df <- result$survival_data
        
        km_fit <- survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df)
        
        # Store the KM fit
        result$km_fit <- km_fit
        
        # Generate ggsurvplot (includes both plot and table)
        km_plot_obj <- survminer::ggsurvplot(
          km_fit,
          data = surv_df,
          palette = input$surv_palette,
          risk.table = input$surv_include_risk_table,
          risk.table.y.text.col = TRUE,
          conf.int = FALSE,
          xlab = "Time (Days)",
          ylab = "Survival Probability",
          title = "Kaplan-Meier Survival Curves by Treatment"
        )
        
        # Store both plot and table
        result$km_plot <- km_plot_obj$plot
        result$km_table <- km_plot_obj$table
        result$km_plot_obj <- km_plot_obj
        
        cat("KM plot generated successfully\n")
      }, error = function(e) {
        cat("ERROR generating KM plot:", e$message, "\n")
      })
    }
    
    rv$surv_results <- result
    rv$surv_running <- FALSE
    
    # Debug: Log the structure of the result
    cat("\n=== RESULT OBJECT STRUCTURE ===\n")
    cat("Class:", class(result), "\n")
    cat("Fields:", paste(names(result), collapse = ", "), "\n")
    cat("Has 'error' field:", "error" %in% names(result), "\n")
    cat("Has 'results' field:", "results" %in% names(result), "\n")
    cat("Has 'model' field:", "model" %in% names(result), "\n")
    cat("Has 'km_plot' field:", "km_plot" %in% names(result), "\n")
    if ("results" %in% names(result) && !is.null(result$results)) {
      cat("Results class:", class(result$results), "\n")
      cat("Results dimensions:", nrow(result$results), "rows x", ncol(result$results), "cols\n")
      cat("Results columns:", paste(colnames(result$results), collapse = ", "), "\n")
    }
    cat("====================================\n\n")
    
    if (isTRUE(result$error)) {
      output$surv_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error:", result$message)
      })
    } else {
      output$surv_status_msg <- renderUI({
        div(class = "alert alert-success",
            icon("check"), "Survival analysis completed successfully!")
      })
    }
  })
  
  # Summary output
  output$surv_summary <- renderText({
    req(rv$surv_results)
    cat("\n[surv_summary] Output function called\n")
    cat("[surv_summary] Result class:", class(rv$surv_results), "\n")
    cat("[surv_summary] Result fields:", paste(names(rv$surv_results), collapse = ", "), "\n")
    if (isTRUE(rv$surv_results$error)) {
      return(paste("Error:", rv$surv_results$message))
    }
    
    summary_text <- capture.output({
      cat("=== SURVIVAL ANALYSIS SUMMARY ===\n\n")
      
      # Show time specification method
      if (rv$surv_results$time_method == "from_date") {
        cat("Time Method: Calculated from Date Column\n")
        cat("Date Column:", input$col_date_for_day, "\n")
        cat("Reference Date:", input$reference_date, "\n")
      } else {
        cat("Time Method: Direct Day Column\n")
        cat("Day Column:", input$col_day, "\n")
      }
      
      cat("Event Column:", input$col_survival_event, "\n")
      cat("Treatment Column:", input$col_treatment, "\n")
      cat("Statistical Method:", rv$surv_results$method_used, "\n")
      
      if (!is.null(rv$surv_results$results)) {
        res <- rv$surv_results$results
        cat("Reference Group:", rv$surv_results$reference_group, "\n\n")
        
        cat("--- Median Survival ---\n")
        for (i in seq_len(nrow(res))) {
          med <- if (is.na(res$Median_Survival[i])) "Not reached" else paste0(res$Median_Survival[i], " days")
          note <- if (!is.na(res$Note[i]) && res$Note[i] != "") paste0(" [", res$Note[i], "]") else ""
          cat(sprintf("  %-20s %s%s\n", res$Group[i], med, note))
        }
      }
    })
    
    paste(summary_text, collapse = "\n")
  })
  
  # Kaplan-Meier plot
  output$surv_plot_km <- renderPlot({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) return(NULL)
    
    if (!is.null(rv$surv_results$km_plot)) {
      return(rv$surv_results$km_plot)
    } else {
      plot(1, 1, type="n", main="KM plot not available", xlab="", ylab="")
    }
  })
  
  # Cox model summary
  output$surv_cox_summary <- renderText({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$model)) {
      return("Model not available")
    }
    
    paste(capture.output(summary(rv$surv_results$model)), collapse = "\n")
  })
  
  # Results table
  output$surv_cox_table <- renderDT({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$results)) {
      return(datatable(data.frame(Message = "Results not available")))
    }
    
    results_df <- rv$surv_results$results
    numeric_cols <- sapply(results_df, is.numeric)
    results_df[numeric_cols] <- round(results_df[numeric_cols], 4)
    
    datatable(results_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Risk table output - numbers at risk per time point per group
  output$surv_risk_table <- renderDT({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$km_fit)) {
      return(datatable(data.frame(Message = "Risk table not available")))
    }
    
    km_sum <- summary(rv$surv_results$km_fit)
    risk_df <- data.frame(
      Time      = km_sum$time,
      Group     = sub("^Treatment=", "", as.character(km_sum$strata)),
      N_Risk    = km_sum$n.risk,
      N_Events  = km_sum$n.event,
      N_Censored = km_sum$n.censor,
      Survival  = round(km_sum$surv, 4)
    )
    # Reorder columns for readability
    risk_df <- risk_df[, c("Group", "Time", "N_Risk", "N_Events", "N_Censored", "Survival")]
    risk_df <- risk_df[order(risk_df$Group, risk_df$Time), ]
    
    datatable(risk_df,
              options = list(scrollX = TRUE, pageLength = 15),
              rownames = FALSE)
  })
  
  # Log-rank test (method information)
  # Download handlers
  output$download_surv_plot_km <- downloadHandler(
    filename = "kaplan_meier_plot.png",
    content = function(file) {
      if (!is.null(rv$surv_results$km_plot)) {
        png(file, width = 800, height = 600)
        print(rv$surv_results$km_plot)
        dev.off()
      }
    }
  )
  
  output$download_surv_cox <- downloadHandler(
    filename = "survival_results.csv",
    content = function(file) {
      if (!is.null(rv$surv_results$results)) {
        write.csv(rv$surv_results$results, file, row.names = FALSE)
      }
    }
  )
  
  output$download_surv_risk_table <- downloadHandler(
    filename = "survival_risk_table.csv",
    content = function(file) {
      req(rv$surv_results$km_fit)
      km_sum <- summary(rv$surv_results$km_fit)
      risk_df <- data.frame(
        Group     = sub("^Treatment=", "", as.character(km_sum$strata)),
        Time      = km_sum$time,
        N_Risk    = km_sum$n.risk,
        N_Events  = km_sum$n.event,
        N_Censored = km_sum$n.censor,
        Survival  = round(km_sum$surv, 4)
      )
      write.csv(risk_df, file, row.names = FALSE)
    }
  )
  
}
