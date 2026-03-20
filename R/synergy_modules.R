#' Drug Synergy Analysis UI
#' @keywords internal
synergyUI <- function() {
  fluidPage(
    h2("Drug Synergy Analysis"),
    
    p("Analyze drug interactions and synergy effects using combination index and Bliss model. ",
      "Requires data with multiple drug treatments and dose information."),
    
    fluidRow(
      column(3,
        wellPanel(
          h4("Analysis Parameters"),
          
          selectInput("syn_drug1_column",
                     "Drug 1 Column",
                     choices = c("-- Select column --" = "")),
          
          selectInput("syn_drug2_column",
                     "Drug 2 Column",
                     choices = c("-- Select column --" = "")),
          
          selectInput("syn_dose_column",
                     "Dose Column",
                     choices = c("-- Select column --" = "")),
          
          selectInput("syn_response_column",
                     "Response/Effect Column",
                     choices = c("-- Select column --" = "")),
          
          selectInput("syn_model",
                     "Synergy Model",
                     choices = c("Combination Index" = "combindex",
                                "Bliss Independence" = "bliss",
                                "Loewe Additivity" = "loewe"),
                     selected = "bliss"),
          
          checkboxInput("syn_include_3d",
                       "Generate 3D Surface Plot",
                       value = TRUE),
          
          sliderInput("syn_confidence",
                     "Confidence Level",
                     min = 0.80, max = 0.99, value = 0.95, step = 0.01),
          
          br(),
          
          actionButton("run_syn_analysis", 
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),
          
          uiOutput("syn_status_msg")
        )
      ),
      
      column(9,
        tabsetPanel(
          id = "syn_results_tabs",
          
          tabPanel("Summary",
            shinycssloaders::withSpinner(
              verbatimTextOutput("syn_summary")
            )
          ),
          
          tabPanel("Combination Index",
            shinycssloaders::withSpinner(
              DTOutput("syn_ci_table")
            ),
            downloadButton("download_syn_ci", "Download Results")
          ),
          
          tabPanel("Bliss Plot",
            shinycssloaders::withSpinner(
              plotOutput("syn_plot_bliss", height = "500px")
            ),
            downloadButton("download_syn_plot_bliss", "Download Plot")
          ),
          
          tabPanel("Synergy Map",
            shinycssloaders::withSpinner(
              plotOutput("syn_plot_synergy_map", height = "500px")
            ),
            downloadButton("download_syn_plot_synergy_map", "Download Plot")
          ),
          
          tabPanel("3D Surface",
            shinycssloaders::withSpinner(
              conditionalPanel(
                condition = "input.syn_include_3d",
                plotOutput("syn_plot_3d", height = "600px")
              )
            ),
            downloadButton("download_syn_plot_3d", "Download Plot")
          ),
          
          tabPanel("Dose Response Curves",
            shinycssloaders::withSpinner(
              plotOutput("syn_plot_doseresponse", height = "500px")
            ),
            downloadButton("download_syn_plot_doseresponse", "Download Plot")
          )
        )
      )
    )
  )
}

#' Drug Synergy Analysis Server Logic
#' @keywords internal
synergy_server <- function(input, output, session, rv) {
  
  rv$syn_results <- NULL
  rv$syn_running <- FALSE
  
  # Update column choices when data loads
  observe({
    if (!is.null(rv$uploaded_data) && rv$data_valid) {
      col_names <- colnames(rv$uploaded_data)
      updateSelectInput(session, "syn_drug1_column", choices = c("-- Select column --" = "", col_names))
      updateSelectInput(session, "syn_drug2_column", choices = c("-- Select column --" = "", col_names))
      updateSelectInput(session, "syn_dose_column", choices = c("-- Select column --" = "", col_names))
      updateSelectInput(session, "syn_response_column", choices = c("-- Select column --" = "", col_names))
    }
  })
  
  # Run synergy analysis
  observeEvent(input$run_syn_analysis, {
    req(rv$uploaded_data, rv$data_valid)
    
    # Validate selections
    if (input$syn_drug1_column == "" || input$syn_drug2_column == "" || 
        input$syn_dose_column == "" || input$syn_response_column == "") {
      output$syn_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error: Please select all required columns")
      })
      return()
    }
    
    rv$syn_running <- TRUE
    output$syn_status_msg <- renderUI({
      div(class = "alert alert-info",
          icon("spinner", class = "fa-spin"), "Running synergy analysis...")
    })
    
    # Run analysis
    result <- safe_analysis({
      mouseExperiment::analyze_drug_synergy(
        df = rv$uploaded_data,
        drug1_column = input$syn_drug1_column,
        drug2_column = input$syn_drug2_column,
        dose_column = input$syn_dose_column,
        response_column = input$syn_response_column,
        model = input$syn_model,
        confidence_level = input$syn_confidence
      )
    })
    
    # Generate Bliss plot
    if (!isTRUE(result$error) && !is.null(result$bliss_data)) {
      tryCatch({
        result$plot_bliss <- mouseExperiment::plot_bliss(result$bliss_data)
      }, error = function(e) {
        cat("WARNING: Error generating Bliss plot:", e$message, "\n")
      })
    }
    
    # Generate combination index plot
    if (!isTRUE(result$error) && !is.null(result$ci_data)) {
      tryCatch({
        result$plot_ci <- mouseExperiment::plot_combination_index(result$ci_data)
      }, error = function(e) {
        cat("WARNING: Error generating CI plot:", e$message, "\n")
      })
    }
    
    # Generate synergy combined visualization
    if (!isTRUE(result$error)) {
      tryCatch({
        result$plot_synergy_combined <- mouseExperiment::plot_synergy_combined(result)
      }, error = function(e) {
        cat("WARNING: Error generating synergy plot:", e$message, "\n")
      })
    }
    
    rv$syn_results <- result
    rv$syn_running <- FALSE
    
    if (isTRUE(result$error)) {
      output$syn_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error:", result$message)
      })
    } else {
      output$syn_status_msg <- renderUI({
        div(class = "alert alert-success",
            icon("check"), "Synergy analysis completed successfully!")
      })
    }
  })
  
  # Summary output
  output$syn_summary <- renderText({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) {
      return(paste("Error:", rv$syn_results$message))
    }
    
    summary_text <- capture.output({
      cat("=== DRUG SYNERGY ANALYSIS SUMMARY ===\n\n")
      cat("Drug 1:", input$syn_drug1_column, "\n")
      cat("Drug 2:", input$syn_drug2_column, "\n")
      cat("Dose Column:", input$syn_dose_column, "\n")
      cat("Response Column:", input$syn_response_column, "\n\n")
      cat("Model:", input$syn_model, "\n")
      cat("Confidence Level:", round(input$syn_confidence * 100), "%\n\n")
      
      if (!is.null(rv$syn_results$synergy_assessment)) {
        cat("Overall Assessment:\n")
        print(rv$syn_results$synergy_assessment)
      }
    })
    
    paste(summary_text, collapse = "\n")
  })
  
  # Combination Index table
  output$syn_ci_table <- renderDT({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || is.null(rv$syn_results$ci_data)) {
      return(datatable(data.frame(Message = "CI data not available")))
    }
    
    ci_df <- as.data.frame(rv$syn_results$ci_data)
    numeric_cols <- sapply(ci_df, is.numeric)
    ci_df[numeric_cols] <- round(ci_df[numeric_cols], 4)
    
    datatable(ci_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Bliss plot
  output$syn_plot_bliss <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    
    if (!is.null(rv$syn_results$plot_bliss)) {
      return(rv$syn_results$plot_bliss)
    } else {
      plot(1, 1, type="n", main="Bliss plot not available", xlab="", ylab="")
    }
  })
  
  # Synergy map/heatmap
  output$syn_plot_synergy_map <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    
    if (!is.null(rv$syn_results$plot_synergy_combined)) {
      return(rv$syn_results$plot_synergy_combined)
    } else {
      plot(1, 1, type="n", main="Synergy map not available", xlab="", ylab="")
    }
  })
  
  # 3D surface plot
  output$syn_plot_3d <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    
    if (!is.null(rv$syn_results$plot_3d)) {
      return(rv$syn_results$plot_3d)
    } else if (!input$syn_include_3d) {
      plot(1, 1, type="n", main="3D plot disabled", xlab="", ylab="")
    } else {
      plot(1, 1, type="n", main="3D plot not available", xlab="", ylab="")
    }
  })
  
  # Dose-response plot
  output$syn_plot_doseresponse <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    
    if (!is.null(rv$syn_results$plot_doseresponse)) {
      return(rv$syn_results$plot_doseresponse)
    } else {
      plot(1, 1, type="n", main="Dose-response plot not available", xlab="", ylab="")
    }
  })
  
  # Download handlers
  output$download_syn_ci <- downloadHandler(
    filename = "combination_index_results.csv",
    content = function(file) {
      if (!is.null(rv$syn_results$ci_data)) {
        write.csv(as.data.frame(rv$syn_results$ci_data), file, row.names = FALSE)
      }
    }
  )
  
  output$download_syn_plot_bliss <- downloadHandler(
    filename = "bliss_plot.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_bliss)) {
        png(file, width = 800, height = 600)
        print(rv$syn_results$plot_bliss)
        dev.off()
      }
    }
  )
  
  output$download_syn_plot_synergy_map <- downloadHandler(
    filename = "synergy_map.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_synergy_combined)) {
        png(file, width = 800, height = 600)
        print(rv$syn_results$plot_synergy_combined)
        dev.off()
      }
    }
  )
  
  output$download_syn_plot_3d <- downloadHandler(
    filename = "synergy_3d_surface.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_3d)) {
        png(file, width = 900, height = 600)
        print(rv$syn_results$plot_3d)
        dev.off()
      }
    }
  )
  
  output$download_syn_plot_doseresponse <- downloadHandler(
    filename = "dose_response_curves.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_doseresponse)) {
        png(file, width = 800, height = 600)
        print(rv$syn_results$plot_doseresponse)
        dev.off()
      }
    }
  )
}
