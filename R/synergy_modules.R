#' Drug Synergy Analysis UI
#' @keywords internal
synergyUI <- function() {
  fluidPage(
    h2("Drug Synergy Analysis"),
    
    p("Test whether a drug combination shows synergistic, additive, or antagonistic effects ",
      "using Bliss independence and Loewe additivity models. ",
      "Ensure tumor growth data is uploaded and configured in the Data Upload tab."),
    
    fluidRow(
      column(3,
        wellPanel(
          h4("Treatment Groups"),
          
          p(tags$small("Select which treatment group in your data corresponds to each role:")),
          
          uiOutput("syn_group_selectors_ui"),
          
          hr(),
          
          h4("Analysis Options"),
          
          numericInput("syn_eval_timepoint",
                      "Evaluation Timepoint (leave blank for last)",
                      value = NA, min = 0, step = 1),
          
          checkboxInput("syn_over_time",
                       "Analyze synergy across all timepoints",
                       value = FALSE),
          
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
            tab_info(strong("Synergy summary: "), "Overall assessment combining Bliss independence and Combination Index (CI) metrics. Reports whether the drug combination is synergistic, additive, or antagonistic at the evaluated timepoint, and lists the key quantitative scores."),
            shinycssloaders::withSpinner(
              verbatimTextOutput("syn_summary")
            )
          ),
          
          tabPanel("Synergy Metrics",
            tab_info(strong("Quantitative synergy scores. "),
              strong("Bliss excess: "), "observed\u202f\u2212\u202fexpected effect (positive\u202f=\u202fsynergy, negative\u202f=\u202fantagonism). ",
              strong("Combination Index (CI): "), "CI\u202f<\u202f1 synergistic, CI\u202f=\u202f1 additive, CI\u202f>\u202f1 antagonistic. ",
              "When \u2018Analyze over time\u2019 is enabled, this table shows metrics at every measured timepoint."),
            shinycssloaders::withSpinner(
              DTOutput("syn_metrics_table")
            ),
            downloadButton("download_syn_metrics", "Download Results")
          ),
          
          tabPanel("TGI Table",
            tab_info(strong("Tumor Growth Inhibition (TGI): "), "percentage reduction in tumor volume relative to the control group at the evaluation timepoint. Higher TGI\u202f% indicates greater anti-tumor activity. The table also shows mean volumes and standard errors per group."),
            shinycssloaders::withSpinner(
              DTOutput("syn_tgi_table")
            ),
            downloadButton("download_syn_tgi", "Download TGI")
          ),
          
          tabPanel("Synergy Plot",
            tab_info(strong("Visual synergy comparison: "), "bars showing observed combination effect alongside the expected effect if the two drugs acted independently (Bliss model). A combination bar taller than the expected bar indicates synergy at the evaluation timepoint."),
            shinycssloaders::withSpinner(
              plotOutput("syn_plot_main", height = "500px")
            ),
            downloadButton("download_syn_plot_main", "Download Plot")
          ),
          
          tabPanel("Trend Over Time",
            tab_info(strong("Synergy trend: "), "synergy/antagonism scores plotted across all measured timepoints. Shows whether the interaction between the two drugs strengthens, weakens, or reverses over the course of the experiment. Requires \u2018Analyze synergy across all timepoints\u2019 to be enabled."),
            conditionalPanel(
              condition = "input.syn_over_time",
              shinycssloaders::withSpinner(
                plotOutput("syn_plot_trend", height = "500px")
              ),
              downloadButton("download_syn_plot_trend", "Download Plot")
            ),
            conditionalPanel(
              condition = "!input.syn_over_time",
              br(),
              p("Enable 'Analyze synergy across all timepoints' to see trend plots.")
            )
          ),
          
          tabPanel("Combination Index",
            tab_info(strong("Combination Index (CI) over time: "), "CI values at each measured timepoint. The dashed reference line at CI\u202f=\u202f1 represents additivity. Points below 1 indicate synergism at that timepoint; points above 1 indicate antagonism. Requires \u2018Analyze synergy across all timepoints\u2019 to be enabled."),
            conditionalPanel(
              condition = "input.syn_over_time",
              shinycssloaders::withSpinner(
                plotOutput("syn_plot_ci", height = "500px")
              ),
              downloadButton("download_syn_plot_ci", "Download Plot")
            ),
            conditionalPanel(
              condition = "!input.syn_over_time",
              br(),
              p("Enable 'Analyze synergy across all timepoints' to see combination index over time.")
            )
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
  rv$syn_over_time <- FALSE
  
  # Dynamically render group selector dropdowns from treatment column
  output$syn_group_selectors_ui <- renderUI({
    req(rv$uploaded_data, rv$data_valid)
    
    treatment_col <- input$col_treatment
    if (is.null(treatment_col) || treatment_col == "") {
      return(p(class = "text-warning",
               icon("exclamation-triangle"),
               "Select Treatment column in Data Upload first."))
    }
    
    groups <- sort(unique(as.character(na.omit(rv$uploaded_data[[treatment_col]]))))
    choices <- setNames(groups, groups)
    
    tagList(
      selectInput("syn_drug_a", "Drug A (single agent)",
                 choices = choices,
                 selected = if (length(groups) >= 2) groups[2] else groups[1]),
      
      selectInput("syn_drug_b", "Drug B (single agent)",
                 choices = choices,
                 selected = if (length(groups) >= 3) groups[3] else groups[1]),
      
      selectInput("syn_combo", "Combination Treatment",
                 choices = choices,
                 selected = if (length(groups) >= 4) groups[4] else groups[1]),
      
      selectInput("syn_control", "Control Group",
                 choices = choices,
                 selected = pick_reference_group(groups))
    )
  })
  
  # Run synergy analysis
  observeEvent(input$run_syn_analysis, {
    req(rv$uploaded_data, rv$data_valid)
    
    needed <- c(input$syn_drug_a, input$syn_drug_b, input$syn_combo, input$syn_control)
    if (any(is.null(needed)) || any(needed == "")) {
      output$syn_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"),
            "Please select all four treatment groups.")
      })
      return()
    }
    
    if (length(unique(needed)) < 4) {
      output$syn_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"),
            "Each role (Drug A, Drug B, Combination, Control) must be a different group.")
      })
      return()
    }
    
    output$syn_status_msg <- renderUI({
      div(class = "alert alert-info",
          icon("spinner", class = "fa-spin"), "Running synergy analysis...")
    })
    
    # Use preprocessed data (Volume + Day resolved) if TG analysis has already
    # been run; otherwise preprocess inline from the raw uploaded data.
    if (!is.null(rv$processed_data)) {
      analysis_data <- rv$processed_data
      vol_col  <- rv$vol_col
      time_col <- rv$time_col
      tx_col   <- rv$tx_col
    } else {
      # Inline preprocessing mirroring the tumor_growth_server logic
      analysis_data <- rv$uploaded_data
      tx_col   <- if (!is.null(input$col_treatment) && input$col_treatment != "") input$col_treatment else "Treatment"
      id_col   <- if (!is.null(input$col_id)        && input$col_id        != "") input$col_id        else "ID"
      
      if (!is.null(input$time_method) && input$time_method == "from_date" &&
          !is.null(input$col_date_for_day) && input$col_date_for_day != "") {
        first_val <- as.character(analysis_data[1, input$col_date_for_day])
        fmt <- if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_val)) "%m/%d/%Y" else
               if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", first_val)) "%Y-%m-%d" else
               if (grepl("^\\d{1,2}-[A-Za-z]{3}$", first_val)) "%d-%b" else NULL
        tryCatch({
          analysis_data <- mouseExperiment::calculate_dates(analysis_data,
            start_date = input$reference_date,
            date_column = input$col_date_for_day,
            date_format = fmt)
        }, error = function(e) NULL)
        time_col <- "Day"
      } else {
        time_col <- if (!is.null(input$col_day) && input$col_day != "") input$col_day else "Day"
      }
      
      if (!is.null(input$volume_method) && input$volume_method == "calculate" &&
          !is.null(input$col_length) && input$col_length != "" &&
          !is.null(input$col_width)  && input$col_width  != "") {
        tryCatch({
          analysis_data <- mouseExperiment::calculate_volume(analysis_data,
            length_column = input$col_length, width_column = input$col_width,
            formula = input$volume_formula)
        }, error = function(e) NULL)
        vol_col <- "Volume"
      } else {
        vol_col <- if (!is.null(input$col_volume) && input$col_volume != "") input$col_volume else "Volume"
      }
    }
    
    # Coerce time and volume to numeric (guard against CSV character reads)
    if (time_col %in% colnames(analysis_data))
      analysis_data[[time_col]] <- suppressWarnings(as.numeric(analysis_data[[time_col]]))
    if (vol_col %in% colnames(analysis_data))
      analysis_data[[vol_col]]  <- suppressWarnings(as.numeric(analysis_data[[vol_col]]))
    
    eval_tp  <- if (!is.null(input$syn_eval_timepoint) && !is.na(input$syn_eval_timepoint)) {
      input$syn_eval_timepoint
    } else {
      NULL
    }
    
    rv$syn_over_time <- isTRUE(input$syn_over_time)
    
    if (rv$syn_over_time) {
      result <- safe_analysis({
        mouseExperiment::analyze_drug_synergy_over_time(
          df               = analysis_data,
          treatment_column = tx_col,
          volume_column    = vol_col,
          time_column      = time_col,
          drug_a_name      = input$syn_drug_a,
          drug_b_name      = input$syn_drug_b,
          combo_name       = input$syn_combo,
          control_name     = input$syn_control
        )
      })
    } else {
      result <- safe_analysis({
        mouseExperiment::analyze_drug_synergy(
          df               = analysis_data,
          treatment_column = tx_col,
          volume_column    = vol_col,
          time_column      = time_col,
          drug_a_name      = input$syn_drug_a,
          drug_b_name      = input$syn_drug_b,
          combo_name       = input$syn_combo,
          control_name     = input$syn_control,
          eval_time_point  = eval_tp
        )
      })
    }
    
    if (!isTRUE(result$error)) {
      if (rv$syn_over_time) {
        tryCatch({
          syn_df <- result$synergy_summary

          # Bliss Independence plot — matches demo Fig 4 style
          result$plot_main <- ggplot2::ggplot(syn_df, ggplot2::aes(x = Time_Point)) +
            ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_A, color = "Drug A")) +
            ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_B, color = "Drug B")) +
            ggplot2::geom_line(ggplot2::aes(y = TGI_Combo, color = "Combination")) +
            ggplot2::geom_line(ggplot2::aes(y = Bliss_Expected_TGI, color = "Bliss Value"),
                               alpha = 0.5, linewidth = 3) +
            ggplot2::scale_color_manual(
              values = c("Drug A" = "blue", "Drug B" = "red",
                         "Combination" = "purple", "Bliss Value" = "grey")) +
            ggplot2::labs(x = "Days", y = "Percent Tumor Growth Inhibition", color = NULL) +
            ggplot2::ylim(0, 100) +
            ggplot2::theme_classic() +
            ggplot2::theme(legend.position = c(0.85, 0.2))

          # Combination Index plot — matches demo Fig 5 style with shaded regions
          mid_x <- mean(syn_df$Time_Point)
          result$plot_ci <- ggplot2::ggplot(syn_df,
                                            ggplot2::aes(x = Time_Point, y = Combination_Index)) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin = 0,    ymax = 0.85), fill = "blue",  alpha = 0.2) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin = 0.85, ymax = 1.15), fill = "grey",  alpha = 0.4) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin = 1.15, ymax = 2),    fill = "red",   alpha = 0.2) +
            ggplot2::geom_line(color = "purple", linewidth = 1.5) +
            ggplot2::annotate("text", x = mid_x, y = 0.20, label = "Synergistic",
                              color = "blue", size = 5) +
            ggplot2::annotate("text", x = mid_x, y = 1.05, label = "Additive",    size = 5) +
            ggplot2::annotate("text", x = mid_x, y = 1.80, label = "Antagonistic",
                              color = "red",  size = 5) +
            ggplot2::labs(x = "Days", y = "Combination Index", title = "") +
            ggplot2::ylim(0, 2) +
            ggplot2::theme_classic()

        }, error = function(e) cat("WARNING: Synergy over-time plot error:", e$message, "\n"))
      } else {
        tryCatch({
          result$plot_main  <- mouseExperiment::plot_drug_synergy(result)
        }, error = function(e) cat("WARNING: Synergy plot error:", e$message, "\n"))
      }
    }
    
    rv$syn_results <- result
    
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
  
  output$syn_summary <- renderText({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(paste("Error:", rv$syn_results$message))
    
    res <- rv$syn_results
    paste(capture.output({
      cat("=== DRUG SYNERGY ANALYSIS SUMMARY ===\n\n")
      cat("Drug A (single agent): ", input$syn_drug_a, "\n")
      cat("Drug B (single agent): ", input$syn_drug_b, "\n")
      cat("Combination:           ", input$syn_combo, "\n")
      cat("Control:               ", input$syn_control, "\n\n")
      
      if (isTRUE(rv$syn_over_time)) {
        cat("Analysis type: Over all timepoints\n\n")
        if (!is.null(res$peak_bliss_synergy)) cat("Peak Bliss synergy at timepoint:", res$peak_bliss_synergy$time_point, "\n")
        if (!is.null(res$peak_ci_synergy))    cat("Strongest CI synergy at timepoint:", res$peak_ci_synergy$time_point, "\n")
      } else {
        cat("Analysis type: Single timepoint\n")
        if (!is.null(res$evaluation_time_point)) cat("Evaluated at Day:", res$evaluation_time_point, "\n")
        cat("\nOverall Assessment:", res$overall_assessment, "\n\n")
        if (!is.null(res$combination_index)) {
          cat("Combination Index (CI):", round(res$combination_index$ci, 4), "\n")
          cat("Interpretation:", res$combination_index$interpretation, "\n")
          cat("  CI < 1 = Synergistic, CI = 1 = Additive, CI > 1 = Antagonistic\n\n")
        }
        if (!is.null(res$bliss_independence)) {
          bi <- res$bliss_independence
          cat("Bliss Independence:\n")
          cat("  Expected effect:", round(bi$expected_effect, 4), "\n")
          cat("  Observed effect:", round(bi$observed_effect, 4), "\n")
          cat("  Synergistic:    ", bi$synergy, "\n")
        }
      }
    }), collapse = "\n")
  })
  
  output$syn_metrics_table <- renderDT({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(datatable(data.frame(Message = "Results not available")))
    
    df <- if (isTRUE(rv$syn_over_time) && !is.null(rv$syn_results$synergy_summary)) {
      rv$syn_results$synergy_summary
    } else if (!is.null(rv$syn_results$synergy_metrics)) {
      as.data.frame(rv$syn_results$synergy_metrics)
    } else {
      data.frame(Message = "No metrics available")
    }
    df <- format_df_for_display(df)
    datatable(df, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
  })
  
  output$syn_tgi_table <- renderDT({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || is.null(rv$syn_results$summary)) {
      return(datatable(data.frame(Message = "TGI data not available")))
    }
    df <- as.data.frame(rv$syn_results$summary)
    df <- format_df_for_display(df)
    datatable(df, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
  })
  
  output$syn_plot_main <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    if (!is.null(rv$syn_results$plot_main)) print(rv$syn_results$plot_main)
    else plot(1, 1, type = "n", main = "Plot not available", xlab = "", ylab = "")
  })
  
  output$syn_plot_trend <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || !isTRUE(rv$syn_over_time)) return(NULL)
    if (!is.null(rv$syn_results$plot_trend)) print(rv$syn_results$plot_trend)
    else plot(1, 1, type = "n", main = "Trend plot not available", xlab = "", ylab = "")
  })
  
  output$syn_plot_ci <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || !isTRUE(rv$syn_over_time)) return(NULL)
    if (!is.null(rv$syn_results$plot_ci)) print(rv$syn_results$plot_ci)
    else plot(1, 1, type = "n", main = "CI plot not available", xlab = "", ylab = "")
  })
  
  output$download_syn_metrics <- downloadHandler(
    filename = "synergy_metrics.csv",
    content = function(file) {
      df <- if (isTRUE(rv$syn_over_time) && !is.null(rv$syn_results$synergy_summary)) {
        rv$syn_results$synergy_summary
      } else if (!is.null(rv$syn_results$synergy_metrics)) {
        as.data.frame(rv$syn_results$synergy_metrics)
      } else NULL
      if (!is.null(df)) write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$download_syn_tgi <- downloadHandler(
    filename = "synergy_tgi.csv",
    content = function(file) {
      if (!is.null(rv$syn_results$summary))
        write.csv(as.data.frame(rv$syn_results$summary), file, row.names = FALSE)
    }
  )
  
  output$download_syn_plot_main <- downloadHandler(
    filename = "synergy_plot.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_main)) {
        png(file, width = 900, height = 600); print(rv$syn_results$plot_main); dev.off()
      }
    }
  )
  
  output$download_syn_plot_trend <- downloadHandler(
    filename = "synergy_trend.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_trend)) {
        png(file, width = 900, height = 600); print(rv$syn_results$plot_trend); dev.off()
      }
    }
  )
  
  output$download_syn_plot_ci <- downloadHandler(
    filename = "combination_index.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_ci)) {
        png(file, width = 900, height = 600); print(rv$syn_results$plot_ci); dev.off()
      }
    }
  )
}
