#' Dose Response Analysis UI
#' @keywords internal
doseResponseUI <- function() {
  fluidPage(
    h2("Dose-Response Analysis"),

    p("Analyze tumor volume responses across dose levels using linear regression, ",
      "ANOVA, and non-linear dose-response models. ",
      "Upload and configure your data (including a Dose column) in the Data Upload tab first, ",
      "then run Tumor Growth Analysis before proceeding here."),

    fluidRow(
      column(3,
        wellPanel(
          h4("Analysis Parameters"),

          uiOutput("dr_dose_status"),

          conditionalPanel(
            condition = "output.dr_dose_available == true",

            selectInput("dr_control_group",
                       "Control Group Name",
                       choices = "",
                       selected = ""),

            uiOutput("dr_time_point_ui"),

            br(),
            actionButton("run_dr_analysis",
                        "Run Dose-Response Analysis",
                        class = "btn-lg btn-primary btn-block")
          ),

          br(),
          uiOutput("dr_status_msg")
        )
      ),

      column(9,
        tabsetPanel(
          id = "dr_results_tabs",

          tabPanel("Summary",
            tab_info(
              strong("Dose-response summary: "),
              "Results of linear regression (slope, R\u00b2, p-value) and one-way ANOVA across dose levels. ",
              "A significant linear slope indicates a systematic increase or decrease in tumor volume with dose. ",
              "A significant ANOVA F-test indicates that at least one dose group differs from the others."
            ),
            shinycssloaders::withSpinner(
              verbatimTextOutput("dr_summary")
            )
          ),

          tabPanel("Dose-Response Table",
            tab_info(
              strong("Per-dose summary statistics: "),
              "Mean, SD, and 95\u00a0% CI of tumor volume at each dose level at the selected time point. ",
              "Use this to quickly compare mean responses across doses."
            ),
            shinycssloaders::withSpinner(
              DTOutput("dr_summary_table")
            ),
            downloadButton("download_dr_table", "Download CSV")
          ),

          tabPanel("Scatter Plot",
            tab_info(
              strong("Dose vs tumor volume scatter plot: "),
              "Each point is one animal. The blue line is the linear regression fit. ",
              "A steep slope with a narrow confidence band suggests a strong dose-response relationship."
            ),
            shinycssloaders::withSpinner(
              plotOutput("dr_plot_scatter", height = "500px")
            ),
            downloadButton("download_dr_scatter", "Download Plot")
          ),

          tabPanel("Boxplot by Dose",
            tab_info(
              strong("Tumor volume distribution by dose level: "),
              "Each box shows the median, IQR, and individual data points per dose group. ",
              "Useful for assessing variance across dose levels and identifying outliers."
            ),
            shinycssloaders::withSpinner(
              plotOutput("dr_plot_boxplot", height = "500px")
            ),
            downloadButton("download_dr_boxplot", "Download Plot")
          ),

          tabPanel("Mean by Dose",
            tab_info(
              strong("Mean tumor volume with 95\u00a0% confidence intervals: "),
              "Bar heights represent group means; error bars are 95\u00a0% CIs. ",
              "Non-overlapping error bars between adjacent dose groups are a rough indicator of significance."
            ),
            shinycssloaders::withSpinner(
              plotOutput("dr_plot_barplot", height = "500px")
            ),
            downloadButton("download_dr_barplot", "Download Plot")
          ),

          tabPanel("Statistical Tests",
            tab_info(
              strong("Full statistical test results: "),
              "Includes linear regression coefficients (slope, intercept, p-value, R\u00b2), ",
              "one-way ANOVA table, and Tukey HSD post-hoc comparisons between dose groups (shown when ANOVA p\u00a0<\u00a00.05). ",
              "Polynomial trend p-values indicate whether a quadratic or higher-order relationship fits better than a line. ",
              "Non-linear model AIC/BIC are shown when the drc package fit converges."
            ),
            shinycssloaders::withSpinner(
              verbatimTextOutput("dr_stat_tests")
            ),
            shinycssloaders::withSpinner(
              DTOutput("dr_tukey_table")
            ),
            downloadButton("download_dr_stats", "Download Results CSV")
          )
        )
      )
    )
  )
}

#' Dose Response Server
#' @keywords internal
dose_response_server <- function(input, output, session, rv) {

  rv$dr_results <- NULL

  # Reactive: is a dose column configured?
  dr_dose_available <- reactive({
    !is.null(rv$dose_col) && rv$dose_col != ""
  })

  # Export this to JS for conditionalPanel
  output$dr_dose_available <- reactive({ dr_dose_available() })
  outputOptions(output, "dr_dose_available", suspendWhenHidden = FALSE)

  # Status / warning banner
  output$dr_dose_status <- renderUI({
    if (!rv$data_valid) {
      return(div(class = "alert alert-warning", icon("exclamation-triangle"),
                 "Please upload and configure data in the Data Upload tab first."))
    }
    if (is.null(rv$processed_data)) {
      return(div(class = "alert alert-warning", icon("exclamation-triangle"),
                 "Please run Tumor Growth Analysis first to resolve Volume and Day columns."))
    }
    if (!dr_dose_available()) {
      return(div(class = "alert alert-warning", icon("exclamation-triangle"),
                 "No Dose column selected. Please configure it in the Data Upload tab (Optional Columns \u2192 Dose Level)."))
    }
    div(class = "alert alert-success", icon("check"),
        paste("Dose column:", rv$dose_col))
  })

  # Populate control group choices from treatment column
  observe({
    req(rv$processed_data, rv$tx_col)
    groups <- sort(unique(as.character(na.omit(rv$processed_data[[rv$tx_col]]))))
    default_ctrl <- pick_reference_group(groups)
    updateSelectInput(session, "dr_control_group",
                      choices = groups,
                      selected = default_ctrl)
  })

  # Dynamic time-point selector: "All dates" (default) + individual Day values
  output$dr_time_point_ui <- renderUI({
    req(rv$processed_data, rv$time_col)
    days <- sort(unique(na.omit(rv$processed_data[[rv$time_col]])))
    day_choices <- c("All dates" = "all",
                     setNames(as.character(days), paste("Day", days)))
    selectInput("dr_time_point",
                "Time Point for Analysis",
                choices  = day_choices,
                selected = "all")
  })

  observeEvent(input$run_dr_analysis, {
    req(rv$processed_data, dr_dose_available())

    output$dr_status_msg <- renderUI({
      div(class = "alert alert-info", icon("spinner", class = "fa-spin"),
          "Running dose-response analysis...")
    })

    data     <- rv$processed_data
    dose_col <- rv$dose_col
    vol_col  <- rv$vol_col
    time_col <- rv$time_col
    tx_col   <- rv$tx_col
    id_col   <- rv$id_col

    # Coerce dose column to numeric
    data[[dose_col]] <- suppressWarnings(as.numeric(data[[dose_col]]))
    n_na_dose <- sum(is.na(data[[dose_col]]))
    if (n_na_dose == nrow(data)) {
      output$dr_status_msg <- renderUI({
        div(class = "alert alert-danger", icon("exclamation-circle"),
            paste("Dose column", dose_col, "contains no numeric values. Please check your data."))
      })
      return()
    }

    time_point_val <- if (!is.null(input$dr_time_point) &&
                           input$dr_time_point != "all" &&
                           nzchar(input$dr_time_point)) {
      as.numeric(input$dr_time_point)
    } else {
      NULL
    }

    result <- safe_analysis({
      mouseExperiment::dose_response_statistics(
        df                 = data,
        dose_column        = dose_col,
        treatment_column   = tx_col,
        volume_column      = vol_col,
        day_column         = time_col,
        id_column          = id_col,
        time_point         = time_point_val,
        control_group_name = input$dr_control_group
      )
    })

    rv$dr_results <- result

    if (isTRUE(result$error)) {
      output$dr_status_msg <- renderUI({
        div(class = "alert alert-danger", icon("exclamation-circle"),
            "Error: ", result$message)
      })
    } else {
      output$dr_status_msg <- renderUI({
        div(class = "alert alert-success", icon("check"),
            "Dose-response analysis completed successfully!")
      })
    }
  })

  # ---- Summary ----
  output$dr_summary <- renderText({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(paste("Error:", rv$dr_results$message))

    res   <- rv$dr_results
    stats <- res$statistics

    lines <- c(
      "=== DOSE-RESPONSE ANALYSIS SUMMARY ===",
      "",
      "-- Linear Regression --",
      sprintf("  Slope:   %.4f", stats$linear_slope %||% NA),
      sprintf("  R\u00b2:      %.4f", stats$linear_r_squared %||% NA),
      sprintf("  p-value: %.4g", stats$linear_p_value %||% NA),
      "",
      "-- ANOVA (Dose as Factor) --",
      sprintf("  p-value: %.4g", stats$anova_p_value %||% NA),
      ""
    )

    if (!is.null(stats$ec50)) {
      lines <- c(lines,
        "-- Non-linear Dose-Response Model --",
        sprintf("  EC50:         %.4g", stats$ec50),
        sprintf("  Hill Slope:   %.4f", stats$hill_slope),
        sprintf("  Lower Limit:  %.4f", stats$lower_limit),
        sprintf("  Upper Limit:  %.4f", stats$upper_limit),
        sprintf("  Linear AIC:   %.2f  |  Non-linear AIC: %.2f",
                stats$linear_aic %||% NA, stats$nonlinear_aic %||% NA),
        ""
      )
    }

    if (!is.null(stats$linear_trend_pvalue)) {
      lines <- c(lines,
        "-- Polynomial Trends --",
        sprintf("  Linear trend p:    %.4g", stats$linear_trend_pvalue),
        sprintf("  Quadratic trend p: %.4g", stats$quadratic_trend_pvalue %||% NA),
        ""
      )
    }

    paste(lines, collapse = "\n")
  })

  # ---- Dose-Response Table ----
  output$dr_summary_table <- renderDT({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(datatable(data.frame(Message = "Results not available")))
    tbl <- as.data.frame(rv$dr_results$summary_table)
    tbl <- format_df_for_display(tbl)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
  })

  # ---- Plots ----
  output$dr_plot_scatter <- renderPlot({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error) || is.null(rv$dr_results$plots$scatter)) return(NULL)
    rv$dr_results$plots$scatter + ggplot2::theme_classic()
  })

  output$dr_plot_boxplot <- renderPlot({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error) || is.null(rv$dr_results$plots$boxplot)) return(NULL)
    rv$dr_results$plots$boxplot + ggplot2::theme_classic()
  })

  output$dr_plot_barplot <- renderPlot({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error) || is.null(rv$dr_results$plots$barplot)) return(NULL)
    rv$dr_results$plots$barplot + ggplot2::theme_classic()
  })

  # ---- Statistical Tests ----
  output$dr_stat_tests <- renderText({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(paste("Error:", rv$dr_results$message))
    stats <- rv$dr_results$statistics

    lines <- capture.output({
      cat("=== LINEAR REGRESSION ===\n")
      if (!is.null(rv$dr_results$linear_model)) print(summary(rv$dr_results$linear_model))
      cat("\n=== ANOVA TABLE ===\n")
      if (!is.null(rv$dr_results$anova_model)) print(summary(rv$dr_results$anova_model))
      if (!is.null(stats$growth_dose_p_value)) {
        cat("\n=== GROWTH RATE vs DOSE ===\n")
        cat(sprintf("  p-value: %.4g,  R\u00b2: %.4f\n",
                    stats$growth_dose_p_value, stats$growth_dose_r_squared))
      }
    })
    paste(lines, collapse = "\n")
  })

  output$dr_tukey_table <- renderDT({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(datatable(data.frame(Message = "No results")))
    tukey <- rv$dr_results$statistics$tukey_results
    if (is.null(tukey)) return(datatable(data.frame(
      Message = "Tukey HSD not computed (ANOVA p \u2265 0.05 or not available)")))
    # TukeyHSD returns a list; extract the factor component
    tbl <- as.data.frame(tukey[[1]])
    tbl$Comparison <- rownames(tbl)
    tbl <- tbl[, c("Comparison", setdiff(names(tbl), "Comparison"))]
    tbl <- format_df_for_display(tbl)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  # ---- Download handlers ----
  output$download_dr_table <- downloadHandler(
    filename = "dose_response_summary.csv",
    content  = function(file) {
      req(rv$dr_results, !isTRUE(rv$dr_results$error))
      write.csv(as.data.frame(rv$dr_results$summary_table), file, row.names = FALSE)
    }
  )

  output$download_dr_scatter <- downloadHandler(
    filename = "dose_scatter.png",
    content  = function(file) {
      req(rv$dr_results$plots$scatter)
      png(file, width = 800, height = 600)
      print(rv$dr_results$plots$scatter)
      dev.off()
    }
  )

  output$download_dr_boxplot <- downloadHandler(
    filename = "dose_boxplot.png",
    content  = function(file) {
      req(rv$dr_results$plots$boxplot)
      png(file, width = 800, height = 600)
      print(rv$dr_results$plots$boxplot)
      dev.off()
    }
  )

  output$download_dr_barplot <- downloadHandler(
    filename = "dose_barplot.png",
    content  = function(file) {
      req(rv$dr_results$plots$barplot)
      png(file, width = 800, height = 600)
      print(rv$dr_results$plots$barplot)
      dev.off()
    }
  )

  output$download_dr_stats <- downloadHandler(
    filename = "dose_response_stats.csv",
    content  = function(file) {
      req(rv$dr_results, !isTRUE(rv$dr_results$error))
      stats <- rv$dr_results$statistics
      out <- data.frame(
        Metric = c("Linear Slope", "Linear R2", "Linear P-value",
                   "ANOVA P-value",
                   "EC50", "Hill Slope", "Lower Limit", "Upper Limit",
                   "Linear AIC", "Nonlinear AIC",
                   "Growth Rate P-value", "Growth Rate R2"),
        Value  = c(
          stats$linear_slope %||% NA,
          stats$linear_r_squared %||% NA,
          stats$linear_p_value %||% NA,
          stats$anova_p_value %||% NA,
          stats$ec50 %||% NA,
          stats$hill_slope %||% NA,
          stats$lower_limit %||% NA,
          stats$upper_limit %||% NA,
          stats$linear_aic %||% NA,
          stats$nonlinear_aic %||% NA,
          stats$growth_dose_p_value %||% NA,
          stats$growth_dose_r_squared %||% NA
        )
      )
      write.csv(out, file, row.names = FALSE)
    }
  )
}

# Null-coalescing helper for safe field access
`%||%` <- function(a, b) if (!is.null(a)) a else b
