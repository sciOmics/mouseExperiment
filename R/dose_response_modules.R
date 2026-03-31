#' Dose Response Analysis UI
#' @keywords internal
doseResponseUI <- function() {
  fluidPage(
    h2("Dose-Response Analysis"),

    p("Analyze tumor volume responses across dose levels using linear regression, ",
      "ANOVA, and non-linear dose-response models. ",
      "Upload and configure your data (including a Dose column) in the Data Upload tab first, ",
      "then run Tumor Growth Analysis before proceeding here."),

    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
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

      tags$div(class = "col-sm-9 main-panel-col",
        tags$button(
          class = "btn btn-default btn-xs sidebar-toggle-btn",
          type  = "button",
          title = "Toggle options panel",
          tags$i(class = "fa fa-chevron-left fa-fw"), " Collapse"
        ),
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
              strong("Per-dose summary statistics "), "of tumor volume at the selected time point, one row per dose level.",
              tags$ul(
                tags$li(strong("Dose"), " \u2014 dose level (concentration or quantity as provided in your data)."),
                tags$li(strong("n"), " \u2014 number of animals at this dose."),
                tags$li(strong("mean_volume"), " \u2014 arithmetic mean tumor volume across all animals at this dose."),
                tags$li(strong("sd_volume"), " \u2014 standard deviation of tumor volume; reflects within-dose variability."),
                tags$li(strong("sem_volume"), " \u2014 standard error of the mean (SD \u00f7 \u221an); indicates precision of the mean estimate."),
                tags$li(strong("ci95_lower / ci95_upper"), " \u2014 lower and upper bounds of the 95\u202f% confidence interval for the mean (using the t-distribution). Non-overlapping intervals between adjacent doses suggest a significant dose step.")
              )
            ),
            shinycssloaders::withSpinner(
              DTOutput("dr_summary_table")
            ),
            downloadButton("download_dr_table", "Download CSV")
          ),

          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  selectInput("dr_plot_selector", "Select Plot",
                    choices = c(
                      "Scatter Plot"       = "scatter",
                      "Boxplot by Dose"    = "boxplot",
                      "Mean by Dose (Bar)" = "barplot"
                    )
                  ),
                  hr(),
                  h6(icon("palette"), " Appearance"),
                  conditionalPanel("input.dr_plot_selector == 'scatter'",
                    colourpicker::colourInput("dr_scatter_pt_color",
                                             "Point colour",
                                             value = "#555555"),
                    colourpicker::colourInput("dr_scatter_line_color",
                                             "Regression line colour",
                                             value = "#2166AC")
                  ),
                  conditionalPanel("input.dr_plot_selector == 'boxplot'",
                    colourpicker::colourInput("dr_box_fill_color",
                                             "Box fill colour",
                                             value = "#4292C6")
                  ),
                  conditionalPanel("input.dr_plot_selector == 'barplot'",
                    colourpicker::colourInput("dr_bar_fill_color",
                                             "Bar fill colour",
                                             value = "#4292C6")
                  ),
                  hr(),
                  h6(style = "font-size: 1.15rem; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                     icon("sliders-h"), " Options"),
                  sliderInput("dr_opt_text_size", "Text size (pt)",
                              min = 10, max = 24, value = 16, step = 1),
                  hr(),
                  h6(style = "font-size: 1.15rem; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                     icon("heading"), " Title & Subtitle"),
                  textInput("dr_plot_title", "Title",
                            placeholder = "Leave blank for default"),
                  selectInput("dr_plot_title_align", "Title alignment",
                    choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                    selected = "center"),
                  textInput("dr_plot_subtitle", "Subtitle",
                            placeholder = "Optional subtitle"),
                  hr(),
                  h6(style = "font-size: 1.15rem; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                     icon("arrows-alt"), " Plot Size / Export"),
                  fluidRow(
                    column(6, numericInput("dr_plot_width",  "Width (px)",  800, min = 200, max = 5000, step = 50)),
                    column(6, numericInput("dr_plot_height", "Height (px)", 600, min = 200, max = 5000, step = 50))
                  ),
                  numericInput("dr_plot_dpi", "Resolution (DPI)", 100, min = 72, max = 600, step = 10),
                  tags$div(style = "margin-top: 10px;", uiOutput("dr_download_plot_btn"))
                )
              ),
              column(9,
                conditionalPanel("input.dr_plot_selector == 'scatter'",
                  h4("Scatter Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput("dr_plot_scatter", width = "auto", height = "auto")
                  )
                ),
                conditionalPanel("input.dr_plot_selector == 'boxplot'",
                  h4("Boxplot by Dose"),
                  shinycssloaders::withSpinner(
                    plotOutput("dr_plot_boxplot", width = "auto", height = "auto")
                  )
                ),
                conditionalPanel("input.dr_plot_selector == 'barplot'",
                  h4("Mean by Dose"),
                  shinycssloaders::withSpinner(
                    plotOutput("dr_plot_barplot", width = "auto", height = "auto")
                  )
                )
              )
            )
          ),

          tabPanel("Statistical Tests",
            tab_info(
              strong("Full statistical test results. "),
              "The text output includes linear regression (slope, intercept, p-value, R\u00b2), one-way ANOVA table, and polynomial trend p-values. Non-linear model AIC/BIC are shown when the drc package fit converges. ",
              "When ANOVA p\u202f<\u202f0.05 the Tukey HSD post-hoc table is shown below:",
              tags$ul(
                tags$li(strong("Comparison"), " \u2014 the two dose levels being compared (\u2018DoseA\u202f-\u202fDoseB\u2019)."),
                tags$li(strong("diff"), " \u2014 estimated mean difference in tumor volume (Dose A minus Dose B)."),
                tags$li(strong("lwr / upr"), " \u2014 lower and upper bounds of the Tukey-adjusted 95\u202f% confidence interval for the difference."),
                tags$li(strong("p adj"), " \u2014 Tukey family-wise error rate adjusted p-value. Values below 0.05 indicate a significant difference between that dose pair.")
              )
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

    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Running dose-response analysis...", value = 0.2)
    on.exit(progress$close(), add = TRUE)

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

    # Attach analysis_data and summary_stats for reactive plot rebuilding
    if (!isTRUE(result$error)) {
      tryCatch({
        # Re-prepare the analysis_data subset used by plots
        ad <- mouseExperiment::prepare_dose_data(
          data,
          dose_column      = dose_col,
          treatment_column = tx_col,
          volume_column    = vol_col,
          day_column       = time_col,
          id_column        = id_col,
          time_point       = time_point_val
        )
        ss <- mouseExperiment::generate_summary_statistics(ad,
          dose_column = dose_col, volume_column = vol_col)
        result$analysis_data <- ad
        result$summary_stats <- ss
      }, error = function(e) {\n        message(\"[DoseResponse] Data preparation failed: \", conditionMessage(e))\n        result$dr_prep_error <<- conditionMessage(e)\n      })
    }

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

  output$dr_download_plot_btn <- renderUI({
    btn_id <- switch(input$dr_plot_selector %||% "scatter",
      scatter = "download_dr_scatter",
      boxplot = "download_dr_boxplot",
      barplot = "download_dr_barplot",
      "download_dr_scatter"
    )
    downloadButton(btn_id, "Download Plot")
  })

  # ---- Plots (rebuilt reactively so colour options take effect immediately) ----
  output$dr_plot_scatter <- renderPlot({

    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(NULL)
    ad <- rv$dr_results$analysis_data
    req(ad)
    dose_col <- rv$dose_col %||% "Dose"
    vol_col  <- rv$vol_col  %||% "Volume"
    pt_col   <- input$dr_scatter_pt_color   %||% "#555555"
    ln_col   <- input$dr_scatter_line_color %||% "#2166AC"
    text_sz  <- input$dr_opt_text_size %||% 16
    p <- ggplot2::ggplot(ad, ggplot2::aes(x = .data[[dose_col]], y = .data[[vol_col]])) +
      ggplot2::geom_point(alpha = 0.7, color = pt_col) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = ln_col) +
      ggplot2::labs(title = "Linear Dose-Response Relationship",
                   x = "Dose", y = "Tumor Volume") +
      ggplot2::theme_classic()
    apply_plot_labels(p, input$dr_plot_title %||% "", input$dr_plot_subtitle %||% "",
                      input$dr_plot_title_align %||% "center", text_sz)
  }, width = function() input$dr_plot_width %||% 800,
     height = function() input$dr_plot_height %||% 600)

  output$dr_plot_boxplot <- renderPlot({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(NULL)
    ad <- rv$dr_results$analysis_data
    req(ad)
    dose_col    <- rv$dose_col %||% "Dose"
    vol_col     <- rv$vol_col  %||% "Volume"
    box_color   <- input$dr_box_fill_color %||% "#4292C6"
    text_sz     <- input$dr_opt_text_size %||% 16
    p <- ggplot2::ggplot(ad, ggplot2::aes(x = factor(.data[[dose_col]]), y = .data[[vol_col]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, fill = box_color, alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
      ggplot2::labs(title = "Tumor Volume by Dose Level",
                   x = "Dose", y = "Tumor Volume") +
      ggplot2::theme_classic()
    apply_plot_labels(p, input$dr_plot_title %||% "", input$dr_plot_subtitle %||% "",
                      input$dr_plot_title_align %||% "center", text_sz)
  }, width = function() input$dr_plot_width %||% 800,
     height = function() input$dr_plot_height %||% 600)

  output$dr_plot_barplot <- renderPlot({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(NULL)
    ss <- rv$dr_results$summary_stats
    req(ss)
    dose_col  <- rv$dose_col %||% "Dose"
    bar_color <- input$dr_bar_fill_color %||% "#4292C6"
    text_sz   <- input$dr_opt_text_size %||% 16
    p <- ggplot2::ggplot(ss, ggplot2::aes(x = .data[[dose_col]], y = .data[["mean_volume"]])) +
      ggplot2::geom_bar(stat = "identity", fill = bar_color, alpha = 0.8) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[["ci95_lower"]], ymax = .data[["ci95_upper"]]), width = 0.2) +
      ggplot2::labs(title = "Mean Tumor Volume by Dose Level (with 95% CI)",
                   x = "Dose", y = "Mean Tumor Volume") +
      ggplot2::theme_classic()
    apply_plot_labels(p, input$dr_plot_title %||% "", input$dr_plot_subtitle %||% "",
                      input$dr_plot_title_align %||% "center", text_sz)
  }, width = function() input$dr_plot_width %||% 800,
     height = function() input$dr_plot_height %||% 600)

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
      req(rv$dr_results)
      ad       <- rv$dr_results$analysis_data %||% rv$dr_results$plots$scatter$data
      dose_col <- rv$dose_col %||% "Dose"
      vol_col  <- rv$vol_col  %||% "Volume"
      pt_col   <- input$dr_scatter_pt_color   %||% "#555555"
      ln_col   <- input$dr_scatter_line_color %||% "#2166AC"
      p <- if (!is.null(ad)) {
        ggplot2::ggplot(ad, ggplot2::aes(x = .data[[dose_col]], y = .data[[vol_col]])) +
          ggplot2::geom_point(alpha = 0.7, color = pt_col) +
          ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = ln_col) +
          ggplot2::labs(title = "Linear Dose-Response Relationship",
                        x = "Dose", y = "Tumor Volume") +
          ggplot2::theme_classic()
      } else rv$dr_results$plots$scatter
      req(!is.null(p))
      png(file, width = input$dr_plot_width %||% 800,
                height = input$dr_plot_height %||% 600,
                units = "px", res = input$dr_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_dr_boxplot <- downloadHandler(
    filename = "dose_boxplot.png",
    content  = function(file) {
      req(rv$dr_results)
      ad       <- rv$dr_results$analysis_data
      dose_col <- rv$dose_col %||% "Dose"
      vol_col  <- rv$vol_col  %||% "Volume"
      bc       <- input$dr_box_fill_color %||% "#4292C6"
      p <- if (!is.null(ad)) {
        ggplot2::ggplot(ad, ggplot2::aes(x = factor(.data[[dose_col]]),
                                         y = .data[[vol_col]])) +
          ggplot2::geom_boxplot(outlier.shape = NA, fill = bc, alpha = 0.7) +
          ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
          ggplot2::labs(title = "Tumor Volume by Dose Level",
                        x = "Dose", y = "Tumor Volume") +
          ggplot2::theme_classic()
      } else rv$dr_results$plots$boxplot
      req(!is.null(p))
      png(file, width = input$dr_plot_width %||% 800,
                height = input$dr_plot_height %||% 600,
                units = "px", res = input$dr_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_dr_barplot <- downloadHandler(
    filename = "dose_barplot.png",
    content  = function(file) {
      req(rv$dr_results)
      ss       <- rv$dr_results$summary_stats
      dose_col <- rv$dose_col %||% "Dose"
      bc       <- input$dr_bar_fill_color %||% "#4292C6"
      p <- if (!is.null(ss)) {
        ggplot2::ggplot(ss, ggplot2::aes(x = .data[[dose_col]], y = .data[["mean_volume"]])) +
          ggplot2::geom_bar(stat = "identity", fill = bc, alpha = 0.8) +
          ggplot2::geom_errorbar(
            ggplot2::aes(ymin = .data[["ci95_lower"]], ymax = .data[["ci95_upper"]]), width = 0.2) +
          ggplot2::labs(title = "Mean Tumor Volume by Dose Level (with 95% CI)",
                        x = "Dose", y = "Mean Tumor Volume") +
          ggplot2::theme_classic()
      } else rv$dr_results$plots$barplot
      req(!is.null(p))
      png(file, width = input$dr_plot_width %||% 800,
                height = input$dr_plot_height %||% 600,
                units = "px", res = input$dr_plot_dpi %||% 100)
      print(p); dev.off()
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
