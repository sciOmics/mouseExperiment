#' Power Analysis UI
#' @keywords internal
powerAnalysisUI <- function() {
  fluidPage(
    h2("Post-hoc Power Analysis"),
    
    p("Estimate statistical power and required sample sizes based on your uploaded tumor growth data. ",
      "Upload and configure data in the Data Upload tab first."),
    
    fluidRow(
      column(3,
        wellPanel(
          h4("Analysis Method"),
          
          selectInput("pw_method",
                     "Method",
                     choices = c(
                       "Parametric"  = "parametric",
                       "Simulation"  = "simulation",
                       "AUC-based"   = "auc"
                     ),
                     selected = "parametric"),
          
          conditionalPanel(
            condition = "input.pw_method == 'simulation'",
            numericInput("pw_n_simulations",
                        "Number of Simulations",
                        value = 1000, min = 100, max = 10000, step = 100)
          ),
          
          hr(),
          
          h4("Alpha Levels"),
          checkboxGroupInput("pw_alpha",
                            "Significance Levels",
                            choices  = c("0.05" = "0.05",
                                        "0.01" = "0.01",
                                        "0.10" = "0.10"),
                            selected = c("0.05", "0.01")),
          
          hr(),
          
          h4("Target Power"),
          checkboxGroupInput("pw_power",
                            "Power Targets",
                            choices  = c("0.80" = "0.80",
                                        "0.85" = "0.85",
                                        "0.90" = "0.90",
                                        "0.95" = "0.95"),
                            selected = c("0.80", "0.90")),

          hr(),

          h4("Custom Effect Sizes"),
          textInput("pw_custom_effect_sizes",
                   "Effect Sizes (comma-separated, optional)",
                   value = "",
                   placeholder = "e.g. 0.2, 0.5, 0.8, 1.0"),
          p(tags$small("Overrides the auto-detected effect sizes. Leave blank to use the values estimated from the data.")),

          br(),
          
          actionButton("run_pw_analysis",
                      "Run Power Analysis",
                      class = "btn-lg btn-primary btn-block"),
          
          uiOutput("pw_status_msg")
        )
      ),
      
      column(9,
        tabsetPanel(
          id = "pw_results_tabs",
          
          tabPanel("Summary",
            tab_info(strong("Power analysis summary: "), "overall observed power for the current sample size and key recommendations. Post-hoc power is calculated from the effect sizes observed in your data \u2014 interpret cautiously as the same data is used to both estimate effects and compute power."),
            shinycssloaders::withSpinner(
              verbatimTextOutput("pw_summary")
            )
          ),
          
          tabPanel("Effect Sizes",
            tab_info(
              strong("Cohen\u2019s\u202fd effect sizes "), "between each treatment group pair estimated from the observed data. Conventional benchmarks: ",
              strong("small\u202f\u2248\u202f0.2"), ", ",
              strong("medium\u202f\u2248\u202f0.5"), ", ",
              strong("large\u202f\u2248\u202f0.8"), ". Larger effect sizes require smaller samples to detect reliably.",
              tags$ul(
                tags$li(strong("Treatment"), " \u2014 the treatment group being compared."),
                tags$li(strong("Reference"), " \u2014 the group it is compared against."),
                tags$li(strong("Raw_Difference"), " \u2014 unadjusted mean difference in the outcome measure between Treatment and Reference."),
                tags$li(strong("Pooled_SD"), " \u2014 pooled within-group standard deviation used as the denominator for Cohen\u2019s d."),
                tags$li(strong("Standardized_Effect"), " \u2014 Cohen\u2019s d = Raw_Difference \u00f7 Pooled_SD. Positive values mean Treatment > Reference.")
              )
            ),
            shinycssloaders::withSpinner(
              DTOutput("pw_effect_sizes_table")
            ),
            downloadButton("download_pw_effect_sizes", "Download")
          ),
          
          tabPanel("Power Analysis",
            tab_info(
              strong("Achieved power (1\u202f\u2212\u202f\u03b2) "), "at the current observed N per group, for each combination of comparison, alpha level, and effect size scenario.",
              tags$ul(
                tags$li(strong("Treatment / Reference"), " \u2014 the pair of groups being evaluated."),
                tags$li(strong("Effect_Size"), " \u2014 Cohen\u2019s d used for this row (drawn from the observed effect sizes above or from a standard scenario grid)."),
                tags$li(strong("Alpha"), " \u2014 significance threshold (typically 0.05)."),
                tags$li(strong("Power"), " \u2014 probability of correctly detecting an effect of this size at this alpha with the current N. Values below 0.80 suggest the study may be underpowered.")
              ),
              em("Post-hoc power should be interpreted cautiously: the same data used to fit the model is used to estimate the effect size, which inflates apparent power.")
            ),
            shinycssloaders::withSpinner(
              DTOutput("pw_power_table")
            ),
            downloadButton("download_pw_power_table", "Download")
          ),
          
          tabPanel("Sample Size Recommendations",
            tab_info(
              strong("N per group needed "), "to achieve each target power level, for each combination of comparison, alpha, and effect size. Use these to plan future experiments.",
              tags$ul(
                tags$li(strong("Treatment / Reference"), " \u2014 the pair being evaluated."),
                tags$li(strong("Effect_Size"), " \u2014 Cohen\u2019s d for this scenario."),
                tags$li(strong("Alpha"), " \u2014 significance threshold."),
                tags$li(strong("Target_Power"), " \u2014 the desired power level (e.g. 0.80 or 0.90)."),
                tags$li(strong("Sample_Size"), " \u2014 N per group required to achieve the target power. Estimates for small effects or stringent targets may be impractically large.")
              )
            ),
            shinycssloaders::withSpinner(
              DTOutput("pw_sample_sizes_table")
            ),
            downloadButton("download_pw_sample_sizes", "Download")
          ),
          
          tabPanel("Power Curves",
            tab_info(strong("Power vs N per group: "), "shows how power increases with sample size for each detected effect size. Horizontal dashed lines mark conventional power thresholds (0.80, 0.90). Use these curves to identify the minimum N that achieves adequate power for your effect of interest."),
            shinycssloaders::withSpinner(
              plotOutput("pw_plot_curves", height = "500px")
            ),
            downloadButton("download_pw_plot", "Download Plot")
          )
        )
      )
    )
  )
}

#' Power Analysis Server Logic
#' @keywords internal
power_server <- function(input, output, session, rv) {
  
  rv$pw_results <- NULL
  
  observeEvent(input$run_pw_analysis, {
    req(rv$uploaded_data, rv$data_valid)
    
    if (length(input$pw_alpha) == 0) {
      output$pw_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), "Please select at least one alpha level.")
      })
      return()
    }
    
    if (length(input$pw_power) == 0) {
      output$pw_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), "Please select at least one power target.")
      })
      return()
    }
    
    output$pw_status_msg <- renderUI({
      div(class = "alert alert-info",
          icon("spinner", class = "fa-spin"), "Running power analysis...")
    })
    
    # Prefer preprocessed data (Volume + Day already resolved); fall back to raw
    if (is.null(rv$processed_data)) {
      output$pw_status_msg <- renderUI({
        div(class = "alert alert-warning",
            icon("exclamation-triangle"),
            "Please run the Tumor Growth Analysis first so that Volume and Day columns are resolved.")
      })
      return()
    }
    
    time_col <- rv$time_col
    vol_col  <- rv$vol_col
    tx_col   <- rv$tx_col
    id_col   <- rv$id_col
    
    alpha_vals  <- as.numeric(input$pw_alpha)
    power_vals  <- as.numeric(input$pw_power)
    n_sims      <- if (input$pw_method == "simulation") input$pw_n_simulations else 1000

    # Parse optional custom effect sizes
    custom_effect_sizes <- NULL
    if (!is.null(input$pw_custom_effect_sizes) && nzchar(trimws(input$pw_custom_effect_sizes))) {
      parsed <- suppressWarnings(as.numeric(strsplit(trimws(input$pw_custom_effect_sizes), "[,\\s]+")[[1]]))
      parsed <- parsed[!is.na(parsed) & parsed > 0]
      if (length(parsed) > 0) custom_effect_sizes <- parsed
    }

    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Running power analysis...", value = 0.2)
    on.exit(progress$close(), add = TRUE)

    result <- safe_analysis({
      mouseExperiment::post_power_analysis(
        data             = rv$processed_data,
        alpha            = alpha_vals,
        power            = power_vals,
        method           = input$pw_method,
        n_simulations    = n_sims,
        time_column      = time_col,
        volume_column    = vol_col,
        treatment_column = tx_col,
        id_column        = id_col,
        effect_sizes     = custom_effect_sizes
      )
    })

    # Inject run-time parameters into summary so renderText can display them
    if (!isTRUE(result$error)) {
      if (is.null(result$summary)) result$summary <- list()
      result$summary$method        <- input$pw_method
      result$summary$alpha_levels  <- alpha_vals
      result$summary$power_targets <- power_vals
    }

    rv$pw_results <- result
    
    if (isTRUE(result$error)) {
      output$pw_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"),
            "Error:", result$message)
      })
    } else {
      output$pw_status_msg <- renderUI({
        div(class = "alert alert-success",
            icon("check"), "Power analysis completed successfully!")
      })
    }
  })
  
  output$pw_summary <- renderText({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error)) return(paste("Error:", rv$pw_results$message))
    
    res <- rv$pw_results
    paste(capture.output({
      cat("=== POST-HOC POWER ANALYSIS SUMMARY ===\n\n")
      if (!is.null(res$summary)) {
        if (is.character(res$summary)) {
          cat(res$summary)
        } else {
          method_label <- switch(res$summary$method,
            parametric  = "Parametric (t-test based)",
            simulation  = "Simulation",
            auc         = "AUC-based",
            res$summary$method
          )
          cat("Method:        ", method_label, "\n")
          cat("Alpha levels:  ", paste(res$summary$alpha_levels,  collapse = ", "), "\n")
          cat("Power targets: ", paste(res$summary$power_targets, collapse = ", "), "\n")
          if (!is.null(res$summary$sample_sizes) && length(res$summary$sample_sizes) > 0) {
            cat("\nSample sizes per group:\n")
            for (nm in names(res$summary$sample_sizes)) {
              cat("  ", nm, ": ", res$summary$sample_sizes[[nm]], "\n", sep = "")
            }
          }
          cat("\n")
          if (!is.null(res$summary$overall_power)) {
            cat("Overall observed power:", round(res$summary$overall_power, 4), "\n\n")
          }
          if (!is.null(res$summary$recommendations)) {
            cat("Recommendations:\n")
            for (r in res$summary$recommendations) cat(" -", r, "\n")
          }
        }
      } else {
        cat("Power analysis completed. See individual tabs for detailed results.\n")
      }
    }), collapse = "\n")
  })
  
  output$pw_effect_sizes_table <- renderDT({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error) || is.null(rv$pw_results$effect_sizes)) {
      return(datatable(data.frame(Message = "Effect sizes not available")))
    }
    df <- as.data.frame(rv$pw_results$effect_sizes)
    df <- format_df_for_display(df)
    datatable(df, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
  })
  
  output$pw_power_table <- renderDT({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error) || is.null(rv$pw_results$post_power_analysis)) {
      return(datatable(data.frame(Message = "Power analysis table not available")))
    }
    df <- as.data.frame(rv$pw_results$post_power_analysis)
    df <- format_df_for_display(df)
    datatable(df, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })
  
  output$pw_sample_sizes_table <- renderDT({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error) || is.null(rv$pw_results$sample_size_estimates)) {
      return(datatable(data.frame(Message = "Sample size estimates not available")))
    }
    df <- as.data.frame(rv$pw_results$sample_size_estimates)
    df <- format_df_for_display(df)
    datatable(df, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })
  
  output$pw_plot_curves <- renderPlot({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error)) return(NULL)

    # Build power-vs-N curves from observed effect sizes
    es_df <- rv$pw_results$effect_sizes
    if (!is.null(es_df) && has_col(es_df, \"Standardized_Effect\")) {
      unique_effects <- sort(unique(round(
        abs(es_df$Standardized_Effect[!is.na(es_df$Standardized_Effect)]), 2
      )))
      unique_effects <- unique_effects[unique_effects > 0]
    } else {
      unique_effects <- c(0.2, 0.5, 0.8, 1.0, 1.5)
    }

    # Use the smallest user-selected alpha for the power curve visualization
    alpha_for_curve <- min(as.numeric(input$pw_alpha %||% 0.05))

    if (length(unique_effects) > 0) {
      n_seq <- seq(2, 60, by = 1)
      rows <- lapply(unique_effects, function(es) {
        data.frame(
          N            = n_seq,
          Power        = vapply(n_seq, function(n)
            tryCatch(
              stats::power.t.test(n = n, delta = es, sd = 1, sig.level = alpha_for_curve)$power,
              error = function(e) NA_real_
            ), numeric(1)),
          Effect_Size  = paste0("d\u202f=\u202f", round(es, 2)),
          stringsAsFactors = FALSE
        )
      })
      pw_df <- do.call(rbind, rows)

      print(
        ggplot2::ggplot(pw_df,
          ggplot2::aes(x = N, y = Power,
                       color = Effect_Size, group = Effect_Size)) +
          ggplot2::geom_line(linewidth = 0.9) +
          ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed",
                              colour = "grey40", linewidth = 0.6) +
          ggplot2::geom_hline(yintercept = 0.90, linetype = "dashed",
                              colour = "grey60", linewidth = 0.6) +
          ggplot2::annotate("text", x = 60, y = 0.81, label = "0.80",
                            hjust = 1, size = 3, colour = "grey40") +
          ggplot2::annotate("text", x = 60, y = 0.91, label = "0.90",
                            hjust = 1, size = 3, colour = "grey60") +
          ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
          ggplot2::labs(
            title  = "Power vs N per Group",
            subtitle = paste0("\u03b1 = ", alpha_for_curve,
                              ";  dashed lines: 0.80 and 0.90 power thresholds"),
            x      = "N per Group",
            y      = "Power",
            colour = "Effect Size\n(Cohen\u2019s d)"
          ) +
          ggplot2::theme_classic() +
          ggplot2::theme(legend.position = "right")
      )
    } else {
      plots <- rv$pw_results$plots
      if (!is.null(plots$power_curves)) {
        print(plots$power_curves)
      } else {
        plot(1, 1, type = "n", main = "Power curves not available", xlab = "", ylab = "")
      }
    }
  })
  
  # Download handlers
  output$download_pw_effect_sizes <- downloadHandler(
    filename = "power_effect_sizes.csv",
    content = function(file) {
      if (!is.null(rv$pw_results$effect_sizes))
        write.csv(as.data.frame(rv$pw_results$effect_sizes), file, row.names = FALSE)
    }
  )
  
  output$download_pw_power_table <- downloadHandler(
    filename = "power_analysis.csv",
    content = function(file) {
      if (!is.null(rv$pw_results$post_power_analysis))
        write.csv(as.data.frame(rv$pw_results$post_power_analysis), file, row.names = FALSE)
    }
  )
  
  output$download_pw_sample_sizes <- downloadHandler(
    filename = "sample_size_estimates.csv",
    content = function(file) {
      if (!is.null(rv$pw_results$sample_size_estimates))
        write.csv(as.data.frame(rv$pw_results$sample_size_estimates), file, row.names = FALSE)
    }
  )
  
  output$download_pw_plot <- downloadHandler(
    filename = "power_curves.png",
    content = function(file) {
      if (!is.null(rv$pw_results$plots$power_curves)) {
        png(file, width = 900, height = 600)
        print(rv$pw_results$plots$power_curves)
        dev.off()
      }
    }
  )
}
