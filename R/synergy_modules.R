#' Drug Synergy Analysis UI
#' @keywords internal
synergyUI <- function() {
  fluidPage(
    h2("Drug Synergy Analysis"),
    
    p("Test whether a drug combination shows synergistic, additive, or antagonistic effects ",
      "using Bliss independence and Loewe additivity models. ",
      "Ensure tumor growth data is uploaded and configured in the Data Upload tab."),
    
    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Treatment Groups"),
          
          p(tags$small("Select which treatment group in your data corresponds to each role:")),
          
          uiOutput("syn_group_selectors_ui"),
          
          hr(),
          
          h4("Analysis Options"),
          
          uiOutput("syn_time_point_ui"),
          
          br(),
          
          actionButton("run_syn_analysis",
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),
          
          uiOutput("syn_status_msg"),

          plot_export_ui("syn")
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
          id = "syn_results_tabs",
          
          tabPanel("Summary",
            tab_info(strong("Synergy summary: "), "Overall assessment combining Bliss independence and Combination Index (CI) metrics. Reports whether the drug combination is synergistic, additive, or antagonistic at the evaluated timepoint, and lists the key quantitative scores."),
            shinycssloaders::withSpinner(
              verbatimTextOutput("syn_summary")
            )
          ),
          
          tabPanel("Synergy Metrics",
            tab_info(
              strong("Quantitative synergy scores. "),
              "When \\u2018Analyze over time\\u2019 is disabled, the table has two columns (",
              strong("Metric"), " / ", strong("Value"), ") summarising scores at the selected timepoint: ",
              tags$ul(
                tags$li(strong("Bliss Difference"), " \\u2014 observed combination effect minus the Bliss-independent expected effect (TGI\\u2009combo \\u2212 TGI\\u2009Bliss). Positive = synergy; negative = antagonism."),
                tags$li(strong("Loewe Difference"), " \\u2014 same concept using the Loewe additivity model."),
                tags$li(strong("Combination Index (CI)"), " \\u2014 CI\\u202f<\\u202f0.85 synergistic; 0.85\\u202f\\u2264\\u202fCI\\u202f\\u2264\\u202f1.15 additive; CI\\u202f>\\u202f1.15 antagonistic."),
                tags$li(strong("Interpretation"), " \\u2014 plain-language summary based on the CI threshold.")
              ),
              "When \\u2018Analyze over time\\u2019 is enabled, each row is a measured timepoint with columns:",
              tags$ul(
                tags$li(strong("Time_Point"), " \\u2014 study day."),
                tags$li(strong("TGI_Drug_A / TGI_Drug_B / TGI_Combo"), " \\u2014 % tumor growth inhibition for Drug A, Drug B, and the combination relative to control."),
                tags$li(strong("Bliss_Expected_TGI / Loewe_Expected_TGI"), " \\u2014 predicted TGI if the drugs were purely independent under each model."),
                tags$li(strong("Bliss_Difference / Loewe_Difference"), " \\u2014 TGI_Combo minus the expected TGI; positive values indicate synergy."),
                tags$li(strong("Combination_Index"), " \\u2014 CI at this timepoint (see thresholds above)."),
                tags$li(strong("P_Value_vs_Drug_A / P_Value_vs_Drug_B"), " \\u2014 two-sample t-test p-value comparing the combination group to each mono-therapy; values below 0.05 indicate a statistically significant difference."),
                tags$li(strong("Synergy_Assessment"), " \\u2014 plain-language classification at each timepoint.")
              )
            ),
            shinycssloaders::withSpinner(
              DTOutput("syn_metrics_table")
            ),
            downloadButton("download_syn_metrics", "Download Results")
          ),
          
          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  uiOutput("syn_plot_selector_ui"),
                  hr(),
                  uiOutput("syn_color_pickers_ui"),
                  hr(),
                  h6(icon("sliders-h"), " Options"),
                  sliderInput("syn_opt_text_size", "Text size (pt)",
                              min = 8, max = 20, value = 12, step = 1),
                  selectInput("syn_opt_legend_pos", "Legend position",
                    choices = c("Auto" = "auto", "None" = "none",
                                "Top" = "top", "Bottom" = "bottom",
                                "Left" = "left", "Right" = "right",
                                "In plot - top left"  = "inplot_tl",
                                "In plot - top right" = "inplot_tr"),
                    selected = "inplot_tr"
                  ),
                  conditionalPanel(
                    "input.syn_plot_selector == 'trend'",
                    checkboxInput("syn_opt_show_shading",
                                  "Show Synergy / Antagonism shading",
                                  value = TRUE)
                  )
                )
              ),
              column(9,
                conditionalPanel("input.syn_plot_selector == 'main'",
                  h4("Synergy Plot"),
                  tab_info(strong("Visual synergy comparison: "),
                    "bars showing observed combination effect alongside the expected effect ",
                    "if the two drugs acted independently (Bliss model). A combination bar ",
                    "taller than the expected bar indicates synergy at the evaluation timepoint."),
                  shinycssloaders::withSpinner(
                    plotOutput("syn_plot_main", height = "500px")
                  ),
                  downloadButton("download_syn_plot_main", "Download Plot")
                ),
                conditionalPanel("input.syn_plot_selector == 'trend'",
                  h4("Synergy Plot"),
                  tab_info(strong("Tumor Growth Inhibition over time: "),
                    "percent tumor growth inhibition for each treatment group at every measured timepoint. ",
                    "Synergy and antagonism annotations are shown where the combination significantly ",
                    "exceeds or falls short of the expected Bliss-independent effect."),
                  shinycssloaders::withSpinner(
                    plotOutput("syn_plot_trend", height = "500px")
                  ),
                  downloadButton("download_syn_plot_trend", "Download Plot")
                ),
                conditionalPanel("input.syn_plot_selector == 'ci'",
                  h4("Combination Index"),
                  tab_info(strong("Combination Index (CI) over time: "),
                    "CI values at each measured timepoint. The dashed reference line at ",
                    "CI\u202f=\u202f1 represents additivity. Points below 1 indicate synergism; ",
                    "points above 1 indicate antagonism."),
                  shinycssloaders::withSpinner(
                    plotOutput("syn_plot_ci", height = "500px")
                  ),
                  downloadButton("download_syn_plot_ci", "Download Plot")
                )
              )
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

  # Dynamic time-point selector: "All dates" (default) + individual Day values
  output$syn_time_point_ui <- renderUI({
    days <- NULL
    if (!is.null(rv$processed_data) && !is.null(rv$time_col) &&
        rv$time_col %in% colnames(rv$processed_data)) {
      days <- sort(unique(na.omit(rv$processed_data[[rv$time_col]])))
    } else if (!is.null(rv$uploaded_data)) {
      day_candidate <- if (!is.null(input$col_day) && nzchar(input$col_day)) input$col_day else "Day"
      if (day_candidate %in% colnames(rv$uploaded_data))
        days <- sort(unique(na.omit(suppressWarnings(
          as.numeric(rv$uploaded_data[[day_candidate]])))))
    }
    day_choices <- c("All dates" = "all")
    if (!is.null(days) && length(days) > 0)
      day_choices <- c(day_choices, setNames(as.character(days), paste("Day", days)))
    selectInput("syn_time_point",
                "Time Point for Analysis",
                choices  = day_choices,
                selected = "all")
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
    
    eval_tp <- if (!is.null(input$syn_time_point) && input$syn_time_point != "all") {
      as.numeric(input$syn_time_point)
    } else {
      NULL
    }

    rv$syn_over_time <- is.null(eval_tp)
    
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

    # Store the 4 group names so the Plots tab can build the color-named vector
    rv$syn_plot_inputs <- list(
      drug_a  = input$syn_drug_a,
      drug_b  = input$syn_drug_b,
      combo   = input$syn_combo,
      control = input$syn_control
    )
    
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
  
  # Colour pickers with real drug names once analysis has been run
  output$syn_color_pickers_ui <- renderUI({
    grps <- rv$syn_plot_inputs
    label_a     <- if (!is.null(grps)) grps$drug_a  else "Drug A"
    label_b     <- if (!is.null(grps)) grps$drug_b  else "Drug B"
    label_combo <- if (!is.null(grps)) grps$combo   else "Combination"
    label_ctrl  <- if (!is.null(grps)) paste0(grps$control, " / Bliss expected") else "Control / Bliss expected"
    tagList(
      h6(icon("palette"), " Appearance"),
      colourpicker::colourInput("syn_color_A",     label_a,     value = "#2166AC"),
      colourpicker::colourInput("syn_color_B",     label_b,     value = "#D6604D"),
      colourpicker::colourInput("syn_color_combo", label_combo, value = "#762A83"),
      colourpicker::colourInput("syn_color_ctrl",  label_ctrl,  value = "#999999")
    )
  })

  # Dynamic plot selector (Trend / CI only available for all-dates analysis)
  output$syn_plot_selector_ui <- renderUI({
    if (isTRUE(rv$syn_over_time)) {
      choices <- c("Synergy Plot" = "trend", "Combination Index" = "ci")
    } else {
      choices <- c("Synergy Plot" = "main")
    }
    selectInput("syn_plot_selector", "Select Plot", choices = choices)
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
        if (!is.null(res$peak_bliss_synergy)) {
          pb <- res$peak_bliss_synergy
          cat("Peak Bliss synergy at timepoint:", pb$Time_Point,
              " (Bliss Difference =", round(pb$Bliss_Difference, 1), "%)\n")
        }
        if (!is.null(res$peak_ci_synergy)) {
          pc <- res$peak_ci_synergy
          cat("Strongest CI synergy at timepoint:", pc$Time_Point,
              " (CI =", round(pc$Combination_Index, 2), ")\n")
        }
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
    datatable(df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE), rownames = FALSE)
  })
  
  # TGI Table removed (duplicate of Synergy Metrics)
  
  output$syn_plot_main <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)

    # Build colour vector matching the 4 synergy roles
    col_A     <- input$syn_color_A     %||% "#2166AC"
    col_B     <- input$syn_color_B     %||% "#D6604D"
    col_combo <- input$syn_color_combo %||% "#762A83"
    col_ctrl  <- input$syn_color_ctrl  %||% "#999999"

    if (isTRUE(rv$syn_over_time)) {
      # Over-time Bliss plot stored in rv uses "Drug A" / "Drug B" etc. as colour labels
      p <- rv$syn_results$plot_main
      req(p)
      print(p + ggplot2::scale_color_manual(
        values = c("Drug A" = col_A, "Drug B" = col_B,
                   "Combination" = col_combo, "Bliss Value" = col_ctrl)))
    } else {
      # Single timepoint - re-call plot_drug_synergy with custom_colors mapped to actual group names
      grps <- rv$syn_plot_inputs
      if (!is.null(grps)) {
        custom_colors <- c(grps$drug_a, grps$drug_b, grps$combo, grps$control)
        names(custom_colors) <- c(grps$drug_a, grps$drug_b, grps$combo, grps$control)
        custom_colors[grps$drug_a]  <- col_A
        custom_colors[grps$drug_b]  <- col_B
        custom_colors[grps$combo]   <- col_combo
        custom_colors[grps$control] <- col_ctrl
      } else {
        custom_colors <- NULL
      }
      tryCatch({
        p <- mouseExperiment::plot_drug_synergy(rv$syn_results,
                                                custom_colors = custom_colors)
        print(p)
      }, error = function(e) {
        if (!is.null(rv$syn_results$plot_main)) print(rv$syn_results$plot_main)
        else plot(1, 1, type = "n", main = "Synergy plot not available",
                  xlab = "", ylab = "")
      })
    }
  })

  output$syn_plot_trend <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || !isTRUE(rv$syn_over_time)) return(NULL)
    req(rv$syn_results$synergy_summary)

    syn_df      <- rv$syn_results$synergy_summary
    grps        <- rv$syn_plot_inputs %||%
                   list(drug_a = "Drug A", drug_b = "Drug B",
                        combo = "Combination", control = "Control")
    col_A       <- input$syn_color_A     %||% "#2166AC"
    col_B       <- input$syn_color_B     %||% "#D6604D"
    col_combo   <- input$syn_color_combo %||% "#762A83"
    col_ctrl    <- input$syn_color_ctrl  %||% "#999999"
    text_sz     <- input$syn_opt_text_size  %||% 12
    leg_pos     <- input$syn_opt_legend_pos %||% "inplot_tr"
    show_shade  <- isTRUE(input$syn_opt_show_shading)

    # Map internal color keys -> real drug name labels for the legend
    color_vals   <- c("Drug A"      = col_A,
                      "Drug B"      = col_B,
                      "Combination" = col_combo,
                      "Bliss Value" = col_ctrl)
    color_labels <- c("Drug A"      = grps$drug_a,
                      "Drug B"      = grps$drug_b,
                      "Combination" = grps$combo,
                      "Bliss Value" = paste0(grps$control, " (Bliss)"))

    p <- ggplot2::ggplot(syn_df, ggplot2::aes(x = Time_Point))

    # Synergy / antagonism shading between Combo TGI and Bliss expected TGI
    if (show_shade &&
        all(c("TGI_Combo", "Bliss_Expected_TGI") %in% colnames(syn_df))) {
      p <- p +
        # synergy region: combo > bliss (blue)
        ggplot2::geom_ribbon(
          ggplot2::aes(
            ymin = Bliss_Expected_TGI,
            ymax = ifelse(TGI_Combo >= Bliss_Expected_TGI,
                          TGI_Combo, Bliss_Expected_TGI)
          ),
          fill = "steelblue", alpha = 0.18
        ) +
        # antagonism region: combo < bliss (red)
        ggplot2::geom_ribbon(
          ggplot2::aes(
            ymin = ifelse(TGI_Combo < Bliss_Expected_TGI,
                          TGI_Combo, Bliss_Expected_TGI),
            ymax = Bliss_Expected_TGI
          ),
          fill = "tomato", alpha = 0.18
        )
    }

    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_A,         color = "Drug A")) +
      ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_B,         color = "Drug B")) +
      ggplot2::geom_line(ggplot2::aes(y = TGI_Combo,          color = "Combination")) +
      ggplot2::geom_line(ggplot2::aes(y = Bliss_Expected_TGI, color = "Bliss Value"),
                         alpha = 0.5, linewidth = 3) +
      ggplot2::scale_color_manual(values = color_vals, labels = color_labels) +
      ggplot2::labs(x = "Days", y = "Percent Tumor Growth Inhibition", color = NULL) +
      ggplot2::ylim(0, 100) +
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = text_sz))

    p <- apply_legend_pos(p, leg_pos)
    print(p)
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
  
  output$download_syn_plot_main <- downloadHandler(
    filename = "synergy_plot.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_main)) {
        png(file, width = input$syn_plot_width %||% 900,
                  height = input$syn_plot_height %||% 600,
                  units = "px", res = input$syn_plot_dpi %||% 100)
        print(rv$syn_results$plot_main); dev.off()
      }
    }
  )
  
  output$download_syn_plot_trend <- downloadHandler(
    filename = "synergy_trend.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_trend)) {
        png(file, width = input$syn_plot_width %||% 900,
                  height = input$syn_plot_height %||% 600,
                  units = "px", res = input$syn_plot_dpi %||% 100)
        print(rv$syn_results$plot_trend); dev.off()
      }
    }
  )
  
  output$download_syn_plot_ci <- downloadHandler(
    filename = "combination_index.png",
    content = function(file) {
      if (!is.null(rv$syn_results$plot_ci)) {
        png(file, width = input$syn_plot_width %||% 900,
                  height = input$syn_plot_height %||% 600,
                  units = "px", res = input$syn_plot_dpi %||% 100)
        print(rv$syn_results$plot_ci); dev.off()
      }
    }
  )
}
