#' Efficacy-Toxicity Analysis UI
#' @keywords internal
efficacyToxicityUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Efficacy-Toxicity Analysis"),

    p("Combine tumor growth inhibition with body weight toxicity data to assess the therapeutic window. ",
      "Requires both tumor volume and body weight columns."),

    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Analysis Parameters"),

          uiOutput(ns("et_data_status")),

          conditionalPanel(ns = ns,
            condition = "output.et_data_available == true",

            uiOutput(ns("et_reference_group_ui")),

            hr(),
            h5(strong("Efficacy Metric")),
            radioButtons(ns("et_efficacy_metric"), "Metric for scatter/bivariate plots",
              choices = c("Final TGI (%)" = "tgi",
                          "Tumor AUC" = "tumor_auc",
                          "Log Cell Kill" = "log_cell_kill"),
              selected = "tgi"
            ),

            hr(),
            h5(strong("Total Benefit")),
            sliderInput(ns("et_lambda"), "\u03bb (toxicity penalty weight)",
                        min = 0, max = 5, value = 1, step = 0.1),

            hr(),
            h5(strong("Weight-Corrected TGI")),
            sliderInput(ns("et_safety_threshold"), "Safety threshold (%)",
                        min = 5, max = 40, value = 20, step = 1),

            hr(),
            h5(strong("Tumor Weight Adjustment")),
            checkboxInput(ns("et_adjust_tumor"), "Subtract estimated tumor weight", value = TRUE),
            conditionalPanel(ns = ns, "input.et_adjust_tumor == true",
              numericInput(ns("et_tumor_density"), "Tumor density (g/cm\u00b3)",
                           value = 1.0, min = 0.1, max = 3.0, step = 0.1)
            ),

            br(),
            actionButton(ns("run_et_analysis"),
                        "Run Efficacy-Toxicity Analysis",
                        class = "btn-lg btn-primary btn-block")
          ),

          br(),
          uiOutput(ns("et_status_msg"))
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
          id = ns("et_results_tabs"),

          tabPanel("Therapeutic Window",
            tab_info(
              strong("Therapeutic Window Metric (TWM): "),
              "Ratio of efficacy (TGI%) to maximum body weight loss (%). ",
              "Higher TWM indicates a more favorable therapeutic index."
            ),
            fluidRow(
              column(3,
                wellPanel(
                  h6(icon("arrows-alt"), " Export"),
                  fluidRow(
                    column(6, numericInput(ns("et_twm_width"), "Width", 700, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("et_twm_height"), "Height", 500, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("et_twm_dpi"), "DPI", 100, min = 72, max = 600, step = 10),
                  downloadButton(ns("download_et_twm_plot"), "Download Plot")
                )
              ),
              column(9,
                shinycssloaders::withSpinner(
                  plotOutput(ns("et_twm_plot"), width = "auto", height = "auto")
                )
              )
            ),
            hr(),
            h4("TWM Table"),
            shinycssloaders::withSpinner(DTOutput(ns("et_twm_table"))),
            downloadButton(ns("download_et_twm_csv"), "Download CSV")
          ),

          tabPanel("Safety-Efficacy Plot",
            tab_info(
              strong("Per-mouse bivariate plot "),
              "of efficacy (x-axis) versus toxicity (y-axis). ",
              "Ideal treatments appear in the upper-right quadrant (high efficacy, low toxicity). ",
              "Group ellipses and centroids summarize each treatment arm."
            ),
            fluidRow(
              column(3,
                wellPanel(
                  checkboxInput(ns("et_show_ellipses"), "Show group ellipses", value = TRUE),
                  checkboxInput(ns("et_show_centroids"), "Show group centroids", value = TRUE),
                  sliderInput(ns("et_scatter_pt_size"), "Point size", 1, 6, 3, step = 0.5),
                  hr(),
                  h6(icon("arrows-alt"), " Export"),
                  fluidRow(
                    column(6, numericInput(ns("et_scatter_width"), "Width", 700, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("et_scatter_height"), "Height", 500, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("et_scatter_dpi"), "DPI", 100, min = 72, max = 600, step = 10),
                  downloadButton(ns("download_et_scatter_plot"), "Download Plot")
                )
              ),
              column(9,
                shinycssloaders::withSpinner(
                  plotOutput(ns("et_scatter_plot"), width = "auto", height = "auto")
                )
              )
            ),
            hr(),
            h4("Per-Mouse Data"),
            shinycssloaders::withSpinner(DTOutput(ns("et_per_mouse_table"))),
            h4("Group Summary"),
            shinycssloaders::withSpinner(DTOutput(ns("et_per_group_table"))),
            downloadButton(ns("download_et_bivariate_csv"), "Download Per-Mouse CSV")
          ),

          tabPanel("Total Benefit",
            tab_info(
              strong("Total Benefit Score: "),
              "B = AUC(efficacy) \u2212 \u03bb \u00d7 AUC(toxicity). ",
              "Positive values indicate net benefit. Adjust \u03bb in the sidebar to change the toxicity penalty."
            ),
            fluidRow(
              column(3,
                wellPanel(
                  h6(icon("arrows-alt"), " Export"),
                  fluidRow(
                    column(6, numericInput(ns("et_benefit_width"), "Width", 700, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("et_benefit_height"), "Height", 500, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("et_benefit_dpi"), "DPI", 100, min = 72, max = 600, step = 10),
                  downloadButton(ns("download_et_benefit_plot"), "Download Plot")
                )
              ),
              column(9,
                shinycssloaders::withSpinner(
                  plotOutput(ns("et_benefit_plot"), width = "auto", height = "auto")
                )
              )
            ),
            hr(),
            h4("Benefit Table"),
            shinycssloaders::withSpinner(DTOutput(ns("et_benefit_table"))),
            downloadButton(ns("download_et_benefit_csv"), "Download CSV")
          ),

          tabPanel("Weight-Corrected TGI",
            tab_info(
              strong("TGI recalculated "),
              "after excluding mice that exceeded the safety threshold for body weight loss. ",
              "Side-by-side comparison with standard TGI shows how toxicity exclusions affect efficacy estimates."
            ),
            fluidRow(
              column(3,
                wellPanel(
                  h6(icon("arrows-alt"), " Export"),
                  fluidRow(
                    column(6, numericInput(ns("et_wctgi_width"), "Width", 700, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("et_wctgi_height"), "Height", 500, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("et_wctgi_dpi"), "DPI", 100, min = 72, max = 600, step = 10),
                  downloadButton(ns("download_et_wctgi_plot"), "Download Plot")
                )
              ),
              column(9,
                shinycssloaders::withSpinner(
                  plotOutput(ns("et_wctgi_plot"), width = "auto", height = "auto")
                )
              )
            ),
            hr(),
            h4("Corrected vs Uncorrected TGI"),
            shinycssloaders::withSpinner(DTOutput(ns("et_wctgi_table"))),
            h4("Excluded Mice"),
            shinycssloaders::withSpinner(DTOutput(ns("et_excluded_table"))),
            downloadButton(ns("download_et_wctgi_csv"), "Download CSV")
          ),

          tabPanel("Correlation Heatmap",
            tab_info(
              strong("Correlation matrix "),
              "between efficacy and toxicity metrics across all mice. ",
              "Helps identify which toxicity measures most strongly track with efficacy."
            ),
            fluidRow(
              column(3,
                wellPanel(
                  h6(icon("arrows-alt"), " Export"),
                  fluidRow(
                    column(6, numericInput(ns("et_heatmap_width"), "Width", 600, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("et_heatmap_height"), "Height", 600, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("et_heatmap_dpi"), "DPI", 100, min = 72, max = 600, step = 10),
                  downloadButton(ns("download_et_heatmap"), "Download Heatmap")
                )
              ),
              column(9,
                shinycssloaders::withSpinner(
                  plotOutput(ns("et_heatmap_plot"), width = "auto", height = "auto")
                )
              )
            )
          ),

          tabPanel("Combined Report",
            tab_info(
              strong("Download a combined report "),
              "with all efficacy-toxicity metrics for every treatment group."
            ),
            shinycssloaders::withSpinner(DTOutput(ns("et_combined_table"))),
            downloadButton(ns("download_et_report"), "Download Full Report CSV")
          )
        )
      )
    )
  )
}

#' Efficacy-Toxicity Analysis Server
#' @keywords internal
efficacy_toxicity_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

  rv$et_results <- NULL
  current_et_twm_plot     <- reactiveVal(NULL)
  current_et_scatter_plot  <- reactiveVal(NULL)
  current_et_benefit_plot  <- reactiveVal(NULL)
  current_et_wctgi_plot    <- reactiveVal(NULL)
  current_et_heatmap_plot  <- reactiveVal(NULL)

  # Reactive: both weight and volume must be available

  et_data_available <- reactive({
    has_wt  <- !is.null(rv$weight_col) && rv$weight_col != ""
    has_vol <- !is.null(rv$vol_col) && rv$vol_col != ""
    has_wt && has_vol
  })

  output$et_data_available <- reactive({ et_data_available() })
  outputOptions(output, "et_data_available", suspendWhenHidden = FALSE)

  # Status banner
  output$et_data_status <- renderUI({
    if (!rv$data_valid) {
      return(status_alert("Please upload and configure data in the Data Upload tab first.", "warning"))
    }
    if (is.null(rv$processed_data)) {
      tryCatch(ensure_preprocessed(rv), error = function(e) NULL)
    }
    has_wt  <- !is.null(rv$weight_col) && rv$weight_col != ""
    has_vol <- !is.null(rv$vol_col) && rv$vol_col != ""
    if (!has_wt && !has_vol) {
      return(status_alert("Both Weight/Mass and Volume columns are required. Configure them in Data Upload.", "warning"))
    }
    if (!has_wt) {
      return(status_alert("No Weight/Mass column selected. Configure it in Data Upload.", "warning"))
    }
    if (!has_vol) {
      return(status_alert("No Volume column selected. Configure it in Data Upload.", "warning"))
    }
    status_alert("Ready for efficacy-toxicity analysis.", "success")
  })

  # Reference group selector
  output$et_reference_group_ui <- renderUI({
    req(rv$processed_data, rv$tx_col)
    groups <- sort(unique(as.character(na.omit(rv$processed_data[[rv$tx_col]]))))
    default_ctrl <- pick_reference_group(groups)
    selectInput(session$ns("et_reference_group"),
                "Control / Reference Group",
                choices = groups,
                selected = default_ctrl)
  })

  # ---- Run Analysis ----
  observeEvent(input$run_et_analysis, {
    req(rv$processed_data, et_data_available())

    output$et_status_msg <- renderUI({
      status_alert("Running efficacy-toxicity analysis...", "info")
    })

    data     <- rv$processed_data
    wt_col   <- rv$weight_col
    vol_col  <- rv$vol_col
    time_col <- rv$time_col
    tx_col   <- rv$tx_col
    id_col   <- rv$id_col
    ref_group <- input$et_reference_group %||% NULL

    adjust_tumor <- isTRUE(input$et_adjust_tumor)
    density_val  <- input$et_tumor_density %||% 1.0
    lambda       <- input$et_lambda %||% 1
    threshold    <- (input$et_safety_threshold %||% 20) / 100
    eff_metric   <- input$et_efficacy_metric %||% "tgi"

    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Running efficacy-toxicity analysis...", value = 0.1)
    on.exit(progress$close(), add = TRUE)

    results <- list()

    # 1. TWM
    progress$set(value = 0.2, message = "Computing Therapeutic Window Metric...")
    results$twm <- safe_analysis(module = "et_twm", {
      mouseExperiment::therapeutic_window_metric(
        df                  = data,
        volume_column       = vol_col,
        weight_column       = wt_col,
        time_column         = time_col,
        treatment_column    = tx_col,
        id_column           = id_col,
        adjust_tumor_weight = adjust_tumor,
        tumor_density       = density_val,
        reference_group     = ref_group
      )
    })

    # 2. Bivariate
    progress$set(value = 0.4, message = "Computing bivariate safety-efficacy data...")
    results$bivariate <- safe_analysis(module = "et_bivariate", {
      mouseExperiment::efficacy_toxicity_bivariate(
        df                  = data,
        volume_column       = vol_col,
        weight_column       = wt_col,
        time_column         = time_col,
        treatment_column    = tx_col,
        id_column           = id_col,
        adjust_tumor_weight = adjust_tumor,
        tumor_density       = density_val,
        efficacy_metric     = eff_metric,
        reference_group     = ref_group
      )
    })

    # 3. Total Benefit
    progress$set(value = 0.6, message = "Computing Total Benefit Score...")
    results$benefit <- safe_analysis(module = "et_benefit", {
      mouseExperiment::total_benefit_area(
        df                  = data,
        volume_column       = vol_col,
        weight_column       = wt_col,
        time_column         = time_col,
        treatment_column    = tx_col,
        id_column           = id_col,
        adjust_tumor_weight = adjust_tumor,
        tumor_density       = density_val,
        lambda              = lambda,
        reference_group     = ref_group
      )
    })

    # 4. Weight-Corrected TGI
    progress$set(value = 0.8, message = "Computing Weight-Corrected TGI...")
    results$wctgi <- safe_analysis(module = "et_wctgi", {
      mouseExperiment::weight_corrected_tgi(
        df                  = data,
        volume_column       = vol_col,
        weight_column       = wt_col,
        time_column         = time_col,
        treatment_column    = tx_col,
        id_column           = id_col,
        adjust_tumor_weight = adjust_tumor,
        tumor_density       = density_val,
        safety_threshold    = threshold,
        reference_group     = ref_group
      )
    })

    rv$et_results <- results

    progress$set(value = 1.0, message = "Done!")

    n_errors <- sum(vapply(results, function(x) isTRUE(x$error), logical(1)))
    if (n_errors == 4) {
      output$et_status_msg <- renderUI({
        status_alert("All analyses failed. Check your data and settings.", "danger")
      })
    } else if (n_errors > 0) {
      output$et_status_msg <- renderUI({
        status_alert(sprintf("Completed with %d warning(s). Some analyses may have failed.", n_errors), "warning")
      })
    } else {
      output$et_status_msg <- renderUI({
        status_alert("Efficacy-toxicity analysis completed successfully!", "success")
      })
    }
  })

  # ---- TWM Plot ----
  output$et_twm_plot <- renderPlot({
    req(rv$et_results, rv$et_results$twm)
    twm <- rv$et_results$twm
    if (isTRUE(twm$error)) return(NULL)

    twm_tbl <- twm$twm_table
    twm_tbl$Treatment <- factor(twm_tbl$Treatment, levels = twm_tbl$Treatment[order(twm_tbl$TWM)])

    p <- ggplot2::ggplot(twm_tbl, ggplot2::aes(x = .data[["Treatment"]], y = .data[["TWM"]],
                                                fill = .data[["Treatment"]])) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["TWM"]], 2)),
                         vjust = -0.5, size = 4) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = "Therapeutic Window Metric by Treatment",
                   subtitle = "TWM = TGI(%) / Max Weight Loss(%). Dashed line = TWM of 1.",
                   x = "Treatment", y = "TWM") +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

    current_et_twm_plot(p)
    p
  }, width = function() input$et_twm_width %||% 700,
     height = function() input$et_twm_height %||% 500)

  output$et_twm_table <- renderDT({
    req(rv$et_results, rv$et_results$twm)
    twm <- rv$et_results$twm
    if (isTRUE(twm$error)) return(datatable(data.frame(Message = twm$message)))
    tbl <- format_df_for_display(twm$twm_table)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ---- Safety-Efficacy Plot ----
  output$et_scatter_plot <- renderPlot({
    req(rv$et_results, rv$et_results$bivariate)
    biv <- rv$et_results$bivariate
    if (isTRUE(biv$error)) return(NULL)

    pm        <- biv$per_mouse
    pt_size   <- input$et_scatter_pt_size %||% 3
    show_ell  <- isTRUE(input$et_show_ellipses)
    show_cent <- isTRUE(input$et_show_centroids)

    eff_label <- switch(biv$efficacy_metric,
      tgi           = "TGI (%)",
      tumor_auc     = "Tumor AUC",
      log_cell_kill = "Log Cell Kill"
    )

    p <- ggplot2::ggplot(pm, ggplot2::aes(x = .data[["Efficacy"]],
                                           y = .data[["Max_Weight_Loss_Pct"]],
                                           color = .data[["Treatment"]])) +
      ggplot2::geom_point(size = pt_size, alpha = 0.7)

    if (show_ell && nrow(pm) > 3) {
      p <- p + ggplot2::stat_ellipse(level = 0.68, linewidth = 0.8,
                                     show.legend = FALSE)
    }

    if (show_cent) {
      pg <- biv$per_group
      if (!is.null(pg)) {
        p <- p + ggplot2::geom_point(
          data = pg,
          ggplot2::aes(x = .data[["Mean_Efficacy"]], y = .data[["Mean_Max_Weight_Loss_Pct"]]),
          shape = 18, size = 6
        )
      }
    }

    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Safety-Efficacy Scatter",
                   subtitle = paste0("Each point is one mouse. Efficacy: ", eff_label),
                   x = eff_label,
                   y = "Max Weight Loss (%)",
                   color = "Treatment") +
      ggplot2::theme_classic(base_size = 14)

    current_et_scatter_plot(p)
    p
  }, width = function() input$et_scatter_width %||% 700,
     height = function() input$et_scatter_height %||% 500)

  output$et_per_mouse_table <- renderDT({
    req(rv$et_results, rv$et_results$bivariate)
    biv <- rv$et_results$bivariate
    if (isTRUE(biv$error)) return(datatable(data.frame(Message = biv$message)))
    tbl <- format_df_for_display(biv$per_mouse)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  output$et_per_group_table <- renderDT({
    req(rv$et_results, rv$et_results$bivariate)
    biv <- rv$et_results$bivariate
    if (isTRUE(biv$error) || is.null(biv$per_group)) {
      return(datatable(data.frame(Message = "Not available")))
    }
    tbl <- format_df_for_display(biv$per_group)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ---- Total Benefit ----
  output$et_benefit_plot <- renderPlot({
    req(rv$et_results, rv$et_results$benefit)
    ben <- rv$et_results$benefit
    if (isTRUE(ben$error)) return(NULL)

    bt <- ben$benefit_table
    bt$Treatment <- factor(bt$Treatment, levels = bt$Treatment[order(bt$Benefit)])

    p <- ggplot2::ggplot(bt, ggplot2::aes(x = .data[["Treatment"]], y = .data[["Benefit"]],
                                           fill = .data[["Benefit"]] > 0)) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["Benefit"]], 1)),
                         vjust = ifelse(bt$Benefit >= 0, -0.5, 1.5), size = 4) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "#d7191c"),
                                 labels = c("TRUE" = "Net Benefit", "FALSE" = "Net Harm"),
                                 name = "") +
      ggplot2::labs(title = "Total Benefit Score by Treatment",
                   subtitle = sprintf("B = AUC(efficacy) \u2212 \u03bb \u00d7 AUC(toxicity), \u03bb = %.1f", ben$lambda),
                   x = "Treatment", y = "Benefit Score") +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

    current_et_benefit_plot(p)
    p
  }, width = function() input$et_benefit_width %||% 700,
     height = function() input$et_benefit_height %||% 500)

  output$et_benefit_table <- renderDT({
    req(rv$et_results, rv$et_results$benefit)
    ben <- rv$et_results$benefit
    if (isTRUE(ben$error)) return(datatable(data.frame(Message = ben$message)))
    tbl <- format_df_for_display(ben$benefit_table)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ---- Weight-Corrected TGI ----
  output$et_wctgi_plot <- renderPlot({
    req(rv$et_results, rv$et_results$wctgi)
    wc <- rv$et_results$wctgi
    if (isTRUE(wc$error)) return(NULL)

    comp <- wc$comparison
    # Pivot for side-by-side bar chart
    plot_df <- data.frame(
      Treatment = rep(comp$Treatment, 2),
      Type      = rep(c("Uncorrected TGI", "Weight-Corrected TGI"), each = nrow(comp)),
      TGI       = c(comp$Uncorrected_TGI, comp$Corrected_TGI),
      stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["Treatment"]],
                                                y = .data[["TGI"]],
                                                fill = .data[["Type"]])) +
      ggplot2::geom_col(position = ggplot2::position_dodge(0.7), width = 0.6, alpha = 0.8) +
      ggplot2::labs(title = "Weight-Corrected vs Standard TGI",
                   subtitle = sprintf("Safety threshold: %d%% weight loss",
                                      round(wc$safety_threshold * 100)),
                   x = "Treatment", y = "TGI (%)", fill = "") +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

    current_et_wctgi_plot(p)
    p
  }, width = function() input$et_wctgi_width %||% 700,
     height = function() input$et_wctgi_height %||% 500)

  output$et_wctgi_table <- renderDT({
    req(rv$et_results, rv$et_results$wctgi)
    wc <- rv$et_results$wctgi
    if (isTRUE(wc$error)) return(datatable(data.frame(Message = wc$message)))
    tbl <- format_df_for_display(wc$comparison)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$et_excluded_table <- renderDT({
    req(rv$et_results, rv$et_results$wctgi)
    wc <- rv$et_results$wctgi
    if (isTRUE(wc$error) || is.null(wc$excluded_summary)) {
      return(datatable(data.frame(Message = "Not available")))
    }
    tbl <- format_df_for_display(wc$excluded_summary)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ---- Correlation Heatmap ----
  output$et_heatmap_plot <- renderPlot({
    req(rv$et_results, rv$et_results$bivariate, rv$et_results$twm)
    biv <- rv$et_results$bivariate
    twm_r <- rv$et_results$twm
    if (isTRUE(biv$error)) return(NULL)

    pm <- biv$per_mouse
    # Build a combined metrics matrix
    metrics_df <- pm[, c("Efficacy", "Max_Weight_Loss_Pct"), drop = FALSE]
    colnames(metrics_df) <- c("Efficacy", "Max Weight Loss (%)")

    # Add AUC if available
    if (!isTRUE(rv$et_results$benefit$error) && !is.null(rv$et_results$benefit$efficacy_auc)) {
      eff_auc <- rv$et_results$benefit$efficacy_auc
      tox_auc <- rv$et_results$benefit$toxicity_auc
      if (nrow(eff_auc) == nrow(pm)) {
        metrics_df[["Efficacy AUC"]] <- eff_auc$AUC
      }
      if (nrow(tox_auc) == nrow(pm)) {
        metrics_df[["Toxicity AUC"]] <- tox_auc$AUC
      }
    }

    # Need at least 2 columns for correlation
    if (ncol(metrics_df) < 2) return(NULL)

    cor_mat <- stats::cor(metrics_df, use = "pairwise.complete.obs")

    # Simple heatmap with ggplot2
    cor_df <- expand.grid(
      Var1 = colnames(cor_mat),
      Var2 = colnames(cor_mat),
      stringsAsFactors = FALSE
    )
    cor_df$Correlation <- as.vector(cor_mat)

    p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = .data[["Var1"]], y = .data[["Var2"]],
                                               fill = .data[["Correlation"]])) +
      ggplot2::geom_tile(color = "white", linewidth = 1) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data[["Correlation"]])),
                         size = 5, color = "black") +
      ggplot2::scale_fill_gradient2(low = "#d7191c", mid = "white", high = "#2c7bb6",
                                    midpoint = 0, limits = c(-1, 1)) +
      ggplot2::labs(title = "Correlation Heatmap: Efficacy vs Toxicity Metrics",
                   x = "", y = "") +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    current_et_heatmap_plot(p)
    p
  }, width = function() input$et_heatmap_width %||% 600,
     height = function() input$et_heatmap_height %||% 600)

  # ---- Combined Report ----
  output$et_combined_table <- renderDT({
    req(rv$et_results)
    res <- rv$et_results

    # Start with TWM table
    report <- NULL
    if (!isTRUE(res$twm$error)) {
      twm_tbl <- res$twm$twm_table
      report <- twm_tbl
    }

    # Add benefit if available
    if (!isTRUE(res$benefit$error) && !is.null(report)) {
      bt <- res$benefit$benefit_table
      report <- merge(report, bt[, c("Treatment", "Benefit")], by = "Treatment", all.x = TRUE)
    }

    # Add corrected TGI if available
    if (!isTRUE(res$wctgi$error) && !is.null(report)) {
      comp <- res$wctgi$comparison
      report <- merge(report, comp[, c("Treatment", "Corrected_TGI", "N_Excluded")],
                      by = "Treatment", all.x = TRUE)
    }

    if (is.null(report)) {
      return(datatable(data.frame(Message = "No results available. Run the analysis first.")))
    }

    tbl <- format_df_for_display(report)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ---- Download Handlers ----
  output$download_et_twm_plot <- downloadHandler(
    filename = "therapeutic_window.png",
    content  = function(file) {
      p <- current_et_twm_plot(); req(!is.null(p))
      png(file, width = input$et_twm_width %||% 700,
          height = input$et_twm_height %||% 500,
          units = "px", res = input$et_twm_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_et_twm_csv <- downloadHandler(
    filename = "therapeutic_window.csv",
    content  = function(file) {
      req(rv$et_results, rv$et_results$twm)
      if (isTRUE(rv$et_results$twm$error)) return()
      utils::write.csv(rv$et_results$twm$twm_table, file, row.names = FALSE)
    }
  )

  output$download_et_scatter_plot <- downloadHandler(
    filename = "safety_efficacy_scatter.png",
    content  = function(file) {
      p <- current_et_scatter_plot(); req(!is.null(p))
      png(file, width = input$et_scatter_width %||% 700,
          height = input$et_scatter_height %||% 500,
          units = "px", res = input$et_scatter_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_et_bivariate_csv <- downloadHandler(
    filename = "safety_efficacy_per_mouse.csv",
    content  = function(file) {
      req(rv$et_results, rv$et_results$bivariate)
      if (isTRUE(rv$et_results$bivariate$error)) return()
      utils::write.csv(rv$et_results$bivariate$per_mouse, file, row.names = FALSE)
    }
  )

  output$download_et_benefit_plot <- downloadHandler(
    filename = "total_benefit.png",
    content  = function(file) {
      p <- current_et_benefit_plot(); req(!is.null(p))
      png(file, width = input$et_benefit_width %||% 700,
          height = input$et_benefit_height %||% 500,
          units = "px", res = input$et_benefit_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_et_benefit_csv <- downloadHandler(
    filename = "total_benefit.csv",
    content  = function(file) {
      req(rv$et_results, rv$et_results$benefit)
      if (isTRUE(rv$et_results$benefit$error)) return()
      utils::write.csv(rv$et_results$benefit$benefit_table, file, row.names = FALSE)
    }
  )

  output$download_et_wctgi_plot <- downloadHandler(
    filename = "weight_corrected_tgi.png",
    content  = function(file) {
      p <- current_et_wctgi_plot(); req(!is.null(p))
      png(file, width = input$et_wctgi_width %||% 700,
          height = input$et_wctgi_height %||% 500,
          units = "px", res = input$et_wctgi_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_et_wctgi_csv <- downloadHandler(
    filename = "weight_corrected_tgi.csv",
    content  = function(file) {
      req(rv$et_results, rv$et_results$wctgi)
      if (isTRUE(rv$et_results$wctgi$error)) return()
      utils::write.csv(rv$et_results$wctgi$comparison, file, row.names = FALSE)
    }
  )

  output$download_et_heatmap <- downloadHandler(
    filename = "correlation_heatmap.png",
    content  = function(file) {
      p <- current_et_heatmap_plot(); req(!is.null(p))
      png(file, width = input$et_heatmap_width %||% 600,
          height = input$et_heatmap_height %||% 600,
          units = "px", res = input$et_heatmap_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_et_report <- downloadHandler(
    filename = "efficacy_toxicity_report.csv",
    content  = function(file) {
      req(rv$et_results)
      res <- rv$et_results
      report <- NULL
      if (!isTRUE(res$twm$error)) report <- res$twm$twm_table
      if (!isTRUE(res$benefit$error) && !is.null(report)) {
        bt <- res$benefit$benefit_table
        report <- merge(report, bt[, c("Treatment", "Benefit")], by = "Treatment", all.x = TRUE)
      }
      if (!isTRUE(res$wctgi$error) && !is.null(report)) {
        comp <- res$wctgi$comparison
        report <- merge(report, comp[, c("Treatment", "Corrected_TGI", "N_Excluded")],
                        by = "Treatment", all.x = TRUE)
      }
      if (is.null(report)) return()
      utils::write.csv(report, file, row.names = FALSE)
    }
  )

  }) # end moduleServer
}
