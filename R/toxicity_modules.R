#' Toxicity Analysis UI
#' @keywords internal
toxicityUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Toxicity Analysis"),

    p("Analyze body weight toxicity, or combine with tumor efficacy to assess the therapeutic window. ",
      "Requires a Weight/Mass column configured in the Data Upload tab."),

    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Analysis Parameters"),

          uiOutput(ns("tox_weight_status")),

          conditionalPanel(ns = ns,
            condition = "output.tox_weight_available == true",

            selectInput(ns("tox_model_type"),
              "Statistical Model",
              choices = c(
                "Body Weight Analysis"       = "bw",
                "Efficacy-Toxicity Analysis" = "et"
              ),
              selected = "bw"
            ),
            helpText(
              "Body Weight Analysis models weight trajectories with mixed-effects models, ",
              "AUC, and time-to-threshold survival. ",
              "Efficacy-Toxicity Analysis also requires a tumor volume column and computes ",
              "the therapeutic window, safety-efficacy scatter, and total benefit score."
            ),

            uiOutput(ns("tox_reference_group_ui")),

            hr(),
            h5(strong("Tumor Weight Adjustment")),
            checkboxInput(ns("tox_adjust_tumor"), "Subtract estimated tumor weight", value = TRUE),
            conditionalPanel(ns = ns, "input.tox_adjust_tumor == true",
              numericInput(ns("tox_tumor_density"), "Tumor density (g/cm\u00b3)",
                           value = 1.0, min = 0.1, max = 3.0, step = 0.1),
              tags$p(tags$small("Tumor weight = Volume (mm\u00b3) \u00d7 density / 1000"))
            ),

            hr(),

            # Body Weight-specific options
            conditionalPanel(ns = ns, "input.tox_model_type == 'bw'",
              h5(strong("Mixed Model Options")),
              checkboxGroupInput(ns("tox_covariates"), "Covariates",
                choices = c("Tumor Volume" = "volume",
                            "Sex" = "sex",
                            "Initial Mass" = "initial_mass"),
                selected = "volume"
              ),
              radioButtons(ns("tox_estimation"), "Estimation method",
                choices = c("REML" = "REML", "ML" = "ML"),
                selected = "REML", inline = TRUE
              ),
              hr(),
              h5(strong("Weight Loss Threshold")),
              sliderInput(ns("tox_threshold"), "Threshold (%)",
                          min = 5, max = 40, value = 20, step = 1),
              uiOutput(ns("tox_baseline_day_ui"))
            ),

            # Efficacy-Toxicity-specific options
            conditionalPanel(ns = ns, "input.tox_model_type == 'et'",
              uiOutput(ns("tox_et_vol_status")),
              h5(strong("Efficacy Metric")),
              radioButtons(ns("tox_et_efficacy_metric"), NULL,
                choices = c("Final TGI (%)" = "tgi",
                            "Tumor AUC" = "tumor_auc",
                            "Log Cell Kill" = "log_cell_kill"),
                selected = "tgi"
              ),
              hr(),
              h5(strong("Total Benefit")),
              sliderInput(ns("tox_et_lambda"), "\u03bb (toxicity penalty)",
                          min = 0, max = 5, value = 1, step = 0.1),
              hr(),
              h5(strong("Safety Threshold")),
              sliderInput(ns("tox_et_safety_threshold"), "Safety threshold (%)",
                          min = 5, max = 40, value = 20, step = 1)
            ),

            br(),
            actionButton(ns("run_tox_analysis"),
                        "Run Analysis",
                        class = "btn-lg btn-primary btn-block")
          ),

          br(),
          uiOutput(ns("tox_status_msg"))
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
          id = ns("tox_results_tabs"),

          tabPanel("Summary",
            tab_info(
              strong("Analysis summary. "),
              "Body Weight Analysis: mixed-effects model output, fixed effects, and marginal means. ",
              "Efficacy-Toxicity Analysis: therapeutic window, benefit scores, and corrected TGI."
            ),
            shinycssloaders::withSpinner(verbatimTextOutput(ns("tox_summary")))
          ),

          tabPanel("Results",

            # BW results
            conditionalPanel(ns = ns, "output.tox_model_type_val == 'bw'",
              tab_info(
                strong("Body Weight Mixed Model. "),
                "Fixed effects test whether weight trajectories differ by treatment group over time. ",
                "Estimated marginal means give group-level estimates adjusted for covariates."
              ),
              h4("Fixed Effects"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_fixed_effects"))),
              downloadButton(ns("download_tox_model"), "Download Fixed Effects CSV"),
              hr(),
              h4("Estimated Marginal Means by Treatment"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_emmeans"))),
              hr(),
              h4("AUC Summary by Group"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_auc_summary"))),
              downloadButton(ns("download_tox_auc_csv"), "Download AUC CSV"),
              hr(),
              h4("Nadir Weight Analysis"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_nadir_table"))),
              hr(),
              h4("Pairwise AUC Comparisons"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_auc_comparisons"))),
              hr(),
              h4("Time-to-Threshold Event Data"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_event_table"))),
              downloadButton(ns("download_tox_events"), "Download Events CSV")
            ),

            # ET results
            conditionalPanel(ns = ns, "output.tox_model_type_val == 'et'",
              tab_info(
                strong("Efficacy-Toxicity Results. "),
                "TWM = TGI(%) / Max Body Weight Loss(%). ",
                "Total Benefit B = AUC(efficacy) \u2212 \u03bb \u00d7 AUC(toxicity). ",
                "Weight-corrected TGI excludes mice exceeding the safety threshold."
              ),
              h4("Therapeutic Window Metric"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_et_twm_table"))),
              downloadButton(ns("download_et_twm_csv"), "Download TWM CSV"),
              hr(),
              h4("Per-Mouse Safety-Efficacy Data"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_et_per_mouse_table"))),
              downloadButton(ns("download_et_bivariate_csv"), "Download Per-Mouse CSV"),
              h4("Group Summary"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_et_per_group_table"))),
              hr(),
              h4("Total Benefit Score"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_et_benefit_table"))),
              downloadButton(ns("download_et_benefit_csv"), "Download Benefit CSV"),
              hr(),
              h4("Weight-Corrected vs Standard TGI"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_et_wctgi_table"))),
              h4("Excluded Mice"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_et_excluded_table"))),
              downloadButton(ns("download_et_wctgi_csv"), "Download Corrected TGI CSV"),
              hr(),
              h4("Combined Report"),
              shinycssloaders::withSpinner(DTOutput(ns("tox_et_combined_table"))),
              downloadButton(ns("download_et_report"), "Download Combined Report CSV")
            )
          ),

          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  conditionalPanel(ns = ns, "output.tox_model_type_val == 'bw'",
                    selectInput(ns("tox_bw_plot_selector"), "Plot Type",
                      choices = c(
                        "Body Weight Trajectory" = "trajectory",
                        "Body Weight AUC"        = "auc",
                        "Weight Loss KM Curve"   = "km"
                      ),
                      selected = "trajectory"
                    )
                  ),
                  conditionalPanel(ns = ns, "output.tox_model_type_val == 'et'",
                    selectInput(ns("tox_et_plot_selector"), "Plot Type",
                      choices = c(
                        "Therapeutic Window"      = "twm",
                        "Safety-Efficacy Scatter" = "scatter",
                        "Total Benefit Score"     = "benefit",
                        "Weight-Corrected TGI"    = "wctgi",
                        "Correlation Heatmap"     = "heatmap"
                      ),
                      selected = "twm"
                    )
                  ),
                  conditionalPanel(ns = ns,
                    "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'trajectory'",
                    hr(),
                    radioButtons(ns("tox_traj_type"), "Display style",
                      choices = c("Mean \u00b1 SEM" = "mean_sem",
                                  "Mean \u00b1 SD"  = "mean_sd",
                                  "Individual mice" = "spaghetti"),
                      selected = "mean_sem"
                    )
                  ),
                  conditionalPanel(ns = ns,
                    "output.tox_model_type_val == 'et' && input.tox_et_plot_selector == 'scatter'",
                    hr(),
                    checkboxInput(ns("tox_et_show_ellipses"),  "Show group ellipses",  value = TRUE),
                    checkboxInput(ns("tox_et_show_centroids"), "Show group centroids", value = TRUE),
                    sliderInput(ns("tox_et_scatter_pt_size"), "Point size", 1, 6, 3, step = 0.5)
                  ),
                  ctrl_section("tox-scaling", "expand", "Scaling",
                    sliderInput(ns("tox_plot_text_size"), "Text size (pt)",
                                min = 10, max = 24, value = 14, step = 1),
                    conditionalPanel(ns = ns,
                      "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'trajectory'",
                      sliderInput(ns("tox_traj_line_width"), "Line width",
                                  min = 0.5, max = 3, value = 1, step = 0.25)
                    )
                  ),
                  ctrl_section("tox-export", "arrows-alt", "Plot Size / Export",
                    fluidRow(
                      column(6, numericInput(ns("tox_plot_width"),  "Width (px)",  800, min = 200, max = 5000, step = 50)),
                      column(6, numericInput(ns("tox_plot_height"), "Height (px)", 500, min = 200, max = 5000, step = 50))
                    ),
                    fluidRow(
                      column(6, numericInput(ns("tox_plot_dpi"), "DPI", 100, min = 72, max = 600, step = 10)),
                      column(6, tags$div(style = "margin-top: 25px;",
                               downloadButton(ns("download_tox_active_plot"), "Download")))
                    )
                  )
                )
              ),
              column(9,
                conditionalPanel(ns = ns,
                  "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'trajectory'",
                  h4("Body Weight Trajectory"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_trajectory_plot"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,
                  "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'auc'",
                  h4("Body Weight AUC (% Change from Baseline)"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_auc_plot"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,
                  "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'km'",
                  h4("Weight Loss Kaplan-Meier Curve"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_km_plot"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,
                  "output.tox_model_type_val == 'et' && input.tox_et_plot_selector == 'twm'",
                  h4("Therapeutic Window Metric"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_et_twm_plot"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,
                  "output.tox_model_type_val == 'et' && input.tox_et_plot_selector == 'scatter'",
                  h4("Safety-Efficacy Scatter"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_et_scatter_plot"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,
                  "output.tox_model_type_val == 'et' && input.tox_et_plot_selector == 'benefit'",
                  h4("Total Benefit Score"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_et_benefit_plot"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,
                  "output.tox_model_type_val == 'et' && input.tox_et_plot_selector == 'wctgi'",
                  h4("Weight-Corrected vs Standard TGI"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_et_wctgi_plot"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,
                  "output.tox_model_type_val == 'et' && input.tox_et_plot_selector == 'heatmap'",
                  h4("Correlation Heatmap"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_et_heatmap_plot"), width = "auto", height = "auto")
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

#' Toxicity Analysis Server
#' @keywords internal
toxicity_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

  rv$tox_results <- NULL

  tox_weight_available <- reactive({
    !is.null(rv$weight_col) && rv$weight_col != ""
  })
  output$tox_weight_available <- reactive({ tox_weight_available() })
  outputOptions(output, "tox_weight_available", suspendWhenHidden = FALSE)

  output$tox_model_type_val <- reactive({ input$tox_model_type %||% "bw" })
  outputOptions(output, "tox_model_type_val", suspendWhenHidden = FALSE)

  output$tox_weight_status <- renderUI({
    if (!rv$data_valid) {
      return(status_alert("Please upload and configure data in the Data Upload tab first.", "warning"))
    }
    if (is.null(rv$processed_data)) {
      tryCatch(ensure_preprocessed(rv), error = function(e) NULL)
    }
    if (!tox_weight_available()) {
      return(status_alert(
        "No Weight/Mass column selected. Configure it in Data Upload \u2192 Optional Columns.",
        "warning"))
    }
    status_alert(paste("Weight column:", rv$weight_col), "success")
  })

  output$tox_et_vol_status <- renderUI({
    has_vol <- !is.null(rv$vol_col) && rv$vol_col != ""
    if (!has_vol)
      return(status_alert("Efficacy-Toxicity Analysis also requires a tumor volume column.", "warning"))
    invisible(NULL)
  })

  output$tox_reference_group_ui <- renderUI({
    req(rv$processed_data, rv$tx_col)
    groups <- sort(unique(as.character(na.omit(rv$processed_data[[rv$tx_col]]))))
    selectInput(session$ns("tox_reference_group"), "Control / Reference Group",
                choices = groups, selected = pick_reference_group(groups))
  })

  output$tox_baseline_day_ui <- renderUI({
    req(rv$processed_data, rv$time_col)
    days <- sort(unique(na.omit(rv$processed_data[[rv$time_col]])))
    day_choices <- c("First measurement" = "first",
                     setNames(as.character(days), paste("Day", days)))
    selectInput(session$ns("tox_baseline_day"), "Baseline day",
                choices = day_choices, selected = "first")
  })

  observe({
    req(rv$processed_data)
    has_vol <- !is.null(rv$vol_col)  && rv$vol_col  != ""
    has_sex <- !is.null(rv$sex_col) && rv$sex_col != ""
    available <- "initial_mass"
    if (has_vol) available <- c("volume", available)
    if (has_sex) available <- c(available, "sex")
    choices  <- c("Tumor Volume" = "volume", "Sex" = "sex", "Initial Mass" = "initial_mass")
    choices  <- choices[choices %in% available]
    selected <- intersect(input$tox_covariates %||% "volume", available)
    updateCheckboxGroupInput(session, "tox_covariates", choices = choices, selected = selected)
  })

  # ---- Run Analysis ----
  observeEvent(input$run_tox_analysis, {
    req(rv$processed_data, tox_weight_available())

    model_type <- input$tox_model_type %||% "bw"

    if (model_type == "et" && (is.null(rv$vol_col) || rv$vol_col == "")) {
      output$tox_status_msg <- renderUI({
        status_alert("Efficacy-Toxicity Analysis requires a tumor volume column.", "danger")
      })
      return()
    }

    output$tox_status_msg <- renderUI({ status_alert("Running analysis...", "info") })

    data      <- rv$processed_data
    wt_col    <- rv$weight_col
    vol_col   <- rv$vol_col
    time_col  <- rv$time_col
    tx_col    <- rv$tx_col
    id_col    <- rv$id_col
    sex_col   <- rv$sex_col
    cage_col  <- rv$cage_col
    ref_group <- input$tox_reference_group %||% NULL
    has_vol   <- !is.null(vol_col) && vol_col %in% names(data)
    adjust    <- isTRUE(input$tox_adjust_tumor) && has_vol
    density   <- input$tox_tumor_density %||% 1.0

    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Running analysis...", value = 0.05)
    on.exit(progress$close(), add = TRUE)

    results            <- list()
    results$model_type <- model_type

    if (model_type == "bw") {
      covs         <- input$tox_covariates   %||% character(0)
      est_method   <- input$tox_estimation   %||% "REML"
      threshold    <- (input$tox_threshold   %||% 20) / 100
      baseline_raw <- input$tox_baseline_day %||% "first"
      baseline_day <- if (baseline_raw == "first") NULL else as.numeric(baseline_raw)

      progress$set(value = 0.2, message = "Fitting mixed model...")
      results$mixed_model <- safe_analysis(module = "tox_mixed_model", {
        mouseExperiment::analyze_body_weight(
          df = data, weight_column = wt_col, time_column = time_col,
          treatment_column = tx_col, id_column = id_col,
          volume_column = if (has_vol) vol_col else NULL,
          sex_column = sex_col, cage_column = cage_col,
          adjust_tumor_weight = adjust, tumor_density = density,
          covariates = covs, estimation = est_method, reference_group = ref_group
        )
      })

      progress$set(value = 0.5, message = "Computing AUC...")
      results$bw_auc <- safe_analysis(module = "tox_auc", {
        mouseExperiment::body_weight_auc(
          df = data, weight_column = wt_col, time_column = time_col,
          treatment_column = tx_col, id_column = id_col,
          volume_column = if (has_vol) vol_col else NULL,
          adjust_tumor_weight = adjust, tumor_density = density, reference_group = ref_group
        )
      })

      progress$set(value = 0.75, message = "Time-to-threshold analysis...")
      results$weight_loss <- safe_analysis(module = "tox_threshold", {
        mouseExperiment::weight_loss_threshold(
          df = data, weight_column = wt_col, time_column = time_col,
          treatment_column = tx_col, id_column = id_col,
          volume_column = if (has_vol) vol_col else NULL,
          adjust_tumor_weight = adjust, tumor_density = density,
          threshold = threshold, baseline_day = baseline_day, reference_group = ref_group
        )
      })

      results$trajectory_data <- tryCatch({
        wd <- data.frame(
          ID        = as.character(data[[id_col]]),
          Treatment = as.character(data[[tx_col]]),
          Day       = as.numeric(data[[time_col]]),
          Weight    = as.numeric(data[[wt_col]]),
          stringsAsFactors = FALSE
        )
        if (has_vol && adjust)
          wd$Weight <- wd$Weight - (as.numeric(data[[vol_col]]) / 1000 * density)
        wd[!is.na(wd$Weight) & !is.na(wd$Day), ]
      }, error = function(e) NULL)

      n_err   <- sum(vapply(results[c("mixed_model", "bw_auc", "weight_loss")],
                             function(x) isTRUE(x$error), logical(1)))
      n_total <- 3L

    } else {
      eff_metric <- input$tox_et_efficacy_metric  %||% "tgi"
      lambda     <- input$tox_et_lambda            %||% 1
      threshold  <- (input$tox_et_safety_threshold %||% 20) / 100

      progress$set(value = 0.2, message = "Computing Therapeutic Window Metric...")
      results$twm <- safe_analysis(module = "et_twm", {
        mouseExperiment::therapeutic_window_metric(
          df = data, volume_column = vol_col, weight_column = wt_col,
          time_column = time_col, treatment_column = tx_col, id_column = id_col,
          adjust_tumor_weight = adjust, tumor_density = density, reference_group = ref_group
        )
      })

      progress$set(value = 0.4, message = "Computing bivariate data...")
      results$bivariate <- safe_analysis(module = "et_bivariate", {
        mouseExperiment::efficacy_toxicity_bivariate(
          df = data, volume_column = vol_col, weight_column = wt_col,
          time_column = time_col, treatment_column = tx_col, id_column = id_col,
          adjust_tumor_weight = adjust, tumor_density = density,
          efficacy_metric = eff_metric, reference_group = ref_group
        )
      })

      progress$set(value = 0.65, message = "Computing Total Benefit Score...")
      results$benefit <- safe_analysis(module = "et_benefit", {
        mouseExperiment::total_benefit_area(
          df = data, volume_column = vol_col, weight_column = wt_col,
          time_column = time_col, treatment_column = tx_col, id_column = id_col,
          adjust_tumor_weight = adjust, tumor_density = density,
          lambda = lambda, reference_group = ref_group
        )
      })

      progress$set(value = 0.85, message = "Computing Weight-Corrected TGI...")
      results$wctgi <- safe_analysis(module = "et_wctgi", {
        mouseExperiment::weight_corrected_tgi(
          df = data, volume_column = vol_col, weight_column = wt_col,
          time_column = time_col, treatment_column = tx_col, id_column = id_col,
          adjust_tumor_weight = adjust, tumor_density = density,
          safety_threshold = threshold, reference_group = ref_group
        )
      })

      n_err   <- sum(vapply(results[c("twm", "bivariate", "benefit", "wctgi")],
                             function(x) isTRUE(x$error), logical(1)))
      n_total <- 4L
    }

    rv$tox_results <- results
    progress$set(value = 1)

    if (n_err == n_total) {
      output$tox_status_msg <- renderUI({
        status_alert("All analyses failed. Check your data and settings.", "danger")
      })
    } else if (n_err > 0) {
      output$tox_status_msg <- renderUI({
        status_alert(sprintf("Completed with %d warning(s). Some results may be unavailable.", n_err), "warning")
      })
    } else {
      output$tox_status_msg <- renderUI({
        status_alert("Analysis completed successfully!", "success")
      })
    }
  })

  # ---- Summary ----
  output$tox_summary <- renderText({
    req(rv$tox_results)
    res <- rv$tox_results
    if (res$model_type == "bw") {
      mm <- res$mixed_model
      if (!is.null(mm) && !isTRUE(mm$error)) return(mm$summary_text)
      if (!is.null(mm) && isTRUE(mm$error))  return(paste("Mixed model error:", mm$message))
      return("No results. Run the analysis first.")
    }
    lines <- c("=== EFFICACY-TOXICITY ANALYSIS SUMMARY ===", "")
    if (!is.null(res$twm) && !isTRUE(res$twm$error)) {
      tbl   <- res$twm$twm_table
      lines <- c(lines, "Therapeutic Window Metric (TWM = TGI / Max Weight Loss%):", "")
      for (i in seq_len(nrow(tbl)))
        lines <- c(lines, sprintf("  %s:  TWM = %.2f  |  TGI = %.1f%%  |  Max Weight Loss = %.1f%%",
                                  tbl$Treatment[i], tbl$TWM[i],
                                  tbl$TGI_Pct[i], tbl$Max_Weight_Loss_Pct[i]))
      lines <- c(lines, "")
    }
    if (!is.null(res$benefit) && !isTRUE(res$benefit$error)) {
      bt    <- res$benefit$benefit_table
      lines <- c(lines, sprintf("Total Benefit (lambda = %.1f):", res$benefit$lambda), "")
      for (i in seq_len(nrow(bt)))
        lines <- c(lines, sprintf("  %s:  Benefit = %.2f", bt$Treatment[i], bt$Benefit[i]))
    }
    paste(lines, collapse = "\n")
  })

  # ---- BW: Result tables ----
  output$tox_fixed_effects <- renderDT({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error)) return(datatable(data.frame(Message = mm$message)))
    datatable(format_df_for_display(mm$fixed_effects),
              options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  output$tox_emmeans <- renderDT({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$emmeans_table))
      return(datatable(data.frame(Message = "Not available")))
    datatable(format_df_for_display(mm$emmeans_table),
              options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_auc_summary <- renderDT({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error)) return(datatable(data.frame(Message = auc$message)))
    datatable(format_df_for_display(auc$auc_summary),
              options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_nadir_table <- renderDT({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error)) return(datatable(data.frame(Message = "Not available")))
    datatable(format_df_for_display(auc$nadir_data),
              options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  output$tox_auc_comparisons <- renderDT({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error) || is.null(auc$comparisons))
      return(datatable(data.frame(Message = "Not available")))
    datatable(format_df_for_display(auc$comparisons),
              options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_event_table <- renderDT({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(datatable(data.frame(Message = wl$message)))
    datatable(format_df_for_display(wl$event_data),
              options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  # ---- ET: Result tables ----
  output$tox_et_twm_table <- renderDT({
    req(rv$tox_results, rv$tox_results$twm)
    twm <- rv$tox_results$twm
    if (isTRUE(twm$error)) return(datatable(data.frame(Message = twm$message)))
    datatable(format_df_for_display(twm$twm_table),
              options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_et_per_mouse_table <- renderDT({
    req(rv$tox_results, rv$tox_results$bivariate)
    biv <- rv$tox_results$bivariate
    if (isTRUE(biv$error)) return(datatable(data.frame(Message = biv$message)))
    datatable(format_df_for_display(biv$per_mouse),
              options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  output$tox_et_per_group_table <- renderDT({
    req(rv$tox_results, rv$tox_results$bivariate)
    biv <- rv$tox_results$bivariate
    if (isTRUE(biv$error) || is.null(biv$per_group))
      return(datatable(data.frame(Message = "Not available")))
    datatable(format_df_for_display(biv$per_group),
              options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_et_benefit_table <- renderDT({
    req(rv$tox_results, rv$tox_results$benefit)
    ben <- rv$tox_results$benefit
    if (isTRUE(ben$error)) return(datatable(data.frame(Message = ben$message)))
    datatable(format_df_for_display(ben$benefit_table),
              options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_et_wctgi_table <- renderDT({
    req(rv$tox_results, rv$tox_results$wctgi)
    wc <- rv$tox_results$wctgi
    if (isTRUE(wc$error)) return(datatable(data.frame(Message = wc$message)))
    datatable(format_df_for_display(wc$comparison),
              options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_et_excluded_table <- renderDT({
    req(rv$tox_results, rv$tox_results$wctgi)
    wc <- rv$tox_results$wctgi
    if (isTRUE(wc$error) || is.null(wc$excluded_summary))
      return(datatable(data.frame(Message = "Not available")))
    datatable(format_df_for_display(wc$excluded_summary),
              options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_et_combined_table <- renderDT({
    req(rv$tox_results)
    res    <- rv$tox_results
    report <- NULL
    if (!isTRUE(res$twm$error))    report <- res$twm$twm_table
    if (!isTRUE(res$benefit$error) && !is.null(report))
      report <- merge(report, res$benefit$benefit_table[, c("Treatment","Benefit")],
                      by = "Treatment", all.x = TRUE)
    if (!isTRUE(res$wctgi$error)   && !is.null(report))
      report <- merge(report,
                      res$wctgi$comparison[, c("Treatment","Corrected_TGI","N_Excluded")],
                      by = "Treatment", all.x = TRUE)
    if (is.null(report))
      return(datatable(data.frame(Message = "No results. Run the analysis first.")))
    datatable(format_df_for_display(report), options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ---- BW Plots ----
  output$tox_trajectory_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$trajectory_data)
    td <- rv$tox_results$trajectory_data
    req(nrow(td) > 0)
    traj_type <- input$tox_traj_type %||% "mean_sem"
    text_sz   <- input$tox_plot_text_size %||% 14
    line_w    <- input$tox_traj_line_width %||% 1

    if (traj_type == "spaghetti") {
      td$.__LineGroup <- paste(td[["Treatment"]], td[["ID"]], sep = "\u00b7")
      ggplot2::ggplot(td, ggplot2::aes(.data[["Day"]], .data[["Weight"]],
                                        group = .data[[".__LineGroup"]], color = .data[["Treatment"]])) +
        ggplot2::geom_line(alpha = 0.6, linewidth = line_w) +
        ggplot2::geom_point(size = 1.5, alpha = 0.6) +
        ggplot2::labs(title = "Body Weight Trajectories",
                     x = "Day", y = "Weight (g)", color = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
    } else {
      agg <- stats::aggregate(Weight ~ Treatment + Day, data = td, FUN = function(x)
        c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = length(x)))
      af  <- do.call(data.frame, agg)
      names(af) <- c("Treatment", "Day", "Mean", "SD", "N")
      af$SEM   <- af$SD / sqrt(af$N)
      af$Error <- if (traj_type == "mean_sd") af$SD else af$SEM
      err_lbl  <- if (traj_type == "mean_sd") "SD" else "SEM"
      ggplot2::ggplot(af, ggplot2::aes(.data[["Day"]], .data[["Mean"]],
                                        color = .data[["Treatment"]], fill = .data[["Treatment"]])) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["Mean"]] - .data[["Error"]],
                                          ymax = .data[["Mean"]] + .data[["Error"]]),
                             alpha = 0.2, color = NA) +
        ggplot2::geom_line(linewidth = line_w) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(title = paste0("Body Weight Over Time (Mean \u00b1 ", err_lbl, ")"),
                     x = "Day", y = "Weight (g)", color = "Treatment", fill = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
    }
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_auc_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error)) return(NULL)
    text_sz <- input$tox_plot_text_size %||% 14
    ggplot2::ggplot(auc$auc_per_mouse,
      ggplot2::aes(.data[["Treatment"]], .data[["AUC_Pct_Change"]], fill = .data[["Treatment"]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Body Weight AUC (% Change from Baseline)",
                   x = "Treatment", y = "AUC of % Weight Change") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "none")
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_km_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(NULL)
    text_sz <- input$tox_plot_text_size %||% 14
    km_plot <- tryCatch(
      survminer::ggsurvplot(wl$km_fit, data = wl$event_data,
        risk.table = TRUE, pval = TRUE, conf.int = TRUE,
        title = paste0("Time to ", round(wl$threshold * 100), "% Body Weight Loss"),
        xlab = "Day", ylab = "Proportion Without Threshold Weight Loss",
        ggtheme = ggplot2::theme_classic(base_size = text_sz)),
      error = function(e) NULL)
    req(!is.null(km_plot))
    print(km_plot)
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  # ---- ET Plots ----
  output$tox_et_twm_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$twm)
    twm <- rv$tox_results$twm
    if (isTRUE(twm$error)) return(NULL)
    text_sz <- input$tox_plot_text_size %||% 14
    tbl <- twm$twm_table
    tbl$Treatment <- factor(tbl$Treatment, levels = tbl$Treatment[order(tbl$TWM)])
    ggplot2::ggplot(tbl, ggplot2::aes(.data[["Treatment"]], .data[["TWM"]], fill = .data[["Treatment"]])) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["TWM"]], 2)), vjust = -0.5, size = 4) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = "Therapeutic Window Metric by Treatment",
                   subtitle = "TWM = TGI(%) / Max Weight Loss(%). Dashed line = TWM of 1.",
                   x = "Treatment", y = "TWM") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_scatter_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bivariate)
    biv <- rv$tox_results$bivariate
    if (isTRUE(biv$error)) return(NULL)
    pm        <- biv$per_mouse
    text_sz   <- input$tox_plot_text_size %||% 14
    pt_size   <- input$tox_et_scatter_pt_size %||% 3
    show_ell  <- isTRUE(input$tox_et_show_ellipses)
    show_cent <- isTRUE(input$tox_et_show_centroids)
    eff_label <- switch(biv$efficacy_metric,
      tgi = "TGI (%)", tumor_auc = "Tumor AUC", log_cell_kill = "Log Cell Kill")
    p <- ggplot2::ggplot(pm,
           ggplot2::aes(.data[["Efficacy"]], .data[["Max_Weight_Loss_Pct"]],
                        color = .data[["Treatment"]])) +
      ggplot2::geom_point(size = pt_size, alpha = 0.7)
    if (show_ell && nrow(pm) > 3)
      p <- p + ggplot2::stat_ellipse(level = 0.68, linewidth = 0.8, show.legend = FALSE)
    if (show_cent && !is.null(biv$per_group))
      p <- p + ggplot2::geom_point(data = biv$per_group,
                 ggplot2::aes(.data[["Mean_Efficacy"]], .data[["Mean_Max_Weight_Loss_Pct"]]),
                 shape = 18, size = 6)
    p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Safety-Efficacy Scatter",
                   subtitle = paste0("Each point is one mouse. Efficacy: ", eff_label),
                   x = eff_label, y = "Max Weight Loss (%)", color = "Treatment") +
      ggplot2::theme_classic(base_size = text_sz)
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_benefit_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$benefit)
    ben <- rv$tox_results$benefit
    if (isTRUE(ben$error)) return(NULL)
    text_sz <- input$tox_plot_text_size %||% 14
    bt <- ben$benefit_table
    bt$Treatment <- factor(bt$Treatment, levels = bt$Treatment[order(bt$Benefit)])
    ggplot2::ggplot(bt, ggplot2::aes(.data[["Treatment"]], .data[["Benefit"]],
                                      fill = .data[["Benefit"]] > 0)) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["Benefit"]], 1)),
                         vjust = ifelse(bt$Benefit >= 0, -0.5, 1.5), size = 4) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "#d7191c"),
                                 labels = c("TRUE" = "Net Benefit", "FALSE" = "Net Harm"), name = "") +
      ggplot2::labs(title = "Total Benefit Score by Treatment",
                   subtitle = sprintf("B = AUC(efficacy) - lambda x AUC(toxicity), lambda = %.1f",
                                      ben$lambda),
                   x = "Treatment", y = "Benefit Score") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_wctgi_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$wctgi)
    wc <- rv$tox_results$wctgi
    if (isTRUE(wc$error)) return(NULL)
    text_sz <- input$tox_plot_text_size %||% 14
    comp <- wc$comparison
    plot_df <- data.frame(
      Treatment = rep(comp$Treatment, 2),
      Type      = rep(c("Uncorrected TGI", "Weight-Corrected TGI"), each = nrow(comp)),
      TGI       = c(comp$Uncorrected_TGI, comp$Corrected_TGI),
      stringsAsFactors = FALSE
    )
    ggplot2::ggplot(plot_df,
      ggplot2::aes(.data[["Treatment"]], .data[["TGI"]], fill = .data[["Type"]])) +
      ggplot2::geom_col(position = ggplot2::position_dodge(0.7), width = 0.6, alpha = 0.8) +
      ggplot2::labs(title = "Weight-Corrected vs Standard TGI",
                   subtitle = sprintf("Safety threshold: %d%% weight loss",
                                      round(wc$safety_threshold * 100)),
                   x = "Treatment", y = "TGI (%)", fill = "") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_heatmap_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bivariate)
    biv <- rv$tox_results$bivariate
    if (isTRUE(biv$error)) return(NULL)
    text_sz <- input$tox_plot_text_size %||% 14
    pm <- biv$per_mouse
    md <- pm[, c("Efficacy", "Max_Weight_Loss_Pct"), drop = FALSE]
    colnames(md) <- c("Efficacy", "Max Weight Loss (%)")
    if (!is.null(rv$tox_results$benefit) && !isTRUE(rv$tox_results$benefit$error)) {
      ea <- rv$tox_results$benefit$efficacy_auc
      ta <- rv$tox_results$benefit$toxicity_auc
      if (!is.null(ea) && nrow(ea) == nrow(pm)) md[["Efficacy AUC"]] <- ea$AUC
      if (!is.null(ta) && nrow(ta) == nrow(pm)) md[["Toxicity AUC"]] <- ta$AUC
    }
    req(ncol(md) >= 2)
    cm     <- stats::cor(md, use = "pairwise.complete.obs")
    cor_df <- expand.grid(Var1 = colnames(cm), Var2 = colnames(cm), stringsAsFactors = FALSE)
    cor_df$Correlation <- as.vector(cm)
    ggplot2::ggplot(cor_df, ggplot2::aes(.data[["Var1"]], .data[["Var2"]],
                                          fill = .data[["Correlation"]])) +
      ggplot2::geom_tile(color = "white", linewidth = 1) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data[["Correlation"]])),
                         size = 5, color = "black") +
      ggplot2::scale_fill_gradient2(low = "#d7191c", mid = "white", high = "#2c7bb6",
                                    midpoint = 0, limits = c(-1, 1)) +
      ggplot2::labs(title = "Correlation Heatmap: Efficacy vs Toxicity Metrics", x = "", y = "") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }, width  = function() input$tox_plot_width  %||% 600,
     height = function() input$tox_plot_height %||% 600)

  # ---- Active plot download ----
  output$download_tox_active_plot <- downloadHandler(
    filename = function() {
      if ((input$tox_model_type %||% "bw") == "bw") {
        switch(input$tox_bw_plot_selector %||% "trajectory",
          trajectory = "body_weight_trajectory.png", auc = "body_weight_auc.png",
          km = "weight_loss_km.png", "toxicity_plot.png")
      } else {
        switch(input$tox_et_plot_selector %||% "twm",
          twm = "therapeutic_window.png", scatter = "safety_efficacy_scatter.png",
          benefit = "total_benefit.png", wctgi = "weight_corrected_tgi.png",
          heatmap = "correlation_heatmap.png", "efficacy_toxicity_plot.png")
      }
    },
    content = function(file) {
      w   <- input$tox_plot_width  %||% 800
      h   <- input$tox_plot_height %||% 500
      dpi <- input$tox_plot_dpi    %||% 100
      model <- input$tox_model_type %||% "bw"
      if (model == "bw" && (input$tox_bw_plot_selector %||% "") == "km") {
        req(rv$tox_results, rv$tox_results$weight_loss)
        wl <- rv$tox_results$weight_loss
        if (isTRUE(wl$error)) return()
        png(file, width = w, height = h, units = "px", res = dpi)
        km_p <- survminer::ggsurvplot(wl$km_fit, data = wl$event_data,
          risk.table = TRUE, pval = TRUE, conf.int = TRUE,
          title = paste0("Time to ", round(wl$threshold*100), "% Weight Loss"),
          ggtheme = ggplot2::theme_classic())
        print(km_p); dev.off(); return()
      }
      plot_id <- if (model == "bw") {
        switch(input$tox_bw_plot_selector %||% "trajectory",
          trajectory = "tox_trajectory_plot", auc = "tox_auc_plot", "tox_trajectory_plot")
      } else {
        switch(input$tox_et_plot_selector %||% "twm",
          twm = "tox_et_twm_plot", scatter = "tox_et_scatter_plot",
          benefit = "tox_et_benefit_plot", wctgi = "tox_et_wctgi_plot",
          heatmap = "tox_et_heatmap_plot", "tox_et_twm_plot")
      }
      tryCatch({
        png(file, width = w, height = h, units = "px", res = dpi)
        print(output[[plot_id]]())
        dev.off()
      }, error = function(e) dev.off())
    }
  )

  # ---- CSV downloads ----
  output$download_tox_model <- downloadHandler("body_weight_model.csv", function(file) {
    req(rv$tox_results, rv$tox_results$mixed_model)
    if (isTRUE(rv$tox_results$mixed_model$error)) return()
    utils::write.csv(rv$tox_results$mixed_model$fixed_effects, file, row.names = FALSE)
  })
  output$download_tox_auc_csv <- downloadHandler("body_weight_auc.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bw_auc)
    if (isTRUE(rv$tox_results$bw_auc$error)) return()
    utils::write.csv(rv$tox_results$bw_auc$auc_per_mouse, file, row.names = FALSE)
  })
  output$download_tox_events <- downloadHandler("weight_loss_events.csv", function(file) {
    req(rv$tox_results, rv$tox_results$weight_loss)
    if (isTRUE(rv$tox_results$weight_loss$error)) return()
    utils::write.csv(rv$tox_results$weight_loss$event_data, file, row.names = FALSE)
  })
  output$download_et_twm_csv <- downloadHandler("therapeutic_window.csv", function(file) {
    req(rv$tox_results, rv$tox_results$twm)
    if (isTRUE(rv$tox_results$twm$error)) return()
    utils::write.csv(rv$tox_results$twm$twm_table, file, row.names = FALSE)
  })
  output$download_et_bivariate_csv <- downloadHandler("safety_efficacy_per_mouse.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bivariate)
    if (isTRUE(rv$tox_results$bivariate$error)) return()
    utils::write.csv(rv$tox_results$bivariate$per_mouse, file, row.names = FALSE)
  })
  output$download_et_benefit_csv <- downloadHandler("total_benefit.csv", function(file) {
    req(rv$tox_results, rv$tox_results$benefit)
    if (isTRUE(rv$tox_results$benefit$error)) return()
    utils::write.csv(rv$tox_results$benefit$benefit_table, file, row.names = FALSE)
  })
  output$download_et_wctgi_csv <- downloadHandler("weight_corrected_tgi.csv", function(file) {
    req(rv$tox_results, rv$tox_results$wctgi)
    if (isTRUE(rv$tox_results$wctgi$error)) return()
    utils::write.csv(rv$tox_results$wctgi$comparison, file, row.names = FALSE)
  })
  output$download_et_report <- downloadHandler("efficacy_toxicity_report.csv", function(file) {
    req(rv$tox_results)
    res    <- rv$tox_results
    report <- NULL
    if (!isTRUE(res$twm$error))    report <- res$twm$twm_table
    if (!isTRUE(res$benefit$error) && !is.null(report))
      report <- merge(report, res$benefit$benefit_table[, c("Treatment","Benefit")],
                      by = "Treatment", all.x = TRUE)
    if (!isTRUE(res$wctgi$error)   && !is.null(report))
      report <- merge(report,
                      res$wctgi$comparison[, c("Treatment","Corrected_TGI","N_Excluded")],
                      by = "Treatment", all.x = TRUE)
    if (is.null(report)) return()
    utils::write.csv(report, file, row.names = FALSE)
  })

  }) # end moduleServer
}
