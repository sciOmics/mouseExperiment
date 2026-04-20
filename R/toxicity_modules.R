#' Toxicity Analysis UI
#' @keywords internal
toxicityUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Toxicity Analysis"),

    p("Analyze body weight changes over time using mixed-effects models, AUC analysis, ",
      "and time-to-threshold survival curves. ",
      "Requires a Weight/Mass column configured in the Data Upload tab."),

    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Analysis Parameters"),

          uiOutput(ns("tox_weight_status")),

          conditionalPanel(ns = ns,
            condition = "output.tox_weight_available == true",

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
            sliderInput(ns("tox_threshold"), "Weight loss threshold (%)",
                        min = 5, max = 40, value = 20, step = 1),
            uiOutput(ns("tox_baseline_day_ui")),

            br(),
            actionButton(ns("run_tox_analysis"),
                        "Run Toxicity Analysis",
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

          tabPanel("Body Weight Trajectory",
            tab_info(
              strong("Body weight over time "),
              "by treatment group. Shows individual mouse trajectories (spaghetti) or group ",
              "mean \u00b1 error. Tumor-adjusted net weight is shown when enabled."
            ),
            fluidRow(
              column(3,
                wellPanel(
                  radioButtons(ns("tox_traj_type"), "Display",
                    choices = c("Mean \u00b1 SEM" = "mean_sem",
                                "Mean \u00b1 SD" = "mean_sd",
                                "Individual mice" = "spaghetti"),
                    selected = "mean_sem"
                  ),
                  hr(),
                  h6(icon("palette"), " Appearance"),
                  sliderInput(ns("tox_traj_text_size"), "Text size (pt)",
                              min = 10, max = 24, value = 14, step = 1),
                  sliderInput(ns("tox_traj_line_width"), "Line width",
                              min = 0.5, max = 3, value = 1, step = 0.25),
                  hr(),
                  h6(icon("arrows-alt"), " Export"),
                  fluidRow(
                    column(6, numericInput(ns("tox_traj_width"), "Width", 800, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("tox_traj_height"), "Height", 500, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("tox_traj_dpi"), "DPI", 100, min = 72, max = 600, step = 10),
                  downloadButton(ns("download_tox_trajectory"), "Download Plot")
                )
              ),
              column(9,
                shinycssloaders::withSpinner(
                  plotOutput(ns("tox_trajectory_plot"), width = "auto", height = "auto")
                )
              )
            )
          ),

          tabPanel("Mixed Model",
            tab_info(
              strong("Linear mixed-effects model "),
              "of body weight over time. Treatment \u00d7 Day interaction tests whether weight ",
              "trajectories differ between groups. Random effects account for repeated measures within mice."
            ),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("tox_model_summary"))
            ),
            h4("Fixed Effects"),
            shinycssloaders::withSpinner(
              DTOutput(ns("tox_fixed_effects"))
            ),
            h4("Estimated Marginal Means by Treatment"),
            shinycssloaders::withSpinner(
              DTOutput(ns("tox_emmeans"))
            ),
            downloadButton(ns("download_tox_model"), "Download Results CSV")
          ),

          tabPanel("Body Weight AUC",
            tab_info(
              strong("Area under the curve "),
              "of body weight change (% from baseline). Negative AUC indicates net weight loss. ",
              "Includes nadir (minimum) weight analysis per mouse."
            ),
            fluidRow(
              column(3,
                wellPanel(
                  h6(icon("arrows-alt"), " Export"),
                  fluidRow(
                    column(6, numericInput(ns("tox_auc_width"), "Width", 800, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("tox_auc_height"), "Height", 500, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("tox_auc_dpi"), "DPI", 100, min = 72, max = 600, step = 10),
                  downloadButton(ns("download_tox_auc_plot"), "Download Plot")
                )
              ),
              column(9,
                shinycssloaders::withSpinner(
                  plotOutput(ns("tox_auc_plot"), width = "auto", height = "auto")
                )
              )
            ),
            hr(),
            h4("AUC Summary by Group"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_auc_summary"))),
            h4("Nadir Weight Analysis"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_nadir_table"))),
            h4("Pairwise Comparisons"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_auc_comparisons"))),
            downloadButton(ns("download_tox_auc_csv"), "Download AUC CSV")
          ),

          tabPanel("Weight Loss Survival",
            tab_info(
              strong("Time-to-threshold analysis: "),
              "Kaplan-Meier curves for time to reach the specified body weight loss threshold. ",
              "Mice that never reach the threshold are right-censored at their last observation."
            ),
            fluidRow(
              column(3,
                wellPanel(
                  h6(icon("arrows-alt"), " Export"),
                  fluidRow(
                    column(6, numericInput(ns("tox_wls_width"), "Width", 800, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("tox_wls_height"), "Height", 500, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("tox_wls_dpi"), "DPI", 100, min = 72, max = 600, step = 10),
                  downloadButton(ns("download_tox_km_plot"), "Download KM Plot")
                )
              ),
              column(9,
                shinycssloaders::withSpinner(
                  plotOutput(ns("tox_km_plot"), width = "auto", height = "auto")
                )
              )
            ),
            hr(),
            h4("Survival Summary"),
            shinycssloaders::withSpinner(verbatimTextOutput(ns("tox_km_summary"))),
            h4("Event Data"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_event_table"))),
            downloadButton(ns("download_tox_events"), "Download Events CSV")
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
  current_tox_trajectory <- reactiveVal(NULL)
  current_tox_auc_plot   <- reactiveVal(NULL)
  current_tox_km_plot    <- reactiveVal(NULL)

  # Reactive: is a weight column configured?
  tox_weight_available <- reactive({
    !is.null(rv$weight_col) && rv$weight_col != ""
  })

  output$tox_weight_available <- reactive({ tox_weight_available() })
  outputOptions(output, "tox_weight_available", suspendWhenHidden = FALSE)

  # Status banner
  output$tox_weight_status <- renderUI({
    if (!rv$data_valid) {
      return(status_alert("Please upload and configure data in the Data Upload tab first.", "warning"))
    }
    if (is.null(rv$processed_data)) {
      tryCatch(ensure_preprocessed(rv), error = function(e) NULL)
    }
    if (!tox_weight_available()) {
      return(status_alert("No Weight/Mass column selected. Please configure it in the Data Upload tab (Optional Columns \u2192 Weight / Mass).", "warning"))
    }
    status_alert(paste("Weight column:", rv$weight_col), "success")
  })

  # Reference group selector
  output$tox_reference_group_ui <- renderUI({
    req(rv$processed_data, rv$tx_col)
    groups <- sort(unique(as.character(na.omit(rv$processed_data[[rv$tx_col]]))))
    default_ctrl <- pick_reference_group(groups)
    selectInput(session$ns("tox_reference_group"),
                "Control / Reference Group",
                choices = groups,
                selected = default_ctrl)
  })

  # Baseline day selector
  output$tox_baseline_day_ui <- renderUI({
    req(rv$processed_data, rv$time_col)
    days <- sort(unique(na.omit(rv$processed_data[[rv$time_col]])))
    day_choices <- c("First measurement" = "first",
                     setNames(as.character(days), paste("Day", days)))
    selectInput(session$ns("tox_baseline_day"),
                tagList("Baseline day", tags$small("(day of first treatment is best)")),
                choices  = day_choices,
                selected = "first")
  })

  # Disable covariates that aren't available
  observe({
    has_vol <- !is.null(rv$vol_col) && rv$vol_col != ""
    has_sex <- !is.null(rv$sex_col) && rv$sex_col != ""
    available <- character(0)
    if (has_vol) available <- c(available, "volume")
    available <- c(available, "initial_mass")
    if (has_sex) available <- c(available, "sex")

    choices <- c("Tumor Volume" = "volume",
                 "Sex" = "sex",
                 "Initial Mass" = "initial_mass")
    # Only show available choices
    choices <- choices[choices %in% available]
    selected <- intersect(input$tox_covariates %||% "volume", available)
    updateCheckboxGroupInput(session, "tox_covariates",
                             choices = choices, selected = selected)
  })

  # ---- Run Analysis ----
  observeEvent(input$run_tox_analysis, {
    req(rv$processed_data, tox_weight_available())

    output$tox_status_msg <- renderUI({
      status_alert("Running toxicity analysis...", "info")
    })

    data     <- rv$processed_data
    wt_col   <- rv$weight_col
    vol_col  <- rv$vol_col
    time_col <- rv$time_col
    tx_col   <- rv$tx_col
    id_col   <- rv$id_col
    sex_col  <- rv$sex_col
    cage_col <- rv$cage_col
    ref_group <- input$tox_reference_group %||% NULL

    has_vol <- !is.null(vol_col) && vol_col %in% names(data)
    adjust_tumor <- isTRUE(input$tox_adjust_tumor) && has_vol
    density_val  <- input$tox_tumor_density %||% 1.0
    covs         <- input$tox_covariates %||% character(0)
    est_method   <- input$tox_estimation %||% "REML"
    threshold    <- (input$tox_threshold %||% 20) / 100
    baseline_day_raw <- input$tox_baseline_day %||% "first"
    baseline_day <- if (baseline_day_raw == "first") NULL else as.numeric(baseline_day_raw)

    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Running toxicity analysis...", value = 0.1)
    on.exit(progress$close(), add = TRUE)

    results <- list()

    # 1. Mixed model
    progress$set(value = 0.3, message = "Fitting mixed model...")
    results$mixed_model <- safe_analysis(module = "toxicity_mixed_model", {
      mouseExperiment::analyze_body_weight(
        df                  = data,
        weight_column       = wt_col,
        time_column         = time_col,
        treatment_column    = tx_col,
        id_column           = id_col,
        volume_column       = if (has_vol) vol_col else NULL,
        sex_column          = sex_col,
        cage_column         = cage_col,
        adjust_tumor_weight = adjust_tumor,
        tumor_density       = density_val,
        covariates          = covs,
        estimation          = est_method,
        reference_group     = ref_group
      )
    })

    # 2. Body weight AUC
    progress$set(value = 0.5, message = "Computing AUC...")
    results$bw_auc <- safe_analysis(module = "toxicity_auc", {
      mouseExperiment::body_weight_auc(
        df                  = data,
        weight_column       = wt_col,
        time_column         = time_col,
        treatment_column    = tx_col,
        id_column           = id_col,
        volume_column       = if (has_vol) vol_col else NULL,
        adjust_tumor_weight = adjust_tumor,
        tumor_density       = density_val,
        reference_group     = ref_group
      )
    })

    # 3. Weight loss threshold (KM)
    progress$set(value = 0.7, message = "Time-to-threshold analysis...")
    results$weight_loss <- safe_analysis(module = "toxicity_threshold", {
      mouseExperiment::weight_loss_threshold(
        df                  = data,
        weight_column       = wt_col,
        time_column         = time_col,
        treatment_column    = tx_col,
        id_column           = id_col,
        volume_column       = if (has_vol) vol_col else NULL,
        adjust_tumor_weight = adjust_tumor,
        tumor_density       = density_val,
        threshold           = threshold,
        baseline_day        = baseline_day,
        reference_group     = ref_group
      )
    })

    # Store trajectory data for plotting
    results$trajectory_data <- tryCatch({
      wd <- data.frame(
        ID        = as.character(data[[id_col]]),
        Treatment = as.character(data[[tx_col]]),
        Day       = as.numeric(data[[time_col]]),
        Weight    = as.numeric(data[[wt_col]]),
        stringsAsFactors = FALSE
      )
      if (has_vol && adjust_tumor) {
        vol <- as.numeric(data[[vol_col]])
        wd$Weight <- wd$Weight - (vol / 1000 * density_val)
      }
      wd <- wd[!is.na(wd$Weight) & !is.na(wd$Day), ]
      wd
    }, error = function(e) NULL)

    rv$tox_results <- results

    progress$set(value = 1.0, message = "Done!")

    n_errors <- sum(vapply(results[c("mixed_model", "bw_auc", "weight_loss")],
                           function(x) isTRUE(x$error), logical(1)))
    if (n_errors == 3) {
      output$tox_status_msg <- renderUI({
        status_alert("All analyses failed. Check your data and settings.", "danger")
      })
    } else if (n_errors > 0) {
      output$tox_status_msg <- renderUI({
        status_alert(sprintf("Completed with %d warning(s). Some analyses may have failed.", n_errors), "warning")
      })
    } else {
      output$tox_status_msg <- renderUI({
        status_alert("Toxicity analysis completed successfully!", "success")
      })
    }
  })

  # ---- Body Weight Trajectory Plot ----
  output$tox_trajectory_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$trajectory_data)
    td <- rv$tox_results$trajectory_data
    req(nrow(td) > 0)

    traj_type <- input$tox_traj_type %||% "mean_sem"
    text_sz   <- input$tox_traj_text_size %||% 14
    line_w    <- input$tox_traj_line_width %||% 1

    if (traj_type == "spaghetti") {
      p <- ggplot2::ggplot(td, ggplot2::aes(x = .data[["Day"]], y = .data[["Weight"]],
                                             group = .data[["ID"]], color = .data[["Treatment"]])) +
        ggplot2::geom_line(alpha = 0.6, linewidth = line_w) +
        ggplot2::geom_point(size = 1.5, alpha = 0.6) +
        ggplot2::labs(title = "Body Weight Trajectories",
                     x = "Day", y = "Weight (g)", color = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
    } else {
      # Compute group summaries
      agg <- stats::aggregate(Weight ~ Treatment + Day, data = td, FUN = function(x) {
        c(mean = mean(x, na.rm = TRUE),
          sd   = stats::sd(x, na.rm = TRUE),
          n    = length(x))
      })
      agg_flat <- do.call(data.frame, agg)
      names(agg_flat) <- c("Treatment", "Day", "Mean", "SD", "N")
      agg_flat$SEM <- agg_flat$SD / sqrt(agg_flat$N)
      agg_flat$Error <- if (traj_type == "mean_sd") agg_flat$SD else agg_flat$SEM

      p <- ggplot2::ggplot(agg_flat, ggplot2::aes(x = .data[["Day"]], y = .data[["Mean"]],
                                                    color = .data[["Treatment"]],
                                                    fill  = .data[["Treatment"]])) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["Mean"]] - .data[["Error"]],
                                          ymax = .data[["Mean"]] + .data[["Error"]]),
                             alpha = 0.2, color = NA) +
        ggplot2::geom_line(linewidth = line_w) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(title = paste0("Body Weight Over Time (Mean \u00b1 ",
                                     if (traj_type == "mean_sd") "SD" else "SEM", ")"),
                     x = "Day", y = "Weight (g)", color = "Treatment", fill = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
    }

    current_tox_trajectory(p)
    p
  }, width = function() input$tox_traj_width %||% 800,
     height = function() input$tox_traj_height %||% 500)

  # ---- Mixed Model Outputs ----
  output$tox_model_summary <- renderText({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error)) return(paste("Error:", mm$message))
    mm$summary_text
  })

  output$tox_fixed_effects <- renderDT({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error)) return(datatable(data.frame(Message = mm$message)))
    tbl <- format_df_for_display(mm$fixed_effects)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  output$tox_emmeans <- renderDT({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$emmeans_table)) {
      return(datatable(data.frame(Message = "Not available")))
    }
    tbl <- format_df_for_display(mm$emmeans_table)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ---- AUC Outputs ----
  output$tox_auc_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error)) return(NULL)

    auc_df <- auc$auc_per_mouse
    p <- ggplot2::ggplot(auc_df, ggplot2::aes(x = .data[["Treatment"]],
                                               y = .data[["AUC_Pct_Change"]],
                                               fill = .data[["Treatment"]])) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Body Weight AUC (% Change from Baseline)",
                   x = "Treatment", y = "AUC of % Weight Change") +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(legend.position = "none")

    current_tox_auc_plot(p)
    p
  }, width = function() input$tox_auc_width %||% 800,
     height = function() input$tox_auc_height %||% 500)

  output$tox_auc_summary <- renderDT({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error)) return(datatable(data.frame(Message = auc$message)))
    tbl <- format_df_for_display(auc$auc_summary)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  output$tox_nadir_table <- renderDT({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error)) return(datatable(data.frame(Message = "Not available")))
    tbl <- format_df_for_display(auc$nadir_data)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  output$tox_auc_comparisons <- renderDT({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error) || is.null(auc$comparisons)) {
      return(datatable(data.frame(Message = "Not available")))
    }
    tbl <- format_df_for_display(auc$comparisons)
    datatable(tbl, options = list(scrollX = TRUE), rownames = FALSE)
  })

  # ---- Weight Loss KM ----
  output$tox_km_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(NULL)

    threshold_pct <- round(wl$threshold * 100)
    km_plot <- tryCatch({
      survminer::ggsurvplot(
        wl$km_fit,
        data       = wl$event_data,
        risk.table = TRUE,
        pval       = TRUE,
        conf.int   = TRUE,
        title      = paste0("Time to ", threshold_pct, "% Body Weight Loss"),
        xlab       = "Day",
        ylab       = "Proportion Without Threshold Weight Loss",
        ggtheme    = ggplot2::theme_classic(base_size = 14)
      )
    }, error = function(e) NULL)

    if (!is.null(km_plot)) {
      current_tox_km_plot(km_plot)
      print(km_plot)
    }
  }, width = function() input$tox_wls_width %||% 800,
     height = function() input$tox_wls_height %||% 500)

  output$tox_km_summary <- renderText({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(paste("Error:", wl$message))

    lines <- c(
      sprintf("=== WEIGHT LOSS SURVIVAL (Threshold: %d%%) ===", round(wl$threshold * 100)),
      ""
    )

    if (!is.null(wl$log_rank)) {
      lines <- c(lines,
        sprintf("Log-rank test p-value: %.4g", wl$log_rank$p.value %||% NA),
        "")
    }

    if (!is.null(wl$cox_summary)) {
      lines <- c(lines, "-- Cox PH Model --")
      cox_lines <- utils::capture.output(print(wl$cox_summary))
      lines <- c(lines, cox_lines)
    }

    paste(lines, collapse = "\n")
  })

  output$tox_event_table <- renderDT({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(datatable(data.frame(Message = wl$message)))
    tbl <- format_df_for_display(wl$event_data)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  # ---- Download Handlers ----
  output$download_tox_trajectory <- downloadHandler(
    filename = "body_weight_trajectory.png",
    content  = function(file) {
      p <- current_tox_trajectory()
      req(!is.null(p))
      png(file, width = input$tox_traj_width %||% 800,
                height = input$tox_traj_height %||% 500,
                units = "px", res = input$tox_traj_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_tox_model <- downloadHandler(
    filename = "toxicity_mixed_model.csv",
    content  = function(file) {
      req(rv$tox_results, rv$tox_results$mixed_model)
      mm <- rv$tox_results$mixed_model
      if (isTRUE(mm$error)) return()
      utils::write.csv(mm$fixed_effects, file, row.names = FALSE)
    }
  )

  output$download_tox_auc_plot <- downloadHandler(
    filename = "body_weight_auc.png",
    content  = function(file) {
      p <- current_tox_auc_plot()
      req(!is.null(p))
      png(file, width = input$tox_auc_width %||% 800,
                height = input$tox_auc_height %||% 500,
                units = "px", res = input$tox_auc_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_tox_auc_csv <- downloadHandler(
    filename = "body_weight_auc.csv",
    content  = function(file) {
      req(rv$tox_results, rv$tox_results$bw_auc)
      auc <- rv$tox_results$bw_auc
      if (isTRUE(auc$error)) return()
      utils::write.csv(auc$auc_per_mouse, file, row.names = FALSE)
    }
  )

  output$download_tox_km_plot <- downloadHandler(
    filename = "weight_loss_km.png",
    content  = function(file) {
      req(rv$tox_results, rv$tox_results$weight_loss)
      wl <- rv$tox_results$weight_loss
      if (isTRUE(wl$error)) return()
      png(file, width = input$tox_wls_width %||% 800,
                height = input$tox_wls_height %||% 500,
                units = "px", res = input$tox_wls_dpi %||% 100)
      km_plot <- survminer::ggsurvplot(wl$km_fit, data = wl$event_data,
        risk.table = TRUE, pval = TRUE, conf.int = TRUE,
        title = paste0("Time to ", round(wl$threshold * 100), "% Weight Loss"),
        ggtheme = ggplot2::theme_classic())
      print(km_plot); dev.off()
    }
  )

  output$download_tox_events <- downloadHandler(
    filename = "weight_loss_events.csv",
    content  = function(file) {
      req(rv$tox_results, rv$tox_results$weight_loss)
      wl <- rv$tox_results$weight_loss
      if (isTRUE(wl$error)) return()
      utils::write.csv(wl$event_data, file, row.names = FALSE)
    }
  )

  }) # end moduleServer
}
