#' Dose Response Analysis UI
#' @keywords internal
doseResponseUI <- function(id) {
  ns <- NS(id)
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

          uiOutput(ns("dr_dose_status")),

          conditionalPanel(ns = ns,
            condition = "output.dr_dose_available == true",

            selectInput(ns("dr_control_group"),
                       "Reference Group",
                       choices = "",
                       selected = ""),

            selectInput(ns("dr_model_type"), "Statistical Model",
              choices = c(
                "Frequentist (ANOVA/regression)" = "frequentist",
                "Bayesian Hill/Emax (brms)"      = "bayes"
              ),
              selected = "frequentist"
            ),
            helpText("Frequentist fits ANOVA with Tukey HSD and linear/non-linear ",
                     "regression. Bayesian fits a Hill/Emax concentration-response model ",
                     "with brms; computation takes 3–12 minutes."),

            conditionalPanel(ns = ns, "input.dr_model_type == 'frequentist'",
              selectInput(ns("dr_comparison_type"),
                         "Pairwise Comparison Type",
                         choices = c(
                           "vs. Reference Group" = "reference",
                           "All Pairwise (Tukey HSD)" = "all"
                         ),
                         selected = "reference")
            ),

            conditionalPanel(ns = ns, "input.dr_model_type == 'bayes'",
              hr(),
              h5(strong("Bayesian Settings")),
              selectInput(ns("dr_prior_strength"), "Prior Strength",
                choices = c(
                  "Skeptical (recommended)" = "skeptical",
                  "Weakly informative"      = "weakly_informative",
                  "Informative"             = "informative",
                  "Diffuse"                 = "diffuse"
                ),
                selected = "skeptical"
              ),
              numericInput(ns("dr_bayes_chains"), "MCMC chains",
                           value = 4, min = 1, max = 8, step = 1),
              numericInput(ns("dr_bayes_iter"), "Iterations per chain",
                           value = 2000, min = 500, max = 10000, step = 500),
              div(class = "alert alert-warning",
                  style = "font-size:0.85em; padding:8px; margin-top:6px;",
                icon("clock"), " Bayesian fitting takes ~3–12 minutes. ",
                "The app stays usable while the model runs in the background."
              )
            ),

            uiOutput(ns("dr_time_point_ui")),

            br(),
            actionButton(ns("run_dr_analysis"),
                        "Run Dose-Response Analysis",
                        class = "btn-lg btn-primary btn-block")
          ),

          br(),
          uiOutput(ns("dr_status_msg"))
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
          id = ns("dr_results_tabs"),

          tabPanel("Info",
            fluidRow(
              column(10, offset = 1,
                br(),

                # ── Overview ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Overview")),
                  div(class = "panel-body",
                    p("The Dose-Response tab models how tumor volume (or TGI) varies across dose levels. Three complementary analyses are performed simultaneously:"),
                    tags$ul(
                      tags$li(strong("Linear regression"), " — tests whether there is a significant linear trend in tumor volume with increasing dose."),
                      tags$li(strong("One-way ANOVA"), " — tests whether any dose group differs from the others, followed by Tukey HSD post-hoc comparisons if significant."),
                      tags$li(strong("Non-linear dose-response model (4PL)"), " — fits a four-parameter log-logistic (Hill) model to quantify ED50 and curve steepness.")
                    ),
                    p(tags$small(class = "text-muted",
                      "Tumor Growth Analysis must be run before using this tab — dose-response uses volume data preprocessed there. ",
                      "A Dose column must also be mapped in Data Upload."
                    ))
                  )
                ),

                # ── Linear Regression ────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-info", "Linear Regression"))),
                  div(class = "panel-body",
                    p("Fits ", tags$code("Volume ~ Dose"), " by ordinary least squares. The slope quantifies the change in mean tumor volume per unit increase in dose."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Slope"),
                      tags$dd("Change in volume per dose unit. A negative slope indicates dose-dependent tumor suppression."),
                      tags$dt("R²"),
                      tags$dd("Proportion of variance in volume explained by dose. Low R² with a significant slope can occur when individual variability is high relative to group differences."),
                      tags$dt("p-value"),
                      tags$dd("Significance of the slope. A non-significant p-value does not exclude a non-linear relationship — inspect the 4PL model output and plots as well.")
                    )
                  )
                ),

                # ── ANOVA ────────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "One-Way ANOVA & Tukey HSD"))),
                  div(class = "panel-body",
                    p("Treats dose as a categorical factor and tests whether group means differ. A significant F-test (p < 0.05) triggers Tukey HSD post-hoc comparisons to identify which specific dose pairs differ."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("F-statistic"),
                      tags$dd("Ratio of between-group to within-group variance. Larger values indicate greater separation between dose groups relative to within-group noise."),
                      tags$dt("Tukey HSD"),
                      tags$dd("Family-wise error rate corrected pairwise differences. Values in 'p adj' below 0.05 indicate a significant difference at the 5 % level after correction for all comparisons.")
                    ),
                    p(tags$small(class = "text-muted",
                      "ANOVA assumes equal variances across groups (homoscedasticity). ",
                      "When this is violated, interpret alongside the linear regression, which is more robust to mild heteroscedasticity."
                    ))
                  )
                ),

                # ── 4PL ──────────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-success", "Non-Linear Model (4-Parameter Log-Logistic)"))),
                  div(class = "panel-body",
                    p("A four-parameter log-logistic (4PL / Hill) model is fitted via the ", tags$code("drc"), " package:"),
                    tags$blockquote(tags$code(
                      "Volume = c + (d − c) / (1 + (Dose / ED50)^b)"
                    )),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("b"), tags$dd("Hill coefficient: steepness of the dose-response curve. Larger absolute values indicate a sharper transition."),
                      tags$dt("c"), tags$dd("Lower asymptote (response at saturating inhibitory dose)."),
                      tags$dt("d"), tags$dd("Upper asymptote (response in the absence of drug)."),
                      tags$dt("ED50"), tags$dd("Dose producing 50 % of the maximal effect. The primary pharmacodynamic summary statistic.")
                    ),
                    div(class = "alert alert-info", style = "margin-top:10px;",
                      icon("circle-info"), " ",
                      strong("AIC / BIC note: "),
                      "AIC and BIC from the linear model (", tags$code("stats::lm"), ") and the non-linear model (",
                      tags$code("drc::drm"), ") are ", strong("not directly comparable"), " for model selection — ",
                      "the two packages use different likelihood parameterisations. ",
                      "Prefer visual inspection of the fitted curve and residuals when choosing between linear and non-linear models."
                    )
                  )
                ),

                br()
              )
            )
          ),

          tabPanel("Summary",
            tab_info(
              strong("Dose-response summary: "),
              "Results of linear regression (slope, R\u00b2, p-value) and one-way ANOVA across dose levels. ",
              "A significant linear slope indicates a systematic increase or decrease in tumor volume with dose. ",
              "A significant ANOVA F-test indicates that at least one dose group differs from the others."
            ),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("dr_summary"))
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
            div(class = "text-right", style = "margin-bottom: 6px;",
              downloadButton(ns("download_dr_table"), "Download CSV")
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("dr_summary_table"))
            )
          ),

          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  selectInput(ns("dr_plot_selector"), "Select Plot",
                    choices = c(
                      "Tumor Volume vs Dose" = "scatter",
                      "Tumor Growth Inhibition vs Dose" = "tgi"
                    )
                  ),
                  hr(),
                  h6(icon("palette"), " Appearance"),
                  colourpicker::colourInput(ns("dr_scatter_pt_color"),
                                           "Point colour",
                                           value = "#555555"),
                  colourpicker::colourInput(ns("dr_scatter_line_color"),
                                           "Trend line colour",
                                           value = "#2166AC"),
                  sliderInput(ns("dr_scatter_pt_size"), "Point size",
                              min = 0.5, max = 6, value = 2, step = 0.5),
                  checkboxInput(ns("dr_scatter_show_ci"),
                                "Show trend line confidence interval",
                                value = TRUE),
                  radioButtons(ns("dr_trend_method"), "Trend line fit",
                    choices = c("Linear" = "linear", "Non-linear (LOESS)" = "loess"),
                    selected = "linear", inline = TRUE
                  ),
                  hr(),
                  h6(style = "font-size: 1.15rem; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                     icon("sliders-h"), " Options"),
                  sliderInput(ns("dr_opt_text_size"), "Text size (pt)",
                              min = 10, max = 24, value = 16, step = 1),
                  hr(),
                  h6(style = "font-size: 1.15rem; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                     icon("heading"), " Title & Subtitle"),
                  textInput(ns("dr_plot_title"), "Title",
                            placeholder = "Leave blank for default"),
                  selectInput(ns("dr_plot_title_align"), "Title alignment",
                    choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                    selected = "center"),
                  textInput(ns("dr_plot_subtitle"), "Subtitle",
                            placeholder = "Optional subtitle"),
                  hr(),
                  h6(style = "font-size: 1.15rem; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                     icon("arrows-alt"), " Plot Size / Export"),
                  fluidRow(
                    column(6, numericInput(ns("dr_plot_width"),  "Width (px)",  800, min = 200, max = 5000, step = 50)),
                    column(6, numericInput(ns("dr_plot_height"), "Height (px)", 600, min = 200, max = 5000, step = 50))
                  ),
                  numericInput(ns("dr_plot_dpi"), "Resolution (DPI)", 100, min = 72, max = 600, step = 10),
                  tags$div(style = "margin-top: 10px;", uiOutput(ns("dr_download_plot_btn")))
                )
              ),
              column(9,
                conditionalPanel(ns = ns,"input.dr_plot_selector == 'scatter'",
                  h4("Tumor Volume vs Dose"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("dr_plot_scatter"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.dr_plot_selector == 'tgi'",
                  h4("Tumor Growth Inhibition vs Dose"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("dr_plot_tgi"), width = "auto", height = "auto")
                  )
                )
              )
            )
          ),

          tabPanel("Model Results",
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
            div(class = "text-right", style = "margin-bottom: 6px;",
              downloadButton(ns("download_dr_stats"), "Download Results CSV")
            ),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("dr_stat_tests"))
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("dr_tukey_table"))
            )
          )
        )
      )
    )
  )
}

#' Dose Response Server
#' @keywords internal
dose_response_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

  current_dr_scatter <- reactiveVal(NULL)
  current_dr_tgi     <- reactiveVal(NULL)

  rv$dr_results <- NULL

  # ── Bayesian DR tab management ─────────────────────────────────────────────
  bayes_dr_tabs_inserted <- reactiveVal(FALSE)

  observeEvent(input$dr_model_type, {
    tabs_id <- "dr_results_tabs"
    if (isTRUE(input$dr_model_type == "bayes") && !isTRUE(bayes_dr_tabs_inserted())) {
      appendTab(inputId = tabs_id,
        tabPanel("Posterior Check",
          tab_info(
            strong("Posterior predictive check (PPC)"), " — model-simulated dose-response overlaid on observed data.",
            em("Good fit: simulated densities envelope the observed density.")
          ),
          div(class = "text-right", style = "margin-bottom: 6px;",
            downloadButton(session$ns("download_dr_pp_check"), "Download Plot")
          ),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("dr_pp_check_plot"), height = "450px")
          )
        ),
        select = FALSE
      )
      appendTab(inputId = tabs_id,
        tabPanel("MCMC Diagnostics",
          tab_info(
            strong("MCMC convergence diagnostics. "),
            "All Rhat values should be ≤ 1.01 before interpreting results."
          ),
          shinycssloaders::withSpinner(
            DTOutput(session$ns("dr_mcmc_diag_table"))
          ),
          br(),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("dr_mcmc_trace_plot"), height = "500px")
          )
        ),
        select = FALSE
      )
      appendTab(inputId = tabs_id,
        tabPanel("Prior vs Posterior",
          tab_info(
            strong("Prior vs Posterior density overlay"), " — shows how much the data have updated the model's prior distributions."
          ),
          div(class = "text-right", style = "margin-bottom: 6px;",
            downloadButton(session$ns("download_dr_prior_posterior"), "Download Plot")
          ),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("dr_prior_posterior_plot"), height = "400px")
          )
        ),
        select = FALSE
      )
      bayes_dr_tabs_inserted(TRUE)
    } else if (!isTRUE(input$dr_model_type == "bayes") && isTRUE(bayes_dr_tabs_inserted())) {
      removeTab(inputId = tabs_id, target = "Posterior Check")
      removeTab(inputId = tabs_id, target = "MCMC Diagnostics")
      removeTab(inputId = tabs_id, target = "Prior vs Posterior")
      bayes_dr_tabs_inserted(FALSE)
    }
  }, ignoreNULL = FALSE)

  output$dr_pp_check_plot <- renderPlot({
    req(rv$dr_results)
    req(isTRUE(rv$dr_results$model_type_used == "bayes_dr"))
    req(!isTRUE(rv$dr_results$error), !is.null(rv$bayes_models$dr))
    brms::pp_check(rv$bayes_models$dr, type = "dens_overlay", ndraws = 50)
  })

  output$dr_mcmc_diag_table <- renderDT({
    req(rv$dr_results)
    req(isTRUE(rv$dr_results$model_type_used == "bayes_dr"))
    req(!isTRUE(rv$dr_results$error))
    diag <- rv$dr_results$mcmc_diagnostics
    req(!is.null(diag))
    datatable(format_df_for_display(as.data.frame(diag)),
              options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })

  output$dr_mcmc_trace_plot <- renderPlot({
    req(rv$dr_results)
    req(isTRUE(rv$dr_results$model_type_used == "bayes_dr"))
    req(!isTRUE(rv$dr_results$error), !is.null(rv$bayes_models$dr))
    bayesplot::mcmc_trace(brms::as.array(rv$bayes_models$dr))
  })

  output$dr_prior_posterior_plot <- renderPlot({
    req(rv$dr_results)
    req(isTRUE(rv$dr_results$model_type_used == "bayes_dr"))
    req(!isTRUE(rv$dr_results$error))
    p <- rv$dr_results$prior_posterior_plot
    req(!is.null(p))
    print(p)
  })

  output$download_dr_pp_check <- downloadHandler(
    filename = function() paste0("dr_posterior_check_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$bayes_models$dr))
      p <- brms::pp_check(rv$bayes_models$dr, type = "dens_overlay", ndraws = 50)
      ggplot2::ggsave(file, plot = p, width = 8, height = 5, dpi = 150)
    }
  )

  output$download_dr_prior_posterior <- downloadHandler(
    filename = function() paste0("dr_prior_posterior_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$dr_results$prior_posterior_plot))
      ggplot2::ggsave(file, plot = rv$dr_results$prior_posterior_plot,
                      width = 8, height = 5, dpi = 150)
    }
  )

  observe({
    toggle_downloads(
      session,
      c("download_dr_table", "download_dr_stats"),
      enabled = !is.null(rv$dr_results) && !isTRUE(rv$dr_results$error)
    )
  })

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
      return(status_alert("Please upload and configure data in the Data Upload tab first.", "warning"))
    }
    if (is.null(rv$processed_data)) {
      # Attempt to preprocess inline so TG analysis isn't required first
      pp_chk <- tryCatch(ensure_preprocessed(rv), error = function(e) NULL)
      if (!is.null(pp_chk) && !is.null(pp_chk$date_fmt_warning)) {
        showNotification(pp_chk$date_fmt_warning, type = "warning", duration = 8)
      }
    }
    if (is.null(rv$processed_data)) {
      return(status_alert("Could not resolve Volume and Day columns. Please check your column mappings in the Data Upload tab.", "warning"))
    }
    if (!dr_dose_available()) {
      return(status_alert("No Dose column selected. Please configure it in the Data Upload tab (Optional Columns \u2192 Dose Level).", "warning"))
    }
    status_alert(paste("Dose column:", rv$dose_col), "success")
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
    selectInput(session$ns("dr_time_point"),
                "Time Point for Analysis",
                choices  = day_choices,
                selected = "all")
  })

  observeEvent(input$run_dr_analysis, {
    req(rv$processed_data, dr_dose_available())

    output$dr_status_msg <- renderUI({
      status_alert("Running dose-response analysis...", "info")
    })

    # Track analysis start with Umami
    shinyjs::runjs("window.trackAnalysisStart('dose_response', 'linear_regression+anova', 'control');")

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
        status_alert(paste("Dose column", dose_col, "contains no numeric values. Please check your data."), "danger")
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

    # ── Bayesian async path ────────────────────────────────────────────────────
    if (isTRUE(input$dr_model_type == "bayes")) {
      if (!all(c(requireNamespace("brms",     quietly = TRUE),
                 requireNamespace("future",   quietly = TRUE),
                 requireNamespace("promises", quietly = TRUE)))) {
        output$dr_status_msg <- renderUI({
          status_alert(paste0("Bayesian mode requires 'brms', 'future', and 'promises'. ",
                              "Install with: install.packages(c('brms', 'future', 'promises'))"),
                       "danger")
        })
        return()
      }

      dr_bayes_params <- list(
        data          = data,
        dose_col      = dose_col,
        vol_col       = vol_col,
        tx_col        = tx_col,
        time_col      = time_col,
        id_col        = id_col,
        time_point    = time_point_val,
        ref_group     = input$dr_control_group,
        prior_str     = input$dr_prior_strength %||% "skeptical",
        n_chains      = as.integer(input$dr_bayes_chains %||% 4L),
        n_iter        = as.integer(input$dr_bayes_iter   %||% 2000L)
      )

      progress_dr <- shiny::Progress$new(session, min = 0, max = 1)
      progress_dr$set(message = "Starting Bayesian dose-response fit…", value = 0.05)

      output$dr_status_msg <- renderUI({
        tagList(
          h4("Bayesian model running…",
             style = "margin-top: 18px; font-weight: 600;"),
          p("Fitting Bayesian Hill/Emax model with brms. This may take 3–12 minutes.",
            style = "color:#666;"),
          shinycssloaders::withSpinner(div())
        )
      })

      shinyjs::runjs("window.trackAnalysisStart('dose_response', 'bayes', 'hill_emax');")

      promises::future_promise({
        mouseExperiment::bayesian_dose_response(
          df               = dr_bayes_params$data,
          dose_column      = dr_bayes_params$dose_col,
          volume_column    = dr_bayes_params$vol_col,
          treatment_column = dr_bayes_params$tx_col,
          day_column       = dr_bayes_params$time_col,
          id_column        = dr_bayes_params$id_col,
          endpoint_day     = dr_bayes_params$time_point,
          reference_group  = dr_bayes_params$ref_group,
          prior_strength   = dr_bayes_params$prior_str,
          n_chains         = dr_bayes_params$n_chains,
          n_iter           = dr_bayes_params$n_iter,
          return_model     = TRUE,
          plots            = TRUE
        )
      }) %...>% (function(result) {
        progress_dr$close()
        result$model_type_used <- "bayes_dr"
        if (!is.null(rv$bayes_models)) rv$bayes_models$dr <- result$model
        else                           rv$bayes_models <- list(dr = result$model)
        result$model <- NULL
        gc()
        rv$dr_results <- result
        if (isTRUE(result$error)) {
          err_msg <- result$message
          shinyjs::runjs(sprintf("window.trackAnalysisComplete('dose_response', false, %s);",
                                 js_string(err_msg)))
          output$dr_status_msg <- renderUI({
            status_alert(paste("Bayesian DR failed:", err_msg), "danger")
          })
        } else {
          shinyjs::runjs("window.trackAnalysisComplete('dose_response', true);")
          output$dr_status_msg <- renderUI({
            status_alert("Bayesian dose-response analysis completed!", "success")
          })
        }
      }) %...!% (function(err) {
        progress_dr$close()
        err_msg <- translate_error(conditionMessage(err))
        shinyjs::runjs(sprintf("window.trackAnalysisComplete('dose_response', false, %s);",
                               js_string(err_msg)))
        output$dr_status_msg <- renderUI({
          status_alert(paste("Bayesian DR failed:", err_msg), "danger")
        })
      })
      return()
    }

    # ── Frequentist path ───────────────────────────────────────────────────────
    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Validating data...", value = 0.1)
    on.exit(progress$close(), add = TRUE)

    result <- safe_analysis(module = "dose_response", {
      progress$set(value = 0.6, message = "Fitting dose-response model...")
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
      }, error = function(e) {
        message("[DoseResponse] Data preparation failed: ", conditionMessage(e))
        result$dr_prep_error <<- conditionMessage(e)
      })
    }

    progress$set(value = 0.95, message = "Preparing results...")
    rv$dr_results <- result

    # Track analysis completion with Umami
    if (isTRUE(result$error)) {
      shinyjs::runjs(sprintf("window.trackAnalysisComplete('dose_response', false, '%s');",
                             gsub("'", "\\\\'", result$message)))
    } else {
      shinyjs::runjs("window.trackAnalysisComplete('dose_response', true);")
    }

    if (isTRUE(result$error)) {
      output$dr_status_msg <- renderUI({
        status_alert(paste0("Error: ", result$message), "danger")
      })
    } else {
      output$dr_status_msg <- renderUI({
        status_alert("Dose-response analysis completed successfully!", "success")
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
      tgi     = "download_dr_tgi",
      "download_dr_scatter"
    )
    downloadButton(session$ns(btn_id), "Download Plot")
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
    pt_sz    <- input$dr_scatter_pt_size     %||% 2
    show_ci  <- isTRUE(input$dr_scatter_show_ci)
    trend    <- input$dr_trend_method %||% "linear"
    text_sz  <- input$dr_opt_text_size %||% 16
    smooth_args <- if (trend == "loess") {
      list(method = "loess", formula = y ~ x, color = ln_col, se = show_ci)
    } else {
      list(method = "lm", formula = y ~ x, color = ln_col, se = show_ci)
    }
    p <- ggplot2::ggplot(ad, ggplot2::aes(x = .data[[dose_col]], y = .data[[vol_col]])) +
      ggplot2::geom_point(alpha = 0.7, color = pt_col, size = pt_sz) +
      do.call(ggplot2::geom_smooth, smooth_args) +
      ggplot2::labs(title = "Tumor Volume vs Dose",
                   x = "Dose", y = "Tumor Volume") +
      ggplot2::theme_classic()
    p <- apply_plot_labels(p, input$dr_plot_title %||% "", input$dr_plot_subtitle %||% "",
                      input$dr_plot_title_align %||% "center", text_sz)
    current_dr_scatter(p)
    p
  }, width = function() input$dr_plot_width %||% 800,
     height = function() input$dr_plot_height %||% 600)

  output$dr_plot_tgi <- renderPlot({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(NULL)
    ad <- rv$dr_results$analysis_data
    req(ad)
    dose_col <- rv$dose_col %||% "Dose"
    vol_col  <- rv$vol_col  %||% "Volume"
    tx_col   <- rv$tx_col   %||% "Treatment"
    ctrl     <- input$dr_control_group %||% ""

    # Compute mean control volume
    ctrl_mean <- mean(ad[[vol_col]][ad[[tx_col]] == ctrl], na.rm = TRUE)
    if (is.na(ctrl_mean) || ctrl_mean == 0) return(NULL)

    # Compute per-observation %TGI
    tgi_data <- ad
    tgi_data$pct_tgi <- (1 - tgi_data[[vol_col]] / ctrl_mean) * 100

    pt_col   <- input$dr_scatter_pt_color   %||% "#555555"
    ln_col   <- input$dr_scatter_line_color %||% "#2166AC"
    pt_sz    <- input$dr_scatter_pt_size     %||% 2
    show_ci  <- isTRUE(input$dr_scatter_show_ci)
    trend    <- input$dr_trend_method %||% "linear"
    text_sz  <- input$dr_opt_text_size %||% 16
    smooth_args <- if (trend == "loess") {
      list(method = "loess", formula = y ~ x, color = ln_col, se = show_ci)
    } else {
      list(method = "lm", formula = y ~ x, color = ln_col, se = show_ci)
    }
    p <- ggplot2::ggplot(tgi_data, ggplot2::aes(x = .data[[dose_col]], y = .data[["pct_tgi"]])) +
      ggplot2::geom_point(alpha = 0.7, color = pt_col, size = pt_sz) +
      do.call(ggplot2::geom_smooth, smooth_args) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Tumor Growth Inhibition vs Dose",
                   x = "Dose", y = "% Tumor Growth Inhibition") +
      ggplot2::theme_classic()
    p <- apply_plot_labels(p, input$dr_plot_title %||% "", input$dr_plot_subtitle %||% "",
                      input$dr_plot_title_align %||% "center", text_sz)
    current_dr_tgi(p)
    p
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
    tbl <- as.data.frame(tukey[[1]])

    comp_type <- input$dr_comparison_type %||% "reference"
    if (comp_type == "reference") {
      ad  <- rv$dr_results$analysis_data
      if (!is.null(ad)) {
        dose_col <- rv$dose_col %||% "Dose"
        tx_col   <- rv$tx_col   %||% "Treatment"
        ctrl     <- input$dr_control_group %||% ""
        ref_doses <- unique(ad[ad[[tx_col]] == ctrl, dose_col])
        if (length(ref_doses) > 0) {
          ref_str <- as.character(ref_doses[[1]])
          keep <- grepl(paste0("^", ref_str, "-"), rownames(tbl)) |
                  grepl(paste0("-", ref_str, "$"), rownames(tbl))
          if (any(keep)) tbl <- tbl[keep, , drop = FALSE]
        }
      }
    }

    tbl$Comparison <- rownames(tbl)
    tbl <- tbl[, c("Comparison", setdiff(names(tbl), "Comparison"))]
    tbl <- format_df_for_display(tbl)
    datatable(tbl, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  # ---- Download handlers ----
  output$download_dr_table <- downloadHandler(
    filename = function() paste0("dose_response_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(file) {
      req(rv$dr_results, !isTRUE(rv$dr_results$error))
      write.csv(as.data.frame(rv$dr_results$summary_table), file, row.names = FALSE)
    }
  )

  output$download_dr_scatter <- downloadHandler(
    filename = function() paste0("dose_scatter_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      p <- current_dr_scatter()
      req(!is.null(p))
      png(file, width = input$dr_plot_width %||% 800,
                height = input$dr_plot_height %||% 600,
                units = "px", res = input$dr_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_dr_tgi <- downloadHandler(
    filename = function() paste0("dose_tgi_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content  = function(file) {
      p <- current_dr_tgi()
      req(!is.null(p))
      png(file, width = input$dr_plot_width %||% 800,
                height = input$dr_plot_height %||% 600,
                units = "px", res = input$dr_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_dr_stats <- downloadHandler(
    filename = function() paste0("dose_response_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
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
  })
}
