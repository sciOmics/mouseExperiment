#' Drug Synergy Analysis UI
#' @keywords internal
synergyUI <- function(id) {
  ns <- NS(id)
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
          
          uiOutput(ns("syn_group_selectors_ui")),
          
          hr(),
          
          h4("Analysis Options"),

          selectInput(ns("syn_model_type"), "Statistical Model",
            choices = c(
              "Frequentist (Bliss/Loewe)" = "frequentist",
              "Bayesian LMM (brms)"       = "bayes"
            ),
            selected = "frequentist"
          ),
          helpText("Frequentist computes Bliss/Loewe synergy scores. Bayesian fits a ",
                   "linear mixed-effects model with brms and reports posterior credible ",
                   "intervals; computation takes 3–12 minutes."),

          conditionalPanel(ns = ns, "input.syn_model_type == 'frequentist'",
            uiOutput(ns("syn_time_point_ui"))
          ),

          conditionalPanel(ns = ns, "input.syn_model_type == 'bayes'",
            hr(),
            h5(strong("Bayesian Settings")),
            selectInput(ns("syn_prior_strength"), "Prior Strength",
              choices = c(
                "Skeptical (recommended)" = "skeptical",
                "Weakly informative"      = "weakly_informative",
                "Informative"             = "informative",
                "Diffuse"                 = "diffuse"
              ),
              selected = "skeptical"
            ),
            selectInput(ns("syn_random_effects"), "Random Effects",
              choices = c(
                "Intercept only"    = "intercept_only",
                "Intercept + Slope" = "slope"
              ),
              selected = "intercept_only"
            ),
            numericInput(ns("syn_bayes_chains"), "MCMC chains",
                         value = 4, min = 1, max = 8, step = 1),
            numericInput(ns("syn_bayes_warmup"), "Warmup per chain",
                         value = 1000, min = 200, max = 5000, step = 100),
            numericInput(ns("syn_bayes_iter"), "Post-warmup draws per chain",
                         value = 500, min = 200, max = 5000, step = 100),
            selectInput(ns("syn_bayes_backend"), "Stan backend",
                        choices = c("rstan", "cmdstanr"),
                        selected = "rstan"),
            div(class = "alert alert-warning",
                style = "font-size:0.85em; padding:8px; margin-top:6px;",
              icon("clock"), " Bayesian fitting takes ~3–12 minutes. ",
              "A progress modal will block other tabs while the model runs, ",
              "but you can safely close the browser and return later — the ",
              "fit continues on the server and the results will be waiting ",
              "when you reconnect."
            )
          ),

          br(),

          actionButton(ns("run_syn_analysis"),
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),
          
          uiOutput(ns("syn_status_msg"))
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
          id = ns("syn_results_tabs"),

          tabPanel("Info",
            fluidRow(
              column(10, offset = 1,
                br(),

                # ── Overview ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Overview")),
                  div(class = "panel-body",
                    p("The Drug Synergy tab assesses whether a drug combination produces greater tumor suppression than predicted from the individual drugs acting independently. Two reference models are evaluated simultaneously:"),
                    tags$ul(
                      tags$li(strong("Bliss Independence"), " — treats the two drugs as probabilistically independent. The expected combined TGI is the sum of individual TGIs, capped at 100 %."),
                      tags$li(strong("Loewe Additivity"), " — treats the drugs as dose-equivalent. The Combination Index (CI) measures how the observed combination effect compares to the expected additive effect.")
                    ),
                    p("Four treatment groups must be assigned in the sidebar: Control, Drug A only, Drug B only, and the Combination (A + B). Analysis can be run at a single timepoint or tracked over the full study duration."),
                    p(tags$small(class = "text-muted",
                      "Tumor Growth Analysis must be run first — Drug Synergy uses preprocessed volume data from that tab."
                    ))
                  )
                ),

                # ── Bliss ────────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-info", "Bliss Independence"))),
                  div(class = "panel-body",
                    p("Under Bliss independence, drugs act on independent targets and their joint fractional effect (FE) is:"),
                    tags$blockquote(tags$code("Expected TGI = min(FE_A + FE_B, 1.0) × 100")),
                    p("where FE = TGI / 100 for each single agent. The cap at 1.0 prevents the expected value from exceeding the biological ceiling."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Bliss Difference"),
                      tags$dd("TGI_combo − Expected TGI. Positive = synergy; negative = antagonism; near zero = independence.")
                    ),
                    div(class = "alert alert-warning", style = "margin-top:10px;",
                      icon("triangle-exclamation"), " ",
                      strong("Ceiling effect: "),
                      "When either individual drug achieves TGI > 50 %, the Bliss expected TGI approaches or reaches 100 %, ",
                      "leaving no room to detect synergy even if it is real. ",
                      "In this regime Bliss Difference will be near zero or negative regardless of combination performance. ",
                      "Interpret Bliss results cautiously when single-agent TGIs are high."
                    )
                  )
                ),

                # ── Loewe ────────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "Loewe Additivity & Combination Index"))),
                  div(class = "panel-body",
                    p("Loewe additivity assumes that drugs acting on the same pathway are dose-equivalent and interchangeable. The Combination Index (CI) is:"),
                    tags$blockquote(tags$code("CI = min(FE_A + FE_B, 1.0) / FE_combo")),
                    p("Interpretation thresholds:"),
                    tags$ul(
                      tags$li(tags$b("CI < 0.85"), " — synergistic: the combination outperforms the additive expectation."),
                      tags$li(tags$b("0.85 ≤ CI ≤ 1.15"), " — additive: the combination matches the expected additive effect."),
                      tags$li(tags$b("CI > 1.15"), " — antagonistic: the combination underperforms the additive expectation.")
                    ),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Loewe Difference"),
                      tags$dd("TGI_combo − Loewe Expected TGI. Positive = synergy; negative = antagonism.")
                    ),
                    div(class = "alert alert-info", style = "margin-top:10px;",
                      icon("circle-info"), " ",
                      strong("Single-dose approximation: "),
                      "This implementation derives fractional effects from single-dose measurements rather than full dose-response curves for each drug. ",
                      "The CI is therefore an approximation of the true Loewe CI. ",
                      "A rigorous Loewe analysis requires IC50 values from complete dose-response curves — collect full dose-response data if a definitive CI is needed."
                    )
                  )
                ),

                # ── Over-time ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Over-Time Analysis")),
                  div(class = "panel-body",
                    p("When 'Analyze over time' is enabled in the sidebar, Bliss and Loewe scores are computed at every measured timepoint and plotted as a synergy trend. This reveals whether synergy emerges, fades, or reverses over the course of treatment."),
                    p("Per-timepoint p-values (Combination vs Drug A and vs Drug B, two-sample t-test) are included in the Synergy Metrics table. These are unadjusted for multiple testing across timepoints — interpret overall patterns rather than individual p-values at isolated days.")
                  )
                ),

                br()
              )
            )
          ),

          tabPanel("Summary",
            tab_info(strong("Synergy summary: "), "Overall assessment combining Bliss independence and Combination Index (CI) metrics. Reports whether the drug combination is synergistic, additive, or antagonistic at the evaluated timepoint, and lists the key quantitative scores."),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("syn_summary"))
            )
          ),
          
          tabPanel("Synergy Metrics",
            tab_info(
              strong("Quantitative synergy scores. "),
              "When \u2018Analyze over time\u2019 is disabled, the table has two columns (",
              strong("Metric"), " / ", strong("Value"), ") summarising scores at the selected timepoint: ",
              tags$ul(
                tags$li(strong("Bliss Difference"), " \u2014 observed combination effect minus the Bliss-independent expected effect (TGI\u2009combo \u2212 TGI\u2009Bliss). Positive = synergy; negative = antagonism."),
                tags$li(strong("Loewe Difference"), " \u2014 observed combination effect minus the Loewe additivity expected effect (sum of individual fractional effects, capped at 100%). Positive = synergy; negative = antagonism."),
                tags$li(strong("Combination Index (CI)"), " \u2014 Loewe-based CI: CI\u202f<\u202f0.85 synergistic; 0.85\u202f\u2264\u202fCI\u202f\u2264\u202f1.15 additive; CI\u202f>\u202f1.15 antagonistic."),
                tags$li(strong("Interpretation"), " \u2014 plain-language summary based on the CI threshold.")
              ),
              "When \u2018Analyze over time\u2019 is enabled, each row is a measured timepoint with columns:",
              tags$ul(
                tags$li(strong("Time_Point"), " \u2014 study day."),
                tags$li(strong("TGI_Drug_A / TGI_Drug_B / TGI_Combo"), " \u2014 % tumor growth inhibition for Drug A, Drug B, and the combination relative to control."),
                tags$li(strong("Bliss_Expected_TGI / Loewe_Expected_TGI"), " \u2014 predicted TGI under the Bliss independence and Loewe additivity models, respectively."),
                tags$li(strong("Bliss_Difference / Loewe_Difference"), " \u2014 TGI_Combo minus the expected TGI; positive values indicate synergy."),
                tags$li(strong("Combination_Index"), " \u2014 Loewe-based CI at this timepoint (see thresholds above)."),
                tags$li(strong("P_Value_vs_Drug_A / P_Value_vs_Drug_B"), " \u2014 two-sample t-test p-value comparing the combination group to each mono-therapy; values below 0.05 indicate a statistically significant difference."),
                tags$li(strong("Synergy_Assessment"), " \u2014 plain-language classification at each timepoint.")
              )
            ),
            div(class = "text-right", style = "margin-bottom: 6px;",
              downloadButton(ns("download_syn_metrics"), "Download Results")
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("syn_metrics_table"))
            )
          ),
          
          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  uiOutput(ns("syn_plot_selector_ui")),
                  ctrl_section("syn-appearance", "palette", "Appearance",
                    uiOutput(ns("syn_color_pickers_ui"))
                  ),
                  ctrl_section("syn-options", "sliders-h", "Options",
                    sliderInput(ns("syn_opt_text_size"), "Text size (pt)",
                                min = 10, max = 24, value = 16, step = 1),
                    selectInput(ns("syn_opt_legend_pos"), "Legend position",
                      choices = c("None" = "none",
                                  "Top" = "top", "Bottom" = "bottom",
                                  "Left" = "left", "Right" = "right",
                                  "In plot - top left"  = "inplot_tl",
                                  "In plot - top right" = "inplot_tr",
                                  "In plot - bottom left" = "inplot_bl",
                                  "In plot - bottom right" = "inplot_br"),
                      selected = "right"
                    ),
                    conditionalPanel(ns = ns,
                      "input.syn_plot_selector == 'trend'",
                      checkboxInput(ns("syn_opt_show_shading"),
                                    "Show Synergy / Antagonism shading",
                                    value = TRUE)
                    ),
                    conditionalPanel(ns = ns,
                      "input.syn_plot_selector == 'ci'",
                      checkboxInput(ns("syn_opt_ci_show_regions"),
                                    "Show region shading & labels",
                                    value = TRUE)
                    ),
                    conditionalPanel(ns = ns,
                      "input.syn_plot_selector == 'tgi'",
                      selectInput(ns("syn_opt_error_bar"), "Error bars",
                        choices  = c("SEM" = "SEM", "SD" = "SD",
                                     "95% CI" = "CI", "None" = "none"),
                        selected = "SEM"
                      )
                    )
                  ),
                  ctrl_section("syn-title", "heading", "Title & Subtitle",
                    textInput(ns("syn_plot_title"), "Title",
                              placeholder = "Leave blank for default"),
                    selectInput(ns("syn_plot_title_align"), "Title alignment",
                      choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                      selected = "center"),
                    textInput(ns("syn_plot_subtitle"), "Subtitle",
                              placeholder = "Optional subtitle")
                  ),
                  ctrl_section("syn-export", "arrows-alt", "Plot Size / Export",
                    fluidRow(
                      column(6, numericInput(ns("syn_plot_width"),  "Width (px)",  800, min = 200, max = 5000, step = 50)),
                      column(6, numericInput(ns("syn_plot_height"), "Height (px)", 600, min = 200, max = 5000, step = 50))
                    ),
                    numericInput(ns("syn_plot_dpi"), "Resolution (DPI)", 100, min = 72, max = 600, step = 10),
                    tags$div(style = "margin-top: 10px;", uiOutput(ns("syn_download_plot_btn")))
                  )
                )
              ),
              column(9,
                conditionalPanel(ns = ns,"input.syn_plot_selector == 'main'",
                  h4("Synergy Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("syn_plot_main"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.syn_plot_selector == 'trend'",
                  h4("Synergy Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("syn_plot_trend"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.syn_plot_selector == 'ci'",
                  h4("Combination Index"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("syn_plot_ci"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.syn_plot_selector == 'tgi'",
                  h4("Tumor Growth Inhibition"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("syn_plot_tgi"), width = "auto", height = "auto")
                  )
                )
              )
            )
          ),

          # CODE_REVIEW.md DIAGNOSTICS gap (6) — frequentist synergy
          # had no checkable t-test assumptions. New per-group Q-Q + boxplot.
          tabPanel("Diagnostics",
            tab_info(
              strong("Model-fit diagnostics. "),
              "Frequentist (t-test path): per-group normality + outlier ",
              "checks. Bayesian: MCMC convergence + posterior diagnostics. ",
              "The view switches based on the selected model type."
            ),
            # ── Freq diagnostics ───────────────────────────────────────────
            conditionalPanel(
              ns = ns,
              condition = "input.syn_model_type != 'bayes'",
              fluidRow(
                column(6,
                  h4("Per-group Q-Q of tumor volumes"),
                  tags$details(
                    tags$summary(tags$small("How to read this plot")),
                    tags$small(
                      "One panel per treatment group; points on the red dashed ",
                      "line = normally distributed. ",
                      strong("S-shape: "), "heavy tails — Welch is robust but ",
                      "consider a Mann-Whitney sanity check. ",
                      strong("Far outliers: "), "investigate before reporting."
                    )
                  ),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("syn_diag_group_qq_plot"), height = "350px")
                  )
                ),
                column(6,
                  h4("Per-group tumor volume boxplot"),
                  tags$details(
                    tags$summary(tags$small("How to read this plot")),
                    tags$small(
                      "Boxplots with individual mice jittered. Outliers in red. ",
                      "Use to confirm group means are dominated by the bulk of ",
                      "mice (not 1-2 extreme cases). ",
                      strong("Wide IQR + an outlier: "), "the synergy / ",
                      "antagonism call may be driven by that one mouse."
                    )
                  ),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("syn_diag_group_boxplot"), height = "350px")
                  )
                )
              )
            ),
            # ── Bayes MCMC diagnostics ─────────────────────────────────────
            conditionalPanel(
              ns = ns,
              condition = "input.syn_model_type == 'bayes'",
              tab_info(
                strong("MCMC convergence diagnostics. "),
                "All Rhat values should be ≤ 1.01 before interpreting results."
              ),
              uiOutput(ns("syn_bayes_diagnostics")),
              shinycssloaders::withSpinner(DTOutput(ns("syn_mcmc_diag_table"))),
              br(),
              h5(strong("Influential observations (Pareto-k > 0.7)")),
              em(class = "text-muted",
                 "Observations the LOO approximation can't reliably estimate."),
              shinycssloaders::withSpinner(DTOutput(ns("syn_bayes_influential"))),
              br(),
              shinycssloaders::withSpinner(
                plotOutput(ns("syn_mcmc_trace_plot"), height = "500px"))
            )
          )
        )
      )
    )
  )
}

#' Drug Synergy Analysis Server Logic
#' @keywords internal
synergy_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
  
  rv$syn_results <- NULL
  rv$syn_over_time <- FALSE

  observe({
    toggle_downloads(
      session,
      c("download_syn_metrics"),
      enabled = !is.null(rv$syn_results) && !isTRUE(rv$syn_results$error)
    )
  })

  # ReactiveVals to hold the current rendered plot (for download fidelity — D4, D5)
  current_syn_main_plot  <- reactiveVal(NULL)
  current_syn_trend_plot <- reactiveVal(NULL)
  current_syn_ci_plot    <- reactiveVal(NULL)
  current_syn_tgi_plot   <- reactiveVal(NULL)

  # ── Bayesian synergy tab management ───────────────────────────────────────
  bayes_syn_tabs_inserted <- reactiveVal(FALSE)

  observeEvent(input$syn_model_type, {
    tabs_id <- "syn_results_tabs"
    if (isTRUE(input$syn_model_type == "bayes") && !isTRUE(bayes_syn_tabs_inserted())) {
      appendTab(inputId = tabs_id,
        tabPanel("Posterior Check",
          tab_info(
            strong("Posterior predictive check (PPC)"), " — model-simulated distributions overlaid on observed data.",
            em("Good fit: simulated densities envelope the observed density.")
          ),
          div(class = "text-right", style = "margin-bottom: 6px;",
            downloadButton(session$ns("download_syn_pp_check"), "Download Plot")
          ),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("syn_pp_check_plot"), height = "450px")
          )
        ),
        select = FALSE
      )
      bayes_syn_tabs_inserted(TRUE)
    } else if (!isTRUE(input$syn_model_type == "bayes") && isTRUE(bayes_syn_tabs_inserted())) {
      removeTab(inputId = tabs_id, target = "Posterior Check")
      bayes_syn_tabs_inserted(FALSE)
    }
  }, ignoreNULL = FALSE)

  output$syn_pp_check_plot <- renderPlot({
    req(rv$syn_results)
    req(isTRUE(rv$syn_results$model_type_used == "bayes_synergy"))
    req(!isTRUE(rv$syn_results$error), !is.null(rv$bayes_models$syn))
    brms::pp_check(rv$bayes_models$syn, type = "dens_overlay", ndraws = 50)
  })

  output$syn_mcmc_diag_table <- renderDT({
    req(rv$syn_results)
    req(isTRUE(rv$syn_results$model_type_used == "bayes_synergy"))
    req(!isTRUE(rv$syn_results$error))
    diag <- rv$syn_results$mcmc_diagnostics
    req(!is.null(diag))
    datatable(format_df_for_display(as.data.frame(diag)),
              options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })

  output$syn_bayes_diagnostics <- renderUI({
    req(rv$syn_results)
    req(isTRUE(rv$syn_results$model_type_used == "bayes_synergy"))
    req(!isTRUE(rv$syn_results$error))
    bayes_diagnostics_panel(rv$syn_results)
  })

  output$syn_bayes_influential <- renderDT({
    req(rv$syn_results)
    req(isTRUE(rv$syn_results$model_type_used == "bayes_synergy"))
    req(!isTRUE(rv$syn_results$error))
    inf_df <- bayes_influential_obs(rv$syn_results)
    if (nrow(inf_df) == 0L) {
      return(datatable(
        data.frame(Note = "No observations exceed the Pareto-k = 0.7 threshold."),
        options = list(dom = "t", paging = FALSE, searching = FALSE), rownames = FALSE
      ))
    }
    datatable(inf_df,
              options = list(scrollX = TRUE, paging = FALSE, searching = FALSE),
              rownames = FALSE)
  })

  output$syn_mcmc_trace_plot <- renderPlot({
    req(rv$syn_results)
    req(isTRUE(rv$syn_results$model_type_used == "bayes_synergy"))
    req(!isTRUE(rv$syn_results$error), !is.null(rv$bayes_models$syn))
    bayesplot::mcmc_trace(posterior::as_draws_array(rv$bayes_models$syn))
  })

  syn_colors <- reactive({
    list(
      A     = input$syn_color_A     %||% "#2166AC",
      B     = input$syn_color_B     %||% "#D6604D",
      combo = input$syn_color_combo %||% "#762A83",
      ctrl  = input$syn_color_ctrl  %||% "#999999"
    )
  })
  
  # Dynamically render group selector dropdowns from treatment column
  output$syn_group_selectors_ui <- renderUI({
    req(rv$uploaded_data, rv$data_valid)
    
    treatment_col <- rv$settings$col_treatment
    if (is.null(treatment_col) || treatment_col == "") {
      return(p(class = "text-warning",
               icon("exclamation-triangle"),
               "Select Treatment column in Data Upload first."))
    }
    
    groups <- sort(unique(as.character(na.omit(rv$uploaded_data[[treatment_col]]))))
    choices <- setNames(groups, groups)
    
    tagList(
      selectInput(session$ns("syn_drug_a"), "Drug A (single agent)",
                 choices = choices,
                 selected = if (length(groups) >= 2) groups[2] else groups[1]),
      
      selectInput(session$ns("syn_drug_b"), "Drug B (single agent)",
                 choices = choices,
                 selected = if (length(groups) >= 3) groups[3] else groups[1]),
      
      selectInput(session$ns("syn_combo"), "Combination Treatment",
                 choices = choices,
                 selected = if (length(groups) >= 4) groups[4] else groups[1]),
      
      selectInput(session$ns("syn_control"), "Reference Group",
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
      day_candidate <- if (!is.null(rv$settings$col_day) && nzchar(rv$settings$col_day)) rv$settings$col_day else "Day"
      if (day_candidate %in% colnames(rv$uploaded_data))
        days <- sort(unique(na.omit(suppressWarnings(
          as.numeric(rv$uploaded_data[[day_candidate]])))))
    }
    day_choices <- c("All dates" = "all")
    if (!is.null(days) && length(days) > 0)
      day_choices <- c(day_choices, setNames(as.character(days), paste("Day", days)))
    selectInput(session$ns("syn_time_point"),
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
        status_alert("Please select all four treatment groups.", "danger")
      })
      return()
    }
    
    if (length(unique(needed)) < 4) {
      output$syn_status_msg <- renderUI({
        status_alert("Each role (Drug A, Drug B, Combination, Control) must be a different group.", "danger")
      })
      return()
    }
    
    output$syn_status_msg <- renderUI({
      status_alert("Running synergy analysis...", "info")
    })
    
    # Track analysis start with Umami
    shinyjs::runjs("window.trackAnalysisStart('drug_synergy', 'bliss_loewe', 'combo_vs_monotherapy');")
    
    # Use preprocessed data (Volume + Day resolved) if TG analysis has already
    # been run; otherwise preprocess inline from the raw uploaded data.
    pp <- tryCatch(ensure_preprocessed(rv), error = function(e) {
      log_caught_error("synergy", "preprocessing", e)
      output$syn_status_msg <- renderUI({
        status_alert(paste0("Preprocessing failed: ", conditionMessage(e)), "danger")
      })
      rv$syn_running <- FALSE
      return(NULL)
    })
    if (is.null(pp)) return()
    if (!is.null(pp$date_fmt_warning)) {
      showNotification(pp$date_fmt_warning, type = "warning", duration = 8)
    }
    analysis_data <- pp$data
    vol_col  <- pp$vol_col
    time_col <- pp$time_col
    tx_col   <- pp$tx_col
    
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

    # ── Bayesian async path ────────────────────────────────────────────────────
    if (isTRUE(input$syn_model_type == "bayes")) {

      syn_bayes_params <- list(
        data       = analysis_data,
        vol_col    = vol_col,
        time_col   = time_col,
        tx_col     = tx_col,
        drug_a     = input$syn_drug_a,
        drug_b     = input$syn_drug_b,
        combo      = input$syn_combo,
        control    = input$syn_control,
        endpoint_day = eval_tp,
        id_col     = rv$id_col,
        cage_col   = rv$cage_col,
        prior_str  = input$syn_prior_strength  %||% "skeptical",
        re_spec    = input$syn_random_effects   %||% "intercept_only",
        n_chains   = as.integer(input$syn_bayes_chains  %||% 4L),
        n_warmup   = as.integer(input$syn_bayes_warmup  %||% 1000L),
        n_iter     = as.integer(input$syn_bayes_iter    %||% 500L),
        backend    = input$syn_bayes_backend %||% "rstan"
      )

      # F.7 — route through the job registry.
      launch_syn_bayes_job(rv$session_id, syn_bayes_params,
                           input, output, session, rv)
      return()
    }

    # ── Frequentist path ───────────────────────────────────────────────────────
    show_analysis_modal("Fitting synergy model…", "Results will appear automatically when complete.")
    on.exit(shiny::removeModal(), add = TRUE)

    if (rv$syn_over_time) {
      result <- safe_analysis(module = "synergy_over_time", {
        mouseExperiment::analyze_drug_synergy_over_time(
          df               = analysis_data,
          treatment_column = tx_col,
          volume_column    = vol_col,
          time_column      = time_col,
          drug_a_name      = input$syn_drug_a,
          drug_b_name      = input$syn_drug_b,
          combo_name       = input$syn_combo,
          control_name     = input$syn_control,
          verbose          = FALSE
        )
      })
    } else {
      result <- safe_analysis(module = "synergy", {
        mouseExperiment::analyze_drug_synergy(
          df               = analysis_data,
          treatment_column = tx_col,
          volume_column    = vol_col,
          time_column      = time_col,
          drug_a_name      = input$syn_drug_a,
          drug_b_name      = input$syn_drug_b,
          combo_name       = input$syn_combo,
          control_name     = input$syn_control,
          eval_time_point  = eval_tp,
          verbose          = FALSE
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
            ggplot2::theme_classic()

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

        }, error = function(e) {
          log_caught_error("synergy", "trend plot generation", e)
        })
      } else {
        tryCatch({
          result$plot_main  <- mouseExperiment::plot_drug_synergy(result)
        }, error = function(e) {
          log_caught_error("synergy", "main plot generation", e)
        })
      }
    }
    
    rv$syn_results <- result
    updateTabsetPanel(session, "syn_results_tabs", selected = "Summary")

    # Track analysis completion with Umami
    if (isTRUE(result$error)) {
      shinyjs::runjs(sprintf("window.trackAnalysisComplete('drug_synergy', false, %s);",
                             js_string(result$message)))
    } else {
      shinyjs::runjs("window.trackAnalysisComplete('drug_synergy', true);")
    }

    # Store the 4 group names so the Plots tab can build the color-named vector
    rv$syn_plot_inputs <- list(
      drug_a  = input$syn_drug_a,
      drug_b  = input$syn_drug_b,
      combo   = input$syn_combo,
      control = input$syn_control
    )
    
    if (isTRUE(result$error)) {
      output$syn_status_msg <- renderUI({
        status_alert(paste("Error:", result$message), "danger")
      })
    } else {
      output$syn_status_msg <- renderUI({
        status_alert("Synergy analysis completed successfully!", "success")
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
      colourpicker::colourInput(session$ns("syn_color_A"),     label_a,     value = "#2166AC"),
      colourpicker::colourInput(session$ns("syn_color_B"),     label_b,     value = "#D6604D"),
      colourpicker::colourInput(session$ns("syn_color_combo"), label_combo, value = "#762A83"),
      colourpicker::colourInput(session$ns("syn_color_ctrl"),  label_ctrl,  value = "#999999")
    )
  })

  # Dynamic plot selector (Trend / CI only available for all-dates analysis)
  output$syn_plot_selector_ui <- renderUI({
    if (isTRUE(rv$syn_over_time)) {
      choices <- c("Synergy Plot" = "trend", "Combination Index" = "ci",
                   "Tumor Growth Inhibition" = "tgi")
    } else {
      choices <- c("Synergy Plot" = "main")
    }
    selectInput(session$ns("syn_plot_selector"), "Select Plot", choices = choices)
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
      elapsed_str <- format_elapsed(rv$syn_results$elapsed_seconds)
      if (!is.null(elapsed_str))
        cat("\nCOMPUTATION TIME:", elapsed_str, "\n")
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

  output$syn_download_plot_btn <- renderUI({
    btn_id <- switch(input$syn_plot_selector %||% "main",
      main   = "download_syn_plot_main",
      trend  = "download_syn_plot_trend",
      ci     = "download_syn_plot_ci",
      tgi    = "download_syn_plot_tgi",
      "download_syn_plot_main"
    )
    downloadButton(session$ns(btn_id), "Download Plot")
  })

  output$syn_plot_main <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)

    # Build colour vector matching the 4 synergy roles
    cols      <- syn_colors()
    col_A     <- cols$A
    col_B     <- cols$B
    col_combo <- cols$combo
    col_ctrl  <- cols$ctrl
    text_sz   <- input$syn_opt_text_size %||% 16

    final_p <- if (isTRUE(rv$syn_over_time)) {
      # Over-time Bliss plot stored in rv uses "Drug A" / "Drug B" etc. as colour labels
      p <- rv$syn_results$plot_main
      req(p)
      p + ggplot2::scale_color_manual(
        values = c("Drug A" = col_A, "Drug B" = col_B,
                   "Combination" = col_combo, "Bliss Value" = col_ctrl))
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
      tryCatch(
        mouseExperiment::plot_drug_synergy(rv$syn_results, custom_colors = custom_colors),
        error = function(e) {
          log_caught_error("synergy", "plot_drug_synergy re-render", e)
          rv$syn_results$plot_main
        }
      )
    }
    req(!is.null(final_p))
    final_p <- apply_plot_labels(final_p,
                                 input$syn_plot_title       %||% "",
                                 input$syn_plot_subtitle    %||% "",
                                 input$syn_plot_title_align %||% "center",
                                 text_sz)
    current_syn_main_plot(final_p)
    print(final_p)
  }, width = function() input$syn_plot_width %||% 800,
     height = function() input$syn_plot_height %||% 600)

  output$syn_plot_trend <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || !isTRUE(rv$syn_over_time)) return(NULL)
    req(rv$syn_results$synergy_summary)

    syn_df      <- rv$syn_results$synergy_summary
    grps        <- rv$syn_plot_inputs %||%
                   list(drug_a = "Drug A", drug_b = "Drug B",
                        combo = "Combination", control = "Control")
    cols        <- syn_colors()
    col_A       <- cols$A
    col_B       <- cols$B
    col_combo   <- cols$combo
    col_ctrl    <- cols$ctrl
    text_sz     <- input$syn_opt_text_size  %||% 16
    leg_pos     <- input$syn_opt_legend_pos %||% "right"
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
    p <- apply_plot_labels(p,
                           input$syn_plot_title       %||% "",
                           input$syn_plot_subtitle    %||% "",
                           input$syn_plot_title_align %||% "center",
                           text_sz)
    current_syn_trend_plot(p)
    print(p)
  }, width = function() input$syn_plot_width %||% 800,
     height = function() input$syn_plot_height %||% 600)
  
  output$syn_plot_ci <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || !isTRUE(rv$syn_over_time)) return(NULL)
    req(rv$syn_results$synergy_summary)

    syn_df      <- rv$syn_results$synergy_summary
    text_sz     <- input$syn_opt_text_size  %||% 16
    leg_pos     <- input$syn_opt_legend_pos %||% "right"
    show_regions <- isTRUE(input$syn_opt_ci_show_regions)

    mid_x <- mean(syn_df$Time_Point)
    p <- ggplot2::ggplot(syn_df,
                         ggplot2::aes(x = Time_Point, y = Combination_Index))

    if (show_regions) {
      p <- p +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0,    ymax = 0.85), fill = "blue",  alpha = 0.2) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0.85, ymax = 1.15), fill = "grey",  alpha = 0.4) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 1.15, ymax = 2),    fill = "red",   alpha = 0.2) +
        ggplot2::annotate("text", x = mid_x, y = 0.20, label = "Synergistic",
                          color = "blue", size = 5) +
        ggplot2::annotate("text", x = mid_x, y = 1.05, label = "Additive",    size = 5) +
        ggplot2::annotate("text", x = mid_x, y = 1.80, label = "Antagonistic",
                          color = "red",  size = 5)
    }

    p <- p +
      ggplot2::geom_line(color = "purple", linewidth = 1.5) +
      ggplot2::labs(x = "Days", y = "Combination Index", title = "") +
      ggplot2::ylim(0, 2) +
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = text_sz))

    p <- apply_legend_pos(p, leg_pos)
    p <- apply_plot_labels(p,
                           input$syn_plot_title       %||% "",
                           input$syn_plot_subtitle    %||% "",
                           input$syn_plot_title_align %||% "center",
                           text_sz)
    current_syn_ci_plot(p)
    print(p)
  }, width = function() input$syn_plot_width %||% 800,
     height = function() input$syn_plot_height %||% 600)
  
  # ------------------------------------------------------------------
  # TGI over time plot: per-mouse TGI with group lines + error bars
  # ------------------------------------------------------------------
  output$syn_plot_tgi <- renderPlot({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || !isTRUE(rv$syn_over_time)) return(NULL)
    req(rv$syn_results$synergy_summary)

    grps <- rv$syn_plot_inputs %||%
            list(drug_a = "Drug A", drug_b = "Drug B",
                 combo = "Combination", control = "Control")

    # Retrieve the preprocessed data and column names
    pp <- tryCatch(ensure_preprocessed(rv), error = function(e) NULL)
    if (!is.null(pp) && !is.null(pp$date_fmt_warning)) {
      showNotification(pp$date_fmt_warning, type = "warning", duration = 8)
    }
    req(pp)
    analysis_data <- pp$data
    vol_col  <- pp$vol_col
    time_col <- pp$time_col
    tx_col   <- pp$tx_col

    # Filter to the 4 synergy groups
    syn_groups <- c(grps$control, grps$drug_a, grps$drug_b, grps$combo)
    df <- analysis_data[analysis_data[[tx_col]] %in% syn_groups, ]
    req(nrow(df) > 0)

    # Coerce types
    df[[time_col]] <- suppressWarnings(as.numeric(df[[time_col]]))
    df[[vol_col]]  <- suppressWarnings(as.numeric(df[[vol_col]]))
    df <- df[!is.na(df[[time_col]]) & !is.na(df[[vol_col]]), ]

    # Compute per-mouse TGI at each timepoint
    # TGI = 100 * (1 - mouse_vol / control_mean_vol)
    days <- sort(unique(df[[time_col]]))
    tgi_rows <- list()
    for (d in days) {
      day_df <- df[df[[time_col]] == d, ]
      ctrl_vols <- day_df[[vol_col]][day_df[[tx_col]] == grps$control]
      ctrl_mean <- mean(ctrl_vols, na.rm = TRUE)
      if (is.na(ctrl_mean) || ctrl_mean <= 0) next
      # Compute TGI for each non-control mouse
      for (g in c(grps$drug_a, grps$drug_b, grps$combo)) {
        g_rows <- day_df[day_df[[tx_col]] == g, ]
        if (nrow(g_rows) == 0) next
        tgi_vals <- 100 * (1 - g_rows[[vol_col]] / ctrl_mean)
        tgi_rows[[length(tgi_rows) + 1L]] <- data.frame(
          Day       = d,
          Treatment = g,
          TGI       = tgi_vals,
          stringsAsFactors = FALSE
        )
      }
    }
    req(length(tgi_rows) > 0)
    tgi_df <- do.call(rbind, tgi_rows)
    tgi_df$Treatment <- factor(tgi_df$Treatment,
                               levels = c(grps$drug_a, grps$drug_b, grps$combo))

    # Summarise per group per day
    err_type <- input$syn_opt_error_bar %||% "SEM"
    summ <- do.call(rbind, lapply(split(tgi_df, list(tgi_df$Day, tgi_df$Treatment)), function(chunk) {
      if (nrow(chunk) == 0) return(NULL)
      x <- chunk$TGI[!is.na(chunk$TGI)]
      n <- length(x)
      if (n == 0) return(NULL)
      m <- mean(x)
      half <- switch(err_type,
        SEM = sd(x) / sqrt(n), SD = sd(x),
        CI = qt(0.975, df = max(n - 1L, 1L)) * sd(x) / sqrt(n), 0)
      data.frame(Day = chunk$Day[1], Treatment = chunk$Treatment[1],
                 mean_tgi = m, lo = m - half, hi = m + half,
                 stringsAsFactors = FALSE)
    }))
    summ$Treatment <- factor(summ$Treatment,
                             levels = c(grps$drug_a, grps$drug_b, grps$combo))

    # Colours and linetypes
    cols      <- syn_colors()
    col_A     <- cols$A
    col_B     <- cols$B
    col_combo <- cols$combo
    text_sz   <- input$syn_opt_text_size  %||% 16
    leg_pos   <- input$syn_opt_legend_pos %||% "right"

    color_vals <- stats::setNames(c(col_A, col_B, col_combo),
                                  c(grps$drug_a, grps$drug_b, grps$combo))
    lty_vals   <- stats::setNames(c("solid", "solid", "dashed"),
                                  c(grps$drug_a, grps$drug_b, grps$combo))

    p <- ggplot2::ggplot(summ, ggplot2::aes(x = Day, y = mean_tgi,
                                            color = Treatment, linetype = Treatment,
                                            group = Treatment)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2.5)

    if (err_type != "none") {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lo, ymax = hi),
        width = diff(range(summ$Day, na.rm = TRUE)) * 0.02,
        linewidth = 0.6)
    }

    p <- p +
      ggplot2::scale_color_manual(values = color_vals) +
      ggplot2::scale_linetype_manual(values = lty_vals) +
      ggplot2::labs(x = "Days", y = "Tumor Growth Inhibition (%)",
                    color = NULL, linetype = NULL) +
      ggplot2::theme_classic() +
      ggplot2::theme(text = ggplot2::element_text(size = text_sz))

    p <- apply_legend_pos(p, leg_pos)
    p <- apply_plot_labels(p,
                           input$syn_plot_title       %||% "",
                           input$syn_plot_subtitle    %||% "",
                           input$syn_plot_title_align %||% "center",
                           text_sz)
    current_syn_tgi_plot(p)
    print(p)
  }, width = function() input$syn_plot_width %||% 800,
     height = function() input$syn_plot_height %||% 600)
  
  # All 6 download handlers live in `mod_synergy_downloads.R`
  # (CODE_REVIEW.md C.1/C.3). The four plot PNG handlers each pull from a
  # `current_syn_*_plot` reactiveVal that the renderPlot bindings update;
  # those reactiveVals are passed in explicitly.
  add_syn_download_handlers(input, output, session, rv,
                            current_syn_main_plot, current_syn_trend_plot,
                            current_syn_ci_plot,   current_syn_tgi_plot)

  # CODE_REVIEW.md DIAGNOSTICS gap (6) — frequentist synergy diagnostics.
  output$syn_diag_group_qq_plot <- renderPlot({
    req(rv$syn_results)
    p <- rv$syn_results$diag_group_qq_plot
    if (is.null(p)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Per-group Q-Q not available (Bayesian or over-time path)")
      return(invisible(NULL))
    }
    p
  })

  output$syn_diag_group_boxplot <- renderPlot({
    req(rv$syn_results)
    p <- rv$syn_results$diag_group_boxplot
    if (is.null(p)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Per-group boxplot not available (Bayesian or over-time path)")
      return(invisible(NULL))
    }
    p
  })

  # F.7 — register the resume handler for bayes_syn jobs.
  register_syn_resume(input, output, session, rv)
  })
}
