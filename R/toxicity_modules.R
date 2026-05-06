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
              uiOutput(ns("tox_baseline_day_ui")),

              hr(),
              h5(strong("Post-hoc Comparisons")),
              selectInput(ns("tox_comparison_type"), NULL,
                choices = c(
                  "vs Reference Group"          = "reference",
                  "All Pairwise (BH-corrected)" = "pairwise"
                ),
                selected = "reference"
              ),
              helpText("'vs Reference' shows comparisons between each treatment and ",
                       "the reference/control group only. ",
                       "'All Pairwise' shows every group pair. Both apply ",
                       "Benjamini-Hochberg p-value correction.")
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

          tabPanel("Info",
            fluidRow(
              column(10, offset = 1,
                br(),

                # ── Overview ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Overview")),
                  div(class = "panel-body",
                    p("The Toxicity tab provides two complementary analysis modes selectable from the ",
                      strong("Statistical Model"), " dropdown:"),
                    tags$ul(
                      tags$li(strong("Body Weight Analysis"), " — models longitudinal weight trajectories using a linear mixed-effects model, and derives AUC-based and time-to-event (KM) toxicity summaries."),
                      tags$li(strong("Efficacy-Toxicity Analysis"), " — integrates body weight toxicity with tumor efficacy data to assess the overall therapeutic window. Requires a tumor volume column.")
                    ),
                    p(tags$small(class = "text-muted",
                      "Body weight is used as the toxicity surrogate throughout. ",
                      "When 'Subtract estimated tumor weight' is enabled, the gross weight is corrected by ",
                      "subtracting Volume(mm\u00b3) \u00d7 density(g/cm\u00b3) / 1000, giving a net host weight."))
                  )
                ),

                # ── Body Weight Analysis ────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-info", "Body Weight Analysis"), " \u2014 Analyses & Metrics")),
                  div(class = "panel-body",

                    h5(strong("Mixed-Effects Model")),
                    p("Fits a linear mixed-effects model to all repeated weight measurements:"),
                    tags$blockquote(
                      tags$code("Weight ~ Treatment \u00d7 Day + [covariates] + (Day | Mouse)")
                    ),
                    p("Random slopes allow each mouse's weight trend to deviate from the group average. ",
                      "Fixed effects test whether trajectories differ by treatment group over time."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("REML"), tags$dd("Restricted Maximum Likelihood — preferred for final parameter estimates (default)."),
                      tags$dt("ML"), tags$dd("Maximum Likelihood — use when comparing models with different fixed-effect structures.")
                    ),
                    tags$p(tags$b("Covariates (optional):")),
                    tags$ul(
                      tags$li(strong("Tumor Volume"), " — controls for tumor-driven weight changes independent of drug toxicity."),
                      tags$li(strong("Sex"), " — adjusts for sex-specific baseline weight differences."),
                      tags$li(strong("Initial Mass"), " — adjusts for pre-treatment weight variation across animals.")
                    ),

                    hr(),
                    h5(strong("Body Weight AUC")),
                    p("For each mouse, body weight is expressed as % change from baseline at each timepoint, then a trapezoidal AUC is computed over the study. ",
                      "A negative AUC indicates sustained weight loss; a positive AUC indicates net weight gain."),
                    tags$ul(
                      tags$li("Individual AUC values are shown as a jitter plot (± boxplot)."),
                      tags$li("Pairwise group comparisons use t-tests with Benjamini-Hochberg correction.")
                    ),

                    hr(),
                    h5(strong("Nadir Analysis")),
                    p("The nadir (minimum weight) is identified per mouse. Reported metrics include nadir day, nadir weight, % loss at nadir, and days to recovery."),

                    hr(),
                    h5(strong("Weight Loss KM Curve")),
                    p("Treats each mouse as an event-time observation where the event is crossing the weight-loss threshold (default 20%). ",
                      "Mice that never reach the threshold are censored at their last observation. ",
                      "A Kaplan-Meier curve is fitted per group; group differences are tested with the log-rank test."),
                    p(tags$b("Weight Loss Threshold:"),
                      " The % body weight loss that constitutes a toxicity event. The default (20%) is a common IACUC endpoint; adjust to match study-specific welfare criteria.")
                  )
                ),

                # ── Efficacy-Toxicity Analysis ──────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "Efficacy-Toxicity Analysis"), " \u2014 Analyses & Metrics")),
                  div(class = "panel-body",

                    h5(strong("Efficacy Metrics")),
                    p("Three metrics quantify anti-tumor efficacy per mouse relative to the control group:"),
                    tags$dl(
                      tags$dt("Final TGI (%)"),
                      tags$dd(
                        "Tumor Growth Inhibition at the final timepoint: ",
                        tags$code("TGI = (1 \u2212 V_treated / V\u0305_control) \u00d7 100"),
                        ". Values > 0 indicate tumor suppression; > 100 indicates tumor regression."
                      ),
                      tags$dt("Tumor AUC"),
                      tags$dd(
                        "Trapezoidal AUC of tumor volume over time, normalised to the control AUC: ",
                        tags$code("(1 \u2212 AUC_treated / AUC_control) \u00d7 100"),
                        ". Captures the full trajectory rather than a single endpoint."
                      ),
                      tags$dt("Log Cell Kill (LCK)"),
                      tags$dd(
                        "Estimates tumor cell kill based on growth delay: ",
                        tags$code("LCK = (T \u2212 C) / (3.32 \u00d7 T\u2082)"),
                        " where T and C are the days for treated and control tumors to reach a target volume, and T\u2082 is the control doubling time. Requires adequate tumor growth in controls."
                      )
                    ),

                    hr(),
                    h5(strong("Therapeutic Window Metric (TWM)")),
                    p("A summary ratio combining efficacy and toxicity per group:"),
                    tags$blockquote(tags$code("TWM = TGI(%) / Max Body Weight Loss(%)")),
                    p("Higher TWM indicates a more favourable balance. A noise floor prevents division by near-zero weight changes. Displayed as a ranked bar chart."),

                    hr(),
                    h5(strong("Safety-Efficacy Scatter")),
                    p("A per-mouse bivariate plot where each point represents one animal. ",
                      "X-axis shows the chosen efficacy metric; Y-axis shows maximum % body weight loss. ",
                      "Points in the top-right quadrant (high efficacy, high toxicity) represent the least favourable profile."),
                    tags$ul(
                      tags$li(strong("Group ellipses"), " — 68% confidence ellipses (1 SD) per treatment group."),
                      tags$li(strong("Group centroids"), " — mean efficacy and mean max weight loss per group, shown as large diamonds.")
                    ),

                    hr(),
                    h5(strong("Total Benefit Score")),
                    p("An integrated scalar combining efficacy and toxicity AUCs over time:"),
                    tags$blockquote(tags$code("B = AUC(efficacy) \u2212 \u03bb \u00d7 AUC(toxicity)")),
                    p("A positive score indicates net benefit; negative indicates net harm. The trade-off parameter ",
                      tags$b("\u03bb"), " (lambda) controls how strongly toxicity is penalised:"),
                    tags$ul(
                      tags$li(tags$b("\u03bb = 1.0"), " (default) — equal weight on efficacy and toxicity."),
                      tags$li(tags$b("\u03bb > 1"), " — penalise toxicity more heavily; use for fragile populations or curative intent."),
                      tags$li(tags$b("\u03bb < 1"), " — tolerate more toxicity; use when efficacy improvement is the primary goal.")
                    ),

                    hr(),
                    h5(strong("Weight-Corrected TGI")),
                    p("Recalculates TGI using only mice that stayed ", em("within"), " the safety weight-loss threshold throughout the study. ",
                      "Mice exceeding the threshold are excluded from the corrected calculation."),
                    p("This prevents a drug that arrested tumor growth ", em("due to severe toxicity"), " from appearing falsely effective. ",
                      "Comparing the standard and corrected TGI reveals the extent of this confound."),
                    p(tags$b("Safety threshold:"),
                      " Mirrors the Weight Loss Threshold in Body Weight Analysis. Mice whose maximum weight loss exceeds this value are excluded from the corrected TGI."),

                    hr(),
                    h5(strong("Correlation Heatmap")),
                    p("Shows Pearson correlations between all available per-mouse metrics (Efficacy, Max Weight Loss, Efficacy AUC, Toxicity AUC). ",
                      "Strong positive correlation between efficacy and toxicity suggests the drug cannot separate the two effects at the doses tested.")
                  )
                ),

                # ── Shared Parameters ────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Shared Parameters")),
                  div(class = "panel-body",
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Reference Group"),
                      tags$dd("The vehicle/control group used as the denominator for TGI, LCK, and AUC ratios. Auto-detected from common names (DMSO, Vehicle, Control, PBS, Saline) or defaults to the first group alphabetically."),
                      tags$dt("Subtract Tumor Weight"),
                      tags$dd("Removes the estimated tumor mass (Volume \u00d7 density / 1000) from gross body weight before all calculations, so that body weight reflects host tissue rather than tumor burden."),
                      tags$dt("Tumor Density"),
                      tags$dd("Used only when tumor weight adjustment is enabled. Default 1.0 g/cm\u00b3 (water density). Typical soft-tissue tumor density is 0.9\u20131.1 g/cm\u00b3."),
                      tags$dt("Baseline Day"),
                      tags$dd("The study day used to compute % change from baseline. Defaults to the earliest available day.")
                    )
                  )
                ),

                br()
              )
            )
          ),

          tabPanel("Summary",
            tab_info(
              strong("Analysis summary. "),
              "Body Weight Analysis: mixed-effects model output, fixed effects, and marginal means. ",
              "Efficacy-Toxicity Analysis: therapeutic window, benefit scores, and corrected TGI."
            ),
            shinycssloaders::withSpinner(verbatimTextOutput(ns("tox_summary")))
          ),

          tabPanel("Mixed Model",
            uiOutput(ns("tox_bw_mixed_model_tab_info")),
            h4("Fixed Effects"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_fixed_effects"))),
            downloadButton(ns("download_tox_model"), "Download Fixed Effects CSV"),
            hr(),
            h4("Estimated Marginal Means by Treatment"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_emmeans"))),
            downloadButton(ns("download_tox_emmeans_csv"), "Download EMMs CSV")
          ),

          tabPanel("AUC",
            uiOutput(ns("tox_bw_auc_tab_info")),
            h4("AUC Summary by Group"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_auc_summary"))),
            downloadButton(ns("download_tox_auc_csv"), "Download AUC CSV"),
            hr(),
            h4("Nadir Weight Analysis"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_nadir_table"))),
            downloadButton(ns("download_tox_nadir_csv"), "Download Nadir CSV"),
            hr(),
            uiOutput(ns("tox_comparisons_title")),
            shinycssloaders::withSpinner(DTOutput(ns("tox_auc_comparisons"))),
            downloadButton(ns("download_tox_comparisons_csv"), "Download Comparisons CSV")
          ),

          tabPanel("Weight Loss KM",
            tab_info(
              strong("Time-to-threshold event data. "),
              "Each row is one mouse. Event = 1 when the weight-loss threshold was crossed; ",
              "Event = 0 (censored) when the mouse never reached the threshold. ",
              "This table feeds the Kaplan-Meier plot in the Plots tab."
            ),
            h4("Time-to-Threshold Event Data"),
            shinycssloaders::withSpinner(DTOutput(ns("tox_event_table"))),
            downloadButton(ns("download_tox_events"), "Download Events CSV")
          ),

          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  conditionalPanel(ns = ns, "output.tox_model_type_val == 'bw'",
                    selectInput(ns("tox_bw_plot_selector"), "Plot Type",
                      choices = c(
                        "Body Weight Trajectory"  = "trajectory",
                        "Body Weight AUC"         = "auc",
                        "Weight Loss KM Curve"    = "km",
                        "AUC Comparison Forest"   = "bw_forest"
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
                    "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'auc'",
                    hr(),
                    checkboxInput(ns("tox_auc_show_boxplot"), "Show box plot", value = TRUE),
                    sliderInput(ns("tox_auc_pt_size"), "Point size",
                                min = 0.5, max = 6, value = 2, step = 0.5)
                  ),
                  conditionalPanel(ns = ns, "output.tox_model_type_val == 'bw'",
                    ctrl_section("tox-groups", "sort", "Groups",
                      uiOutput(ns("tox_group_order_color_ui"))
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
                  "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'bw_forest'",
                  h4("AUC Comparison Forest Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_bw_forest_plot"), width = "auto", height = "auto")
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
          ),

          tabPanel("Diagnostics",
            conditionalPanel(ns = ns, "output.tox_model_type_val == 'bw'",
              tab_info(
                strong("Mixed-effects model diagnostics "), "for the body weight trajectory model. ",
                "Inspect residual normality (Q-Q plot) and homoscedasticity (residuals vs fitted). ",
                em("Diagnostics are only available in Body Weight Analysis mode.")
              ),
              fluidRow(
                column(6,
                  h4("Q-Q Plot of Residuals"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_diag_qq"), height = "400px")
                  )
                ),
                column(6,
                  h4("Residuals vs Fitted"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_diag_resid_fitted"), height = "400px")
                  )
                )
              )
            ),
            conditionalPanel(ns = ns, "output.tox_model_type_val == 'et'",
              tab_info(
                "Model diagnostics are only available in ",
                strong("Body Weight Analysis"), " mode. ",
                "Switch the Statistical Model in the sidebar to access them."
              )
            )
          ),

          tabPanel(
            title = tags$span(icon("file-zipper"), " Download All"),
            value = "__tox_download_all__",
            wellPanel(
              style = "text-align:center; padding:30px;",
              h4("Export All Results"),
              p("Download every available result table (CSV) as a single ZIP file."),
              downloadButton(ns("download_tox_all"), "Download All Results",
                             class = "btn-lg btn-primary",
                             icon = icon("download"))
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

  # ---- Dynamic tab insertion based on model type ----
  bw_tabs_inserted <- reactiveVal(TRUE)   # BW tabs present in initial UI
  et_tabs_inserted <- reactiveVal(FALSE)

  observeEvent(input$tox_model_type, {
    mt      <- input$tox_model_type %||% "bw"
    tabs_id <- "tox_results_tabs"
    ns      <- session$ns

    if (mt == "et" && !isTRUE(et_tabs_inserted())) {
      if (isTRUE(bw_tabs_inserted())) {
        removeTab(inputId = tabs_id, target = "Mixed Model")
        removeTab(inputId = tabs_id, target = "AUC")
        removeTab(inputId = tabs_id, target = "Weight Loss KM")
        bw_tabs_inserted(FALSE)
      }
      insertTab(inputId = tabs_id, target = "Summary", position = "after",
        tabPanel("Therapeutic Window",
          tab_info(
            strong("Therapeutic Window Metric (TWM). "),
            "TWM = TGI(%) / Max Body Weight Loss(%). Higher values indicate a more ",
            "favourable efficacy-to-toxicity ratio."
          ),
          shinycssloaders::withSpinner(DTOutput(ns("tox_et_twm_table"))),
          downloadButton(ns("download_et_twm_csv"), "Download TWM CSV")
        )
      )
      insertTab(inputId = tabs_id, target = "Therapeutic Window", position = "after",
        tabPanel("Safety-Efficacy",
          tab_info(
            strong("Per-mouse and group-level safety-efficacy data. "),
            "Each row is one mouse (or group). Efficacy is the chosen metric relative to control; ",
            "toxicity is maximum % body weight loss."
          ),
          h4("Per-Mouse Data"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_et_per_mouse_table"))),
          downloadButton(ns("download_et_bivariate_csv"), "Download Per-Mouse CSV"),
          hr(),
          h4("Group Summary"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_et_per_group_table")))
        )
      )
      insertTab(inputId = tabs_id, target = "Safety-Efficacy", position = "after",
        tabPanel("Benefit Score",
          tab_info(
            strong("Total Benefit Score. "),
            "B = AUC(efficacy) \u2212 \u03bb \u00d7 AUC(toxicity). Positive = net benefit; ",
            "negative = net harm. Adjust \u03bb in the sidebar to change the toxicity penalty."
          ),
          shinycssloaders::withSpinner(DTOutput(ns("tox_et_benefit_table"))),
          downloadButton(ns("download_et_benefit_csv"), "Download Benefit CSV")
        )
      )
      insertTab(inputId = tabs_id, target = "Benefit Score", position = "after",
        tabPanel("Corrected TGI",
          tab_info(
            strong("Weight-corrected vs standard TGI. "),
            "Mice exceeding the safety weight-loss threshold are excluded from the corrected TGI. ",
            "A large gap between the two values indicates that efficacy was partly driven by toxicity."
          ),
          h4("Weight-Corrected vs Standard TGI"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_et_wctgi_table"))),
          hr(),
          h4("Excluded Mice"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_et_excluded_table"))),
          downloadButton(ns("download_et_wctgi_csv"), "Download Corrected TGI CSV")
        )
      )
      insertTab(inputId = tabs_id, target = "Corrected TGI", position = "after",
        tabPanel("Combined Report",
          tab_info(
            strong("Combined results table. "),
            "Merges TWM, Benefit Score, Corrected TGI, and exclusion count per group into a ",
            "single summary for export."
          ),
          shinycssloaders::withSpinner(DTOutput(ns("tox_et_combined_table"))),
          downloadButton(ns("download_et_report"), "Download Combined Report CSV")
        )
      )
      et_tabs_inserted(TRUE)

    } else if (mt == "bw" && !isTRUE(bw_tabs_inserted())) {
      if (isTRUE(et_tabs_inserted())) {
        removeTab(inputId = tabs_id, target = "Therapeutic Window")
        removeTab(inputId = tabs_id, target = "Safety-Efficacy")
        removeTab(inputId = tabs_id, target = "Benefit Score")
        removeTab(inputId = tabs_id, target = "Corrected TGI")
        removeTab(inputId = tabs_id, target = "Combined Report")
        et_tabs_inserted(FALSE)
      }
      insertTab(inputId = tabs_id, target = "Summary", position = "after",
        tabPanel("Mixed Model",
          uiOutput(ns("tox_bw_mixed_model_tab_info")),
          h4("Fixed Effects"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_fixed_effects"))),
          downloadButton(ns("download_tox_model"), "Download Fixed Effects CSV"),
          hr(),
          h4("Estimated Marginal Means by Treatment"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_emmeans"))),
          downloadButton(ns("download_tox_emmeans_csv"), "Download EMMs CSV")
        )
      )
      insertTab(inputId = tabs_id, target = "Mixed Model", position = "after",
        tabPanel("AUC",
          uiOutput(ns("tox_bw_auc_tab_info")),
          h4("AUC Summary by Group"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_auc_summary"))),
          downloadButton(ns("download_tox_auc_csv"), "Download AUC CSV"),
          hr(),
          h4("Nadir Weight Analysis"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_nadir_table"))),
          downloadButton(ns("download_tox_nadir_csv"), "Download Nadir CSV"),
          hr(),
          uiOutput(ns("tox_comparisons_title")),
          shinycssloaders::withSpinner(DTOutput(ns("tox_auc_comparisons"))),
          downloadButton(ns("download_tox_comparisons_csv"), "Download Comparisons CSV")
        )
      )
      insertTab(inputId = tabs_id, target = "AUC", position = "after",
        tabPanel("Weight Loss KM",
          tab_info(
            strong("Time-to-threshold event data. "),
            "Each row is one mouse. Event = 1 when the weight-loss threshold was crossed; ",
            "Event = 0 (censored) when the mouse never reached the threshold. ",
            "This table feeds the Kaplan-Meier plot in the Plots tab."
          ),
          h4("Time-to-Threshold Event Data"),
          shinycssloaders::withSpinner(DTOutput(ns("tox_event_table"))),
          downloadButton(ns("download_tox_events"), "Download Events CSV")
        )
      )
      bw_tabs_inserted(TRUE)
    }
  }, ignoreNULL = FALSE)

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

    # Store reference group and treatment groups for colour/order pickers
    rv$tox_ref_group   <- ref_group
    rv$tox_plot_groups <- if (!is.null(results$trajectory_data) && nrow(results$trajectory_data) > 0)
      unique(as.character(results$trajectory_data$Treatment))
    else if (!is.null(tx_col) && has_col(data, tx_col))
      unique(as.character(data[[tx_col]]))
    else NULL

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

  # ---- BW: Tab info outputs ----

  output$tox_bw_mixed_model_tab_info <- renderUI({
    tab_info(
      strong("Body Weight Mixed Model. "),
      "Fixed effects test whether weight trajectories differ by treatment group over time. ",
      "The ", tags$code("Treatment \u00d7 Day"), " interaction is the primary term of interest. ",
      "Estimated marginal means give group-level weight estimates adjusted for covariates."
    )
  })

  output$tox_bw_auc_tab_info <- renderUI({
    ctype <- input$tox_comparison_type %||% "reference"
    if (ctype == "reference") {
      tab_info(
        strong("Body Weight AUC summaries and group comparisons. "),
        "AUC is the trapezoidal integral of % weight change from baseline over the study. ",
        "Comparisons shown are vs the reference/control group with BH p-value correction. ",
        tags$small(class = "text-muted", "Change comparison type in the sidebar.")
      )
    } else {
      tab_info(
        strong("Body Weight AUC summaries and all-pairwise comparisons. "),
        "AUC is the trapezoidal integral of % weight change from baseline over the study. ",
        "All group pairs are shown with Benjamini-Hochberg p-value correction. ",
        tags$small(class = "text-muted", "Change comparison type in the sidebar.")
      )
    }
  })

  output$tox_comparisons_title <- renderUI({
    ctype <- input$tox_comparison_type %||% "reference"
    lbl   <- if (ctype == "reference") "vs Reference Group" else "All Pairwise (BH-corrected)"
    tagList(
      h4(paste0("AUC Comparisons \u2014 ", lbl)),
      tags$p(class = "text-muted", tags$small("Change comparison type in the sidebar."))
    )
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
    auc      <- rv$tox_results$bw_auc
    cmp_type <- input$tox_comparison_type %||% "reference"
    if (isTRUE(auc$error) || is.null(auc$comparisons))
      return(datatable(data.frame(Message = "Not available")))
    df  <- auc$comparisons
    ref <- rv$tox_ref_group %||% ""
    if (cmp_type == "reference" && ref != "")
      df <- df[df$Group1 == ref | df$Group2 == ref, , drop = FALSE]
    if (nrow(df) > 0) df$P_Adjusted_BH <- stats::p.adjust(df$P_Value, method = "BH")
    datatable(format_df_for_display(df), options = list(scrollX = TRUE), rownames = FALSE)
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

  # ---- Diagnostics ----

  output$tox_diag_qq <- renderPlot({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$model)) return(NULL)
    r <- tryCatch(stats::residuals(mm$model), error = function(e) NULL)
    req(!is.null(r))
    stats::qqnorm(r, main = "Q-Q Plot of Mixed Model Residuals", pch = 16,
                  col = "#2c7bb6")
    stats::qqline(r, col = "#d7191c", lwd = 2)
  })

  output$tox_diag_resid_fitted <- renderPlot({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$model)) return(NULL)
    r <- tryCatch(stats::residuals(mm$model), error = function(e) NULL)
    f <- tryCatch(stats::fitted(mm$model),    error = function(e) NULL)
    req(!is.null(r), !is.null(f))
    plot(f, r, main = "Residuals vs Fitted",
         xlab = "Fitted Values", ylab = "Residuals",
         pch = 16, col = "#2c7bb6")
    graphics::abline(h = 0, col = "#d7191c", lty = 2, lwd = 2)
  })

  # ---- BW AUC Forest Plot ----

  output$tox_bw_forest_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc      <- rv$tox_results$bw_auc
    if (isTRUE(auc$error) || is.null(auc$comparisons)) return(NULL)
    cmp_type <- input$tox_comparison_type %||% "reference"
    text_sz  <- input$tox_plot_text_size %||% 14
    df  <- auc$comparisons
    ref <- rv$tox_ref_group %||% ""
    if (cmp_type == "reference" && ref != "")
      df <- df[df$Group1 == ref | df$Group2 == ref, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df$P_Adjusted_BH <- stats::p.adjust(df$P_Value, method = "BH")
    df$Comparison    <- paste(df$Group1, "vs", df$Group2)
    df$Comparison    <- factor(df$Comparison, levels = rev(df$Comparison))
    df$Significant   <- df$P_Adjusted_BH < 0.05
    lbl <- if (cmp_type == "reference") "vs Reference Group" else "All Pairwise"
    ggplot2::ggplot(df,
      ggplot2::aes(x = .data[["Mean_Diff"]], y = .data[["Comparison"]],
                   color = .data[["Significant"]])) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data[["CI_Lower"]], xmax = .data[["CI_Upper"]]),
        height = 0.25, linewidth = 0.8
      ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "#d7191c", "FALSE" = "#999999"),
        labels = c("TRUE" = "BH p\u202f<\u202f0.05", "FALSE" = "n.s."),
        name = ""
      ) +
      ggplot2::labs(
        title    = "BW AUC Comparison Forest Plot",
        subtitle = sprintf("%s comparisons; Benjamini-Hochberg correction", lbl),
        x = "Mean Difference in AUC (% weight change)", y = NULL
      ) +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "top")
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  # ---- BW Plots ----

  output$tox_group_order_color_ui <- renderUI({
    req(rv$tox_plot_groups)
    group_order_and_color_ui("tox", rv$tox_plot_groups, ns = session$ns)
  })

  observe_palette_preset(input, session, "tox", reactive(rv$tox_plot_groups))

  output$tox_trajectory_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$trajectory_data)
    td <- rv$tox_results$trajectory_data
    req(nrow(td) > 0)
    traj_type    <- input$tox_traj_type %||% "mean_sem"
    text_sz      <- input$tox_plot_text_size %||% 14
    line_w       <- input$tox_traj_line_width %||% 1
    ordered_grps <- get_ordered_groups(input, "tox", rv$tox_plot_groups)
    colors       <- build_color_vector(input, "tox_color", ordered_grps)

    if (traj_type == "spaghetti") {
      td$.__LineGroup <- paste(td[["Treatment"]], td[["ID"]], sep = "\u00b7")
      if (!is.null(ordered_grps))
        td$Treatment <- factor(td$Treatment, levels = ordered_grps)
      p <- ggplot2::ggplot(td, ggplot2::aes(.data[["Day"]], .data[["Weight"]],
                                            group = .data[[".__LineGroup"]], color = .data[["Treatment"]])) +
        ggplot2::geom_line(alpha = 0.6, linewidth = line_w) +
        ggplot2::geom_point(size = 1.5, alpha = 0.6) +
        ggplot2::labs(title = "Body Weight Trajectories",
                     x = "Day", y = "Weight (g)", color = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
      if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)
      p
    } else {
      agg <- stats::aggregate(Weight ~ Treatment + Day, data = td, FUN = function(x)
        c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = length(x)))
      af  <- do.call(data.frame, agg)
      names(af) <- c("Treatment", "Day", "Mean", "SD", "N")
      af$SEM   <- af$SD / sqrt(af$N)
      af$Error <- if (traj_type == "mean_sd") af$SD else af$SEM
      err_lbl  <- if (traj_type == "mean_sd") "SD" else "SEM"
      if (!is.null(ordered_grps))
        af$Treatment <- factor(af$Treatment, levels = ordered_grps)
      p <- ggplot2::ggplot(af, ggplot2::aes(.data[["Day"]], .data[["Mean"]],
                                            color = .data[["Treatment"]], fill = .data[["Treatment"]])) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["Mean"]] - .data[["Error"]],
                                          ymax = .data[["Mean"]] + .data[["Error"]]),
                             alpha = 0.2, color = NA) +
        ggplot2::geom_line(linewidth = line_w) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(title = paste0("Body Weight Over Time (Mean \u00b1 ", err_lbl, ")"),
                     x = "Day", y = "Weight (g)", color = "Treatment", fill = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
      if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors) +
                                     ggplot2::scale_fill_manual(values = colors)
      p
    }
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_auc_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error)) return(NULL)
    text_sz       <- input$tox_plot_text_size %||% 14
    show_box      <- isTRUE(input$tox_auc_show_boxplot %||% TRUE)
    pt_sz         <- input$tox_auc_pt_size %||% 2
    ordered_grps  <- get_ordered_groups(input, "tox", rv$tox_plot_groups)
    colors        <- build_color_vector(input, "tox_color", ordered_grps)
    auc_df        <- auc$auc_per_mouse
    if (!is.null(ordered_grps))
      auc_df$Treatment <- factor(auc_df$Treatment, levels = ordered_grps)
    p <- ggplot2::ggplot(auc_df,
      ggplot2::aes(.data[["Treatment"]], .data[["AUC_Pct_Change"]], fill = .data[["Treatment"]]))
    if (show_box)
      p <- p + ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6)
    p <- p +
      ggplot2::geom_jitter(ggplot2::aes(color = .data[["Treatment"]]),
                           width = 0.2, alpha = 0.8, size = pt_sz) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Body Weight AUC (% Change from Baseline)",
                   x = "Treatment", y = "AUC of % Weight Change") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "none")
    if (!is.null(colors)) p <- p +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_color_manual(values = colors)
    p
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_km_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(NULL)
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- get_ordered_groups(input, "tox", rv$tox_plot_groups)
    colors       <- build_color_vector(input, "tox_color", ordered_grps)
    # Apply group order to event_data so ggsurvplot respects it
    event_df <- wl$event_data
    if (!is.null(ordered_grps) && has_col(event_df, "Treatment"))
      event_df$Treatment <- factor(event_df$Treatment, levels = ordered_grps)
    km_fit <- tryCatch(
      survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = event_df),
      error = function(e) wl$km_fit
    )
    pal <- if (!is.null(colors)) unname(colors) else "hue"
    km_plot <- tryCatch(
      survminer::ggsurvplot(km_fit, data = event_df,
        palette    = pal,
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
           ggplot2::aes(.data[["Efficacy"]], .data[["Max_Pct_Weight_Loss"]],
                        color = .data[["Treatment"]])) +
      ggplot2::geom_point(size = pt_size, alpha = 0.7)
    if (show_ell && nrow(pm) > 3)
      p <- p + ggplot2::stat_ellipse(level = 0.68, linewidth = 0.8, show.legend = FALSE)
    if (show_cent && !is.null(biv$per_group))
      p <- p + ggplot2::geom_point(data = biv$per_group,
                 ggplot2::aes(.data[["Efficacy_Mean"]], .data[["Toxicity_Mean"]]),
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
    bt$Treatment <- factor(bt$Treatment, levels = bt$Treatment[order(bt$Benefit_Score)])
    ggplot2::ggplot(bt, ggplot2::aes(.data[["Treatment"]], .data[["Benefit_Score"]],
                                      fill = .data[["Benefit_Score"]] > 0)) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["Benefit_Score"]], 1)),
                         vjust = ifelse(bt$Benefit_Score >= 0, -0.5, 1.5), size = 4) +
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
      TGI       = c(comp$TGI_Uncorrected, comp$TGI_Corrected),
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
    md <- pm[, c("Efficacy", "Max_Pct_Weight_Loss"), drop = FALSE]
    colnames(md) <- c("Efficacy", "Max Weight Loss (%)")
    if (!is.null(rv$tox_results$benefit) && !isTRUE(rv$tox_results$benefit$error)) {
      ea <- rv$tox_results$benefit$efficacy_auc
      ta <- rv$tox_results$benefit$toxicity_auc
      if (!is.null(ea) && nrow(ea) == nrow(pm)) md[["Efficacy AUC"]] <- ea$Efficacy_AUC
      if (!is.null(ta) && nrow(ta) == nrow(pm)) md[["Toxicity AUC"]] <- ta$Toxicity_AUC
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
  output$download_tox_emmeans_csv <- downloadHandler("body_weight_emmeans.csv", function(file) {
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$emmeans_table)) return()
    utils::write.csv(mm$emmeans_table, file, row.names = FALSE)
  })
  output$download_tox_auc_csv <- downloadHandler("body_weight_auc.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bw_auc)
    if (isTRUE(rv$tox_results$bw_auc$error)) return()
    utils::write.csv(rv$tox_results$bw_auc$auc_per_mouse, file, row.names = FALSE)
  })
  output$download_tox_nadir_csv <- downloadHandler("body_weight_nadir.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bw_auc)
    if (isTRUE(rv$tox_results$bw_auc$error)) return()
    utils::write.csv(rv$tox_results$bw_auc$nadir_data, file, row.names = FALSE)
  })
  output$download_tox_comparisons_csv <- downloadHandler("bw_auc_comparisons.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc      <- rv$tox_results$bw_auc
    cmp_type <- isolate(input$tox_comparison_type) %||% "reference"
    if (isTRUE(auc$error) || is.null(auc$comparisons)) return()
    df  <- auc$comparisons
    ref <- rv$tox_ref_group %||% ""
    if (cmp_type == "reference" && ref != "")
      df <- df[df$Group1 == ref | df$Group2 == ref, , drop = FALSE]
    if (nrow(df) > 0) df$P_Adjusted_BH <- stats::p.adjust(df$P_Value, method = "BH")
    utils::write.csv(df, file, row.names = FALSE)
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

  # ---- Download All (ZIP) ----
  output$download_tox_all <- downloadHandler(
    filename = "toxicity_results.zip",
    content  = function(file) {
      req(rv$tox_results)
      res    <- rv$tox_results
      tmpdir <- tempfile(pattern = "tox_export_")
      dir.create(tmpdir, recursive = TRUE)
      on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

      write_if <- function(df, fname) {
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0)
          utils::write.csv(df, file.path(tmpdir, fname), row.names = FALSE)
      }

      if (res$model_type == "bw") {
        mm <- res$mixed_model
        if (!isTRUE(mm$error)) {
          write_if(mm$fixed_effects,  "bw_fixed_effects.csv")
          write_if(mm$emmeans_table,  "bw_emmeans.csv")
        }
        bwa <- res$bw_auc
        if (!isTRUE(bwa$error)) {
          write_if(bwa$auc_per_mouse, "bw_auc_per_mouse.csv")
          write_if(bwa$auc_summary,   "bw_auc_summary.csv")
          write_if(bwa$nadir_data,    "bw_nadir.csv")
          if (!is.null(bwa$comparisons)) {
            df  <- bwa$comparisons
            df$P_Adjusted_BH <- stats::p.adjust(df$P_Value, method = "BH")
            write_if(df, "bw_auc_comparisons.csv")
          }
        }
        wl <- res$weight_loss
        if (!isTRUE(wl$error)) write_if(wl$event_data, "bw_weight_loss_events.csv")
      } else {
        if (!isTRUE(res$twm$error))      write_if(res$twm$twm_table,          "et_twm.csv")
        if (!isTRUE(res$bivariate$error)) {
          write_if(res$bivariate$per_mouse, "et_per_mouse.csv")
          write_if(res$bivariate$per_group, "et_per_group.csv")
        }
        if (!isTRUE(res$benefit$error))  write_if(res$benefit$benefit_table,  "et_benefit.csv")
        if (!isTRUE(res$wctgi$error)) {
          write_if(res$wctgi$comparison,     "et_corrected_tgi.csv")
          write_if(res$wctgi$excluded_summary, "et_excluded_mice.csv")
        }
      }

      old_wd <- getwd()
      setwd(tmpdir)
      on.exit(setwd(old_wd), add = TRUE)
      csv_files <- list.files(".", pattern = "\\.csv$")
      if (length(csv_files) == 0) return()
      utils::zip(zipfile = file, files = csv_files)
    }
  )

  }) # end moduleServer
}
