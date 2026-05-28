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
                "Body Weight Analysis"            = "bw",
                "Efficacy-Toxicity Analysis"      = "et",
                "Bayesian Therapeutic Window (brms)" = "bayes_twm"
              ),
              selected = "bw"
            ),
            helpText(
              "Body Weight Analysis models weight trajectories with mixed-effects models, ",
              "AUC, and time-to-threshold survival. ",
              "Efficacy-Toxicity Analysis also requires a tumor volume column and computes ",
              "the therapeutic window, safety-efficacy scatter, and total benefit score. ",
              "Bayesian Therapeutic Window requires Bayesian TG and Bayesian BW analyses ",
              "to have been run first."
            ),

            conditionalPanel(ns = ns, "input.tox_model_type == 'bayes_twm'",
              uiOutput(ns("bayes_twm_readiness_ui"))
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
              selectInput(ns("tox_bw_method"), "Statistical Method",
                choices = c(
                  "Frequentist (LME4)" = "lme4",
                  "Bayesian (brms)"    = "bayes"
                ),
                selected = "lme4"
              ),
              helpText("LME4 fits a linear mixed-effects model. Bayesian reports full ",
                       "posterior distributions and credible intervals; computation ",
                       "takes 3–12 minutes."),
              conditionalPanel(ns = ns,
                "input.tox_model_type == 'bw' && input.tox_bw_method == 'bayes'",
                hr(),
                h5(strong("Bayesian Settings")),
                selectInput(ns("tox_bw_prior_strength"), "Prior Strength",
                  choices = c(
                    "Skeptical (recommended)" = "skeptical",
                    "Weakly informative"      = "weakly_informative",
                    "Informative"             = "informative",
                    "Diffuse"                 = "diffuse"
                  ),
                  selected = "skeptical"
                ),
                selectInput(ns("tox_bw_random_effects"), "Random Effects",
                  choices = c(
                    "Intercept only" = "intercept_only",
                    "Intercept + Slope" = "slope"
                  ),
                  selected = "intercept_only"
                ),
                numericInput(ns("tox_bw_bayes_chains"), "MCMC chains",
                             value = 4, min = 1, max = 8, step = 1),
                numericInput(ns("tox_bw_bayes_warmup"), "Warmup per chain",
                             value = 1000, min = 200, max = 5000, step = 100),
                numericInput(ns("tox_bw_bayes_iter"), "Post-warmup draws per chain",
                             value = 500, min = 200, max = 5000, step = 100),
                div(class = "alert alert-warning",
                    style = "font-size:0.85em; padding:8px; margin-top:6px;",
                  icon("clock"), " Bayesian fitting takes ~3–12 minutes. ",
                  "The app stays usable while the model runs in the background."
                )
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
                            "Tumor AUC" = "tumor_auc"),
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
                      tags$li("Individual AUC values are shown as a jitter plot with a configurable centre mark (mean or median) and error bars (SEM, SD, or 95% CI)."),
                      tags$li("Units are ", tags$b("% \u00b7 days"), ": e.g. a value of \u2212100 means an average 5% weight loss sustained over 20 days."),
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
                    ),

                    hr(),
                    h5(strong("Therapeutic Window Metric (TWM)")),
                    p("A summary ratio combining efficacy and toxicity per group:"),
                    tags$blockquote(tags$code("TWM = TGI(%) / Mean Body Weight Loss(%)")),
                    p("Higher TWM indicates a more favourable balance. The denominator is the group mean of per-mouse maximum weight loss (not the single worst-case mouse). A noise floor prevents division by near-zero weight changes. Displayed as a ranked bar chart."),

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

                    br()
                  )
                ),

                # ── Shared Parameters ────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Shared Parameters")),
                  div(class = "panel-body",
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Reference Group"),
                      tags$dd("The vehicle/control group used as the denominator for TGI and AUC ratios. Auto-detected from common names (DMSO, Vehicle, Control, PBS, Saline) or defaults to the first group alphabetically."),
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

          tabPanel("Model Results",
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
                        "Weight-Corrected TGI"    = "wctgi"
                      ),
                      selected = "twm"
                    )
                  ),
                  conditionalPanel(ns = ns,
                    "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'trajectory'",
                    hr(),
                    radioButtons(ns("tox_traj_type"), "Display style",
                      choices = c("Individual mice"            = "spaghetti",
                                  "Individual + Group Mean"    = "mean_spaghetti",
                                  "Mean \u00b1 SEM"            = "mean_sem",
                                  "Mean \u00b1 SD"             = "mean_sd"),
                      selected = "spaghetti"
                    )
                  ),
                  conditionalPanel(ns = ns,
                    "output.tox_model_type_val == 'bw' && input.tox_bw_plot_selector == 'auc'",
                    hr(),
                    selectInput(ns("tox_auc_center_bar"), "Center bar",
                      choices  = c("Mean" = "mean", "Median" = "median", "None" = "none"),
                      selected = "mean"),
                    selectInput(ns("tox_auc_error_bar"), "Error bars",
                      choices  = c("SEM" = "SEM", "SD" = "SD", "95% CI" = "CI", "None" = "none"),
                      selected = "SEM"),
                    sliderInput(ns("tox_auc_pt_size"), "Point size",
                                min = 0.5, max = 6, value = 2, step = 0.5),
                    sliderInput(ns("tox_auc_jitter_width"), "Jitter width",
                                min = 0, max = 0.4, value = 0.15, step = 0.02)
                  ),
                  conditionalPanel(ns = ns, "output.tox_model_type_val == 'bw'",
                    ctrl_section("tox-groups", "sort", "Groups",
                      uiOutput(ns("tox_group_order_color_ui"))
                    )
                  ),
                  conditionalPanel(ns = ns, "output.tox_model_type_val == 'et'",
                    ctrl_section("tox-et-groups", "sort", "Groups",
                      uiOutput(ns("tox_et_group_order_color_ui"))
                    )
                  ),
                  conditionalPanel(ns = ns,
                    "output.tox_model_type_val == 'et' && input.tox_et_plot_selector == 'scatter'",
                    hr(),
                    checkboxInput(ns("tox_et_show_ellipses"),  "Show group ellipses",  value = TRUE),
                    checkboxInput(ns("tox_et_show_centroids"), "Show group centroids", value = TRUE),
                    sliderInput(ns("tox_et_scatter_pt_size"),  "Point size",       1,   6,  3,   step = 0.5),
                    sliderInput(ns("tox_et_centroid_size"),    "Centroid size",    2,  12,  6,   step = 0.5),
                    sliderInput(ns("tox_et_ellipse_lwd"),      "Ellipse line width", 0.3, 3, 0.8, step = 0.1),
                    selectInput(ns("tox_et_ellipse_lty"),      "Ellipse line style",
                      choices = c("Solid" = "solid", "Dashed" = "dashed", "Dotted" = "dotted"),
                      selected = "solid"),
                    checkboxInput(ns("tox_et_shape_by_group"), "Shape by group", value = FALSE)
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

  ordered_tox_grps <- reactive(
    get_ordered_groups(input, "tox", rv$tox_plot_groups)
  )
  ordered_et_grps <- reactive(
    get_ordered_groups(input, "tox_et", rv$tox_et_plot_groups)
  )

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
        removeTab(inputId = tabs_id, target = "Model Results")
        removeTab(inputId = tabs_id, target = "AUC")
        removeTab(inputId = tabs_id, target = "Weight Loss KM")
        bw_tabs_inserted(FALSE)
      }
      insertTab(inputId = tabs_id, target = "Summary", position = "after",
        tabPanel("Therapeutic Window",
          tab_info(
            strong("Therapeutic Window Metric (TWM). "),
            "TWM = TGI(%) / Mean Body Weight Loss(%). Higher values indicate a more ",
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
        tabPanel("Model Results",
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
      insertTab(inputId = tabs_id, target = "Model Results", position = "after",
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
    } else if (mt == "bayes_twm") {
      if (isTRUE(bw_tabs_inserted())) {
        removeTab(inputId = tabs_id, target = "Model Results")
        removeTab(inputId = tabs_id, target = "AUC")
        removeTab(inputId = tabs_id, target = "Weight Loss KM")
        bw_tabs_inserted(FALSE)
      }
      if (isTRUE(et_tabs_inserted())) {
        removeTab(inputId = tabs_id, target = "Therapeutic Window")
        removeTab(inputId = tabs_id, target = "Safety-Efficacy")
        removeTab(inputId = tabs_id, target = "Benefit Score")
        removeTab(inputId = tabs_id, target = "Corrected TGI")
        removeTab(inputId = tabs_id, target = "Combined Report")
        et_tabs_inserted(FALSE)
      }
      insertTab(inputId = tabs_id, target = "Summary", position = "after",
        tabPanel("Bayesian TWM Results",
          tab_info(
            strong("Bayesian Therapeutic Window Metric. "),
            "Combines posterior distributions of tumor growth inhibition and body weight change ",
            "to estimate a credible interval for the therapeutic window."
          ),
          div(class = "text-right", style = "margin-bottom: 6px;",
            downloadButton(ns("download_bayes_twm_csv"), "Download Results CSV")
          ),
          shinycssloaders::withSpinner(DTOutput(ns("bayes_twm_table"))),
          br(),
          shinycssloaders::withSpinner(
            plotOutput(ns("bayes_twm_plot"), height = "450px")
          ),
          downloadButton(ns("download_bayes_twm_plot"), "Download Plot")
        )
      )
    }
  }, ignoreNULL = FALSE)

  # Bayesian TWM readiness banner
  output$bayes_twm_readiness_ui <- renderUI({
    has_tg <- !is.null(rv$bayes_models$tg)
    has_bw <- !is.null(rv$bayes_models$bw)
    if (has_tg && has_bw) {
      status_alert("Bayesian TG and BW results available. Click Run Analysis.", "success")
    } else {
      missing <- c(if (!has_tg) "Bayesian TG (run in Tumor Growth tab)", if (!has_bw) "Bayesian BW (run above with Bayesian method)")
      status_alert(paste0("Required: ", paste(missing, collapse = "; ")), "warning")
    }
  })

  # ── Bayesian BW tab management ─────────────────────────────────────────────
  bayes_bw_tabs_inserted <- reactiveVal(FALSE)

  observeEvent(list(input$tox_model_type, input$tox_bw_method), {
    tabs_id  <- "tox_results_tabs"
    ns_here  <- session$ns
    is_bayes <- isTRUE(input$tox_model_type == "bw") &&
                isTRUE(input$tox_bw_method   == "bayes")
    if (is_bayes && !isTRUE(bayes_bw_tabs_inserted())) {
      appendTab(inputId = tabs_id,
        tabPanel("Posterior Check",
          tab_info(
            strong("Posterior predictive check (PPC)"), " — model-simulated weight trajectories overlaid on observed data.",
            em("Good fit: simulated densities envelope the observed density.")
          ),
          div(class = "text-right", style = "margin-bottom: 6px;",
            downloadButton(ns_here("download_tox_bw_pp_check"), "Download Plot")
          ),
          shinycssloaders::withSpinner(
            plotOutput(ns_here("tox_bw_pp_check_plot"), height = "450px")
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
            DTOutput(ns_here("tox_bw_mcmc_diag_table"))
          ),
          br(),
          shinycssloaders::withSpinner(
            plotOutput(ns_here("tox_bw_mcmc_trace_plot"), height = "500px")
          )
        ),
        select = FALSE
      )
      appendTab(inputId = tabs_id,
        tabPanel("Prior vs Posterior",
          tab_info(
            strong("Prior vs Posterior density overlay"), " — shows how much the data have updated each parameter's prior distribution."
          ),
          div(class = "text-right", style = "margin-bottom: 6px;",
            downloadButton(ns_here("download_tox_bw_prior_posterior"), "Download Plot")
          ),
          shinycssloaders::withSpinner(
            plotOutput(ns_here("tox_bw_prior_posterior_plot"), height = "400px")
          )
        ),
        select = FALSE
      )
      bayes_bw_tabs_inserted(TRUE)
    } else if (!is_bayes && isTRUE(bayes_bw_tabs_inserted())) {
      removeTab(inputId = tabs_id, target = "Posterior Check")
      removeTab(inputId = tabs_id, target = "MCMC Diagnostics")
      removeTab(inputId = tabs_id, target = "Prior vs Posterior")
      bayes_bw_tabs_inserted(FALSE)
    }
  }, ignoreNULL = FALSE)

  output$tox_bw_pp_check_plot <- renderPlot({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bw_bayes"))
    req(!isTRUE(rv$tox_results$mixed_model$error))
    req(!is.null(rv$bayes_models$bw))
    brms::pp_check(rv$bayes_models$bw, type = "dens_overlay", ndraws = 50)
  })

  output$tox_bw_mcmc_diag_table <- renderDT({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bw_bayes"))
    diag <- rv$tox_results$mixed_model$mcmc_diagnostics
    req(!is.null(diag))
    datatable(format_df_for_display(as.data.frame(diag)),
              options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })

  output$tox_bw_mcmc_trace_plot <- renderPlot({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bw_bayes"))
    req(!is.null(rv$bayes_models$bw))
    bayesplot::mcmc_trace(brms::as.array(rv$bayes_models$bw))
  })

  output$tox_bw_prior_posterior_plot <- renderPlot({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bw_bayes"))
    p <- rv$tox_results$mixed_model$prior_posterior_plot
    req(!is.null(p))
    print(p)
  })

  output$download_tox_bw_pp_check <- downloadHandler(
    filename = function() paste0("bw_posterior_check_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$bayes_models$bw))
      p <- brms::pp_check(rv$bayes_models$bw, type = "dens_overlay", ndraws = 50)
      ggplot2::ggsave(file, plot = p, width = 8, height = 5, dpi = 150)
    }
  )

  output$download_tox_bw_prior_posterior <- downloadHandler(
    filename = function() paste0("bw_prior_posterior_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$tox_results$mixed_model$prior_posterior_plot))
      ggplot2::ggsave(file, plot = rv$tox_results$mixed_model$prior_posterior_plot,
                      width = 8, height = 5, dpi = 150)
    }
  )

  # ── Bayesian TWM renders ───────────────────────────────────────────────────
  output$bayes_twm_table <- renderDT({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bayes_twm"))
    twm <- rv$tox_results$twm
    req(!is.null(twm), !isTRUE(twm$error))
    tbl <- tryCatch(as.data.frame(twm$twm_table %||% twm$treatment_effects), error = function(e) NULL)
    req(!is.null(tbl))
    datatable(format_df_for_display(tbl),
              options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })

  output$bayes_twm_plot <- renderPlot({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bayes_twm"))
    twm <- rv$tox_results$twm
    req(!is.null(twm), !isTRUE(twm$error))
    p <- twm$twm_plot %||% twm$plot
    req(!is.null(p))
    print(p)
  })

  output$download_bayes_twm_csv <- downloadHandler(
    filename = function() paste0("bayesian_twm_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      req(!is.null(rv$tox_results$twm))
      tbl <- rv$tox_results$twm$twm_table %||% rv$tox_results$twm$treatment_effects
      req(!is.null(tbl))
      write.csv(as.data.frame(tbl), file, row.names = FALSE)
    }
  )

  output$download_bayes_twm_plot <- downloadHandler(
    filename = function() paste0("bayesian_twm_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$tox_results$twm))
      p <- rv$tox_results$twm$twm_plot %||% rv$tox_results$twm$plot
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 8, height = 5, dpi = 150)
    }
  )

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

    if (model_type == "bayes_twm") {
      has_tg_bayes <- !is.null(rv$bayes_models$tg)
      has_bw_bayes <- !is.null(rv$bayes_models$bw)
      if (!has_tg_bayes || !has_bw_bayes) {
        missing <- c(if (!has_tg_bayes) "Bayesian Tumor Growth", if (!has_bw_bayes) "Bayesian Body Weight")
        output$tox_status_msg <- renderUI({
          status_alert(paste0("Run first: ", paste(missing, collapse = " and "), "."), "danger")
        })
        return()
      }
      if (!all(c(requireNamespace("brms",     quietly = TRUE),
                 requireNamespace("future",   quietly = TRUE),
                 requireNamespace("promises", quietly = TRUE)))) {
        output$tox_status_msg <- renderUI({
          status_alert(paste0("Bayesian mode requires 'brms', 'future', and 'promises'. ",
                              "Install with: install.packages(c('brms', 'future', 'promises'))"),
                       "danger")
        })
        return()
      }

      tg_for_twm <- list(model = rv$bayes_models$tg)
      bw_for_twm <- list(model = rv$bayes_models$bw)
      twm_params <- list(
        tg_result  = tg_for_twm,
        bw_result  = bw_for_twm,
        tx_col     = rv$tx_col,
        time_col   = rv$time_col,
        ref_group  = input$tox_reference_group %||% NULL
      )

      progress_twm <- shiny::Progress$new(session, min = 0, max = 1, style = "old")
      progress_twm$set(message = "Starting Bayesian TWM fit…", value = 0.05)

      output$tox_status_msg <- renderUI({
        tagList(
          h4("Bayesian TWM running…", style = "margin-top: 18px; font-weight: 600;"),
          p("Combining TG and BW posterior samples. This may take a few minutes.",
            style = "color:#666;"),
          shinycssloaders::withSpinner(div())
        )
      })

      promises::future_promise({
        mouseExperiment::bayesian_therapeutic_window(
          tg_result        = twm_params$tg_result,
          bw_result        = twm_params$bw_result,
          treatment_column = twm_params$tx_col,
          time_column      = twm_params$time_col,
          reference_group  = twm_params$ref_group,
          plots            = TRUE
        )
      }) %...>% (function(result) {
        progress_twm$close()
        result$model_type <- "bayes_twm"
        rv$tox_results <- list(
          model_type = "bayes_twm",
          twm        = result
        )
        if (isTRUE(result$error)) {
          output$tox_status_msg <- renderUI({
            status_alert(paste("Bayesian TWM failed:", result$message), "danger")
          })
        } else {
          output$tox_status_msg <- renderUI({
            status_alert("Bayesian Therapeutic Window analysis completed!", "success")
          })
        }
      }) %...!% (function(err) {
        progress_twm$close()
        err_msg <- translate_error(conditionMessage(err))
        output$tox_status_msg <- renderUI({
          status_alert(paste("Bayesian TWM failed:", err_msg), "danger")
        })
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

    progress <- shiny::Progress$new(session, min = 0, max = 1, style = "old")
    progress$set(message = "Validating data...", value = 0.1)
    on.exit(progress$close(), add = TRUE)

    results            <- list()
    results$model_type <- model_type

    if (model_type == "bw" && isTRUE(input$tox_bw_method == "bayes")) {
      if (!all(c(requireNamespace("brms",     quietly = TRUE),
                 requireNamespace("future",   quietly = TRUE),
                 requireNamespace("promises", quietly = TRUE)))) {
        output$tox_status_msg <- renderUI({
          status_alert(paste0("Bayesian mode requires 'brms', 'future', and 'promises'. ",
                              "Install with: install.packages(c('brms', 'future', 'promises'))"),
                       "danger")
        })
        return()
      }
      progress$close()

      bw_bayes_params <- list(
        data          = data,
        wt_col        = wt_col,
        time_col      = time_col,
        tx_col        = tx_col,
        id_col        = id_col,
        cage_col      = cage_col,
        vol_col       = if (has_vol) vol_col else NULL,
        ref_group     = ref_group,
        adjust        = adjust,
        density       = density,
        prior_str     = input$tox_bw_prior_strength %||% "skeptical",
        re_spec       = input$tox_bw_random_effects %||% "intercept_only",
        n_chains      = as.integer(input$tox_bw_bayes_chains  %||% 4L),
        n_warmup      = as.integer(input$tox_bw_bayes_warmup  %||% 1000L),
        n_iter        = as.integer(input$tox_bw_bayes_iter    %||% 500L)
      )

      progress_bw <- shiny::Progress$new(session, min = 0, max = 1, style = "old")
      progress_bw$set(message = "Starting Bayesian body weight fit…", value = 0.05)

      output$tox_status_msg <- renderUI({
        tagList(
          h4("Bayesian model running…",
             style = "margin-top: 18px; font-weight: 600;"),
          p("Fitting Bayesian body weight model with brms. This may take 3–12 minutes.",
            style = "color:#666;"),
          shinycssloaders::withSpinner(div())
        )
      })

      promises::future_promise({
        mouseExperiment::bayesian_body_weight(
          df                           = bw_bayes_params$data,
          weight_column                = bw_bayes_params$wt_col,
          time_column                  = bw_bayes_params$time_col,
          treatment_column             = bw_bayes_params$tx_col,
          id_column                    = bw_bayes_params$id_col,
          cage_column                  = bw_bayes_params$cage_col,
          volume_column                = bw_bayes_params$vol_col,
          adjust_tumor_weight          = bw_bayes_params$adjust,
          tumor_density                = bw_bayes_params$density,
          random_effects_specification = bw_bayes_params$re_spec,
          reference_group              = bw_bayes_params$ref_group,
          prior_strength               = bw_bayes_params$prior_str,
          n_chains                     = bw_bayes_params$n_chains,
          n_warmup                     = bw_bayes_params$n_warmup,
          n_iter                       = bw_bayes_params$n_iter,
          return_model                 = TRUE,
          plots                        = TRUE
        )
      }) %...>% (function(result) {
        progress_bw$close()
        result$model_type  <- "bw_bayes"
        if (!is.null(rv$bayes_models)) rv$bayes_models$bw <- result$model
        else                           rv$bayes_models <- list(bw = result$model)
        result$model <- NULL
        gc()
        rv$tox_results <- list(
          model_type    = "bw_bayes",
          mixed_model   = result,
          bw_auc        = NULL,
          weight_loss   = NULL
        )
        rv$tox_plot_groups <- if (!is.null(result$treatment_effects) &&
                                   has_col(as.data.frame(result$treatment_effects), "Group")) {
          unique(as.character(as.data.frame(result$treatment_effects)$Group))
        } else NULL
        if (isTRUE(result$error)) {
          err_msg <- result$message
          output$tox_status_msg <- renderUI({
            status_alert(paste("Bayesian BW failed:", err_msg), "danger")
          })
        } else {
          output$tox_status_msg <- renderUI({
            status_alert("Bayesian body weight analysis completed!", "success")
          })
        }
      }) %...!% (function(err) {
        progress_bw$close()
        err_msg <- translate_error(conditionMessage(err))
        output$tox_status_msg <- renderUI({
          status_alert(paste("Bayesian BW failed:", err_msg), "danger")
        })
      })
      return()
    }

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

    rv$tox_et_plot_groups <- if (model_type == "et" &&
                                 !is.null(results$bivariate) &&
                                 !isTRUE(results$bivariate$error) &&
                                 !is.null(results$bivariate$per_mouse))
      unique(as.character(results$bivariate$per_mouse$Treatment))
    else if (model_type == "et" && !is.null(tx_col) && has_col(data, tx_col))
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
      lines <- c(lines, "Therapeutic Window Metric (TWM = TGI / Mean Weight Loss%):", "")
      for (i in seq_len(nrow(tbl)))
        lines <- c(lines, sprintf("  %s:  TWM = %.2f  |  TGI = %.1f%%  |  Mean Weight Loss = %.1f%%",
                                  tbl$Treatment[i], tbl$TWM[i],
                                  tbl$TGI[i], tbl$Mean_Pct_Weight_Loss[i]))
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

  output$tox_et_group_order_color_ui <- renderUI({
    req(rv$tox_et_plot_groups)
    group_order_and_color_ui("tox_et", rv$tox_et_plot_groups, ns = session$ns)
  })
  observe_palette_preset(input, session, "tox_et", reactive(rv$tox_et_plot_groups))

  output$tox_trajectory_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$trajectory_data)
    td <- rv$tox_results$trajectory_data
    req(nrow(td) > 0)
    traj_type    <- input$tox_traj_type %||% "spaghetti"
    text_sz      <- input$tox_plot_text_size %||% 14
    line_w       <- input$tox_traj_line_width %||% 1
    ordered_grps <- ordered_tox_grps()
    colors       <- build_color_vector(input, "tox_color", ordered_grps)

    if (traj_type %in% c("spaghetti", "mean_spaghetti")) {
      td$.__LineGroup <- paste(td[["Treatment"]], td[["ID"]], sep = "\u00b7")
      if (!is.null(ordered_grps))
        td$Treatment <- factor(td$Treatment, levels = ordered_grps)
      p <- ggplot2::ggplot(td, ggplot2::aes(.data[["Day"]], .data[["Weight"]],
                                            group = .data[[".__LineGroup"]], color = .data[["Treatment"]])) +
        ggplot2::geom_line(alpha = if (traj_type == "mean_spaghetti") 0.3 else 0.6,
                           linewidth = line_w) +
        ggplot2::geom_point(size = 1.5, alpha = if (traj_type == "mean_spaghetti") 0.3 else 0.6) +
        ggplot2::labs(title = "Body Weight Trajectories",
                     x = "Day", y = "Weight (g)", color = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
      if (traj_type == "mean_spaghetti") {
        grp_mean <- stats::aggregate(Weight ~ Treatment + Day, data = td, FUN = mean, na.rm = TRUE)
        if (!is.null(ordered_grps))
          grp_mean$Treatment <- factor(grp_mean$Treatment, levels = ordered_grps)
        p <- p + ggplot2::geom_line(
          data = grp_mean,
          ggplot2::aes(.data[["Day"]], .data[["Weight"]], group = .data[["Treatment"]],
                       color = .data[["Treatment"]]),
          linewidth = line_w * 2.2, inherit.aes = FALSE)
      }
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
    text_sz      <- input$tox_plot_text_size %||% 14
    center       <- input$tox_auc_center_bar %||% "mean"
    err_bars     <- input$tox_auc_error_bar  %||% "SEM"
    pt_sz        <- input$tox_auc_pt_size    %||% 2
    jitter_w     <- input$tox_auc_jitter_width %||% 0.15
    ordered_grps <- ordered_tox_grps()
    colors       <- build_color_vector(input, "tox_color", ordered_grps)
    auc_df       <- auc$auc_per_mouse
    if (!is.null(ordered_grps))
      auc_df$Treatment <- factor(auc_df$Treatment, levels = ordered_grps)

    # Per-group summary for center bar and error bars
    gs <- stats::aggregate(AUC_Pct_Change ~ Treatment, data = auc_df,
      FUN = function(x) c(mean   = mean(x, na.rm = TRUE),
                          median = stats::median(x, na.rm = TRUE),
                          sd     = stats::sd(x, na.rm = TRUE),
                          n      = sum(!is.na(x))))
    gs <- do.call(data.frame, gs)
    names(gs) <- c("Treatment", "Mean", "Median", "SD", "N")
    gs$SEM  <- gs$SD / sqrt(gs$N)
    gs$CI95 <- gs$SEM * stats::qt(0.975, pmax(gs$N - 1, 1))
    if (!is.null(ordered_grps))
      gs$Treatment <- factor(gs$Treatment, levels = ordered_grps)
    gs$center_y <- if (center == "median") gs$Median else gs$Mean

    p <- ggplot2::ggplot(auc_df,
      ggplot2::aes(x = .data[["Treatment"]], y = .data[["AUC_Pct_Change"]])) +
      ggplot2::geom_jitter(ggplot2::aes(color = .data[["Treatment"]]),
                           width = jitter_w, alpha = 0.8, size = pt_sz)

    if (center != "none") {
      p <- p + ggplot2::geom_crossbar(
        data = gs,
        ggplot2::aes(x = .data[["Treatment"]], y = .data[["center_y"]],
                     ymin = .data[["center_y"]], ymax = .data[["center_y"]]),
        width = 0.45, linewidth = 0.8, color = "black", fill = NA,
        inherit.aes = FALSE)
      if (err_bars != "none") {
        err_half <- switch(err_bars,
          SEM = gs$SEM, SD = gs$SD, CI = gs$CI95, NULL)
        if (!is.null(err_half)) {
          gs$ymin_err <- gs$center_y - err_half
          gs$ymax_err <- gs$center_y + err_half
          p <- p + ggplot2::geom_errorbar(
            data = gs,
            ggplot2::aes(x = .data[["Treatment"]], ymin = .data[["ymin_err"]],
                         ymax = .data[["ymax_err"]]),
            width = 0.2, linewidth = 0.8, color = "black",
            inherit.aes = FALSE)
        }
      }
    }

    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Body Weight AUC (% Change from Baseline)",
                   x = "Treatment", y = "AUC of % Weight Change (\u00b7days)") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "none")
    if (!is.null(colors)) p <- p +
      ggplot2::scale_color_manual(values = colors)
    p
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_km_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(NULL)
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- ordered_tox_grps()
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
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- ordered_et_grps()
    colors       <- build_color_vector(input, "tox_et_color", ordered_grps)
    tbl <- twm$twm_table
    if (!is.null(ordered_grps))
      tbl$Treatment <- factor(tbl$Treatment, levels = ordered_grps)
    else
      tbl$Treatment <- factor(tbl$Treatment, levels = tbl$Treatment[order(tbl$TWM)])
    p <- ggplot2::ggplot(tbl,
           ggplot2::aes(.data[["Treatment"]], .data[["TWM"]], fill = .data[["Treatment"]])) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["TWM"]], 2)), vjust = -0.5, size = 4) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = "Therapeutic Window Metric by Treatment",
                   subtitle = "TWM = TGI(%) / Mean Weight Loss(%). Dashed line = TWM of 1.",
                   x = "Treatment", y = "TWM") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
    p
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_scatter_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bivariate)
    biv <- rv$tox_results$bivariate
    if (isTRUE(biv$error)) return(NULL)
    pm           <- biv$per_mouse
    text_sz      <- input$tox_plot_text_size %||% 14
    pt_size      <- input$tox_et_scatter_pt_size %||% 3
    cent_size    <- input$tox_et_centroid_size   %||% 6
    ell_lwd      <- input$tox_et_ellipse_lwd     %||% 0.8
    ell_lty      <- input$tox_et_ellipse_lty     %||% "solid"
    show_ell     <- isTRUE(input$tox_et_show_ellipses)
    show_cent    <- isTRUE(input$tox_et_show_centroids)
    shape_by_grp <- isTRUE(input$tox_et_shape_by_group)
    ordered_grps <- ordered_et_grps()
    colors       <- build_color_vector(input, "tox_et_color", ordered_grps)
    if (!is.null(ordered_grps))
      pm$Treatment <- factor(pm$Treatment, levels = ordered_grps)
    eff_label <- switch(biv$efficacy_metric,
      tgi = "TGI (%)", tumor_auc = "Tumor AUC", "Unknown metric")
    aes_base <- if (shape_by_grp)
      ggplot2::aes(.data[["Efficacy"]], .data[["Max_Pct_Weight_Loss"]],
                   color = .data[["Treatment"]], shape = .data[["Treatment"]])
    else
      ggplot2::aes(.data[["Efficacy"]], .data[["Max_Pct_Weight_Loss"]],
                   color = .data[["Treatment"]])
    p <- ggplot2::ggplot(pm, aes_base) +
      ggplot2::geom_point(size = pt_size, alpha = 0.7)
    if (show_ell && nrow(pm) > 3)
      p <- p + ggplot2::stat_ellipse(level = 0.68, linewidth = ell_lwd,
                                     linetype = ell_lty, show.legend = FALSE)
    if (show_cent && !is.null(biv$per_group)) {
      pg <- biv$per_group
      if (!is.null(ordered_grps))
        pg$Treatment <- factor(pg$Treatment, levels = ordered_grps)
      p <- p + ggplot2::geom_point(data = pg,
                 ggplot2::aes(.data[["Efficacy_Mean"]], .data[["Toxicity_Mean"]],
                              color = .data[["Treatment"]]),
                 shape = 18, size = cent_size, inherit.aes = FALSE)
    }
    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Safety-Efficacy Scatter",
                   subtitle = paste0("Each point is one mouse. Efficacy: ", eff_label),
                   x = eff_label, y = "Max Weight Loss (%)", color = "Treatment") +
      ggplot2::theme_classic(base_size = text_sz)
    if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)
    p
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_benefit_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$benefit)
    ben <- rv$tox_results$benefit
    if (isTRUE(ben$error)) return(NULL)
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- ordered_et_grps()
    bt <- ben$benefit_table
    if (!is.null(ordered_grps))
      bt$Treatment <- factor(bt$Treatment, levels = ordered_grps)
    else
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
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- ordered_et_grps()
    comp <- wc$comparison
    if (!is.null(ordered_grps))
      comp$Treatment <- factor(comp$Treatment, levels = ordered_grps)
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
          "efficacy_toxicity_plot.png")
      }
    },
    content = function(file) {
      w   <- input$tox_plot_width  %||% 800
      h   <- input$tox_plot_height %||% 500
      dpi <- input$tox_plot_dpi    %||% 100
      model <- input$tox_model_type %||% "bw"
      if (model == "bw" && (input$tox_bw_plot_selector %||% "") == "km") {
        p <- tryCatch(output$tox_km_plot(), error = function(e) NULL)
        if (is.null(p)) { showNotification("Plot is not available.", type = "warning"); return(invisible(NULL)) }
        png(file, width = w, height = h, units = "px", res = dpi)
        print(p); dev.off(); return()
      }
      plot_id <- if (model == "bw") {
        switch(input$tox_bw_plot_selector %||% "trajectory",
          trajectory = "tox_trajectory_plot", auc = "tox_auc_plot", "tox_trajectory_plot")
      } else {
        switch(input$tox_et_plot_selector %||% "twm",
          twm = "tox_et_twm_plot", scatter = "tox_et_scatter_plot",
          benefit = "tox_et_benefit_plot", wctgi = "tox_et_wctgi_plot",
          "tox_et_twm_plot")
      }
      p <- tryCatch(output[[plot_id]](), error = function(e) NULL)
      if (is.null(p)) {
        showNotification("Plot is not available.", type = "warning")
        return(invisible(NULL))
      }
      png(file, width = w, height = h, units = "px", res = dpi)
      print(p)
      dev.off()
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
