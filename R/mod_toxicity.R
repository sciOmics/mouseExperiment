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
                selectInput(ns("tox_bw_bayes_backend"), "Stan backend",
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
                          min = 5, max = 40, value = 20, step = 1),
              hr(),
              # CODE_REVIEW.md backend R3.5 / G.3 + R3.6 / R3.7.
              h5(strong("Endpoint & Uncertainty")),
              selectInput(ns("tox_endpoint_method"), "Endpoint estimand",
                choices = c(
                  "Model-based (recommended)" = "model",
                  "Each animal's last observation" = "last_obs",
                  "Survivors at the endpoint day (legacy)" = "survivors"),
                selected = "model"),
              helpText(
                "Animals leave these studies because their tumours grew, so ",
                "restricting to survivors selects the slower growers \u2014 hardest ",
                "in the control arm, which biases TGI downward. ",
                strong("Model-based"), " instead reads each arm's geometric mean ",
                "at the endpoint day from a log-scale mixed model fitted to every ",
                "observation, so nothing is discarded. Use ", strong("Survivors"),
                " only to reproduce a pre-0.8.0 figure."),
              numericInput(ns("tox_n_boot"), "Bootstrap resamples (0 = skip)",
                          value = 2000, min = 0, max = 20000, step = 500),
              helpText("Mouse-level bootstrap giving 95% intervals for TGI, mean ",
                       "weight loss and TWM. The TWM table is sorted, so the ",
                       "intervals are what tell you whether the ranking is real.")
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
                  ),
                  # CODE_REVIEW.md backend R3.26 -- censoring here is informative:
                  # an animal whose record ends early usually ended because it was
                  # euthanised for tumour burden, which shares a cause with weight
                  # loss. KM assumes a censored animal remains at risk, so 1 - KM
                  # overstates weight-loss incidence whenever that happens.
                  uiOutput(ns("tox_km_competing_note")),
                  DTOutput(ns("tox_km_censoring_table"))
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
                strong("Body-weight model diagnostics. "),
                "Frequentist (LME4): residual diagnostics. ",
                "Bayesian: MCMC convergence + posterior diagnostics. ",
                "The view switches based on the Statistical Model selected."
              ),
              # ── Freq diagnostics (LME4 BW) ─────────────────────────────
              conditionalPanel(ns = ns,
                              "input.tox_bw_method != 'bayes'",
              fluidRow(
                column(6,
                  h4("Q-Q Plot of Residuals"),
                  div(style = "margin-bottom: 24px;",
                    shinycssloaders::withSpinner(
                      plotOutput(ns("tox_diag_qq"), height = "380px")
                    )
                  ),
                  tags$details(
                    tags$summary(tags$small("How to read this plot")),
                    tags$small(
                      "Points should fall along the red reference line. ",
                      strong("S-shape: "), "heavy-tailed residuals — outliers or a ",
                      "transform issue. ",
                      strong("Bow: "), "skewed residuals — usually fixed by changing the ",
                      "weight transform."
                    )
                  )
                ),
                column(6,
                  h4("Residuals vs Fitted"),
                  div(style = "margin-bottom: 24px;",
                    shinycssloaders::withSpinner(
                      plotOutput(ns("tox_diag_resid_fitted"), height = "380px")
                    )
                  ),
                  tags$details(
                    tags$summary(tags$small("How to read this plot")),
                    tags$small(
                      "Random scatter around zero = LMM assumptions met. ",
                      strong("Curved blue line: "), "the linear-in-time model misses ",
                      "non-linearity — consider a polynomial term or a GAM. ",
                      strong("Fan shape: "), "heteroscedasticity — consider log-",
                      "transforming or weighting the response."
                    )
                  )
                )
              ),
              br(),
              fluidRow(
                column(6,
                  h4("Scale-Location"),
                  div(style = "margin-bottom: 24px;",
                    shinycssloaders::withSpinner(
                      plotOutput(ns("tox_diag_scale_location"), height = "380px")
                    )
                  ),
                  tags$details(
                    tags$summary(tags$small("How to read this plot")),
                    tags$small(
                      "sqrt(|residuals|) vs fitted — magnifies variance trends. ",
                      "Flat blue line = constant variance. Upward slope = ",
                      "heteroscedasticity."
                    )
                  )
                ),
                column(6,
                  h4("Random-effects Q-Q"),
                  div(style = "margin-bottom: 24px;",
                    shinycssloaders::withSpinner(
                      plotOutput(ns("tox_diag_re_qq"), height = "380px")
                    )
                  ),
                  tags$details(
                    tags$summary(tags$small("How to read this plot")),
                    tags$small(
                      "Per-animal random intercepts should be normally ",
                      "distributed. Heavy tails are rare; safe to ignore at ",
                      "typical preclinical N."
                    )
                  )
                )
              ),
              # CODE_REVIEW.md DIAGNOSTICS gap (13) — LMM influence
              # diagnostics for BW.
              conditionalPanel(
                ns = ns,
                condition = "output.tox_show_influence === true",
                hr(),
                h4("Influence diagnostics (case-deletion)"),
                tab_info(
                  strong("Cook's distance"), " measures the overall ",
                  "effect of deleting each observation. Above the dashed ",
                  "threshold (4/n) = an observation that moves the fit. ",
                  strong("DFBETAS"), " shows per-coefficient influence; ",
                  "values | > 1 | are flagged."
                ),
                fluidRow(
                  column(6,
                    h5("Cook's distance"),
                    shinycssloaders::withSpinner(
                      DTOutput(ns("tox_cooks_distance"))
                    )
                  ),
                  column(6,
                    h5("DFBETAS"),
                    shinycssloaders::withSpinner(
                      DTOutput(ns("tox_dfbetas"))
                    )
                  )
                )
              )
              ),  # /freq BW conditionalPanel
              # ── Bayes MCMC diagnostics (BW Bayesian) ──────────────────
              conditionalPanel(ns = ns,
                              "input.tox_bw_method == 'bayes'",
                tab_info(
                  strong("MCMC convergence diagnostics. "),
                  "All Rhat values should be ≤ 1.01 before interpreting results."
                ),
                uiOutput(ns("tox_bw_bayes_diagnostics")),
                shinycssloaders::withSpinner(DTOutput(ns("tox_bw_mcmc_diag_table"))),
                br(),
                h5(strong("Influential observations (Pareto-k > 0.7)")),
                em(class = "text-muted",
                   "Observations the LOO approximation can't reliably estimate."),
                shinycssloaders::withSpinner(DTOutput(ns("tox_bw_bayes_influential"))),
                br(),
                shinycssloaders::withSpinner(
                  plotOutput(ns("tox_bw_mcmc_trace_plot"), height = "500px"))
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

  # D.7: this was an unconditional `rv$tox_results <- NULL`. `rv` is created empty
  # in app.R and each module mounts once, so the assignment never cleared anything
  # in production -- but it wiped the slot at module init, which made these five
  # modules impossible to drive through testServer() with a pre-populated result.
  # That is why I.2 stayed open across four rounds: every attempt to assert render
  # behaviour got the "run an analysis first" placeholder back. Removed; the slot
  # reads NULL until an analysis writes it either way.

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
          # CODE_REVIEW.md R3.D9 -- the table is sorted by TWM, so it reads as a
          # ranking. Backend v0.7.0 attaches bootstrap bounds; without them the
          # order looks more resolved than the data supports.
          uiOutput(ns("tox_twm_ranking_note")),
          shinycssloaders::withSpinner(DTOutput(ns("tox_et_twm_table"))),
          uiOutput(ns("tox_twm_attrition_ui")),
          DTOutput(ns("tox_twm_attrition_table")),
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
      # CODE_REVIEW.md Round 2 H — mE v0.4.6 G.5 added aliased diagnostic
      # plots from the two input models (pp_check_tg/bw, prior_posterior_tg/bw,
      # mcmc_trace_tg/bw). Surface them in a single Diagnostics tab so users
      # don't have to dig back into the TG/BW results to inspect convergence
      # of the models that fed into TWM.
      insertTab(inputId = tabs_id, target = "Bayesian TWM Results",
                position = "after",
        tabPanel("TWM Diagnostics",
          tab_info(
            strong("Diagnostics for the two input models. "),
            "TWM does not fit its own brms model — it combines posterior ",
            "draws from a Bayesian TG fit and a Bayesian BW fit. Each pair ",
            "of plots below shows the same standard diagnostic (PPC, prior ",
            "vs posterior, MCMC trace) for both inputs side-by-side."
          ),
          h4("Posterior Predictive Check"),
          fluidRow(
            column(6, h5("TG model"),
              shinycssloaders::withSpinner(plotOutput(ns("twm_pp_check_tg"),
                                                      height = "300px"))),
            column(6, h5("BW model"),
              shinycssloaders::withSpinner(plotOutput(ns("twm_pp_check_bw"),
                                                      height = "300px")))
          ),
          br(),
          h4("Prior vs Posterior"),
          fluidRow(
            column(6, h5("TG model"),
              shinycssloaders::withSpinner(plotOutput(ns("twm_prior_post_tg"),
                                                      height = "300px"))),
            column(6, h5("BW model"),
              shinycssloaders::withSpinner(plotOutput(ns("twm_prior_post_bw"),
                                                      height = "300px")))
          ),
          br(),
          h4("MCMC Trace"),
          fluidRow(
            column(6, h5("TG model"),
              shinycssloaders::withSpinner(plotOutput(ns("twm_trace_tg"),
                                                      height = "300px"))),
            column(6, h5("BW model"),
              shinycssloaders::withSpinner(plotOutput(ns("twm_trace_bw"),
                                                      height = "300px")))
          )
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

  output$tox_bw_bayes_diagnostics <- renderUI({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bw_bayes"))
    bayes_diagnostics_panel(rv$tox_results$mixed_model)
  })

  output$tox_bw_bayes_influential <- renderDT({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bw_bayes"))
    inf_df <- bayes_influential_obs(rv$tox_results$mixed_model)
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

  output$tox_bw_mcmc_trace_plot <- renderPlot({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bw_bayes"))
    req(!is.null(rv$bayes_models$bw))
    bayesplot::mcmc_trace(posterior::as_draws_array(rv$bayes_models$bw))
  })

  output$tox_bw_prior_posterior_plot <- renderPlot({
    req(rv$tox_results)
    req(isTRUE(rv$tox_results$model_type == "bw_bayes"))
    p <- rv$tox_results$mixed_model$prior_posterior_plot
    req(!is.null(p))
    print(p)
  })

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

  # CODE_REVIEW.md Round 2 H — surface the six aliased diagnostic plots
  # that mE v0.4.6 G.5 added to bayesian_therapeutic_window() output.
  twm_plot_renderer <- function(field_name) {
    renderPlot({
      req(rv$tox_results)
      req(isTRUE(rv$tox_results$model_type == "bayes_twm"))
      req(!is.null(rv$tox_results$twm))
      p <- rv$tox_results$twm[[field_name]]
      if (is.null(p)) {
        plot.new()
        title(main = paste0(field_name, " not available"),
              sub = "Re-run the upstream TG/BW Bayesian fits with plots = TRUE.",
              cex.sub = 0.8)
        return(invisible(NULL))
      }
      print(p)
    })
  }
  output$twm_pp_check_tg    <- twm_plot_renderer("pp_check_tg")
  output$twm_pp_check_bw    <- twm_plot_renderer("pp_check_bw")
  output$twm_prior_post_tg  <- twm_plot_renderer("prior_posterior_tg")
  output$twm_prior_post_bw  <- twm_plot_renderer("prior_posterior_bw")
  output$twm_trace_tg       <- twm_plot_renderer("mcmc_trace_tg")
  output$twm_trace_bw       <- twm_plot_renderer("mcmc_trace_bw")

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

      tg_for_twm <- list(model = rv$bayes_models$tg)
      bw_for_twm <- list(model = rv$bayes_models$bw)
      twm_params <- list(
        tg_result  = tg_for_twm,
        bw_result  = bw_for_twm,
        tx_col     = rv$tx_col,
        time_col   = rv$time_col,
        ref_group  = input$tox_reference_group %||% NULL
      )

      # F.7 — route through the job registry.
      launch_twm_bayes_job(rv$session_id, twm_params,
                           input, output, session, rv)
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

    show_analysis_modal("Running toxicity analysis…", "Results will appear automatically when complete.")
    on.exit(shiny::removeModal(), add = TRUE)

    results            <- list()
    results$model_type <- model_type

    if (model_type == "bw" && isTRUE(input$tox_bw_method == "bayes")) {
      # Cancel the frequentist on.exit(removeModal) registered above — otherwise it
      # tears down the Bayesian modal shown below the moment this branch returns.
      # The future_promise callbacks remove the modal when the fit completes.
      on.exit()

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
        n_iter        = as.integer(input$tox_bw_bayes_iter    %||% 500L),
        backend       = input$tox_bw_bayes_backend %||% "rstan"
      )

      # F.7 — route through the job registry.
      launch_bw_bayes_job(rv$session_id, bw_bayes_params,
                          input, output, session, rv)
      return()
    }

    if (model_type == "bw") {
      covs         <- input$tox_covariates   %||% character(0)
      est_method   <- input$tox_estimation   %||% "REML"
      threshold    <- (input$tox_threshold   %||% 20) / 100
      baseline_raw <- input$tox_baseline_day %||% "first"
      baseline_day <- if (baseline_raw == "first") NULL else as.numeric(baseline_raw)

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

      results$bw_auc <- safe_analysis(module = "tox_auc", {
        mouseExperiment::body_weight_auc(
          df = data, weight_column = wt_col, time_column = time_col,
          treatment_column = tx_col, id_column = id_col,
          volume_column = if (has_vol) vol_col else NULL,
          adjust_tumor_weight = adjust, tumor_density = density, reference_group = ref_group
        )
      })

      results$weight_loss <- safe_analysis(module = "tox_threshold", {
        mouseExperiment::weight_loss_threshold(
          df = data, weight_column = wt_col, time_column = time_col,
          treatment_column = tx_col, id_column = id_col,
          # Backend R3.25: cage is part of the mouse key, so reused IDs across
          # cages are not collapsed into one subject.
          cage_column = cage_col,
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

      # Backend v0.8.0+ / v0.7.0+: the endpoint estimand and the mouse-level
      # bootstrap. Defaults are the statistically correct ones; "survivors"
      # exists only to reproduce pre-0.8.0 figures. CODE_REVIEW.md R3.D9.
      ep_method <- input$tox_endpoint_method %||% "model"
      n_boot_tw <- as.integer(input$tox_n_boot %||% 2000L)

      results$twm <- safe_analysis(module = "et_twm", {
        mouseExperiment::therapeutic_window_metric(
          df = data, volume_column = vol_col, weight_column = wt_col,
          time_column = time_col, treatment_column = tx_col, id_column = id_col,
          adjust_tumor_weight = adjust, tumor_density = density, reference_group = ref_group,
          endpoint_method = ep_method, n_boot = n_boot_tw, boot_seed = 20260730L
        )
      })

      results$bivariate <- safe_analysis(module = "et_bivariate", {
        mouseExperiment::efficacy_toxicity_bivariate(
          df = data, volume_column = vol_col, weight_column = wt_col,
          time_column = time_col, treatment_column = tx_col, id_column = id_col,
          adjust_tumor_weight = adjust, tumor_density = density,
          efficacy_metric = eff_metric, reference_group = ref_group,
          endpoint_method = ep_method
        )
      })

      results$benefit <- safe_analysis(module = "et_benefit", {
        mouseExperiment::total_benefit_area(
          df = data, volume_column = vol_col, weight_column = wt_col,
          time_column = time_col, treatment_column = tx_col, id_column = id_col,
          adjust_tumor_weight = adjust, tumor_density = density,
          lambda = lambda, reference_group = ref_group
        )
      })

      results$wctgi <- safe_analysis(module = "et_wctgi", {
        mouseExperiment::weight_corrected_tgi(
          df = data, volume_column = vol_col, weight_column = wt_col,
          time_column = time_col, treatment_column = tx_col, id_column = id_col,
          adjust_tumor_weight = adjust, tumor_density = density,
          safety_threshold = threshold, reference_group = ref_group,
          endpoint_method = ep_method
        )
      })

      n_err   <- sum(vapply(results[c("twm", "bivariate", "benefit", "wctgi")],
                             function(x) isTRUE(x$error), logical(1)))
      n_total <- 4L
    }

    rv$tox_results <- results
    updateTabsetPanel(session, "tox_results_tabs", selected = "Summary")

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
    elapsed_str <- format_elapsed(res$elapsed_seconds)
    if (!is.null(elapsed_str)) {
      lines <- c(lines, "", paste("COMPUTATION TIME:", elapsed_str))
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

  # Does the TWM ranking survive its own uncertainty? With typical group sizes the
  # intervals overlap heavily and the order of adjacent arms is usually not
  # resolvable -- worth saying before anyone acts on the top row.
  output$tox_twm_ranking_note <- renderUI({
    req(rv$tox_results, rv$tox_results$twm)
    twm <- rv$tox_results$twm
    if (isTRUE(twm$error)) return(NULL)
    tbl <- twm$twm_table
    if (is.null(tbl) || !all(c("TWM_Lower", "TWM_Upper") %in% names(tbl)) ||
        nrow(tbl) < 2L) return(NULL)

    ord <- tbl[order(-tbl$TWM), , drop = FALSE]
    top <- ord[1, ]; second <- ord[2, ]
    separable <- is.finite(top$TWM_Lower) && is.finite(second$TWM_Upper) &&
      top$TWM_Lower > second$TWM_Upper

    status_alert(
      if (separable) {
        sprintf(paste0("%s ranks highest and its 95%% interval does not overlap ",
                       "the next arm's, so the top of this ranking is resolved."),
                top$Treatment)
      } else {
        sprintf(paste0("The 95%% intervals for %s and %s overlap, so this ranking ",
                       "does not distinguish them. Read the order as indicative ",
                       "and check TWM_Lower / TWM_Upper before acting on it."),
                top$Treatment, second$Treatment)
      },
      if (separable) "success" else "warning")
  })

  output$tox_twm_attrition_ui <- renderUI({
    req(rv$tox_results, rv$tox_results$twm)
    twm <- rv$tox_results$twm
    if (isTRUE(twm$error) || is.null(twm$attrition)) return(NULL)
    a <- twm$attrition
    lost <- if (!is.null(a$Pct_Lost)) suppressWarnings(max(a$Pct_Lost, na.rm = TRUE)) else 0
    tagList(
      tags$hr(),
      h4("Animals contributing at the endpoint day"),
      p(class = "text-muted",
        "Endpoint day: ", strong(twm$endpoint_day %||% "\u2014"),
        ". Estimand: ", strong(twm$endpoint_method %||% "model"), "."),
      if (is.finite(lost) && lost > 0 &&
          identical(twm$endpoint_method, "survivors")) {
        status_alert(sprintf(paste0(
          "Up to %.0f%% of animals in one arm did not reach the endpoint day, and ",
          "the survivor estimand conditions on reaching it. Because animals leave ",
          "when their tumours grow, this biases TGI downward -- most in the ",
          "control arm. Switch to the model-based estimand unless you are ",
          "reproducing a pre-0.8.0 figure."), lost), "warning")
      }
    )
  })

  output$tox_twm_attrition_table <- renderDT({
    req(rv$tox_results, rv$tox_results$twm)
    twm <- rv$tox_results$twm
    if (isTRUE(twm$error) || is.null(twm$attrition)) return(NULL)
    datatable(format_df_for_display(as.data.frame(twm$attrition)),
              options = list(paging = FALSE, searching = FALSE), rownames = FALSE)
  })

  output$tox_km_competing_note <- renderUI({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(NULL)
    n_cr <- wl$n_competing_risk %||% 0L
    if (!is.finite(n_cr) || n_cr == 0L) {
      return(tagList(
        tags$hr(),
        p(class = "text-muted",
          "All censoring here is administrative (animals still on study at the ",
          "end), so the Kaplan-Meier curve is unbiased for this endpoint.")
      ))
    }
    tagList(
      tags$hr(),
      h4("Competing risks"),
      status_alert(sprintf(paste0(
        "%d animal(s) left the study before the end without reaching the ",
        "weight-loss threshold. That censoring is NOT independent \u2014 animals are ",
        "usually removed for tumour burden, which shares a cause with weight loss. ",
        "Kaplan-Meier assumes a censored animal stays at risk, so the curve above ",
        "OVERSTATES weight-loss incidence. Use the Aalen-Johansen cumulative ",
        "incidence in the result object (`cuminc`) for incidence statements."),
        n_cr), "warning"),
      p(class = "text-muted",
        "Breakdown of how each animal left the study:")
    )
  })

  output$tox_km_censoring_table <- renderDT({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error) || is.null(wl$censoring_summary)) return(NULL)
    cs <- as.data.frame(wl$censoring_summary)
    if (!nrow(cs)) return(NULL)
    datatable(format_df_for_display(cs),
              options = list(paging = FALSE, searching = FALSE), rownames = FALSE)
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

  # Frequentist plot renderers + diagnostics + group-order UI live in
  # `mod_toxicity_plots.R` (CODE_REVIEW.md C.1/C.3).
  add_tox_plots(input, output, session, rv, ordered_tox_grps, ordered_et_grps)

  # All 17 download handlers (Bayesian BW + Bayesian TWM + frequentist CSV
  # + active plot + ZIP) live in `mod_toxicity_downloads.R` (CODE_REVIEW.md
  # C.1/C.3). Active-plot handler reaches into `output[[plot_id]]()` which
  # still works because every handler binds to the same `output` env.
  add_tox_download_handlers(input, output, session, rv)

  # F.7 — register resume handlers for both Bayesian fits in this module.
  register_bw_resume(input, output, session, rv)
  register_twm_resume(input, output, session, rv)

  }) # end moduleServer
}
