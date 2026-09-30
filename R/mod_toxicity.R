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
                "Body Weight Analysis" = "bw",
                "Therapeutic Window"   = "et"
              ),
              selected = "bw"
            ),
            helpText(
              "Body Weight Analysis models weight trajectories with a mixed-effects ",
              "model and times each animal to a weight-loss threshold. ",
              "Therapeutic Window also requires a tumor volume column and sets each ",
              "arm's tumor growth inhibition against its weight loss."
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

            # Therapeutic-window options
            conditionalPanel(ns = ns, "input.tox_model_type == 'et'",
              uiOutput(ns("tox_et_vol_status")),
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
                      tags$li(strong("Body Weight Analysis"), " — models longitudinal weight trajectories using a linear mixed-effects model, and times each animal to a weight-loss threshold (Kaplan-Meier)."),
                      tags$li(strong("Therapeutic Window"), " — sets each arm's tumor growth inhibition against its body-weight loss. Requires a tumor volume column.")
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
                    h5(strong("Weight Loss KM Curve")),
                    p("Treats each mouse as an event-time observation where the event is crossing the weight-loss threshold (default 20%). ",
                      "Mice that never reach the threshold are censored at their last observation. ",
                      "A Kaplan-Meier curve is fitted per group; group differences are tested with the log-rank test."),
                    div(class = "alert alert-warning", style = "font-size:0.9em;",
                      icon("triangle-exclamation"), " ",
                      strong("Competing risks."), " A mouse euthanised for tumour burden can no longer experience a ",
                      "weight-loss event, so treating it as an ordinary censoring overstates the weight-loss risk — ",
                      "censoring assumes the animal remained at risk, and it did not. Where tumour-driven removals ",
                      "are present the panel reports how many there were and provides an Aalen-Johansen cumulative ",
                      "incidence estimate alongside the Kaplan-Meier curve. The censoring summary breaks down, per ",
                      "arm, how many observations were true censorings versus competing events; a large competing ",
                      "count means the KM curve should not be read as a probability of weight loss."),
                    p(tags$b("Weight Loss Threshold:"),
                      " The % body weight loss that constitutes a toxicity event. The default (20%) is a common IACUC endpoint; adjust to match study-specific welfare criteria.")
                  )
                ),

                # ── Therapeutic Window ──────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "Therapeutic Window"), " \u2014 Analyses & Metrics")),
                  div(class = "panel-body",

                    h5(strong("Efficacy")),
                    p("Efficacy is tumor growth inhibition at the endpoint day, relative to the reference group: ",
                      tags$code("TGI = (1 \u2212 V_treated / V_control) \u00d7 100"),
                      ". Values > 0 indicate tumor suppression; > 100 indicates regression. ",
                      "Each arm's endpoint volume comes from the endpoint estimand described below."),

                    hr(),
                    h5(strong("Therapeutic Window Metric (TWM)")),
                    p("A summary ratio combining efficacy and toxicity per group:"),
                    tags$blockquote(tags$code("TWM = TGI(%) / Mean Body Weight Loss(%)")),
                    p("Higher TWM indicates a more favourable balance. The denominator is the group mean of per-mouse maximum weight loss (not the single worst-case mouse). A noise floor prevents division by near-zero weight changes. Displayed as a ranked bar chart."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Ranking bounds"),
                      tags$dd("TWM is a ratio of two estimated quantities, so its uncertainty is wider than either ",
                              "component's and does not follow a normal distribution. Bootstrap intervals are ",
                              "reported and the ranking note states whether the top arm is genuinely separated from ",
                              "the second or whether their intervals overlap. A bar chart ordered by point estimate ",
                              "looks like a ranking whether or not the data support one."),
                      tags$dt("Bootstrap resamples"),
                      tags$dd("Number of mouse-level resamples used for those intervals. The control arm is resampled ",
                              "too, so uncertainty in the TGI denominator propagates rather than being treated as ",
                              "known. Set to 0 to skip, which returns point estimates with no interval."),
                      tags$dt("Endpoint estimand"),
                      tags$dd(strong("Model"), " (default) estimates group means at the endpoint from the fitted ",
                              "model, using every observation each animal contributed. ", strong("Observed"),
                              " uses only animals still present at the endpoint. The two differ whenever animals are ",
                              "removed early, and the difference is not small: because the fastest-growing controls ",
                              "are the ones that leave, the observed estimand flatters the control arm and understates ",
                              "TGI. Per-arm attrition is reported beneath the table so the size of the issue is visible.")
                    ),

                    br()
                  )
                ),

                # ── Shared Parameters ────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Shared Parameters")),
                  div(class = "panel-body",
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Reference Group"),
                      tags$dd("The vehicle/control group used as the denominator for TGI. Auto-detected from common names (DMSO, Vehicle, Control, PBS, Saline) or defaults to the first group alphabetically."),
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
              "Therapeutic Window: TGI, weight loss and TWM per arm."
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
                        "Weight Loss KM Curve"    = "km"
                      ),
                      selected = "trajectory"
                    )
                  ),
                  conditionalPanel(ns = ns, "output.tox_model_type_val == 'et'",
                    selectInput(ns("tox_et_plot_selector"), "Plot Type",
                      choices = c(
                        "Therapeutic Window"      = "twm"
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
                  "output.tox_model_type_val == 'et' && input.tox_et_plot_selector == 'twm'",
                  h4("Therapeutic Window Metric"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tox_et_twm_plot"), width = "auto", height = "auto")
                  )
                )
              )
            )
          ),

          tabPanel("Diagnostics",
            conditionalPanel(ns = ns, "output.tox_model_type_val == 'bw'",
              tab_info(
                strong("Body-weight model diagnostics. "),
                "Residual, random-effects and influence diagnostics for the ",
                "body-weight mixed model."
              ),
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
                      strong("Curved blue line: "), "weight does not change linearly ",
                      "with time, so the Day slope averages over a curved ",
                      "trajectory — read it alongside the trajectory plot. ",
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
      et_tabs_inserted(TRUE)

    } else if (mt == "bw" && !isTRUE(bw_tabs_inserted())) {
      if (isTRUE(et_tabs_inserted())) {
        removeTab(inputId = tabs_id, target = "Therapeutic Window")
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
      return(status_alert("Therapeutic Window also requires a tumor volume column.", "warning"))
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
        status_alert("Therapeutic Window requires a tumor volume column.", "danger")
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

    show_analysis_modal("Running toxicity analysis…", "Results will appear automatically when complete.")
    on.exit(shiny::removeModal(), add = TRUE)

    results            <- list()
    results$model_type <- model_type

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

      n_err   <- sum(vapply(results[c("mixed_model", "weight_loss")],
                             function(x) isTRUE(x$error), logical(1)))
      n_total <- 2L

    } else {
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

      n_err   <- as.integer(isTRUE(results$twm$error))
      n_total <- 1L
    }

    rv$tox_results <- results
    updateTabsetPanel(session, "tox_results_tabs", selected = "Summary")

    # Store treatment groups for colour/order pickers
    rv$tox_plot_groups <- if (!is.null(results$trajectory_data) && nrow(results$trajectory_data) > 0)
      unique(as.character(results$trajectory_data$Treatment))
    else if (!is.null(tx_col) && has_col(data, tx_col))
      unique(as.character(data[[tx_col]]))
    else NULL

    # Data order (so colours match the other tabs), limited to the arms the
    # TWM table reports.
    rv$tox_et_plot_groups <- if (model_type == "et" && !is.null(tx_col) &&
                                 has_col(data, tx_col)) {
      grps    <- unique(as.character(data[[tx_col]]))
      twm_tbl <- results$twm$twm_table
      if (!isTRUE(results$twm$error) && !is.null(twm_tbl))
        grps <- grps[grps %in% twm_tbl$Treatment]
      grps
    } else NULL

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
    lines <- c("=== THERAPEUTIC WINDOW SUMMARY ===", "")
    if (!is.null(res$twm) && !isTRUE(res$twm$error)) {
      tbl   <- res$twm$twm_table
      lines <- c(lines, "Therapeutic Window Metric (TWM = TGI / Mean Weight Loss%):", "")
      for (i in seq_len(nrow(tbl)))
        lines <- c(lines, sprintf("  %s:  TWM = %.2f  |  TGI = %.1f%%  |  Mean Weight Loss = %.1f%%",
                                  tbl$Treatment[i], tbl$TWM[i],
                                  tbl$TGI[i], tbl$Mean_Pct_Weight_Loss[i]))
      lines <- c(lines, "")
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

  # Plot renderers + diagnostics + group-order UI live in
  # `mod_toxicity_plots.R` (CODE_REVIEW.md C.1/C.3).
  add_tox_plots(input, output, session, rv, ordered_tox_grps, ordered_et_grps)

  # The six download handlers (CSV + active plot + ZIP) live in
  # `mod_toxicity_downloads.R` (CODE_REVIEW.md C.1/C.3).
  add_tox_download_handlers(input, output, session, rv)

  }) # end moduleServer
}
