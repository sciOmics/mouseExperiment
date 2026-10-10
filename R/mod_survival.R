#' Survival Analysis UI
#' @keywords internal
survivalUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Survival Analysis",
       tags$small(class = "text-muted", style = "font-size:0.55em; margin-left:0.6em;",
                  "(Time-to-Endpoint)")),


    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Analysis Parameters"),

          selectInput(ns("surv_model_type"), "Statistical Model",
            choices = c(
              "Frequentist (Cox PH)"    = "frequentist",
              "Bayesian (Weibull/brms)" = "bayes"
            ),
            selected = "frequentist"
          ),
          helpText("Cox PH fits a proportional-hazards model. Bayesian Weibull ",
                   "models the full survival time distribution and reports posterior ",
                   "credible intervals; it runs in the background and takes from ",
                   "seconds to a few minutes."),

          # The reference is the study's control group (R20-D.H).
          uiOutput(ns("surv_reference_note")),
          uiOutput(ns("surv_design_note")),

          conditionalPanel(ns = ns, "input.surv_model_type == 'bayes'",
            hr(),
            h5(strong("Bayesian Settings")),
            selectInput(ns("surv_bayes_family"), "Survival family",
              choices = c(
                "Weibull"     = "weibull",
                "Log-normal"  = "lognormal"
              ),
              selected = "weibull"
            ),
            selectInput(ns("surv_prior_strength"), "Prior Strength",
              choices = c(
                "Skeptical (recommended)" = "skeptical",
                "Weakly informative"      = "weakly_informative",
                "Informative"             = "informative",
                "Diffuse"                 = "diffuse"
              ),
              selected = "skeptical"
            ),
            numericInput(ns("surv_bayes_chains"), "MCMC chains",
                         value = 4, min = 1, max = 8, step = 1),
            numericInput(ns("surv_bayes_warmup"), "Warmup per chain",
                         value = 1000, min = 200, max = 5000, step = 100),
            numericInput(ns("surv_bayes_iter"), "Post-warmup draws per chain",
                         value = 500, min = 200, max = 5000, step = 100),
            selectInput(ns("surv_bayes_backend"), "Stan backend",
                        choices = c("rstan", "cmdstanr"),
                        selected = "rstan"),
            div(class = "alert alert-warning",
                style = "font-size:0.85em; padding:8px; margin-top:6px;",
              icon("clock"), " The first Bayesian fit of a model compiles it (about 20 s); ",
              "a later fit with the same structure reuses it and usually takes under a minute. ",
              "Fits run one at a time on a background worker, so other users and tabs stay ",
              "responsive, and yours may wait for one ahead of it. You can close the browser ",
              "and return later; the result waits for you, and the banner at the top lets ",
              "you resume or cancel the fit."
            )
          ),

          hr(),
          h5(strong("Endpoint")),
          uiOutput(ns("surv_endpoint_source_ui")),

          hr(),
          # R20.D13: a "Treat natural death as event" checkbox sat here and
          # nothing read it. It works now, and appears only when a
          # removal-reason column says which animals died (R20-D.H);
          # otherwise the tab states what it assumes.
          uiOutput(ns("surv_death_option_ui")),
          uiOutput(ns("surv_event_assumption")),

          br(),
          actionButton(ns("run_surv_analysis"),
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),

          uiOutput(ns("surv_status_msg"))
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
          id = ns("surv_results_tabs"),

          tabPanel("Info",
            fluidRow(
              column(10, offset = 1,
                br(),

                # ── Note on "survival" semantics ────────────────────────────
                div(class = "alert alert-warning",
                    style = "font-size:0.88em; margin-bottom:12px; padding:10px 14px;",
                  icon("triangle-exclamation"), " ",
                  strong("Note on ‘survival’ in mouse tumour experiments: "),
                  "Most IACUC protocols require euthanasia once tumour volume exceeds a set limit ",
                  "(e.g. 2000 mm³). The recorded endpoint is therefore ",
                  strong("time until tumour volume threshold"),
                  ", not true death. This is mathematically identical to survival analysis, but results ",
                  "should be labelled as ‘Time-to-Endpoint’ or ‘Time-to-Tumour-Volume-Threshold’ ",
                  "rather than ‘survival’ in publications.",
                  tags$br(),
                  tags$span(style = "color:#555; font-size:0.95em;",
                    "Occasionally a mouse dies before reaching the volume limit. ",
                    "Competing-risks models (e.g. Fine–Gray) are the rigorous solution, but they require ",
                    "large per-group sample sizes rarely achievable in mouse experiments. ",
                    "The recommended pragmatic approach is to treat any death as an event (‘endpoint reached’), ",
                    "because the animal did not survive the study regardless of cause. ",
                    "Map a removal-reason column on Data Upload and the sidebar offers exactly that. ",
                    "Without one, with an Event/Censor column events are exactly the rows coded 1, so code such ",
                    "deaths 1 to follow this approach; when events are derived from a volume threshold, only a ",
                    "tumour reaching the limit is an event and other removals are censored. The sidebar says ",
                    "which applies."
                  )
                ),

                # ── Overview ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Overview")),
                  div(class = "panel-body",
                    p("The Survival tab fits Kaplan–Meier (KM) survival curves and a Cox proportional hazards model to time-to-endpoint data. In mouse tumor experiments the 'endpoint' is typically the day a tumor exceeds the IACUC volume limit, not true animal death."),
                    p("Three outputs are produced:"),
                    tags$ul(
                      tags$li(strong("Kaplan–Meier curves"), ", with the log-rank test across all arms (not adjusted for multiple comparisons)."),
                      tags$li(strong("Cox proportional hazards model"), " reporting hazard ratios (HR) per group relative to the reference."),
                      tags$li(strong("Risk table"), " showing animals at risk, events, and KM estimates at each event time.")
                    )
                  )
                ),

                # ── KM ───────────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-info", "Kaplan–Meier Estimator"))),
                  div(class = "panel-body",
                    p("The KM estimator computes a non-parametric step-function survival probability S(t) = P(T > t) from event times. Each downward step corresponds to one or more events on that day."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Log-rank test"),
                      tags$dd("Tests whether the survival curves of two or more groups are the same. The plot shows the test across all arms at once, unadjusted. The comparisons of each arm with the reference are in Model Results: Cox (or Firth) p-values, adjusted for the number of comparisons; when the Cox model cannot be fitted, they come from pairwise log-rank tests."),
                      tags$dt("Median survival"),
                      tags$dd("Day by which 50 % of animals reached the endpoint. 'Not reached' if fewer than half had events before study end."),
                      tags$dt("Censoring"),
                      tags$dd("Animals that never reach the endpoint contribute their last observation day as a censored time. This correctly uses partial information rather than excluding them.")
                    ),
                    p(tags$small(class = "text-muted",
                      "When an arm has no events, its hazard ratio is estimated with Firth's correction, which stays finite. ",
                      "If that fit fails, the comparisons come from pairwise log-rank tests."
                    ))
                  )
                ),

                # ── Cox ───────────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "Cox Proportional Hazards Model"))),
                  div(class = "panel-body",
                    p("The Cox model estimates the hazard ratio (HR) for each treatment group relative to the reference:"),
                    tags$blockquote(tags$code("h(t) = h₀(t) × exp(β₁ × Treatment)")),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("HR < 1"),
                      tags$dd("Lower hazard than reference — animals reach the endpoint more slowly (better outcome)."),
                      tags$dt("HR = 1"),
                      tags$dd("Identical hazard to reference."),
                      tags$dt("HR > 1"),
                      tags$dd("Higher hazard than reference — animals reach the endpoint faster (worse outcome).")
                    ),
                    p(strong("Firth's bias-reduced estimation"), " is applied automatically when any group has zero events. Standard Cox partial likelihood produces infinite coefficient estimates in this case; Firth's penalisation regularises them."),
                    p(tags$small(class = "text-muted",
                      "The proportional hazards assumption requires that the HR is constant over time. ",
                      "Crossing KM curves indicate a time-varying HR, which violates this assumption — interpret Cox results with caution in that case."
                    ))
                  )
                ),

                # ── Endpoint ─────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Endpoint Configuration")),
                  div(class = "panel-body",
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Mapped censor column"),
                      tags$dd("When an Event/Censor column is mapped in Data Upload (0 = censored, 1 = event), it is used directly."),
                      tags$dt("No censor column"),
                      tags$dd("Without a mapped Event/Censor column, choose 'Derive from tumour volume threshold' as the event source. An animal's disappearance is not treated as an event by itself."),
                      tags$dt("Volume threshold"),
                      tags$dd("Optionally generate an event indicator by flagging animals whose tumor volume exceeds a specified threshold at any timepoint."),
                      tags$dt("Natural death"),
                      tags$dd("With a removal-reason column mapped on Data Upload, the sidebar offers to count deaths from other causes as events (recommended): an animal whose reason is marked 'Death or moribund' under Study Design is an event at its last measurement, with either event source. Without one, or with the box unticked: Event/Censor column, an animal that died of other causes is an event only if its row is coded 1; volume threshold, it is censored at its last measurement, which assumes its later risk was like that of the animals still on study; if such deaths are treatment-related, this favours the treatment."),
                      tags$dt("Event time"),
                      tags$dd("The first day an animal's row is flagged 1 (or its tumour reaches the limit). An animal whose every row is flagged has an animal-level status, and its last day is used. Animals are identified by treatment, ID and cage, so an ear tag reused in another arm is a different animal.")
                    )
                  )
                ),

                br()
              )
            )
          ),

          tabPanel("Summary",
            tab_info(strong("Time-to-Endpoint summary: "), "Per-group event counts, total subjects, event rate, and median time-to-endpoint (the day by which 50\u202f% of animals reached the endpoint \u2014 typically tumour volume exceeding the IACUC limit). \u2018Not reached\u2019 means fewer than half the group experienced an event by study end."),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("surv_summary"))
            )
          ),
          
          tabPanel("Model Results",
            tab_info(
              strong("Cox proportional hazards model. "),
              "Hazard Ratios (HR) are computed relative to the reference group. Firth\u2019s bias-reduced method is applied automatically when a group has zero events.",
              tags$ul(
                tags$li(strong("Group"), " \u2014 treatment group name."),
                tags$li(strong("HR"), " \u2014 Hazard Ratio: the instantaneous risk of the event in this group relative to the reference. HR\u202f<\u202f1 means lower hazard (better survival); HR\u202f>\u202f1 means higher hazard."),
                tags$li(strong("CI_Lower / CI_Upper"), " \u2014 95\u202f% confidence interval for the HR, not adjusted for multiple comparisons. It can exclude 1.0 while the adjusted P_Value is above 0.05; judge significance by P_Value."),
                tags$li(strong("P_Value"), " \u2014 p-value for the comparison with the reference, adjusted for the number of comparisons (P_Adjust_Method; unadjusted in P_Value_Unadjusted)."),
                tags$li(strong("P_Method"), " \u2014 shown when whole cages were randomised (Study Design on Data Upload): the p-value comes from a cage-level permutation log-rank, and the hazard ratio has no interval."),
                tags$li(strong("Events"), " \u2014 number of animals that reached the endpoint (tumour volume exceeded the IACUC threshold, or died before reaching it if natural death is treated as an event)."),
                tags$li(strong("Total"), " \u2014 total number of animals in the group."),
                tags$li(strong("Event_Rate"), " \u2014 Events \u00f7 Total. Values close to 1.0 mean most animals reached the endpoint."),
                tags$li(strong("Median_Survival"), " \u2014 day by which 50\u202f% of the group reached the endpoint. \u2018Not reached\u2019 if fewer than half had events by study end."),
                tags$li(strong("Note"), " \u2014 flags the reference group row.")
              )
            ),
            div(class = "text-right", style = "margin-bottom: 6px;",
              downloadButton(ns("download_surv_cox"), "Download Results")
            ),
            # C-index (mE v0.4.5 E.2) + Cox PH check summary (mE v0.4.5 A.1)
            uiOutput(ns("surv_model_quality_panel")),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("surv_cox_summary"))
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("surv_cox_table"))
            )
          ),

          tabPanel("Diagnostics",
            tab_info(
              strong("Model-fit diagnostics. "),
              "Cox PH: proportional-hazards assumption check (Schoenfeld). ",
              "Bayesian (brms): MCMC convergence + posterior diagnostics. ",
              "The view switches based on the selected Statistical Model."
            ),
            # ── Freq diagnostics (Cox PH) ──────────────────────────────────
            conditionalPanel(ns = ns,
                            "input.surv_model_type != 'bayes'",
            tab_info(
              strong("Cox proportional-hazards assumption check via Schoenfeld residuals (cox.zph). "),
              "The Cox model assumes the hazard ratio is constant over time. A small global p-value indicates the assumption is violated and the reported HRs / CIs may be misleading (a treatment whose effect accelerates or attenuates over time gives inconsistent HRs depending on when you look).",
              tags$ul(
                tags$li(strong("chisq"), " \u2014 Schoenfeld-residual test statistic per covariate, plus a GLOBAL row."),
                tags$li(strong("df"), " \u2014 degrees of freedom for the test."),
                tags$li(strong("p"), " \u2014 p-value. < 0.05 in the GLOBAL row means the PH assumption is rejected at the 5 % level.")
              ),
              em("Only available for the standard Cox path. NA for Firth (coxphf) and log-rank fallbacks.")
            ),
            uiOutput(ns("surv_ph_test_banner")),
            shinycssloaders::withSpinner(
              DTOutput(ns("surv_ph_test_table"))
            ),
            hr(),
            h4("Schoenfeld residuals plot"),
            tab_info(
              strong("How to read this plot. "),
              "Each panel shows the smoothed Schoenfeld residuals for one ",
              "covariate against transformed survival time. ",
              strong("Flat horizontal smoother around 0 (within the shaded CI band)"),
              " = PH assumption met for that term. ",
              strong("Monotonic trend"),
              " = the hazard ratio drifts over time (treatment effect ",
              "strengthens or weakens). ",
              strong("Sinusoidal / non-monotonic shape"),
              " = a more complex time-varying effect — usually means a ",
              "Cox model with a time-varying coefficient (or a parametric ",
              "AFT) fits better. The dashed line at zero is the null."
            ),
            shinycssloaders::withSpinner(
              plotOutput(ns("surv_schoenfeld_plot"), height = "500px")
            )
            ),  # /freq conditionalPanel
            # ── Bayes MCMC diagnostics ─────────────────────────────────────
            conditionalPanel(ns = ns,
                            "input.surv_model_type == 'bayes'",
              tab_info(
                strong("MCMC convergence diagnostics"),
                " for all fixed-effect parameters. All chains should have ",
                "converged (R̂ ≈ 1.00) and have sufficient effective ",
                "samples before interpreting results.",
                tags$ul(
                  tags$li(strong("Rhat"), " — potential scale-reduction factor; values > 1.01 indicate non-convergence."),
                  tags$li(strong("Bulk_ESS"), " — bulk effective sample size; aim for > 400 per chain."),
                  tags$li(strong("Tail_ESS"), " — tail ESS; important for credible interval accuracy."),
                  tags$li(strong("Converged"), " — TRUE when R̂ ≤ 1.01.")
                ),
                em("If Rhat > 1.01, try increasing iterations or chains before reporting results.")
              ),
              uiOutput(ns("surv_bayes_diagnostics")),
              shinycssloaders::withSpinner(DTOutput(ns("surv_mcmc_diag_table"))),
              br(),
              h5(strong("Influential observations (Pareto-k > 0.7)")),
              em(class = "text-muted",
                 "Observations the LOO approximation can't reliably estimate. Investigate these mice individually."),
              shinycssloaders::withSpinner(DTOutput(ns("surv_bayes_influential"))),
              br(),
              shinycssloaders::withSpinner(
                plotOutput(ns("surv_mcmc_trace_plot"), height = "500px"))
            )
          ),

          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  selectInput(ns("surv_plot_selector"), "Select Plot",
                    choices = c(
                      "Kaplan-Meier Curve"       = "km",
                      "Hazard Ratio Forest Plot" = "forest"
                    )
                  ),
                  ctrl_section("surv-groups", "sort", "Groups",
                    conditionalPanel(ns = ns,"input.surv_plot_selector == 'km'",
                      uiOutput(ns("surv_group_order_color_ui"))
                    ),
                    conditionalPanel(ns = ns,"input.surv_plot_selector == 'forest'",
                      tags$p(class = "text-muted", style = "font-size:12px; margin:4px 0;",
                             "Group colours are mapped automatically from the analysis groups.")
                    )
                  ),
                  ctrl_section("surv-elements", "layer-group", "Elements",
                    conditionalPanel(ns = ns,"input.surv_plot_selector == 'km'",
                      checkboxInput(ns("surv_opt_conf_int"),
                                    tagList("Show confidence interval",
                                            help_tip("Displays 95% confidence bands around the Kaplan-Meier survival curves.")),
                                    value = FALSE),
                      checkboxInput(ns("surv_opt_risk_table"),
                                    tagList("Show risk table on plot",
                                            help_tip("Number of subjects still at risk at each timepoint, displayed beneath the KM curve.")),
                                    value = TRUE),
                      numericInput(ns("surv_km_xstart"),
                                   "X-axis start (Day)",
                                   value = 0, min = 0, step = 1)
                    ),
                    conditionalPanel(ns = ns,"input.surv_plot_selector == 'forest'",
                      colourpicker::colourInput(ns("surv_forest_pt_color"),
                                               "Point colour", value = "steelblue"),
                      h6("Annotation"),
                      checkboxInput(ns("surv_forest_show_hr"), "Show HR",      value = TRUE),
                      checkboxInput(ns("surv_forest_show_ci"), "Show 95% CI",  value = TRUE),
                      checkboxInput(ns("surv_forest_show_p"),  "Show p-value", value = TRUE)
                    ),
                    selectInput(ns("surv_opt_legend_pos"), "Legend position",
                      choices = c("None"   = "none",
                                  "Top"    = "top",    "Bottom" = "bottom",
                                  "Left"   = "left",   "Right"  = "right",
                                  "In plot \u2013 top left"     = "inplot_tl",
                                  "In plot \u2013 top right"    = "inplot_tr",
                                  "In plot \u2013 bottom left"  = "inplot_bl",
                                  "In plot \u2013 bottom right" = "inplot_br"),
                      selected = "right"
                    ),
                    textInput(ns("surv_plot_title"), "Title",
                              placeholder = "Leave blank for default"),
                    selectInput(ns("surv_plot_title_align"), "Title alignment",
                      choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                      selected = "center"),
                    textInput(ns("surv_plot_subtitle"), "Subtitle",
                              placeholder = "Optional subtitle"),
                    textInput(ns("surv_plot_xaxis_label"), "X-axis label",
                              placeholder = "Leave blank for default"),
                    textInput(ns("surv_plot_yaxis_label"), "Y-axis label",
                              placeholder = "Leave blank for default"),
                    textInput(ns("surv_plot_legend_label"), "Legend label",
                              placeholder = "Leave blank for default")
                  ),
                  ctrl_section("surv-scaling", "expand", "Scaling",
                    sliderInput(ns("surv_opt_text_size"), "Text size (pt)",
                                min = 10, max = 24, value = 16, step = 1)
                  ),
                  ctrl_section("surv-export", "arrows-alt", "Plot Size / Export",
                    fluidRow(
                      column(6, numericInput(ns("surv_plot_width"),  "Width (px)",  800, min = 200, max = 5000, step = 50)),
                      column(6, numericInput(ns("surv_plot_height"), "Height (px)", 600, min = 200, max = 5000, step = 50))
                    ),
                    numericInput(ns("surv_plot_dpi"), "Resolution (DPI)", 100, min = 72, max = 600, step = 10),
                    tags$div(style = "margin-top: 10px;", uiOutput(ns("surv_download_plot_btn")))
                  )
                )
              ),
              column(9,
                conditionalPanel(ns = ns,"input.surv_plot_selector == 'km'",
                  h4("Kaplan-Meier Survival Curves"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("surv_plot_km"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.surv_plot_selector == 'forest'",
                  h4("Hazard Ratio Forest Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("surv_forest_plot"), width = "auto", height = "auto")
                  )
                )
              )
            )
          ),

          tabPanel("Risk Table",
            tab_info(
              strong("Numbers at risk table: "), "one row per group per event time point, summarising what happened to the cohort at that moment.",
              tags$ul(
                tags$li(strong("Group"), " \u2014 treatment group name."),
                tags$li(strong("Time"), " \u2014 the day on which one or more events occurred."),
                tags$li(strong("N_Risk"), " \u2014 number of animals still alive and uncensored at the start of this time interval."),
                tags$li(strong("N_Events"), " \u2014 number of animals that experienced the endpoint at exactly this time."),
                tags$li(strong("N_Censored"), " \u2014 number of animals whose last observation was at this time without reaching the endpoint (e.g. still within tumour volume limit at study end, or removed from study)."),
                tags$li(strong("Survival"), " \u2014 Kaplan\u2013Meier estimated survival probability at this time point (1.0\u202f=\u202f100\u202f% surviving).")
              ),
              "N_Risk decreases over time as events and censoring accumulate; a rapid drop indicates many events in a short window."
            ),
            div(class = "text-right", style = "margin-bottom: 6px;",
              downloadButton(ns("download_surv_risk_table"), "Download Table")
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("surv_risk_table"))
            )
          )
        )
      )
    )
  )
}

#' Survival Analysis Server Logic
#' @keywords internal
survival_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

  
  # D.7: this was an unconditional `rv$surv_results <- NULL`. `rv` is created empty
  # in app.R and each module mounts once, so the assignment never cleared anything
  # in production -- but it wiped the slot at module init, which made these five
  # modules impossible to drive through testServer() with a pre-populated result.
  # That is why I.2 stayed open across four rounds: every attempt to assert render
  # behaviour got the "run an analysis first" placeholder back. Removed; the slot
  # reads NULL until an analysis writes it either way.
  rv$surv_running <- FALSE

  ordered_surv_grps <- reactive(
    get_ordered_groups(input, "surv", rv$surv_plot_groups)
  )

  observe({
    toggle_downloads(
      session,
      c("download_surv_cox", "download_surv_risk_table"),
      enabled = !is.null(rv$surv_results) && !isTRUE(rv$surv_results$error)
    )
  })

  # The reference is the study's control group in this tab's labels
  # ("Treatment Dose" when a dose column is mapped).
  surv_ref <- reactive({
    req(rv$uploaded_data, rv$data_valid)
    study_reference(rv, analysis_group_labels(rv))
  })

  output$surv_reference_note <- renderUI({
    reference_note(rv, tryCatch(surv_ref(), error = function(e) NULL))
  })

  output$surv_design_note <- renderUI({
    unit <- rv$settings$randomisation_unit %||% "mouse"
    if (identical(unit, "cage")) {
      tags$p(style = "font-size:12px;", "Unit of randomisation: ", strong("whole cages"),
             ". Each arm is compared with the reference by moving whole cages ",
             "between them (a cage-level permutation log-rank). With 2 cages per ",
             "arm no p-value can fall below 1/3, and with 3 below 0.1.")
    } else {
      tags$p(style = "font-size:12px;", "Unit of randomisation: ",
             strong("individual animals"), " (Study Design on Data Upload).")
    }
  })

  output$surv_death_option_ui <- renderUI({
    spec <- removal_reason_spec(rv)
    if (is.null(spec)) return(NULL)
    deaths <- reasons_in(spec, "death")
    tagList(
      checkboxInput(session$ns("surv_death_as_event"),
                    "Count deaths from other causes as events (recommended)",
                    value = isolate(input$surv_death_as_event) %||% TRUE),
      tags$p(tags$small(class = "text-muted",
        "Deaths are the removal reasons marked 'Death or moribund' under Study ",
        "Design on Data Upload: ",
        if (length(deaths)) paste(sQuote(deaths, FALSE), collapse = ", ")
        else "none at present", "."))
    )
  })

  output$surv_event_assumption <- renderUI({
    src <- input$surv_event_source %||% "column"
    spec <- removal_reason_spec(rv)
    deaths_counted <- !is.null(spec) && isTRUE(input$surv_death_as_event)
    div(class = "alert alert-info", style = "font-size:0.85em; padding:8px;",
      icon("circle-info"), " ",
      if (deaths_counted) {
        tagList(strong("Deaths from other causes are events. "),
          if (identical(src, "volume"))
            "A tumour reaching the volume limit is an event, and so is an animal whose removal reason is a death, at its last measurement. Other removals are censored."
          else
            "Events are the rows coded 1, and any animal whose removal reason is a death, whatever the Event/Censor column says.")
      } else if (identical(src, "volume")) {
        tagList(strong("Deaths from other causes are censored. "),
          "Only a tumour reaching the volume limit is an event; an animal that ",
          "leaves the study for another reason is censored at its last ",
          "measurement. If such removals are treatment-related, this favours ",
          "the treatment.",
          if (is.null(spec)) " Map a removal-reason column on Data Upload to count deaths as events.")
      } else {
        tagList(strong("Events are the rows coded 1. "),
          "An animal that died of other causes counts as an event only if the ",
          "Event/Censor column codes it 1.",
          if (is.null(spec)) " Map a removal-reason column on Data Upload to count deaths as events.")
      })
  })

  # Endpoint source UI: always shows both options; warns inline if volume column not configured
  output$surv_endpoint_source_ui <- renderUI({
    has_vol <- nzchar(rv$settings$col_volume %||% "") ||
               (isTRUE(rv$settings$volume_method == "calculate") &&
                nzchar(rv$settings$col_length %||% "") &&
                nzchar(rv$settings$col_width  %||% ""))
    tagList(
      radioButtons(session$ns("surv_event_source"), NULL,
        choices = c(
          "Use Event/Censor column"             = "column",
          "Derive from tumour volume threshold" = "volume"
        ),
        selected = isolate(input$surv_event_source) %||% "column"
      ),
      conditionalPanel(ns = session$ns, "input.surv_event_source == 'column'",
        tags$p(tags$small(class = "text-muted",
          "Using the Event/Censor column configured in Data Upload. ",
          "Coding: 1\u202f=\u202fendpoint reached, 0\u202f=\u202fcensored."
        ))
      ),
      conditionalPanel(ns = session$ns, "input.surv_event_source == 'volume'",
        if (!has_vol)
          div(class = "alert alert-warning", style = "font-size:0.85em; padding:6px 10px; margin-bottom:6px;",
            icon("triangle-exclamation"), " ",
            "No Tumour Volume column is configured. Go to ",
            strong("Data Upload"), " and map the Volume column, then return here."
          ),
        numericInput(session$ns("surv_vol_threshold"),
                     paste0("IACUC volume limit (", volume_unit_label(rv), ")"),
                     value = isolate(input$surv_vol_threshold) %||%
                       if (identical(rv$settings$volume_units, "cm3")) 2 else 2000,
                     min = 0, step = if (identical(rv$settings$volume_units, "cm3")) 0.05 else 50),
        tags$p(tags$small(
          "First day each mouse\u2019s tumour reaches or exceeds this volume = Event\u202f=\u202f1 ",
          "(Time = that day). Mice that never exceed the limit are censored at their ",
          "last observation (Event\u202f=\u202f0)."
        ))
      )
    )
  })

  # Run survival analysis
  observeEvent(input$run_surv_analysis, {
    req(rv$uploaded_data, rv$data_valid)

    event_src <- input$surv_event_source %||% "column"

    # Validate that required columns are configured from data upload
    if (rv$settings$col_treatment == "" || rv$settings$col_treatment == "-- Select a reference group --") {
      output$surv_status_msg <- renderUI({
        status_alert("Error: Please select a Treatment Group on the Data Upload tab", "danger")
      })
      shinyjs::runjs("window.trackError('validation', 'Treatment column not selected', 'survival');")
      return()
    }

    if (event_src == "column") {
      if (is.null(rv$settings$col_survival_event) || rv$settings$col_survival_event == "" || rv$settings$col_survival_event == "None") {
        output$surv_status_msg <- renderUI({
          status_alert("Error: Please select an Event/Censor Status column in Data Upload, or switch to \u2018Derive from tumour volume threshold\u2019.", "danger")
        })
        shinyjs::runjs("window.trackError('validation', 'Survival event column not selected', 'survival');")
        return()
      }
    } else {
      has_vol_surv <- nzchar(rv$settings$col_volume %||% "") ||
                      (isTRUE(rv$settings$volume_method == "calculate") &&
                       nzchar(rv$settings$col_length %||% "") &&
                       nzchar(rv$settings$col_width  %||% ""))
      if (!has_vol_surv) {
        output$surv_status_msg <- renderUI({
          status_alert("Error: No Tumour Volume column configured. Configure it in Data Upload or switch to \u2018Use Event/Censor column\u2019.", "danger")
        })
        return()
      }
    }
    
    ref_group <- surv_ref()
    if (is.null(ref_group) || !nzchar(ref_group)) {
      output$surv_status_msg <- renderUI(status_alert(paste(
        "Error: no reference group. Choose the control group under Study",
        "Design on Data Upload."), "danger"))
      return()
    }
    if (!nzchar(rv$settings$col_id %||% "")) {
      output$surv_status_msg <- renderUI(status_alert(paste(
        "Error: map the Mouse/Subject ID column on Data Upload; survival",
        "needs one record per animal."), "danger"))
      return()
    }

    # Track analysis start (no study strings: R20.D12)
    shinyjs::runjs("window.trackAnalysisStart('survival', 'cox', 'study_control');")
    
    rv$surv_running <- TRUE
    output$surv_status_msg <- renderUI({
      status_alert("Running survival analysis...", "info")
    })
    
    # Prepare data — reuse canonical preprocessing if available.
    # R20.D4: every other module wraps this; unwrapped, a date error reached
    # the user as a raw message and left the run flag set.
    pp <- tryCatch(ensure_preprocessed(rv), error = function(e) {
      log_caught_error("survival", "preprocessing", e)
      output$surv_status_msg <- renderUI({
        status_alert(paste0("Preprocessing failed: ", conditionMessage(e)), "danger")
      })
      rv$surv_running <- FALSE
      NULL
    })
    if (is.null(pp)) return()
    if (!is.null(pp$date_fmt_warning)) {
      showNotification(pp$date_fmt_warning, type = "warning", duration = 8)
    }
    data <- pp$data

    # Resolve dose column from UI input directly — do NOT use rv$dose_col which
    # is only set after TG analysis runs. rv$settings$col_dose always reflects the
    # current column mapping regardless of analysis order.
    dose_col_surv <- if (!is.null(rv$settings$col_dose) && nzchar(rv$settings$col_dose) &&
                         rv$settings$col_dose != "None") rv$settings$col_dose else NULL

    # When a dose column is mapped, build a combined "Treatment Dose" group label
    # so each dose level is treated as its own group.
    has_dose_surv <- has_col(data, dose_col_surv) && has_col(data, rv$settings$col_treatment)
    if (has_dose_surv) {
      data[["TreatmentGroup"]] <- paste(data[[rv$settings$col_treatment]],
                                        data[[dose_col_surv]], sep = " ")
    }

    # Handle time column — reuse from preprocessed data
    time_col  <- pp$time_col
    full_data <- data  # preserve full longitudinal data for volume-threshold event derivation

    # One row per animal (R20.D2): keyed by arm + ID (+ cage), at the first
    # flagged event day. See survival_animal_rows().
    id_col_surv   <- rv$settings$col_id
    cage_col_surv <- if (!is.null(rv$settings$col_cage) && nzchar(rv$settings$col_cage) &&
                         rv$settings$col_cage != "None") rv$settings$col_cage else NULL
    if (!has_col(data, cage_col_surv)) cage_col_surv <- NULL
    tx_col_surv <- if (has_dose_surv) "TreatmentGroup" else rv$settings$col_treatment

    if (!has_col(data, id_col_surv) || nrow(data) == 0L) {
      output$surv_status_msg <- renderUI(status_alert(
        if (nrow(data) == 0L) "Error: No data available for analysis"
        else paste0("Error: ID column '", id_col_surv, "' is not in the data."),
        "danger"))
      rv$surv_running <- FALSE
      return()
    }
    data <- survival_animal_rows(
      full_data, time_col, tx_col_surv, id_col_surv, cage_col_surv,
      event_col = if (event_src == "column") rv$settings$col_survival_event)

    # Resolve event/time columns: either from user-supplied censor column or derived
    # by finding the first day each animal's tumour exceeded the volume threshold
    if (event_src == "volume") {
      if (is.null(pp$vol_col) || !pp$vol_col %in% names(full_data)) {
        showNotification(
          "Volume column not mapped. Cannot derive survival events from volume threshold.",
          type = "error"
        )
        rv$surv_running <- FALSE
        return(NULL)
      }
      vol_threshold <- as.numeric(input$surv_vol_threshold %||% 2000)
      ev <- volume_threshold_events(full_data, time_col, tx_col_surv, id_col_surv,
                                    cage_col_surv, vol_col = pp$vol_col,
                                    threshold = vol_threshold)
      m <- match(data$.animal_key, ev$.animal_key)
      data[[".__DerivedTime"]]  <- ev$Time[m]
      data[[".__DerivedEvent"]] <- ev$Event[m]

      censor_col_use <- ".__DerivedEvent"
      time_col_use   <- ".__DerivedTime"
    } else {
      censor_col_use <- rv$settings$col_survival_event
      time_col_use   <- time_col
    }

    # R20.D13 / R20-D.H: with a removal-reason column, deaths from other
    # causes can count as events, at the animal's last measurement.
    death_spec <- removal_reason_spec(rv)
    n_death_events <- NA_integer_
    if (!is.null(death_spec) && isTRUE(input$surv_death_as_event) &&
        has_col(full_data, death_spec$column) && has_col(data, censor_col_use)) {
      reasons <- animal_removal_reasons(full_data, tx_col_surv, id_col_surv,
                                        cage_col_surv, death_spec$column, time_col)
      died <- names(reasons)[!is.na(reasons) &
                               reasons %in% reasons_in(death_spec, "death")]
      ev_now <- suppressWarnings(as.numeric(data[[censor_col_use]]))
      hit <- data$.animal_key %in% died & (is.na(ev_now) | ev_now == 0)
      data[[censor_col_use]] <- ifelse(hit, 1, ev_now)
      n_death_events <- sum(hit)
    }

    # Validate all required columns exist before analysis
    required <- c(Treatment = tx_col_surv,
                  Time      = time_col_use,
                  Event     = censor_col_use)
    missing_cols <- check_required_cols(data, required)
    if (length(missing_cols) > 0) {
      output$surv_status_msg <- renderUI({
        status_alert(paste0("Required column(s) not found: ",
                   paste(names(missing_cols), "=", sQuote(missing_cols), collapse = ", ")), "danger")
      })
      rv$surv_running <- FALSE
      return()
    }

    # ── Bayesian async path ────────────────────────────────────────────────────
    if (isTRUE(input$surv_model_type == "bayes")) {

      bayes_params <- list(
        data           = data,
        time_col       = time_col_use,
        event_col      = censor_col_use,
        tx_col         = tx_col_surv,
        id_col         = id_col_surv,
        cage_col       = cage_col_surv,
        ref_group      = ref_group,
        family         = input$surv_bayes_family    %||% "weibull",
        prior_str      = input$surv_prior_strength  %||% "skeptical",
        n_chains       = as.integer(input$surv_bayes_chains  %||% 4L),
        n_warmup       = as.integer(input$surv_bayes_warmup  %||% 1000L),
        n_iter         = as.integer(input$surv_bayes_iter    %||% 500L),
        backend        = input$surv_bayes_backend %||% "rstan",
        event_src      = event_src,
        death_events_added = n_death_events,
        vol_threshold  = if (event_src == "volume") as.numeric(input$surv_vol_threshold %||% 2000) else NA_real_,
        time_method    = rv$settings$time_method,
        time_col_orig  = time_col,
        ref_date       = rv$settings$reference_date
      )

      # F.7 — route through the job registry so the fit survives tab
      # disconnect / NPM proxy timeout.
      launch_surv_bayes_job(rv$session_id, bayes_params,
                            input, output, session, rv)
      return()
    }

    # ── Frequentist (Cox) path ─────────────────────────────────────────────────
    show_analysis_modal("Fitting survival model…", "Results will appear automatically when complete.")
    on.exit(shiny::removeModal(), add = TRUE)

    result <- safe_analysis(module = "survival", {
      # R20.D14: the ID column is passed (it used to have to be named "ID").
      # R20.6: the declared unit of randomisation decides the test.
      mouseExperiment::survival_statistics(
        df = data,
        time_column = time_col_use,
        censor_column = censor_col_use,
        treatment_column = tx_col_surv,
        cage_column = cage_col_surv,
        id_column = id_col_surv,
        randomisation_unit = if (!is.null(cage_col_surv))
                               rv$settings$randomisation_unit %||% "mouse" else "mouse",
        dose_column = dose_col_surv,
        reference_group = ref_group,
        firth_correction = TRUE,
        verbose = FALSE
      )
    })

    # Store event source for downstream use (KM plot titles, summary text)
    result$event_src_used     <- event_src
    result$vol_threshold_used <- if (event_src == "volume") as.numeric(input$surv_vol_threshold %||% 2000) else NA_real_
    # Deaths counted as events through the removal-reason column (NA: not used).
    result$death_events_added <- n_death_events

    # Store time calculation info for later display
    result$time_method <- rv$settings$time_method
    result$time_column_used <- time_col
    result$reference_date <- rv$settings$reference_date
    
    # Generate Kaplan-Meier plot if analysis succeeded
    if (!isTRUE(result$error) && !is.null(result$survival_data)) {
      tryCatch({
        # Use the survival_data returned by survival_statistics (has Time, Event, Treatment columns)
        surv_df <- result$survival_data

        # Keep a survfit object for the risk table tab
        km_fit <- survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df)
        result$km_fit <- km_fit

        # R17.3: this called mouseExperiment::plot_survival(), which was removed
        # from the backend in ace5b94 ("Remove plot_survival and plot_forest
        # functions"). The dashboard was never updated, so the call threw on
        # every survival run, the tryCatch below logged it, and the fallback
        # ggsurvplot rendered instead. Users saw a plot, so nothing looked wrong
        # -- but the intended plot has never once been drawn. Build it here
        # directly rather than reaching for a function that no longer exists.
        pal <- build_color_vector(input, "surv_color", levels(factor(surv_df$Treatment)))
        km_plot_obj <- survminer::ggsurvplot(
          km_fit,
          data         = surv_df,
          pval         = km_pval_label(surv_df),   # R20.D29
          conf.int     = TRUE,
          risk.table   = TRUE,
          palette      = if (!is.null(pal)) unname(pal) else NULL,
          ggtheme      = ggplot2::theme_classic(),
          title        = "Kaplan-Meier Survival Curves by Treatment",
          legend.title = "Treatment"
        )
        result$km_plot <- km_plot_obj

      }, error = function(e) {
        # Fallback to manual ggsurvplot if plot_survival() fails
        log_caught_error("survival", "plot_survival()", e)
        tryCatch({
          surv_df <- result$survival_data
          km_fit <- survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df)
          result$km_fit <<- km_fit
          km_plot_obj <- survminer::ggsurvplot(
            km_fit,
            data = surv_df,
            risk.table = TRUE,
            risk.table.y.text.col = TRUE,
            conf.int = FALSE,
            xlab = if (event_src == "volume") "Time to Endpoint (Days)" else "Time (Days)",
            ylab = if (event_src == "volume") "Proportion Below Volume Threshold" else "Survival Probability",
            title = if (event_src == "volume") "Time-to-Tumour-Volume-Threshold by Treatment" else "Kaplan-Meier Survival Curves by Treatment"
          )
          result$km_plot <<- km_plot_obj
        }, error = function(e2) {
          log_caught_error("survival", "ggsurvplot fallback", e2)
          NULL
        })
      })
    }
    
    rv$surv_results <- result
    updateTabsetPanel(session, "surv_results_tabs", selected = "Summary")
    rv$surv_running <- FALSE

    # Track analysis completion with Umami
    if (isTRUE(result$error)) {
      shinyjs::runjs(sprintf("window.trackAnalysisComplete('survival', false, %s);",
                             js_string(result$message)))
    } else {
      shinyjs::runjs("window.trackAnalysisComplete('survival', true);")
    }

    # Store treatment group names for the per-group colour pickers in the Plots tab
    rv$surv_plot_groups <- if (!is.null(result$survival_data) && has_col(result$survival_data, "Treatment")) {
      unique(as.character(result$survival_data$Treatment))
    } else NULL
    
    if (isTRUE(result$error)) {
      output$surv_status_msg <- renderUI({
        status_alert(paste("Error:", result$message), "danger")
      })
    } else {
      # R20.D34: the backend's warnings are shown, not dropped.
      output$surv_status_msg <- renderUI(
        completion_status(list(result), "Survival analysis completed successfully!"))
    }
  })
  
  # ── Bayesian tab management ────────────────────────────────────────────────
  bayes_surv_tabs_inserted <- reactiveVal(FALSE)

  observeEvent(input$surv_model_type, {
    tabs_id <- "surv_results_tabs"
    if (isTRUE(input$surv_model_type == "bayes") && !isTRUE(bayes_surv_tabs_inserted())) {
      appendTab(inputId = tabs_id,
        tabPanel("Posterior Check",
          tab_info(
            strong("Posterior predictive check (PPC)"), " — the model generates new datasets from the posterior and overlays their densities on the observed data.",
            tags$ul(
              tags$li(strong("y"), " — observed survival times."),
              tags$li(strong("y_rep"), " — simulated datasets from the posterior predictive distribution.")
            ),
            em("Good model fit: simulated densities closely envelope the observed density.")
          ),
          div(class = "text-right", style = "margin-bottom: 6px;",
            downloadButton(session$ns("download_surv_pp_check"), "Download Plot")
          ),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("surv_pp_check_plot"), height = "450px")
          )
        ),
        select = FALSE
      )
      appendTab(inputId = tabs_id,
        tabPanel("Prior vs Posterior",
          tab_info(
            strong("Prior vs Posterior density overlay"), " — shows how much the data have updated each parameter's prior distribution.",
            em("Parameters whose posterior is similar to the prior are not well identified by the data; consider collecting more data or using an informative prior.")
          ),
          div(class = "text-right", style = "margin-bottom: 6px;",
            downloadButton(session$ns("download_surv_prior_posterior"), "Download Plot")
          ),
          shinycssloaders::withSpinner(
            plotOutput(session$ns("surv_prior_posterior_plot"), height = "400px")
          )
        ),
        select = FALSE
      )
      bayes_surv_tabs_inserted(TRUE)
    } else if (!isTRUE(input$surv_model_type == "bayes") && isTRUE(bayes_surv_tabs_inserted())) {
      removeTab(inputId = tabs_id, target = "Posterior Check")
      removeTab(inputId = tabs_id, target = "Prior vs Posterior")
      bayes_surv_tabs_inserted(FALSE)
    }
  }, ignoreNULL = FALSE)

  # ── Bayesian render outputs ────────────────────────────────────────────────
  output$surv_pp_check_plot <- renderPlot({
    req(rv$surv_results)
    req(isTRUE(rv$surv_results$model_type_used == "bayes"))
    # R20.D18: from the result; the fitted model is no longer kept.
    req(!isTRUE(rv$surv_results$error), !is.null(rv$surv_results$pp_check_plot))
    rv$surv_results$pp_check_plot
  })

  output$surv_mcmc_diag_table <- renderDT({
    req(rv$surv_results)
    req(isTRUE(rv$surv_results$model_type_used == "bayes"))
    req(!isTRUE(rv$surv_results$error))
    diag <- rv$surv_results$mcmc_diagnostics
    req(!is.null(diag))
    datatable(
      format_df_for_display(as.data.frame(diag)),
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE)
    )
  })

  # Round 2 D — Bayesian diagnostics panel + influential observations.
  output$surv_bayes_diagnostics <- renderUI({
    req(rv$surv_results)
    req(isTRUE(rv$surv_results$model_type_used == "bayes"))
    req(!isTRUE(rv$surv_results$error))
    bayes_diagnostics_panel(rv$surv_results)
  })

  output$surv_bayes_influential <- renderDT({
    req(rv$surv_results)
    req(isTRUE(rv$surv_results$model_type_used == "bayes"))
    req(!isTRUE(rv$surv_results$error))
    inf_df <- bayes_influential_obs(rv$surv_results)
    if (nrow(inf_df) == 0L) {
      return(datatable(
        data.frame(Note = sprintf("No observations exceed the Pareto-k threshold (%s).",
                                  format(rv$surv_results$loo_diagnostics$k_threshold %||% 0.7))),
        options = list(dom = "t", paging = FALSE, searching = FALSE), rownames = FALSE
      ))
    }
    datatable(
      inf_df,
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE), rownames = FALSE
    )
  })

  output$surv_mcmc_trace_plot <- renderPlot({
    req(rv$surv_results)
    req(isTRUE(rv$surv_results$model_type_used == "bayes"))
    req(!isTRUE(rv$surv_results$error), !is.null(rv$surv_results$posterior_draws))
    bayesplot::mcmc_trace(rv$surv_results$posterior_draws)
  })

  output$surv_prior_posterior_plot <- renderPlot({
    req(rv$surv_results)
    req(isTRUE(rv$surv_results$model_type_used == "bayes"))
    req(!isTRUE(rv$surv_results$error))
    p <- rv$surv_results$prior_posterior_plot
    req(!is.null(p))
    print(p)
  })

  # Summary output
  output$surv_summary <- renderText({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) {
      return(paste("Error:", rv$surv_results$message))
    }

    summary_text <- capture.output({
      cat("=== SURVIVAL ANALYSIS SUMMARY ===\n\n")
      
      # Show time specification method
      if (rv$surv_results$time_method == "from_date") {
        cat("Time Method: Calculated from Date Column\n")
        cat("Date Column:", rv$settings$col_date_for_day, "\n")
        cat("Reference Date:", rv$settings$reference_date, "\n")
      } else {
        cat("Time Method: Direct Day Column\n")
        cat("Day Column:", rv$settings$col_day, "\n")
      }
      
      if (isTRUE(rv$surv_results$event_src_used == "volume")) {
        cat("Event Source: Derived from tumour volume threshold (>=",
            rv$surv_results$vol_threshold_used, volume_unit_label(rv), ")\n")
      } else {
        cat("Event Column:", rv$settings$col_survival_event, "\n")
      }
      nd <- rv$surv_results$death_events_added
      if (length(nd) == 1L && !is.na(nd)) {
        cat("Deaths from other causes counted as events:", nd,
            "animal(s), from the removal-reason column\n")
      }
      cat("Treatment Column:", rv$settings$col_treatment, "\n")
      cat("Statistical Method:", rv$surv_results$method_used, "\n")
      unit <- rv$surv_results$randomisation_unit %||% "mouse"
      cat("Unit of Randomisation:",
          if (identical(unit, "cage")) "whole cages (p-values from a cage-level permutation log-rank)"
          else "individual animals", "\n")
      
      if (!is.null(rv$surv_results$results)) {
        res <- rv$surv_results$results
        cat("Reference Group:", rv$surv_results$reference_group, "\n\n")
        
        cat("--- Median Survival ---\n")
        for (i in seq_len(nrow(res))) {
          med <- if (is.na(res$Median_Survival[i])) "Not reached" else paste0(res$Median_Survival[i], " days")
          note <- if (!is.na(res$Note[i]) && res$Note[i] != "") paste0(" [", res$Note[i], "]") else ""
          cat(sprintf("  %-20s %s%s\n", res$Group[i], med, note))
        }
      }
      elapsed_str <- format_elapsed(rv$surv_results$elapsed_seconds)
      if (!is.null(elapsed_str))
        cat("\nCOMPUTATION TIME:", elapsed_str, "\n")
    })

    paste(summary_text, collapse = "\n")
  })
  
  # Group order and color combined control for KM/Survival plots
  output$surv_group_order_color_ui <- renderUI({
    req(rv$surv_plot_groups)
    group_order_and_color_ui("surv", rv$surv_plot_groups, ns = session$ns)
  })

  # Apply colorblind-friendly palette presets when selector changes
  observe_palette_preset(input, session, "surv", reactive(rv$surv_plot_groups))

  output$surv_download_plot_btn <- renderUI({
    btn_id <- switch(input$surv_plot_selector %||% "km",
      km     = "download_surv_plot_km",
      forest = "download_surv_forest_plot",
      "download_surv_plot_km"
    )
    downloadButton(session$ns(btn_id), "Download Plot")
  })

  # KM plot builder — shared by renderPlot and download handler
  km_plot_r <- reactive({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) return(NULL)
    req(rv$surv_results$km_fit, rv$surv_results$survival_data)

    surv_df      <- rv$surv_results$survival_data
    ordered_grps <- ordered_surv_grps()
    if (!is.null(ordered_grps) && has_col(surv_df, "Treatment"))
      surv_df$Treatment <- factor(surv_df$Treatment, levels = ordered_grps)
    km_fit <- tryCatch(
      survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df),
      error = function(e) {
        log_caught_error("survival", "KM refit", e)
        rv$surv_results$km_fit
      }
    )
    groups  <- if (!is.null(ordered_grps)) ordered_grps else rv$surv_plot_groups
    pal     <- build_color_vector(input, "surv_color", groups)
    if (!is.null(pal)) pal <- unname(pal)
    if (is.null(pal)) pal <- "hue"

    x_start       <- if (!is.null(input$surv_km_xstart) && !is.na(input$surv_km_xstart))
                       as.numeric(input$surv_km_xstart) else 0
    text_size     <- input$surv_opt_text_size %||% 16
    event_src_km  <- rv$surv_results$event_src_used %||% "column"
    default_title <- if (event_src_km == "volume") "Time-to-Tumour-Volume-Threshold by Treatment"
                     else                          "Kaplan-Meier Survival Curves by Treatment"
    plot_title    <- if (nzchar(trimws(input$surv_plot_title %||% "")))
                       input$surv_plot_title else default_title
    plot_sub      <- input$surv_plot_subtitle %||% ""
    title_hjust   <- switch(input$surv_plot_title_align %||% "center",
                            left = 0, center = 0.5, right = 1, 0.5)
    xlab_custom   <- trimws(input$surv_plot_xaxis_label  %||% "")
    ylab_custom   <- trimws(input$surv_plot_yaxis_label  %||% "")
    leg_label     <- trimws(input$surv_plot_legend_label %||% "")
    legend_pos    <- input$surv_opt_legend_pos %||% "right"

    tryCatch({
      km_obj <- survminer::ggsurvplot(
        km_fit,
        data                  = surv_df,
        palette               = pal,
        risk.table            = isTRUE(input$surv_opt_risk_table),
        risk.table.y.text.col = TRUE,
        conf.int              = isTRUE(input$surv_opt_conf_int),
        xlim                  = c(x_start, NA),
        xlab                  = if (nzchar(xlab_custom)) xlab_custom else
                                  if (event_src_km == "volume") "Time to Endpoint (Days)" else "Time (Days)",
        ylab                  = if (nzchar(ylab_custom)) ylab_custom else
                                  if (event_src_km == "volume") "Proportion Below Volume Threshold" else "Survival Probability",
        legend.title          = if (nzchar(leg_label)) leg_label else "Treatment",
        title                 = plot_title,
        ggtheme               = ggplot2::theme_classic(base_size = text_size)
      )
      km_obj$plot <- km_obj$plot +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = title_hjust))
      if (nzchar(trimws(plot_sub))) {
        km_obj$plot <- km_obj$plot +
          ggplot2::labs(subtitle = plot_sub) +
          ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = title_hjust))
      }
      if (startsWith(legend_pos, "inplot_")) {
        coord <- switch(legend_pos,
          inplot_tl = c(0.02, 0.98), inplot_tr = c(0.98, 0.98),
          inplot_bl = c(0.02, 0.02), inplot_br = c(0.98, 0.02), c(0.5, 0.5))
        just <- switch(legend_pos,
          inplot_tl = c("left","top"),  inplot_tr = c("right","top"),
          inplot_bl = c("left","bottom"), inplot_br = c("right","bottom"),
          c("center","center"))
        km_obj$plot <- km_obj$plot +
          ggplot2::theme(
            legend.position      = coord,
            legend.justification = just,
            legend.background    = ggplot2::element_rect(
              fill = scales::alpha("white", 0.8), colour = NA))
      } else {
        km_obj$plot <- km_obj$plot + ggplot2::theme(legend.position = legend_pos)
      }
      km_obj
    }, error = function(e) rv$surv_results$km_plot)
  })

  output$surv_plot_km <- renderPlot({
    km_obj <- km_plot_r()
    req(!is.null(km_obj))
    print(km_obj)
  }, width = function() input$surv_plot_width %||% 800,
     height = function() input$surv_plot_height %||% 600)
  
  # Cox model summary
  output$surv_cox_summary <- renderText({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) return("Analysis error — see Status.")
    if (isTRUE(rv$surv_results$model_type_used == "bayes")) {
      ps <- rv$surv_results$posterior_summary
      if (is.null(ps)) return("Posterior summary not available.")
      return(paste(capture.output(print(as.data.frame(ps))), collapse = "\n"))
    }
    if (is.null(rv$surv_results$model)) return("Model not available")
    paste(capture.output(summary(rv$surv_results$model)), collapse = "\n")
  })

  # C-index + PH-check summary panel (frequentist Cox path only).
  # mE v0.4.5 E.2 (c_index) and A.1 (ph_test); rendered alongside the Cox table
  # so users can read model quality at a glance.
  output$surv_model_quality_panel <- renderUI({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) return(NULL)
    if (isTRUE(rv$surv_results$model_type_used == "bayes")) return(NULL)

    method <- rv$surv_results$method_used %||% NA_character_
    c_idx  <- rv$surv_results$c_index
    ph     <- rv$surv_results$ph_test

    parts <- list()

    # Concordance / C-index (Cox path only)
    if (!is.null(c_idx) && !is.null(c_idx$concordance)) {
      val <- as.numeric(c_idx$concordance)
      se  <- tryCatch(sqrt(as.numeric(c_idx$var)), error = function(e) NA_real_)
      ci_lo <- if (is.finite(se)) val - 1.96 * se else NA_real_
      ci_hi <- if (is.finite(se)) val + 1.96 * se else NA_real_
      txt <- if (is.finite(se)) {
        sprintf("Concordance (C-index): %.3f (95%% CI %.3f – %.3f). 0.5 = no discrimination, 1.0 = perfect.",
                val, ci_lo, ci_hi)
      } else {
        sprintf("Concordance (C-index): %.3f. 0.5 = no discrimination, 1.0 = perfect.", val)
      }
      parts$c_index <- shiny::div(class = "alert alert-info",
        shiny::icon("scale-balanced"), " ", txt)
    }

    # PH check (cox.zph) global summary
    if (!is.null(ph) && !is.null(ph$table)) {
      tab <- as.data.frame(ph$table)
      if ("GLOBAL" %in% rownames(tab) && "p" %in% colnames(tab)) {
        gp <- as.numeric(tab["GLOBAL", "p"])
        if (is.finite(gp) && gp < 0.05) {
          parts$ph <- shiny::div(class = "alert alert-warning",
            shiny::icon("triangle-exclamation"),
            shiny::strong(" Proportional-hazards assumption may be violated: "),
            sprintf("global cox.zph p = %.4f (< 0.05). ", gp),
            "Inspect the PH Assumption tab for details; reported HRs / CIs may be misleading.")
        } else if (is.finite(gp)) {
          parts$ph <- shiny::div(class = "alert alert-success",
            shiny::icon("check-circle"),
            shiny::strong(" PH assumption OK: "),
            sprintf("global cox.zph p = %.4f (≥ 0.05).", gp))
        }
      }
    }

    # Unit of randomisation (package R20.6)
    cp <- rv$surv_results$cage_permutation
    if (identical(rv$surv_results$randomisation_unit, "cage") && !is.null(cp) && nrow(cp)) {
      parts$cage <- shiny::div(class = "alert alert-info",
        shiny::icon("layer-group"), shiny::strong(" Whole cages were randomised. "),
        "Each p-value comes from a cage-level permutation log-rank, which moves whole ",
        "cages between the arm and the reference; the hazard ratios have no interval, ",
        "because no valid cage-level interval is available. With few cages the ",
        "smallest attainable p-value can exceed 0.05, whatever the effect.",
        shiny::tags$table(class = "table table-condensed",
          style = "width:auto; margin:8px 0 0 0; background:#fff;",
          shiny::tags$thead(shiny::tags$tr(
            shiny::tags$th("Arm"), shiny::tags$th("Cages (arm / reference)"),
            shiny::tags$th("Cage assignments"), shiny::tags$th("p (unadjusted)"),
            shiny::tags$th("Smallest attainable p"))),
          shiny::tags$tbody(lapply(seq_len(nrow(cp)), function(i) shiny::tags$tr(
            shiny::tags$td(cp$Group[i]),
            shiny::tags$td(sprintf("%d / %d", cp$Cages_Treated[i], cp$Cages_Reference[i])),
            shiny::tags$td(paste0(format(cp$Assignments[i], big.mark = ","),
                                  if (isTRUE(cp$Exact[i])) " (all enumerated)" else " (sampled)")),
            shiny::tags$td(format_pval(cp$P_Value[i])),
            shiny::tags$td(signif(cp$Min_Attainable_P[i], 3))))))
      )
    }
    if (!is.null(rv$surv_results$cage_caveat)) {
      parts$cage_caveat <- shiny::div(class = "alert alert-warning",
        shiny::icon("triangle-exclamation"), " ", rv$surv_results$cage_caveat)
    }

    if (length(parts) == 0L) return(NULL)
    do.call(shiny::tagList, unname(parts))
  })

  # PH Assumption tab — Schoenfeld test table (cox.zph)
  output$surv_ph_test_banner <- renderUI({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) return(NULL)
    ph <- rv$surv_results$ph_test
    if (is.null(ph)) {
      return(shiny::div(class = "alert alert-info",
        shiny::icon("info-circle"),
        " The PH check is only computed for the standard Cox path. ",
        sprintf("This analysis used method_used = '%s'; no Schoenfeld test available.",
                rv$surv_results$method_used %||% "unknown")))
    }
    NULL
  })

  # Schoenfeld residuals plot — the `cox.zph` object supports plot()
  # natively (base graphics). We wrap in a guarded block so the panel
  # gracefully shows a message when ph_test isn't available (Firth /
  # log-rank fallback) or when the user is in Bayesian mode.
  output$surv_schoenfeld_plot <- renderPlot({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) return(NULL)
    ph <- rv$surv_results$ph_test
    if (is.null(ph)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Schoenfeld plot not available for this fit path")
      return(invisible(NULL))
    }
    # cox.zph's plot method does the right thing — one panel per
    # covariate plus a global panel; smoothers with 95% CI ribbons.
    tryCatch(plot(ph), error = function(e) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = sprintf("Schoenfeld plot failed: %s",
                          conditionMessage(e)))
    })
  })

  output$surv_ph_test_table <- renderDT({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) {
      return(datatable(data.frame(Message = "Results not available")))
    }
    ph <- rv$surv_results$ph_test
    if (is.null(ph) || is.null(ph$table)) {
      return(datatable(
        data.frame(Note = "cox.zph not available — non-Cox method was used."),
        options = list(dom = "t", paging = FALSE, searching = FALSE), rownames = FALSE
      ))
    }
    tab <- as.data.frame(ph$table)
    tab <- cbind(Covariate = rownames(tab), tab)
    rownames(tab) <- NULL
    tab <- format_df_for_display(tab)
    datatable(tab,
              options = list(scrollX = TRUE, paging = FALSE, searching = FALSE),
              rownames = FALSE)
  })

  # Results table — Cox HR table (frequentist) or treatment effects (Bayesian)
  output$surv_cox_table <- renderDT({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) {
      return(datatable(data.frame(Message = "Results not available")))
    }
    if (isTRUE(rv$surv_results$model_type_used == "bayes")) {
      te <- rv$surv_results$treatment_effects
      if (is.null(te)) return(datatable(data.frame(Message = "Treatment effects not available")))
      return(datatable(format_df_for_display(as.data.frame(te)),
                       options = list(scrollX = TRUE, paging = FALSE, searching = FALSE)))
    }
    if (is.null(rv$surv_results$results)) {
      return(datatable(data.frame(Message = "Results not available")))
    }
    results_df <- format_df_for_display(rv$surv_results$results)
    datatable(results_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Risk table output - numbers at risk per time point per group
  output$surv_risk_table <- renderDT({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$km_fit)) {
      return(datatable(data.frame(Message = "Risk table not available")))
    }
    
    km_sum <- summary(rv$surv_results$km_fit)
    risk_df <- data.frame(
      Time      = km_sum$time,
      Group     = sub("^Treatment=", "", as.character(km_sum$strata)),
      N_Risk    = km_sum$n.risk,
      N_Events  = km_sum$n.event,
      N_Censored = km_sum$n.censor,
      Survival  = round(km_sum$surv, 4)
    )
    # Reorder columns for readability
    risk_df <- risk_df[, c("Group", "Time", "N_Risk", "N_Events", "N_Censored", "Survival")]
    risk_df <- risk_df[order(risk_df$Group, risk_df$Time), ]
    
    datatable(risk_df,
              options = list(scrollX = TRUE, paging = FALSE, searching = FALSE),
              rownames = FALSE)
  })
  
  # build_surv_forest_plot() lives in helpers_plot.R (R20.D19): the download
  # handler in mod_survival_downloads.R could not see it as a local here.
  # Forest plot
  output$surv_forest_plot <- renderPlot({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$results)) return(NULL)
    
    res <- rv$surv_results$results
    req(all(c("Group", "HR", "CI_Lower", "CI_Upper", "P_Value") %in% colnames(res)))
    
    ref_grp <- rv$surv_results$reference_group
    
    # Verify at least one non-reference group with valid HR
    plot_df <- res[!is.na(res$HR), ]
    if (!is.null(ref_grp)) plot_df <- plot_df[plot_df$Group != ref_grp, ]
    req(nrow(plot_df) > 0)

    pt_color <- input$surv_forest_pt_color %||% "steelblue"
    show_hr  <- isTRUE(input$surv_forest_show_hr)
    show_ci  <- isTRUE(input$surv_forest_show_ci)
    show_p   <- isTRUE(input$surv_forest_show_p)

    p <- build_surv_forest_plot(res, ref_grp, pt_color, show_hr, show_ci, show_p)
    req(!is.null(p))
    text_size   <- input$surv_opt_text_size %||% 16
    legend_pos  <- input$surv_opt_legend_pos %||% "right"
    p <- apply_plot_labels(p,
                           input$surv_plot_title        %||% "",
                           input$surv_plot_subtitle     %||% "",
                           input$surv_plot_title_align  %||% "center",
                           text_size,
                           xlab       = input$surv_plot_xaxis_label  %||% "",
                           ylab       = input$surv_plot_yaxis_label  %||% "",
                           legend_lab = input$surv_plot_legend_label %||% "")
    apply_legend_pos(p, legend_pos)
  }, width = function() input$surv_plot_width %||% 800,
     height = function() input$surv_plot_height %||% 600)

  # All 6 download handlers (2 Bayesian + 4 frequentist) live in
  # `mod_survival_downloads.R` (CODE_REVIEW.md C.1/C.3). km_plot_r is
  # passed explicitly so the KM-PNG handler's reactive dependency stays
  # visible at the call site.
  add_surv_download_handlers(input, output, session, rv, km_plot_r)

  # F.7 — register the resume handler for bayes_surv jobs.
  register_surv_resume(input, output, session, rv)

  })
}
