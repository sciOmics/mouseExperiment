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

            # The reference is the study's control group (R20-D.H).
            uiOutput(ns("dr_reference_note")),

            # R20.D47: one agent's dose series.
            uiOutput(ns("dr_series_ui")),

            selectInput(ns("dr_comparison_type"),
                       tagList("Pairwise Comparison Type", help_tip(paste(
                         "vs. Reference Group compares each dose with the control.",
                         "All Pairwise compares every pair of doses (Tukey HSD)."))),
                       choices = c(
                         "vs. Reference Group" = "reference",
                         "All Pairwise (Tukey HSD)" = "all"
                       ),
                       selected = "reference"),

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
                      tags$li(strong("Linear regression"), " \u2014 tests whether there is a significant linear trend in tumor volume with increasing dose."),
                      tags$li(strong("One-way ANOVA"), " \u2014 tests whether any dose group differs from the others, followed by Tukey HSD post-hoc comparisons if significant."),
                      tags$li(strong("Non-linear dose-response model (4PL)"), " \u2014 fits a four-parameter log-logistic (Hill) model to quantify ED50 and curve steepness.")
                    ),
                    p(strong("One agent at a time. "),
                      "A dose-response curve describes one agent at several doses, against the study's control. ",
                      "Choose the agent's arms under ", strong("Dose series"), " in the sidebar; the dashboard ",
                      "suggests the largest set of single-agent arms that share a name and differ in dose. Arms ",
                      "that share a dose, or other agents and combinations, cannot join the same curve, and the ",
                      "analysis stops rather than pool them. When one treatment name is used at every dose, the ",
                      "rows with dose 0 are the control."),
                    p(tags$small(class = "text-muted",
                      "A Dose column must be mapped in Data Upload."
                    ))
                  )
                ),

                # ── Evaluable days (R20-K) ───────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Evaluable days")),
                  div(class = "panel-body",
                    p(evaluable_rule_text()),
                    p("The analysis uses one day on which every dose group is evaluable: by default the ",
                      "last such day, with the animals measured on it. The day selector lists only ",
                      "evaluable days. Before, the default was each animal's own last measurement, so ",
                      "controls removed at the volume limit contributed capped, earlier volumes and the ",
                      "EC50 was biased. The results name the day and list the days and dose groups ",
                      "the rule excluded."),
                    p(strong("TGI by dose"), " is read at the same day from the model that the Drug ",
                      "Synergy and Toxicity tabs use, so the three tabs agree.")
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
                      tags$dt("R\u00b2"),
                      tags$dd("Proportion of variance in volume explained by dose. Low R\u00b2 with a significant slope can occur when individual variability is high relative to group differences."),
                      tags$dt("p-value"),
                      tags$dd("Significance of the slope. A non-significant p-value does not exclude a non-linear relationship \u2014 inspect the 4PL model output and plots as well.")
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
                      "Volume = c + (d \u2212 c) / (1 + (Dose / ED50)^b)"
                    )),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("b"), tags$dd("Hill coefficient: steepness of the dose-response curve. Larger absolute values indicate a sharper transition."),
                      tags$dt("c"), tags$dd("Lower asymptote (response at saturating inhibitory dose)."),
                      tags$dt("d"), tags$dd("Upper asymptote (response in the absence of drug)."),
                      tags$dt("ED50"), tags$dd("Dose producing 50 % of the maximal effect. The primary pharmacodynamic summary statistic.")
                    ),
                    p("The EC50 reported is the dose giving half of the fitted response range, with a 95 % interval ",
                      "computed on log dose, so it stays positive. The lower asymptote cannot go below zero, and the ",
                      "five-parameter (asymmetric) curve is considered only with at least six dose levels. The summary ",
                      "flags an EC50 outside the tested doses, and a curve with as many parameters as dose levels, which ",
                      "fits the dose means exactly and so says nothing about its own shape; read either as descriptive."),
                    div(class = "alert alert-info", style = "margin-top:10px;",
                      icon("circle-info"), " ",
                      strong("AIC / BIC note: "),
                      "AIC and BIC from the linear model (", tags$code("stats::lm"), ") and the non-linear model (",
                      tags$code("drc::drm"), ") are ", strong("not directly comparable"), " for model selection \u2014 ",
                      "the two packages use different likelihood parameterisations. ",
                      "Prefer visual inspection of the fitted curve and residuals when choosing between linear and non-linear models."
                    )
                  )
                ),

                # ── Trend tests ──────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-danger", "Trend Tests"))),
                  div(class = "panel-body",
                    p("A dose-response study usually asks whether response is ", em("monotone in dose"),
                      " \u2014 not merely whether the groups differ somewhere. Three tests address that directly:"),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Jonckheere\u2013Terpstra"),
                      tags$dd("The test that actually respects dose ordering. It is sensitive to a consistent ",
                              "monotone trend and, unlike ANOVA, gains power from the ordering rather than ",
                              "discarding it. This is the primary trend result. It is two-sided; the summary ",
                              "also says whether volumes fall or rise with dose. (It was one-sided in the ",
                              "direction the data suggested, which doubled its false-positive rate.)"),
                      tags$dt("Linear trend"),
                      tags$dd("p-value for a linear term in dose. Significant means a straight-line component ",
                              "is present; it says nothing about whether the relationship is ", em("only"), " linear."),
                      tags$dt("Quadratic trend"),
                      tags$dd("Detects curvature \u2014 a plateau or a reversal at high dose. A significant quadratic ",
                              "term alongside a significant linear one usually means the response saturates, which ",
                              "matters when extrapolating to untested doses.")
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
            uiOutput(ns("dr_evaluability_note")),
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
              DT::DTOutput(ns("dr_summary_table"))
            ),
            h4("TGI by dose"),
            p(class = "text-muted",
              "At the analysis day, from the endpoint model (the TGI the Drug Synergy and Toxicity tabs ",
              "report), with 95% intervals from draws of the model."),
            DT::DTOutput(ns("dr_tgi_table"))
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
                  sliderInput(ns("dr_opt_jitter_width"), "Jitter width",
                              min = 0, max = 1, value = 0, step = 0.05),
                  checkboxInput(ns("dr_scatter_show_ci"),
                                "Show trend line confidence interval",
                                value = TRUE),
                  radioButtons(ns("dr_trend_method"), "Trend line fit",
                    choices = c("Linear" = "linear", "Non-linear (LOESS)" = "loess"),
                    selected = "linear", inline = TRUE
                  ),
                  selectInput(ns("dr_opt_center_bar"), "Centre overlay (per dose)",
                    choices = c("None" = "none",
                                "Mean"   = "mean",
                                "Median" = "median"),
                    selected = "none"
                  ),
                  selectInput(ns("dr_opt_error_bar"), "Error bar overlay",
                    choices = c("None" = "none",
                                "SEM" = "SEM",
                                "SD"  = "SD",
                                "95% CI" = "CI"),
                    selected = "none"
                  ),
                  hr(),
                  h6(style = "font-size: 1.15rem; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                     icon("sliders-h"), " Options"),
                  sliderInput(ns("dr_opt_text_size"), "Text size (pt)",
                              min = 10, max = 24, value = 16, step = 1),
                  selectInput(ns("dr_opt_legend_pos"), "Legend position",
                    choices = c("None" = "none",
                                "Top" = "top", "Bottom" = "bottom",
                                "Left" = "left", "Right" = "right",
                                "In plot - top left"  = "inplot_tl",
                                "In plot - top right" = "inplot_tr",
                                "In plot - bottom left" = "inplot_bl",
                                "In plot - bottom right" = "inplot_br"),
                    selected = "right"
                  ),
                  hr(),
                  h6(style = "font-size: 1.15rem; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                     icon("heading"), " Title, Labels & Subtitle"),
                  textInput(ns("dr_plot_title"), "Title",
                            placeholder = "Leave blank for default"),
                  selectInput(ns("dr_plot_title_align"), "Title alignment",
                    choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                    selected = "center"),
                  textInput(ns("dr_plot_subtitle"), "Subtitle",
                            placeholder = "Optional subtitle"),
                  textInput(ns("dr_plot_xaxis_label"), "X-axis label",
                            placeholder = "Leave blank for default"),
                  textInput(ns("dr_plot_yaxis_label"), "Y-axis label",
                            placeholder = "Leave blank for default"),
                  textInput(ns("dr_plot_legend_label"), "Legend label",
                            placeholder = "Leave blank for default"),
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
              DT::DTOutput(ns("dr_tukey_table"))
            )
          ),

          # CODE_REVIEW.md DIAGNOSTICS gap (5) \u2014 frequentist DR had no
          # goodness-of-fit. v0.4.8 surfaces drc::modelFit + residuals plot.
          tabPanel("Diagnostics",
            tab_info(
              strong("Model-fit diagnostics. "),
              "Hill / Emax lack-of-fit and residuals."
            ),
            em("Only available when the drc Hill / Emax fit converged; ",
               "in pure-linear mode this section will be empty."),
            fluidRow(
              column(6,
                h4("Lack-of-fit test"),
                plot_help(
                    "drc::modelFit compares the Hill fit against a one-mean-",
                    "per-dose smoother via an F-test. ",
                    strong("Large p (\u2265 0.05): "),
                    "the Hill model is as good as the smoother \u2014 ",
                    "functional form is adequate. ",
                    strong("Small p (< 0.05): "),
                    "the smoother fits significantly better \u2014 the Hill ",
                    "equation may be the wrong shape; inspect the residuals ",
                    "plot for clues.",
                    summary = "How to read this table"),
                shinycssloaders::withSpinner(
                  DT::DTOutput(ns("dr_lack_of_fit_table"))
                )
              ),
              column(6,
                h4("Residuals vs Dose"),
                plot_help(
                    "Points should scatter around the red zero line with no ",
                    "trend. ",
                    strong("Curved blue smoother: "),
                    "the Hill function misses curvature at some dose range ",
                    "\u2014 try a different model. ",
                    strong("Fan shape (variance grows with dose): "),
                    "heteroscedasticity; consider log-transforming the ",
                    "volume axis or weighting."),
                shinycssloaders::withSpinner(
                  plotOutput(ns("dr_residuals_plot"), height = "400px")
                )
              )
            )
          ),
          # R20-D.Y: choose which results go into one ZIP.
          tabPanel(
            title = tags$span(icon("file-zipper"), " Download Results"),
            value = "__download_results__",
            wellPanel(download_results_ui(ns, "dr"))
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

  # D.7: this was an unconditional `rv$dr_results <- NULL`. `rv` is created empty
  # in app.R and each module mounts once, so the assignment never cleared anything
  # in production -- but it wiped the slot at module init, which made these five
  # modules impossible to drive through testServer() with a pre-populated result.
  # That is why I.2 stayed open across four rounds: every attempt to assert render
  # behaviour got the "run an analysis first" placeholder back. Removed; the slot
  # reads NULL until an analysis writes it either way.

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

  # The study's control group, in this tab's treatment labels.
  dr_ref <- reactive({
    req(rv$processed_data, rv$tx_col)
    study_reference(rv, rv$processed_data[[rv$tx_col]])
  })

  output$dr_reference_note <- renderUI({
    si <- tryCatch(dr_series_info(), error = function(e) NULL)
    if (!is.null(si) && is.null(si$control)) {
      # One treatment name at every dose: its dose-0 rows are the control.
      return(tags$p(style = "font-size:12px;", "Control: ", strong("the rows with dose 0"),
        ", because the study's control group spans several doses."))
    }
    reference_note(rv, tryCatch(dr_ref(), error = function(e) NULL))
  })

  # R20.D47 / package R20.16: the arms that can form one agent's dose series,
  # a suggested series, and the control. The study's control arm is the
  # control when it is a dose-0 arm; when one treatment name spans every dose,
  # the rows with dose 0 are.
  dr_series_info <- reactive({
    req(rv$processed_data, rv$dose_col)
    pp <- ensure_preprocessed(rv)
    d  <- pp$data
    tx <- pp$tx_col
    dose <- suppressWarnings(as.numeric(d[[rv$dose_col]]))
    ctrl <- study_reference(rv, d[[tx]])
    ctrl_rows <- !is.null(ctrl) & !is.na(d[[tx]]) & d[[tx]] %in% ctrl
    ctrl_is_arm <- any(ctrl_rows) && all(is.na(dose[ctrl_rows]) | dose[ctrl_rows] == 0)
    arms <- treatment_values(d, tx)
    choices <- if (ctrl_is_arm) setdiff(arms, ctrl) else arms
    list(data = d, tx_col = tx, time_col = pp$time_col, vol_col = pp$vol_col,
         control = if (ctrl_is_arm) ctrl, choices = choices,
         suggested = guess_dose_series(d, tx, rv$dose_col, choices))
  })

  output$dr_series_ui <- renderUI({
    si <- tryCatch(dr_series_info(), error = function(e) NULL)
    if (is.null(si)) return(NULL)
    cur <- isolate(input$dr_treatments)
    sel <- if (length(cur) && all(cur %in% si$choices)) cur else si$suggested
    tagList(
      selectizeInput(session$ns("dr_treatments"),
        tagList("Dose series", help_tip("The arms of one agent at its different doses. Other agents, combinations and arms that share a dose cannot join one curve.")),
        choices = si$choices, selected = sel, multiple = TRUE,
        options = list(plugins = list("remove_button"),
                       placeholder = "Choose one agent's arms")),
      if (!length(sel)) tags$p(tags$small(class = "text-warning",
        "Choose the arms of one agent; they could not be guessed from the names."))
    )
  })

  # The rows the analysis will use: the chosen series and its control.
  dr_series_rows <- function(si, series) {
    d <- si$data
    dose <- suppressWarnings(as.numeric(d[[rv$dose_col]]))
    ctrl_rows <- if (!is.null(si$control)) !is.na(d[[si$tx_col]]) & d[[si$tx_col]] == si$control
                 else is.finite(dose) & dose == 0
    keep <- ctrl_rows | (!is.na(d[[si$tx_col]]) & d[[si$tx_col]] %in% series & is.finite(dose))
    d <- d[keep, , drop = FALSE]
    d[[".dose_group"]] <- ifelse(ctrl_rows[keep], 0, dose[keep])
    d
  }

  # Day selector (R20-K, R20.D48): the last evaluable day (default) or another
  # evaluable day of the chosen series. "All dates" meant each animal's own
  # last observation.
  output$dr_time_point_ui <- renderUI({
    req(rv$processed_data, rv$time_col, rv$dose_col)
    days <- tryCatch({
      si <- dr_series_info()
      d  <- dr_series_rows(si, input$dr_treatments %||% si$suggested)
      mouseExperiment::evaluable_days(
        d, treatment_column = ".dose_group", id_column = rv$settings$col_id,
        time_column = si$time_col,
        cage_column = if (has_col(d, rv$settings$col_cage)) rv$settings$col_cage,
        volume_column = si$vol_col)$days
    }, error = function(e) NULL)
    day_choices <- c("Last evaluable day" = "default")
    if (length(days)) {
      day_choices <- c(day_choices, setNames(as.character(days), paste("Day", days)))
    }
    selectInput(session$ns("dr_time_point"),
                tagList("Time Point for Analysis", help_tip(paste(
                  "Only days on which every dose of the series is evaluable are offered:",
                  "at least 50% of a group's enrolled animals, and at least 3, still on",
                  "study. The default is the last of them."))),
                choices  = day_choices,
                selected = "default")
  })

  observeEvent(input$run_dr_analysis, {
    req(rv$processed_data, dr_dose_available())

    output$dr_status_msg <- renderUI({
      status_alert("Running dose-response analysis...", "info")
    })

    # Track analysis start with Umami
    shinyjs::runjs("window.trackAnalysisStart('dose_response', 'linear_regression+anova', 'control');")

    # R20.D29: ensure_preprocessed() applies the current subset filter;
    # rv$processed_data carried the filter as it was when it was cached.
    pp <- tryCatch(ensure_preprocessed(rv), error = function(e) {
      output$dr_status_msg <- renderUI(status_alert(
        paste0("Preprocessing failed: ", conditionMessage(e)), "danger"))
      NULL
    })
    if (is.null(pp)) return()
    data     <- pp$data
    dose_col <- rv$dose_col
    vol_col  <- pp$vol_col
    time_col <- pp$time_col
    tx_col   <- pp$tx_col
    id_col   <- rv$id_col
    ctrl     <- study_reference(rv, data[[tx_col]])

    # Coerce dose column to numeric
    data[[dose_col]] <- suppressWarnings(as.numeric(data[[dose_col]]))
    prep_steps <- c(pp$steps, list(record_code(sprintf(   # R20-AC: for the R script
      "{d}[[%s]] <- suppressWarnings(as.numeric({d}[[%s]]))",
      script_value(dose_col), script_value(dose_col)))))
    n_na_dose <- sum(is.na(data[[dose_col]]))
    if (n_na_dose == nrow(data)) {
      output$dr_status_msg <- renderUI({
        status_alert(paste("Dose column", dose_col, "contains no numeric values. Please check your data."), "danger")
      })
      return()
    }

    # R20.D47: one agent's arms, and its control: the study's control arm when
    # it is a dose-0 arm, otherwise the rows with dose 0.
    si <- tryCatch(dr_series_info(), error = function(e) NULL)
    series <- input$dr_treatments %||% character(0)
    if (is.null(si) || !length(series)) {
      output$dr_status_msg <- renderUI(status_alert(paste(
        "Choose the arms of one agent's dose series under 'Dose series' in the",
        "sidebar."), "warning"))
      return()
    }

    time_point_val <- if (!is.null(input$dr_time_point) &&
                           !input$dr_time_point %in% c("all", "default") &&
                           nzchar(input$dr_time_point)) {
      as.numeric(input$dr_time_point)
    } else {
      NULL
    }

    # ── Frequentist path ───────────────────────────────────────────────────────
    show_analysis_modal("Fitting dose-response model\u2026", "Results will appear automatically when complete.")
    on.exit(shiny::removeModal(), add = TRUE)

    dr_args <- list(
      df                 = data,
      dose_column        = dose_col,
      treatment_column   = tx_col,
      volume_column      = vol_col,
      time_column        = time_col,
      id_column          = id_col,
      cage_column        = rv$cage_col,  # part of the animal key
      endpoint_day       = time_point_val,
      reference_group    = si$control,
      treatments         = series,
      verbose            = FALSE,  # the report went to the server log
      volume_units       = backend_volume_units(rv)   # backend 0.33.0 zero floor
    )
    result <- safe_analysis(module = "dose_response",
                            do.call(mouseExperiment::dose_response_statistics, dr_args))
    # R20-AA: the call exactly as made, for the R script in Download Results.
    result <- with_provenance(result, record_call("dose_response_statistics", dr_args),
                              prep_steps, rv$source_file)
    # The plots and tables use the control of the displayed run (R20.D48).
    if (is.list(result)) result$control_group_used <- si$control %||% "dose 0"
    # The backend names its argument; here the choice is the sidebar selector.
    if (isTRUE(result$error) && is.character(result$message)) {
      result$message <- gsub("with `treatments =`", "under 'Dose series' in the sidebar",
                             result$message, fixed = TRUE)
    }

    # R17.3: this used to call mouseExperiment::prepare_dose_data() and
    # ::generate_summary_statistics() to re-derive the plotting frame. Neither is
    # exported -- both are internal helpers -- so the call threw on every run, the
    # tryCatch swallowed it, and the scatter, box and trend plots all rendered
    # blank behind a silent req(ad). The backend now returns the frame it
    # actually analysed, so read that instead of re-deriving it.
    if (!isTRUE(result$error)) {
      result$summary_stats <- result$summary_table
      if (is.null(result$analysis_data)) {
        result$dr_prep_error <-
          "The dose-response result has no analysis data; plots unavailable."
      }
    }

    rv$dr_results <- result
    updateTabsetPanel(session, "dr_results_tabs", selected = "Summary")

    # Track analysis completion with Umami
    if (isTRUE(result$error)) {
      shinyjs::runjs("window.trackAnalysisComplete('dose_response', false, 'error');")  # R20.D12/D22
    } else {
      shinyjs::runjs("window.trackAnalysisComplete('dose_response', true);")
    }

    if (isTRUE(result$error)) {
      output$dr_status_msg <- renderUI({
        status_alert(paste0("Error: ", result$message), "danger")
      })
    } else {
      # R20.D34: the backend's warnings are shown, not dropped.
      output$dr_status_msg <- renderUI(completion_status(
        list(result), "Dose-response analysis completed successfully!"))
    }
  })

  # ---- Summary ----
  output$dr_summary <- renderText({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(paste("Error:", rv$dr_results$message))

    res   <- rv$dr_results
    stats <- res$statistics

    ser <- res$series
    jt  <- res$trend_test$jonckheere_test
    lines <- c(
      "=== DOSE-RESPONSE ANALYSIS SUMMARY ===",
      "",
      if (!is.null(ser)) c(
        sprintf("Dose series: %s", paste(ser$treatments, collapse = ", ")),
        sprintf("Control:     %s (dose 0)", paste(ser$control, collapse = ", ")),
        sprintf("Day:         %s", format(res$endpoint_day)),
        ""),
      "-- Trend (Jonckheere-Terpstra, two-sided) --",
      if (!is.null(jt)) c(
        sprintf("  p-value: %.4g  (volumes %s with dose)", jt$p.value,
                jt$direction %||% "change"))
      else "  Not run (needs at least three dose levels).",
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
      ci <- stats$ec50_ci
      lines <- c(lines,
        "-- Non-linear Dose-Response Model --",
        sprintf("  EC50:         %.4g  (95%% interval %s, on log dose)", stats$ec50,
                if (!is.null(ci)) sprintf("%.4g to %.4g", ci[["lower"]], ci[["upper"]])
                else "not available"),
        if (nzchar(stats$ec50_note %||% "")) paste("  Note:        ", stats$ec50_note),
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

    elapsed_str <- format_elapsed(rv$dr_results$elapsed_seconds)
    if (!is.null(elapsed_str)) {
      lines <- c(lines, paste("COMPUTATION TIME:", elapsed_str))
    }

    paste(lines, collapse = "\n")
  })

  # ---- Dose-Response Table ----
  output$dr_summary_table <- DT::renderDT({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(dt_table(data.frame(Message = "Results not available")))
    tbl <- as.data.frame(rv$dr_results$summary_table)
    tbl <- format_df_for_display(tbl)
    DT::datatable(tbl, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
  })

  output$dr_download_plot_btn <- renderUI({
    btn_id <- switch(input$dr_plot_selector %||% "scatter",
      scatter = "download_dr_scatter",
      tgi     = "download_dr_tgi",
      "download_dr_scatter"
    )
    downloadButton(session$ns(btn_id), "Download Plot")
  })

  # Helper: build per-dose summary df (center + half-width) for overlay layers.
  # Returns NULL when no overlay was requested. Aggregates `y_col` over `x_col`.
  dr_dose_summary <- function(df, x_col, y_col, center, err) {
    if (center == "none" && err == "none") return(NULL)
    doses <- if (is.factor(df[[x_col]])) levels(df[[x_col]]) else unique(df[[x_col]])
    do.call(rbind, lapply(doses, function(d) {
      y <- df[[y_col]][df[[x_col]] == d]
      y <- y[is.finite(y)]
      n <- length(y)
      if (n == 0L) return(NULL)
      ctr <- switch(center,
        mean   = mean(y),
        median = stats::median(y),
        mean(y))
      half <- switch(err,
        SEM = stats::sd(y) / sqrt(n),
        SD  = stats::sd(y),
        CI  = stats::qt(0.975, df = max(n - 1L, 1L)) * stats::sd(y) / sqrt(n),
        0)
      data.frame(x = d, ctr = ctr, lo = ctr - half, hi = ctr + half,
                 stringsAsFactors = FALSE)
    }))
  }

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
    jit_w    <- input$dr_opt_jitter_width    %||% 0
    show_ci  <- isTRUE(input$dr_scatter_show_ci)
    trend    <- input$dr_trend_method %||% "linear"
    text_sz  <- input$dr_opt_text_size %||% 16
    center   <- input$dr_opt_center_bar %||% "none"
    err_bar  <- input$dr_opt_error_bar  %||% "none"
    smooth_args <- if (trend == "loess") {
      list(method = "loess", formula = y ~ x, color = ln_col, se = show_ci)
    } else {
      list(method = "lm", formula = y ~ x, color = ln_col, se = show_ci)
    }
    p <- ggplot2::ggplot(ad, ggplot2::aes(x = .data[[dose_col]], y = .data[[vol_col]])) +
      ggplot2::geom_jitter(alpha = 0.7, color = pt_col, size = pt_sz,
                           width = jit_w, height = 0) +
      do.call(ggplot2::geom_smooth, smooth_args) +
      ggplot2::labs(title = "Tumor Volume vs Dose",
                   x = "Dose", y = "Tumor Volume") +
      ggplot2::theme_classic()
    summ <- dr_dose_summary(ad, dose_col, vol_col, center, err_bar)
    if (!is.null(summ)) {
      if (err_bar != "none")
        p <- p + ggplot2::geom_errorbar(data = summ,
                ggplot2::aes(x = x, ymin = lo, ymax = hi),
                inherit.aes = FALSE, colour = "black",
                width = 0.05 * diff(range(ad[[dose_col]], na.rm = TRUE)),
                linewidth = 0.9)
      if (center != "none")
        p <- p + ggplot2::geom_crossbar(data = summ,
                ggplot2::aes(x = x, y = ctr, ymin = ctr, ymax = ctr),
                inherit.aes = FALSE, colour = "black",
                width = 0.07 * diff(range(ad[[dose_col]], na.rm = TRUE)),
                linewidth = 0.9)
    }
    p <- apply_plot_labels(p, input$dr_plot_title %||% "", input$dr_plot_subtitle %||% "",
                      input$dr_plot_title_align %||% "center", text_sz,
                      input$dr_plot_xaxis_label %||% "",
                      input$dr_plot_yaxis_label %||% "",
                      input$dr_plot_legend_label %||% "")
    p <- apply_legend_pos(p, input$dr_opt_legend_pos %||% "right")
    current_dr_scatter(p)
    p
  }, width = function() input$dr_plot_width %||% 800,
     height = function() input$dr_plot_height %||% 600)

  output$dr_plot_tgi <- renderPlot({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(NULL)
    # R20.D48: plot the TGI the backend reported -- from the endpoint model at
    # the analysis day, the same TGI as the Synergy and Toxicity tabs. It was
    # recomputed here per animal against the live Reference dropdown and
    # contradicted the Toxicity tab.
    tt <- rv$dr_results$tgi_table
    if (is.null(tt) || !nrow(tt)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "TGI needs a control (dose 0, or the study's control group)")
      return(invisible(NULL))
    }
    ln_col  <- input$dr_scatter_line_color %||% "#2166AC"
    pt_sz   <- input$dr_scatter_pt_size     %||% 2
    text_sz <- input$dr_opt_text_size %||% 16
    p <- ggplot2::ggplot(tt, ggplot2::aes(x = .data[["Dose"]], y = .data[["TGI"]])) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_line(colour = ln_col) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[["TGI_Lower"]], ymax = .data[["TGI_Upper"]]),
                             width = 0.03 * diff(range(tt$Dose)), colour = "black") +
      ggplot2::geom_point(size = pt_sz + 1, colour = ln_col) +
      ggplot2::labs(title = "Tumor Growth Inhibition vs Dose",
                    subtitle = sprintf("Day %s, from the endpoint model; bars are 95%% intervals",
                                       format(rv$dr_results$endpoint_day)),
                    x = "Dose", y = "% Tumor Growth Inhibition") +
      ggplot2::theme_classic()
    p <- apply_plot_labels(p, input$dr_plot_title %||% "", input$dr_plot_subtitle %||% "",
                      input$dr_plot_title_align %||% "center", text_sz,
                      input$dr_plot_xaxis_label %||% "",
                      input$dr_plot_yaxis_label %||% "",
                      input$dr_plot_legend_label %||% "")
    p <- apply_legend_pos(p, input$dr_opt_legend_pos %||% "right")
    current_dr_tgi(p)
    p
  }, width = function() input$dr_plot_width %||% 800,
     height = function() input$dr_plot_height %||% 600)

  output$dr_tgi_table <- DT::renderDT({
    req(rv$dr_results)
    tt <- rv$dr_results$tgi_table
    if (isTRUE(rv$dr_results$error) || is.null(tt)) return(NULL)
    DT::datatable(format_df_for_display(tt), rownames = FALSE,
              options = list(paging = FALSE, searching = FALSE, scrollX = TRUE))
  })

  output$dr_evaluability_note <- renderUI({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(NULL)
    evaluability_note(rv$dr_results$evaluability, day = rv$dr_results$endpoint_day,
                      arms_label = "all dose groups")
  })

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
      jt <- rv$dr_results$trend_test$jonckheere_test
      if (!is.null(jt)) {
        cat("\n=== JONCKHEERE-TERPSTRA TREND TEST (two-sided) ===\n")
        cat(sprintf("  JT = %s,  p-value: %.4g;  volumes %s with dose\n",
                    format(unname(jt$statistic)), jt$p.value, jt$direction %||% "change"))
      }
      if (!is.null(stats$ec50)) {
        cat("\n=== EC50 (ED50 of the fitted curve) ===\n")
        ci <- stats$ec50_ci
        cat(sprintf("  EC50: %.4g;  95%% interval on log dose: %s\n", stats$ec50,
                    if (!is.null(ci)) sprintf("%.4g to %.4g", ci[["lower"]], ci[["upper"]])
                    else "not available"))
        cat(sprintf("  Curve: %s, %s parameters for %s dose levels\n",
                    stats$dr_model_type %||% "?", stats$n_curve_parameters %||% "?",
                    stats$n_dose_levels %||% "?"))
        if (nzchar(stats$ec50_note %||% "")) cat("  Note:", stats$ec50_note, "\n")
      }
      if (!is.null(stats$growth_dose_p_value)) {
        cat("\n=== GROWTH RATE vs DOSE ===\n")
        cat(sprintf("  p-value: %.4g,  R\u00b2: %.4f\n",
                    stats$growth_dose_p_value, stats$growth_dose_r_squared))
        n_out <- stats$growth_rate_animals_left_out %||% 0
        if (n_out > 0) cat(sprintf(paste0("  %d animal(s) left out: fewer than three days ",
                                          "with a positive measured volume.\n"), n_out))
      }
    })
    paste(lines, collapse = "\n")
  })

  output$dr_tukey_table <- DT::renderDT({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(dt_table(data.frame(Message = "No results")))
    tukey <- rv$dr_results$statistics$tukey_results
    if (is.null(tukey)) return(dt_table(data.frame(
      Message = "Tukey HSD not computed (ANOVA p \u2265 0.05 or not available)")))
    tbl <- as.data.frame(tukey[[1]])

    comp_type <- input$dr_comparison_type %||% "reference"
    if (comp_type == "reference") {
      ad  <- rv$dr_results$analysis_data
      if (!is.null(ad)) {
        dose_col <- rv$dose_col %||% "Dose"
        tx_col   <- rv$tx_col   %||% "Treatment"
        # The control is dose 0 by construction (package R20.16).
        ref_doses <- rv$dr_results$control_dose %||% "0"
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
    DT::datatable(tbl, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })

  # CODE_REVIEW.md DIAGNOSTICS gap (5) — frequentist DR Diagnostics tab.
  output$dr_lack_of_fit_table <- DT::renderDT({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(dt_table(data.frame(Message = "No results")))
    lof <- rv$dr_results$statistics$dr_lack_of_fit
    if (is.null(lof) || nrow(lof) == 0L)
      return(dt_table(data.frame(
        Message = "Lack-of-fit not available (Hill / Emax fit did not converge)")))
    dt_table(format_df_for_display(as.data.frame(lof)),
              options = list(dom = "t", paging = FALSE, searching = FALSE))
  })

  output$dr_residuals_plot <- renderPlot({
    req(rv$dr_results)
    if (isTRUE(rv$dr_results$error)) return(NULL)
    p <- rv$dr_results$statistics$dr_residuals_plot
    if (is.null(p)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Residuals plot not available (Hill / Emax fit did not converge)")
      return(invisible(NULL))
    }
    p
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
      out <- dr_stats_table(rv$dr_results)
      write.csv(out, file, row.names = FALSE)
    }
  )


  # "Download Results" (CODE_REVIEW.md R20-D.Y).
  add_download_results(
    input, output, session, prefix = "dr", module = "dose_response",
    items = reactive(dr_dl_items(rv$dr_results,
                                 scatter_plot = function() current_dr_scatter(),
                                 tgi_plot     = function() current_dr_tgi())),
    readme_header = function() dl_readme_header("Dose Response", rv$dr_results, function(res) list(
      dose_series = paste(input$dr_treatments, collapse = ", "),
      evaluation_day = res$endpoint_day, comparisons = input$dr_comparison_type)),
    plot_size = function() list(width  = input$dr_plot_width  %||% 800,
                                height = input$dr_plot_height %||% 600,
                                dpi    = input$dr_plot_dpi    %||% 100))
  })
}

#' The statistics table as downloaded: trend tests, EC50 and the curve
#' @keywords internal
dr_stats_table <- function(res) {
  stats <- res$statistics
  out <- data.frame(
    Metric = c("Linear Slope", "Linear R2", "Linear P-value",
               "ANOVA P-value",
               "Jonckheere-Terpstra P-value (two-sided)",
               "EC50", "EC50 95% Lower (log dose)", "EC50 95% Upper (log dose)",
               "EC50 Within Tested Doses",
               "Hill Slope", "Lower Limit", "Upper Limit",
               "Linear AIC", "Nonlinear AIC",
               "Growth Rate P-value", "Growth Rate R2"),
    Value  = c(
      stats$linear_slope %||% NA,
      stats$linear_r_squared %||% NA,
      stats$linear_p_value %||% NA,
      stats$anova_p_value %||% NA,
      res$trend_test$jonckheere_test$p.value %||% NA,
      stats$ec50 %||% NA,
      if (!is.null(stats$ec50_ci)) stats$ec50_ci[["lower"]] else NA,
      if (!is.null(stats$ec50_ci)) stats$ec50_ci[["upper"]] else NA,
      as.numeric(stats$ec50_in_range %||% NA),
      stats$hill_slope %||% NA,
      stats$lower_limit %||% NA,
      stats$upper_limit %||% NA,
      stats$linear_aic %||% NA,
      stats$nonlinear_aic %||% NA,
      stats$growth_dose_p_value %||% NA,
      stats$growth_dose_r_squared %||% NA
    )
  )
  out
}

