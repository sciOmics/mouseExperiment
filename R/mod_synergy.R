#' Drug Synergy Analysis UI
#' @keywords internal
synergyUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Drug Synergy Analysis"),
    
    p("Test whether a drug combination shows synergistic, additive, or antagonistic effects ",
      "using the Bliss independence model. ",
      "Ensure tumor growth data is uploaded and configured in the Data Upload tab."),
    
    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Treatment Groups"),
          
          p(tags$small("Select which treatment group in your data corresponds to each role:")),
          
          uiOutput(ns("syn_group_selectors_ui")),
          
          hr(),
          
          h4("Analysis Options"),

          uiOutput(ns("syn_time_point_ui")),

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

                # ── Uncertainty & estimand ───────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-danger", "Uncertainty"), " — how far to trust the verdict")),
                  div(class = "panel-body",
                    p("The synergy verdict is a categorical label derived from three group means. Without an ",
                      "interval, a “Strong Synergy” resting on noise is indistinguishable on screen from one ",
                      "resting on a real effect."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Bootstrap resamples"),
                      tags$dd("Mouse-level resampling produces confidence intervals for the Bliss ",
                              "metrics. The control arm is resampled too, so uncertainty in the denominator ",
                              "propagates rather than being treated as known — which matters because every metric ",
                              "here is a ratio against control. Set to 0 to skip; the verdict is then a point ",
                              "estimate with no stated uncertainty."),
                      tags$dt("Interval callout"),
                      tags$dd("Beneath the verdict, an explicit statement of whether the 95 % interval for Bliss ",
                              "excess excludes zero. This is the single most useful line in the tab: a verdict ",
                              "whose interval spans zero is not evidence of synergy, whatever the label says."),
                      tags$dt("Endpoint estimand"),
                      tags$dd(strong("Model"), " (default) estimates endpoint group means from the fitted model, ",
                              "using every observation each animal contributed. ", strong("Observed"), " uses only ",
                              "animals still present at the endpoint. Because the fastest-growing controls are the ",
                              "ones removed early, the observed estimand shrinks the control mean and inflates every ",
                              "fractional effect derived from it — which propagates directly into both Bliss and ",
                              "the Bliss expectation. Per-arm attrition is shown so the exposure is visible.")
                    )
                  )
                ),

                # ── Over-time ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Over-Time Analysis")),
                  div(class = "panel-body",
                    p("When 'Analyze over time' is enabled in the sidebar, the Bliss score is computed at every measured timepoint and plotted as a synergy trend. This reveals whether synergy emerges, fades, or reverses over the course of treatment."),
                    p("Per-timepoint p-values (Combination vs Drug A and vs Drug B, two-sample t-test) are included in the Synergy Metrics table. These are unadjusted for multiple testing across timepoints — interpret overall patterns rather than individual p-values at isolated days.")
                  )
                ),

                br()
              )
            )
          ),

          tabPanel("Summary",
            tab_info(strong("Synergy summary: "), "Overall assessment from the Bliss independence model. Reports whether the drug combination is synergistic, additive, or antagonistic at the evaluated timepoint, and lists the key quantitative scores."),
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
                tags$li(strong("Interpretation"), " \u2014 plain-language summary based on the CI threshold.")
              ),
              "When \u2018Analyze over time\u2019 is enabled, each row is a measured timepoint with columns:",
              tags$ul(
                tags$li(strong("Time_Point"), " \u2014 study day."),
                tags$li(strong("TGI_Drug_A / TGI_Drug_B / TGI_Combo"), " \u2014 % tumor growth inhibition for Drug A, Drug B, and the combination relative to control."),
                tags$li(strong("Bliss_Expected_TGI"), " \u2014 predicted TGI under the Bliss independence model."),
                tags$li(strong("Bliss_Difference"), " \u2014 TGI_Combo minus the expected TGI; positive values indicate synergy."),
                tags$li(strong("P_Value_vs_Drug_A / P_Value_vs_Drug_B"), " \u2014 two-sample t-test p-value comparing the combination group to each mono-therapy; values below 0.05 indicate a statistically significant difference."),
                tags$li(strong("Synergy_Assessment"), " \u2014 plain-language classification at each timepoint."),
                tags$li(strong("Evaluable"), " \u2014 FALSE where a single agent did not inhibit growth that day. Bliss independence is defined only for inhibitory agents, so the Bliss columns are blank on those days and they cannot be the peak.")
              )
            ),
            div(class = "text-right", style = "margin-bottom: 6px;",
              downloadButton(ns("download_syn_metrics"), "Download Results")
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("syn_metrics_table"))
            ),

            # CODE_REVIEW.md R3.D9 / backend R3.6-R3.7 -- the verdict above is
            # derived from three group means. Without an interval a "Strong
            # Synergy" label whose Bliss-excess CI spans zero is
            # indistinguishable on screen from one that does not.
            # R14.5 / R17.2 -- the Combination Index was removed from the
            # backend: it computed response additivity rather than Loewe
            # additivity and labelled Bliss-additive combinations antagonistic
            # across most of the range where active single agents sit. The
            # verdict below is now the Bliss assessment, read against its
            # bootstrap interval.
            uiOutput(ns("syn_consensus_note")),
            uiOutput(ns("syn_ci_heading")),
            uiOutput(ns("syn_ci_verdict_note")),
            DTOutput(ns("syn_ci_table")),

            # Backend R3.5 -- animals leave these studies because their tumours
            # grew, so how many reached the evaluation day is part of reading
            # any TGI figure.
            uiOutput(ns("syn_attrition_heading")),
            DTOutput(ns("syn_attrition_table"))
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
              "Per-group normality and outlier checks for the t-test path."
            ),
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
  
  # D.7: this was an unconditional `rv$syn_results <- NULL`. `rv` is created empty
  # in app.R and each module mounts once, so the assignment never cleared anything
  # in production -- but it wiped the slot at module init, which made these five
  # modules impossible to drive through testServer() with a pre-populated result.
  # That is why I.2 stayed open across four rounds: every attempt to assert render
  # behaviour got the "run an analysis first" placeholder back. Removed; the slot
  # reads NULL until an analysis writes it either way.
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
  current_syn_tgi_plot   <- reactiveVal(NULL)

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
    
    # R20.D36: the choices are the arms the analysis uses ("Treatment Dose"
    # when a dose column is mapped), whether or not Tumour Growth has run.
    # They used to come from the raw column while the analysis switched to
    # "Treatment Dose" labels after a Tumour Growth run, so every run failed.
    groups  <- analysis_group_labels(rv)
    choices <- setNames(groups, groups)
    ctrl    <- study_reference(rv, groups)
    others  <- setdiff(groups, ctrl)
    combo   <- others[grepl("\\+|combo", others, ignore.case = TRUE)][1]
    if (is.na(combo)) combo <- others[3]
    singles <- setdiff(others, combo)

    tagList(
      selectInput(session$ns("syn_drug_a"), "Drug A (single agent)",
                 choices = choices, selected = singles[1]),

      selectInput(session$ns("syn_drug_b"), "Drug B (single agent)",
                 choices = choices, selected = singles[2]),

      selectInput(session$ns("syn_combo"), "Combination Treatment",
                 choices = choices, selected = combo),

      # The reference is the study's control group (R20-D.H).
      reference_note(rv, ctrl),

      # CODE_REVIEW.md backend R3.5 / G.3 -- the estimand, exposed. The default
      # is the statistically correct one; "survivors" exists to reproduce
      # pre-0.8.0 figures.
      selectInput(session$ns("syn_endpoint_method"), "Endpoint estimand",
                 choices = c(
                   "Model-based (recommended)" = "model",
                   "Each animal's last observation" = "last_obs",
                   "Survivors at the evaluation day (legacy)" = "survivors"),
                 selected = "model"),
      helpText(
        strong("Model-based"), " reads each arm's geometric mean at the ",
        "evaluation day from a log-scale mixed model fitted to every ",
        "observation, so animals removed earlier still contribute. ",
        strong("Survivors"), " uses only animals present on that day \u2014 which ",
        "conditions on survival and biases TGI downward, because animals leave ",
        "when their tumours grow. On a simulated study with volume-triggered ",
        "euthanasia the survivor estimand understated TGI by ~23 percentage ",
        "points where the model-based one was off by ~5."),

      numericInput(session$ns("syn_n_boot"),
                  "Bootstrap resamples (0 = skip)",
                  value = 2000, min = 0, max = 20000, step = 500),
      helpText("Mouse-level bootstrap giving 95% intervals for every synergy ",
               "metric. Set 0 for a faster point-estimate-only run.")
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
    
    syn_ctrl <- study_reference(rv, analysis_group_labels(rv))
    needed <- list(input$syn_drug_a, input$syn_drug_b, input$syn_combo, syn_ctrl)
    if (any(vapply(needed, function(x) is.null(x) || is.na(x) || !nzchar(x), logical(1L)))) {
      output$syn_status_msg <- renderUI({
        status_alert(paste("Please select Drug A, Drug B and the combination, and the",
                           "control group under Study Design on Data Upload."), "danger")
      })
      return()
    }
    needed <- unlist(needed)

    if (length(unique(needed)) < 4) {
      output$syn_status_msg <- renderUI({
        status_alert(paste("Each role (Drug A, Drug B, Combination, Control) must be a different",
                           "group; the control is the study's control group,", sQuote(syn_ctrl), "."),
                     "danger")
      })
      return()
    }
    
    output$syn_status_msg <- renderUI({
      status_alert("Running synergy analysis...", "info")
    })
    
    # Track analysis start with Umami
    shinyjs::runjs("window.trackAnalysisStart('drug_synergy', 'bliss', 'combo_vs_monotherapy');")
    
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
    # R20.D36: analyse the same arm labels the selectors offer.
    analysis_data[[".arm"]] <- arm_labels(analysis_data, rv$settings$col_treatment,
                                          rv$settings$col_dose)
    tx_col   <- ".arm"
    cage_col_syn <- if (has_col(analysis_data, rv$settings$col_cage)) rv$settings$col_cage else NULL
    # Animal identifier: needed by the backend's mouse-level bootstrap and by its
    # endpoint helper. Fall back to the mapped ID column, then to "ID".
    id_col_syn <- if (!is.null(rv$settings$col_id) && nzchar(rv$settings$col_id)) {
      rv$settings$col_id
    } else if ("ID" %in% colnames(analysis_data)) "ID" else NULL
    if (is.null(id_col_syn)) {
      showNotification(
        paste0("No animal ID column is mapped, so the synergy bootstrap and the ",
               "model-based endpoint cannot resample animals. Map an ID column in ",
               "Data Upload to enable the confidence intervals."),
        type = "warning", duration = 10)
      id_col_syn <- "ID"
    }
    
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
          control_name     = syn_ctrl,
          # R20.D14: the ID column (it had to be named "ID"), cage, and the
          # estimand and bootstrap settings, which "All dates" used to ignore.
          id_column        = id_col_syn,
          cage_column      = cage_col_syn,
          endpoint_method  = input$syn_endpoint_method %||% "model",
          n_boot           = as.integer(input$syn_n_boot %||% 2000L),
          boot_seed        = 20260730L,
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
          id_column        = id_col_syn,
          cage_column      = cage_col_syn,  # part of the animal key
          drug_a_name      = input$syn_drug_a,
          drug_b_name      = input$syn_drug_b,
          combo_name       = input$syn_combo,
          control_name     = syn_ctrl,
          eval_time_point  = eval_tp,
          # Backend v0.8.0+: "model" reads each arm's geometric mean at the
          # evaluation day off a log-scale LMM fitted to every observation, so
          # animals removed earlier still contribute. "survivors" reproduces the
          # pre-0.8.0 raw survivor mean.
          endpoint_method  = input$syn_endpoint_method %||% "model",
          # Backend v0.7.0+: mouse-level bootstrap so the verdict carries an
          # interval instead of standing on a bare point estimate.
          n_boot           = as.integer(input$syn_n_boot %||% 2000L),
          boot_seed        = 20260730L,
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
      control = syn_ctrl
    )
    
    if (isTRUE(result$error)) {
      output$syn_status_msg <- renderUI({
        status_alert(paste("Error:", result$message), "danger")
      })
    } else {
      # R20.D34: the backend's warnings are shown, not dropped.
      output$syn_status_msg <- renderUI(completion_status(
        list(result), "Synergy analysis completed successfully!"))
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
      choices <- c("Synergy Plot" = "trend",
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
      # The roles of the displayed run, not the current selectors.
      roles <- rv$syn_plot_inputs %||% list()
      cat("Drug A (single agent): ", roles$drug_a %||% "", "\n")
      cat("Drug B (single agent): ", roles$drug_b %||% "", "\n")
      cat("Combination:           ", roles$combo %||% "", "\n")
      cat("Control:               ", roles$control %||% "", "\n\n")
      
      if (isTRUE(rv$syn_over_time)) {
        cat("Analysis type: Over all timepoints\n\n")
        # R20.8: the peak is taken over evaluable days only, and there is none
        # when a single agent never inhibited growth.
        pb <- res$peak_bliss_synergy
        n_days <- NROW(res$synergy_summary)
        n_eval <- sum(res$synergy_summary$Evaluable %||% rep(TRUE, n_days))
        if (n_eval < n_days)
          cat("Evaluable timepoints:", n_eval, "of", n_days,
              "(Bliss is not evaluated where a single agent did not inhibit growth)\n")
        if (!is.null(pb) && NROW(pb) > 0L) {
          cat("Peak Bliss synergy at timepoint:", pb$Time_Point,
              " (Bliss Difference =", round(pb$Bliss_Difference, 1), "%)\n")
        } else if (!is.null(pb)) {
          cat("No Bliss peak: no timepoint was evaluable.\n")
        }
      } else {
        cat("Analysis type: Single timepoint\n")
        if (!is.null(res$evaluation_time_point)) cat("Evaluated at Day:", res$evaluation_time_point, "\n")
        cat("\nOverall Assessment:", res$overall_assessment, "\n\n")
        if (!is.null(res$bliss_independence)) {
          bi <- res$bliss_independence
          cat("Bliss Independence:\n")
          # R20.8: when an agent did not inhibit growth the backend reports every
          # Bliss quantity as NA; say so rather than printing "NA" as a verdict.
          if (isFALSE(res$evaluable) || is.na(bi$synergy)) {
            cat("  Not evaluable: a single agent did not inhibit growth, and Bliss\n",
                "  independence is defined only for inhibitory agents.\n", sep = "")
            cat("  Observed effect:", round(bi$observed_effect, 4), "\n")
          } else {
            cat("  Expected effect:", round(bi$expected_effect, 4), "\n")
            cat("  Observed effect:", round(bi$observed_effect, 4), "\n")
            cat("  Synergistic:    ", bi$synergy, "\n")
          }
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
  
  # ── Bootstrap intervals + attrition (CODE_REVIEW.md R3.D9) ─────────────────
  #
  # The Synergy Metrics table above states a categorical verdict derived from
  # three group means. Backend v0.7.0 attaches a mouse-level bootstrap so that
  # verdict can be read against its uncertainty, and v0.8.0 reports how many
  # animals actually reached the evaluation day. Neither was surfaced.

  output$syn_ci_heading <- renderUI({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || is.null(rv$syn_results$synergy_ci)) return(NULL)
    tagList(
      tags$hr(),
      h4("Uncertainty \u2014 mouse-level bootstrap 95% intervals"),
      p(class = "text-muted",
        "Animals are resampled within group, ", strong("including the control arm"),
        ", and every metric recomputed. Resampling the control propagates the ",
        "uncertainty of the TGI denominator, which a point estimate treats as a ",
        "known constant.")
    )
  })

  output$syn_ci_table <- renderDT({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    ci <- rv$syn_results$synergy_ci
    if (is.null(ci)) return(NULL)
    df <- format_df_for_display(as.data.frame(ci))
    datatable(df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE),
              rownames = FALSE)
  })

  # The single most important reading aid: does the Bliss-excess interval exclude
  # zero? If it does not, the verdict above is not evidence of synergy.
  # R14.5 / R17.2: the headline verdict. It used to sit above a separate
  # Combination Index verdict and flag disagreement between the two models; the
  # CI was removed from the backend because it computed response additivity
  # rather than Loewe additivity and labelled Bliss-additive combinations
  # antagonistic across most of the range where active single agents sit. One
  # model now, so there is nothing to reconcile -- but the label is still a point
  # estimate, which is why the interval note below it matters more than this box.
  output$syn_consensus_note <- renderUI({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    oa <- rv$syn_results$overall_assessment
    if (is.null(oa) || !nzchar(oa)) return(NULL)

    style <- if (grepl("Not evaluable", oa)) "warning"
             else if (grepl("Synergy", oa)) "success"
             else if (grepl("Antagonism", oa)) "danger" else "info"

    div(class = paste0("alert alert-", style),
      tags$strong("Synergy verdict: "), oa,
      tags$br(),
      tags$small(
        "From the Bliss independence model: expected effect = ",
        tags$code("FE_A + FE_B \u2212 FE_A \u00d7 FE_B"),
        ", so the drugs are assumed to act independently on the surviving ",
        "fraction. This is a threshold applied to three group means \u2014 read it ",
        "against the interval below, not on its own."
      )
    )
  })

  output$syn_ci_verdict_note <- renderUI({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    ci <- rv$syn_results$synergy_ci
    if (is.null(ci) || !"Metric" %in% names(ci)) return(NULL)
    row <- ci[ci$Metric == "Bliss_Excess_FE", , drop = FALSE]
    if (!nrow(row) || !is.finite(row$CI_Lower[1]) || !is.finite(row$CI_Upper[1])) {
      return(NULL)
    }
    spans_zero <- row$CI_Lower[1] <= 0 && row$CI_Upper[1] >= 0
    if (spans_zero) {
      status_alert(paste0(
        "The Bliss-excess 95% interval [",
        signif(row$CI_Lower[1], 3), ", ", signif(row$CI_Upper[1], 3),
        "] includes zero, so this data does not distinguish the combination from ",
        "Bliss additivity. Read the synergy label above as descriptive only."),
        "warning")
    } else {
      status_alert(paste0(
        "The Bliss-excess 95% interval [",
        signif(row$CI_Lower[1], 3), ", ", signif(row$CI_Upper[1], 3),
        "] excludes zero, consistent with a departure from Bliss additivity in ",
        "the ", if (row$CI_Lower[1] > 0) "synergistic" else "antagonistic",
        " direction."),
        "success")
    }
  })

  output$syn_attrition_heading <- renderUI({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error) || is.null(rv$syn_results$attrition)) return(NULL)
    a <- rv$syn_results$attrition
    lost <- if (!is.null(a$Pct_Lost)) max(a$Pct_Lost, na.rm = TRUE) else 0
    tagList(
      tags$hr(),
      h4("Animals contributing at the evaluation day"),
      p(class = "text-muted",
        "Animals leave these studies because their tumours grew, so survivors are ",
        "systematically the slower growers \u2014 most severely in the control arm, ",
        "which loses animals first. The endpoint method shown in the sidebar ",
        "controls whether the analysis conditions on survival or models through it."),
      if (is.finite(lost) && lost > 0) {
        status_alert(sprintf(
          paste0("Up to %.0f%% of animals in one arm did not reach the evaluation ",
                 "day. With endpoint method '%s' this %s."),
          lost, rv$syn_results$endpoint_method %||% "model",
          if (identical(rv$syn_results$endpoint_method, "survivors")) {
            "biases TGI downward -- switch to 'model' unless you are reproducing an older figure"
          } else {
            "is handled by fitting through the dropout rather than discarding it"
          }),
          if (identical(rv$syn_results$endpoint_method, "survivors")) "warning" else "info")
      }
    )
  })

  output$syn_attrition_table <- renderDT({
    req(rv$syn_results)
    if (isTRUE(rv$syn_results$error)) return(NULL)
    a <- rv$syn_results$attrition
    if (is.null(a)) return(NULL)
    datatable(format_df_for_display(as.data.frame(a)),
              options = list(paging = FALSE, searching = FALSE), rownames = FALSE)
  })

  # TGI Table removed (duplicate of Synergy Metrics)

  output$syn_download_plot_btn <- renderUI({
    btn_id <- switch(input$syn_plot_selector %||% "main",
      main   = "download_syn_plot_main",
      trend  = "download_syn_plot_trend",
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
    # The run's roles are arm labels (R20.D36).
    analysis_data[[".arm"]] <- arm_labels(analysis_data, rv$settings$col_treatment,
                                          rv$settings$col_dose)
    tx_col   <- ".arm"

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
                            current_syn_tgi_plot)

  # CODE_REVIEW.md DIAGNOSTICS gap (6) — frequentist synergy diagnostics.
  output$syn_diag_group_qq_plot <- renderPlot({
    req(rv$syn_results)
    p <- rv$syn_results$diag_group_qq_plot
    if (is.null(p)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Per-group Q-Q not available (over-time path)")
      return(invisible(NULL))
    }
    p
  })

  output$syn_diag_group_boxplot <- renderPlot({
    req(rv$syn_results)
    p <- rv$syn_results$diag_group_boxplot
    if (is.null(p)) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Per-group boxplot not available (over-time path)")
      return(invisible(NULL))
    }
    p
  })

  })
}
