#' Power Analysis UI
#' @keywords internal
powerAnalysisUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Power Analysis"),

    p("Estimate statistical power and sample sizes from your data or from user-supplied parameters.",
      " Upload and configure data in the Data Upload tab first (required for post-hoc mode)."),

    fluidRow(
      # ---- Sidebar -----------------------------------------------------------
      column(3,
        wellPanel(

          # ---- Mode selector -------------------------------------------------
          radioButtons(ns("pw_mode"),
                       label    = "Analysis Mode",
                       choices  = c("Post-hoc (observed data)"     = "posthoc",
                                    "A priori \u2014 Analytic"      = "apriori_analytic",
                                    "A priori \u2014 LMM Simulation" = "apriori_lmm"),
                       selected = "posthoc"),

          hr(),

          # ---- POST-HOC controls ---------------------------------------------
          conditionalPanel(ns = ns, condition = "input.pw_mode == 'posthoc'",
            h4("Analysis Method"),

            selectInput(ns("pw_method"),
                        "Method",
                        choices = c(
                          "Parametric" = "parametric",
                          "Simulation" = "simulation",
                          "AUC-based"  = "auc"
                        ),
                        selected = "parametric"),

            conditionalPanel(ns = ns,
              condition = "input.pw_method == 'simulation'",
              numericInput(ns("pw_n_simulations"),
                           "Number of Simulations",
                           value = 1000, min = 100, max = 10000, step = 100)
            ),

            hr(),

            h4("Custom Effect Sizes"),
            textInput(ns("pw_custom_effect_sizes"),
                      "Effect sizes (comma-separated, optional)",
                      value       = "",
                      placeholder = "e.g.\u202f0.2,\u202f0.5,\u202f0.8,\u202f1.0"),
            p(tags$small("Overrides auto-detected effect sizes. Leave blank to use values estimated from the data.")),

            hr()
          ),

          # ---- A PRIORI ANALYTIC controls ------------------------------------
          conditionalPanel(ns = ns, condition = "input.pw_mode == 'apriori_analytic'",

            h4("Effect Size"),

            radioButtons(ns("pw_ap_es_mode"),
                         label    = NULL,
                         choices  = c("Cohen\u2019s d directly"          = "cohens_d",
                                      "Mean difference + Pooled SD"       = "delta_sd"),
                         selected = "cohens_d"),

            conditionalPanel(ns = ns, condition = "input.pw_ap_es_mode == 'cohens_d'",
              textInput(ns("pw_ap_cohens_d"),
                        "Cohen\u2019s d values (comma-separated)",
                        value       = "0.5, 0.8, 1.0, 1.5",
                        placeholder = "e.g.\u202f0.5,\u202f0.8,\u202f1.0"),
              p(tags$small("Benchmarks: small\u202f\u2248\u202f0.2, medium\u202f\u2248\u202f0.5, large\u202f\u2248\u202f0.8"))
            ),

            conditionalPanel(ns = ns, condition = "input.pw_ap_es_mode == 'delta_sd'",
              numericInput(ns("pw_ap_delta"),
                           "Expected mean difference (\u0394)",
                           value = 100, min = 0, step = 1),
              numericInput(ns("pw_ap_pooled_sd"),
                           "Pooled within-group SD",
                           value = 50, min = 0.1, step = 0.1),
              actionButton(ns("pw_seed_sd"),
                           label = tagList(icon("database"), "Seed SD from loaded data"),
                           class = "btn btn-sm btn-default btn-block")
            ),

            hr(),

            selectInput(ns("pw_ap_n_groups"),
                        "Number of treatment groups",
                        choices  = as.character(2:6),
                        selected = "2"),

            hr(),

            radioButtons(ns("pw_ap_submode"),
                         label    = "Solve for:",
                         choices  = c("Required N (given target power)" = "find_n",
                                      "Achieved power (given N)"        = "find_power"),
                         selected = "find_n"),

            conditionalPanel(ns = ns, condition = "input.pw_ap_submode == 'find_power'",
              textInput(ns("pw_ap_n_values"),
                        "N per group to evaluate (comma-separated)",
                        value       = "5, 8, 10, 12, 15",
                        placeholder = "e.g.\u202f5,\u202f8,\u202f10")
            ),

            hr()
          ),

          # ---- A PRIORI LMM SIMULATION controls ------------------------------
          conditionalPanel(ns = ns, condition = "input.pw_mode == 'apriori_lmm'",

            h4("Tumor Growth Parameters"),

            numericInput(ns("pw_lmm_baseline_mean"),
                         "Baseline tumor volume (mean)",
                         value = 100, min = 1, step = 1),
            numericInput(ns("pw_lmm_baseline_sd"),
                         "Baseline tumor volume (SD)",
                         value = 20, min = 0.1, step = 0.1),
            numericInput(ns("pw_lmm_ctrl_rate"),
                         "Control growth rate (volume units/day)",
                         value = 15, min = 0, step = 0.5),
            numericInput(ns("pw_lmm_tx_effect"),
                         "Treatment rate reduction (units/day)",
                         value = 10, min = 0, step = 0.5),
            p(tags$small("Treated groups grow at [control rate \u2212 this value] per day.")),

            hr(),

            h4("Variance Components"),

            numericInput(ns("pw_lmm_rand_int_sd"),
                         "Random intercept SD",
                         value = 10, min = 0, step = 0.1),
            numericInput(ns("pw_lmm_rand_slope_sd"),
                         "Random slope SD",
                         value = 3, min = 0, step = 0.1),
            numericInput(ns("pw_lmm_resid_sd"),
                         "Residual (measurement) SD",
                         value = 10, min = 0.1, step = 0.1),

            hr(),

            h4("Study Design"),

            numericInput(ns("pw_lmm_n_groups"),
                         "Number of groups (incl. control)",
                         value = 2, min = 2, max = 8, step = 1),
            textInput(ns("pw_lmm_n_per_group"),
                      "N per group to evaluate (comma-separated)",
                      value       = "5, 8, 10, 12, 15",
                      placeholder = "e.g.\u202f5,\u202f8,\u202f10,\u202f12,\u202f15"),
            numericInput(ns("pw_lmm_n_timepoints"),
                         "Number of measurement days",
                         value = 15, min = 2, max = 100, step = 1),
            numericInput(ns("pw_lmm_day_start"),
                         "First measurement day",
                         value = 0, min = 0, step = 1),
            numericInput(ns("pw_lmm_day_interval"),
                         "Interval between days",
                         value = 1, min = 1, step = 1),

            hr(),

            h4("Simulation"),

            numericInput(ns("pw_lmm_n_sims"),
                         "Number of simulations",
                         value = 500, min = 50, max = 2000, step = 50),
            p(tags$small(icon("exclamation-triangle"),
                         tags$strong(" Computationally intensive."),
                         " Use \u2264500 simulations for quick estimates.")),

            hr()
          ),

          # ---- Shared: Alpha (all modes) ------------------------------------
          h4("Alpha Levels"),
          checkboxGroupInput(ns("pw_alpha"),
                             "Significance Levels",
                             choices  = c("0.05" = "0.05",
                                          "0.01" = "0.01",
                                          "0.10" = "0.10"),
                             selected = c("0.05", "0.01")),

          # ---- Shared: Target Power (post-hoc + analytic only) --------------
          conditionalPanel(ns = ns, condition = "input.pw_mode != 'apriori_lmm'",
            hr(),
            h4("Target Power"),
            checkboxGroupInput(ns("pw_power"),
                               "Power Targets",
                               choices  = c("0.80" = "0.80",
                                            "0.85" = "0.85",
                                            "0.90" = "0.90",
                                            "0.95" = "0.95"),
                               selected = c("0.80", "0.90"))
          ),

          br(),

          actionButton(ns("run_pw_analysis"),
                       "Run Analysis",
                       class = "btn-lg btn-primary btn-block"),

          uiOutput(ns("pw_status_msg"))
        )
      ),

      # ---- Results area ------------------------------------------------------
      column(9,

        # Post-hoc results panel
        conditionalPanel(ns = ns, condition = "input.pw_mode == 'posthoc'",
          tabsetPanel(
            id = ns("pw_results_tabs"),

            tabPanel("Summary",
              tab_info(strong("Power analysis summary: "), "overall observed power for the current sample size and key recommendations. Post-hoc power is calculated from the effect sizes observed in your data \u2014 interpret cautiously as the same data is used to both estimate effects and compute power."),
              shinycssloaders::withSpinner(
                verbatimTextOutput(ns("pw_summary"))
              )
            ),

            tabPanel("Effect Sizes",
              tab_info(
                strong("Cohen\u2019s\u202fd effect sizes "), "between each treatment group pair estimated from the observed data. Conventional benchmarks: ",
                strong("small\u202f\u2248\u202f0.2"), ", ",
                strong("medium\u202f\u2248\u202f0.5"), ", ",
                strong("large\u202f\u2248\u202f0.8"), ". Larger effect sizes require smaller samples to detect reliably.",
                tags$ul(
                  tags$li(strong("Treatment"), " \u2014 the treatment group being compared."),
                  tags$li(strong("Reference"), " \u2014 the group it is compared against."),
                  tags$li(strong("Raw_Difference"), " \u2014 unadjusted mean difference in the outcome measure between Treatment and Reference."),
                  tags$li(strong("Pooled_SD"), " \u2014 pooled within-group standard deviation used as the denominator for Cohen\u2019s d."),
                  tags$li(strong("Standardized_Effect"), " \u2014 Cohen\u2019s d = Raw_Difference \u00f7 Pooled_SD. Positive values mean Treatment > Reference.")
                )
              ),
              shinycssloaders::withSpinner(
                DTOutput(ns("pw_effect_sizes_table"))
              ),
              downloadButton(ns("download_pw_effect_sizes"), "Download")
            ),

            tabPanel("Power Analysis",
              tab_info(
                strong("Achieved power (1\u202f\u2212\u202f\u03b2) "), "at the current observed N per group, for each combination of comparison, alpha level, and effect size scenario.",
                tags$ul(
                  tags$li(strong("Treatment / Reference"), " \u2014 the pair of groups being evaluated."),
                  tags$li(strong("Effect_Size"), " \u2014 Cohen\u2019s d used for this row (drawn from the observed effect sizes above or from a standard scenario grid)."),
                  tags$li(strong("Alpha"), " \u2014 significance threshold (typically 0.05)."),
                  tags$li(strong("Power"), " \u2014 probability of correctly detecting an effect of this size at this alpha with the current N. Values below 0.80 suggest the study may be underpowered.")
                ),
                em("Post-hoc power should be interpreted cautiously: the same data used to fit the model is used to estimate the effect size, which inflates apparent power.")
              ),
              shinycssloaders::withSpinner(
                DTOutput(ns("pw_power_table"))
              ),
              downloadButton(ns("download_pw_power_table"), "Download")
            ),

            tabPanel("Sample Size Recommendations",
              tab_info(
                strong("N per group needed "), "to achieve each target power level, for each combination of comparison, alpha, and effect size. Use these to plan future experiments.",
                tags$ul(
                  tags$li(strong("Treatment / Reference"), " \u2014 the pair being evaluated."),
                  tags$li(strong("Effect_Size"), " \u2014 Cohen\u2019s d for this scenario."),
                  tags$li(strong("Alpha"), " \u2014 significance threshold."),
                  tags$li(strong("Target_Power"), " \u2014 the desired power level (e.g. 0.80 or 0.90)."),
                  tags$li(strong("Sample_Size"), " \u2014 N per group required to achieve the target power. Estimates for small effects or stringent targets may be impractically large.")
                )
              ),
              shinycssloaders::withSpinner(
                DTOutput(ns("pw_sample_sizes_table"))
              ),
              downloadButton(ns("download_pw_sample_sizes"), "Download")
            ),

            tabPanel("Power Curves",
              tab_info(strong("Power vs N per group: "), "shows how power increases with sample size for each detected effect size. Horizontal dashed lines mark conventional power thresholds (0.80, 0.90). Use these curves to identify the minimum N that achieves adequate power for your effect of interest."),
              shinycssloaders::withSpinner(
                plotOutput(ns("pw_plot_curves"), height = "500px")
              ),
              downloadButton(ns("download_pw_plot"), "Download Plot")
            )
          )
        ),

        # A priori results panel (analytic + LMM)
        conditionalPanel(ns = ns, condition = "input.pw_mode != 'posthoc'",
          tabsetPanel(
            id = ns("pw_ap_results_tabs"),

            tabPanel("Scenario Table",
              tab_info(
                strong("A priori power / sample size table: "),
                "For each combination of effect size, alpha, and target power, shows the required N per group (Find N mode) or achieved power (Find Power / LMM mode).",
                tags$ul(
                  tags$li(strong("Effect_Size"), " \u2014 Cohen\u2019s d (analytic) or the simulated scenario (LMM)."),
                  tags$li(strong("N_Groups"), " \u2014 number of treatment groups in the design."),
                  tags$li(strong("Alpha"), " \u2014 significance threshold."),
                  tags$li(strong("Required_N"), " \u2014 N per group needed to reach the target power (Find N mode)."),
                  tags$li(strong("Achieved_Power / Power"), " \u2014 power reached at the specified N (Find Power / LMM mode).")
                )
              ),
              shinycssloaders::withSpinner(
                DTOutput(ns("pw_ap_scenario_table"))
              ),
              downloadButton(ns("download_pw_ap_scenario"), "Download")
            ),

            tabPanel("Power Curves",
              tab_info(strong("Power vs N per group: "), "prospective power curves for the specified effect size scenario(s). Each line represents one Cohen\u2019s d / alpha combination. Dashed lines mark 0.80 and 0.90 power thresholds."),
              shinycssloaders::withSpinner(
                plotOutput(ns("pw_ap_plot_curves"), height = "500px")
              ),
              downloadButton(ns("download_pw_ap_plot"), "Download Plot")
            ),

            tabPanel("SD Sensitivity",
              tab_info(
                strong("Sensitivity to the assumed SD: "), "shows how the required N changes when the pooled within-group SD used to convert a raw mean difference into Cohen\u2019s d is perturbed by \u00b120\u202f% and \u00b140\u202f%. Only available when using the Mean\u202fdifference\u202f+\u202fPooled\u202fSD input mode.",
                tags$ul(
                  tags$li(strong("SD_Change"), " \u2014 percentage change applied to the assumed SD."),
                  tags$li(strong("Assumed_SD"), " \u2014 resulting SD value."),
                  tags$li(strong("Cohens_d"), " \u2014 Cohen\u2019s d after adjusting for the new SD."),
                  tags$li(strong("Required_N"), " \u2014 sample size required at the reference alpha / target power.")
                )
              ),
              shinycssloaders::withSpinner(
                DTOutput(ns("pw_ap_sensitivity_table"))
              ),
              downloadButton(ns("download_pw_ap_sensitivity"), "Download")
            )
          )
        )
      )
    )
  )
}

#' Power Analysis Server Logic
#' @keywords internal
power_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
  
  current_pw_plot    <- reactiveVal(NULL)
  current_ap_plot    <- reactiveVal(NULL)

  rv$pw_results <- NULL
  rv$ap_results <- NULL

  # ---- Seed pooled SD from loaded data (analytic a priori) ------------------
  observeEvent(input$pw_seed_sd, {
    req(rv$processed_data)
    df      <- rv$processed_data
    tx_col  <- rv$tx_col
    vol_col <- rv$vol_col
    if (!has_col(df, tx_col) || !has_col(df, vol_col)) {
      showNotification("Could not find treatment or volume column in processed data.",
                       type = "warning", duration = 4)
      return()
    }
    group_sds  <- tapply(df[[vol_col]], df[[tx_col]], stats::sd,  na.rm = TRUE)
    group_ns   <- tapply(df[[vol_col]], df[[tx_col]], function(x) sum(!is.na(x)))
    wvar_sum   <- sum((group_ns - 1) * group_sds^2, na.rm = TRUE)
    df_total   <- sum(group_ns - 1, na.rm = TRUE)
    pooled_sd_val <- if (df_total > 0) sqrt(wvar_sum / df_total) else mean(group_sds, na.rm = TRUE)
    if (is.na(pooled_sd_val) || pooled_sd_val <= 0) {
      showNotification("Could not compute a valid pooled SD from the data.", type = "warning", duration = 4)
      return()
    }
    updateNumericInput(session, "pw_ap_pooled_sd", value = round(pooled_sd_val, 2))
    showNotification(
      paste0("Pooled SD seeded from data: ", round(pooled_sd_val, 2)),
      type = "message", duration = 3
    )
  })

  observeEvent(input$run_pw_analysis, {

    # ---- Shared validation ---------------------------------------------------
    if (length(input$pw_alpha) == 0) {
      output$pw_status_msg <- renderUI({
        status_alert("Please select at least one alpha level.", "danger")
      })
      return()
    }

    alpha_vals <- as.numeric(input$pw_alpha)

    output$pw_status_msg <- renderUI({
      status_alert("Running analysis\u2026", "info")
    })

    # =========================================================================
    # BRANCH A — POST-HOC (observed data)
    # =========================================================================
    if (input$pw_mode == "posthoc") {

      req(rv$uploaded_data, rv$data_valid)

      if (length(input$pw_power) == 0) {
        output$pw_status_msg <- renderUI({
          status_alert("Please select at least one power target.", "danger")
        })
        return()
      }

      shinyjs::runjs(sprintf("window.trackAnalysisStart('power', '%s', 'post_hoc');",
                             input$pw_method))

      if (is.null(rv$processed_data)) {
        output$pw_status_msg <- renderUI({
          status_alert("Please run the Tumor Growth Analysis first so that Volume and Day columns are resolved.", "warning")
        })
        return()
      }

      time_col   <- rv$time_col
      vol_col    <- rv$vol_col
      tx_col     <- rv$tx_col
      id_col     <- rv$id_col
      power_vals <- as.numeric(input$pw_power)
      n_sims     <- if (input$pw_method == "simulation") input$pw_n_simulations else 1000

      custom_effect_sizes <- NULL
      if (!is.null(input$pw_custom_effect_sizes) && nzchar(trimws(input$pw_custom_effect_sizes))) {
        parsed <- suppressWarnings(as.numeric(strsplit(trimws(input$pw_custom_effect_sizes), "[,\\s]+")[[1]]))
        parsed <- parsed[!is.na(parsed) & parsed > 0]
        if (length(parsed) > 0) custom_effect_sizes <- parsed
      }

      progress <- shiny::Progress$new(session, min = 0, max = 1)
      progress$set(message = "Running post-hoc power analysis\u2026", value = 0.2)
      on.exit(progress$close(), add = TRUE)

      result <- safe_analysis(module = "power", {
        mouseExperiment::post_power_analysis(
          data             = rv$processed_data,
          alpha            = alpha_vals,
          power            = power_vals,
          method           = input$pw_method,
          n_simulations    = n_sims,
          time_column      = time_col,
          volume_column    = vol_col,
          treatment_column = tx_col,
          id_column        = id_col,
          effect_sizes     = custom_effect_sizes
        )
      })

      if (!isTRUE(result$error)) {
        if (is.null(result$summary)) result$summary <- list()
        result$summary$method        <- input$pw_method
        result$summary$alpha_levels  <- alpha_vals
        result$summary$power_targets <- power_vals
      }

      rv$pw_results <- result

      if (isTRUE(result$error)) {
        shinyjs::runjs(sprintf("window.trackAnalysisComplete('power', false, '%s');",
                               gsub("'", "\\\\'", result$message)))
        output$pw_status_msg <- renderUI({
          status_alert(paste("Error:", result$message), "danger")
        })
      } else {
        shinyjs::runjs("window.trackAnalysisComplete('power', true);")
        output$pw_status_msg <- renderUI({
          status_alert("Post-hoc power analysis completed.", "success")
        })
      }

    # =========================================================================
    # BRANCH B — A PRIORI ANALYTIC
    # =========================================================================
    } else if (input$pw_mode == "apriori_analytic") {

      if (length(input$pw_power) == 0 && input$pw_ap_submode == "find_n") {
        output$pw_status_msg <- renderUI({
          status_alert("Please select at least one power target.", "danger")
        })
        return()
      }

      shinyjs::runjs("window.trackAnalysisStart('power', 'analytic', 'apriori');")

      # Resolve effect sizes
      es_vals <- NULL
      if (input$pw_ap_es_mode == "cohens_d") {
        raw_txt <- trimws(input$pw_ap_cohens_d %||% "")
        if (!nzchar(raw_txt)) {
          output$pw_status_msg <- renderUI({
            status_alert("Please enter at least one Cohen\u2019s d value.", "danger")
          })
          return()
        }
        parsed <- suppressWarnings(as.numeric(strsplit(raw_txt, "[,\\s]+")[[1]]))
        es_vals <- parsed[!is.na(parsed) & parsed > 0]
        if (length(es_vals) == 0) {
          output$pw_status_msg <- renderUI({
            status_alert("No valid positive Cohen\u2019s d values found.", "danger")
          })
          return()
        }
      }

      power_vals <- if (input$pw_ap_submode == "find_n") as.numeric(input$pw_power) else c(0.80, 0.90)
      n_vals     <- if (input$pw_ap_submode == "find_power") {
        raw_n <- suppressWarnings(as.integer(strsplit(trimws(input$pw_ap_n_values %||% ""), "[,\\s]+")[[1]]))
        raw_n[!is.na(raw_n) & raw_n >= 2]
      } else NULL

      n_groups <- as.integer(input$pw_ap_n_groups %||% 2)

      progress <- shiny::Progress$new(session, min = 0, max = 1)
      progress$set(message = "Computing a priori power\u2026", value = 0.1)
      on.exit(progress$close(), add = TRUE)

      # Call apriori_power_analysis() once per effect size; combine results
      if (input$pw_ap_es_mode == "cohens_d") {
        ap_runs <- lapply(es_vals, function(d) {
          tryCatch(
            mouseExperiment::apriori_power_analysis(
              effect_size          = d,
              n_groups             = n_groups,
              alpha                = alpha_vals,
              target_power         = power_vals,
              n_per_group          = n_vals,
              mode                 = input$pw_ap_submode,
              sd_sensitivity_range = 0   # no SD to vary
            ),
            error = function(e) list(error = TRUE, message = conditionMessage(e))
          )
        })
      } else {
        # delta + SD path: single call, sensitivity table available
        delta_val <- input$pw_ap_delta %||% 0
        sd_val    <- input$pw_ap_pooled_sd %||% 1
        if (sd_val <= 0) {
          output$pw_status_msg <- renderUI({
            status_alert("Pooled SD must be > 0.", "danger")
          })
          return()
        }
        ap_runs <- list(
          tryCatch(
            mouseExperiment::apriori_power_analysis(
              delta                = delta_val,
              pooled_sd            = sd_val,
              n_groups             = n_groups,
              alpha                = alpha_vals,
              target_power         = power_vals,
              n_per_group          = n_vals,
              mode                 = input$pw_ap_submode,
              sd_sensitivity_range = 0.4
            ),
            error = function(e) list(error = TRUE, message = conditionMessage(e))
          )
        )
      }

      # Check for errors
      errors <- Filter(function(x) isTRUE(x$error), ap_runs)
      if (length(errors) > 0) {
        shinyjs::runjs("window.trackAnalysisComplete('power', false, 'analytic error');")
        output$pw_status_msg <- renderUI({
          status_alert(paste("Error:", errors[[1]]$message), "danger")
        })
        return()
      }

      scenario_df  <- do.call(rbind, lapply(ap_runs, `[[`, "scenario_table"))
      curve_df     <- do.call(rbind, lapply(ap_runs, `[[`, "power_curve_data"))
      sens_list    <- Filter(Negate(is.null), lapply(ap_runs, `[[`, "sensitivity_table"))
      sens_df      <- if (length(sens_list) > 0) do.call(rbind, sens_list) else NULL

      rv$ap_results <- list(
        scenario_table    = scenario_df,
        power_curve_data  = curve_df,
        sensitivity_table = sens_df,
        mode              = "apriori_analytic",
        ap_submode        = input$pw_ap_submode
      )

      progress$set(value = 1)
      shinyjs::runjs("window.trackAnalysisComplete('power', true);")
      output$pw_status_msg <- renderUI({
        status_alert("A priori power analysis completed.", "success")
      })

    # =========================================================================
    # BRANCH C — A PRIORI LMM SIMULATION
    # =========================================================================
    } else if (input$pw_mode == "apriori_lmm") {

      shinyjs::runjs("window.trackAnalysisStart('power', 'lmm_sim', 'apriori');")

      n_per_grp_lmm <- suppressWarnings(
        as.integer(strsplit(trimws(input$pw_lmm_n_per_group %||% ""), "[,\\s]+")[[1]])
      )
      n_per_grp_lmm <- n_per_grp_lmm[!is.na(n_per_grp_lmm) & n_per_grp_lmm >= 2]
      if (length(n_per_grp_lmm) == 0) {
        output$pw_status_msg <- renderUI({
          status_alert("Please enter at least one valid N per group value.", "danger")
        })
        return()
      }

      tp_start    <- input$pw_lmm_day_start    %||% 0
      tp_interval <- input$pw_lmm_day_interval %||% 1
      n_tp        <- input$pw_lmm_n_timepoints %||% 15
      timepoints  <- seq(tp_start, by = tp_interval, length.out = n_tp)

      progress <- shiny::Progress$new(session, min = 0, max = 1)
      progress$set(message = "Running LMM simulation\u2026 (this may take a while)", value = 0.05)
      on.exit(progress$close(), add = TRUE)

      total_scenarios  <- length(n_per_grp_lmm) * length(alpha_vals)
      completed_scens  <- 0L
      progress_fn_shiny <- function(prop) {
        completed_scens <<- completed_scens + 1L
        progress$set(value = 0.05 + 0.90 * prop,
                     message = sprintf("LMM simulation: scenario %d / %d\u2026",
                                       completed_scens, total_scenarios))
      }

      result_lmm <- safe_analysis(module = "power_lmm", {
        mouseExperiment::apriori_power_simulation(
          n_per_group         = n_per_grp_lmm,
          n_groups            = as.integer(input$pw_lmm_n_groups %||% 2),
          baseline_volume     = input$pw_lmm_baseline_mean %||% 100,
          baseline_sd         = input$pw_lmm_baseline_sd   %||% 20,
          control_growth_rate = input$pw_lmm_ctrl_rate     %||% 15,
          treatment_effect    = input$pw_lmm_tx_effect     %||% 10,
          random_slope_sd     = input$pw_lmm_rand_slope_sd %||% 3,
          random_intercept_sd = input$pw_lmm_rand_int_sd   %||% 10,
          residual_sd         = input$pw_lmm_resid_sd      %||% 10,
          timepoints          = timepoints,
          alpha               = alpha_vals,
          n_simulations       = as.integer(input$pw_lmm_n_sims %||% 500),
          progress_fn         = progress_fn_shiny
        )
      })

      progress$set(value = 1)

      if (isTRUE(result_lmm$error)) {
        shinyjs::runjs(sprintf("window.trackAnalysisComplete('power', false, '%s');",
                               gsub("'", "\\\\'", result_lmm$message)))
        output$pw_status_msg <- renderUI({
          status_alert(paste("Error:", result_lmm$message), "danger")
        })
        return()
      }

      rv$ap_results <- list(
        scenario_table    = result_lmm$power_table,
        power_curve_data  = result_lmm$power_curve_data,
        sensitivity_table = NULL,
        mode              = "apriori_lmm",
        ap_submode        = "find_power",
        lmm_params        = result_lmm$params
      )

      shinyjs::runjs("window.trackAnalysisComplete('power', true);")
      output$pw_status_msg <- renderUI({
        status_alert("LMM simulation completed.", "success")
      })
    }
  })
  
  output$pw_summary <- renderText({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error)) return(paste("Error:", rv$pw_results$message))
    
    res <- rv$pw_results
    paste(capture.output({
      cat("=== POST-HOC POWER ANALYSIS SUMMARY ===\n\n")
      if (!is.null(res$summary)) {
        if (is.character(res$summary)) {
          cat(res$summary)
        } else {
          method_label <- switch(res$summary$method,
            parametric  = "Parametric (t-test based)",
            simulation  = "Simulation",
            auc         = "AUC-based",
            res$summary$method
          )
          cat("Method:        ", method_label, "\n")
          cat("Alpha levels:  ", paste(res$summary$alpha_levels,  collapse = ", "), "\n")
          cat("Power targets: ", paste(res$summary$power_targets, collapse = ", "), "\n")
          if (!is.null(res$summary$sample_sizes) && length(res$summary$sample_sizes) > 0) {
            cat("\nSample sizes per group:\n")
            for (nm in names(res$summary$sample_sizes)) {
              cat("  ", nm, ": ", res$summary$sample_sizes[[nm]], "\n", sep = "")
            }
          }
          cat("\n")
          if (!is.null(res$summary$overall_power)) {
            cat("Overall observed power:", round(res$summary$overall_power, 4), "\n\n")
          }
          if (!is.null(res$summary$recommendations)) {
            cat("Recommendations:\n")
            for (r in res$summary$recommendations) cat(" -", r, "\n")
          }
        }
      } else {
        cat("Power analysis completed. See individual tabs for detailed results.\n")
      }
    }), collapse = "\n")
  })
  
  output$pw_effect_sizes_table <- renderDT({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error) || is.null(rv$pw_results$effect_sizes)) {
      return(datatable(data.frame(Message = "Effect sizes not available")))
    }
    df <- as.data.frame(rv$pw_results$effect_sizes)
    df <- format_df_for_display(df)
    datatable(df, options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
  })
  
  output$pw_power_table <- renderDT({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error) || is.null(rv$pw_results$post_power_analysis)) {
      return(datatable(data.frame(Message = "Power analysis table not available")))
    }
    df <- as.data.frame(rv$pw_results$post_power_analysis)
    df <- format_df_for_display(df)
    datatable(df, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })
  
  output$pw_sample_sizes_table <- renderDT({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error) || is.null(rv$pw_results$sample_size_estimates)) {
      return(datatable(data.frame(Message = "Sample size estimates not available")))
    }
    df <- as.data.frame(rv$pw_results$sample_size_estimates)
    df <- format_df_for_display(df)
    datatable(df, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  })
  
  output$pw_plot_curves <- renderPlot({
    req(rv$pw_results)
    if (isTRUE(rv$pw_results$error)) return(NULL)

    # Build power-vs-N curves from observed effect sizes
    es_df <- rv$pw_results$effect_sizes
    if (!is.null(es_df) && has_col(es_df, "Standardized_Effect")) {
      unique_effects <- sort(unique(round(
        abs(es_df$Standardized_Effect[!is.na(es_df$Standardized_Effect)]), 2
      )))
      unique_effects <- unique_effects[unique_effects > 0]
    } else {
      unique_effects <- c(0.2, 0.5, 0.8, 1.0, 1.5)
    }

    # Use the smallest user-selected alpha for the power curve visualization
    alpha_for_curve <- min(as.numeric(input$pw_alpha %||% 0.05))

    if (length(unique_effects) > 0) {
      n_seq <- seq(2, 60, by = 1)
      rows <- lapply(unique_effects, function(es) {
        data.frame(
          N            = n_seq,
          Power        = vapply(n_seq, function(n)
            tryCatch(
              stats::power.t.test(n = n, delta = es, sd = 1, sig.level = alpha_for_curve)$power,
              error = function(e) NA_real_
            ), numeric(1)),
          Effect_Size  = paste0("d\u202f=\u202f", round(es, 2)),
          stringsAsFactors = FALSE
        )
      })
      pw_df <- do.call(rbind, rows)

      p <- ggplot2::ggplot(pw_df,
          ggplot2::aes(x = N, y = Power,
                       color = Effect_Size, group = Effect_Size)) +
          ggplot2::geom_line(linewidth = 0.9) +
          ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed",
                              colour = "grey40", linewidth = 0.6) +
          ggplot2::geom_hline(yintercept = 0.90, linetype = "dashed",
                              colour = "grey60", linewidth = 0.6) +
          ggplot2::annotate("text", x = 60, y = 0.81, label = "0.80",
                            hjust = 1, size = 3, colour = "grey40") +
          ggplot2::annotate("text", x = 60, y = 0.91, label = "0.90",
                            hjust = 1, size = 3, colour = "grey60") +
          ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
          ggplot2::labs(
            title  = "Power vs N per Group",
            subtitle = paste0("\u03b1 = ", alpha_for_curve,
                              ";  dashed lines: 0.80 and 0.90 power thresholds"),
            x      = "N per Group",
            y      = "Power",
            colour = "Effect Size\n(Cohen\u2019s d)"
          ) +
          ggplot2::theme_classic() +
          ggplot2::theme(legend.position = "right")
      current_pw_plot(p)
      p
    } else {
      plots <- rv$pw_results$plots
      if (!is.null(plots$power_curves)) {
        p <- plots$power_curves
        current_pw_plot(p)
        p
      } else {
        plot(1, 1, type = "n", main = "Power curves not available", xlab = "", ylab = "")
      }
    }
  })
  
  output$download_pw_effect_sizes <- downloadHandler(
    filename = "power_effect_sizes.csv",
    content = function(file) {
      if (!is.null(rv$pw_results$effect_sizes))
        write.csv(as.data.frame(rv$pw_results$effect_sizes), file, row.names = FALSE)
    }
  )
  
  output$download_pw_power_table <- downloadHandler(
    filename = "power_analysis.csv",
    content = function(file) {
      if (!is.null(rv$pw_results$post_power_analysis))
        write.csv(as.data.frame(rv$pw_results$post_power_analysis), file, row.names = FALSE)
    }
  )
  
  output$download_pw_sample_sizes <- downloadHandler(
    filename = "sample_size_estimates.csv",
    content = function(file) {
      if (!is.null(rv$pw_results$sample_size_estimates))
        write.csv(as.data.frame(rv$pw_results$sample_size_estimates), file, row.names = FALSE)
    }
  )
  
  output$download_pw_plot <- downloadHandler(
    filename = "power_curves.png",
    content = function(file) {
      p <- current_pw_plot()
      req(!is.null(p))
      png(file, width = 900, height = 600)
      print(p); dev.off()
    }
  )

  # ===========================================================================
  # A PRIORI OUTPUTS
  # ===========================================================================

  output$pw_ap_scenario_table <- renderDT({
    req(rv$ap_results)
    df <- rv$ap_results$scenario_table
    if (is.null(df) || nrow(df) == 0)
      return(datatable(data.frame(Message = "No scenario results available")))
    df <- format_df_for_display(as.data.frame(df))
    datatable(df, options = list(scrollX = TRUE, pageLength = 25), rownames = FALSE)
  })

  output$pw_ap_plot_curves <- renderPlot({
    req(rv$ap_results)
    curve_df <- rv$ap_results$power_curve_data
    if (is.null(curve_df) || nrow(curve_df) == 0) {
      plot(1, 1, type = "n", main = "A priori power curves not available", xlab = "", ylab = "")
      return(invisible(NULL))
    }

    mode_ap <- rv$ap_results$mode %||% "apriori_analytic"

    if (mode_ap == "apriori_lmm") {
      # LMM: x = N_Per_Group, y = Power, facet by Alpha
      p <- ggplot2::ggplot(curve_df,
              ggplot2::aes(x = N_Per_Group, y = Power,
                           colour = factor(Alpha), group = factor(Alpha))) +
            ggplot2::geom_line(linewidth = 0.9) +
            ggplot2::geom_point(size = 2) +
            ggplot2::geom_ribbon(
              ggplot2::aes(ymin = pmax(0, Power - 1.96 * SE_Power),
                           ymax = pmin(1, Power + 1.96 * SE_Power),
                           fill = factor(Alpha)),
              alpha = 0.15, colour = NA) +
            ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed",
                                colour = "grey40", linewidth = 0.5) +
            ggplot2::geom_hline(yintercept = 0.90, linetype = "dashed",
                                colour = "grey60", linewidth = 0.5) +
            ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
            ggplot2::labs(
              title    = "A Priori Power \u2014 LMM Simulation",
              subtitle = "Shaded band: \u00b11.96 \u00d7 SE (Monte Carlo uncertainty)",
              x        = "N per Group",
              y        = "Estimated Power",
              colour   = "Alpha",
              fill     = "Alpha"
            ) +
            ggplot2::theme_classic() +
            ggplot2::theme(legend.position = "right")
    } else {
      # Analytic: x = N_Per_Group, y = Power, colour by Effect_Size, facet by Alpha
      curve_df$Comparison <- paste0("d = ", round(curve_df$Effect_Size, 2),
                                    ",  \u03b1 = ", curve_df$Alpha)
      p <- ggplot2::ggplot(curve_df,
              ggplot2::aes(x = N_Per_Group, y = Power,
                           colour = factor(round(Effect_Size, 2)),
                           linetype = factor(Alpha),
                           group = interaction(Effect_Size, Alpha))) +
            ggplot2::geom_line(linewidth = 0.9) +
            ggplot2::geom_hline(yintercept = 0.80, linetype = "dotted",
                                colour = "grey40", linewidth = 0.5) +
            ggplot2::geom_hline(yintercept = 0.90, linetype = "dotted",
                                colour = "grey60", linewidth = 0.5) +
            ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
            ggplot2::labs(
              title    = "A Priori Power Curves",
              subtitle = "Dotted lines: 0.80 and 0.90 power thresholds",
              x        = "N per Group",
              y        = "Power (1 \u2212 \u03b2)",
              colour   = "Cohen\u2019s d",
              linetype = "Alpha (\u03b1)"
            ) +
            ggplot2::theme_classic() +
            ggplot2::theme(legend.position = "right")
    }

    current_ap_plot(p)
    p
  })

  output$pw_ap_sensitivity_table <- renderDT({
    req(rv$ap_results)
    df <- rv$ap_results$sensitivity_table
    if (is.null(df) || nrow(df) == 0)
      return(datatable(data.frame(
        Message = "SD sensitivity analysis is only available in the Mean difference + Pooled SD input mode."
      )))
    df <- format_df_for_display(as.data.frame(df))
    datatable(df, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
  })

  output$download_pw_ap_scenario <- downloadHandler(
    filename = function() {
      if (!is.null(rv$ap_results) && rv$ap_results$mode == "apriori_lmm")
        "lmm_power_table.csv" else "apriori_scenario_table.csv"
    },
    content = function(file) {
      df <- rv$ap_results$scenario_table
      if (!is.null(df)) write.csv(as.data.frame(df), file, row.names = FALSE)
    }
  )

  output$download_pw_ap_plot <- downloadHandler(
    filename = "apriori_power_curves.png",
    content = function(file) {
      p <- current_ap_plot()
      req(!is.null(p))
      png(file, width = 900, height = 600)
      print(p); dev.off()
    }
  )

  output$download_pw_ap_sensitivity <- downloadHandler(
    filename = "apriori_sd_sensitivity.csv",
    content = function(file) {
      df <- rv$ap_results$sensitivity_table
      if (!is.null(df)) write.csv(as.data.frame(df), file, row.names = FALSE)
    }
  )

  })
}
