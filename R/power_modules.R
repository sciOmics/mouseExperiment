#' Power Analysis UI
#' @keywords internal
powerAnalysisUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Power Analysis"),

    p("Estimate statistical power and sample sizes using a priori methods.",
      " Configure parameters in the sidebar and click Run Analysis."),

    fluidRow(
      # ---- Sidebar -----------------------------------------------------------
      column(3,
        wellPanel(

          # ---- Mode selector -------------------------------------------------
          radioButtons(ns("pw_mode"),
                       label    = "Analysis Mode",
                       choices  = c("A priori — Analytic"      = "apriori_analytic",
                                    "A priori — LMM Simulation" = "apriori_lmm"),
                       selected = "apriori_analytic"),

          hr(),

          # ---- A PRIORI ANALYTIC controls ------------------------------------
          conditionalPanel(ns = ns, condition = "input.pw_mode == 'apriori_analytic'",

            h4("Effect Size"),

            radioButtons(ns("pw_ap_es_mode"),
                         label    = NULL,
                         choices  = c("Cohen’s d directly"          = "cohens_d",
                                      "Mean difference + Pooled SD"       = "delta_sd"),
                         selected = "cohens_d"),

            conditionalPanel(ns = ns, condition = "input.pw_ap_es_mode == 'cohens_d'",
              textInput(ns("pw_ap_cohens_d"),
                        "Cohen’s d values (comma-separated)",
                        value       = "0.5, 0.8, 1.0, 1.5",
                        placeholder = "e.g. 0.5, 0.8, 1.0"),
              p(tags$small("Benchmarks: small ≈ 0.2, medium ≈ 0.5, large ≈ 0.8"))
            ),

            conditionalPanel(ns = ns, condition = "input.pw_ap_es_mode == 'delta_sd'",
              numericInput(ns("pw_ap_delta"),
                           "Expected mean difference (Δ)",
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
                        placeholder = "e.g. 5, 8, 10")
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
                         "Control log-growth rate (/day)",
                         value = 0.15, min = 0, step = 0.01),
            numericInput(ns("pw_lmm_tx_effect"),
                         "Treatment rate reduction (/day, log-scale)",
                         value = 0.10, min = 0, step = 0.01),
            p(tags$small("All growth rates are on the log scale. Rate 0.15 ≈ 16% daily volume increase.")),

            hr(),

            h4("Variance Components (log scale)"),

            numericInput(ns("pw_lmm_rand_int_sd"),
                         "Random intercept SD (log-scale)",
                         value = 0.20, min = 0, step = 0.01),
            numericInput(ns("pw_lmm_rand_slope_sd"),
                         "Random slope SD (log-scale)",
                         value = 0.05, min = 0, step = 0.01),
            numericInput(ns("pw_lmm_resid_sd"),
                         "Residual SD (log-scale)",
                         value = 0.10, min = 0.01, step = 0.01),

            hr(),

            h4("Study Design"),

            numericInput(ns("pw_lmm_n_groups"),
                         "Number of groups (incl. control)",
                         value = 2, min = 2, max = 8, step = 1),
            textInput(ns("pw_lmm_n_per_group"),
                      "N per group to evaluate (comma-separated)",
                      value       = "5, 8, 10, 12, 15",
                      placeholder = "e.g. 5, 8, 10, 12, 15"),
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
                         " Use ≤500 simulations for quick estimates.")),

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

          # ---- Target Power (analytic a priori only) ------------------------
          conditionalPanel(ns = ns, condition = "input.pw_mode == 'apriori_analytic'",
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
        tabsetPanel(
          id = ns("pw_ap_results_tabs"),

          tabPanel("Info",
            fluidRow(
              column(10, offset = 1,
                br(),

                # ── Overview ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Overview")),
                  div(class = "panel-body",
                    p("The Power Analysis tab supports prospective (a priori) sample-size planning before an experiment is run. Both methods estimate the probability of detecting a true effect of a specified magnitude at a chosen significance level."),
                    tags$ul(
                      tags$li(strong("A priori — Analytic"), " — closed-form power calculation using a two-sample t-test (two groups) or one-way ANOVA (three or more groups). Fast; no simulation required."),
                      tags$li(strong("A priori — LMM Simulation"), " — Monte Carlo simulation of longitudinal tumor-growth data. Fits the same mixed-effects model used in the Tumor Growth tab and extracts power from the Treatment × Day interaction. More realistic for repeated-measures designs.")
                    ),
                    p(tags$small(class = "text-muted",
                      "Post-hoc power analysis (computed from observed data after the experiment) is not provided. ",
                      "Post-hoc power is a monotone function of the p-value and adds no information beyond the test result itself ",
                      "(Hoenig & Heisey, 2001)."
                    ))
                  )
                ),

                # ── Analytic Mode ────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-info", "Analytic Mode"), " — Method & Parameters")),
                  div(class = "panel-body",

                    h5(strong("Method")),
                    p("For two groups, a two-sample t-test power calculation is used (non-central t-distribution). For three or more groups, one-way ANOVA power is computed via the non-central F distribution."),
                    tags$blockquote(
                      p(tags$b("Two groups:"), " Cohen's d = Δ / SD_pooled; power from the non-central t-distribution."),
                      p(tags$b("Three or more groups:"), " Cohen's d is converted to Cohen's f via ", tags$code("f = d / √2"), ". ",
                        "This assumes the two most extreme groups span the full effect and overestimates f — and therefore power — when all treated groups are uniformly shifted vs control. ",
                        "A method note is shown in the Scenario Table output when ≥ 3 groups are selected.")
                    ),

                    hr(),
                    h5(strong("Effect Size Input")),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Cohen's d"),
                      tags$dd("Enter one or more d values directly (comma-separated). Conventional benchmarks: small ≈ 0.2, medium ≈ 0.5, large ≈ 0.8."),
                      tags$dt("Mean difference + SD"),
                      tags$dd("Enter the expected mean difference (Δ) and pooled within-group SD. Cohen's d is computed internally. Use 'Seed SD from loaded data' to populate the SD from your uploaded dataset.")
                    ),

                    hr(),
                    h5(strong("Solve For")),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Find N"),
                      tags$dd("Returns the required N per group to achieve each target power × alpha combination. Results shown in the Scenario Table."),
                      tags$dt("Find Power"),
                      tags$dd("Returns achieved power at specified N values. Useful when sample size is constrained by practical considerations.")
                    ),

                    hr(),
                    h5(strong("SD Sensitivity Tab")),
                    p("Shows how the required N changes when the assumed SD is perturbed by ±20 % and ±40 %. Only available in Mean Difference + Pooled SD mode. Use this to gauge sensitivity of your sample size estimate to uncertainty in the expected variance.")
                  )
                ),

                # ── LMM Simulation Mode ───────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "LMM Simulation Mode"), " — Method & Parameters")),
                  div(class = "panel-body",

                    h5(strong("Method")),
                    p("Synthetic tumor-growth datasets are generated on the log scale:"),
                    tags$blockquote(tags$code(
                      "log(Volume) = log(baseline) + b₀ + (rate + b₁) × t + ε"
                    )),
                    p("where ", tags$code("b₀"), " ~ N(0, random intercept SD), ", tags$code("b₁"), " ~ N(0, random slope SD), and ",
                      tags$code("ε"), " ~ N(0, residual SD). Volumes are exponentiated back to the original scale."),
                    p("For each simulated dataset, ", tags$code("log(Volume) ~ Treatment × Day + (Day | ID)"),
                      " is fitted via lme4 and significance assessed by a likelihood-ratio test (LRT). Power = proportion of simulations where p < α."),

                    hr(),
                    h5(strong("Key Parameters")),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Control log-growth rate"),
                      tags$dd("Exponential growth constant for the control group on the log scale. A rate of 0.15 corresponds to ~16 % daily volume increase."),
                      tags$dt("Treatment rate reduction"),
                      tags$dd("Reduction in growth rate for treated vs control, on the log scale."),
                      tags$dt("Random intercept SD"),
                      tags$dd("Between-animal variability in baseline volume (log scale). Higher values require larger N."),
                      tags$dt("Random slope SD"),
                      tags$dd("Between-animal variability in growth rate (log scale)."),
                      tags$dt("Residual SD"),
                      tags$dd("Within-animal measurement noise (log scale)."),
                      tags$dt("Simulations"),
                      tags$dd("Number of Monte Carlo draws. Use ≥ 1000 for publication-quality estimates; ≤ 500 for quick exploration. Power estimates are shown with ±1.96 × SE uncertainty bands in the Power Curves tab.")
                    )
                  )
                ),

                # ── Interpretation ───────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Interpretation")),
                  div(class = "panel-body",
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Alpha (α)"),
                      tags$dd("Significance threshold (Type I error rate). Conventional choices are 0.05 (5 %) and 0.01 (1 %). Lower alpha requires a larger N to achieve the same power."),
                      tags$dt("Power (1−β)"),
                      tags$dd("Probability of detecting a true effect. 0.80 is a common minimum; 0.90 is recommended when false negatives are costly."),
                      tags$dt("N per group"),
                      tags$dd("Sample size per treatment arm. Total study animals = N per group × number of groups.")
                    ),
                    p(tags$small(class = "text-muted",
                      "Power estimates assume the specified effect size is exactly correct. ",
                      "If the true effect is smaller, actual power will be lower than estimated. ",
                      "Use the SD Sensitivity tab to understand robustness to variance assumptions."
                    ))
                  )
                ),

                br()
              )
            )
          ),

          tabPanel("Scenario Table",
            tab_info(
              strong("A priori power / sample size table: "),
              "For each combination of effect size, alpha, and target power, shows the required N per group (Find N mode) or achieved power (Find Power / LMM mode).",
              tags$ul(
                tags$li(strong("Effect_Size"), " — Cohen’s d (analytic) or the simulated scenario (LMM)."),
                tags$li(strong("N_Groups"), " — number of treatment groups in the design."),
                tags$li(strong("Alpha"), " — significance threshold."),
                tags$li(strong("Required_N"), " — N per group needed to reach the target power (Find N mode)."),
                tags$li(strong("Achieved_Power / Power"), " — power reached at the specified N (Find Power / LMM mode).")
              )
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("pw_ap_scenario_table"))
            ),
            downloadButton(ns("download_pw_ap_scenario"), "Download")
          ),

          tabPanel("Power Curves",
            tab_info(strong("Power vs N per group: "), "prospective power curves for the specified effect size scenario(s). Each line represents one Cohen’s d / alpha combination. Dashed lines mark 0.80 and 0.90 power thresholds."),
            shinycssloaders::withSpinner(
              plotOutput(ns("pw_ap_plot_curves"), height = "500px")
            ),
            downloadButton(ns("download_pw_ap_plot"), "Download Plot")
          ),

          tabPanel("SD Sensitivity",
            tab_info(
              strong("Sensitivity to the assumed SD: "), "shows how the required N changes when the pooled within-group SD used to convert a raw mean difference into Cohen’s d is perturbed by ±20 % and ±40 %. Only available when using the Mean difference + Pooled SD input mode.",
              tags$ul(
                tags$li(strong("SD_Change"), " — percentage change applied to the assumed SD."),
                tags$li(strong("Assumed_SD"), " — resulting SD value."),
                tags$li(strong("Cohens_d"), " — Cohen’s d after adjusting for the new SD."),
                tags$li(strong("Required_N"), " — sample size required at the reference alpha / target power.")
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
}

#' Power Analysis Server Logic
#' @keywords internal
power_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

  current_ap_plot    <- reactiveVal(NULL)

  rv$ap_results <- NULL

  # ---- Seed pooled SD from loaded data (analytic a priori) ------------------
  observeEvent(input$pw_seed_sd, {
    prep <- tryCatch(ensure_preprocessed(rv), error = function(e) NULL)
    if (is.null(prep)) {
      showNotification("No data loaded yet — upload data first.", type = "warning", duration = 4)
      return()
    }
    df      <- prep$data
    tx_col  <- prep$tx_col
    vol_col <- prep$vol_col
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
      status_alert("Running analysis…", "info")
    })

    # =========================================================================
    # BRANCH A — A PRIORI ANALYTIC
    # =========================================================================
    if (input$pw_mode == "apriori_analytic") {

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
            status_alert("Please enter at least one Cohen’s d value.", "danger")
          })
          return()
        }
        parsed <- suppressWarnings(as.numeric(strsplit(raw_txt, "[,\\s]+")[[1]]))
        es_vals <- parsed[!is.na(parsed) & parsed > 0]
        if (length(es_vals) == 0) {
          output$pw_status_msg <- renderUI({
            status_alert("No valid positive Cohen’s d values found.", "danger")
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

      progress <- shiny::Progress$new(session, min = 0, max = 1, style = "old")
      progress$set(message = "Computing a priori power…", value = 0.1)
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
    # BRANCH B — A PRIORI LMM SIMULATION
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

      progress <- shiny::Progress$new(session, min = 0, max = 1, style = "old")
      progress$set(message = "Running LMM simulation… (this may take a while)", value = 0.05)
      on.exit(progress$close(), add = TRUE)

      total_scenarios  <- length(n_per_grp_lmm) * length(alpha_vals)
      completed_scens  <- 0L
      progress_fn_shiny <- function(prop) {
        completed_scens <<- completed_scens + 1L
        progress$set(value = 0.05 + 0.90 * prop,
                     message = sprintf("LMM simulation: scenario %d / %d…",
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
        shinyjs::runjs(sprintf("window.trackAnalysisComplete('power', false, %s);",
                               js_string(result_lmm$message)))
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
              title    = "A Priori Power — LMM Simulation",
              subtitle = "Shaded band: ±1.96 × SE (Monte Carlo uncertainty)",
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
                                    ",  α = ", curve_df$Alpha)
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
              y        = "Power (1 − β)",
              colour   = "Cohen’s d",
              linetype = "Alpha (α)"
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
