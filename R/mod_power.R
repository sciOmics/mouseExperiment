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

            # R20.D17 -- the analytic call passed none of these, so the backend's
            # multiplicity and attrition handling (R3.15, K.15) was unreachable.
            radioButtons(ns("pw_ap_p_adjust"), "Multiplicity correction",
                         choices = c("Bonferroni over the comparisons (recommended)" = "bonferroni",
                                     "None (unadjusted)" = "none"),
                         selected = "bonferroni"),
            numericInput(ns("pw_ap_n_comparisons"), "Treated-vs-control comparisons",
                         value = 1, min = 1, max = 20, step = 1),
            numericInput(ns("pw_ap_dropout"), "Expected attrition (%)",
                         value = 0, min = 0, max = 90, step = 5),
            helpText("Power is for each treated-vs-control comparison. A study with ",
                     "k groups makes k \u2212 1 of them, and Bonferroni divides alpha ",
                     "among them. Dunnett, the Tumor Growth default, is slightly less ",
                     "conservative, so the N here is a small overestimate for it. ",
                     "Attrition is the share of animals you expect to lose before the ",
                     "endpoint; Enroll_N is the number to start with."),

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
            p(tags$small(tags$em(
              "Baseline SD does not change the estimated power. The test is on ",
              "Treatment × Day — a contrast of growth ", tags$em("rates"), " — and ",
              "the per-mouse random intercept absorbs baseline spread exactly. ",
              "It is here because it shapes the simulated volumes, not the ",
              "sample-size answer."))),
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
            numericInput(ns("pw_lmm_dropout_limit"),
                         "Euthanasia volume limit (blank = none)",
                         value = NA, min = 1, step = 50),
            p(tags$small(tags$em(
              "Set your IACUC limit to plan against a study that loses animals. ",
              "Leaving it blank assumes every animal is measured at every ",
              "timepoint, which overstates power — by roughly 6 points at 30% ",
              "attrition and 12 at 40%."))),

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

          tabPanel("Summary",
            fluidRow(
              column(10, offset = 1,
                br(),
                shinycssloaders::withSpinner(
                  uiOutput(ns("pw_ap_summary_ui"))
                )
              )
            )
          ),

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
                      tags$li(strong("A priori — Analytic"), " — closed-form power for each treated-vs-control comparison (two-sample t-test), with the alpha split across the k \u2212 1 comparisons a k-group study makes. Fast; no simulation required."),
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
                    p("Power is computed for one treated-vs-control comparison: a two-sample t-test (non-central t-distribution) with Cohen's d = \u0394 / SD_pooled. That is how these studies are analysed, as k \u2212 1 comparisons against the control arm."),
                    tags$blockquote(
                      p(tags$b("Two groups:"), " one comparison at the chosen alpha."),
                      p(tags$b("Three or more groups:"), " each comparison is tested at alpha / (number of comparisons) under the default Bonferroni correction, so each treated arm that differs from control by d is detected with the reported power. ",
                        "Dunnett, the Tumor Growth default, is slightly less conservative, so the N shown is a small overestimate for it. ",
                        "Versions before 26.09.002 used a one-way ANOVA with f = d / \u221a2 here, which gave sample sizes 2.6\u20134\u00d7 too small (CODE_REVIEW.md R20.7).")
                    ),
                    p(tags$b("Attrition:"), " Required_N is the number that must reach the endpoint; with an expected attrition above 0 %, Enroll_N is the number to start with."),

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
                      " is fitted as in the Tumor Growth analysis (REML; uncorrelated slopes, then a random intercept, when a fit fails). ",
                      "Each treated arm's growth rate is compared with the control's by a Satterthwaite t-test at \u03b1/(k\u22121). ",
                      "Power = proportion of simulations in which a given arm is significant; the Scenario Table also gives the chance that ",
                      "any arm, or every arm, is significant."),

                    hr(),
                    h5(strong("Key Parameters")),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Control log-growth rate"),
                      tags$dd("Exponential growth constant for the control group on the log scale. A rate of 0.15 corresponds to ~16 % daily volume increase."),
                      tags$dt("Treatment rate reduction"),
                      tags$dd("Reduction in growth rate for treated vs control, on the log scale."),
                      tags$dt("Random intercept SD"),
                      tags$dd("Between-animal variability in baseline volume (log scale). ",
                              tags$strong("Does not affect the power estimate"),
                              " — the test is on Treatment × Day, a contrast of growth rates, and the per-animal random intercept absorbs baseline spread exactly. Varying it from 0.05 to 1.50 leaves power unchanged. The same is true of Baseline SD, with which it is aliased."),
                      tags$dt("Random slope SD"),
                      tags$dd("Between-animal variability in growth ", tags$em("rate"), " (log scale). ",
                              tags$strong("This is the dominant driver of power"),
                              " — it is the noise the treatment effect must be detected against. In a worked 2-arm, n = 8 design, moving it from 0.01 to 0.15 dropped power from 1.00 to 0.15. Estimate it from your own historical control data before trusting a sample-size answer."),
                      tags$dt("Residual SD"),
                      tags$dd("Within-animal measurement noise (log scale)."),
                      tags$dt("Euthanasia volume limit"),
                      tags$dd("Removes an animal's observations once its simulated volume crosses the limit, mimicking an IACUC removal. Leaving it blank assumes complete data and overstates power — about 6 points at 30 % attrition, 12 at 40 %. Realised attrition is reported with the results."),
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
#' What the LMM simulation's power refers to
#'
#' Backend 0.30.0 (R20.33): `Power` is per treated arm against the control,
#' tested at alpha/(k-1). It used to be an omnibus likelihood-ratio test that
#' rejected 9 % of the time under the null at n = 5.
#'
#' @param st The simulation's `power_table`.
#' @return A `tags$p`, or NULL.
#' @keywords internal
lmm_power_note <- function(st) {
  if (is.null(st) || !nrow(st) || is.null(st$Alpha_Per_Comparison)) return(NULL)
  k <- st$N_Groups[1]
  failed <- sum(st$N_Failed, na.rm = TRUE)
  tags$p(
    tags$strong("Power:"), " the chance that a given treated arm is significant against ",
    "the control, using the tumour-growth analysis's test of growth rates (Satterthwaite t)",
    if (isTRUE(k > 2)) sprintf(" at α/(k−1) = %s for %d treated arms",
                               format(signif(st$Alpha_Per_Comparison[1], 3)), k - 1L),
    ". In the Scenario Table, ", tags$em("Power_Any"), " is the chance at least one arm ",
    "is significant, ", tags$em("Power_All"), " that every arm is, and ",
    tags$em("Power_Omnibus"), " that the joint test of all arms is.",
    if (failed > 0) sprintf(" %d simulated stud%s could not be fitted and are excluded.",
                            failed, if (failed == 1) "y" else "ies")
  )
}

power_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {

  current_ap_plot    <- reactiveVal(NULL)

  # D.7: this was an unconditional `rv$ap_results <- NULL`. `rv` is created empty
  # in app.R and each module mounts once, so the assignment never cleared anything
  # in production -- but it wiped the slot at module init, which made these five
  # modules impossible to drive through testServer() with a pre-populated result.
  # That is why I.2 stayed open across four rounds: every attempt to assert render
  # behaviour got the "run an analysis first" placeholder back. Removed; the slot
  # reads NULL until an analysis writes it either way.

  # ---- Seed pooled SD from loaded data (analytic a priori) ------------------
  # R20.D17: the comparison family follows the number of groups (k - 1
  # treated-vs-control comparisons) unless the user sets it afterwards.
  observeEvent(input$pw_ap_n_groups, {
    k <- as.integer(input$pw_ap_n_groups %||% 2)
    updateNumericInput(session, "pw_ap_n_comparisons", value = max(1L, k - 1L))
  })

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
      # R20.D17: plan against the analysis that will be run -- k - 1
      # vs-control comparisons with a correction -- and against attrition.
      p_adj   <- input$pw_ap_p_adjust %||% "bonferroni"
      n_comp  <- max(1L, as.integer(input$pw_ap_n_comparisons %||% (n_groups - 1L)))
      dropout <- min(max(as.numeric(input$pw_ap_dropout %||% 0), 0), 90) / 100

      show_analysis_modal("Computing a priori power…", "Results will appear automatically when complete.")
      on.exit(shiny::removeModal(), add = TRUE)

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
              p_adjust_method      = p_adj,
              n_comparisons        = n_comp,
              dropout_rate         = dropout,
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
              p_adjust_method      = p_adj,
              n_comparisons        = n_comp,
              dropout_rate         = dropout,
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
        ap_submode        = input$pw_ap_submode,
        method_note       = ap_runs[[1]]$method_note,
        p_adjust_method   = p_adj,
        n_comparisons     = n_comp,
        dropout_rate      = dropout
      )

      shinyjs::runjs("window.trackAnalysisComplete('power', true);")
      output$pw_status_msg <- renderUI({
        status_alert("A priori power analysis completed.", "success")
      })
      updateTabsetPanel(session, "pw_ap_results_tabs", selected = "Summary")

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

      show_analysis_modal("Running LMM power simulation…", "This may take several minutes depending on simulation count.")
      on.exit(shiny::removeModal(), add = TRUE)

      result_lmm <- safe_analysis(module = "power_lmm", {
        mouseExperiment::apriori_power_simulation(
          n_per_group         = n_per_grp_lmm,
          n_groups            = as.integer(input$pw_lmm_n_groups %||% 2),
          baseline_volume     = input$pw_lmm_baseline_mean %||% 100,
          baseline_sd         = input$pw_lmm_baseline_sd   %||% 20,
          # These fallbacks must stay on the log scale. They previously read
          # 15 / 10 / 3 / 10 / 10 -- values from the percent-scale parameterisation
          # the backend abandoned when the data-generating model moved to the log
          # scale. Had one fired, control_growth_rate = 15/day would have made
          # exp(15 * t) overflow to Inf and the fit fail with an opaque error.
          control_growth_rate = input$pw_lmm_ctrl_rate     %||% 0.15,
          treatment_effect    = input$pw_lmm_tx_effect     %||% 0.10,
          random_slope_sd     = input$pw_lmm_rand_slope_sd %||% 0.05,
          random_intercept_sd = input$pw_lmm_rand_int_sd   %||% 0.20,
          residual_sd         = input$pw_lmm_resid_sd      %||% 0.10,
          timepoints          = timepoints,
          dropout_limit       = {
            dl <- suppressWarnings(as.numeric(input$pw_lmm_dropout_limit))
            if (length(dl) != 1L || is.na(dl) || dl <= 0) Inf else dl
          },
          alpha               = alpha_vals,
          n_simulations       = as.integer(input$pw_lmm_n_sims %||% 500)
        )
      })

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
        lmm_params        = result_lmm$params,
        lmm_attrition     = result_lmm$attrition
      )

      shinyjs::runjs("window.trackAnalysisComplete('power', true);")
      output$pw_status_msg <- renderUI({
        status_alert("LMM simulation completed.", "success")
      })
      updateTabsetPanel(session, "pw_ap_results_tabs", selected = "Summary")
    }
  })

  # ===========================================================================
  # A PRIORI OUTPUTS
  # ===========================================================================

  # ── Summary tab — concise overview of the most recent run ───────────────
  # CODE_REVIEW.md Round 2 E.1: every other module has a Summary tab, and
  # the auto-switch-to-Summary behaviour added in v26.05.022 was skipping
  # Power Analysis because it had no such tab. Adding it both restores
  # parity and gives users a one-glance "what did this analysis answer".
  output$pw_ap_summary_ui <- renderUI({
    if (is.null(rv$ap_results)) {
      return(div(
        class = "alert alert-info",
        icon("info-circle"), " ",
        "Run an a priori power analysis to see results here."
      ))
    }
    res     <- rv$ap_results
    mode_ap <- res$mode %||% "apriori_analytic"
    submode <- res$ap_submode %||% NA_character_

    # Recommended N: smallest N reaching the target power across scenarios.
    target  <- isolate(input$pw_ap_target_power) %||% 0.80
    rec_n_html <- if (mode_ap == "apriori_lmm") {
      cd <- res$power_curve_data
      ok <- cd[!is.na(cd$Power) & cd$Power >= target, , drop = FALSE]
      if (nrow(ok) > 0L) {
        min_n <- min(ok$N_Per_Group, na.rm = TRUE)
        tags$p(
          tags$strong("Recommended N per group:"), " ",
          tags$span(class = "badge", min_n),
          tags$br(),
          tags$small(class = "text-muted",
            sprintf("Smallest simulated N reaching power ≥ %.0f%%.", 100 * target))
        )
      } else {
        tags$p(
          tags$strong("Recommended N per group:"), " not reached",
          tags$br(),
          tags$small(class = "text-muted",
            sprintf("None of the simulated N values reached power ≥ %.0f%%; extend the range.", 100 * target))
        )
      }
    } else {
      st <- res$scenario_table
      n_col <- intersect(c("Required_N", "N_Required"), colnames(st))[1]
      if (!is.na(n_col) && any(is.finite(st[[n_col]]))) {
        min_n <- min(st[[n_col]][is.finite(st[[n_col]])], na.rm = TRUE)
        tags$p(
          tags$strong("Recommended N per group (smallest scenario):"), " ",
          tags$span(class = "badge", min_n),
          tags$br(),
          tags$small(class = "text-muted",
            sprintf("From the scenario with the most favourable (effect, alpha) combination at target power = %.0f%%.", 100 * target))
        )
      } else {
        tags$p(
          tags$strong("Recommended N per group:"), " see Scenario Table"
        )
      }
    }

    method_str <- switch(mode_ap,
      "apriori_lmm"      = "A priori — LMM simulation (longitudinal mixed model)",
      "apriori_analytic" = sprintf(
        "A priori — analytic (%s)",
        switch(submode %||% "",
          "find_n"     = "find N for target power",
          "find_power" = "find power for fixed N",
          "analytic"
        )
      ),
      paste0("A priori — ", mode_ap)
    )

    # What the reported power was actually planned against. A power figure from a
    # complete-data simulation is an upper bound for a study that loses animals,
    # so the assumption is stated rather than left implicit.
    attrition_html <- local({
      att <- res$lmm_attrition
      if (is.null(att) || !nrow(att)) return(NULL)
      lost <- max(att$Pct_Obs_Lost, na.rm = TRUE)
      if (!is.finite(lost) || lost <= 0) {
        tags$p(tags$span(class = "text-muted",
          tags$strong("Dropout:"), " none modelled — every animal contributes ",
          "every timepoint. This is an upper bound on achievable power; a study ",
          "losing 30 % of its observations runs roughly 6 points lower. Set an ",
          "euthanasia volume limit to plan against your actual attrition."))
      } else {
        tags$p(
          tags$strong("Dropout modelled:"), " ",
          sprintf("limit %s — %.1f %% of observations lost, %.1f %% of animals truncated.",
                  format(att$Dropout_Limit[1]), lost,
                  max(att$Pct_Mice_Truncated, na.rm = TRUE))
        )
      }
    })

    div(
      div(class = "panel panel-default",
        div(class = "panel-heading",
            h4(class = "panel-title", "Power-Analysis Summary")),
        div(class = "panel-body",
          tags$p(tags$strong("Method:"), " ", method_str),
          tags$p(tags$strong("Scenarios:"), " ",
                 nrow(res$scenario_table), " row(s) in Scenario Table"),
          rec_n_html,
          attrition_html,
          # Backend 0.30.0 (R20.33): what the simulated power refers to.
          if (identical(mode_ap, "apriori_lmm")) lmm_power_note(res$scenario_table),
          # R20.D17: say what the analytic power refers to and what it assumes.
          if (identical(mode_ap, "apriori_analytic")) tags$p(
            tags$strong("Comparisons:"), " ",
            if (identical(res$p_adjust_method, "none")) {
              "no multiplicity correction"
            } else {
              sprintf("Bonferroni over %d treated-vs-control comparison(s)",
                      res$n_comparisons %||% 1L)
            },
            ". ", tags$strong("Attrition:"), " ",
            sprintf("%.0f %%", 100 * (res$dropout_rate %||% 0)),
            if ((res$dropout_rate %||% 0) > 0) " (Enroll_N in the Scenario Table is the number to start with)." else
              " \u2014 none assumed; Required_N is the number that must reach the endpoint."),
          if (!is.null(res$method_note)) tags$p(tags$small(class = "text-muted", res$method_note)),
          tags$hr(),
          tags$p(tags$small(class = "text-muted",
            "Full scenario grid is in the ", tags$strong("Scenario Table"),
            " tab; power curves are in ", tags$strong("Power Curves"),
            "; SD sensitivity (analytic mode only) is in ",
            tags$strong("SD Sensitivity"), "."))
        )
      )
    )
  })

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
            # Backend 0.30.0 (R20.33): Wilson limits, which stay inside
            # [0, 1] and do not collapse to zero width at 0 % or 100 %.
            ggplot2::geom_ribbon(
              ggplot2::aes(ymin = Power_Lower, ymax = Power_Upper,
                           fill = factor(Alpha)),
              alpha = 0.15, colour = NA) +
            ggplot2::geom_hline(yintercept = 0.80, linetype = "dashed",
                                colour = "grey40", linewidth = 0.5) +
            ggplot2::geom_hline(yintercept = 0.90, linetype = "dashed",
                                colour = "grey60", linewidth = 0.5) +
            ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
            ggplot2::labs(
              title    = "A Priori Power — LMM Simulation",
              subtitle = "Power for one treated arm vs control. Shaded band: Wilson 95% interval (simulation uncertainty)",
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
