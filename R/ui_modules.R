#' UI Modules for mouseExperimentDashboard

# Shared helper: creates a styled info box for tab descriptions
tab_info <- function(...) {
  div(
    class = "alert alert-info",
    style = "font-size:0.88em; margin:10px 0 8px 0; padding:8px 14px;",
    ...
  )
}

# Shared helper: inline help tooltip (question-circle icon with title hint)
help_tip <- function(text) {
  tags$span(
    class = "text-muted",
    style = "cursor:help; margin-left:4px;",
    title = text,
    icon("circle-question", class = "fa-sm")
  )
}

# Helper: collapsible section for the plot control panel.
# Starts collapsed; click the header to expand.
ctrl_section <- function(id, icon_nm, title, ...) {
  collapse_id <- paste0("ctrl-", id)
  tagList(
    tags$a(
      href = paste0("#", collapse_id),
      `data-toggle` = "collapse",
      class = "collapsed",
      style = paste0(
        "display:block; font-size:1.2rem; font-weight:700; ",
        "color:#333; text-decoration:none; ",
        "padding:5px 2px; margin-top:4px; user-select:none; ",
        "cursor:pointer;"
      ),
      icon(icon_nm, style = "margin-right:5px;"),
      title,
      tags$span(
        class = "pull-right text-muted",
        style = "font-size:0.75rem; padding-top:3px;",
        HTML("&#9660;")
      )
    ),
    tags$div(id = collapse_id, class = "collapse", ...),
    hr(style = "margin:5px 0;")
  )
}
#' Home Page UI
#' @keywords internal
homeUI <- function() {
  fluidPage(
    HTML("<div style='padding: 40px;'>"),
    h1("Welcome to mouseExperimentDashboard"),
    p("Analyze mouse tumor growth experiments without writing R code."),
    
    fluidRow(
      column(6,
        h3("Getting Started"),
        p("1. Navigate to", strong("Data Upload"), "and load your tumor growth data"),
        p("2. Configure column mappings if needed"),
        p("3. Select", strong("Tumor Growth Analysis"), " to run statistical models"),
        p("4. Download plots and tables from your analysis")
      ),
      column(6,
        h3("Quick Links"),
        p(actionButton("goto_upload", "Upload Data", icon("upload"), class = "btn-lg btn-primary"), 
          style = "margin-bottom: 20px;"),
        p(actionButton("goto_analysis", "Run Analysis", icon("calculator"), class = "btn-lg btn-success"))
      )
    ),
    
    hr(),
    h4("Features"),
    tags$ul(
      tags$li("Upload CSV or RDA files"),
      tags$li("Automatic column detection"),
      tags$li("Multiple statistical models (LME4, AUC)"),
      tags$li("Interactive visualizations"),
      tags$li("Download results in multiple formats")
    ),
    
    hr(),
    HTML("</div>")
  )
}

#' Data Upload UI
#' @keywords internal
dataUploadUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Data Upload & Configuration"),

    fluidRow(
      column(12,
        wellPanel(
          fluidRow(
            column(4,
              h5(strong("Try Demo Data")),
              p("Load a built-in demo dataset to explore the dashboard's features."),
              selectInput(ns("demo_dataset"), NULL,
                choices = c(
                  "-- Select a demo dataset --" = "",
                  "Combo Treatment (Control / HDACi / aPD1 / HDACi+PD1, survival data)" = "combo",
                  "Dose Levels (dose-response data with survival)" = "dose_levels",
                  "Necrotic Data (None / Drug, length & width columns)" = "necrotic",
                  "Weight Data (DMSO / Drug A / Drug B / Drug A+B, pre-calculated volume)" = "weight"
                ),
                selected = ""
              ),
              fluidRow(
                column(6,
                  actionButton(ns("load_demo_data"), "Load Demo Data",
                               icon = icon("play-circle"),
                               class = "btn-default btn-block")
                ),
                column(6, uiOutput(ns("demo_load_msg")))
              )
            ),
            column(8,
              tags$label("Choose file (csv)"),
              tags$p(tags$small(
                icon("info-circle"), " Data must be in ",
                tags$strong("long format"), " — one row per observation/measurement."
              )),
              fileInput(ns("data_file"), NULL,
                       accept = c(".csv", ".rda", ".RData")),
              uiOutput(ns("file_validation_msg"))
            )
          )
        )
      )
    ),
    
    # Column mapping section (shown after file upload)
    uiOutput(ns("column_mapping_ui")),
    
    # Dataset summary (shown after successful upload + column config)
    uiOutput(ns("data_summary_ui")),
    
    # Data preview (shown after successful upload)
    uiOutput(ns("data_preview_ui")),

    # Session save / restore
    session_save_restore_ui(ns)
  )
}

#' Column Mapping UI
#' @keywords internal
column_mapping_ui <- function(col_names, detected_cols, ns = NS(NULL)) {
  wellPanel(
    h4(strong("Configure Column Mapping")),
    tags$p("Select which columns contain your data. ", 
      "Suggested columns are pre-selected."),
    
    fluidRow(
      column(4,
        h4(strong("Required Columns")),
        radioButtons(ns("time_method"), "How to specify time?",
                    choices = c(
                      "Direct Day column" = "direct_day",
                      "Calculate from Date column" = "from_date"
                    ),
                    selected = if (!is.null(detected_cols$day) && detected_cols$day != "") "direct_day" else "from_date"),
        
        conditionalPanel(ns = ns,
          condition = "input.time_method == 'direct_day'",
          selectInput(ns("col_day"), "Day Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$day)
        ),
        
        conditionalPanel(ns = ns,
          condition = "input.time_method == 'from_date'",
          selectInput(ns("col_date_for_day"), "Date Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$date),
          textInput(ns("reference_date"), "Reference Date (YYYY-MM-DD)",
                   value = "2025-03-24",
                   placeholder = "First measurement date")
        ),

        selectInput(ns("col_id"),
                   tagList("Mouse/Subject ID", help_tip("Unique identifier for each mouse/subject. Each ID should have multiple rows (one per timepoint).")),
                   choices = c("Select..." = "", col_names),
                   selected = detected_cols$id),
        
        selectInput(ns("col_treatment"),
                   tagList("Treatment Group", help_tip("Column indicating which treatment/group each subject belongs to (e.g. Control, Drug A, Drug B).")),
                   choices = c("Select..." = "", col_names),
                   selected = detected_cols$treatment)
      ),
      
      column(4,
        h4(strong("Volume Measurement")),
        radioButtons(ns("volume_method"), "How to specify volume?",
                    choices = c(
                      "Direct volume column" = "direct",
                      "Calculate from Length & Width" = "calculate"
                    ),
                    selected = if (!is.null(detected_cols$volume)) "direct" else "calculate"),
        
        conditionalPanel(ns = ns,
          condition = "input.volume_method == 'direct'",
          selectInput(ns("col_volume"), "Volume Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$volume)
        ),
        
        conditionalPanel(ns = ns,
          condition = "input.volume_method == 'calculate'",
          selectInput(ns("col_length"), "Length Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$length),
          selectInput(ns("col_width"), "Width Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$width),
          selectInput(ns("col_height"), tagList("Height Column",
                     help_tip("Optional third dimension. Required for Ellipsoid (3-axis) and Box formulas. When omitted, height is approximated as equal to width.")),
                     choices = c("None" = "", col_names),
                     selected = ""),
          selectInput(ns("volume_formula"), "Volume Formula",
                     choices = c("Ellipsoid (default)" = "ellipsoid",
                                "Modified Ellipsoid" = "modified_ellipsoid",
                                "Ellipsoid (3-axis, L\u00d7W\u00d7H)" = "ellipsoid_3axis",
                                "Cylinder" = "cylinder",
                                "Sphere" = "sphere",
                                "Box (L\u00d7W\u00d7H)" = "box"),
                     selected = "ellipsoid"),
          conditionalPanel(ns = ns, "input.volume_formula == 'ellipsoid'",
            tags$p(tags$small(HTML("V = (&pi;/6) &times; L &times; W&sup2;")))
          ),
          conditionalPanel(ns = ns, "input.volume_formula == 'modified_ellipsoid'",
            tags$p(tags$small(HTML("V = (L &times; W&sup2;) / 2")))
          ),
          conditionalPanel(ns = ns, "input.volume_formula == 'ellipsoid_3axis'",
            tags$p(tags$small(HTML("V = (&pi;/6) &times; L &times; W &times; H")))
          ),
          conditionalPanel(ns = ns, "input.volume_formula == 'cylinder'",
            tags$p(tags$small(HTML("V = &pi; &times; (W/2)&sup2; &times; L")))
          ),
          conditionalPanel(ns = ns, "input.volume_formula == 'sphere'",
            tags$p(tags$small(HTML("V = (&pi;/6) &times; ((L + W) / 2)&sup3;")))
          ),
          conditionalPanel(ns = ns, "input.volume_formula == 'box'",
            tags$p(tags$small(HTML("V = L &times; W &times; H")))
          )
        )
      ),
      
      column(4,
        h4(strong("Optional Columns")),
        selectInput(ns("col_cage"),
                   tagList("Cage ID", help_tip("Optional grouping variable. Used as a random effect in mixed models to account for cage-to-cage variability.")),
                   choices = c("None" = "", col_names),
                   selected = detected_cols$cage),
        tags$p(tags$small("Cage ID is used to distinguish mice with same ID in different cages.")),
        
        br(),
        
        selectInput(ns("col_dose"),
                   tagList("Dose Level", help_tip("Numeric dose values. Groups are split into Treatment+Dose combinations for dose-response modelling.")),
                   choices = c("None" = "", col_names),
                   selected = detected_cols$dose),
        tags$p(tags$small("Required for Dose-Response Analysis. Should contain numeric dose concentrations or levels.")),
        
        br(),
        selectInput(ns("col_survival_event"),
                   tagList("Event/Censor Status", help_tip("Binary column: 1 = endpoint reached (tumour exceeded IACUC volume limit, or animal euthanised/died), 0 = censored (did not reach endpoint by study end). In most mouse tumour experiments, event = tumour volume threshold, not natural death. Alternatively, leave unset and use the volume-threshold derivation option in the Survival tab sidebar.")),
                   choices = c("None" = "", col_names),
                   selected = detected_cols$survival_event),
        tags$p(tags$small("Required for Survival Analysis (or use volume-threshold derivation in Survival sidebar). Coding: ", tags$strong("1 = endpoint reached (e.g. tumour exceeded IACUC volume limit)"), ", ", tags$strong("0 = censored"), ".")),

        br(),
        selectInput(ns("col_weight"),
                   tagList("Weight / Mass", help_tip("Body weight or mass measurements. Used for toxicity analysis, tumor weight adjustment, and weight-loss survival.")),
                   choices = c("None" = "", col_names),
                   selected = detected_cols$weight),
        tags$p(tags$small("Required for Toxicity Analysis. Should contain body weight or mass measurements.")),

        br(),
        selectInput(ns("col_sex"),
                   tagList("Sex / Gender", help_tip("Optional covariate for mixed models. Accounts for metabolic differences between male and female mice.")),
                   choices = c("None" = "", col_names),
                   selected = detected_cols$sex),
        tags$p(tags$small("Optional covariate for Toxicity mixed models.")),

        br(),
        selectInput(ns("col_necrotic"),
                   tagList(icon("circle-exclamation"), " Necrotic Flag",
                           help_tip("Binary column indicating whether a tumor was visually necrotic at each measurement (e.g. Y/N, yes/no, 1/0). When mapped, enables necrosis-aware adjustment in Tumor Growth analysis.")),
                   choices = c("None" = "", col_names),
                   selected = detected_cols$necrotic %||% ""),
        tags$p(tags$small("Optional. Enables necrotic-observation handling in Tumor Growth analysis."))
      )
    )
  )
}

#' Data Preview UI
#' @keywords internal
data_preview_ui <- function(ns = NS(NULL)) {
  wellPanel(
    h4(strong("Data Preview")),
    p(tags$small("Shows calculated columns (Volume, Day from start date) where applicable.")),
    fluidRow(
      column(12,
        p(strong(textOutput(ns("preview_summary")))),
        DTOutput(ns("data_table")),
        br(),
        downloadButton(ns("download_preview_data"), "Download Processed Data as CSV",
                      class = "btn-sm btn-default")
      )
    )
  )
}

#' Tumor Growth Analysis UI
#' @keywords internal
tumorGrowthUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Tumor Growth Analysis"),
    
    p("Configure analysis parameters and run statistical modeling. ",
      "Ensure data is uploaded first from the Data Upload tab."),
    
    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Analysis Parameters"),

          # 1. Statistical Model
          selectInput(ns("tg_model_type"),
                     "Statistical Model",
                     choices = c("Linear Mixed-Effects (LME4)" = "lme4",
                                "Area Under Curve (AUC)"       = "auc",
                                "Bayesian LMM (brms)"          = "bayes"),
                     selected = "lme4"),
          helpText("LME4 models tumor volume trajectories over time using a linear ",
                   "mixed-effects model, which handles repeated measures and unequal ",
                   "timepoints. AUC summarizes total tumor burden as the area under ",
                   "each animal's growth curve. Bayesian LMM fits the same mixed-effects ",
                   "model as LME4 but reports full posterior distributions and 95 % ",
                   "credible intervals; computation takes 3–12 minutes."),

          # 2. Reference Group
          selectInput(ns("tg_reference"),
                     "Reference Group (Required)",
                     choices = c("-- Select a reference group --" = "")),
          helpText("The group all comparisons are made against (typically a vehicle or ",
                   "IgG control). Automatically set to the most likely control group ",
                   "when data is uploaded."),

          # 3. Post-hoc Comparisons
          selectInput(ns("tg_comparison_type"),
                     "Post-hoc Comparisons",
                     choices = c("vs Reference Group" = "reference",
                                "All Comparisons" = "tukey"),
                     selected = "reference"),
          helpText("Determines which group comparisons are performed after the primary ",
                   "model. 'vs Reference' tests each treatment against the reference ",
                   "group using Dunnett's correction, which is designed for this ",
                   "comparison structure and is less conservative than Bonferroni. ",
                   "'All Comparisons' tests every pair with Tukey HSD familywise error ",
                   "correction."),

          # 4. Data Transformation (LME4 only)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4'",
            selectInput(ns("tg_transform"),
                       "Data Transformation",
                       choices = c("Log" = "log",
                                  "Square Root" = "sqrt",
                                  "None" = "none"),
                       selected = "log"),
            helpText("Transformation applied to tumor volume before modeling. Log ",
                     "transformation is recommended for most tumor growth data as it ",
                     "stabilizes variance and linearizes exponential growth. Use 'None' ",
                     "for already-linear data.")
          ),

          # 5. Extrapolation Points (always visible)
          numericInput(ns("tg_extrapolation_points"),
                      "Extrapolation Points",
                      value = 3, min = 0, max = 10, step = 1),
          helpText("Number of recent timepoints used to estimate the growth slope for ",
                   "extrapolating tumor volumes of mice that died before the study end. ",
                   "Set to 0 to disable. If an Event/Censor column is mapped in Data ",
                   "Upload it will be used; otherwise mice that disappear before the ",
                   "final timepoint are automatically detected."),

          # 6. Cage Effect Handling (LME4 only)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4'",
            selectInput(ns("tg_cage_effect"),
                       "Cage Effect Handling",
                       choices = c("Include if not collinear" = "include_if_not_collinear",
                                  "Always include" = "always_include",
                                  "Never include" = "never_include"),
                       selected = "include_if_not_collinear"),
            helpText("Controls whether cage is included as a covariate to account for ",
                     "housing effects on tumor growth. 'Include if not collinear' adds ",
                     "cage automatically when it is statistically separable from ",
                     "treatment. Requires a Cage ID column to be mapped in Data Upload.")
          ),

          # 7. Random Effects (LME4 and Bayesian)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4' || input.tg_model_type == 'bayes'",
            selectInput(ns("tg_random_effects"),
                       "Random Effects",
                       choices = c("Intercept only" = "intercept_only",
                                  "Intercept + Slope" = "slope",
                                  "None" = "none"),
                       selected = "intercept_only"),
            helpText("Specifies the random effects structure for individual animals. ",
                     "'Intercept only' allows each animal to have its own baseline ",
                     "volume. 'Intercept + Slope' also allows individual growth rates, ",
                     "which is more flexible but may fail to converge with small samples. ",
                     "('None' is not available in Bayesian mode.)")
          ),

          # 8. Bayesian-specific controls
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'bayes'",
            hr(),
            h5(strong("Bayesian Settings")),
            selectInput(ns("tg_prior_strength"), "Prior Strength",
              choices = c(
                "Weakly informative (recommended)" = "weakly_informative",
                "Informative"                      = "informative",
                "Diffuse"                          = "diffuse"
              ),
              selected = "weakly_informative"
            ),
            helpText("Prior on fixed effects (log-scale slopes/contrasts). ",
                     "Weakly informative: N(0, 1) — regularises estimates without ",
                     "strong assumptions. Informative: N(0, 0.5) — tighter. ",
                     "Diffuse: N(0, 2.5) — minimal regularisation."),
            numericInput(ns("tg_bayes_chains"), "MCMC chains",
                         value = 4, min = 1, max = 8, step = 1),
            numericInput(ns("tg_bayes_iter"), "Iterations per chain",
                         value = 2000, min = 500, max = 10000, step = 500),
            div(class = "alert alert-warning", style = "font-size:0.85em; padding:8px; margin-top:6px;",
              icon("clock"), " Bayesian fitting takes ~3–12 minutes. ",
              "The app stays usable while the model runs in the background."
            )
          ),

          # Necrotic handling — only shown when a necrotic column is mapped
          conditionalPanel(
            ns = ns,
            condition = sprintf("input['%s-col_necrotic'] !== '' && input['%s-col_necrotic'] !== undefined",
                                "data", "data"),
            hr(),
            h5(strong("Necrotic Observations")),
            selectInput(ns("tg_necrotic_handling"), NULL,
              choices = c(
                "Exclude from analysis"       = "exclude",
                "Include as covariate"        = "covariate",
                "No adjustment (raw data)"    = "none"
              ),
              selected = "exclude"
            ),
            helpText(
              strong("Exclude"), " — removes flagged observations before fitting the model.",
              tags$br(),
              strong("Covariate"), " — adds a 0/1 Necrotic term to the fixed effects.",
              tags$br(),
              strong("No adjustment"), " — ignores the flag (keeps all observations)."
            )
          ),

          br(),

          actionButton(ns("run_tg_analysis"),
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),

          uiOutput(ns("tg_status_msg"))
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
          id = ns("tg_results_tabs"),

          tabPanel("Info",
            fluidRow(
              column(10, offset = 1,
                br(),

                # ── Overview ────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Overview")),
                  div(class = "panel-body",
                    p("The Tumor Growth tab provides two statistical modeling approaches selectable from the ",
                      strong("Statistical Model"), " dropdown:"),
                    tags$ul(
                      tags$li(strong("Linear Mixed-Effects Model (LME4)"), " — fits a mixed-effects model to repeated tumor volume measurements. Accounts for inter-animal variability and unequal observation times. Recommended for most longitudinal datasets."),
                      tags$li(strong("Area Under Curve (AUC)"), " — summarises each animal’s total tumor burden as the trapezoidal AUC of its growth curve. Non-parametric and robust to missing timepoints, but discards trajectory shape information.")
                    ),
                    p(tags$small(class = "text-muted",
                      "Both modes compute Tumor Growth Inhibition (TGI) for each treatment group relative to the reference. ",
                      "TGI = (1 − V_treated / V_control) × 100 at the final timepoint."
                    ))
                  )
                ),

                # ── LME4 Mode ───────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-info", "LME4 Mode"), " — Analyses & Outputs")),
                  div(class = "panel-body",

                    h5(strong("Model")),
                    p("A linear mixed-effects model is fitted to transformed tumor volumes:"),
                    tags$blockquote(tags$code("log(Volume) ~ Treatment × Day + (1 | ID)")),
                    p("The Treatment × Day interaction captures whether treatment groups differ in their growth ", em("trajectories"), " over time, not just at a single endpoint."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Transformation"), tags$dd("Log (default) stabilises variance and linearises exponential growth. Square root suits moderate heteroscedasticity. Use ‘None’ only for already-linear data."),
                      tags$dt("Random Effects"), tags$dd("‘Intercept only’ allows each animal to have its own baseline volume. ‘Intercept + Slope’ additionally allows individual growth rates — more flexible but may fail to converge with small groups."),
                      tags$dt("Cage Effect"), tags$dd("When a cage ID is mapped in Data Upload, cage can be included as a covariate to account for housing-level clustering. ‘Include if not collinear’ adds it automatically when it is statistically separable from treatment.")
                    ),

                    hr(),
                    h5(strong("ANOVA Tab")),
                    p("Type III ANOVA F-test for Treatment, Day, and the Treatment × Day interaction. A significant interaction (p < 0.05) indicates that treatment groups differ in their growth trajectories — the primary analysis target."),

                    hr(),
                    h5(strong("Comparisons Tab")),
                    p("Pairwise group comparisons at the mean observation day using estimated marginal means (EMMs):"),
                    tags$ul(
                      tags$li(strong("vs Reference (Dunnett)"), " — each treatment vs the reference group only. Dunnett correction controls the familywise error rate for this one-vs-many structure and is less conservative than Bonferroni."),
                      tags$li(strong("All Comparisons (Tukey HSD)"), " — every pairwise combination. Use when all inter-group differences are of scientific interest.")
                    ),
                    p(tags$small(class = "text-muted",
                      "Comparisons are marginalised at the mean study day. For trajectory differences across time, see the Treatment Effects tab."
                    )),

                    hr(),
                    h5(strong("Treatment Effects Tab")),
                    p("Estimated marginal means (EMMs) of transformed volume per treatment group. EMMs are shown at the mean observation day and at five quantile days (min, Q1, median, Q3, max), capturing the Treatment × Day interaction that a single mean-day estimate would miss."),

                    hr(),
                    h5(strong("Growth Rates Tab")),
                    p("Per-animal exponential growth rates estimated by log-linear regression of each animal’s volume trajectory independently. Reported as growth rate (log-scale slope), fold change, and % change per day with 95 % CIs and R². These are descriptive; for group-level inference use the ANOVA and Comparisons tabs."),

                    hr(),
                    h5(strong("Diagnostics Tab")),
                    p("Q–Q plot of residuals and residuals vs fitted values. Systematic curvature suggests the transformation choice should be revisited. Funnel-shaped spread indicates heteroscedasticity.")
                  )
                ),

                # ── AUC Mode ─────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "AUC Mode"), " — Analyses & Outputs")),
                  div(class = "panel-body",
                    p("Trapezoidal AUC of tumor volume over time is computed per animal. AUC (mm³·day) captures total tumor burden integrated over the study, not just the final timepoint."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Individual AUC"), tags$dd("Per-animal AUC values with group summaries (mean, SD, SEM)."),
                      tags$dt("AUC Summary"), tags$dd("Group comparisons via Welch t-tests with Benjamini–Hochberg correction.")
                    ),
                    p(tags$small(class = "text-muted",
                      "LOCF (last observation carried forward) can extrapolate AUC to a common endpoint for animals censored early, ",
                      "so all animals contribute to the same time window. Set Extrapolation Points to 0 to use only observed data."
                    ))
                  )
                ),

                # ── Bayesian LMM Mode ─────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-success", "Bayesian LMM Mode"), " — Method & Outputs")),
                  div(class = "panel-body",

                    h5(strong("Method")),
                    p("Fits the same model formula as LME4 using full Bayesian inference via ",
                      tags$a("brms", href = "https://paul-buerkner.github.io/brms/", target = "_blank"),
                      " / Stan. Instead of point estimates and p-values, every parameter has a ",
                      strong("posterior distribution"), ": a probability distribution over plausible values ",
                      "given the data and priors."),
                    tags$blockquote(tags$code(
                      "log(Volume) ~ Treatment × Day + (1 | ID)"
                    )),
                    p("Results are summarised as posterior medians with ",
                      strong("95 % highest posterior density (HPD) credible intervals"), ". ",
                      "A 95 % CrI has a direct probability interpretation: there is a 95 % posterior ",
                      "probability that the true parameter lies within the interval — unlike a frequentist ",
                      "confidence interval."),

                    hr(),
                    h5(strong("Prior Strength")),
                    p("Priors regularise the model in the absence of strong data. All presets use the ",
                      "same prior family; only the width changes:"),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Weakly informative"),
                      tags$dd("N(0, 1) on fixed effects (log-scale); Exponential(1) on SDs. ",
                              "Recommended for most datasets — prevents extreme estimates without ",
                              "imposing biological assumptions."),
                      tags$dt("Informative"),
                      tags$dd("N(0, 0.5) — tighter regularisation; use when you have strong prior ",
                              "knowledge that effects are small."),
                      tags$dt("Diffuse"),
                      tags$dd("N(0, 2.5) — minimal regularisation; equivalent to near-flat priors. ",
                              "Use with caution on small datasets.")
                    ),

                    hr(),
                    h5(strong("Result Tabs")),
                    tags$ul(
                      tags$li(strong("Summary"), " — posterior median and 95 % CrI for every model parameter."),
                      tags$li(strong("Treatment Effects"), " — group-level estimated marginal means (posterior) at the mean study day."),
                      tags$li(strong("Comparisons"), " — posterior contrasts of each group vs the reference."),
                      tags$li(strong("Posterior Check"), " — overlay of observed data density vs 50 draws from the posterior predictive distribution. Systematic mismatch indicates model mis-specification."),
                      tags$li(strong("MCMC Diagnostics"), " — Rhat and effective sample size (ESS) per parameter. Rhat > 1.01 indicates non-convergence; increase iterations or inspect trace plots."),
                      tags$li(strong("Growth Rates"), " — per-animal log-linear regression rates (same as LME4, not MCMC-based).")
                    ),

                    div(class = "alert alert-info", style = "margin-top:10px;",
                      icon("clock"), " ",
                      strong("Computation time:"), " typically 3–12 minutes on modern hardware. ",
                      "Stan compiles the model on the first call with a given formula and caches the result — ",
                      "subsequent calls with the same formula are faster. ",
                      "The dashboard runs fitting in the background; the app stays responsive while the model runs."
                    )
                  )
                ),

                # ── Necrotic Tumor Handling ───────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title",
                       tags$span(class = "label label-danger", "Necrotic Tumors"),
                       " — Handling Options")),
                  div(class = "panel-body",
                    p("When a ", strong("Necrotic Flag"), " column is mapped in Data Upload, a handling
                      selector appears in the sidebar. Necrotic tumors have a non-viable core that inflates
                      caliper-based volume estimates without contributing to proliferating tissue."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Exclude (recommended)"),
                      tags$dd("Removes flagged observations from the model input. Conservative and publication-defensible.
                               The Necrosis tab reports how many observations were flagged and which groups were affected,
                               regardless of handling mode."),
                      tags$dt("Covariate"),
                      tags$dd("Adds a 0/1 ", code(".necrotic_flag"), " fixed effect to the LMM formula. Estimates the
                               average volume bias attributable to necrosis. Useful when necrosis is approximately
                               evenly distributed across groups and you want to retain all timepoints."),
                      tags$dt("No adjustment"),
                      tags$dd("All observations are used as-is (current default when no necrotic column is mapped).
                               The Necrosis summary table is still shown for documentation purposes.")
                    ),
                    div(class = "alert alert-info", style = "font-size:0.85em; margin-bottom:0;",
                      icon("circle-info"),
                      " The Necrosis tab in the results panel is inserted automatically when a necrotic column
                        is mapped and an analysis has been run."
                    )
                  )
                ),

                # ── Shared Parameters ─────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading", h4(class = "panel-title", "Shared Parameters")),
                  div(class = "panel-body",
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Reference Group"),
                      tags$dd("The vehicle/control group used as comparator. Auto-detected from common names (Vehicle, Control, DMSO, PBS, Saline) or defaults to the first group alphabetically."),
                      tags$dt("Post-hoc Comparisons"),
                      tags$dd("Dunnett (vs reference, default) or Tukey HSD (all pairs). Applies to the Comparisons tab in LME4 mode."),
                      tags$dt("Extrapolation Points"),
                      tags$dd("Number of recent timepoints used to fit a slope for animals that died before study end, then extrapolate their volume to the final day. Set to 0 to use observed data only."),
                      tags$dt("Time Column"),
                      tags$dd("Numeric study day, or calculated from a date column mapped in Data Upload (days elapsed from first observation).")
                    )
                  )
                ),

                br()
              )
            )
          ),

          tabPanel("Summary",
            uiOutput(ns("tg_summary_tab_info")),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("tg_summary"))
            )
          ),
          
          tabPanel("ANOVA",
            uiOutput(ns("tg_anova_tab_info")),
            shinycssloaders::withSpinner(
              DTOutput(ns("tg_anova_table"))
            ),
            downloadButton(ns("download_tg_anova"), "Download ANOVA")
          ),
          
          tabPanel("Comparisons",
            uiOutput(ns("tg_comparisons_tab_info")),
            shinycssloaders::withSpinner(
              DTOutput(ns("tg_comparisons_table"))
            ),
            downloadButton(ns("download_tg_comparisons"), "Download Comparisons")
          ),
          
          tabPanel("Treatment Effects",
            uiOutput(ns("tg_treatment_effects_tab_info")),
            shinycssloaders::withSpinner(
              DTOutput(ns("tg_effects_table"))
            ),
            downloadButton(ns("download_tg_effects"), "Download Effects")
          ),
          
          tabPanel("Growth Rates",
            tab_info(
              strong("Per-animal exponential growth rates "), "estimated from each subject\u2019s individual tumor volume trajectory using linear regression on log\u2081\u209a(volume) vs time. Requires \u22653 time points per animal. ",
              strong("All rate, fold change, and percent change values are on a per-day basis."),
              tags$ul(
                tags$li(strong("Treatment"), " \u2014 treatment group."),
                tags$li(strong("ID"), " \u2014 individual animal identifier."),
                tags$li(strong("Cage"), " \u2014 cage identifier."),
                tags$li(strong("Growth Rate"), " \u2014 slope of the log\u2081\u209a(volume)\u2013time regression (units: log-volume per day). Represents the instantaneous exponential growth constant: a value of 0.10 means ~10% volume increase per day on a continuous compounding basis."),
                tags$li(strong("GR 95% CI Lower/Upper Bound"), " \u2014 95% confidence interval for the log-scale growth rate."),
                tags$li(strong("Fold Change"), " \u2014 exp(Growth Rate): multiplicative volume change per day. A value of 1.10 means the tumor is 10% larger each day; values below 1.0 indicate regression."),
                tags$li(strong("FC 95% CI Lower/Upper Bound"), " \u2014 exp(GR CI bounds): back-transformed confidence interval for the daily fold change."),
                tags$li(strong("Percent Change"), " \u2014 (Fold Change \u2212 1) \u00d7 100: daily volume change expressed as a percentage (e.g. +10.5%/day or \u22123.2%/day)."),
                tags$li(strong("PC 95% CI Lower/Upper Bound"), " \u2014 (exp(GR CI) \u2212 1) \u00d7 100: back-transformed CI bounds for the percent change."),
                tags$li(strong("R\u00b2"), " \u2014 coefficient of determination from the per-animal log-linear regression. Values close to 1.0 indicate a good exponential fit; lower values suggest irregular growth or too few time points.")
              ),
              em("These are descriptive per-animal values used to visualise within-group spread and for sample size estimation. For group-level inference, use the Treatment Effects tab.")
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("tg_growth_rates_table"))
            ),
            downloadButton(ns("download_tg_rates"), "Download Growth Rates")
          ),
          
          tabPanel("Data Summary",
            tab_info(
              strong("Descriptive statistics for tumor volume "), "by treatment group and timepoint. Shows mean, standard deviation, and standard error across subjects measured at each study day.",
              tags$ul(
                tags$li(strong("Treatment"), " \u2014 treatment group."),
                tags$li(strong("Day"), " \u2014 study day (measurement timepoint)."),
                tags$li(strong("Mean"), " \u2014 mean tumor volume (mm\u00b3) for the group at that timepoint."),
                tags$li(strong("SD"), " \u2014 standard deviation of volume across animals in the group."),
                tags$li(strong("SEM"), " \u2014 standard error of the mean (SD / \u221aN)."),
                tags$li(strong("N"), " \u2014 number of animals with a measurement on that day.")
              ),
              em("Useful for inspecting raw descriptive trends across time and confirming data entry. Statistical comparisons should use the ANOVA, Comparisons, or AUC Summary tabs as appropriate.")
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("tg_data_summary_table"))
            ),
            downloadButton(ns("download_tg_data_summary"), "Download Data Summary")
          ),
          
          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  selectInput(ns("tg_plot_selector"), "Plot Type",
                    choices = c(
                      "Tumor Growth Curves"      = "growth_curves",
                      "Growth Rates Plot"        = "growth_rate",
                      "Area Under the Curve Plot" = "auc",
                      "Waterfall Plot"           = "waterfall",
                      "Effect Size Forest Plot"  = "effect_forest"
                    )
                  ),
                  ctrl_section("group-order", "sort", "Groups",
                    uiOutput(ns("tg_group_order_color_ui"))
                  ),
                  ctrl_section("elements", "layer-group", "Elements",
                    conditionalPanel(ns = ns,
                      "['growth_rate','auc'].indexOf(input.tg_plot_selector) >= 0",
                      selectInput(ns("tg_opt_center_bar"), "Center line / bar",
                        choices = c("Mean" = "mean", "Median" = "median",
                                    "Geometric Mean" = "geometric_mean",
                                    "None" = "none"),
                        selected = "mean"
                      ),
                      selectInput(ns("tg_opt_error_bar"), "Error bars",
                        choices  = c("None" = "none", "SEM" = "SEM",
                                     "SD" = "SD", "95% CI" = "CI"),
                        selected = "SEM"
                      )
                    ),
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'growth_curves'",
                      checkboxInput(ns("tg_opt_group_line"),
                                    "Show group average line", value = TRUE),
                      checkboxInput(ns("tg_opt_show_extrap"),
                                    "Show extrapolated curves/points", value = TRUE)
                    ),
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'effect_forest'",
                      checkboxInput(ns("tg_opt_forest_show_pval"),
                                    "Show p-value",     value = TRUE),
                      checkboxInput(ns("tg_opt_forest_show_effect"),
                                    "Show effect size", value = FALSE),
                      checkboxInput(ns("tg_opt_forest_show_ci"),
                                    "Show 95% CI",      value = FALSE)
                    ),
                    selectInput(ns("tg_opt_legend_pos"), "Legend position",
                      choices = c("None" = "none",
                                  "Top" = "top", "Bottom" = "bottom",
                                  "Left" = "left", "Right" = "right",
                                  "In plot - top left"  = "inplot_tl",
                                  "In plot - top right" = "inplot_tr",
                                  "In plot - bottom left" = "inplot_bl",
                                  "In plot - bottom right" = "inplot_br"),
                      selected = "right"
                    ),
                    textInput(ns("tg_plot_title"), "Title",
                              placeholder = "Leave blank for default"),
                    selectInput(ns("tg_plot_title_align"), "Title alignment",
                      choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                      selected = "center"),
                    textInput(ns("tg_plot_subtitle"), "Subtitle",
                              placeholder = "Optional subtitle"),
                    textInput(ns("tg_plot_xaxis_label"), "X-axis label",
                              placeholder = "Leave blank for default"),
                    textInput(ns("tg_plot_yaxis_label"), "Y-axis label",
                              placeholder = "Leave blank for default"),
                    textInput(ns("tg_plot_legend_label"), "Legend label",
                              placeholder = "Leave blank for default")
                  ),
                  ctrl_section("scaling", "expand", "Scaling",
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'growth_curves'",
                      selectInput(ns("tg_opt_yscale"), "Y-axis scale",
                        choices = c("Linear" = "linear",
                                    "Log10"  = "log10"),
                        selected = "linear"
                      ),
                      numericInput(ns("tg_opt_xstart"), "X-axis start (Day)",
                                   value = 0, min = 0, step = 1)
                    ),
                    sliderInput(ns("tg_opt_text_size"), "Text size (pt)",
                                min = 10, max = 24, value = 16, step = 1),
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'growth_curves'",
                      sliderInput(ns("tg_opt_growth_pt_size"), "Point size",
                                  min = 0.5, max = 10, value = 4, step = 0.5)
                    ),
                    conditionalPanel(ns = ns,
                      "['growth_rate','auc'].indexOf(input.tg_plot_selector) >= 0",
                      sliderInput(ns("tg_opt_point_size"),   "Point size",
                                  0.5, 10, 4, step = 0.5),
                      sliderInput(ns("tg_opt_jitter_width"), "Jitter width",
                                  0, 0.5, 0.1, step = 0.02)
                    ),
                    conditionalPanel(ns = ns,
                      "input.tg_plot_selector == 'effect_forest'",
                      sliderInput(ns("tg_opt_forest_pt_size"), "Point size",
                                  min = 1, max = 10, value = 3, step = 0.5),
                      sliderInput(ns("tg_opt_forest_ann_size"), "Annotation text size",
                                  min = 8, max = 22, value = 12, step = 1)
                    )
                  ),
                  ctrl_section("export", "arrows-alt", "Plot Size / Export",
                    fluidRow(
                      column(6, numericInput(ns("tg_plot_width"),  "Width (px)",  800, min = 200, max = 5000, step = 50)),
                      column(6, numericInput(ns("tg_plot_height"), "Height (px)", 550, min = 200, max = 5000, step = 50))
                    ),
                    fluidRow(
                      column(6, numericInput(ns("tg_plot_dpi"), "DPI", 100, min = 72, max = 300, step = 1)),
                      column(6, tags$div(style = "margin-top: 25px;", uiOutput(ns("tg_download_plot_btn"))))
                    )
                  )
                )
              ),
              column(9,
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'growth_curves'",
                  h4("Tumor Growth Curves"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_growth"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'growth_rate'",
                  h4("Growth Rates Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_rates"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'auc'",
                  h4("Area Under the Curve Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_auc"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'waterfall'",
                  h4("Waterfall Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_waterfall"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.tg_plot_selector == 'effect_forest'",
                  h4("Effect Size Forest Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_plot_effect_forest"), width = "auto", height = "auto")
                  )
                ),

              )
            )
          ),
          
          tabPanel("Diagnostics",
            fluidRow(
              column(6,
                h4("Q-Q Plot of Residuals"),
                shinycssloaders::withSpinner(
                  plotOutput(ns("tg_diag_qq"), height = "400px")
                )
              ),
              column(6,
                h4("Residuals vs Fitted"),
                shinycssloaders::withSpinner(
                  plotOutput(ns("tg_diag_resid_fitted"), height = "400px")
                )
              )
            )
          ),
          
          tabPanel(
            title = tags$span(icon("file-zipper"), " Download All"),
            value = "__download_all__",
            wellPanel(
              style = "text-align:center; padding:30px;",
              h4("Export All Results"),
              p("Download every available result table (CSV) and plot (PNG) ",
                "as a single ZIP file."),
              downloadButton(ns("download_tg_all"), "Download All Results",
                             class = "btn-lg btn-primary",
                             icon = icon("download"))
            )
          )
        )
      )
    )
  )
}

#' Help/About UI
#' @keywords internal
helpUI <- function() {
  fluidPage(
    h2("Help & Documentation"),

    div(
      class = "alert alert-warning",
      style = "border-left: 5px solid #e2a000;",
      h4(icon("triangle-exclamation"), " Disclaimer"),
      p(
        strong("mouseExperiment"), " is provided \"as is\" for research purposes only,",
        " without warranty of any kind, express or implied. The authors and contributors",
        " accept no liability for errors, inaccuracies, or any consequences — including",
        " but not limited to incorrect results, data loss, or decisions made on the basis",
        " of outputs from this software.",
        " All results must be independently validated before use in publications, reports,",
        " regulatory submissions, clinical decisions, or any other application."
      )
    ),

    hr(),

    fluidRow(
      column(6,
        wellPanel(
          h3(icon("book"), " Documentation"),
          p("Full documentation, tutorials, and function reference for the ",
            strong("mouseExperiment"), " package:"),
          tags$a(
            href  = "https://github.com/sciOmics/mouseExperiment",
            target = "_blank",
            class = "btn btn-primary btn-lg",
            icon("github"), " GitHub Repository"
          )
        )
      ),
      column(6,
        wellPanel(
          h3(icon("bug"), " Report Issues"),
          p("Found a bug or want to request a new feature? Open an issue on GitHub:"),
          tags$a(
            href  = "https://github.com/sciOmics/mouseExperiment/issues",
            target = "_blank",
            class = "btn btn-warning btn-lg",
            icon("exclamation-circle"), " Submit Issue / Feature Request"
          )
        )
      )
    ),

    hr(),

    h3("Column Naming"),
    p("The dashboard looks for common column names like:"),
    tags$ul(
      tags$li(strong("ID:"), "ID, Mouse_ID, Animal_ID, Subject"),
      tags$li(strong("Time:"), "Day, Time, Timepoint"),
      tags$li(strong("Volume:"), "Volume, Tumor_Volume"),
      tags$li(strong("Treatment:"), "Treatment, Group"),
      tags$li(strong("Cage:"), "Cage, Cage_ID"),
      tags$li(strong("Dose (optional):"), "Dose, dose, DoseLevel"),
      tags$li(strong("Survival/Censor (optional):"), "Survival_Censor, Censor, Event \u2014 numeric column coded as ", strong("1 = endpoint reached (e.g. tumour exceeded IACUC volume limit, or death)"), " and ", strong("0 = censored (did not reach endpoint by study end)"), ". Alternatively, leave unset and use the volume-threshold derivation option in the Survival tab sidebar.")
    ),

    h3("Statistical Methods"),
    p("The analysis uses linear mixed-effects models (lme4) or area under the curve (AUC) methods."),
    p("Consult the ",
      tags$a("mouseExperiment documentation",
             href = "https://github.com/sciOmics/mouseExperiment", target = "_blank"),
      " for detailed method information."),

    h3("Troubleshooting"),
    p(strong("File won't upload:"), "Ensure it is a CSV or RDA file."),
    p(strong("Columns not detected:"), "Check that column names match expected patterns above."),
    p(strong("Analysis tab greyed out:"), "Upload and configure data first. Survival requires a censor column; Dose Response requires a dose column; Drug Synergy requires at least 3 treatment groups (e.g. control, Drug A, and combination)."),
    p(strong("Analysis fails:"), "Check that required columns are selected and contain valid data.")
  )
}

#' About UI
#' @keywords internal
citationUI <- function() {
  pkg_version   <- tryCatch(as.character(utils::packageVersion("mouseExperiment")),      error = function(e) "(not installed)")
  dash_version  <- tryCatch(as.character(utils::packageVersion("mouseExperimentDashboard")), error = function(e) "(not installed)")
  pkg_date      <- tryCatch({
    d <- utils::packageDescription("mouseExperiment")$Date
    if (is.null(d) || is.na(d)) "2026-04-03" else d
  }, error = function(e) "2026-04-03")
  dash_date     <- tryCatch({
    d <- utils::packageDescription("mouseExperimentDashboard")$Date
    if (is.null(d) || is.na(d)) "2026-04-03" else d
  }, error = function(e) "2026-04-03")

  fluidPage(
    h2("About"),

    fluidRow(
      column(6,
        wellPanel(
          h4(icon("box-open"), " mouseExperiment"),
          tags$table(
            class = "table table-condensed",
            style = "margin-bottom:0;",
            tags$tbody(
              tags$tr(
                tags$td(strong("Version")),
                tags$td(pkg_version)
              ),
              tags$tr(
                tags$td(strong("Last Updated")),
                tags$td(pkg_date)
              ),
              tags$tr(
                tags$td(strong("Repository")),
                tags$td(tags$a("github.com/sciOmics/mouseExperiment",
                               href = "https://github.com/sciOmics/mouseExperiment",
                               target = "_blank"))
              )
            )
          )
        )
      ),
      column(6,
        wellPanel(
          h4(icon("gauge"), " mouseExperimentDashboard"),
          tags$table(
            class = "table table-condensed",
            style = "margin-bottom:0;",
            tags$tbody(
              tags$tr(
                tags$td(strong("Version")),
                tags$td(dash_version)
              ),
              tags$tr(
                tags$td(strong("Last Updated")),
                tags$td(dash_date)
              )
            )
          )
        )
      )
    )
  )
}

#' Analysis Log UI
#'
#' Displays a live panel of analysis log entries emitted during the session.
#' @param id Shiny module id.
#' @keywords internal
analysisLogUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2(icon("clipboard-list"), " Analysis Log"),
    tab_info(
      "Events logged during your session (data uploads, analysis runs, errors).",
      " Logs are stored in memory and cleared when the session ends."
    ),

    fluidRow(
      column(12,
        wellPanel(
          fluidRow(
            column(3, actionButton(ns("refresh_log"), "Refresh",
                                   icon = icon("arrows-rotate"),
                                   class = "btn-sm btn-default")),
            column(3, actionButton(ns("clear_log"), "Clear log",
                                   icon = icon("trash"),
                                   class = "btn-sm btn-danger")),
            column(6, downloadButton(ns("download_log"), "Download log",
                                     class = "btn-sm btn-info"))
          )
        ),
        verbatimTextOutput(ns("log_output"))
      )
    )
  )
}