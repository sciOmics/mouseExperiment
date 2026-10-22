#' UI Modules for mouseExperimentDashboard

# Shared helper: creates a styled info box for tab descriptions
tab_info <- function(...) {
  div(
    class = "alert alert-info",
    style = "font-size:0.88em; margin:10px 0 8px 0; padding:8px 14px;",
    ...
  )
}

# Shared helper: inline help icon. It shows `text` in a tooltip on hover, on
# keyboard focus, and on a click or tap (www/js/help-tips.js). A bare title
# attribute showed nothing on a click or a tap. The icon sits inside input
# labels, and that script also keeps its clicks from reaching the label, which
# passed them on to the input: the Control group menu opened instead.
help_tip <- function(text) {
  tags$span(
    class = "help-tip text-muted",
    style = "cursor:help; margin-left:4px;",
    tabindex = "0", role = "button", `aria-label` = paste("Help:", text),
    title = text,
    icon("circle-question", class = "fa-sm")
  )
}

# Shared helper: a "How to read this plot" note under a plot's title, collapsed
# until clicked. It used to be a small grey line below each plot, often a
# screen further down, and was easy to miss (2026-10-22).
plot_help <- function(..., summary = "How to read this plot") {
  tags$details(
    class = "plot-help",
    tags$summary(icon("circle-info"), " ", summary),
    tags$div(class = "plot-help-body", ...)
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
                  "Combo Treatment" = "combo",
                  "Dose Levels"     = "dose_levels",
                  "Weight"          = "weight",
                  "All Features"    = "all_features"
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
                       accept = ".csv"),
              uiOutput(ns("file_validation_msg"))
            )
          )
        )
      )
    ),
    
    # Column mapping section (shown after file upload)
    uiOutput(ns("column_mapping_ui")),

    # Study design (R20-D.H): set once here, used by every tab.
    uiOutput(ns("study_design_ui")),
    
    # Dataset summary (shown after successful upload + column config)
    uiOutput(ns("data_summary_ui")),

    # Subset filter (CODE_REVIEW.md Round 2 F.5) — exclude mice / cages /
    # treatments. Sits between the Dataset Summary and Data Preview so
    # users see the summary first, choose what to filter, click Apply,
    # then the Data Preview below updates to reflect what they applied.
    subset_filter_ui(ns),

    # Data preview (shown after successful upload)
    uiOutput(ns("data_preview_ui")),

    # Settings preset (column mapping + methods)
    settings_preset_ui(ns)
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
          # R20.D1: dates such as 01/03/2025 are read in the order the
          # column itself requires; this choice applies only when no value
          # decides it.
          selectInput(ns("date_order"),
                      tagList("Date order", help_tip("For dates such as 01/03/2025. With 'Detect from the data', a value such as 15/03/2025 (no month 15) decides the order for the whole column; if no value decides it, you are asked to choose.")),
                      choices = c("Detect from the data" = "auto",
                                  "Day first (DD/MM/YYYY)" = "dmy",
                                  "Month first (MM/DD/YYYY)" = "mdy"),
                      selected = "auto"),
          textInput(ns("reference_date"), "Reference Date (YYYY-MM-DD)",
                   value = "2025-03-24",
                   placeholder = "First measurement date"),
          # R20.D4: a date the preview cannot use is reported here, where it
          # was typed, instead of being swallowed.
          uiOutput(ns("reference_date_status"))
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
        tags$p(tags$small("Optional. Enables necrotic-observation handling in Tumor Growth analysis.")),

        br(),
        # R20-D.H: why each animal left the study.
        selectInput(ns("col_removal_reason"),
                   tagList("Removal Reason",
                           help_tip("Optional. Why each animal left the study, e.g. 'Tumour volume limit', 'Found dead', 'Body condition', 'End of study'. Read from each animal's last non-empty value. Each value is given a meaning under Study Design.")),
                   choices = c("None" = "", col_names),
                   selected = detected_cols$removal_reason %||% ""),
        tags$p(tags$small("Optional. Lets Survival count deaths from other causes as events, and ",
                          "the Toxicity weight-loss threshold separate weight-loss removals from ",
                          "other removals."))
      )
    )
  )
}

#' Study design UI
#'
#' CODE_REVIEW.md R20-D.H. The study's volume units, control group and unit of
#' randomisation are declared once, on Data Upload, and every tab uses them.
#' They replace the per-tab reference and cage selectors, which could
#' disagree with each other.
#'
#' @param treatments Values of the mapped treatment column.
#' @param defaults Named list of initial selections (`volume_units`,
#'   `control_group`, `randomisation_unit`); a missing `volume_units` leaves
#'   the units undeclared.
#' @param ns Namespace function.
#' @keywords internal
study_design_ui <- function(treatments, defaults = list(), ns = NS(NULL)) {
  control <- defaults$control_group %||% ""
  if (!control %in% treatments) control <- pick_reference_group(treatments)
  wellPanel(
    h4(strong("Study Design")),
    tags$p("Set once here; every analysis tab uses these settings."),
    # R20-K: the evaluable-day rule is fixed; it is stated here because it
    # applies to every endpoint result (synergy, therapeutic window,
    # dose-response, AUC).
    tags$p(class = "text-muted", style = "font-size:12px;",
           icon("calendar-check"), " ", strong("Evaluable days: "),
           evaluable_rule_text(),
           " Synergy, the therapeutic window, dose-response and the AUC use only ",
           "days on which every arm they compare is evaluable."),
    fluidRow(
      column(4,
        radioButtons(ns("volume_units"),
          tagList("Tumour volume units", help_tip("The units of the volume column, or of the volume calculated from length and width (millimetres give mm\u00b3). They set the tumour-mass correction on the Toxicity tab and the axis labels.")),
          choices = c("mm\u00b3" = "mm3", "cm\u00b3" = "cm3"),
          selected = defaults$volume_units %||% character(0), inline = TRUE),
        uiOutput(ns("volume_units_status"))
      ),
      column(4,
        selectInput(ns("control_group"),
          tagList("Control group", help_tip("The reference arm for every comparison: tumour growth, survival, body weight, synergy and dose-response.")),
          choices = treatments, selected = control)
      ),
      column(4,
        conditionalPanel(ns = ns, condition = "input.col_cage != ''",
          radioButtons(ns("randomisation_unit"),
            tagList("Treatment was assigned to", help_tip("Whole cages: survival compares arms by moving whole cages between them (a cage-level permutation log-rank), and the tumour-growth permutation test permutes cages. With 2 cages per arm the smallest attainable p-value is 1/3, and with 3 it is 0.1, so a cage-randomised comparison with so few cages cannot reach p < 0.05.")),
            choices = c("Individual animals" = "mouse", "Whole cages" = "cage"),
            selected = defaults$randomisation_unit %||% "mouse")
        )
      )
    ),
    # R20-D.H: the meaning of each removal reason, when that column is mapped.
    uiOutput(ns("removal_reasons_ui"))
  )
}

#' Removal-reason mapping UI
#'
#' One row per distinct value of the removal-reason column, each with its
#' category (R20-D.H). The selectors are numbered because reason values are
#' free text.
#'
#' @param values Distinct reason values.
#' @param defaults Named list, value to category code, from a preset.
#' @param col The column's name.
#' @param ns Namespace function.
#' @keywords internal
removal_reason_map_ui <- function(values, defaults = list(), col = "", ns = NS(NULL)) {
  if (!length(values)) {
    return(tags$p(class = "text-muted", style = "font-size:12px;",
                  icon("circle-info"), " The removal-reason column ",
                  tags$code(col), " has no values."))
  }
  guess <- guess_removal_category(values)
  tagList(
    tags$hr(),
    h5(strong("Removal reasons"), tags$small(class = "text-muted",
       " — what each value of ", tags$code(col), " means")),
    tags$p(class = "text-muted", style = "font-size:12px;",
      "Survival can count ", strong("deaths"), " as events. The Toxicity weight-loss threshold ",
      "counts ", strong("weight-loss"), " removals as events, censors ", strong("planned ends"),
      ", and treats the other removals as competing risks. The categories were guessed ",
      "from the wording; check them."),
    lapply(seq_along(values), function(i) {
      d <- defaults[[values[i]]]
      sel <- if (length(d) == 1L && d %in% removal_reason_categories) d else guess[i]
      fluidRow(
        column(6, tags$div(style = "padding-top:7px;", values[i])),
        column(6, selectInput(ns(paste0("reason_cat_", i)), NULL,
                              choices = removal_reason_categories, selected = sel,
                              width = "100%")))
    })
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
          # Explanations sit in help_tip() tooltips, so the controls and the
          # Run button stay in view (26.10.007); the Info tab has the detail.
          selectInput(ns("tg_model_type"),
                     tagList("Statistical Model", help_tip(paste(
                       "LME4: a linear mixed-effects model of each animal's growth over time.",
                       "AUC: compares total tumour burden, the area under each arm's fitted",
                       "curve, over a window that every arm covers; animals that die early",
                       "still count through all their measurements. Bayesian LMM: the LME4",
                       "model with posterior distributions and credible intervals; it runs in",
                       "the background and takes seconds to a few minutes."))),
                     choices = c("Linear Mixed-Effects (LME4)"           = "lme4",
                                "Area Under Curve (AUC)"                 = "auc",
                                "Bayesian LMM (brms)"                    = "bayes"),
                     selected = "lme4"),

          # 2. Reference Group: the study's control group (R20-D.H)
          uiOutput(ns("tg_reference_note")),

          # 3. Post-hoc Comparisons
          selectInput(ns("tg_comparison_type"),
                     tagList("Post-hoc Comparisons", help_tip(paste(
                       "vs Reference Group compares each arm with the control group.",
                       "All Comparisons compares every pair of arms."))),
                     choices = c("vs Reference Group" = "reference",
                                "All Comparisons" = "tukey"),
                     selected = "reference"),
          # CODE_REVIEW.md R3.D3 — p_adjust_method was documented by the backend
          # with four options and exposed by the UI with none, so users got
          # whichever correction each code path happened to hardcode. The
          # backend now honours the argument (backend R3.1), so surface it.
          # Backend H.3 -- opt-in because each permutation refits the model, and
          # because the unit has to be a deliberate statement about the design
          # rather than a default someone accepts without reading.
          checkboxInput(ns("tg_run_perm"),
                       tagList("Run randomisation test of treatment \u00d7 time", help_tip(paste(
                         "A permutation test of whether the arms' growth rates differ, which",
                         "needs neither the model's degrees-of-freedom approximation nor",
                         "normal random effects. Each permutation refits the model, so it",
                         "takes longer."))),
                       value = FALSE),
          conditionalPanel(ns = ns, condition = "input.tg_run_perm == true",
            uiOutput(ns("tg_perm_unit_note")),
            numericInput(ns("tg_perm_n"), "Permutations",
                        value = 999, min = 99, max = 9999, step = 100)
          ),

          selectInput(ns("tg_p_adjust"),
                     tagList("Multiple-comparison correction", help_tip(paste(
                       "Applied over exactly the comparisons shown. Dunnett (vs reference)",
                       "and Tukey (all pairs) use the correlation between comparisons, so",
                       "they are less conservative than Bonferroni or Holm. The AUC model",
                       "offers Holm, Bonferroni, FDR or none."))),
                     choices = tg_p_adjust_choices("lme4"),
                     selected = "dunnett"),

          # 4. Data Transformation (LME4 only)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4'",
            selectInput(ns("tg_transform"),
                       tagList("Data Transformation", help_tip(paste(
                         "Applied to tumour volume before modelling. Log is recommended for",
                         "most data: it stabilises the variance and makes exponential growth",
                         "linear."))),
                       choices = c("Log" = "log",
                                  "Square Root" = "sqrt",
                                  "None" = "none"),
                       selected = "log")
          ),

          # 5. Extrapolation Points — PLOTTING ONLY as of backend v0.8.0
          numericInput(ns("tg_extrapolation_points"),
                      tagList("Extrapolation Points (plots only)", help_tip(paste(
                        "Affects the growth-trajectory plots only: the number of recent",
                        "timepoints used to project the curve of an animal that left the",
                        "study early, drawn dashed. 0 turns it off. The statistics never use",
                        "projected values; the endpoint metrics extrapolate inside the model."))),
                      value = 3, min = 0, max = 10, step = 1),
          uiOutput(ns("tg_extrap_note")),

          # 6. Cage Effect Handling (LME4 only)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4'",
            selectInput(ns("tg_cage_effect"),
                       tagList("Cage Effect Handling", help_tip(paste(
                         "How cage enters the model. Auto follows the design: a fixed effect",
                         "when cages mix treatments, a random intercept when each arm has two",
                         "or more cages of its own, and left out with a warning when an arm",
                         "has a single cage. The legacy option drops nested cages, which",
                         "understates standard errors."))),
                       choices = c("Auto (match the design)" = "auto",
                                  "Include if not collinear (legacy)" = "include_if_not_collinear",
                                  "Always include" = "always_include",
                                  "Never include" = "never_include",
                                  "As random effect" = "as_random_effect"),
                       selected = "auto")
          ),

          # 7. Random Effects (LME4 and Bayesian)
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'lme4' || input.tg_model_type == 'bayes'",
            # CODE_REVIEW.md R20.5 / R20-D.H: random slopes by default. With a
            # random intercept only, animals' differing growth rates land in
            # the residual and the Treatment x Day test rejected a true null in
            # 73 % of simulated studies.
            selectInput(ns("tg_random_effects"),
                       tagList("Random Effects", help_tip(paste(
                         "Random slopes give each animal its own baseline and growth rate.",
                         "Intercept only assumes every animal in an arm grows at the same",
                         "rate; when they do not, the treatment \u00d7 time test finds",
                         "differences that are not there. If a model cannot be fitted, the",
                         "analysis falls back to simpler random effects and says so. Bayesian",
                         "mode treats Uncorrelated as Random slopes and None as Intercept only."))),
                       choices = c("Random slopes (recommended)" = "slope",
                                   "Uncorrelated random slopes" = "slope_uncorrelated",
                                   "Intercept only (all animals in an arm grow at the same rate)" = "intercept_only",
                                   "None" = "none"),
                       selected = "slope")
          ),

          # 8. Bayesian-specific controls
          conditionalPanel(ns = ns,
            condition = "input.tg_model_type == 'bayes'",
            hr(),
            h5(strong("Bayesian Settings")),
            selectInput(ns("tg_prior_strength"),
              tagList("Prior Strength", help_tip(paste(
                "The prior on the treatment effects. Skeptical, N(0, 0.25), needs strong",
                "data to support a large effect; Weakly informative is N(0, 1);",
                "Informative N(0, 0.5); Diffuse N(0, 2.5). Manual sets each prior below."))),
              choices = c(
                "Skeptical (recommended)" = "skeptical",
                "Weakly informative"      = "weakly_informative",
                "Informative"             = "informative",
                "Diffuse"                 = "diffuse",
                "Manual specification"    = "manual"
              ),
              selected = "skeptical"
            ),
            # Manual prior specification — shown only when "Manual" is selected
            conditionalPanel(ns = ns,
              condition = "input.tg_prior_strength == 'manual'",
              tags$hr(style = "margin:8px 0"),
              tags$p(tags$strong("Fixed effects (b)"), style = "margin-bottom:4px; font-size:0.9em;"),
              fluidRow(
                column(6, selectInput(ns("tg_prior_b_dist"), NULL,
                  choices = c("normal" = "normal", "student_t" = "student_t"),
                  selected = "normal")),
                column(3, numericInput(ns("tg_prior_b_mu"),  "μ",  value = 0,    step = 0.05)),
                column(3, numericInput(ns("tg_prior_b_sd"),  "σ",  value = 0.25, step = 0.05, min = 0.01))
              ),
              conditionalPanel(ns = ns,
                condition = "input.tg_prior_b_dist == 'student_t'",
                numericInput(ns("tg_prior_b_df"), "df (student_t)", value = 3, min = 1, step = 1)
              ),
              tags$p(tags$strong("Intercept"), style = "margin:6px 0 4px; font-size:0.9em;"),
              fluidRow(
                column(6, selectInput(ns("tg_prior_int_dist"), NULL,
                  choices = c("normal" = "normal", "student_t" = "student_t"),
                  selected = "normal")),
                column(3, numericInput(ns("tg_prior_int_mu"), "μ",  value = 0,    step = 0.1)),
                column(3, numericInput(ns("tg_prior_int_sd"), "σ",  value = 0.625, step = 0.05, min = 0.01))
              ),
              conditionalPanel(ns = ns,
                condition = "input.tg_prior_int_dist == 'student_t'",
                numericInput(ns("tg_prior_int_df"), "df (student_t)", value = 3, min = 1, step = 1)
              ),
              tags$p(tags$strong("Random-effect SD (sd)"), style = "margin:6px 0 4px; font-size:0.9em;"),
              fluidRow(
                column(6, selectInput(ns("tg_prior_sd_dist"), NULL,
                  choices = c("exponential" = "exponential", "half-normal" = "half_normal"),
                  selected = "exponential")),
                column(6, numericInput(ns("tg_prior_sd_param"),
                  label = uiOutput(ns("tg_prior_sd_param_label")),
                  value = 2, min = 0.01, step = 0.1))
              ),
              tags$p(tags$strong("Residual SD (sigma)"), style = "margin:6px 0 4px; font-size:0.9em;"),
              fluidRow(
                column(6, selectInput(ns("tg_prior_sigma_dist"), NULL,
                  choices = c("exponential" = "exponential", "half-normal" = "half_normal"),
                  selected = "exponential")),
                column(6, numericInput(ns("tg_prior_sigma_param"),
                  label = uiOutput(ns("tg_prior_sigma_param_label")),
                  value = 2, min = 0.01, step = 0.1))
              )
            ),
            numericInput(ns("tg_bayes_chains"), "MCMC chains",
                         value = 4, min = 1, max = 8, step = 1),
            numericInput(ns("tg_bayes_warmup"), "Warmup per chain",
                         value = 1000, min = 200, max = 5000, step = 100),
            numericInput(ns("tg_bayes_iter"), "Post-warmup draws per chain",
                         value = 500, min = 200, max = 5000, step = 100),
            selectInput(ns("tg_bayes_backend"), "Stan backend",
                        choices  = c("rstan", "cmdstanr"),
                        selected = "rstan"),
            div(class = "alert alert-warning", style = "font-size:0.85em; padding:8px; margin-top:6px;",
              icon("clock"), " The first fit of a model compiles it (about 20 s); later fits ",
              "with the same structure reuse it. Fits run one at a time in the background: ",
              "you can close the browser and resume or cancel from the banner when you return."
            )
          ),

          # Necrotic handling: rendered server-side, shown once a necrotic column
          # is mapped on Data Upload (R20-D.V)
          uiOutput(ns("necrotic_handling_ui")),

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
                    p("The Tumor Growth tab provides three statistical modeling approaches selectable from the ",
                      strong("Statistical Model"), " dropdown:"),
                    tags$ul(
                      tags$li(strong("Linear Mixed-Effects Model (LME4)"), " — fits a mixed-effects model to repeated tumor volume measurements. Accounts for inter-animal variability and unequal observation times. Recommended for most longitudinal datasets."),
                      tags$li(strong("Area Under Curve (AUC)"), " — compares each arm’s total tumor burden as the area under its fitted growth curve, over the window in which every arm still has enough animals on study."),
                      tags$li(strong("Bayesian LMM (brms)"), " — the LME4 model fitted by MCMC, returning full posterior distributions instead of point estimates and p-values.")
                    ),
                    p(tags$small(class = "text-muted",
                      "R14.3: this tab does ", strong("not"), " report Tumor Growth Inhibition. ",
                      "TGI is an endpoint summary — a ratio of two group means on one day — ",
                      "while every model here estimates a growth ", em("trajectory"), ", and ",
                      "collapsing a trajectory to a single-day ratio discards most of what was ",
                      "measured. TGI is computed where it belongs to the question being asked: ",
                      "the ", strong("Toxicity"), " tab (as the efficacy half of the therapeutic ",
                      "window) and the ", strong("Drug Synergy"), " tab (as the fractional effect ",
                      "the Bliss model acts on). Both expose an endpoint-estimand ",
                      "selector, because a ratio of endpoint means is exactly where early ",
                      "removals bias the answer."
                    )),
                    div(class = "alert alert-info", style = "font-size:0.9em; margin-bottom:0;",
                      icon("circle-info"), " ",
                      strong("Which endpoint?"), " Animals removed before study end make “the final timepoint” ",
                      "ambiguous, and taking it literally computes TGI among survivors only — which flatters ",
                      "the control arm, because the fastest-growing controls are exactly the ones that left. ",
                      "The Toxicity and Synergy tabs expose an ", strong("Endpoint estimand"), " selector for this; ",
                      "here TGI is derived from the fitted model, which uses every observation each animal ",
                      "contributed rather than only those present at the end."
                    )
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
                      tags$dt("Transformation"), tags$dd("Log (default) stabilises variance and linearises exponential growth. Square root suits moderate heteroscedasticity. Use 'None' only for already-linear data."),
                      tags$dt("Random Effects"), tags$dd("'Random slopes' (default) gives each animal its own baseline volume and its own growth rate. If that model cannot be fitted, the analysis falls back to uncorrelated slopes and then to an intercept only, and the Summary says so. 'Intercept only' assumes every animal in an arm grows at the same rate, which makes the treatment \u00d7 time test find differences that are not there. 'None' drops the animal-level term entirely and is appropriate only when each animal contributes a single observation."),
                      tags$dt("Cage Effect"),
                      tags$dd(strong("Auto (default)"), " classifies the housing design and applies the matching treatment, because the correct handling depends on the design rather than on preference:",
                              tags$ul(
                                tags$li(em("Crossed"), " \u2014 a cage holds more than one treatment. Cage is separable from treatment and enters as a fixed effect."),
                                tags$li(em("Nested, replicated"), " \u2014 one treatment per cage, several cages per arm. Cage nests within treatment as a random effect; ignoring it treats cagemates as independent and understates standard errors."),
                                tags$li(em("Nested, partly replicated"), " \u2014 some arms have several cages and others one. Cage stays in the model as a random effect, estimated from the replicated arms, and the Summary names the single-cage arms."),
                                tags$li(em("Nested, confounded"), " \u2014 one cage per treatment. Cage and treatment are the same variable and no model can separate them. The result reports the effect as inseparable rather than silently crediting it to treatment.")
                              ),
                              "The remaining options force a specific handling. ", em("Include if not collinear (legacy)"),
                              " is the pre-v0.12 behaviour, kept so earlier analyses can be reproduced.")
                    ),

                    hr(),
                    h5(strong("ANOVA Tab")),
                    p("F-tests with Satterthwaite degrees of freedom, reported as ",
                      code("NumDF"), " / ", code("DenDF"), " / ", code("F value"), " / ", code("Pr(>F)"),
                      ", in this order: whether the arms' growth rates differ (treatment \u00d7 time, the primary ",
                      "analysis target); whether the arms differ at the mean study day; and whether there is growth ",
                      "on average. A ", code("Tests"), " column states each question."),
                    p(tags$small(class = "text-muted",
                      "The treatment row compares the arms at the mean study day, not at day 0, so it does not ",
                      "change when time is counted from implantation instead of randomisation. With few animals ",
                      "per arm, the randomisation test below makes fewer assumptions."
                    )),

                    hr(),
                    h5(strong("Comparisons Tab")),
                    p("Pairwise group comparisons at the mean observation day using estimated marginal means (EMMs). ",
                      strong("Two independent settings"), " control this, and conflating them is the usual mistake:"),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Post-hoc Comparisons"),
                      tags$dd("Which comparisons are made \u2014 ", em("vs Reference Group"), " (each treatment against the reference only) or ", em("All Comparisons"), " (every pair). This defines the family."),
                      tags$dt("Correction"),
                      tags$dd("How that family is corrected: Dunnett (exact, multivariate-t), Tukey HSD, Bonferroni, Holm, Benjamini-Hochberg FDR, or None. Dunnett is the default and is the right partner for a vs-reference family \u2014 it accounts for the correlation between comparisons that share a control arm, which Bonferroni ignores and therefore over-corrects.")
                    ),
                    p(tags$small(class = "text-muted",
                      "The correction is applied over exactly the set of comparisons requested, and the table header ",
                      "states which correction was actually used rather than asserting one. Comparisons are ",
                      "marginalised at the mean study day; for differences in trajectory across time, see the ",
                      "Treatment Effects tab."
                    )),

                    hr(),
                    h5(strong("Randomisation Test")),
                    p("An optional exact test of the Treatment \u00d7 Day interaction that avoids both the ",
                      "denominator-degrees-of-freedom approximation and the assumption of normally distributed ",
                      "random effects. Treatment labels are permuted and the model refitted; the p-value is the ",
                      "proportion of permutations giving a statistic at least as extreme as the observed one."),
                    div(class = "alert alert-warning", style = "font-size:0.9em;",
                      icon("triangle-exclamation"), " ",
                      strong("You must state the unit that was randomised."), " If animals were assigned to arms ",
                      "individually, permute animals. If whole cages were assigned, permute cages \u2014 permuting ",
                      "animals then breaks the exchangeability the test relies on and yields a p-value that is too ",
                      "small. This cannot be inferred from the data; it is a fact about how the experiment was run."),
                    p(tags$small(class = "text-muted",
                      "The panel reports the number of exchangeable units, how many permutations converged, whether ",
                      "the permutation set was enumerated exhaustively, and the smallest p-value the design can ",
                      "attain. With few units that floor can exceed 0.05, in which case no amount of signal can ",
                      "produce a significant result \u2014 and the asymptotic p-value beside it should not be trusted."
                    )),

                    hr(),
                    h5(strong("Treatment Effects Tab")),
                    p("Estimated marginal means (EMMs) of transformed volume per treatment group, at the mean observation day: the day the treatment row of the ANOVA compares the arms at. How the arms diverge over time is in the growth-rate test (ANOVA) and the Growth Rates tab."),

                    hr(),
                    h5(strong("Growth Rates Tab")),
                    p("Per-animal exponential growth rates estimated by log-linear regression of each animal’s volume trajectory independently. Reported as growth rate (log-scale slope), fold change, and % change per day with 95 % CIs and R². These are descriptive; for group-level inference use the ANOVA and Comparisons tabs."),

                    hr(),
                    h5(strong("Diagnostics Tab")),
                    p("Four residual plots, each answering a different question about whether the model's ",
                      "assumptions hold — which is what makes the p-values above mean anything:"),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Q–Q plot"),
                      tags$dd("Are residuals normally distributed? Systematic curvature at the tails suggests the transformation choice should be revisited."),
                      tags$dt("Residuals vs fitted"),
                      tags$dd("Is the mean structure right? Curvature means log-linear growth does not hold over the whole study — read growth-rate comparisons with care."),
                      tags$dt("Scale–location"),
                      tags$dd("Is the variance constant on the modelling scale? A rising trend means heteroscedasticity survived the transformation."),
                      tags$dt("Random-effects Q–Q"),
                      tags$dd("Are the per-animal effects themselves normal? This is a separate assumption from residual normality and the one most often overlooked. If it fails, prefer the randomisation test, which assumes neither.")
                    ),
                    p(tags$small(class = "text-muted",
                      "Influence measures (Cook's distance, DFBETAS) are also available: they identify individual ",
                      "animals whose removal would materially change the conclusion. A result resting on one mouse ",
                      "is worth knowing about before publication rather than after."
                    ))
                  )
                ),

                # ── AUC Mode ─────────────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-warning", "AUC Mode"), " — Analyses & Outputs")),
                  div(class = "panel-body",
                    p("Each arm's AUC is the area under its fitted growth curve, the geometric-mean curve of a ",
                      "mixed model of log volume with a smooth (natural spline) time course per arm and ",
                      "per-animal random slopes. It is integrated over one window common to every arm: from the ",
                      "first study day to the last day on which every arm is evaluable (at least 50 % of its ",
                      "enrolled animals, and at least 3, still on study). AUC is in volume \u00d7 days."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("AUC by arm"), tags$dd("Each arm's AUC with a 95 % interval, from draws of the model."),
                      tags$dt("Comparisons"),
                      tags$dd("The ratio of two arms' AUCs and their difference, with 95 % intervals and a Wald ",
                              "p-value from the same draws, corrected over the family selected in the sidebar. ",
                              "A ratio below 1 means less tumor burden than the comparator."),
                      tags$dt("ANOVA"),
                      tags$dd("A Wald test that every arm has the same AUC (on the log scale)."),
                      tags$dt("Per-animal AUC"),
                      tags$dd("Each animal's own trapezoidal AUC, for description only.")
                    ),
                    div(class = "alert alert-info", style = "font-size:0.9em;",
                      icon("circle-info"), " ",
                      strong("Animals that die or are removed early "), "still count. Every measurement an ",
                      "animal had goes into the model, which estimates its growth rate from them and keeps it in ",
                      "the arm's curve as if the animal had continued. So the AUC is not pulled down when the ",
                      "largest tumours leave the study, as an average of the animals still on study would be. ",
                      "In simulations with 10 animals per arm and about a third of the treated arm lost, removal at ",
                      "the volume limit moved the AUC by \u22120.2 % (the average of the survivors: \u221213 %). ",
                      "Deaths that the earlier measurements only partly foreshadow left a small underestimate, ",
                      "\u22123 % to \u22125 % (survivors: \u221218 % to \u221221 %). When many animals die for ",
                      "reasons their growth does not show, read the Survival analysis alongside the AUC."
                    ),
                    div(class = "alert alert-warning", style = "font-size:0.9em; margin-bottom:0;",
                      icon("triangle-exclamation"), " ",
                      strong("Why not each animal's own AUC? "), "An animal euthanised on day 14 has a smaller ",
                      "trapezoidal AUC than one followed to day 28 because it was measured for less time, not because ",
                      "its tumour was smaller. Controls leave first, so comparing per-animal AUCs made an effective ",
                      "drug look worse than control in about half of simulated studies. The per-animal values are ",
                      "shown under \u2018AUC Individual Data\u2019 for description, not compared.")
                  )
                ),

                # ── Bayesian LMM Mode ─────────────────────────────────────────
                div(class = "panel panel-default",
                  div(class = "panel-heading",
                    h4(class = "panel-title", tags$span(class = "label label-success", "Bayesian Mode"), " — Method & Outputs")),
                  div(class = "panel-body",

                    h5(strong("Method")),
                    p("The Bayesian LMM fits the same formula as LME4 using full Bayesian inference via ",
                      tags$a("brms", href = "https://paul-buerkner.github.io/brms/", target = "_blank"),
                      " / Stan. Instead of point estimates and p-values, every parameter has a ",
                      strong("posterior distribution"), ": a probability distribution over plausible values ",
                      "given the data and priors."),
                    tags$blockquote(tags$code(
                      "log(Volume) ~ Treatment × Day + (1 | ID)"
                    )),
                    p("Results are summarised as ", strong("posterior means"),
                      " with 95 % credible intervals. A 95 % CrI has a direct probability ",
                      "interpretation: there is a 95 % posterior probability that the true ",
                      "parameter lies within the interval — unlike a frequentist confidence ",
                      "interval, which is a statement about the procedure rather than the parameter."),
                    p(tags$small(class = "text-muted",
                      "R14.7: two different interval types appear under the same ",
                      code("Lower_95_CrI"), " / ", code("Upper_95_CrI"), " headings, and it is ",
                      "worth knowing which you are reading. The ", strong("Summary"),
                      " tab reports raw model parameters straight from brms: an ",
                      em("equal-tailed"), " interval, cutting 2.5 % from each tail. The ",
                      strong("Treatment Effects"), " and ", strong("Comparisons"),
                      " tabs come from emmeans, which returns a genuine ", em("highest-density"),
                      " interval — the shortest interval containing 95 % of the posterior. ",
                      "The two coincide when a posterior is symmetric and diverge when it is ",
                      "skewed, which is exactly the case for variance components and for the ",
                      "back-transformed time ratios in the Bayesian survival model."
                    )),

                    hr(),
                    h5(strong("Prior Strength")),
                    p("Priors express what you believe about parameter values ", em("before"), " seeing the data. ",
                      "They regularise estimates and prevent the model from chasing noise in small samples. ",
                      "All presets use Normal priors on fixed effects and Exponential priors on variance components; only the widths differ."),
                    tags$dl(class = "dl-horizontal",
                      tags$dt("Skeptical"),
                      tags$dd(tags$strong("Default."), " N(0, 0.25) on fixed effects; Exponential(2) on SDs. ",
                              "Encodes the reasonable belief that treatment effects are small — effects > ±0.5 log-units per day ",
                              "are considered very unlikely a priori. Recommended for exploratory preclinical studies where ",
                              "large effects have not been previously established."),
                      tags$dt("Weakly informative"),
                      tags$dd("N(0, 1) on fixed effects; Exponential(1) on SDs. ",
                              "Broader than Skeptical — appropriate when prior work suggests the treatment ",
                              "category can produce moderate-to-large effects."),
                      tags$dt("Informative"),
                      tags$dd("N(0, 0.5) on fixed effects; Exponential(2) on SDs. ",
                              "Between Skeptical and Weakly informative. Useful when effect sizes in the ",
                              "0.1–0.5 log-unit range are plausible based on prior data."),
                      tags$dt("Diffuse"),
                      tags$dd("N(0, 2.5) on fixed effects; Exponential(0.5) on SDs. ",
                              "Near-flat priors; use with caution on small datasets as estimates can be ",
                              "driven entirely by noisy data."),
                      tags$dt("Manual"),
                      tags$dd("Specify the prior family and parameters for each model class independently. ",
                              "Fixed effects (b): Normal or Student-t. Intercept: same. ",
                              "Random-effect SD and residual sigma: Exponential or half-Normal.")
                    ),

                    hr(),
                    h5(strong("Result Tabs")),
                    tags$ul(
                      tags$li(strong("Summary"), " — posterior median and 95 % CrI for every model parameter."),
                      tags$li(strong("Treatment Effects"), " — group-level estimated marginal means (posterior) at the mean study day."),
                      tags$li(strong("Comparisons"), " — posterior contrasts of each group vs the reference."),
                      tags$li(strong("Posterior Check"), " — overlay of observed data density vs draws from the posterior predictive distribution. Systematic mismatch indicates model mis-specification. A coverage plot reports whether the stated credible intervals actually contain the observed values at the claimed rate."),
                      tags$li(strong("Diagnostics"), " (Bayesian view) — Rhat and effective sample size (ESS) per parameter, with trace and rank plots. Rhat > 1.01 indicates the chains have not mixed; low ESS means the posterior summaries are themselves noisy. Rank plots detect the mixing failures that trace plots hide once chains overlap."),
                      tags$li(strong("Influence"), " — Pareto-k values from PSIS leave-one-out. Values above the threshold (0.7 at most; lower with fewer posterior draws, e.g. 0.667 at 1,000) mean the importance sampling is unreliable for that observation, which usually flags an animal the model fits badly rather than a coding error."),
                      tags$li(strong("Growth Rates"), " — per-animal log-linear regression rates (same as LME4, not MCMC-based).")
                    ),

                    div(class = "alert alert-info", style = "margin-top:10px;",
                      icon("clock"), " ",
                      strong("Computation time:"), " the first fit of a model compiles its Stan code (about 20 s); ",
                      "a later fit whose model has the same structure (same number of arms, random effects and prior preset) ",
                      "reuses the compiled model and starts sampling at once. ",
                      "Fits run one at a time on a background worker, so the app stays responsive while a model runs; ",
                      "the progress window shows how many fits are ahead of yours."
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
                      tags$dd("The study's control group, chosen under Study Design on Data Upload. With a dose column mapped, arms are labelled 'Treatment Dose' and the control maps to its own dose label."),
                      tags$dt("Post-hoc Comparisons"),
                      tags$dd("Which family of comparisons to make — vs Reference Group (default) or All Comparisons. Paired with the separate Correction setting."),
                      tags$dt("Correction"),
                      tags$dd("How the chosen family is corrected: Dunnett (exact multivariate-t, default), Tukey HSD, Bonferroni, Holm, Benjamini-Hochberg FDR, or None."),
                      tags$dt("Extrapolation Points"),
                      tags$dd(strong("Affects the growth-trajectory plots only."), " Number of recent timepoints used to estimate a slope when projecting the curve for animals that left the study early; projected segments are drawn dashed. Set to 0 to disable. This no longer feeds any statistical model — treating an imputed point as measured data made standard errors ignore the imputation."),
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
            ),
            # CODE_REVIEW.md Round 2 F.1 — reproducible report export.
            br(),
            div(class = "text-right",
              downloadButton(ns("download_tg_report"),
                             "Download Full Report (HTML)",
                             class = "btn-sm btn-primary")
            )
          ),
          
          tabPanel("ANOVA",
            uiOutput(ns("tg_anova_tab_info")),
            shinycssloaders::withSpinner(
              DTOutput(ns("tg_anova_table"))
            ),
            downloadButton(ns("download_tg_anova"), "Download ANOVA"),

            # R20.D30 / R20.D61: the Welch variance panel is gone; the
            # model-based AUC (backend 0.26.0) has no variance test.

            # Backend H.3 -- randomisation test of the interaction. Reported here
            # because it answers the same question as the ANOVA row above but
            # without the denominator-df or normality approximations.
            uiOutput(ns("tg_perm_test_ui"))
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
            downloadButton(ns("download_tg_rates"), "Download Growth Rates"),
            br(), br(),
            # Per-mouse drill-down (CODE_REVIEW.md Round 2 F.4)
            h5(strong(icon("magnifying-glass"), " Per-mouse drill-down"),
               tags$small(style = "color:#6c757d; margin-left:8px;",
                          "Select a row in the table above to see that mouse's trajectory.")),
            shinycssloaders::withSpinner(
              plotOutput(ns("tg_drilldown_plot"), height = "350px")
            )
          ),
          
          tabPanel("Data Summary",
            tab_info(
              strong("Descriptive statistics for tumor volume "), "by treatment group and timepoint, on the volume scale in the declared units, whatever transform the model used.",
              tags$ul(
                tags$li(strong("Treatment"), " \u2014 treatment group."),
                tags$li(strong("Day"), " \u2014 study day (measurement timepoint)."),
                tags$li(strong("N"), " \u2014 number of animals with a volume measured on that day."),
                tags$li(strong("Mean"), " \u2014 arithmetic mean tumor volume for the group at that timepoint."),
                tags$li(strong("SD"), " \u2014 standard deviation of volume across animals in the group."),
                tags$li(strong("SEM"), " \u2014 standard error of the mean (SD / \u221aN)."),
                tags$li(strong("Geometric mean"), " (log transform only) \u2014 the mean on the log scale the model uses, back-transformed; zeros count at the floor stated in the Summary. Tumour volumes are right-skewed, so it is usually below the arithmetic mean.")
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
                    conditionalPanel(ns = ns,
                      "['growth_rate','auc','effect_forest'].indexOf(input.tg_plot_selector) >= 0",
                      checkboxInput(ns("tg_opt_backtransform"),
                                    "Display on original scale",
                                    value = FALSE),
                      helpText("Back-transforms log-scale model output to original units. ",
                               "Growth rates and contrasts become daily fold changes (exp(value)); ",
                               "forest plot uses a log10 axis centred on 1. ",
                               "Inference is still performed on the log scale.")
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
            tab_info(
              strong("Model-fit diagnostics. "),
              "Frequentist fits: residual diagnostics. Bayesian fits: MCMC ",
              "convergence + posterior diagnostics. The view switches based ",
              "on the selected model type."
            ),
            # ── Freq diagnostics (LME4 / AUC) ──────────────────────
            conditionalPanel(
              ns = ns,
              condition = "input.tg_model_type != 'bayes'",
            fluidRow(
              column(6,
                h4("Q-Q Plot of Residuals"),
                plot_help(
                    "The residuals are standardised (centred and divided by their SD), so ",
                    "for a well-fitting model the points fall along the red dashed line. ",
                    strong("S-shape: "), "residuals have heavier tails than normal — ",
                    "consider a robust likelihood or check for outliers. ",
                    strong("Bow: "), "skewed residuals — check the volume transform ",
                    "(log vs sqrt vs none). ",
                    strong("A few outliers at the ends: "), "common; not a problem ",
                    "unless they're clinically obvious errors."),
                div(style = "margin-bottom: 24px;",
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_diag_qq"), height = "380px")
                  )
                )
              ),
              column(6,
                h4("Residuals vs Fitted"),
                plot_help(
                    "Points should scatter randomly around the red zero line. ",
                    strong("Curved blue trend: "), "the model misses non-linearity ",
                    "— growth is not log-linear over the whole study; try a different volume ",
                    "transform. ",
                    strong("Fan shape (wider on the right): "), "heteroscedasticity ",
                    "— variance grows with the mean. Switch to log transform or a ",
                    "model with variance structure."),
                div(style = "margin-bottom: 24px;",
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_diag_resid_fitted"), height = "380px")
                  )
                )
              )
            ),
            br(),
            fluidRow(
              column(6,
                h4("Scale-Location"),
                plot_help(
                    "Plots sqrt(|residuals|) vs fitted to magnify variance ",
                    "trends. ",
                    strong("Flat blue line: "), "constant variance — LMM assumption ",
                    "satisfied. ",
                    strong("Upward slope: "), "heteroscedasticity — see ",
                    "Residuals vs Fitted notes above."),
                div(style = "margin-bottom: 24px;",
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_diag_scale_location"), height = "380px")
                  )
                )
              ),
              column(6,
                h4("Random-effects Q-Q"),
                plot_help(
                    "Per-animal random intercepts (and slopes, if used) should ",
                    "be normally distributed. The values are standardised, so points ",
                    "along the red line = assumption met. A component with no spread ",
                    "(estimated variance zero) is left out. Heavy tails / clear S-shape: the LMM ",
                    "may be giving a few animals too much / too little ",
                    "shrinkage. Rare and usually safe to ignore at N ≤ 30."),
                div(style = "margin-bottom: 24px;",
                  shinycssloaders::withSpinner(
                    plotOutput(ns("tg_diag_re_qq"), height = "380px")
                  )
                )
              )
            ),
            # CODE_REVIEW.md DIAGNOSTICS gap (13) — LMM influence
            # diagnostics. Only shown when diag_cooks_distance was
            # populated (LME4 path with include_diagnostics = TRUE).
            conditionalPanel(
              ns = ns,
              condition = "output.tg_show_influence === true",
              hr(),
              h4("Influence diagnostics (case-deletion)"),
              tab_info(
                strong("Cook's distance"), " measures the overall effect of ",
                "leaving out each animal, one row per animal (arm, ID, cage). ",
                "Values above the threshold (4 / number of animals) flag ",
                "animals that move the fit. Inspect them individually \u2014 ",
                "they are not automatically wrong, but they are driving the ",
                "result. ",
                strong("DFBETAS"), " breaks the influence down per fixed ",
                "effect: each column is one coefficient, each row one ",
                "animal. Values | > 1 | are conventionally flagged."
              ),
              fluidRow(
                column(6,
                  h5("Cook's distance"),
                  shinycssloaders::withSpinner(
                    DTOutput(ns("tg_cooks_distance"))
                  )
                ),
                column(6,
                  h5("DFBETAS"),
                  shinycssloaders::withSpinner(
                    DTOutput(ns("tg_dfbetas"))
                  )
                )
              )
            )
            ),  # /freq conditionalPanel
            # ── Bayes MCMC diagnostics ────────────────────────────────────
            conditionalPanel(
              ns = ns,
              condition = "input.tg_model_type == 'bayes'",
              tab_info(
                strong("MCMC convergence diagnostics"),
                " for all fixed-effect parameters. All chains should have ",
                "converged (R̂ ≈ 1.00) and sufficient effective ",
                "samples before interpreting results.",
                tags$ul(
                  tags$li(strong("Parameter"), " — model parameter name."),
                  tags$li(strong("Rhat"), " — potential scale-reduction factor; values > 1.01 indicate non-convergence."),
                  tags$li(strong("Bulk_ESS"), " — bulk effective sample size; aim for > 400 per chain."),
                  tags$li(strong("Tail_ESS"), " — tail effective sample size; important for credible interval accuracy."),
                  tags$li(strong("Converged"), " — TRUE when R̂ ≤ 1.01.")
                ),
                em("If Rhat > 1.01, try increasing iterations or chains before reporting results.")
              ),
              uiOutput(ns("tg_bayes_diagnostics")),
              shinycssloaders::withSpinner(DTOutput(ns("tg_mcmc_diag_table"))),
              br(),
              h5(strong("Influential observations (Pareto-k > 0.7)")),
              em(class = "text-muted",
                 "Observations the LOO approximation can't reliably estimate. Investigate these mice individually."),
              shinycssloaders::withSpinner(DTOutput(ns("tg_bayes_influential"))),
              br(),
              fluidRow(
                column(6,
                  h5(strong("Pareto-k distribution (all observations)")),
                  plot_help(
                      "Histogram of LOO Pareto-k values across all observations. ",
                      strong("All bars left of 0.5 (orange dotted): "), "LOO is reliable. ",
                      strong("Some bars between 0.5 and 0.7: "), "moderate influence; estimates usable but inspect. ",
                      strong("Bars right of 0.7 (red dashed): "), "those observations are too influential for PSIS-LOO to estimate reliably — see the table above for which."),
                  shinycssloaders::withSpinner(plotOutput(ns("tg_pareto_k_plot"), height = "300px"))
                ),
                column(6,
                  h5(strong("Posterior predictive coverage")),
                  plot_help(
                      "Each bar = empirical coverage of the corresponding nominal posterior-predictive interval. ",
                      strong("Bars near red dashed line: "), "well-calibrated. ",
                      strong("Bars below: "), "under-coverage — model intervals are too narrow. ",
                      strong("Bars above: "), "over-coverage — too wide; model is being conservative."),
                  shinycssloaders::withSpinner(plotOutput(ns("tg_ppc_coverage_plot"), height = "300px"))
                )
              ),
              br(),
              h5(strong("MCMC trace plots")),
              plot_help(
                "Each line is one chain's draws of a parameter, after warmup. Converged chains ",
                "overlap and look like stationary noise around one level (a \u2018fuzzy ",
                "caterpillar\u2019). ", strong("A drift, a trend, or a chain at a different level: "),
                "the chains have not converged and the estimates are unreliable; run more warmup ",
                "iterations or simplify the model."),
              shinycssloaders::withSpinner(plotOutput(ns("tg_mcmc_trace_plot"), height = "500px")),
              br(),
              h5(strong("MCMC rank plots")),
              plot_help(
                  "Modern convergence diagnostic. Each chain's draws are rank-normalised; chains that have converged look like uniform draws (flat histograms) and resemble each other. ",
                  strong("Distinctly non-uniform shapes or chain-specific patterns: "), "the chain hasn't mixed and the posterior estimate is unreliable.",
                  summary = "How to read these plots"),
              shinycssloaders::withSpinner(plotOutput(ns("tg_mcmc_rank_plot"), height = "500px")),
              br(),
              fluidRow(
                column(6,
                  h5(strong("Treatment-parameter posterior densities")),
                  plot_help(
                      "Posterior density (filled area) for each treatment coefficient. ",
                      "Look for: (a) modes shifted from the prior (the data informed the parameter), (b) credible-interval bars that don't include zero (significant posterior direction), (c) heavy tails / multiple modes (potential identifiability issue)."),
                  shinycssloaders::withSpinner(plotOutput(ns("tg_posterior_dist_plot"), height = "400px"))
                ),
                column(6,
                  h5(strong("Bayesian residuals vs day")),
                  plot_help(
                      "Posterior-mean residuals plotted against study day, faceted by treatment, with a loess smoother. Flat horizontal line at zero per group = the linear-time assumption is satisfied. ",
                      strong("Systematic curvature: "), "the log-linear time assumption is violated; read growth-rate comparisons with care."),
                  shinycssloaders::withSpinner(plotOutput(ns("tg_bayes_residuals_plot"), height = "400px"))
                )
              )
            )  # /bayes conditionalPanel
          ),

          # R20-D.X: replaced "Download All"; the user chooses what to export.
          tabPanel(
            title = tags$span(icon("file-zipper"), " Download Results"),
            value = "__download_results__",
            wellPanel(download_results_ui(ns, "tg"))
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

    # \u2500\u2500 Workflow overview \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    h3(icon("arrows-turn-right"), " Basic Workflow"),
    tags$ol(
      tags$li(strong("Data tab"), " \u2014 upload a CSV file and map your columns (ID, time, volume, treatment). Optional columns: cage, dose, body weight, survival/censor indicator."),
      tags$li(strong("Study Design"), " (on the Data tab) \u2014 declare the tumour volume units (mm\u00b3 or cm\u00b3), the control group, and, when a cage column is mapped, whether individual animals or whole cages were randomised. Every analysis tab uses these; none has its own reference selector."),
      tags$li(strong("Analysis tabs"), " \u2014 each tab activates once the required columns are mapped:"),
      tags$ul(
        tags$li(strong("Tumor Growth"), " \u2014 always available after upload"),
        tags$li(strong("Survival"), " \u2014 requires a survival/censor column (or use the volume-threshold derivation)"),
        tags$li(strong("Toxicity"), " \u2014 requires a body weight column"),
        tags$li(strong("Dose Response"), " \u2014 requires a dose column"),
        tags$li(strong("Drug Synergy"), " \u2014 requires control, single-agent A, single-agent B, and combination groups"),
        tags$li(strong("Power Analysis"), " \u2014 always available; no data upload needed for a priori modes")
      ),
      tags$li("Select the statistical ", strong("method"), " in the sidebar, configure parameters, and click ", strong("Run Analysis"), "."),
      tags$li("Explore results across sub-tabs (summary, plots, pairwise comparisons, model output)."),
      tags$li("Download plots (PNG) and tables (CSV) from the per-tab download buttons.")
    ),

    hr(),

    # \u2500\u2500 Data format \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    h3(icon("table"), " Data Format"),
    p("Upload data in ", strong("long format"), ": one row per measurement per animal. Accepted format: CSV (.csv)."),

    h4("Demo data"),
    p("No data of your own? The Data tab ships four demo datasets:"),
    tags$ul(
      tags$li(tags$strong("Combo Treatment"), " \u2014 four arms with survival data, for the Drug Synergy tab."),
      tags$li(tags$strong("Dose Levels"), " \u2014 graded doses with survival, for Dose Response."),
      tags$li(tags$strong("Weight"), " \u2014 four arms with body weight, for Toxicity."),
      tags$li(tags$strong("All Features"), " \u2014 the most complete, and the one that reaches every module: ",
              "six arms (vehicle, a three-dose series, a second drug and a combination), length/width/height, ",
              "body weight, dose, two cages per arm, necrosis flags, survival censoring and a removal-reason ",
              "column (tumour volume limit, weight loss, found dead, end of study).")
    ),

    h4("Column auto-detection"),
    p("The dashboard detects columns by name. Recognised patterns (case-insensitive):"),
    tags$table(
      class = "table table-condensed table-bordered",
      style = "max-width:680px;",
      tags$thead(
        tags$tr(tags$th("Role"), tags$th("Recognised names"), tags$th("Required?"))
      ),
      tags$tbody(
        tags$tr(tags$td("Animal ID"),         tags$td("ID, Mouse_ID, Animal_ID, Subject, Tag, Ear Tag, Mouse Tag"), tags$td("Yes")),
        tags$tr(tags$td("Time / Day"),        tags$td("Day, Time, Timepoint, Day_Number"),  tags$td("Yes \u2021")),
        tags$tr(tags$td("Date"),              tags$td("Date, Measurement_Date"),            tags$td("Yes \u2021")),
        tags$tr(tags$td("Tumor volume"),      tags$td("Volume, Tumor_Volume, TV"),          tags$td("Yes \u2020")),
        tags$tr(tags$td("Tumor length"),      tags$td("Length, Tumor_Length"),              tags$td("Yes \u2020")),
        tags$tr(tags$td("Tumor width"),       tags$td("Width, Tumor_Width"),                tags$td("Yes \u2020")),
        tags$tr(tags$td("Tumor height"),      tags$td("mapped manually in Data Upload"),    tags$td("No")),
        tags$tr(tags$td("Treatment / Group"), tags$td("Treatment, Group, Treatment_Group"), tags$td("Yes")),
        tags$tr(tags$td("Cage"),              tags$td("Cage, Cage_ID"),                     tags$td("No")),
        tags$tr(tags$td("Dose"),              tags$td("Dose, Dose_Level, Concentration, Conc"), tags$td("No")),
        tags$tr(tags$td("Body weight"),       tags$td("Weight, Mass, Body_Weight, Body Mass, BW"), tags$td("No")),
        tags$tr(tags$td("Sex"),               tags$td("Sex, Gender, M/F"),                  tags$td("No")),
        tags$tr(tags$td("Necrotic flag"),     tags$td("Necrotic, Necrosis, Necr, Is_Necrotic, Necrotic_Flag"), tags$td("No")),
        tags$tr(tags$td("Removal reason"),    tags$td("Removal_Reason, Reason, Cause_Of_Death, Fate"), tags$td("No")),
        tags$tr(tags$td("Survival / Censor"), tags$td("Survival_Censor, Censor, Event, Survival_Event"), tags$td("No")),
        tags$tr(tags$td("Survival time"),     tags$td("Survival_Time, Time_To_Event, Days_To_Event"), tags$td("No"))
      )
    ),
    p(tags$em("\u2020 Provide either a pre-calculated Volume column, or Length + Width columns."),
      " When volume is calculated from calipers, six formulas are selectable in Data Upload: ",
      "ellipsoid (default, ", tags$code("\u03c0\u00b7L\u00b7W\u00b2/6"), "), modified ellipsoid (",
      tags$code("L\u00b7W\u00b2/2"), ") \u2014 the common preclinical convention \u2014 plus 3-axis ellipsoid, cylinder, ",
      "sphere and box. The 3-axis and box formulas use a Height column if one is mapped; the others ",
      "approximate height by width. Volumes may be supplied in mm\u00b3 or cm\u00b3."),
    p(tags$em("\u2021 Supply either a numeric Day column or a Date column; elapsed days are derived from the date and a reference start date.")),
    p("Survival/Censor must be numeric: ", strong("1 = endpoint reached"), " (e.g. tumour exceeded IACUC volume limit), ",
      strong("0 = censored"), " (alive or study-end censored). If absent, use the volume-threshold derivation in the Survival tab sidebar."),

    hr(),

    # \u2500\u2500 Analysis modules \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    h3(icon("flask"), " Analysis Modules"),

    h4("Tumor Growth"),
    tags$dl(
      tags$dt("LME4 (linear mixed-effects)"),
      tags$dd("Log-transforms tumor volume and fits a linear mixed model with a Treatment \u00d7 Day interaction. Each animal has its own intercept and growth rate (random slopes), falling back to simpler structures only if the fit fails; the ANOVA reports F-tests with Satterthwaite degrees of freedom. Cage handling is chosen to match the housing design rather than applied uniformly \u2014 see the tab's Info panel. Pairwise contrasts via emmeans; the comparison family and its correction are two separate settings, defaulting to vs-reference with an exact Dunnett correction."),
      tags$dt("AUC (area under the curve)"),
      tags$dd("The area under each arm's fitted growth curve over the window in which every arm is evaluable (at least 50 % of its enrolled animals, and at least 3, on study). Animals that die or are removed early count through all their measurements, so early deaths of large tumours do not pull the AUC down. Arms are compared by the ratio of their AUCs, with intervals and Wald p-values from the model and the selected correction. Each animal's own trapezoidal AUC is shown for description only: an animal removed early accumulates less area from being measured for less time."),
      tags$dt("Bayesian LMM (brms)"),
      tags$dd("The LME4 model fitted via MCMC (Stan). Return posterior distributions for treatment effects with credible intervals, plus convergence diagnostics (Rhat, ESS, trace and rank plots) and PSIS-LOO influence measures. Requires the brms package."),
      tags$dt("Randomisation test (optional, any frequentist mode)"),
      tags$dd("Permutes treatment labels and refits, giving an exact p-value for the Treatment \u00d7 Day interaction that assumes neither a denominator-df approximation nor normal random effects. The unit permuted, individual animals or whole cages, is the unit of randomisation declared under Study Design on the Data tab \u2014 the test is only valid if that matches how the experiment was actually run.")
    ),

    h4("Survival Analysis"),
    tags$dl(
      tags$dt("Kaplan-Meier"),
      tags$dd("Non-parametric survival estimates with 95% confidence bands. Log-rank test for group differences."),
      tags$dt("Cox proportional hazards"),
      tags$dd("Semi-parametric Cox PH model. Hazard ratios with 95% CI for each treatment vs. reference. When whole cages were randomised (Study Design), each arm is instead compared with the reference by a cage-level permutation log-rank, and the hazard ratios have no interval. Firth's bias-reduced estimation is applied automatically when any arm has zero events, where the standard partial likelihood would return infinite coefficients. The proportional-hazards assumption is tested and reported \u2014 if it fails, the hazard ratio is an average over a period in which the true ratio changed, and should not be quoted as a single number."),
      tags$dt("Bayesian AFT (brms)"),
      tags$dd("Accelerated failure time model with a Weibull or log-normal distribution. Returns time ratios and derived hazard ratios with credible intervals.")
    ),

    h4("Toxicity (Body Weight)"),
    tags$dl(
      tags$dt("LME4 body weight model"),
      tags$dd("Linear mixed model on body weight over time. Treatment \u00d7 Day interaction identifies groups with significant weight-change trajectories."),
      tags$dt("Therapeutic window"),
      tags$dd("Sets each arm's tumor growth inhibition beside its mean worst weight loss (each animal's largest drop below its own first weighing), each with a 95% interval, and flags the weight loss as Tolerated, Not tolerated or Unclear against a threshold you set (default 20%). The two axes are read separately; a single efficacy-toxicity ratio was reported before version 26.09.005 and was withdrawn."),
      tags$dt("Weight-loss time-to-event"),
      tags$dd("Kaplan-Meier for time to crossing a weight-loss threshold. With a removal-reason column mapped on Data Upload, removals for weight loss count as events, planned ends as censoring, and other removals (such as tumour burden) as competing events, with an Aalen-Johansen cumulative incidence estimate. Without one, an animal whose record ends early is censored, and the tab says so. A per-arm breakdown of how animals left is reported."),
      tags$dt("Endpoint estimand"),
      tags$dd("Shared by the Toxicity and Synergy tabs. \u2018Model\u2019 estimates endpoint means from the fitted model using every observation; \u2018Observed\u2019 uses only animals still present at the endpoint. They differ whenever animals leave early, and because the fastest-growing controls are the ones that leave, \u2018Observed\u2019 flatters the control arm.")
    ),

    h4("Drug Synergy"),
    tags$dl(
      tags$dt("Bliss Independence"),
      tags$dd("Calculates the expected additive effect from two monotherapies at a chosen time point, then tests whether the observed combination effect exceeds it. Positive excess = synergy; negative excess = antagonism."),
      tags$dt("Why there is no Combination Index"),
      tags$dd("A Loewe-style CI was removed in v0.21.0. What it computed was not ",
              "Loewe additivity \u2014 true Loewe is dose-equivalence and needs a ",
              "dose-response curve per agent, which single-dose designs do not ",
              "provide. The stand-in was ", tags$em("response additivity"),
              ", which fails the sham-combination test: combining a drug with ",
              "itself it predicts twice the effect, so an agent at 50 % inhibition ",
              "should reach 100 %, and when it does not the method calls the drug ",
              "antagonistic with itself. Measured across the effect grid with the ",
              "combination set to exactly Bliss-additive, 42 % of cells were ",
              "labelled antagonistic and none synergistic."),
      tags$dt("Interval on the verdict"),
      tags$dd("The synergy verdict comes from the 95 % interval for the Bliss excess: Synergy when it lies above zero, Antagonism when it lies below, and Additive (no departure from Bliss detected) when it contains zero. Without intervals the label says so.")
    ),

    h4("Dose-Response"),
    tags$dl(
      tags$dt("One agent's dose series"),
      tags$dd("A curve describes one agent at several doses against the study's control. Choose the agent's arms in the sidebar; other agents, combinations, and arms sharing a dose cannot join the same curve."),
      tags$dt("Hill / Emax"),
      tags$dd("Fits a four-parameter log-logistic (4PL) dose-response curve via nonlinear least squares. Reports the EC50 (the dose giving half the fitted response range) with a 95% interval on log dose, the Hill slope, and the fitted curve; an EC50 outside the tested doses is flagged."),
      tags$dt("Trend tests"),
      tags$dd("Alongside the curve fit, the tab reports linear and quadratic trend p-values and a two-sided Jonckheere-Terpstra test. The last is the one that respects dose ordering: it asks whether response is monotone in dose, which is the question a dose-response study is usually posing, rather than merely whether the groups differ.")
    ),

    h4("Power Analysis"),
    tags$dl(
      tags$dt("Analytic (a priori)"),
      tags$dd("Classical power for each treated-vs-control comparison (two-sample t-test), with the alpha split across the comparisons a multi-arm study makes and an optional attrition allowance, across a range of effect sizes and sample sizes."),
      tags$dt("LMM simulation"),
      tags$dd("Simulates repeated-measures data under a linear mixed model and estimates power empirically across a grid of sample sizes. An optional euthanasia volume limit models attrition; leaving it blank assumes every animal is measured at every timepoint, which overstates power by roughly 6 points at 30 % attrition. Note that baseline SD and random-intercept SD do not affect the answer \u2014 the test is on growth rates, which are orthogonal to per-animal baselines. Random-slope SD is the parameter that matters.")
    ),

    hr(),

    # \u2500\u2500 Troubleshooting \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    h3(icon("wrench"), " Troubleshooting"),
    tags$dl(
      tags$dt("File won't upload"),
      tags$dd("Ensure the file is a CSV (.csv) file."),
      tags$dt("Columns not auto-detected"),
      tags$dd("Use the column mapping dropdowns in the Data tab to assign columns manually."),
      tags$dt("Analysis tab is greyed out"),
      tags$dd("Upload and configure data first. Survival requires a censor column (or enable volume-threshold derivation). Toxicity requires a body weight column. Dose Response requires a dose column. Drug Synergy requires at least four treatment groups."),
      tags$dt("Analysis fails or produces no output"),
      tags$dd("Verify that the selected columns contain valid, non-missing data. Tumor volumes must be positive. Time points must be numeric. At least two groups must be present for any comparison."),
      tags$dt("Bayesian model is slow to start"),
      tags$dd("The first run compiles the Stan model (\u223c30\u201360 s). Subsequent runs with the same model family reuse the compiled binary and start faster. Bayesian fits run asynchronously \u2014 the app remains usable while MCMC sampling proceeds."),
      tags$dt("Plot download produces a blank file"),
      tags$dd("Run the analysis first so the plot is rendered, then click the download button.")
    )
  )
}

#' About UI
#' @keywords internal
citationUI <- function() {
  pkg_version   <- tryCatch(as.character(utils::packageVersion("mouseExperiment")),      error = function(e) "(not installed)")
  dash_version  <- tryCatch(as.character(utils::packageVersion("mouseExperimentDashboard")), error = function(e) "(not installed)")
  # Neither package has a Date field any more (it went stale; R20.D26), and
  # the fallback here was a fixed 2026-04-03. The build date is always there.
  built_date <- function(pkg) tryCatch({
    b <- utils::packageDescription(pkg)$Built
    d <- if (is.null(b)) NA_character_ else regmatches(b, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", b))
    if (length(d) && !is.na(d)) d else "unknown"
  }, error = function(e) "unknown")
  pkg_date  <- built_date("mouseExperiment")
  dash_date <- built_date("mouseExperimentDashboard")

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
                tags$td(strong("Built")),
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
                tags$td(strong("Built")),
                tags$td(dash_date)
              )
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        wellPanel(
          h4(icon("envelope"), " Contact"),
          p("Questions, bug reports and feature requests: ",
            tags$a("david@sciomics.bio", href = "mailto:david@sciomics.bio"), ".")
        )
      )
    ),

    # R20.D12: say what is counted.
    fluidRow(
      column(12,
        wellPanel(
          h4(icon("chart-simple"), " Usage analytics"),
          p("mouseexperiment.com counts page views and which features are used, with ",
            "self-hosted Umami analytics: no cookies, and nothing sent to a third party. ",
            "Events carry only the analysis type, the method, counts such as the number ",
            "of rows, and yes/no flags. They never include your data, file names, column ",
            "names, treatment names or error messages. Uploaded data stays in your ",
            "session on this server.")
        )
      )
    )
  )
}
