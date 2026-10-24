# The Data tab's UI: upload, column mapping, Study Design, the data check and
# the preview (split from ui_modules.R, CODE_REVIEW.md R22-K).

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
                tags$strong("long format"), " \u2014 one row per observation/measurement."
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

    # F5: what the data check found, once the columns are mapped.
    uiOutput(ns("data_check_ui")),

    # Subset filter (CODE_REVIEW.md Round 2 F.5) — exclude mice / cages /
    # treatments. Sits between the Dataset Summary and Data Preview so
    # users see the summary first, choose what to filter, click Apply,
    # then the Data Preview below updates to reflect what they applied.
    subset_filter_ui(ns),

    # Data preview (shown after successful upload)
    uiOutput(ns("data_preview_ui")),

    # Settings preset (column mapping + methods)
    settings_preset_ui(ns),

    # F4: one report of every analysis run in this session.
    wellPanel(
      h5(strong(icon("file-lines"), " Study report")),
      tags$p(style = "font-size:12px; color:#6c757d;",
             "Every analysis you have run in this session, with its settings, tables, plots and ",
             "diagnostics, and the data check, in one HTML file to print or keep. Each tab's ",
             "Download Results has its own report too."),
      downloadButton(ns("download_study_report"), "Download the study report (HTML)",
                     class = "btn-sm btn-default"))
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
          uiOutput(ns("cage_labels_ui")),
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

#' Study Design's cage-label question (R22.3)
#'
#' @param reused The labels every arm reuses (see `reused_cage_labels()`);
#'   empty when the labels show no such pattern.
#' @param selected "study" or "arm" from a preset; otherwise guessed: "arm"
#'   when the arms reuse one set of labels.
#' @param ns Namespace function.
#' @keywords internal
cage_labels_ui <- function(reused = character(0), selected = NULL, ns = NS(NULL)) {
  if (length(selected) != 1L || !selected %in% c("study", "arm")) {
    selected <- if (length(reused)) "arm" else "study"
  }
  tagList(
    radioButtons(ns("cage_labels"),
      tagList("Cage labels are", help_tip("Unique across the study: each label is one cage, which may hold several arms. Numbered within each arm: C1 in Vehicle and C1 in Drug A are different cages, and the analyses label them 'Vehicle C1' and 'Drug A C1'. This decides how the cage effect is modelled: one cage holding several arms is estimated across arms.")),
      choices = c("Unique across the study" = "study", "Numbered within each arm" = "arm"),
      selected = selected),
    uiOutput(ns("cage_labels_status"))
  )
}

#' What the cage labels suggest, beside the question
#' @param reused The labels every arm reuses.
#' @param choice The current answer, "study" or "arm".
#' @keywords internal
cage_labels_note <- function(reused, choice) {
  shown <- paste0(paste(utils::head(reused, 5L), collapse = ", "),
                  if (length(reused) > 5L) ", ...")
  if (identical(choice, "arm")) {
    tags$p(class = "text-muted", style = "font-size:12px; margin-top:-6px;",
           icon("circle-info"),
           paste0(" Every arm uses the labels ", shown, ", so they are read as numbered ",
                  "within each arm. Choose 'Unique' only if those cages really held ",
                  "several arms."))
  } else {
    tags$p(class = "text-warning", style = "font-size:12px; margin-top:-6px;",
           icon("triangle-exclamation"),
           paste0(" Every arm uses the labels ", shown, ". Read as unique, each is one ",
                  "cage holding several arms, and the cage effect is estimated across ",
                  "arms. Choose 'Numbered within each arm' unless the cages really were ",
                  "shared."))
  }
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
       " \u2014 what each value of ", tags$code(col), " means")),
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

#' The Data tab's data-check panel (F5)
#'
#' @param qc The result of `check_study_data()`.
#' @param ns Namespace function.
#' @keywords internal
data_check_panel <- function(qc, ns = NS(NULL)) {
  s <- qc$summary
  colour <- c(Problem = "#c62828", Check = "#b26a00", Note = "#6c757d")
  n_prob <- sum(s$Severity == "Problem")
  n_check <- sum(s$Severity == "Check")
  headline <- if (!nrow(s) || (n_prob + n_check) == 0L) {
    status_alert(paste0("No problems found", if (nrow(s)) "; notes on the design below." else "."),
                 "success")
  } else if (n_prob) {
    status_alert(paste0("Possible errors in the data: ", n_prob, " kind", if (n_prob > 1L) "s",
                        " of problem. Correct them in your file, or exclude the rows' animals ",
                        "under Subset Filter, before you rely on an analysis."), "danger")
  } else {
    status_alert("Nothing that is clearly an error; a few things are worth a look.", "warning")
  }
  wellPanel(
    h4(strong("Data check")),
    tags$p(class = "text-muted", style = "font-size:12px;",
           "Your whole file, as uploaded, with the columns mapped above. Nothing is changed. ",
           "\"Line\" is the line in your file, the header being line 1."),
    headline,
    if (nrow(s)) tags$table(class = "table table-condensed", style = "margin-bottom:8px;",
      tags$thead(tags$tr(tags$th(""), tags$th("Check"), tags$th("Animals"), tags$th("Rows"),
                         tags$th("What it means"))),
      tags$tbody(lapply(seq_len(nrow(s)), function(i) tags$tr(
        tags$td(tags$span(style = paste0("color:", colour[[s$Severity[i]]], "; font-weight:bold;"),
                          s$Severity[i])),
        tags$td(s$Check[i]), tags$td(s$Animals[i]), tags$td(s$Rows[i]), tags$td(s$What[i]))))),
    if (nrow(qc$flagged)) tags$details(
      tags$summary("Show the flagged rows"),
      DT::DTOutput(ns("data_check_flagged")),
      downloadButton(ns("download_data_check"), "Download the flagged rows (CSV)",
                     class = "btn-sm btn-default"))
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
        DT::DTOutput(ns("data_table")),
        br(),
        downloadButton(ns("download_preview_data"), "Download Processed Data as CSV",
                      class = "btn-sm btn-default")
      )
    )
  )
}
