#' UI Modules for mouseExperimentDashboard

# Shared helper: creates a styled info box for tab descriptions
tab_info <- function(...) {
  div(
    class = "alert alert-info",
    style = "font-size:0.88em; margin:10px 0 8px 0; padding:8px 14px;",
    ...
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
dataUploadUI <- function() {
  fluidPage(
    h2("Data Upload & Configuration"),

    fluidRow(
      column(12,
        wellPanel(
          fluidRow(
            column(4,
              h5(strong("Try Demo Data")),
              p("Load a built-in demo dataset to explore the dashboard's features."),
              selectInput("demo_dataset", NULL,
                choices = c(
                  "-- Select a demo dataset --" = "",
                  "Combo Treatment (Control / HDACi / aPD1 / HDACi+PD1, survival data)" = "combo",
                  "Dose Levels (dose-response data with survival)" = "dose_levels"
                ),
                selected = ""
              ),
              fluidRow(
                column(6,
                  actionButton("load_demo_data", "Load Demo Data",
                               icon = icon("play-circle"),
                               class = "btn-default btn-block")
                ),
                column(6, uiOutput("demo_load_msg"))
              )
            ),
            column(8,
              fileInput("data_file",
                       "Choose File (CSV or RDA)",
                       accept = c(".csv", ".rda", ".RData")),
              uiOutput("file_validation_msg")
            )
          )
        )
      )
    ),
    
    # Column mapping section (shown after file upload)
    uiOutput("column_mapping_ui"),
    
    # Dataset summary (shown after successful upload + column config)
    uiOutput("data_summary_ui"),
    
    # Data preview (shown after successful upload)
    uiOutput("data_preview_ui"),
    
  )
}

#' Column Mapping UI
#' @keywords internal
column_mapping_ui <- function(col_names, detected_cols) {
  wellPanel(
    h4(strong("Configure Column Mapping")),
    tags$p("Select which columns contain your data. ", 
      "Suggested columns are pre-selected."),
    
    fluidRow(
      column(4,
        h4(strong("Required Columns")),
        selectInput("col_id", "Mouse/Subject ID",
                   choices = c("Select..." = "", col_names),
                   selected = detected_cols$id),
        
        radioButtons("time_method", "How to specify time?",
                    choices = c(
                      "Direct Day column" = "direct_day",
                      "Calculate from Date column" = "from_date"
                    ),
                    selected = if (!is.null(detected_cols$day) && detected_cols$day != "") "direct_day" else "from_date"),
        
        conditionalPanel(
          condition = "input.time_method == 'direct_day'",
          selectInput("col_day", "Day Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$day)
        ),
        
        conditionalPanel(
          condition = "input.time_method == 'from_date'",
          selectInput("col_date_for_day", "Date Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$date),
          textInput("reference_date", "Reference Date (YYYY-MM-DD)",
                   value = "2025-03-24",
                   placeholder = "First measurement date")
        ),
        
        selectInput("col_treatment", "Treatment Group",
                   choices = c("Select..." = "", col_names),
                   selected = detected_cols$treatment)
      ),
      
      column(4,
        h4(strong("Volume Measurement")),
        radioButtons("volume_method", "How to specify volume?",
                    choices = c(
                      "Direct volume column" = "direct",
                      "Calculate from Length & Width" = "calculate"
                    ),
                    selected = if (!is.null(detected_cols$volume)) "direct" else "calculate"),
        
        conditionalPanel(
          condition = "input.volume_method == 'direct'",
          selectInput("col_volume", "Volume Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$volume)
        ),
        
        conditionalPanel(
          condition = "input.volume_method == 'calculate'",
          selectInput("col_length", "Length Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$length),
          selectInput("col_width", "Width Column",
                     choices = c("Select..." = "", col_names),
                     selected = detected_cols$width),
          selectInput("volume_formula", "Volume Formula",
                     choices = c("Ellipsoid (default)" = "ellipsoid",
                                "Modified Ellipsoid" = "modified_ellipsoid",
                                "Cylinder" = "cylinder",
                                "Sphere" = "sphere"),
                     selected = "ellipsoid"),
          conditionalPanel("input.volume_formula == 'ellipsoid'",
            tags$p(tags$small(HTML("V = (&pi;/6) &times; L &times; W&sup2;")))
          ),
          conditionalPanel("input.volume_formula == 'modified_ellipsoid'",
            tags$p(tags$small(HTML("V = (L &times; W&sup2;) / 2")))
          ),
          conditionalPanel("input.volume_formula == 'cylinder'",
            tags$p(tags$small(HTML("V = &pi; &times; (W/2)&sup2; &times; L")))
          ),
          conditionalPanel("input.volume_formula == 'sphere'",
            tags$p(tags$small(HTML("V = (&pi;/6) &times; ((L + W) / 2)&sup3;")))
          )
        )
      ),
      
      column(4,
        h4(strong("Optional Columns")),
        selectInput("col_cage", "Cage ID",
                   choices = c("None" = "", col_names),
                   selected = detected_cols$cage),
        tags$p(tags$small("Cage ID is used to distinguish mice with same ID in different cages.")),
        
        br(),
        
        selectInput("col_dose", "Dose Level",
                   choices = c("None" = "", col_names),
                   selected = detected_cols$dose),
        tags$p(tags$small("Required for Dose-Response Analysis. Should contain numeric dose concentrations or levels.")),
        
        br(),
        selectInput("col_survival_event", "Event/Censor Status Column",
                   choices = c("None" = "", col_names),
                   selected = detected_cols$survival_event),
        tags$p(tags$small("Required for Survival Analysis. Time will use the same Day/Date settings as above."))
      )
    )
  )
}

#' Data Preview UI
#' @keywords internal
data_preview_ui <- function() {
  wellPanel(
    h4(strong("Data Preview")),
    p(tags$small("Shows calculated columns (Volume, Day from start date) where applicable.")),
    fluidRow(
      column(12,
        p(strong(textOutput("preview_summary"))),
        DTOutput("data_table"),
        br(),
        downloadButton("download_preview_data", "Download Processed Data as CSV",
                      class = "btn-sm btn-default")
      )
    )
  )
}

#' Tumor Growth Analysis UI
#' @keywords internal
tumorGrowthUI <- function() {
  fluidPage(
    h2("Tumor Growth Analysis"),
    
    p("Configure analysis parameters and run statistical modeling. ",
      "Ensure data is uploaded first from the Data Upload tab."),
    
    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Analysis Parameters"),

          # 1. Statistical Model
          selectInput("tg_model_type",
                     "Statistical Model",
                     choices = c("Linear Mixed-Effects (LME4)" = "lme4",
                                "Area Under Curve (AUC)" = "auc"),
                     selected = "lme4"),
          helpText("LME4 models tumor volume trajectories over time using a linear ",
                   "mixed-effects model, which handles repeated measures and unequal ",
                   "timepoints. AUC summarizes total tumor burden as the area under ",
                   "each animal's growth curve."),

          # 2. Reference Group
          selectInput("tg_reference",
                     "Reference Group (Required)",
                     choices = c("-- Select a reference group --" = "")),
          helpText("The group all comparisons are made against (typically a vehicle or ",
                   "IgG control). Automatically set to the most likely control group ",
                   "when data is uploaded."),

          # 3. Post-hoc Comparisons
          selectInput("tg_comparison_type",
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
          conditionalPanel(
            condition = "input.tg_model_type == 'lme4'",
            selectInput("tg_transform",
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
          numericInput("tg_extrapolation_points",
                      "Extrapolation Points",
                      value = 3, min = 0, max = 10, step = 1),
          helpText("Number of recent timepoints used to estimate the growth slope for ",
                   "extrapolating tumor volumes of mice that died before the study end. ",
                   "Set to 0 to disable. If an Event/Censor column is mapped in Data ",
                   "Upload it will be used; otherwise mice that disappear before the ",
                   "final timepoint are automatically detected."),

          # 6. Cage Effect Handling (LME4 only)
          conditionalPanel(
            condition = "input.tg_model_type == 'lme4'",
            selectInput("tg_cage_effect",
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

          # 7. Random Effects (LME4 only)
          conditionalPanel(
            condition = "input.tg_model_type == 'lme4'",
            selectInput("tg_random_effects",
                       "Random Effects",
                       choices = c("Intercept only" = "intercept_only",
                                  "Intercept + Slope" = "slope",
                                  "None" = "none"),
                       selected = "intercept_only"),
            helpText("Specifies the random effects structure for individual animals. ",
                     "'Intercept only' allows each animal to have its own baseline ",
                     "volume. 'Intercept + Slope' also allows individual growth rates, ",
                     "which is more flexible but may fail to converge with small samples.")
          ),

          br(),

          actionButton("run_tg_analysis",
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),

          uiOutput("tg_status_msg"),

          plot_export_ui("tg")
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
          id = "tg_results_tabs",
          
          tabPanel("Summary",
            tab_info(strong("Model summary: "), "Reports the model type (LME4 or AUC), data transformation applied, polynomial degree, random effects structure, number of subjects and treatment groups, reference group, and any convergence warnings."),
            shinycssloaders::withSpinner(
              verbatimTextOutput("tg_summary")
            )
          ),
          
          tabPanel("ANOVA",
            tab_info(strong("Type III ANOVA "), "from the linear mixed-effects model. Tests whether treatment, time, and/or their interaction significantly explain variation in tumor volume. ",
              strong("p\u202f<\u202f0.05"), " indicates a statistically significant effect. The ", em("Treatment:Time"), " interaction term is the most biologically informative."),
            shinycssloaders::withSpinner(
              DTOutput("tg_anova_table")
            ),
            downloadButton("download_tg_anova", "Download ANOVA")
          ),
          
          tabPanel("Comparisons",
            tab_info(strong("Post-hoc pairwise contrasts "), "between treatment groups based on estimated marginal means. ",
              br(),
              strong("vs Reference: "), "each group vs the baseline. ",
              strong("All Pairwise: "), "unadjusted contrasts between every pair. ",
              strong("Tukey HSD: "), "all pairs with familywise error correction. ",
              br(),
              "\u2018estimate\u2019 is the mean difference on the transformed scale; ",
              "\u2018p.value\u2019\u202f<\u202f0.05 indicates a significant difference between groups."),
            uiOutput("tg_comparisons_title"),
            shinycssloaders::withSpinner(
              DTOutput("tg_comparisons_table")
            ),
            downloadButton("download_tg_comparisons", "Download Comparisons")
          ),
          
          tabPanel("Treatment Effects",
            tab_info(strong("Estimated marginal means "), "(adjusted group means) from the model. Reflects the model-estimated volume at each timepoint for each treatment group, controlling for cage effects and individual variability. Provides the key effect sizes for biological interpretation."),
            shinycssloaders::withSpinner(
              DTOutput("tg_effects_table")
            ),
            downloadButton("download_tg_effects", "Download Effects")
          ),
          
          tabPanel("Growth Rates",
            tab_info(strong("Per-animal growth rates "), "calculated from individual tumor volume trajectories. ",
              strong("growth_rate"), ": slope of log\u2081\u209a(volume) over time (log-scale units/day). ",
              strong("daily_fold_change"), ": exp(growth_rate) — the multiplicative volume change per day (e.g. 1.10 = 10% growth/day). ",
              "Used to assess within-group variability and for power analysis effect size estimation. Not the primary model output \u2014 use Treatment Effects for group-level inference."),
            shinycssloaders::withSpinner(
              DTOutput("tg_growth_rates_table")
            ),
            downloadButton("download_tg_rates", "Download Growth Rates")
          ),
          
          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  selectInput("tg_plot_selector", "Select Plot",
                    choices = c(
                      "Tumor Growth Curves"      = "growth_curves",
                      "Growth Rate Distribution" = "growth_rate",
                      "AUC Comparison"           = "auc",
                      "Effect Size Forest Plot"  = "effect_forest"
                    )
                  ),
                  hr(),
                  h6(icon("sort"), " Group Order"),
                  uiOutput("tg_group_order_ui"),
                  hr(),
                  h6(icon("palette"), " Appearance"),
                  uiOutput("tg_color_pickers_ui"),
                  hr(),
                  h6(icon("sliders-h"), " Options"),
                  sliderInput("tg_opt_text_size", "Text size (pt)",
                              min = 8, max = 20, value = 12, step = 1),
                  selectInput("tg_opt_legend_pos", "Legend position",
                    choices = c("Auto" = "auto", "None" = "none",
                                "Top" = "top", "Bottom" = "bottom",
                                "Left" = "left", "Right" = "right",
                                "In plot - top left"  = "inplot_tl",
                                "In plot - top right" = "inplot_tr"),
                    selected = "auto"
                  ),
                  conditionalPanel(
                    "input.tg_plot_selector == 'growth_curves'",
                    checkboxInput("tg_opt_group_line",
                                  "Show group average line", value = TRUE),
                    checkboxInput("tg_opt_show_extrap",
                                  "Show extrapolated curves/points", value = FALSE),
                    selectInput("tg_opt_yscale", "Y-axis scale",
                      choices = c("Linear" = "linear",
                                  "Log10"  = "log10"),
                      selected = "linear"
                    ),
                    sliderInput("tg_opt_growth_pt_size", "Point size",
                                min = 0.5, max = 8, value = 2, step = 0.5),
                    numericInput("tg_opt_xstart", "X-axis start (Day)",
                                 value = 0, min = 0, step = 1)
                  ),
                  conditionalPanel(
                    "['growth_rate','auc'].indexOf(input.tg_plot_selector) >= 0",
                    selectInput("tg_opt_center_bar", "Center line / bar",
                      choices = c("Mean" = "mean", "Median" = "median",
                                  "Geometric Mean" = "geometric_mean",
                                  "None" = "none"),
                      selected = "mean"
                    ),
                    selectInput("tg_opt_error_bar", "Error bars",
                      choices  = c("None" = "none", "SEM" = "SEM",
                                   "SD" = "SD", "95% CI" = "CI"),
                      selected = "SEM"
                    ),
                    sliderInput("tg_opt_point_size",   "Point size",
                                0.5, 8, 2.5, step = 0.5),
                    sliderInput("tg_opt_jitter_width", "Jitter width",
                                0, 1, 0.2, step = 0.05)
                  )
                )
              ),
              column(9,
                conditionalPanel("input.tg_plot_selector == 'growth_curves'",
                  h4("Tumor Growth Curves"),
                  tab_info(strong("Tumor Growth Curves: "),
                    "individual mouse trajectories with optional group mean overlay."),
                  shinycssloaders::withSpinner(
                    plotOutput("tg_plot_growth", height = "550px")
                  ),
                  downloadButton("download_tg_plot_growth", "Download Plot")
                ),
                conditionalPanel("input.tg_plot_selector == 'growth_rate'",
                  h4("Growth Rate Distribution"),
                  tab_info(strong("Growth Rate Distribution: "),
                    "spread of per-animal growth rates within each group."),
                  shinycssloaders::withSpinner(
                    plotOutput("tg_plot_rates", height = "550px")
                  ),
                  downloadButton("download_tg_plot_rates", "Download Plot")
                ),
                conditionalPanel("input.tg_plot_selector == 'auc'",
                  h4("AUC Comparison"),
                  tab_info(strong("AUC / Treatment Effects: "),
                    "model-estimated area under the curve, one point per subject."),
                  shinycssloaders::withSpinner(
                    plotOutput("tg_plot_auc", height = "550px")
                  ),
                  downloadButton("download_tg_plot_auc", "Download Plot")
                ),
                conditionalPanel("input.tg_plot_selector == 'effect_forest'",
                  h4("Effect Size Forest Plot"),
                  tab_info(strong("Pairwise effect sizes: "),
                    "mean differences between treatment groups with 95% CIs. ",
                    "Points left of zero indicate lower volume vs. reference; ",
                    "points to the right indicate higher volume."),
                  shinycssloaders::withSpinner(
                    plotOutput("tg_plot_effect_forest", height = "550px")
                  ),
                  downloadButton("download_tg_plot_effect_forest", "Download Plot")
                ),

              )
            )
          ),
          
          tabPanel("Diagnostics",
            tab_info(
              strong("Q\u2011Q Plot: "), "residuals should follow the diagonal for valid inference. ",
              strong("Residuals vs Fitted: "), "random scatter around zero indicates a good model fit; patterns suggest model misspecification or heteroscedasticity."
            ),
            fluidRow(
              column(6,
                h4("Q-Q Plot of Residuals"),
                shinycssloaders::withSpinner(
                  plotOutput("tg_diag_qq", height = "400px")
                )
              ),
              column(6,
                h4("Residuals vs Fitted"),
                shinycssloaders::withSpinner(
                  plotOutput("tg_diag_resid_fitted", height = "400px")
                )
              )
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
      tags$li(strong("Survival/Censor (optional):"), "Survival_Censor, Censor, Event")
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
    p(strong("Analysis tab greyed out:"), "Upload and configure data first. Survival requires a censor column; Dose Response requires a dose column."),
    p(strong("Analysis fails:"), "Check that required columns are selected and contain valid data.")
  )
}

#' Citation UI
#' @keywords internal
citationUI <- function() {
  fluidPage(
    h2("Citation"),

    wellPanel(
      h3("How to Cite"),
      p("If you use mouseExperiment or this dashboard in your research, please cite:"),
      tags$blockquote(
        style = "font-size: 1em; border-left: 4px solid #0073B7; padding-left: 16px;",
        p(strong("mouseExperiment: "), "An R package for statistical analysis of mouse tumor growth experiments."),
        p("Available at: ",
          tags$a("https://github.com/sciOmics/mouseExperiment",
                 href = "https://github.com/sciOmics/mouseExperiment", target = "_blank"))
      )
    ),

    wellPanel(
      h3("R Citation"),
      p("To get the citation information directly in R:"),
      tags$pre(
        style = "background:#f5f5f5; padding:12px; border-radius:4px;",
        'citation("mouseExperiment")'
      )
    ),

    wellPanel(
      h3("License"),
      p("mouseExperiment is released under the ",
        tags$a("MIT License", href = "https://github.com/sciOmics/mouseExperiment/blob/main/LICENSE",
               target = "_blank"), ".")
    )
  )
}
