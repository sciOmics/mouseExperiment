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
        p("Upload your tumor growth data (CSV or RDA format). ",
          "The dashboard will automatically detect column names."),
        br()
      )
    ),
    
    fluidRow(
      column(6,
        wellPanel(
          h4("Step 1: Upload Tumor Growth Data"),
          fileInput("data_file", 
                   "Choose File (CSV or RDA)",
                   accept = c(".csv", ".rda", ".RData"))
        ),
        wellPanel(
          h4("Step 1b: Upload Treatment Schedule (Optional)"),
          p(tags$small("Upload a CSV/RDA file describing the treatment timeline (e.g. dosing days).",
                       "Used to overlay treatment events on tumor growth plots.")),
          fileInput("schedule_file",
                   "Treatment Schedule File (CSV or RDA)",
                   accept = c(".csv", ".rda", ".RData")),
          uiOutput("schedule_validation_msg")
        )
      ),
      
      column(6,
        wellPanel(
          h4("File Information"),
          verbatimTextOutput("file_info"),
          uiOutput("file_validation_msg")
        )
      )
    ),
    
    hr(),
    
    # Column mapping section (shown after file upload)
    uiOutput("column_mapping_ui"),
    
    hr(),
    
    # Dataset summary (shown after successful upload + column config)
    uiOutput("data_summary_ui"),
    
    # Data preview (shown after successful upload)
    uiOutput("data_preview_ui"),
    
    # Treatment schedule preview (shown after schedule upload)
    uiOutput("schedule_preview_ui"),
    
    hr(),
    
    # Action buttons
    fluidRow(
      column(12,
        actionButton("reset_data", "Start Over", class = "btn-warning"),
        actionButton("proceed_analysis", "Proceed to Analysis", 
                    class = "btn-success")
      )
    )
  )
}

#' Column Mapping UI
#' @keywords internal
column_mapping_ui <- function(col_names, detected_cols) {
  wellPanel(
    h4("Step 2: Configure Column Mapping"),
    tags$p("Select which columns contain your data. ", 
      "Suggested columns are pre-selected."),
    
    fluidRow(
      column(4,
        h4(strong("Required Columns")),
        selectInput("col_id", "Mouse/Subject ID",
                   choices = c("Select..." = "", col_names),
                   selected = detected_cols$id),
        
        h5("Time Points"),
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
                     selected = "ellipsoid")
        )
      ),
      
      column(4,
        h4(strong("Optional Columns")),
        selectInput("col_cage", "Cage ID",
                   choices = c("None" = "", col_names),
                   selected = detected_cols$cage),
        tags$p(tags$small("Cage ID is important for distinguishing mice with same ID in different cages.")),
        
        br(),
        
        selectInput("col_dose", "Dose Level (optional)",
                   choices = c("None" = "", col_names),
                   selected = detected_cols$dose),
        tags$p(tags$small("Required for Dose-Response Analysis. Should contain numeric dose concentrations or levels.")),
        
        br(),
        h5("Survival Analysis"),
        selectInput("col_survival_event", "Event/Censor Status Column",
                   choices = c("None" = "", col_names),
                   selected = detected_cols$survival_event),
        tags$p(tags$small("Optional: Only needed if performing survival analysis. Time will use the same Day/Date settings as above."))
      )
    )
  )
}

#' Data Preview UI
#' @keywords internal
data_preview_ui <- function() {
  wellPanel(
    h4("Step 3: Data Preview"),
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
    
    fluidRow(
      column(3,
        wellPanel(
          h4("Analysis Parameters"),
          
          selectInput("tg_model_type",
                     "Statistical Model",
                     choices = c("Linear Mixed-Effects (LME4)" = "lme4",
                                "Area Under Curve (AUC)" = "auc"),
                     selected = "lme4"),
          
          conditionalPanel(
            condition = "input.tg_model_type == 'lme4'",
            selectInput("tg_transform",
                       "Data Transformation",
                       choices = c("Log" = "log",
                                  "Square Root" = "sqrt",
                                  "None" = "none"),
                       selected = "log"),
            
            selectInput("tg_random_effects",
                       "Random Effects",
                       choices = c("Intercept only" = "intercept_only",
                                  "Intercept + Slope" = "slope",
                                  "None" = "none"),
                       selected = "intercept_only"),
            
            selectInput("tg_cage_effect",
                       "Cage Effect Handling",
                       choices = c("Include if not collinear" = "include_if_not_collinear",
                                  "Always include" = "always_include",
                                  "Never include" = "never_include"),
                       selected = "include_if_not_collinear"),

            tags$details(
              style = "margin-top: 10px;",
              tags$summary(
                style = "cursor: pointer; color: #555; font-size: 0.9em;",
                "Advanced Options"
              ),
              div(style = "padding-top: 8px;",
                sliderInput("tg_poly_degree",
                           "Polynomial Degree",
                           min = 1, max = 3, value = 1),
                numericInput("tg_extrapolation_points",
                            "Extrapolation Points",
                            value = 3, min = 0, max = 10, step = 1),
                helpText("Number of final time-points used to extrapolate tumor volumes ",
                         "for mice that died before the last study day. ",
                         "Set to 0 to disable extrapolation. Requires an Event/Censor column ",
                         "to be mapped in Data Upload (for the growth curves plot).")
              )
            )
          ),
          
          selectInput("tg_reference",
                     "Reference Group (Required)",
                     choices = c("-- Select a reference group --" = "")),
          
          selectInput("tg_comparison_type",
                     "Post-hoc Comparisons",
                     choices = c("vs Reference Group" = "reference",
                                "All Pairwise" = "pairwise",
                                "Tukey HSD" = "tukey"),
                     selected = "reference"),
          
          br(),
          
          actionButton("run_tg_analysis", 
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),
          
          uiOutput("tg_status_msg")
        )
      ),
      
      column(9,
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
            tab_info(strong("Per-animal growth rates "), "calculated from individual tumor volume trajectories. Used to assess within-group variability and for power analysis effect size estimation. Not the primary model output \u2014 use Treatment Effects for group-level inference."),
            shinycssloaders::withSpinner(
              DTOutput("tg_growth_rates_table")
            ),
            downloadButton("download_tg_rates", "Download Growth Rates")
          ),
          
          tabPanel("Plots",
            tab_info(
              strong("Tumor Growth Curves: "), "individual mouse trajectories with group mean overlay. ",
              strong("Growth Rate Distribution: "), "spread of per-animal rates within each group. ",
              strong("AUC / Treatment Effects: "), "model-estimated effect sizes per group. ",
              strong("Q\u2011Q Plot: "), "assesses normality of model residuals \u2014 points should follow the diagonal line for valid inference."
            ),
            fluidRow(
              column(6,
                h4("Tumor Growth Curves"),
                shinycssloaders::withSpinner(
                  plotOutput("tg_plot_growth", height = "400px")
                ),
                downloadButton("download_tg_plot_growth", "Download Plot")
              ),
              column(6,
                h4("Growth Rate Distribution"),
                shinycssloaders::withSpinner(
                  plotOutput("tg_plot_rates", height = "400px")
                ),
                downloadButton("download_tg_plot_rates", "Download Plot")
              )
            ),
            
            fluidRow(
              column(6,
                h4("AUC Comparison"),
                shinycssloaders::withSpinner(
                  plotOutput("tg_plot_auc", height = "400px")
                ),
                downloadButton("download_tg_plot_auc", "Download Plot")
              ),
              column(6,
                h4("Random Effects"),
                shinycssloaders::withSpinner(
                  plotOutput("tg_plot_caterpillar", height = "400px")
                ),
                downloadButton("download_tg_plot_caterpillar", "Download Plot")
              )
            ),
            
            hr(),
            
            fluidRow(
              column(12,
                h4("Treatment Schedule"),
                p(tags$small("Requires a treatment schedule file uploaded in Data Upload (Step 1b)")),
                shinycssloaders::withSpinner(
                  plotOutput("tg_plot_treatment_schedule", height = "350px")
                ),
                downloadButton("download_tg_plot_schedule", "Download Plot")
              )
            )
          ),
          
          tabPanel("Diagnostics",
            tab_info(strong("Model diagnostics: "), "convergence status, singularity warnings (random effects structure may be too complex), residual statistics, and variance component estimates. A singular fit warning suggests simplifying the random effects (e.g. removing the slope term)."),
            verbatimTextOutput("tg_diagnostics")
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
    
    h3("Column Naming"),
    p("The dashboard looks for common column names like:"),
    tags$ul(
      tags$li(strong("ID:"), "ID, Mouse_ID, Animal_ID, Subject"),
      tags$li(strong("Time:"), "Day, Time, Timepoint"),
      tags$li(strong("Volume:"), "Volume, Tumor_Volume"),
      tags$li(strong("Treatment:"), "Treatment, Group"),
      tags$li(strong("Cage:"), "Cage, Cage_ID")
    ),
    
    h3("Statistical Methods"),
    p("The analysis uses linear mixed-effects models (lme4) or area under the curve (AUC) methods."),
    p("Consult the", a("mouseExperiment documentation", href = "https://github.com/"), 
      "for detailed method information."),
    
    h3("Troubleshooting"),
    p(strong("File won't upload:"), "Ensure it's a CSV or RDA file."),
    p(strong("Columns not detected:"), "Check that column names match expected patterns."),
    p(strong("Analysis fails:"), "Check that required columns are selected and contain valid data."),
    
    h3("Contact & Support"),
    p("For issues or questions, contact the package maintainers.")
  )
}
