#' Survival Analysis UI
#' @keywords internal
survivalUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Survival Analysis"),
    
    p("Analyze time-to-event data using Kaplan-Meier survival curves and proportional hazards models. ",
      "Ensure data is uploaded and configured in the Data Upload tab."),
    
    tags$div(class = "row collapsible-sidebar-row",
      tags$div(class = "col-sm-3 sidebar-panel-col",
        wellPanel(
          h4("Analysis Parameters"),
          
          uiOutput(ns("surv_reference_group_ui")),
          
          br(),
          
          actionButton(ns("run_surv_analysis"), 
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),
          
          uiOutput(ns("surv_status_msg"))
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
          id = ns("surv_results_tabs"),
          
          tabPanel("Summary",
            tab_info(strong("Survival summary: "), "Per-group event counts, total subjects, event rate, and median survival time (the day by which 50\u202f% of animals have reached the endpoint). \u2018Not reached\u2019 means fewer than half the group experienced an event by the end of the study."),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("surv_summary"))
            )
          ),
          
          tabPanel("Cox Model",
            tab_info(
              strong("Cox proportional hazards model. "),
              "Hazard Ratios (HR) are computed relative to the reference group. Firth\u2019s bias-reduced method is applied automatically when a group has zero events.",
              tags$ul(
                tags$li(strong("Group"), " \u2014 treatment group name."),
                tags$li(strong("HR"), " \u2014 Hazard Ratio: the instantaneous risk of the event in this group relative to the reference. HR\u202f<\u202f1 means lower hazard (better survival); HR\u202f>\u202f1 means higher hazard."),
                tags$li(strong("CI_Lower / CI_Upper"), " \u2014 95\u202f% confidence interval for the HR. If this interval does not include 1.0, the difference is statistically significant at the 5\u202f% level."),
                tags$li(strong("P_Value"), " \u2014 p-value for the individual group coefficient from the Cox model."),
                tags$li(strong("Events"), " \u2014 number of animals that experienced the endpoint (e.g. reached the tumour size threshold or were euthanised)."),
                tags$li(strong("Total"), " \u2014 total number of animals in the group."),
                tags$li(strong("Event_Rate"), " \u2014 Events \u00f7 Total. Values close to 1.0 mean most animals reached the endpoint."),
                tags$li(strong("Median_Survival"), " \u2014 day by which 50\u202f% of the group experienced the event. \u2018Not reached\u2019 if fewer than half had events by study end."),
                tags$li(strong("Note"), " \u2014 flags the reference group row.")
              )
            ),
            shinycssloaders::withSpinner(
              verbatimTextOutput(ns("surv_cox_summary"))
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("surv_cox_table"))
            ),
            downloadButton(ns("download_surv_cox"), "Download Results")
          ),
          
          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  selectInput(ns("surv_plot_selector"), "Select Plot",
                    choices = c(
                      "Kaplan-Meier Curve"       = "km",
                      "Hazard Ratio Forest Plot" = "forest"
                    )
                  ),
                  ctrl_section("surv-groups", "sort", "Groups",
                    conditionalPanel(ns = ns,"input.surv_plot_selector == 'km'",
                      uiOutput(ns("surv_group_order_color_ui"))
                    ),
                    conditionalPanel(ns = ns,"input.surv_plot_selector == 'forest'",
                      tags$p(class = "text-muted", style = "font-size:12px; margin:4px 0;",
                             "Group colours are mapped automatically from the analysis groups.")
                    )
                  ),
                  ctrl_section("surv-elements", "layer-group", "Elements",
                    conditionalPanel(ns = ns,"input.surv_plot_selector == 'km'",
                      checkboxInput(ns("surv_opt_conf_int"),
                                    tagList("Show confidence interval",
                                            help_tip("Displays 95% confidence bands around the Kaplan-Meier survival curves.")),
                                    value = FALSE),
                      checkboxInput(ns("surv_opt_risk_table"),
                                    tagList("Show risk table on plot",
                                            help_tip("Number of subjects still at risk at each timepoint, displayed beneath the KM curve.")),
                                    value = TRUE),
                      numericInput(ns("surv_km_xstart"),
                                   "X-axis start (Day)",
                                   value = 0, min = 0, step = 1)
                    ),
                    conditionalPanel(ns = ns,"input.surv_plot_selector == 'forest'",
                      colourpicker::colourInput(ns("surv_forest_pt_color"),
                                               "Point colour", value = "steelblue"),
                      h6("Annotation"),
                      checkboxInput(ns("surv_forest_show_hr"), "Show HR",      value = TRUE),
                      checkboxInput(ns("surv_forest_show_ci"), "Show 95% CI",  value = TRUE),
                      checkboxInput(ns("surv_forest_show_p"),  "Show p-value", value = TRUE)
                    ),
                    selectInput(ns("surv_opt_legend_pos"), "Legend position",
                      choices = c("None"   = "none",
                                  "Top"    = "top",    "Bottom" = "bottom",
                                  "Left"   = "left",   "Right"  = "right",
                                  "In plot \u2013 top left"     = "inplot_tl",
                                  "In plot \u2013 top right"    = "inplot_tr",
                                  "In plot \u2013 bottom left"  = "inplot_bl",
                                  "In plot \u2013 bottom right" = "inplot_br"),
                      selected = "right"
                    ),
                    textInput(ns("surv_plot_title"), "Title",
                              placeholder = "Leave blank for default"),
                    selectInput(ns("surv_plot_title_align"), "Title alignment",
                      choices = c("Left" = "left", "Center" = "center", "Right" = "right"),
                      selected = "center"),
                    textInput(ns("surv_plot_subtitle"), "Subtitle",
                              placeholder = "Optional subtitle"),
                    textInput(ns("surv_plot_xaxis_label"), "X-axis label",
                              placeholder = "Leave blank for default"),
                    textInput(ns("surv_plot_yaxis_label"), "Y-axis label",
                              placeholder = "Leave blank for default"),
                    textInput(ns("surv_plot_legend_label"), "Legend label",
                              placeholder = "Leave blank for default")
                  ),
                  ctrl_section("surv-scaling", "expand", "Scaling",
                    sliderInput(ns("surv_opt_text_size"), "Text size (pt)",
                                min = 10, max = 24, value = 16, step = 1)
                  ),
                  ctrl_section("surv-export", "arrows-alt", "Plot Size / Export",
                    fluidRow(
                      column(6, numericInput(ns("surv_plot_width"),  "Width (px)",  800, min = 200, max = 5000, step = 50)),
                      column(6, numericInput(ns("surv_plot_height"), "Height (px)", 600, min = 200, max = 5000, step = 50))
                    ),
                    numericInput(ns("surv_plot_dpi"), "Resolution (DPI)", 100, min = 72, max = 600, step = 10),
                    tags$div(style = "margin-top: 10px;", uiOutput(ns("surv_download_plot_btn")))
                  )
                )
              ),
              column(9,
                conditionalPanel(ns = ns,"input.surv_plot_selector == 'km'",
                  h4("Kaplan-Meier Survival Curves"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("surv_plot_km"), width = "auto", height = "auto")
                  )
                ),
                conditionalPanel(ns = ns,"input.surv_plot_selector == 'forest'",
                  h4("Hazard Ratio Forest Plot"),
                  shinycssloaders::withSpinner(
                    plotOutput(ns("surv_forest_plot"), width = "auto", height = "auto")
                  )
                )
              )
            )
          ),

          tabPanel("Risk Table",
            tab_info(
              strong("Numbers at risk table: "), "one row per group per event time point, summarising what happened to the cohort at that moment.",
              tags$ul(
                tags$li(strong("Group"), " \u2014 treatment group name."),
                tags$li(strong("Time"), " \u2014 the day on which one or more events occurred."),
                tags$li(strong("N_Risk"), " \u2014 number of animals still alive and uncensored at the start of this time interval."),
                tags$li(strong("N_Events"), " \u2014 number of animals that experienced the endpoint at exactly this time."),
                tags$li(strong("N_Censored"), " \u2014 number of animals whose last observation was at this time without experiencing the endpoint (e.g. removed from study)."),
                tags$li(strong("Survival"), " \u2014 Kaplan\u2013Meier estimated survival probability at this time point (1.0\u202f=\u202f100\u202f% surviving).")
              ),
              "N_Risk decreases over time as events and censoring accumulate; a rapid drop indicates many events in a short window."
            ),
            shinycssloaders::withSpinner(
              DTOutput(ns("surv_risk_table"))
            ),
            downloadButton(ns("download_surv_risk_table"), "Download Table")
          )
        )
      )
    )
  )
}

#' Survival Analysis Server Logic
#' @keywords internal
survival_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
  
  rv$surv_results <- NULL
  rv$surv_running <- FALSE
  
  # Dynamically generate reference group dropdown based on treatment groups
  output$surv_reference_group_ui <- renderUI({
    req(rv$uploaded_data, rv$data_valid)

    # Determine which column / labels to use.
    # When a dose column is mapped rv$tx_col == "TreatmentGroup" and the
    # actual groups are "Treatment Dose" combinations.  Build that vector
    # on-the-fly from the raw upload so the dropdown is always correct even
    # before a full analysis has been run.
    get_groups <- function() {
      d <- rv$uploaded_data
      # Resolve dose column directly from the UI input, not rv$dose_col which
      # requires TG analysis to have run first.
      dose_col_input <- if (!is.null(rv$settings$col_dose) && nzchar(rv$settings$col_dose) &&
                            rv$settings$col_dose != "None") rv$settings$col_dose else NULL
      has_dose <- has_col(d, dose_col_input)
      if (has_dose &&
          !is.null(rv$settings$col_treatment) && nzchar(rv$settings$col_treatment) &&
          has_col(d, rv$settings$col_treatment)) {
        combined <- paste(d[[rv$settings$col_treatment]], d[[dose_col_input]], sep = " ")
        return(sort(unique(as.character(na.omit(combined)))))
      }
      if (!is.null(rv$settings$col_treatment) && nzchar(rv$settings$col_treatment) &&
          has_col(d, rv$settings$col_treatment)) {
        return(sort(unique(as.character(na.omit(d[[rv$settings$col_treatment]])))))
      }
      character(0)
    }

    treatment_groups <- get_groups()

    if (length(treatment_groups) > 0) {
      default_ref <- pick_reference_group(treatment_groups)
      return(
        selectInput(session$ns("surv_reference_group"),
                   "Reference Group",
                   choices = treatment_groups,
                   selected = default_ref)
      )
    }

    # Fallback if no data
    selectInput(session$ns("surv_reference_group"),
               "Reference Group",
               choices = "",
               selected = "")
  })
  
  # Run survival analysis
  observeEvent(input$run_surv_analysis, {
    req(rv$uploaded_data, rv$data_valid)
    
    # Validate that required columns are configured from data upload
    if (rv$settings$col_treatment == "" || rv$settings$col_treatment == "-- Select a reference group --") {
      output$surv_status_msg <- renderUI({
        status_alert("Error: Please select a Treatment Group on the Data Upload tab", "danger")
      })
      return()
    }
    

    if (is.null(rv$settings$col_survival_event) || rv$settings$col_survival_event == "" || rv$settings$col_survival_event == "None") {
      output$surv_status_msg <- renderUI({
        status_alert("Error: Please select an Event/Censor Status Column on the Data Upload tab", "danger")
      })
      return()
    }
    
    rv$surv_running <- TRUE
    output$surv_status_msg <- renderUI({
      status_alert("Running survival analysis...", "info")
    })
    
    # Prepare data — reuse canonical preprocessing if available
    pp <- ensure_preprocessed(rv)
    data <- pp$data

    # Resolve dose column from UI input directly — do NOT use rv$dose_col which
    # is only set after TG analysis runs. rv$settings$col_dose always reflects the
    # current column mapping regardless of analysis order.
    dose_col_surv <- if (!is.null(rv$settings$col_dose) && nzchar(rv$settings$col_dose) &&
                         rv$settings$col_dose != "None") rv$settings$col_dose else NULL

    # When a dose column is mapped, build a combined "Treatment Dose" group label
    # so each dose level is treated as its own group.
    has_dose_surv <- has_col(data, dose_col_surv) && has_col(data, rv$settings$col_treatment)
    if (has_dose_surv) {
      data[["TreatmentGroup"]] <- paste(data[[rv$settings$col_treatment]],
                                        data[[dose_col_surv]], sep = " ")
    }

    # Handle time column — reuse from preprocessed data
    time_col <- pp$time_col
    
    # For survival analysis, reduce data to one row per animal
    # Use the user-selected ID and Cage columns (not hardcoded names)
    id_col_surv   <- if (!is.null(rv$settings$col_id)   && nzchar(rv$settings$col_id))   rv$settings$col_id   else "ID"
    cage_col_surv <- if (!is.null(rv$settings$col_cage) && nzchar(rv$settings$col_cage) &&
                         rv$settings$col_cage != "None") rv$settings$col_cage else NULL

    has_id   <- has_col(data, id_col_surv)
    has_cage <- has_col(data, cage_col_surv)

    if (has_id) {
      # Build key columns for deduplication
      key_cols <- if (has_cage) c(id_col_surv, cage_col_surv) else id_col_surv
      unique_animals <- unique(data[, key_cols, drop = FALSE])
      survival_rows <- vector("list", nrow(unique_animals))

      for (i in seq_len(nrow(unique_animals))) {
        if (has_cage) {
          mask <- data[[id_col_surv]]   == unique_animals[i, id_col_surv] &
                  data[[cage_col_surv]] == unique_animals[i, cage_col_surv]
        } else {
          mask <- data[[id_col_surv]] == unique_animals[i, id_col_surv]
        }
        animal_rows <- data[mask, , drop = FALSE]
        if (nrow(animal_rows) > 0) {
          survival_rows[[i]] <- animal_rows[nrow(animal_rows), , drop = FALSE]
        }
      }
      survival_data <- do.call(rbind, Filter(Negate(is.null), survival_rows))

      # Validate aggregation
      if (is.null(survival_data) || nrow(survival_data) == 0) {
        result <- list(error = TRUE, message = "Data aggregation resulted in empty dataset. Check ID and Cage columns.")
        rv$surv_results <- result
        rv$surv_running <- FALSE
        output$surv_status_msg <- renderUI({
          status_alert(paste("Error:", result$message), "danger")
        })
        return()
      }

      data <- survival_data
    } else {
      # No ID column mapped — operate on data as-is (assume one row per animal)
      if (nrow(data) == 0) {
        result <- list(error = TRUE, message = "No data available for analysis.")
        rv$surv_results <- result
        rv$surv_running <- FALSE
        output$surv_status_msg <- renderUI({
          status_alert("Error: No data available for analysis", "danger")
        })
        return()
      }
    }
    
    # Run analysis with processed data and calculated/selected time column
    tx_col_surv <- if (has_dose_surv) "TreatmentGroup" else rv$settings$col_treatment

    # Validate all required columns exist before analysis
    required <- c(Treatment = tx_col_surv,
                  Time      = time_col,
                  Event     = rv$settings$col_survival_event)
    missing_cols <- check_required_cols(data, required)
    if (length(missing_cols) > 0) {
      output$surv_status_msg <- renderUI({
        status_alert(paste0("Required column(s) not found: ",
                   paste(names(missing_cols), "=", sQuote(missing_cols), collapse = ", ")), "danger")
      })
      rv$surv_running <- FALSE
      return()
    }

    progress <- shiny::Progress$new(session, min = 0, max = 1)
    progress$set(message = "Running survival analysis...", value = 0.2)
    on.exit(progress$close(), add = TRUE)

    result <- safe_analysis(module = "survival", {
      mouseExperiment::survival_statistics(
        df = data,
        time_column = time_col,
        censor_column = rv$settings$col_survival_event,
        treatment_column = tx_col_surv,
        cage_column = if (!is.null(rv$settings$col_cage) && nzchar(rv$settings$col_cage)) rv$settings$col_cage else "Cage",
        dose_column = dose_col_surv,
        reference_group = input$surv_reference_group,
        firth_correction = TRUE
      )
    })
    progress$set(value = 0.9, detail = "Preparing results")
    
    # Store time calculation info for later display
    result$time_method <- rv$settings$time_method
    result$time_column_used <- time_col
    result$reference_date <- rv$settings$reference_date
    
    # Generate Kaplan-Meier plot if analysis succeeded
    if (!isTRUE(result$error) && !is.null(result$survival_data)) {
      tryCatch({
        # Use the survival_data returned by survival_statistics (has Time, Event, Treatment columns)
        surv_df <- result$survival_data

        # Keep a survfit object for the risk table tab
        km_fit <- survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df)
        result$km_fit <- km_fit

        # Use the exported plot_survival() function from the package to generate the KM plot
        # (matches demo notebook approach; returns a ggsurvplot object)
        result$km_plot <- mouseExperiment::plot_survival(
          data              = data,
          time_column       = time_col,
          censor_column     = rv$settings$col_survival_event,
          treatment_column  = tx_col_surv,
          id_column         = rv$settings$col_id,
          cage_column       = if (!is.null(rv$settings$col_cage) && nzchar(rv$settings$col_cage)) rv$settings$col_cage else NULL,
          dose_column       = dose_col_surv,
          colors            = NULL,
          show_risk_table   = TRUE,
          title             = "Kaplan-Meier Survival Curves by Treatment"
        )

      }, error = function(e) {
        # Fallback to manual ggsurvplot if plot_survival() fails
        log_caught_error("survival", "plot_survival()", e)
        tryCatch({
          surv_df <- result$survival_data
          km_fit <- survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df)
          result$km_fit <<- km_fit
          km_plot_obj <- survminer::ggsurvplot(
            km_fit,
            data = surv_df,
            risk.table = TRUE,
            risk.table.y.text.col = TRUE,
            conf.int = FALSE,
            xlab = "Time (Days)",
            ylab = "Survival Probability",
            title = "Kaplan-Meier Survival Curves by Treatment"
          )
          result$km_plot <<- km_plot_obj
        }, error = function(e2) {
          log_caught_error("survival", "ggsurvplot fallback", e2)
          NULL
        })
      })
    }
    
    rv$surv_results <- result
    rv$surv_running <- FALSE

    # Store treatment group names for the per-group colour pickers in the Plots tab
    rv$surv_plot_groups <- if (!is.null(result$survival_data) && has_col(result$survival_data, "Treatment")) {
      unique(as.character(result$survival_data$Treatment))
    } else NULL
    
    if (isTRUE(result$error)) {
      output$surv_status_msg <- renderUI({
        status_alert(paste("Error:", result$message), "danger")
      })
    } else {
      output$surv_status_msg <- renderUI({
        status_alert("Survival analysis completed successfully!", "success")
      })
    }
  })
  
  # Summary output
  output$surv_summary <- renderText({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) {
      return(paste("Error:", rv$surv_results$message))
    }
    
    summary_text <- capture.output({
      cat("=== SURVIVAL ANALYSIS SUMMARY ===\n\n")
      
      # Show time specification method
      if (rv$surv_results$time_method == "from_date") {
        cat("Time Method: Calculated from Date Column\n")
        cat("Date Column:", rv$settings$col_date_for_day, "\n")
        cat("Reference Date:", rv$settings$reference_date, "\n")
      } else {
        cat("Time Method: Direct Day Column\n")
        cat("Day Column:", rv$settings$col_day, "\n")
      }
      
      cat("Event Column:", rv$settings$col_survival_event, "\n")
      cat("Treatment Column:", rv$settings$col_treatment, "\n")
      cat("Statistical Method:", rv$surv_results$method_used, "\n")
      
      if (!is.null(rv$surv_results$results)) {
        res <- rv$surv_results$results
        cat("Reference Group:", rv$surv_results$reference_group, "\n\n")
        
        cat("--- Median Survival ---\n")
        for (i in seq_len(nrow(res))) {
          med <- if (is.na(res$Median_Survival[i])) "Not reached" else paste0(res$Median_Survival[i], " days")
          note <- if (!is.na(res$Note[i]) && res$Note[i] != "") paste0(" [", res$Note[i], "]") else ""
          cat(sprintf("  %-20s %s%s\n", res$Group[i], med, note))
        }
      }
    })
    
    paste(summary_text, collapse = "\n")
  })
  
  # Group order and color combined control for KM/Survival plots
  output$surv_group_order_color_ui <- renderUI({
    req(rv$surv_plot_groups)
    group_order_and_color_ui("surv", rv$surv_plot_groups, ns = session$ns)
  })

  # Apply colorblind-friendly palette presets when selector changes
  observe_palette_preset(input, session, "surv", reactive(rv$surv_plot_groups))

  output$surv_download_plot_btn <- renderUI({
    btn_id <- switch(input$surv_plot_selector %||% "km",
      km     = "download_surv_plot_km",
      forest = "download_surv_forest_plot",
      "download_surv_plot_km"
    )
    downloadButton(session$ns(btn_id), "Download Plot")
  })

  # Kaplan-Meier plot — re-renders live when colour pickers, order, or conf.int option change
  output$surv_plot_km <- renderPlot({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) return(NULL)
    req(rv$surv_results$km_fit, rv$surv_results$survival_data)

    surv_df <- rv$surv_results$survival_data

    # Apply user group order as factor levels and refit KM so strata follow that order
    ordered_grps <- get_ordered_groups(input, "surv", rv$surv_plot_groups)
    if (!is.null(ordered_grps) && has_col(surv_df, "Treatment")) {
      surv_df$Treatment <- factor(surv_df$Treatment, levels = ordered_grps)
    }
    km_fit <- tryCatch(
      survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df),
      error = function(e) {
        log_caught_error("survival", "KM refit for colour/order", e)
        rv$surv_results$km_fit
      }
    )

    # Build colour vector from per-group colourInputs (ordered to match strata)
    groups  <- if (!is.null(ordered_grps)) ordered_grps else rv$surv_plot_groups
    pal     <- build_color_vector(input, "surv_color", groups)
    # unname so ggsurvplot assigns colours positionally (strata are "Treatment=X", not "X")
    if (!is.null(pal)) pal <- unname(pal)
    if (is.null(pal)) pal <- "hue"

    # X-axis start
    x_start <- if (!is.null(input$surv_km_xstart) && !is.na(input$surv_km_xstart))
                  as.numeric(input$surv_km_xstart) else 0
    text_size   <- input$surv_opt_text_size %||% 16
    plot_title  <- if (nzchar(trimws(input$surv_plot_title %||% "")))
                     input$surv_plot_title else "Kaplan-Meier Survival Curves by Treatment"
    plot_sub    <- input$surv_plot_subtitle %||% ""
    title_hjust <- switch(input$surv_plot_title_align %||% "center",
                          left = 0, center = 0.5, right = 1, 0.5)
    xlab_custom <- trimws(input$surv_plot_xaxis_label  %||% "")
    ylab_custom <- trimws(input$surv_plot_yaxis_label  %||% "")
    leg_label   <- trimws(input$surv_plot_legend_label %||% "")
    legend_pos  <- input$surv_opt_legend_pos %||% "right"

    tryCatch({
      km_obj <- survminer::ggsurvplot(
        km_fit,
        data                  = surv_df,
        palette               = pal,
        risk.table            = isTRUE(input$surv_opt_risk_table),
        risk.table.y.text.col = TRUE,
        conf.int              = isTRUE(input$surv_opt_conf_int),
        xlim                  = c(x_start, NA),
        xlab                  = if (nzchar(xlab_custom)) xlab_custom else "Time (Days)",
        ylab                  = if (nzchar(ylab_custom)) ylab_custom else "Survival Probability",
        legend.title          = if (nzchar(leg_label)) leg_label else "Treatment",
        title                 = plot_title,
        ggtheme               = ggplot2::theme_classic(base_size = text_size)
      )
      # Apply title alignment and optional subtitle
      km_obj$plot <- km_obj$plot +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = title_hjust))
      if (nzchar(trimws(plot_sub))) {
        km_obj$plot <- km_obj$plot +
          ggplot2::labs(subtitle = plot_sub) +
          ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = title_hjust))
      }
      # Apply legend position (ggsurvplot returns a list; position touches the plot component)
      lp <- if (startsWith(legend_pos, "inplot_")) {
        coord <- switch(legend_pos,
          inplot_tl = c(0.02, 0.98), inplot_tr = c(0.98, 0.98),
          inplot_bl = c(0.02, 0.02), inplot_br = c(0.98, 0.02), c(0.5, 0.5))
        just <- switch(legend_pos,
          inplot_tl = c("left","top"),  inplot_tr = c("right","top"),
          inplot_bl = c("left","bottom"), inplot_br = c("right","bottom"),
          c("center","center"))
        km_obj$plot <- km_obj$plot +
          ggplot2::theme(
            legend.position    = coord,
            legend.justification = just,
            legend.background  = ggplot2::element_rect(
              fill = scales::alpha("white", 0.8), colour = NA))
      } else {
        km_obj$plot <- km_obj$plot + ggplot2::theme(legend.position = legend_pos)
      }
      print(km_obj)
    }, error = function(e) {
      # Fallback: print pre-generated stored KM plot
      if (!is.null(rv$surv_results$km_plot)) {
        print(rv$surv_results$km_plot)
      } else {
        plot(1, 1, type = "n", main = "KM plot not available", xlab = "", ylab = "")
      }
    })  
  }, width = function() input$surv_plot_width %||% 800,
     height = function() input$surv_plot_height %||% 600)
  
  # Cox model summary
  output$surv_cox_summary <- renderText({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$model)) {
      return("Model not available")
    }
    
    paste(capture.output(summary(rv$surv_results$model)), collapse = "\n")
  })
  
  # Results table
  output$surv_cox_table <- renderDT({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$results)) {
      return(datatable(data.frame(Message = "Results not available")))
    }
    
    results_df <- rv$surv_results$results
    results_df <- format_df_for_display(results_df)
    
    datatable(results_df, options = list(scrollX = TRUE, paging = FALSE, searching = FALSE))
  })
  
  # Risk table output - numbers at risk per time point per group
  output$surv_risk_table <- renderDT({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$km_fit)) {
      return(datatable(data.frame(Message = "Risk table not available")))
    }
    
    km_sum <- summary(rv$surv_results$km_fit)
    risk_df <- data.frame(
      Time      = km_sum$time,
      Group     = sub("^Treatment=", "", as.character(km_sum$strata)),
      N_Risk    = km_sum$n.risk,
      N_Events  = km_sum$n.event,
      N_Censored = km_sum$n.censor,
      Survival  = round(km_sum$surv, 4)
    )
    # Reorder columns for readability
    risk_df <- risk_df[, c("Group", "Time", "N_Risk", "N_Events", "N_Censored", "Survival")]
    risk_df <- risk_df[order(risk_df$Group, risk_df$Time), ]
    
    datatable(risk_df,
              options = list(scrollX = TRUE, paging = FALSE, searching = FALSE),
              rownames = FALSE)
  })
  
  # Helper: build the HR forest plot (shared by renderPlot + download handler)
  build_surv_forest_plot <- function(res, ref_grp, pt_color, show_hr, show_ci, show_p) {
    plot_df <- res[!is.na(res$HR), ]
    if (!is.null(ref_grp)) plot_df <- plot_df[plot_df$Group != ref_grp, ]
    if (nrow(plot_df) == 0) return(NULL)

    # Add reference group row at HR = 1 (no CI)
    ref_row <- plot_df[1L, , drop = FALSE]
    ref_row[1L, ] <- NA
    ref_row$Group    <- paste0(ref_grp, "\n(Reference)")
    ref_row$HR       <- 1
    ref_row$CI_Lower <- NA_real_
    ref_row$CI_Upper <- NA_real_
    plot_df <- rbind(plot_df, ref_row)

    # Ensure numeric columns stay numeric after rbind
    for (nc in c("HR", "CI_Lower", "CI_Upper", "P_Value")) {
      if (nc %in% colnames(plot_df)) plot_df[[nc]] <- as.numeric(plot_df[[nc]])
    }
    plot_df$is_ref <- plot_df$Group == paste0(ref_grp, "\n(Reference)")

    # Symmetric log10 axis (exclude reference row)
    all_hr_vals <- c(plot_df$HR[!plot_df$is_ref],
                     plot_df$CI_Lower[!plot_df$is_ref],
                     plot_df$CI_Upper[!plot_df$is_ref])
    all_hr_vals <- all_hr_vals[!is.na(all_hr_vals) & all_hr_vals > 0]
    max_log10 <- max(ceiling(max(abs(log10(all_hr_vals)), na.rm = TRUE)), 1L)
    x_lim    <- c(10^(-max_log10), 10^max_log10)
    x_breaks <- 10^seq(-max_log10, max_log10, by = 1)
    x_labels <- ifelse(x_breaks >= 1,
                       as.character(round(x_breaks, 0)),
                       as.character(x_breaks))

    # Build per-row annotation labels (vapply avoids the apply(sapply) dim error)
    hr_vals  <- as.numeric(plot_df$HR)
    ci_lo    <- as.numeric(plot_df$CI_Lower)
    ci_hi    <- as.numeric(plot_df$CI_Upper)
    p_vals   <- as.numeric(plot_df$P_Value)
    plot_df$label_text <- vapply(seq_len(nrow(plot_df)), function(i) {
      if (isTRUE(plot_df$is_ref[i])) return("")
      parts <- character(0)
      if (show_hr && !is.na(hr_vals[i]))  parts <- c(parts, paste0("HR: ", round(hr_vals[i], 2)))
      if (show_ci && !is.na(ci_lo[i]))    parts <- c(parts, paste0("[", round(ci_lo[i], 2), "\u2013", round(ci_hi[i], 2), "]"))
      if (show_p  && !is.na(p_vals[i]))   parts <- c(parts, paste0("p = ", format_pval(p_vals[i])))
      paste(parts, collapse = "  ")
    }, character(1))

    ggplot2::ggplot(plot_df,
                   ggplot2::aes(x = HR, y = reorder(Group, HR),
                                xmin = CI_Lower, xmax = CI_Upper)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
      ggplot2::geom_errorbarh(
        data = ~ subset(., !is_ref),
        height = 0.25, linewidth = 0.8, color = pt_color
      ) +
      ggplot2::geom_point(
        data = ~ subset(., !is_ref), size = 5, color = pt_color
      ) +
      ggplot2::geom_point(
        data = ~ subset(., is_ref), size = 5, color = "grey50", shape = 18
      ) +
      { if (any(nzchar(plot_df$label_text[!plot_df$is_ref])))
          ggplot2::geom_text(
            data = ~ subset(., !is_ref),
            ggplot2::aes(label = label_text),
            vjust = 3, hjust = -0.1, size = 4, color = "grey30")
        else list() } +
      ggplot2::scale_x_continuous(
        trans  = "log10",
        limits = x_lim,
        breaks = x_breaks,
        labels = x_labels
      ) +
      ggplot2::labs(
        title = paste("Forest Plot: Hazard Ratios vs", ref_grp),
        x     = "Hazard Ratio (log\u2081\u2080 scale)",
        y     = NULL
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        axis.text.y    = ggplot2::element_text(size = 11),
        legend.position = "none"
      )
  }

  # Forest plot
  output$surv_forest_plot <- renderPlot({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$results)) return(NULL)
    
    res <- rv$surv_results$results
    req(all(c("Group", "HR", "CI_Lower", "CI_Upper", "P_Value") %in% colnames(res)))
    
    ref_grp <- rv$surv_results$reference_group
    
    # Verify at least one non-reference group with valid HR
    plot_df <- res[!is.na(res$HR), ]
    if (!is.null(ref_grp)) plot_df <- plot_df[plot_df$Group != ref_grp, ]
    req(nrow(plot_df) > 0)

    pt_color <- input$surv_forest_pt_color %||% "steelblue"
    show_hr  <- isTRUE(input$surv_forest_show_hr)
    show_ci  <- isTRUE(input$surv_forest_show_ci)
    show_p   <- isTRUE(input$surv_forest_show_p)

    p <- build_surv_forest_plot(res, ref_grp, pt_color, show_hr, show_ci, show_p)
    req(!is.null(p))
    text_size   <- input$surv_opt_text_size %||% 16
    legend_pos  <- input$surv_opt_legend_pos %||% "right"
    p <- apply_plot_labels(p,
                           input$surv_plot_title        %||% "",
                           input$surv_plot_subtitle     %||% "",
                           input$surv_plot_title_align  %||% "center",
                           text_size,
                           xlab       = input$surv_plot_xaxis_label  %||% "",
                           ylab       = input$surv_plot_yaxis_label  %||% "",
                           legend_lab = input$surv_plot_legend_label %||% "")
    apply_legend_pos(p, legend_pos)
  }, width = function() input$surv_plot_width %||% 800,
     height = function() input$surv_plot_height %||% 600)
  
  # Download handlers
  output$download_surv_plot_km <- downloadHandler(
    filename = "kaplan_meier_plot.png",
    content = function(file) {
      req(rv$surv_results)
      surv_df <- rv$surv_results$survival_data
      ordered_grps <- get_ordered_groups(input, "surv", rv$surv_plot_groups)
      if (!is.null(ordered_grps) && has_col(surv_df, "Treatment"))
        surv_df$Treatment <- factor(surv_df$Treatment, levels = ordered_grps)
      km_fit <- tryCatch(
        survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df),
        error = function(e) {
          log_caught_error("survival", "KM refit for download", e)
          rv$surv_results$km_fit
        }
      )
      groups <- if (!is.null(ordered_grps)) ordered_grps else rv$surv_plot_groups
      pal <- build_color_vector(input, "surv_color", groups)
      if (!is.null(pal)) pal <- unname(pal)
      if (is.null(pal)) pal <- "hue"
      x_start     <- if (!is.null(input$surv_km_xstart) && !is.na(input$surv_km_xstart))
                       as.numeric(input$surv_km_xstart) else 0
      xlab_custom <- trimws(input$surv_plot_xaxis_label  %||% "")
      ylab_custom <- trimws(input$surv_plot_yaxis_label  %||% "")
      leg_label   <- trimws(input$surv_plot_legend_label %||% "")
      legend_pos  <- input$surv_opt_legend_pos %||% "right"
      plot_title  <- if (nzchar(trimws(input$surv_plot_title %||% "")))
                       input$surv_plot_title else "Kaplan-Meier Survival Curves by Treatment"
      plot_sub    <- input$surv_plot_subtitle %||% ""
      title_hjust <- switch(input$surv_plot_title_align %||% "center",
                            left = 0, center = 0.5, right = 1, 0.5)
      km_obj <- tryCatch(
        survminer::ggsurvplot(
          km_fit,
          data                  = surv_df,
          palette               = pal,
          risk.table            = isTRUE(input$surv_opt_risk_table),
          risk.table.y.text.col = TRUE,
          conf.int              = isTRUE(input$surv_opt_conf_int),
          xlim                  = c(x_start, NA),
          xlab                  = if (nzchar(xlab_custom)) xlab_custom else "Time (Days)",
          ylab                  = if (nzchar(ylab_custom)) ylab_custom else "Survival Probability",
          legend.title          = if (nzchar(leg_label)) leg_label else "Treatment",
          title                 = plot_title
        ),
        error = function(e) rv$surv_results$km_plot
      )
      req(!is.null(km_obj))
      if (inherits(km_obj, "ggsurvplot")) {
        km_obj$plot <- km_obj$plot +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = title_hjust))
        if (nzchar(trimws(plot_sub))) {
          km_obj$plot <- km_obj$plot +
            ggplot2::labs(subtitle = plot_sub) +
            ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = title_hjust))
        }
        if (startsWith(legend_pos, "inplot_")) {
          coord <- switch(legend_pos,
            inplot_tl = c(0.02, 0.98), inplot_tr = c(0.98, 0.98),
            inplot_bl = c(0.02, 0.02), inplot_br = c(0.98, 0.02), c(0.5, 0.5))
          just <- switch(legend_pos,
            inplot_tl = c("left","top"),  inplot_tr = c("right","top"),
            inplot_bl = c("left","bottom"), inplot_br = c("right","bottom"),
            c("center","center"))
          km_obj$plot <- km_obj$plot +
            ggplot2::theme(
              legend.position    = coord,
              legend.justification = just,
              legend.background  = ggplot2::element_rect(
                fill = scales::alpha("white", 0.8), colour = NA))
        } else {
          km_obj$plot <- km_obj$plot + ggplot2::theme(legend.position = legend_pos)
        }
      }
      png(file, width  = input$surv_plot_width  %||% 800,
                height = input$surv_plot_height %||% 600,
                units  = "px", res = input$surv_plot_dpi %||% 100)
      print(km_obj)
      dev.off()
    }
  )
  
  output$download_surv_cox <- downloadHandler(
    filename = "survival_results.csv",
    content = function(file) {
      if (!is.null(rv$surv_results$results)) {
        write.csv(rv$surv_results$results, file, row.names = FALSE)
      }
    }
  )
  
  output$download_surv_risk_table <- downloadHandler(
    filename = "survival_risk_table.csv",
    content = function(file) {
      req(rv$surv_results$km_fit)
      km_sum <- summary(rv$surv_results$km_fit)
      risk_df <- data.frame(
        Group     = sub("^Treatment=", "", as.character(km_sum$strata)),
        Time      = km_sum$time,
        N_Risk    = km_sum$n.risk,
        N_Events  = km_sum$n.event,
        N_Censored = km_sum$n.censor,
        Survival  = round(km_sum$surv, 4)
      )
      write.csv(risk_df, file, row.names = FALSE)
    }
  )
  
  output$download_surv_forest_plot <- downloadHandler(
    filename = "forest_plot.png",
    content = function(file) {
      req(rv$surv_results, !isTRUE(rv$surv_results$error), rv$surv_results$results)
      p <- build_surv_forest_plot(
        res      = rv$surv_results$results,
        ref_grp  = rv$surv_results$reference_group,
        pt_color = input$surv_forest_pt_color %||% "steelblue",
        show_hr  = isTRUE(input$surv_forest_show_hr %||% TRUE),
        show_ci  = isTRUE(input$surv_forest_show_ci %||% TRUE),
        show_p   = isTRUE(input$surv_forest_show_p  %||% TRUE)
      )
      req(!is.null(p))
      p <- apply_plot_labels(p,
                             input$surv_plot_title        %||% "",
                             input$surv_plot_subtitle     %||% "",
                             input$surv_plot_title_align  %||% "center",
                             input$surv_opt_text_size     %||% 16,
                             xlab       = input$surv_plot_xaxis_label  %||% "",
                             ylab       = input$surv_plot_yaxis_label  %||% "",
                             legend_lab = input$surv_plot_legend_label %||% "")
      p <- apply_legend_pos(p, input$surv_opt_legend_pos %||% "right")
      png(file, width = input$surv_plot_width %||% 900,
                height = input$surv_plot_height %||% 600,
                units = "px", res = input$surv_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )
  
  })
}
