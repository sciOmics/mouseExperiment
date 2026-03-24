#' Survival Analysis UI
#' @keywords internal
survivalUI <- function() {
  fluidPage(
    h2("Survival Analysis"),
    
    p("Analyze time-to-event data using Kaplan-Meier survival curves and proportional hazards models. ",
      "Ensure data is uploaded and configured in the Data Upload tab."),
    
    fluidRow(
      column(3,
        wellPanel(
          h4("Analysis Parameters"),
          
          selectInput("surv_palette",
                     "Color Palette",
                     choices = c("Default" = "hue",
                                "Colorblind" = "colorblind",
                                "Viridis" = "viridis",
                                "Set1" = "Set1",
                                "Set2" = "Set2"),
                     selected = "hue"),
          
          uiOutput("surv_reference_group_ui"),
          
          br(),
          
          actionButton("run_surv_analysis", 
                      "Run Analysis",
                      class = "btn-lg btn-primary btn-block"),
          
          uiOutput("surv_status_msg")
        )
      ),
      
      column(9,
        tabsetPanel(
          id = "surv_results_tabs",
          
          tabPanel("Summary",
            tab_info(strong("Survival summary: "), "Per-group event counts, total subjects, event rate, and median survival time (the day by which 50\u202f% of animals have reached the endpoint). \u2018Not reached\u2019 means fewer than half the group experienced an event by the end of the study."),
            shinycssloaders::withSpinner(
              verbatimTextOutput("surv_summary")
            )
          ),
          
          tabPanel("Cox Model",
            tab_info(strong("Cox proportional hazards model. "),
              "\u2018Hazard Ratio (HR)\u2019 is relative to the reference group. ",
              strong("HR\u202f<\u202f1"), " = lower hazard (better survival); ",
              strong("HR\u202f>\u202f1"), " = higher hazard. ",
              "95\u202f% CI not crossing 1.0 indicates statistical significance. ",
              "Firth\u2019s bias-reduced method is applied automatically when perfect separation is detected (e.g. a group with zero events)."),
            shinycssloaders::withSpinner(
              verbatimTextOutput("surv_cox_summary")
            ),
            shinycssloaders::withSpinner(
              DTOutput("surv_cox_table")
            ),
            downloadButton("download_surv_cox", "Download Results")
          ),
          
          tabPanel("Plots",
            fluidRow(
              column(3,
                wellPanel(
                  selectInput("surv_plot_selector", "Select Plot",
                    choices = c(
                      "Kaplan-Meier Curve"       = "km",
                      "Hazard Ratio Forest Plot" = "forest"
                    )
                  ),
                  hr(),
                  h6(icon("palette"), " Appearance"),
                  conditionalPanel("input.surv_plot_selector == 'km'",
                    uiOutput("surv_color_pickers_ui"),
                    checkboxInput("surv_opt_conf_int",
                                  "Show confidence interval", value = FALSE)
                  ),
                  conditionalPanel("input.surv_plot_selector == 'forest'",
                    colourpicker::colourInput("surv_forest_pt_color",
                                             "Point colour", value = "steelblue")
                  )
                )
              ),
              column(9,
                conditionalPanel("input.surv_plot_selector == 'km'",
                  h4("Kaplan-Meier Survival Curves"),
                  tab_info(
                    strong("Kaplan-Meier survival curve: "),
                    "non-parametric step-function estimating the probability of ",
                    "surviving past each timepoint. Each drop represents one or more ",
                    "events; tick marks represent censored animals still alive at last ",
                    "observation. Groups that remain higher on the y-axis have better survival."
                  ),
                  shinycssloaders::withSpinner(
                    plotOutput("surv_plot_km", height = "500px")
                  ),
                  downloadButton("download_surv_plot_km", "Download Plot")
                ),
                conditionalPanel("input.surv_plot_selector == 'forest'",
                  h4("Hazard Ratio Forest Plot"),
                  tab_info(
                    strong("Hazard ratio forest plot: "),
                    "each row is a treatment group vs the reference. ",
                    "The x-axis uses a log\u2082 scale, symmetric around HR\u202f=\u202f1. ",
                    "Points to the left indicate a protective effect; points to the right ",
                    "indicate increased hazard relative to the reference group."
                  ),
                  shinycssloaders::withSpinner(
                    plotOutput("surv_forest_plot", height = "500px")
                  ),
                  downloadButton("download_surv_forest_plot", "Download Plot")
                )
              )
            )
          ),

          tabPanel("Risk Table",
            tab_info(strong("Numbers at risk table: "), "at each event time, shows how many animals per group are still at risk (N Risk), how many experienced the event (N Events), and how many were censored (N Censored). Declining N Risk over time reflects events and censoring."),
            shinycssloaders::withSpinner(
              DTOutput("surv_risk_table")
            ),
            downloadButton("download_surv_risk_table", "Download Table")
          )
        )
      )
    )
  )
}

#' Survival Analysis Server Logic
#' @keywords internal
survival_server <- function(input, output, session, rv) {
  
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
      if (!is.null(rv$tx_col) && rv$tx_col == "TreatmentGroup" &&
          !is.null(rv$dose_col) && nzchar(rv$dose_col) &&
          rv$dose_col %in% colnames(d) &&
          !is.null(input$col_treatment) && nzchar(input$col_treatment) &&
          input$col_treatment %in% colnames(d)) {
        combined <- paste(d[[input$col_treatment]], d[[rv$dose_col]], sep = " ")
        return(sort(unique(as.character(na.omit(combined)))))
      }
      if (!is.null(input$col_treatment) && nzchar(input$col_treatment) &&
          input$col_treatment %in% colnames(d)) {
        return(sort(unique(as.character(na.omit(d[[input$col_treatment]])))))
      }
      character(0)
    }

    treatment_groups <- get_groups()

    if (length(treatment_groups) > 0) {
      default_ref <- pick_reference_group(treatment_groups)
      return(
        selectInput("surv_reference_group",
                   "Reference Group",
                   choices = treatment_groups,
                   selected = default_ref)
      )
    }

    # Fallback if no data
    selectInput("surv_reference_group",
               "Reference Group",
               choices = "",
               selected = "")
  })
  
  # Run survival analysis
  observeEvent(input$run_surv_analysis, {
    req(rv$uploaded_data, rv$data_valid)
    
    # Validate that required columns are configured from data upload
    if (input$col_treatment == "" || input$col_treatment == "-- Select a reference group --") {
      output$surv_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error: Please select a Treatment Group on the Data Upload tab")
      })
      return()
    }
    

    if (is.null(input$col_survival_event) || input$col_survival_event == "" || input$col_survival_event == "None") {
      output$surv_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error: Please select an Event/Censor Status Column on the Data Upload tab")
      })
      return()
    }
    
    rv$surv_running <- TRUE
    output$surv_status_msg <- renderUI({
      div(class = "alert alert-info",
          icon("spinner", class = "fa-spin"), "Running survival analysis...")
    })
    
    # Prepare data - apply same time processing as tumor growth analysis
    data <- rv$uploaded_data

    # When a dose column is mapped the treatment column is set to "TreatmentGroup"
    # (a combined "Treatment Dose" label).  rv$uploaded_data never carries that
    # synthetic column, so we must build it here before any deduplication or
    # date calculations so it is preserved throughout the pipeline.
    if (!is.null(rv$dose_col) && nzchar(rv$dose_col) &&
        rv$dose_col %in% colnames(data) &&
        !is.null(rv$tx_col) && rv$tx_col == "TreatmentGroup") {
      data[["TreatmentGroup"]] <- paste(data[[input$col_treatment]],
                                        data[[rv$dose_col]], sep = " ")
    }

    # Handle time column - either use direct day or calculate from date
    # This mirrors the tumor growth analysis logic
    if (input$time_method == "from_date") {
      # Calculate days from reference date
      start_date <- input$reference_date
      
      # Detect date format from the actual data in the date column
      first_date <- as.character(data[1, input$col_date_for_day])
      
      # Determine date format by checking the pattern
      date_format <- NULL
      if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_date)) {
        date_format <- "%m/%d/%Y"  # MM/DD/YYYY format
      } else if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", first_date)) {
        date_format <- "%Y-%m-%d"  # YYYY-MM-DD format
      } else if (grepl("^\\d{1,2}-[A-Za-z]{3}$", first_date)) {
        date_format <- "%d-%b"     # DD-Mon format
      }
      
      data <- mouseExperiment::calculate_dates(
        data,
        start_date = start_date,
        date_column = input$col_date_for_day,
        date_format = date_format
      )
      time_col <- "Day"  # calculate_dates produces a "Day" column
    } else {
      time_col <- input$col_day
    }
    
    # For survival analysis, reduce data to one row per animal
    # Use the user-selected ID and Cage columns (not hardcoded names)
    id_col_surv   <- if (!is.null(input$col_id)   && nzchar(input$col_id))   input$col_id   else "ID"
    cage_col_surv <- if (!is.null(input$col_cage) && nzchar(input$col_cage) &&
                         input$col_cage != "None") input$col_cage else NULL

    has_id   <- !is.null(id_col_surv)   && id_col_surv   %in% colnames(data)
    has_cage <- !is.null(cage_col_surv) && cage_col_surv %in% colnames(data)

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
          div(class = "alert alert-danger",
              icon("exclamation-circle"),
              "Error:", result$message)
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
          div(class = "alert alert-danger",
              icon("exclamation-circle"), "Error: No data available for analysis")
        })
        return()
      }
    }
    
    # Run analysis with processed data and calculated/selected time column
    # Add diagnostic output
    cat("\n=== SURVIVAL ANALYSIS DIAGNOSTICS ===\n")
    cat("Data rows:", nrow(data), "\n")
    cat("Data columns:", paste(colnames(data), collapse = ", "), "\n")
    cat("Time column specified:", time_col, "\n")
    cat("Time column exists:", time_col %in% colnames(data), "\n")
    cat("Censor column specified:", input$col_survival_event, "\n")
    cat("Censor column exists:", input$col_survival_event %in% colnames(data), "\n")
    cat("Treatment column specified:", input$col_treatment, "\n")
    cat("Treatment column exists:", input$col_treatment %in% colnames(data), "\n")
    
    if (time_col %in% colnames(data)) {
      cat("Time values sample:", paste(head(data[[time_col]], 5), collapse = ", "), "\n")
      cat("Time NAs:", sum(is.na(data[[time_col]])), "\n")
    }
    
    if (input$col_survival_event %in% colnames(data)) {
      cat("Censor values sample:", paste(head(data[[input$col_survival_event]], 5), collapse = ", "), "\n")
      cat("Censor NAs:", sum(is.na(data[[input$col_survival_event]])), "\n")
    }
    cat("====================================\n\n")
    
    tx_col_surv <- if (!is.null(rv$tx_col) && nzchar(rv$tx_col)) rv$tx_col else input$col_treatment

    result <- safe_analysis({
      mouseExperiment::survival_statistics(
        df = data,
        time_column = time_col,
        censor_column = input$col_survival_event,
        treatment_column = tx_col_surv,
        cage_column = if (!is.null(input$col_cage) && nzchar(input$col_cage)) input$col_cage else "Cage",
        dose_column = rv$dose_col,
        reference_group = input$surv_reference_group,
        firth_correction = TRUE
      )
    })
    
    # Store time calculation info for later display
    result$time_method <- input$time_method
    result$time_column_used <- time_col
    result$reference_date <- input$reference_date
    
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
          censor_column     = input$col_survival_event,
          treatment_column  = tx_col_surv,
          id_column         = input$col_id,
          cage_column       = if (!is.null(input$col_cage) && nzchar(input$col_cage)) input$col_cage else NULL,
          dose_column       = rv$dose_col,
          colors            = NULL,
          show_risk_table   = TRUE,
          title             = "Kaplan-Meier Survival Curves by Treatment"
        )

        cat("KM plot generated successfully via plot_survival()\n")
      }, error = function(e) {
        cat("ERROR generating KM plot via plot_survival():", e$message, "\n")
        # Fallback to manual ggsurvplot if plot_survival() fails
        tryCatch({
          surv_df <- result$survival_data
          km_fit <- survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = surv_df)
          result$km_fit <<- km_fit
          km_plot_obj <- survminer::ggsurvplot(
            km_fit,
            data = surv_df,
            palette = input$surv_palette,
            risk.table = TRUE,
            risk.table.y.text.col = TRUE,
            conf.int = FALSE,
            xlab = "Time (Days)",
            ylab = "Survival Probability",
            title = "Kaplan-Meier Survival Curves by Treatment"
          )
          result$km_plot <<- km_plot_obj
          cat("KM plot generated via fallback ggsurvplot\n")
        }, error = function(e2) {
          cat("ERROR in KM fallback:", e2$message, "\n")
        })
      })
    }
    
    rv$surv_results <- result
    rv$surv_running <- FALSE

    # Store treatment group names for the per-group colour pickers in the Plots tab
    rv$surv_plot_groups <- if (!is.null(result$survival_data) && "Treatment" %in% colnames(result$survival_data)) {
      unique(as.character(result$survival_data$Treatment))
    } else NULL
    
    # Debug: Log the structure of the result
    cat("\n=== RESULT OBJECT STRUCTURE ===\n")
    cat("Class:", class(result), "\n")
    cat("Fields:", paste(names(result), collapse = ", "), "\n")
    cat("Has 'error' field:", "error" %in% names(result), "\n")
    cat("Has 'results' field:", "results" %in% names(result), "\n")
    cat("Has 'model' field:", "model" %in% names(result), "\n")
    cat("Has 'km_plot' field:", "km_plot" %in% names(result), "\n")
    if ("results" %in% names(result) && !is.null(result$results)) {
      cat("Results class:", class(result$results), "\n")
      cat("Results dimensions:", nrow(result$results), "rows x", ncol(result$results), "cols\n")
      cat("Results columns:", paste(colnames(result$results), collapse = ", "), "\n")
    }
    cat("====================================\n\n")
    
    if (isTRUE(result$error)) {
      output$surv_status_msg <- renderUI({
        div(class = "alert alert-danger",
            icon("exclamation-circle"), 
            "Error:", result$message)
      })
    } else {
      output$surv_status_msg <- renderUI({
        div(class = "alert alert-success",
            icon("check"), "Survival analysis completed successfully!")
      })
    }
  })
  
  # Summary output
  output$surv_summary <- renderText({
    req(rv$surv_results)
    cat("\n[surv_summary] Output function called\n")
    cat("[surv_summary] Result class:", class(rv$surv_results), "\n")
    cat("[surv_summary] Result fields:", paste(names(rv$surv_results), collapse = ", "), "\n")
    if (isTRUE(rv$surv_results$error)) {
      return(paste("Error:", rv$surv_results$message))
    }
    
    summary_text <- capture.output({
      cat("=== SURVIVAL ANALYSIS SUMMARY ===\n\n")
      
      # Show time specification method
      if (rv$surv_results$time_method == "from_date") {
        cat("Time Method: Calculated from Date Column\n")
        cat("Date Column:", input$col_date_for_day, "\n")
        cat("Reference Date:", input$reference_date, "\n")
      } else {
        cat("Time Method: Direct Day Column\n")
        cat("Day Column:", input$col_day, "\n")
      }
      
      cat("Event Column:", input$col_survival_event, "\n")
      cat("Treatment Column:", input$col_treatment, "\n")
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
  
  # Dynamic colour pickers for KM groups
  output$surv_color_pickers_ui <- renderUI({
    req(rv$surv_plot_groups)
    group_color_pickers_ui("surv_color", rv$surv_plot_groups)
  })

  # Kaplan-Meier plot — re-renders live when colour pickers or conf.int option change
  output$surv_plot_km <- renderPlot({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error)) return(NULL)
    req(rv$surv_results$km_fit, rv$surv_results$survival_data)

    surv_df <- rv$surv_results$survival_data
    km_fit  <- rv$surv_results$km_fit

    # Build colour vector from per-group colourInputs
    groups  <- rv$surv_plot_groups
    pal     <- build_color_vector(input, "surv_color", groups)
    if (is.null(pal)) pal <- "hue"

    tryCatch({
      km_obj <- survminer::ggsurvplot(
        km_fit,
        data                  = surv_df,
        palette               = pal,
        risk.table            = TRUE,
        risk.table.y.text.col = TRUE,
        conf.int              = isTRUE(input$surv_opt_conf_int),
        xlab                  = "Time (Days)",
        ylab                  = "Survival Probability",
        title                 = "Kaplan-Meier Survival Curves by Treatment"
      )
      print(km_obj)
    }, error = function(e) {
      # Fallback: print pre-generated stored KM plot
      if (!is.null(rv$surv_results$km_plot)) {
        print(rv$surv_results$km_plot)
      } else {
        plot(1, 1, type = "n", main = "KM plot not available", xlab = "", ylab = "")
      }
    })
  })
  
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
    
    datatable(results_df, options = list(scrollX = TRUE, pageLength = 10))
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
              options = list(scrollX = TRUE, pageLength = 15),
              rownames = FALSE)
  })
  
  # Log-rank test (method information)
  # Forest plot
  output$surv_forest_plot <- renderPlot({
    req(rv$surv_results)
    if (isTRUE(rv$surv_results$error) || is.null(rv$surv_results$results)) return(NULL)
    
    res <- rv$surv_results$results
    req(all(c("Group", "HR", "CI_Lower", "CI_Upper", "P_Value") %in% colnames(res)))
    
    ref_grp <- rv$surv_results$reference_group
    
    # Exclude reference group row if HR == 1 and CI_Lower == CI_Upper
    plot_df <- res[!is.na(res$HR), ]
    if (!is.null(ref_grp)) {
      plot_df <- plot_df[plot_df$Group != ref_grp, ]
    }
    req(nrow(plot_df) > 0)
    
    # Compute symmetric log2 axis limits
    all_hr_vals <- c(plot_df$HR, plot_df$CI_Lower, plot_df$CI_Upper)
    all_hr_vals <- all_hr_vals[!is.na(all_hr_vals) & all_hr_vals > 0]
    max_log2 <- max(abs(log2(all_hr_vals)), na.rm = TRUE)
    max_log2 <- ceiling(max_log2)  # round up to nearest integer
    x_lim <- c(2^(-max_log2), 2^max_log2)
    x_breaks <- 2^seq(-max_log2, max_log2, by = 1)
    x_labels <- ifelse(x_breaks >= 1,
                       as.character(x_breaks),
                       paste0("1/", as.character(1 / x_breaks)))
    
    ggplot2::ggplot(plot_df,
                   ggplot2::aes(x = HR, y = reorder(Group, HR),
                                xmin = CI_Lower, xmax = CI_Upper)) +
      ggplot2::geom_point(size = 5, color = input$surv_forest_pt_color %||% "steelblue") +
      ggplot2::geom_errorbarh(height = 0.25, linewidth = 0.8,
                              color = input$surv_forest_pt_color %||% "steelblue") +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
      ggplot2::geom_text(
        ggplot2::aes(label = paste0("HR: ", round(HR, 2), ", p = ", round(P_Value, 4))),
        vjust = 3, hjust = -0.1, size = 4, color = "grey30") +
      ggplot2::scale_x_continuous(
        trans  = "log2",
        limits = x_lim,
        breaks = x_breaks,
        labels = x_labels
      ) +
      ggplot2::labs(
        title = paste("Forest Plot: Hazard Ratios vs", ref_grp),
        x = "Hazard Ratio (log\u2082 scale)",
        y = NULL
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 11),
        legend.position = "none"
      )
  })
  
  # Download handlers
  output$download_surv_plot_km <- downloadHandler(
    filename = "kaplan_meier_plot.png",
    content = function(file) {
      if (!is.null(rv$surv_results$km_plot)) {
        png(file, width = 800, height = 600)
        print(rv$surv_results$km_plot)
        dev.off()
      }
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
      p <- local({
        res <- rv$surv_results$results
        ref_grp <- rv$surv_results$reference_group
        plot_df <- res[!is.na(res$HR), ]
        if (!is.null(ref_grp)) plot_df <- plot_df[plot_df$Group != ref_grp, ]
        all_hr_vals <- c(plot_df$HR, plot_df$CI_Lower, plot_df$CI_Upper)
        all_hr_vals <- all_hr_vals[!is.na(all_hr_vals) & all_hr_vals > 0]
        max_log2 <- ceiling(max(abs(log2(all_hr_vals)), na.rm = TRUE))
        x_lim    <- c(2^(-max_log2), 2^max_log2)
        x_breaks <- 2^seq(-max_log2, max_log2, by = 1)
        x_labels <- ifelse(x_breaks >= 1,
                           as.character(x_breaks),
                           paste0("1/", as.character(1 / x_breaks)))
        ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = HR, y = reorder(Group, HR),
                                    xmin = CI_Lower, xmax = CI_Upper)) +
          ggplot2::geom_point(size = 5, color = "steelblue") +
          ggplot2::geom_errorbarh(height = 0.25, linewidth = 0.8, color = "steelblue") +
          ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
          ggplot2::geom_text(
            ggplot2::aes(label = paste0("HR: ", round(HR, 2), ", p = ", round(P_Value, 4))),
            vjust = 3, hjust = -0.1, size = 4, color = "grey30") +
          ggplot2::scale_x_continuous(
            trans  = "log2",
            limits = x_lim,
            breaks = x_breaks,
            labels = x_labels
          ) +
          ggplot2::labs(title = paste("Forest Plot: Hazard Ratios vs", ref_grp),
                       x = "Hazard Ratio (log\u2082 scale)", y = NULL) +
          ggplot2::theme_classic(base_size = 14) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 11),
                        legend.position = "none")
      })
      png(file, width = 900, height = 600); print(p); dev.off()
    }
  )
  
}
