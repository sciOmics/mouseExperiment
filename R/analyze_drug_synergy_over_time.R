#' Analyze Drug Combination Synergy Over Time
#'
#' This function tests for synergistic effects of drug combinations in tumor growth data
#' over multiple time points. It extends the analyze_drug_synergy function by evaluating
#' how synergy metrics change over the course of the experiment, providing insights
#' into when synergy may be strongest.
#'
#' @param df A data frame containing tumor growth data.
#' @param treatment_column A character string specifying the column name for treatment groups. Default is "Treatment".
#' @param volume_column A character string specifying the column name for tumor volume measurements. Default is "Volume".
#' @param time_column A character string specifying the column name for time points. Default is "Day".
#' @param drug_a_name A character string specifying the name of the first single agent treatment group.
#' @param drug_b_name A character string specifying the name of the second single agent treatment group.
#' @param combo_name A character string specifying the name of the combination treatment group.
#' @param control_name A character string specifying the name of the control/vehicle group. Default is "Control".
#' @param min_time_point Optional. A numeric value specifying the minimum time point to include in analysis. Default is NULL.
#' @param max_time_point Optional. A numeric value specifying the maximum time point to include in analysis. Default is NULL.
#' @param verbose Logical. If TRUE, prints detailed results to the console.
#'        Default is TRUE for interactive use; set to FALSE for programmatic/dashboard use.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{timepoint_results}{A list of results at each time point, each containing the outputs from analyze_drug_synergy}
#'   \item{synergy_summary}{A data frame summarizing synergy metrics across all time points}
#'   \item{peak_bliss_synergy}{Information about when Bliss synergy was strongest}
#'   \item{drug_a_name, drug_b_name, combo_name}{Names of treatment groups for plotting}
#' }
#'
#' @details
#' This function applies the analyze_drug_synergy approach to each time point in the data,
#' allowing for the assessment of how synergy develops over time. It calculates key metrics
#' including:
#'
#' 1. Bliss Independence effect differences at each time point
#' 2. Combination Index (CI) at each time point
#' 3. Statistical significance of combination advantage over monotherapies
#'
#' To visualize the results, use the plot_synergy_trend
#' function with the output from this function.
#'
#' @section Assumptions and Limitations:
#' \strong{Bliss Independence applied to TGI:} Bliss Independence was formulated for the
#' probability of cell death, not for proportional growth inhibition. Applying it to TGI is a
#' common pragmatic choice but carries a ceiling effect: when individual drug TGIs are large
#' (each > 50\%), the Bliss expected combined TGI approaches 100\%, making it nearly impossible
#' to demonstrate synergy by this criterion regardless of the true biological interaction.
#' Interpret Bliss results cautiously when individual-agent TGIs exceed 50\% at a given time point.
#'
#' \strong{No Combination Index.} The Loewe / CI path was removed in v0.21.0:
#' what it computed was response additivity rather than Loewe additivity, and it
#' labelled Bliss-additive combinations antagonistic across most of the range
#' where active single agents sit. See \code{\link{analyze_drug_synergy}}.
#'
#' @examples
#' # Analyze synergy over all available time points
#' data(combo_treatment_synthetic_data)
#' data_processed <- calculate_volume(combo_treatment_synthetic_data)
#' data_processed <- calculate_dates(data_processed, start_date = "03/24/2025")
#' 
#' synergy_results <- analyze_drug_synergy_over_time(
#'   df = data_processed,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Combo",
#'   control_name = "Control"
#' )
#' 
#' # Access synergy metrics at each time point
#' print(head(synergy_results$synergy_summary))
#' trend_plot <- plot_synergy_trend(synergy_results)
#' print(trend_plot)
#' 
#' # Check when synergy was strongest
#' print(synergy_results$peak_bliss_synergy)
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @export
analyze_drug_synergy_over_time <- function(df, 
                                      treatment_column = "Treatment",
                                      volume_column = "Volume",
                                      time_column = "Day",
                                      drug_a_name,
                                      drug_b_name,
                                      combo_name,
                                      control_name = "Control",
                                      min_time_point = NULL,
                                      max_time_point = NULL,
                                      verbose = TRUE) {
  
  # Input validation
  required_columns <- c(treatment_column, volume_column, time_column)
  missing_cols <- required_columns[!required_columns %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check that the specified groups exist in the data
  all_groups <- c(drug_a_name, drug_b_name, combo_name, control_name)
  missing_groups <- all_groups[!all_groups %in% unique(df[[treatment_column]])]
  
  if (length(missing_groups) > 0) {
    stop("The following specified groups do not exist in the treatment column: ", 
         paste(missing_groups, collapse = ", "))
  }
  
  # Get all time points
  all_timepoints <- sort(unique(df[[time_column]]))
  
  # Apply time point filters if specified
  if (!is.null(min_time_point)) {
    all_timepoints <- all_timepoints[all_timepoints >= min_time_point]
  }
  if (!is.null(max_time_point)) {
    all_timepoints <- all_timepoints[all_timepoints <= max_time_point]
  }
  
  # Check if we have any time points to analyze
  if (length(all_timepoints) == 0) {
    stop("No time points available for analysis after applying filters.")
  }
  
  message(paste("Analyzing drug synergy across", length(all_timepoints), "time points..."))
  
  # Initialize lists to store results
  timepoint_results <- list()
  synergy_rows <- list()
  
  # Loop through each time point and calculate synergy
  for (tp in all_timepoints) {
    # Create subset of data for this time point
    tp_data <- df[df[[time_column]] == tp, ]
    
    # Check if we have data for all required groups at this time point
    groups_at_tp <- unique(tp_data[[treatment_column]])
    missing_at_tp <- all_groups[!all_groups %in% groups_at_tp]
    
    if (length(missing_at_tp) > 0) {
      warning(paste("Time point", tp, "is missing data for groups:", paste(missing_at_tp, collapse = ", "), 
                   "- skipping this time point."))
      next
    }
    
    # Try to calculate synergy for this time point
    tryCatch({
      # Run the synergy analysis for this time point
      synergy_results <- analyze_drug_synergy(
        df = df,
        treatment_column = treatment_column,
        volume_column = volume_column,
        time_column = time_column,
        drug_a_name = drug_a_name,
        drug_b_name = drug_b_name,
        combo_name = combo_name,
        control_name = control_name,
        eval_time_point = tp,
        verbose = FALSE
      )
      
      # Store full results
      timepoint_results[[as.character(tp)]] <- synergy_results
      
      # Extract key metrics for summary
      bliss_result <- synergy_results$bliss_independence
      stat_tests <- synergy_results$statistical_tests
      
      # Get tumor growth inhibition values directly from the results to ensure consistency
      summary_df <- synergy_results$summary
      tgi_drug_a <- summary_df$TGI_Percent[summary_df$Treatment == drug_a_name]
      tgi_drug_b <- summary_df$TGI_Percent[summary_df$Treatment == drug_b_name]
      tgi_combo <- summary_df$TGI_Percent[summary_df$Treatment == combo_name]
      
      # Bliss expected value, read straight from summary_df
      bliss_expected <- summary_df$TGI_Percent[summary_df$Treatment == "Bliss Expected"]

      # Accumulate row for this time point
      synergy_rows[[length(synergy_rows) + 1]] <- data.frame(
        Time_Point = tp,
        TGI_Drug_A = tgi_drug_a,
        TGI_Drug_B = tgi_drug_b,
        TGI_Combo = tgi_combo,
        Bliss_Expected_TGI = bliss_expected,
        Bliss_Difference = bliss_result$difference * 100,
        P_Value_vs_Drug_A = stat_tests$P_Value[1],
        P_Value_vs_Drug_B = stat_tests$P_Value[2],
        Synergy_Assessment = synergy_results$overall_assessment,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      warning(paste("Error analyzing time point", tp, ":", e$message))
    })
  }
  
  # Combine accumulated rows into summary data frame
  synergy_summary <- if (length(synergy_rows) > 0) do.call(rbind, synergy_rows) else data.frame()

  # Check if we have any successful results
  if (nrow(synergy_summary) == 0) {
    stop("Could not calculate synergy for any time points.")
  }
  
  # Order the summary by time point
  synergy_summary <- synergy_summary[order(synergy_summary$Time_Point), ]
  
  # Find when Bliss synergy was strongest
  peak_bliss_synergy <- synergy_summary[which.max(synergy_summary$Bliss_Difference), ]
  
  # Print summary of findings (only when verbose)
  if (isTRUE(verbose)) {
    message("\n=== Drug Combination Synergy Analysis Over Time ===")
    message("Analysis performed across ", nrow(synergy_summary), " time points from ", 
        min(synergy_summary$Time_Point), " to ", max(synergy_summary$Time_Point), "\n")
    
    message("Peak Synergy Findings:")
    message("Strongest Bliss Synergy at Day ", peak_bliss_synergy$Time_Point, 
               " (Difference = ", round(peak_bliss_synergy$Bliss_Difference, 1), "%)\n")
    
    message("Synergy Summary by Time Point:")
    message(paste(utils::capture.output(
      print(synergy_summary[, c("Time_Point", "TGI_Combo", "Bliss_Expected_TGI",
                             "Bliss_Difference", "Synergy_Assessment")])
    ), collapse = "\n"))
  }
  
  # Return comprehensive results
  return(list(
    timepoint_results = timepoint_results,
    synergy_summary = synergy_summary,
    peak_bliss_synergy = peak_bliss_synergy,
    # Add these for plotting functions
    drug_a_name = drug_a_name,
    drug_b_name = drug_b_name,
    combo_name = combo_name
  ))
}

#' Plot Drug Synergy Trend Over Time
#'
#' Creates a line plot visualizing tumor growth inhibition (TGI) trends over time
#' for different treatment groups, highlighting synergy and antagonism regions.
#'
#' @param synergy_results Results object from analyze_drug_synergy_over_time function
#' @param custom_title Optional custom title for the plot
#' @param custom_colors Optional named vector of custom colors for treatment groups
#'
#' @return A ggplot2 object visualizing TGI trends over time
#' @export
#'
#' @examples
#' \dontrun{
#' # First run the analysis
#' results <- analyze_drug_synergy_over_time(
#'   df = tumor_data,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Drug A + Drug B"
#' )
#' 
#' # Then create the trend plot
#' plot_synergy_trend(results)
#' }
plot_synergy_trend <- function(synergy_results, custom_title = NULL, custom_colors = NULL) {
  # Validate input
  if (!is.list(synergy_results) || is.null(synergy_results$synergy_summary)) {
    stop("Input must be a valid result object from analyze_drug_synergy_over_time()")
  }
  
  # Extract necessary data
  synergy_summary <- synergy_results$synergy_summary
  drug_a_name <- synergy_results$drug_a_name
  drug_b_name <- synergy_results$drug_b_name
  combo_name <- synergy_results$combo_name
  
  # Set title
  if (is.null(custom_title)) {
    title <- "Tumor Growth Inhibition and Synergy Over Time"
  } else {
    title <- custom_title
  }
  
  # Create a named color vector to ensure correct assignment
  if (is.null(custom_colors)) {
    color_values <- c("red", "blue", "purple", "gray50")
    names(color_values) <- c(drug_a_name, drug_b_name, combo_name, "Bliss Expected")
  } else {
    color_values <- custom_colors
  }
  
  # Annotation coordinates: additive offsets so they work when Time_Point starts at 0
  .tp_lo   <- min(synergy_summary$Time_Point)
  .tp_span <- max(max(synergy_summary$Time_Point) - .tp_lo, 1)
  .tp_x1   <- .tp_lo + .tp_span * 0.05
  .tp_x2   <- .tp_lo + .tp_span * 0.15
  .tp_x3   <- .tp_lo + .tp_span * 0.20
  .tgi_max <- max(synergy_summary$TGI_Combo, na.rm = TRUE)

  # Create the plot
  trend_plot <- ggplot2::ggplot(synergy_summary, ggplot2::aes(x = Time_Point)) +
    # Treatment TGI lines - use named mapping for consistency
    ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_A, color = drug_a_name)) +
    ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_B, color = drug_b_name)) +
    ggplot2::geom_line(ggplot2::aes(y = TGI_Combo, color = combo_name), linewidth = 1.2) +
    ggplot2::geom_line(ggplot2::aes(y = Bliss_Expected_TGI, color = "Bliss Expected"), linetype = "dashed") +
    # Synergy area (when combo effect > bliss expected)
    ggplot2::geom_ribbon(data = subset(synergy_summary, TGI_Combo > Bliss_Expected_TGI),
                      ggplot2::aes(ymin = Bliss_Expected_TGI, ymax = TGI_Combo),
                      fill = "lightgreen", alpha = 0.4) +
    # Antagonism area (when combo effect < bliss expected)
    ggplot2::geom_ribbon(data = subset(synergy_summary, TGI_Combo < Bliss_Expected_TGI),
                      ggplot2::aes(ymin = TGI_Combo, ymax = Bliss_Expected_TGI),
                      fill = "pink", alpha = 0.4) +
    # Formatting
    ggplot2::scale_color_manual(
      name = "Treatment",
      values = color_values,
      labels = c(drug_a_name, drug_b_name, combo_name, "Bliss Expected")
    ) +
    # Manual legend for ribbon areas
    ggplot2::annotate("rect",
              xmin = .tp_x1, xmax = .tp_x2,
              ymin = .tgi_max * 0.85, ymax = .tgi_max * 0.90,
              fill = "lightgreen", alpha = 0.4) +
    ggplot2::annotate("text",
              x = .tp_x3, y = .tgi_max * 0.875,
              label = "Synergy", hjust = 0) +
    ggplot2::annotate("rect",
              xmin = .tp_x1, xmax = .tp_x2,
              ymin = .tgi_max * 0.75, ymax = .tgi_max * 0.80,
              fill = "pink", alpha = 0.4) +
    ggplot2::annotate("text",
              x = .tp_x3, y = .tgi_max * 0.775,
              label = "Antagonism", hjust = 0) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Green area = synergy (", combo_name, " > Bliss Expected)\n",
                       "Pink area = antagonism (", combo_name, " < Bliss Expected)"),
      x = "Time Point",
      y = "Tumor Growth Inhibition (%)"
    ) +
    ggplot2::theme_minimal()
  
  return(trend_plot)
}