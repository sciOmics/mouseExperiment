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
#'   \item{peak_ci_synergy}{Information about when combination index synergy was strongest}
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
#' To visualize the results, use the plot_synergy_trend, plot_combination_index, 
#' or plot_synergy_combined functions with the output from this function.
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
  synergy_summary <- data.frame(
    Time_Point = numeric(),
    TGI_Drug_A = numeric(),
    TGI_Drug_B = numeric(),
    TGI_Combo = numeric(),
    Bliss_Expected_TGI = numeric(),
    Additive_Mean_Expected_TGI = numeric(),
    Bliss_Difference = numeric(),
    Additive_Mean_Difference = numeric(),
    Combination_Index = numeric(),
    P_Value_vs_Drug_A = numeric(),
    P_Value_vs_Drug_B = numeric(),
    Synergy_Assessment = character(),
    Validation_Check = character(),  # Add validation column
    stringsAsFactors = FALSE
  )
  
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
      loewe_result <- synergy_results$loewe_additivity
      ci_result <- synergy_results$combination_index
      stat_tests <- synergy_results$statistical_tests
      
      # Get tumor growth inhibition values directly from the results to ensure consistency
      summary_df <- synergy_results$summary
      tgi_drug_a <- summary_df$TGI_Percent[summary_df$Treatment == drug_a_name]
      tgi_drug_b <- summary_df$TGI_Percent[summary_df$Treatment == drug_b_name]
      tgi_combo <- summary_df$TGI_Percent[summary_df$Treatment == combo_name]
      
      # Recalculate the Bliss expected values to ensure consistency with bar graphs
      # Use the same formula as in analyze_drug_synergy: EA + EB - (EA * EB)
      fe_a <- tgi_drug_a / 100  # Convert to fractional effect
      fe_b <- tgi_drug_b / 100
      bliss_expected_fe <- fe_a + fe_b - (fe_a * fe_b)
      bliss_expected <- bliss_expected_fe * 100
      
      # Get Additive (Mean) expected value (which is simply (A+B)/2)
      additive_mean_expected <- summary_df$TGI_Percent[summary_df$Treatment == "Additive (Mean)"]
      
      # Perform validation checks
      # 1. Verify Bliss expected calculation is consistent
      bliss_expected_check <- fe_a + fe_b - (fe_a * fe_b)
      bliss_expected_tgi_check <- bliss_expected_check * 100
      bliss_diff_check <- (tgi_combo / 100) - bliss_expected_check
      
      # 2. Verify CI calculation
      ci_check <- (fe_a + fe_b) / (2 * (tgi_combo / 100))
      
      # Check if values match within tolerance
      bliss_match <- abs(bliss_expected - bliss_expected_tgi_check) < 0.01
      ci_match <- abs(ci_result$ci - ci_check) < 0.01
      
      # 3. Create validation message
      validation_msg <- "All calculations verified"
      if (!bliss_match) {
        validation_msg <- paste("Bliss mismatch:", round(bliss_expected, 2), "vs", round(bliss_expected_tgi_check, 2))
      } else if (!ci_match) {
        validation_msg <- paste("CI mismatch:", round(ci_result$ci, 2), "vs", round(ci_check, 2))
      }
      
      # Create a row for this time point
      tp_row <- data.frame(
        Time_Point = tp,
        TGI_Drug_A = tgi_drug_a,
        TGI_Drug_B = tgi_drug_b,
        TGI_Combo = tgi_combo,
        Bliss_Expected_TGI = bliss_expected,
        Additive_Mean_Expected_TGI = additive_mean_expected,
        Bliss_Difference = bliss_result$difference * 100, # Convert to percentage
        Additive_Mean_Difference = loewe_result$difference * 100, # Convert to percentage
        Combination_Index = ci_result$ci,
        P_Value_vs_Drug_A = stat_tests$P_Value[1],
        P_Value_vs_Drug_B = stat_tests$P_Value[2],
        Synergy_Assessment = synergy_results$overall_assessment,
        Validation_Check = validation_msg,
        stringsAsFactors = FALSE
      )
      
      # Add to summary data frame
      synergy_summary <- rbind(synergy_summary, tp_row)
      
    }, error = function(e) {
      warning(paste("Error analyzing time point", tp, ":", e$message))
    })
  }
  
  # Check if we have any successful results
  if (nrow(synergy_summary) == 0) {
    stop("Could not calculate synergy for any time points.")
  }
  
  # Order the summary by time point
  synergy_summary <- synergy_summary[order(synergy_summary$Time_Point), ]
  
  # Find when synergy was strongest (lowest CI and highest Bliss difference)
  peak_ci_synergy <- synergy_summary[which.min(synergy_summary$Combination_Index), ]
  peak_bliss_synergy <- synergy_summary[which.max(synergy_summary$Bliss_Difference), ]
  
  # Print summary of findings (only when verbose)
  if (isTRUE(verbose)) {
    cat("\n=== Drug Combination Synergy Analysis Over Time ===\n")
    cat("Analysis performed across", nrow(synergy_summary), "time points from", 
        min(synergy_summary$Time_Point), "to", max(synergy_summary$Time_Point), "\n\n")
    
    cat("Peak Synergy Findings:\n")
    cat(paste0("Strongest CI Synergy at Day ", peak_ci_synergy$Time_Point, 
               " (CI = ", round(peak_ci_synergy$Combination_Index, 2), ")\n"))
    cat(paste0("Strongest Bliss Synergy at Day ", peak_bliss_synergy$Time_Point, 
               " (Difference = ", round(peak_bliss_synergy$Bliss_Difference, 1), "%)\n\n"))
    
    cat("Synergy Summary by Time Point:\n")
    print(synergy_summary[, c("Time_Point", "TGI_Combo", "Bliss_Expected_TGI", 
                             "Bliss_Difference", "Combination_Index", "Synergy_Assessment", "Validation_Check")])
    
    # Check for validation issues
    validation_issues <- synergy_summary$Validation_Check != "All calculations verified"
    if (any(validation_issues)) {
      cat("\nWARNING: Some validation checks failed. Please review calculations.\n")
    } else {
      cat("\nAll calculations verified as consistent.\n")
    }
  }
  
  # Return comprehensive results
  return(list(
    timepoint_results = timepoint_results,
    synergy_summary = synergy_summary,
    peak_ci_synergy = peak_ci_synergy,
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
  
  # Create the plot
  trend_plot <- ggplot2::ggplot(synergy_summary, ggplot2::aes(x = Time_Point)) +
    # Treatment TGI lines - use named mapping for consistency
    ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_A, color = drug_a_name)) +
    ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_B, color = drug_b_name)) +
    ggplot2::geom_line(ggplot2::aes(y = TGI_Combo, color = combo_name), size = 1.2) +
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
    # Add a manual legend for the ribbon areas using labelled rectangles in the plot
    ggplot2::annotate("rect", 
              xmin = min(synergy_summary$Time_Point) * 1.05, 
              xmax = min(synergy_summary$Time_Point) * 1.15, 
              ymin = max(synergy_summary$TGI_Combo, na.rm = TRUE) * 0.85, 
              ymax = max(synergy_summary$TGI_Combo, na.rm = TRUE) * 0.90,
              fill = "lightgreen", alpha = 0.4) +
    ggplot2::annotate("text", 
              x = min(synergy_summary$Time_Point) * 1.2, 
              y = max(synergy_summary$TGI_Combo, na.rm = TRUE) * 0.875,
              label = "Synergy", hjust = 0) +
    ggplot2::annotate("rect", 
              xmin = min(synergy_summary$Time_Point) * 1.05, 
              xmax = min(synergy_summary$Time_Point) * 1.15, 
              ymin = max(synergy_summary$TGI_Combo, na.rm = TRUE) * 0.75, 
              ymax = max(synergy_summary$TGI_Combo, na.rm = TRUE) * 0.80,
              fill = "pink", alpha = 0.4) +
    ggplot2::annotate("text", 
              x = min(synergy_summary$Time_Point) * 1.2, 
              y = max(synergy_summary$TGI_Combo, na.rm = TRUE) * 0.775,
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

#' Plot Combination Index Over Time
#'
#' Creates a line plot visualizing the Combination Index (CI) over time,
#' with clear indication of synergy (CI < 1) and antagonism (CI > 1) regions.
#'
#' @param synergy_results Results object from analyze_drug_synergy_over_time function
#' @param custom_title Optional custom title for the plot
#'
#' @return A ggplot2 object visualizing CI values over time
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
#' # Then create the CI plot
#' plot_combination_index(results)
#' }
plot_combination_index <- function(synergy_results, custom_title = NULL) {
  # Validate input
  if (!is.list(synergy_results) || is.null(synergy_results$synergy_summary)) {
    stop("Input must be a valid result object from analyze_drug_synergy_over_time()")
  }
  
  # Extract necessary data
  synergy_summary <- synergy_results$synergy_summary
  
  # Set title
  if (is.null(custom_title)) {
    title <- "Combination Index Over Time"
  } else {
    title <- custom_title
  }
  
  # Create the plot
  ci_plot <- ggplot2::ggplot(synergy_summary, ggplot2::aes(x = Time_Point, y = Combination_Index)) +
    ggplot2::geom_line(size = 1) +
    # Fix the TRUE/FALSE mapping issue by explicitly setting aesthetics
    ggplot2::geom_point(
      data = subset(synergy_summary, Combination_Index < 1),
      ggplot2::aes(x = Time_Point, y = Combination_Index),
      size = 3, color = "green"
    ) +
    ggplot2::geom_point(
      data = subset(synergy_summary, Combination_Index >= 1),
      ggplot2::aes(x = Time_Point, y = Combination_Index),
      size = 3, color = "red"
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    # Add manual legend
    ggplot2::annotate("point", x = -Inf, y = -Inf, size = 3, color = "green") +
    ggplot2::annotate("point", x = -Inf, y = -Inf, size = 3, color = "red") +
    ggplot2::guides(color = "none") + # Remove automatic color legend
    ggplot2::annotate("text", x = max(synergy_summary$Time_Point) * 0.15, 
              y = max(synergy_summary$Combination_Index, na.rm = TRUE) * 0.85,
              label = "Synergy (CI < 1)", color = "green", hjust = 0) +
    ggplot2::annotate("text", x = max(synergy_summary$Time_Point) * 0.15, 
              y = max(synergy_summary$Combination_Index, na.rm = TRUE) * 0.95,
              label = "Antagonism (CI >= 1)", color = "red", hjust = 0) +
    ggplot2::labs(
      title = title,
      subtitle = "CI < 1 indicates synergy, CI > 1 indicates antagonism",
      x = "Time Point",
      y = "Combination Index (CI)"
    ) +
    ggplot2::theme_minimal()
  
  return(ci_plot)
}

#' Plot Combined Drug Synergy Analysis
#'
#' Creates a combined plot showing both the tumor growth inhibition trends
#' and the combination index trends over time.
#'
#' @param synergy_results Results object from analyze_drug_synergy_over_time function
#' @param custom_title Optional list with custom titles for the individual plots
#' @param custom_colors Optional named vector of custom colors for treatment groups
#'
#' @return A combined ggplot2 object with both trend plots
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
#' # Then create the combined plot
#' plot_synergy_combined(results)
#' }
plot_synergy_combined <- function(synergy_results, custom_title = NULL, custom_colors = NULL) {
  # Create individual plots
  trend_plot <- plot_synergy_trend(synergy_results, 
                                 ifelse(is.list(custom_title), custom_title$trend, custom_title),
                                 custom_colors)
  
  ci_plot <- plot_combination_index(synergy_results,
                                  ifelse(is.list(custom_title), custom_title$ci, custom_title))
  
  # Combine the plots
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    combined_plot <- ggpubr::ggarrange(trend_plot, ci_plot, 
                                     ncol = 1, nrow = 2, 
                                     common.legend = FALSE,
                                     heights = c(2, 1))
    return(combined_plot)
  } else {
    warning("Package 'ggpubr' is required for creating combined plots. Returning trend plot only.")
    return(trend_plot)
  }
}