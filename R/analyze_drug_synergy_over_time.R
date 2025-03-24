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
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item{timepoint_results}{A list of results at each time point, each containing the outputs from analyze_drug_synergy}
#'   \item{synergy_trend}{A data frame summarizing synergy metrics across all time points}
#'   \item{trend_plot}{A ggplot2 object visualizing how synergy changes over time}
#'   \item{peak_synergy}{Information about when synergy was strongest}
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
#' The function generates visualization of synergy trends and identifies when synergy is
#' strongest during the course of treatment.
#'
#' @examples
#' # Analyze synergy over all available time points
#' synergy_results <- analyze_drug_synergy_over_time(
#'   df = tumor_data,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Drug A + Drug B",
#'   control_name = "Vehicle"
#' )
#' 
#' # Print the trend plot
#' print(synergy_results$trend_plot)
#' 
#' # Check when synergy was strongest
#' print(synergy_results$peak_synergy)
#'
#' @import ggplot2
#' @import dplyr
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
                                      max_time_point = NULL) {
  
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
    Loewe_Expected_TGI = numeric(),
    Bliss_Difference = numeric(),
    Loewe_Difference = numeric(),
    Combination_Index = numeric(),
    P_Value_vs_Drug_A = numeric(),
    P_Value_vs_Drug_B = numeric(),
    Synergy_Assessment = character(),
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
        eval_time_point = tp
      )
      
      # Store full results
      timepoint_results[[as.character(tp)]] <- synergy_results
      
      # Extract key metrics for summary
      bliss_result <- synergy_results$bliss_independence
      loewe_result <- synergy_results$loewe_additivity
      ci_result <- synergy_results$combination_index
      stat_tests <- synergy_results$statistical_tests
      
      # Get tumor growth inhibition values
      summary_df <- synergy_results$summary
      tgi_drug_a <- summary_df$TGI_Percent[summary_df$Treatment == drug_a_name]
      tgi_drug_b <- summary_df$TGI_Percent[summary_df$Treatment == drug_b_name]
      tgi_combo <- summary_df$TGI_Percent[summary_df$Treatment == combo_name]
      bliss_expected <- summary_df$TGI_Percent[summary_df$Treatment == "Bliss Expected"]
      loewe_expected <- summary_df$TGI_Percent[summary_df$Treatment == "Loewe Expected"]
      
      # Create a row for this time point
      tp_row <- data.frame(
        Time_Point = tp,
        TGI_Drug_A = tgi_drug_a,
        TGI_Drug_B = tgi_drug_b,
        TGI_Combo = tgi_combo,
        Bliss_Expected_TGI = bliss_expected,
        Loewe_Expected_TGI = loewe_expected,
        Bliss_Difference = bliss_result$difference * 100, # Convert to percentage
        Loewe_Difference = loewe_result$difference * 100, # Convert to percentage
        Combination_Index = ci_result$ci,
        P_Value_vs_Drug_A = stat_tests$P_Value[1],
        P_Value_vs_Drug_B = stat_tests$P_Value[2],
        Synergy_Assessment = synergy_results$overall_assessment,
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
  
  # Create synergy trend plot
  trend_plot <- ggplot2::ggplot(synergy_summary, ggplot2::aes(x = Time_Point)) +
    # Treatment TGI lines
    ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_A, color = drug_a_name)) +
    ggplot2::geom_line(ggplot2::aes(y = TGI_Drug_B, color = drug_b_name)) +
    ggplot2::geom_line(ggplot2::aes(y = TGI_Combo, color = combo_name), size = 1.2) +
    ggplot2::geom_line(ggplot2::aes(y = Bliss_Expected_TGI, color = "Bliss Expected"), linetype = "dashed") +
    # Synergy area
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Bliss_Expected_TGI, ymax = TGI_Combo), 
                      fill = "lightblue", alpha = 0.3) +
    # Formatting
    ggplot2::scale_color_manual(
      name = "Treatment",
      values = c("red", "blue", "purple", "gray50"),
      labels = c(drug_a_name, drug_b_name, combo_name, "Bliss Expected")
    ) +
    ggplot2::labs(
      title = "Tumor Growth Inhibition and Synergy Over Time",
      subtitle = "Shaded area shows synergy (difference between observed and expected effect)",
      x = "Time Point",
      y = "Tumor Growth Inhibition (%)"
    ) +
    ggplot2::theme_minimal()
  
  # Create Combination Index trend plot
  ci_plot <- ggplot2::ggplot(synergy_summary, ggplot2::aes(x = Time_Point, y = Combination_Index)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(size = 3, ggplot2::aes(color = Combination_Index < 1)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_color_manual(
      name = "Synergy",
      values = c("red", "green"),
      labels = c("Antagonism (CI > 1)", "Synergy (CI < 1)")
    ) +
    ggplot2::labs(
      title = "Combination Index Over Time",
      subtitle = "CI < 1 indicates synergy, CI > 1 indicates antagonism",
      x = "Time Point",
      y = "Combination Index (CI)"
    ) +
    ggplot2::theme_minimal()
  
  # Combine plots
  combined_plots <- ggpubr::ggarrange(trend_plot, ci_plot, 
                                     ncol = 1, nrow = 2, 
                                     common.legend = FALSE,
                                     heights = c(2, 1))
  
  # Find when synergy was strongest (lowest CI and highest Bliss difference)
  peak_ci_synergy <- synergy_summary[which.min(synergy_summary$Combination_Index), ]
  peak_bliss_synergy <- synergy_summary[which.max(synergy_summary$Bliss_Difference), ]
  
  # Print summary of findings
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
                           "Bliss_Difference", "Combination_Index", "Synergy_Assessment")])
  
  # Display the plots
  print(combined_plots)
  
  # Return comprehensive results
  return(list(
    timepoint_results = timepoint_results,
    synergy_summary = synergy_summary,
    trend_plot = trend_plot,
    ci_plot = ci_plot,
    combined_plot = combined_plots,
    peak_ci_synergy = peak_ci_synergy,
    peak_bliss_synergy = peak_bliss_synergy
  ))
}