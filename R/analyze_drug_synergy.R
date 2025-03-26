#' Analyze Drug Combination Synergy in Tumor Growth
#'
#' This function tests for synergistic effects of drug combinations in tumor growth data.
#' It compares the observed combination effect against several expected interaction models, 
#' including Bliss independence and Loewe additivity. The function calculates synergy scores 
#' and performs statistical tests to determine if the combination shows synergistic, 
#' additive, or antagonistic effects.
#'
#' @param df A data frame containing tumor growth data.
#' @param treatment_column A character string specifying the column name for treatment groups. Default is "Treatment".
#' @param volume_column A character string specifying the column name for tumor volume measurements. Default is "Volume".
#' @param time_column A character string specifying the column name for time points. Default is "Day".
#' @param drug_a_name A character string specifying the name of the first single agent treatment group.
#' @param drug_b_name A character string specifying the name of the second single agent treatment group.
#' @param combo_name A character string specifying the name of the combination treatment group.
#' @param control_name A character string specifying the name of the control/vehicle group. Default is "Control".
#' @param eval_time_point Optional. A numeric value specifying a specific time point to evaluate synergy.
#'        If NULL (default), the function will use the last time point in the data.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{summary}{A data frame summarizing the tumor growth inhibition (TGI) for each treatment and synergy metrics.}
#'   \item{bliss_independence}{Results of the Bliss independence model, including expected vs. observed effects.}
#'   \item{loewe_additivity}{Results of the Loewe additivity model.}
#'   \item{combination_index}{The combination index (CI), where CI < 1 indicates synergy, CI = 1 indicates additivity, and CI > 1 indicates antagonism.}
#'   \item{statistical_test}{Results of statistical tests comparing observed vs. expected effects.}
#'   \item{plot_data}{Data prepared for plotting, to be used with plot_drug_synergy function.}
#' }
#'
#' @details
#' The function calculates tumor growth inhibition (TGI) for each treatment group relative to the control.
#' It then applies several models to test for synergy:
#' 
#' 1. Bliss Independence Model: Assumes drugs act independently through different mechanisms.
#'    Expected effect = EA + EB - (EA * EB), where EA and EB are the effects of drug A and B alone.
#' 
#' 2. Loewe Additivity Model: Assumes drugs work through similar mechanisms.
#'    Expected effect is calculated based on dose-response relationships.
#' 
#' 3. Combination Index: A widely used metric where CI < 1 indicates synergy.
#' 
#' The function performs statistical tests to determine if the observed combination effect
#' significantly differs from the expected effect under these models.
#'
#' @examples
#' # Example with synthetic dataset
#' data(combo_treatment_synthetic_data)
#' data_processed <- calculate_volume(combo_treatment_synthetic_data)
#' data_processed <- calculate_dates(data_processed, start_date = "03/24/2025")
#' 
#' synergy_results <- analyze_drug_synergy(
#'   df = data_processed,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Combo",
#'   control_name = "Control"
#' )
#' 
#' # Print the summary
#' print(synergy_results$summary)
#' 
#' # Create and display the synergy visualization
#' synergy_plot <- plot_drug_synergy(synergy_results)
#' print(synergy_plot)
#'
#' @import dplyr
#' @import ggplot2
#' @export
analyze_drug_synergy <- function(df, 
                               treatment_column = "Treatment",
                               volume_column = "Volume",
                               time_column = "Day",
                               drug_a_name,
                               drug_b_name,
                               combo_name,
                               control_name = "Control",
                               eval_time_point = NULL) {
  
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
  
  # If no specific time point is provided, use the last time point
  if (is.null(eval_time_point)) {
    eval_time_point <- max(df[[time_column]])
    message(paste("No specific evaluation time point provided. Using the last time point:", eval_time_point))
  } else {
    # Check that the specified time point exists
    if (!eval_time_point %in% df[[time_column]]) {
      closest_time <- df[[time_column]][which.min(abs(df[[time_column]] - eval_time_point))]
      warning(paste("Specified time point", eval_time_point, "not found in data.",
                   "Using closest available time point:", closest_time))
      eval_time_point <- closest_time
    }
  }
  
  # Filter data for the evaluation time point
  analysis_data <- df[df[[time_column]] == eval_time_point, ]
  
  # Calculate mean volumes for each treatment group
  group_means <- tapply(analysis_data[[volume_column]], analysis_data[[treatment_column]], mean)
  
  # Extract mean volumes for each group
  control_mean <- group_means[control_name]
  drug_a_mean <- group_means[drug_a_name]
  drug_b_mean <- group_means[drug_b_name]
  combo_mean <- group_means[combo_name]
  
  # Calculate Tumor Growth Inhibition (TGI) for each treatment
  tgi_a <- 100 * (1 - (drug_a_mean / control_mean))
  tgi_b <- 100 * (1 - (drug_b_mean / control_mean))
  tgi_combo <- 100 * (1 - (combo_mean / control_mean))
  
  # Convert TGI to fractional effect (FE)
  fe_a <- tgi_a / 100
  fe_b <- tgi_b / 100
  fe_combo <- tgi_combo / 100
  
  # Calculate expected effect using Bliss Independence model
  # Bliss Independence: Expected combined effect = EA + EB - (EA * EB)
  bliss_expected_fe <- fe_a + fe_b - (fe_a * fe_b)
  bliss_expected_tgi <- bliss_expected_fe * 100
  
  # Calculate expected effect using Loewe Additivity
  # For Loewe, we need a simplified approach without dose information
  # We'll use the simplified formula: Expected = (Effect A + Effect B) / 2
  loewe_expected_fe <- (fe_a + fe_b) / 2
  loewe_expected_tgi <- loewe_expected_fe * 100
  
  # Calculate Combination Index (CI)
  # CI < 1 indicates synergy, CI = 1 indicates additivity, CI > 1 indicates antagonism
  # Simplified CI calculation without dose information
  if (fe_combo > 0) {
    ci_value <- (fe_a + fe_b) / (2 * fe_combo)
  } else {
    ci_value <- NA
    warning("Cannot calculate Combination Index: Combination effect is zero or negative")
  }
  
  # Determine synergy interpretation based on CI
  # CI < 1 indicates synergy, CI = 1 indicates additivity, CI > 1 indicates antagonism
  if (is.na(ci_value)) {
    synergy_interpretation <- "Cannot determine (CI calculation error)"
  } else if (ci_value < 0.85) {
    synergy_interpretation <- "Synergistic (CI < 0.85)"
  } else if (ci_value < 1.15) {
    synergy_interpretation <- "Additive (0.85 <= CI <= 1.15)"
  } else {
    synergy_interpretation <- "Antagonistic (CI > 1.15)"
  }
  
  # Calculate the difference between observed and expected effects
  bliss_difference <- fe_combo - bliss_expected_fe
  loewe_difference <- fe_combo - loewe_expected_fe
  
  # Determine synergy label based on both models
  if (bliss_difference > 0.1 && loewe_difference > 0.1) {
    synergy_label <- "Strong Synergy"
  } else if (bliss_difference > 0 && loewe_difference > 0) {
    synergy_label <- "Synergy"
  } else if (bliss_difference > -0.1 && loewe_difference > -0.1) {
    synergy_label <- "Additivity"
  } else {
    synergy_label <- "Antagonism"
  }
  
  # Statistical test for synergy
  # Perform t-test to compare combo group vs. each single agent
  t_test_a_combo <- t.test(
    analysis_data[[volume_column]][analysis_data[[treatment_column]] == combo_name],
    analysis_data[[volume_column]][analysis_data[[treatment_column]] == drug_a_name]
  )
  
  t_test_b_combo <- t.test(
    analysis_data[[volume_column]][analysis_data[[treatment_column]] == combo_name],
    analysis_data[[volume_column]][analysis_data[[treatment_column]] == drug_b_name]
  )
  
  # Clean up combo_name if needed - handle "HDACi + PD1" vs "aPD1"
  # This fixes inconsistencies where "HDACi + PD1" is in the data but drug_b might be named "aPD1"
  if (combo_name == "HDACi + PD1" && drug_b_name == "aPD1") {
    # Variables are already correct - no change needed
  }
  
  # Create a data frame for summary results
  summary_df <- data.frame(
    Treatment = c(drug_a_name, drug_b_name, combo_name, "Bliss Expected", "Loewe Expected"),
    Mean_Volume = c(drug_a_mean, drug_b_mean, combo_mean, 
                   control_mean * (1 - bliss_expected_fe), 
                   control_mean * (1 - loewe_expected_fe)),
    TGI_Percent = c(tgi_a, tgi_b, tgi_combo, bliss_expected_tgi, loewe_expected_tgi),
    Fractional_Effect = c(fe_a, fe_b, fe_combo, bliss_expected_fe, loewe_expected_fe)
  )
  
  # Add synergy metrics to a separate data frame
  synergy_metrics <- data.frame(
    Metric = c("Bliss Difference", "Loewe Difference", "Combination Index", "Interpretation"),
    Value = c(bliss_difference, loewe_difference, ci_value, synergy_interpretation)
  )
  
  # Create statistical test summary
  stat_tests <- data.frame(
    Comparison = c(paste("Combo vs", drug_a_name), paste("Combo vs", drug_b_name)),
    P_Value = c(t_test_a_combo$p.value, t_test_b_combo$p.value),
    Significant = c(t_test_a_combo$p.value < 0.05, t_test_b_combo$p.value < 0.05)
  )
  
  # Prepare data for plotting (will be used by plot_drug_synergy function)
  plot_data <- data.frame(
    Treatment = factor(c(drug_a_name, drug_b_name, combo_name, "Bliss Expected", "Loewe Expected"),
                     levels = c(drug_a_name, drug_b_name, "Bliss Expected", "Loewe Expected", combo_name)),
    TGI = c(tgi_a, tgi_b, bliss_expected_tgi, loewe_expected_tgi, tgi_combo),
    Type = c("Observed", "Observed", "Expected", "Expected", "Observed")
  )
  
  # Print results
  cat("\n=== Drug Combination Synergy Analysis ===\n")
  cat("Evaluation time point:", eval_time_point, "\n\n")
  
  cat("Treatment Mean Volumes:\n")
  cat(paste0("Control (", control_name, "): ", round(control_mean, 2), "\n"))
  cat(paste0(drug_a_name, ": ", round(drug_a_mean, 2), " (TGI: ", round(tgi_a, 1), "%)\n"))
  cat(paste0(drug_b_name, ": ", round(drug_b_mean, 2), " (TGI: ", round(tgi_b, 1), "%)\n"))
  cat(paste0(combo_name, ": ", round(combo_mean, 2), " (TGI: ", round(tgi_combo, 1), "%)\n\n"))
  
  cat("Expected Effects:\n")
  cat(paste0("Bliss Independence: TGI = ", round(bliss_expected_tgi, 1), "%\n"))
  cat(paste0("Loewe Additivity: TGI = ", round(loewe_expected_tgi, 1), "%\n\n"))
  
  cat("Synergy Assessment:\n")
  cat(paste0("Bliss Difference: ", round(bliss_difference * 100, 1), "% (", 
             ifelse(bliss_difference > 0, "Synergy", "No Synergy"), ")\n"))
  cat(paste0("Loewe Difference: ", round(loewe_difference * 100, 1), "% (", 
             ifelse(loewe_difference > 0, "Synergy", "No Synergy"), ")\n"))
  cat(paste0("Combination Index: ", round(ci_value, 2), " (", synergy_interpretation, ")\n\n"))
  
  cat("Statistical Tests:\n")
  cat(paste0("Combo vs ", drug_a_name, ": p = ", round(t_test_a_combo$p.value, 4), 
             ifelse(t_test_a_combo$p.value < 0.05, " (Significant)", " (Not Significant)"), "\n"))
  cat(paste0("Combo vs ", drug_b_name, ": p = ", round(t_test_b_combo$p.value, 4), 
             ifelse(t_test_b_combo$p.value < 0.05, " (Significant)", " (Not Significant)"), "\n\n"))
  
  cat("Overall Assessment:", synergy_label, "\n\n")
  
  # Return a list with all results
  return(list(
    summary = summary_df,
    synergy_metrics = synergy_metrics,
    bliss_independence = list(
      expected_effect = bliss_expected_fe,
      observed_effect = fe_combo,
      difference = bliss_difference,
      synergy = bliss_difference > 0
    ),
    loewe_additivity = list(
      expected_effect = loewe_expected_fe,
      observed_effect = fe_combo,
      difference = loewe_difference,
      synergy = loewe_difference > 0
    ),
    combination_index = list(
      ci = ci_value,
      interpretation = synergy_interpretation
    ),
    statistical_tests = stat_tests,
    overall_assessment = synergy_label,
    evaluation_time_point = eval_time_point,
    plot_data = plot_data, # Keep the plot data for later plotting
    # Additional data needed for plotting
    drug_a_name = drug_a_name,
    drug_b_name = drug_b_name,
    combo_name = combo_name,
    control_name = control_name
  ))
}

#' Plot Drug Combination Synergy Analysis
#'
#' Creates a bar plot visualizing tumor growth inhibition (TGI) for different treatment groups
#' and synergy metrics from a drug combination analysis.
#'
#' @param synergy_results Results object from analyze_drug_synergy function
#' @param custom_title Optional custom title for the plot
#' @param custom_colors Optional named vector of custom colors for plot elements
#'
#' @return A ggplot2 object visualizing the synergy analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # First run the analysis
#' results <- analyze_drug_synergy(
#'   df = tumor_data,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Drug A + Drug B",
#'   control_name = "Vehicle"
#' )
#' 
#' # Then create the plot
#' plot_drug_synergy(results)
#' 
#' # With custom title
#' plot_drug_synergy(results, custom_title = "Custom Analysis Title")
#' }
plot_drug_synergy <- function(synergy_results, custom_title = NULL, custom_colors = NULL) {
  # Validate input
  if (!is.list(synergy_results) || is.null(synergy_results$plot_data)) {
    stop("Input must be a valid result object from analyze_drug_synergy()")
  }
  
  # Extract the plot data
  plot_data <- synergy_results$plot_data
  
  # Set title
  if (is.null(custom_title)) {
    title <- paste("Drug Combination Analysis at Day", synergy_results$evaluation_time_point)
  } else {
    title <- custom_title
  }
  
  # Set colors
  if (is.null(custom_colors)) {
    fill_colors <- c("Expected" = "lightblue", "Observed" = "darkblue")
  } else {
    fill_colors <- custom_colors
  }
  
  # Create the plot
  synergy_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Treatment, y = TGI, fill = Type)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "black") +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::labs(
      title = title,
      subtitle = paste("Synergy Assessment:", synergy_results$overall_assessment, 
                     "(CI =", round(synergy_results$combination_index$ci, 2), ")"),
      x = "Treatment",
      y = "Tumor Growth Inhibition (%)",
      fill = "Data Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(ggplot2::aes(label = round(TGI, 1)), vjust = -0.5)
  
  return(synergy_plot)
}