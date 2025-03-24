#' Perform Survival Analysis and Return Hazard Ratios, Confidence Intervals, and P-Values
#'
#' This function performs survival analysis using the Cox Proportional Hazards Model.
#' It calculates the hazard ratio for different treatment groups, along with the 95% confidence intervals
#' and p-values to assess the statistical significance of the treatment on survival outcomes.
#' The function uses the `survival` package to fit the Cox model and returns a summary of the results.
#'
#' @param df A data frame containing the data to be analyzed.
#' @param time_column A character string specifying the name of the column representing the time to event (default: "Day").
#' @param censor_column A character string specifying the name of the column representing the censoring indicator (1 = event, 0 = censored) (default: "Survival_Censor").
#' @param treatment_column A character string specifying the name of the column representing the treatment variable (default: "Treatment").
#' @param cage_column A character string specifying the name of the column representing the cage identifier (default: "Cage"). The function will test for collinearity between Cage and Treatment variables. If there's no collinearity, Cage will be included as a nested variable within Treatment in the statistical model.
#' @param id_column A character string specifying the name of the column representing the individual identifier for each mouse (default: "ID").
#' @param dose_column Optional. A character string specifying the name of the column containing dose information (default: NULL).
#' @param reference_group A character string specifying the treatment group to use as the reference/comparison group (default: NULL, which uses the first level alphabetically). The reference group will have a hazard ratio of exactly 1.0, and all other groups' hazard ratios will be calculated relative to it.
#'
#' @return A list containing:
#' \item{model}{The fitted Cox Proportional Hazards model object.}
#' \item{results}{A data frame containing the following columns:
#' \itemize{
#'   \item{Group}{The treatment or group category.}
#'   \item{Hazard_Ratio}{The hazard ratio (exp(coef)) for the group.}
#'   \item{CI_Lower}{The lower bound of the 95% confidence interval for the hazard ratio.}
#'   \item{CI_Upper}{The upper bound of the 95% confidence interval for the hazard ratio.}
#'   \item{P_Value}{The p-value associated with the group effect in the Cox model.}
#' }}
#' \item{forest_plot}{A ggplot2 object visualizing the hazard ratios with confidence intervals as a forest plot.
#' The plot displays each group's hazard ratio on a logarithmic scale, with confidence intervals and exact values.}
#'
#' @details The function fits a Cox Proportional Hazards model using the `survival` package,
#' and the output includes hazard ratios, confidence intervals, and p-values for the specified treatment groups.
#' The function assumes that the time-to-event data and censoring indicators are appropriately represented
#' in the input data frame. It ensures each mouse is counted only once by using the last observation for each subject.
#' 
#' The reference group is included in the results with a hazard ratio of 1.0, since all comparisons are made
#' relative to this group. In the forest plot, the reference group is clearly labeled as such.
#' 
#' The function handles cage effects and tests for collinearity between Cage and Treatment variables:
#' \itemize{
#'   \item The function prints a contingency table showing the distribution of cages across treatment groups.
#'   \item If Cage and Treatment are collinear (e.g., each cage contains only one treatment group), 
#'         only the Treatment effect is included in the model.
#'   \item If there's no collinearity, Cage is included as a clustering variable using `cluster(cage)`. 
#'         This accounts for correlation within cages by providing robust standard errors, while still
#'         estimating treatment effects.
#'   \item The function creates a nested cage variable (`cage_in_group`) that explicitly represents
#'         the hierarchical structure of cages within treatment groups.
#'   \item The model formula used is reported as a message when the function runs.
#' }
#'
#' @examples
#' # Example usage with default reference group (first alphabetically)
#' results <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", 
#'                               treatment_column = "Treatment", cage_column = "Cage", id_column = "ID")
#' print(results$results)
#' 
#' # Example usage with specified reference group ("Control")
#' results_control_ref <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", 
#'                                          treatment_column = "Treatment", cage_column = "Cage", id_column = "ID", 
#'                                          reference_group = "Control")
#' print(results_control_ref$results)
#'
#' # Example usage with "Drug A" as reference group
#' results_drug_a_ref <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", 
#'                                         treatment_column = "Treatment", cage_column = "Cage", id_column = "ID", 
#'                                         reference_group = "Drug A")
#' print(results_drug_a_ref$results)
#'
#' # Example usage with dose information
#' results_with_dose <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", 
#'                                         treatment_column = "Treatment", cage_column = "Cage", id_column = "ID", 
#'                                         dose_column = "Dose")
#' print(results_with_dose$results)
#'
#' # Example with dose information and specific reference group
#' # Note: When using dose, the reference group must be the full string with dose, e.g., "Drug A - Dose: 10"
#' results_dose_ref <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", 
#'                                        treatment_column = "Treatment", cage_column = "Cage", id_column = "ID", 
#'                                        dose_column = "Dose", reference_group = "Control - Dose: 0")
#' print(results_dose_ref$results)
#'
#' @import survival
#' @import ggplot2
#'
#' @export
survival_statistics <- function(df, time_column = "Day", censor_column = "Survival_Censor", 
                          treatment_column = "Treatment", cage_column = "Cage", id_column = "ID", 
                          dose_column = NULL, reference_group = NULL) {

  # Ensure required columns exist in the dataframe
  required_columns <- c(time_column, censor_column, treatment_column, cage_column, id_column)
  missing_cols <- required_columns[!required_columns %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for dose column if specified
  if (!is.null(dose_column) && !(dose_column %in% colnames(df))) {
    warning(paste("Dose column", dose_column, "not found in data frame, proceeding without dose information"))
    dose_column <- NULL
  }
  
  # Create a composite group identifier based on Treatment (and Dose if available)
  if (!is.null(dose_column)) {
    # Create a group identifier combining Treatment and Dose
    df$group <- paste(df[[treatment_column]], df[[dose_column]], sep = " - Dose: ")
  } else {
    # Use Treatment as the group identifier
    df$group <- df[[treatment_column]]
  }
  
  # Create a unique subject identifier to ensure each mouse is only counted once
  df$subject_id <- paste(df[[cage_column]], df[[id_column]], sep = "_")
  
  # Create a unique identifier for each mouse (combining cage, treatment, ID, and optionally dose)
  if (!is.null(dose_column)) {
    df$Mouse_ID <- paste(df[[cage_column]], df[[treatment_column]], 
                        df[[dose_column]], df[[id_column]], sep = "_")
  } else {
    df$Mouse_ID <- paste(df[[cage_column]], df[[treatment_column]], 
                        df[[id_column]], sep = "_")
  }
  df$Mouse_ID <- factor(df$Mouse_ID)
  
  # Aggregate data to subject level - keep only the last time point per subject
  # First sort by subject_id and time to ensure we get the latest record per subject
  df <- df[order(df$subject_id, df[[time_column]]), ]
  
  # Get the last entry for each subject
  subjects <- unique(df$subject_id)
  last_records <- list()
  
  for (subject in subjects) {
    subject_data <- df[df$subject_id == subject, ]
    last_records[[length(last_records) + 1]] <- subject_data[nrow(subject_data), ]
  }
  
  # Combine into a single dataframe with one row per subject
  df_aggregated <- do.call(rbind, last_records)
  
  # Create a new data frame with standardized column names for analysis
  analysis_df <- data.frame(
    time = df_aggregated[[time_column]],
    status = df_aggregated[[censor_column]],
    group = df_aggregated$group,
    Mouse_ID = df_aggregated$Mouse_ID
  )
  
  # Set the reference group if specified
  if (!is.null(reference_group)) {
    # Check if the reference group exists
    if (!reference_group %in% unique(analysis_df$group)) {
      stop(paste("Reference group", reference_group, "not found in the", treatment_column, "column."))
    }
    
    # Convert group to factor with specified reference level
    analysis_df$group <- stats::relevel(factor(analysis_df$group), ref = reference_group)
    
    # Print message about reference group
    message(paste("Using", reference_group, "as the reference group for hazard ratios."))
  } else {
    # Convert to factor without changing the reference level (will use alphabetical first)
    analysis_df$group <- factor(analysis_df$group)
    
    # Print message about which group is being used as reference
    ref_group <- levels(analysis_df$group)[1]
    message(paste("Using", ref_group, "as the reference group for hazard ratios (first alphabetically)."))
  }
  
  # Check for collinearity between Cage and Treatment
  # Create a temporary variable for the cage
  analysis_df$cage <- df_aggregated[[cage_column]]
  
  # Test for collinearity between cage and treatment group
  has_collinearity <- FALSE
  
  # Only test for collinearity if we have multiple cages
  if (length(unique(analysis_df$cage)) > 1) {
    # Create a simple contingency table of cage vs treatment
    cont_table <- table(df_aggregated[[cage_column]], df_aggregated[[treatment_column]])
    
    # Check if each cage only has one treatment (perfect collinearity)
    cage_has_one_treatment <- all(rowSums(cont_table > 0) == 1)
    
    # Check if each treatment only has one cage (perfect collinearity)
    treatment_has_one_cage <- all(colSums(cont_table > 0) == 1)
    
    # If either condition is true, there's perfect collinearity
    has_collinearity <- cage_has_one_treatment || treatment_has_one_cage
    
    # Log the result
    if (has_collinearity) {
      message("Detected collinearity between Cage and Treatment. Using Treatment only in the model.")
    } else {
      message("No perfect collinearity detected between Cage and Treatment. Including Cage as a nested variable.")
    }
  } else {
    message("Only one cage detected. Using Treatment only in the model.")
    has_collinearity <- TRUE  # Treat single cage as collinear case
  }
  
  # First, check the cage distribution across treatment groups
  analysis_df$cage <- factor(analysis_df$cage)
  analysis_df$group <- factor(analysis_df$group)
  
  # Create contingency table of cages by treatment groups
  cage_group_table <- table(analysis_df$cage, analysis_df$group)
  print("Cage distribution across treatment groups:")
  print(cage_group_table)
  
  # Fit the Cox Proportional Hazards Model with proper formula
  # Either include cage as nested variable, or exclude it if collinearity exists
  if (has_collinearity) {
    # Use simple model with just treatment groups
    cox_model <- survival::coxph(survival::Surv(time, status) ~ group, data = analysis_df)
    message("Using model without cage effect: Surv(time, status) ~ group")
  } else {
    # Check if we have multiple cages per treatment
    cages_per_group <- colSums(cage_group_table > 0)
    
    # Create a properly defined nested cage variable
    analysis_df$cage_in_group <- interaction(analysis_df$group, analysis_df$cage, sep = ":")
    
    # Approach 1: Include cage as a clustering variable (robust standard errors)
    # This accounts for correlations within cages
    cox_model <- survival::coxph(survival::Surv(time, status) ~ group + survival::cluster(cage), data = analysis_df)
    message("Using model with cage clustering: Surv(time, status) ~ group + cluster(cage)")
    
    # Approach 2 (alternative if needed): Use frailty term for cage
    # Uncomment the line below to use frailty instead of clustering
    # cox_model <- survival::coxph(survival::Surv(time, status) ~ group + survival::frailty(cage), data = analysis_df)
    # message("Using model with cage frailty: Surv(time, status) ~ group + frailty(cage)")
  }

  # Summary of the Cox model
  cox_summary = summary(cox_model)
  print(cox_summary)

  # Extract model results - need to handle the coefficients and their confidence intervals
  coef_table <- cox_summary$coefficients
  conf_table <- cox_summary$conf.int
  
  # Get variable names
  var_names <- rownames(coef_table)
  
  # Extract the results for each group
  results_list <- list()
  
  # Handle columns with different possible names
  p_value_col <- ifelse("Pr(>|z|)" %in% colnames(coef_table), "Pr(>|z|)", 
                 ifelse("p" %in% colnames(coef_table), "p", NA))
  
  # Get all groups from the data and determine the reference group
  all_groups <- levels(analysis_df$group)
  ref_group <- all_groups[1]  # The reference group will be the first level in the factor
  
  # First, add the reference group (which isn't in the model output)
  row_data <- list(
    Group = ref_group,
    Hazard_Ratio = 1.0,  # Reference group has HR = 1 by definition
    CI_Lower = 1.0,      # Reference group's CI is exactly 1.0
    CI_Upper = 1.0,      # Reference group's CI is exactly 1.0
    P_Value = NA         # P-value doesn't apply to reference group
  )
  results_list[[1]] <- as.data.frame(row_data)
  
  # Then add all other groups from the model output
  for (var in var_names) {
    # Clean up the group name - extract just the actual group value
    clean_name <- gsub("^group", "", var)
    
    # Get the index for this variable in both tables
    idx <- which(rownames(coef_table) == var)
    conf_idx <- which(rownames(conf_table) == var)
    
    # Create a row for this group with safe extraction
    row_data <- list(
      Group = clean_name,
      Hazard_Ratio = if(length(conf_idx) > 0) conf_table[conf_idx, "exp(coef)"] else NA,
      CI_Lower = if(length(conf_idx) > 0 && "lower .95" %in% colnames(conf_table)) conf_table[conf_idx, "lower .95"] else NA,
      CI_Upper = if(length(conf_idx) > 0 && "upper .95" %in% colnames(conf_table)) conf_table[conf_idx, "upper .95"] else NA,
      P_Value = if(length(idx) > 0 && !is.na(p_value_col)) coef_table[idx, p_value_col] else NA
    )
    
    results_list[[length(results_list) + 1]] <- as.data.frame(row_data)
  }
  
  # We've removed the frailty term, so this is no longer needed
  
  # Combine all results
  results_df <- do.call(rbind, results_list)
  
  # Create a forest plot for hazard ratios
  if (nrow(results_df) > 0 && !all(is.na(results_df$Hazard_Ratio))) {
    # For plotting, we need to handle the reference group differently
    plot_data <- results_df
    
    # For the reference group, we don't show p-value since it's NA
    plot_data$label <- ifelse(plot_data$Group == ref_group,
                             sprintf("%.2f (Reference)", plot_data$Hazard_Ratio),
                             ifelse(is.na(plot_data$P_Value),
                                   sprintf("%.2f (%.2f-%.2f)", 
                                          plot_data$Hazard_Ratio, 
                                          plot_data$CI_Lower, 
                                          plot_data$CI_Upper),
                                   sprintf("%.2f (%.2f-%.2f), p=%.3f", 
                                          plot_data$Hazard_Ratio, 
                                          plot_data$CI_Lower, 
                                          plot_data$CI_Upper, 
                                          plot_data$P_Value)))
    
    if (nrow(plot_data) > 0) {
      # Create the forest plot using ggplot2
      forest_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Hazard_Ratio, y = Group)) +
        # Reference line at HR = 1
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
        # Confidence intervals
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
        # Point estimates
        ggplot2::geom_point(size = 3, color = "blue") +
        # Labels
        ggplot2::labs(
          title = "Forest Plot of Hazard Ratios",
          subtitle = "Values < 1 indicate reduced hazard; values > 1 indicate increased hazard",
          x = "Hazard Ratio (95% CI)",
          y = NULL
        ) +
        # Log scale for better visualization of hazard ratios
        ggplot2::scale_x_log10() +
        # Add values as text
        ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.1) +
        # Theme adjustments
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank()
        )
      
      # Print the forest plot
      print(forest_plot)
    } else {
      warning("Cannot create forest plot: No valid hazard ratio data with confidence intervals")
      forest_plot <- NULL
    }
  } else {
    warning("Cannot create forest plot: No hazard ratio data available")
    forest_plot <- NULL
  }

  # Return the results as a list
  return(list(
    model = cox_model,
    results = results_df,
    forest_plot = forest_plot
  ))
}
