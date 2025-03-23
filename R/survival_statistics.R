#' Perform Survival Analysis and Return Hazard Ratios, Confidence Intervals, and P-Values
#'
#' This function performs survival analysis using the Cox Proportional Hazards Model.
#' It calculates the hazard ratio for different groups, along with the 95% confidence intervals
#' and p-values to assess the statistical significance of the group variable on survival outcomes.
#' The function uses the `survival` package to fit the Cox model and returns a summary of the results.
#'
#' @param df A data frame containing the data to be analyzed.
#' @param time_column A character string specifying the name of the column representing the time to event (default: "Day").
#' @param censor_column A character string specifying the name of the column representing the censoring indicator (1 = event, 0 = censored) (default: "Survival_Censor").
#' @param group_column A character string specifying the name of the column representing the treatment or group variable (default: "Group").
#' @param id_column A character string specifying the name of the column representing the individual identifier for each mouse (default: "ID").
#' @param reference_group A character string specifying the level to use as the reference/comparison group (default: NULL, which uses the first level alphabetically).
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
#' in the input data frame, and the model accounts for individual mice as random effects.
#'
#' @examples
#' # Example usage with default reference group (first alphabetically)
#' results <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", 
#'                               group_column = "Group", id_column = "ID")
#' print(results$results)
#' 
#' # Example usage with specified reference group
#' results_b_ref <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", 
#'                                     group_column = "Group", id_column = "ID", reference_group = "B")
#' print(results_b_ref$results)
#'
#' @import survival
#' @import ggplot2
#'
#' @export
survival_statistics <- function(df, time_column = "Day", censor_column = "Survival_Censor", 
                          group_column = "Group", id_column = "ID", reference_group = NULL) {

  # Ensure required columns exist in the dataframe
  if (!all(c(time_column, censor_column, group_column, id_column) %in% colnames(df))) {
    stop("One or more specified columns are not found in the dataframe.")
  }

  # Create a unique identifier for each mouse (based on ID and Group)
  df$Mouse_ID <- factor(paste(df[[group_column]], df[[id_column]], sep = "_"))

  # Create a new data frame with standardized column names for analysis
  analysis_df <- data.frame(
    time = df[[time_column]],
    status = df[[censor_column]],
    group = df[[group_column]],
    Mouse_ID = df$Mouse_ID
  )
  
  # Set the reference group if specified
  if (!is.null(reference_group)) {
    # Check if the reference group exists
    if (!reference_group %in% unique(analysis_df$group)) {
      stop(paste("Reference group", reference_group, "not found in the", group_column, "column."))
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
  
  # Fit the Cox Proportional Hazards Model with proper formula
  cox_model = survival::coxph(survival::Surv(time, status) ~ group + survival::frailty(Mouse_ID), data = analysis_df)

  # Summary of the Cox model
  cox_summary = summary(cox_model)
  print(cox_summary)

  # Extract model results - need to handle the coefficients and their confidence intervals
  coef_table <- cox_summary$coefficients
  conf_table <- cox_summary$conf.int
  
  # Get variable names - only from the fixed effects (exclude frailty term)
  var_names <- rownames(coef_table)
  var_names <- var_names[!grepl("frailty", var_names)]
  
  # Extract the results for each group
  results_list <- list()
  
  # Handle columns with different possible names
  p_value_col <- ifelse("Pr(>|z|)" %in% colnames(coef_table), "Pr(>|z|)", 
                 ifelse("p" %in% colnames(coef_table), "p", NA))
  
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
  
  # Add frailty term if it exists
  frailty_idx <- grep("frailty", rownames(coef_table))
  if (length(frailty_idx) > 0) {
    row_data <- list(
      Group = "Mouse ID (random effect)",
      Hazard_Ratio = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      P_Value = if(!is.na(p_value_col)) coef_table[frailty_idx, p_value_col] else NA
    )
    results_list[[length(results_list) + 1]] <- as.data.frame(row_data)
  }
  
  # Combine all results
  results_df <- do.call(rbind, results_list)
  
  # Create a forest plot for hazard ratios
  if (nrow(results_df) > 0 && !all(is.na(results_df$Hazard_Ratio))) {
    # Filter out rows with NA values for plotting
    plot_data <- results_df[!is.na(results_df$Hazard_Ratio) & 
                           !is.na(results_df$CI_Lower) &
                           !is.na(results_df$CI_Upper), ]
    
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
        ggplot2::geom_text(ggplot2::aes(
          label = sprintf("%.2f (%.2f-%.2f), p=%.3f", 
                         Hazard_Ratio, CI_Lower, CI_Upper, P_Value)
        ), hjust = -0.1) +
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
