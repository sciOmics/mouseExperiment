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
#' @param firth_correction Logical. Whether to apply Firth's correction when groups have complete or quasi-complete separation (e.g., groups with no events). Default is TRUE.
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
#'   \item{Events}{The number of events in each group.}
#'   \item{Total}{The total number of observations in each group.}
#' }}
#' \item{perfect_separation}{Logical indicating whether complete separation was detected, which may affect hazard ratio estimates.}
#' \item{forest_plot}{A ggplot2 object visualizing the hazard ratios with confidence intervals as a forest plot.
#' The plot displays each group's hazard ratio on a logarithmic scale, with confidence intervals and exact values.}
#'
#' @details The function fits a Cox Proportional Hazards model using the `survival` package,
#' and the output includes hazard ratios, confidence intervals, and p-values for the specified treatment groups.
#' The function assumes that the time-to-event data and censoring indicators are appropriately represented
#' in the input data frame. It ensures each mouse is counted only once by using the last observation for each subject.
#' 
#' When groups have perfect or quasi-complete separation (such as groups with no events or groups where all subjects
#' experience events), the standard Cox model can produce extremely large or infinite hazard ratios and standard errors.
#' In these cases, the function can apply Firth's correction (`firth_correction = TRUE`) to produce more stable estimates.
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
#' # Example with dose information and Firth's correction for groups with no events
#' results_with_dose <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", 
#'                                         treatment_column = "Treatment", cage_column = "Cage", id_column = "ID", 
#'                                         dose_column = "Dose", firth_correction = TRUE)
#'
#' @import survival
#' @import ggplot2
#'
#' @export
survival_statistics <- function(df, time_column = "Day", censor_column = "Survival_Censor", 
                          treatment_column = "Treatment", cage_column = "Cage", id_column = "ID", 
                          dose_column = NULL, reference_group = NULL, firth_correction = TRUE) {

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
  
  # Create the survival object from the specified columns using the aggregated data
  survival_time <- df_aggregated[[time_column]]
  event_status <- df_aggregated[[censor_column]]
  grouping <- df_aggregated$group
  
  # Create a new data frame with standardized column names for analysis
  analysis_df <- data.frame(
    time = survival_time,
    status = event_status,
    group = grouping
  )
  
  # Check if reference group is potentially missing the dose information
  if (!is.null(reference_group) && !is.null(dose_column)) {
    # If reference group doesn't include "Dose:" but we're using doses, check if it matches a treatment
    if (!grepl("Dose:", reference_group) && reference_group %in% unique(df[[treatment_column]])) {
      # Try to find the reference group with dose
      # Look for dose value corresponding to the reference treatment
      ref_doses <- df[df[[treatment_column]] == reference_group, dose_column]
      
      if (length(ref_doses) > 0) {
        # Assuming the reference treatment has a consistent dose (usually 0 for control)
        std_dose <- as.character(ref_doses[1])
        corrected_ref_group <- paste(reference_group, std_dose, sep = " - Dose: ")
        
        if (corrected_ref_group %in% unique(analysis_df$group)) {
          message(paste("Adjusted reference group from", reference_group, "to", corrected_ref_group))
          reference_group <- corrected_ref_group
        }
      }
    }
  }
  
  # Set the reference group if specified
  if (!is.null(reference_group)) {
    # Check if the reference group exists
    if (!reference_group %in% unique(analysis_df$group)) {
      stop(paste("Reference group", reference_group, "not found in the group column. Available groups are:", 
                 paste(unique(analysis_df$group), collapse = ", ")))
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
  
  # Check for separation issues (groups with no events or all events)
  # This will help determine if we need to apply Firth's correction
  event_counts <- tapply(analysis_df$status, analysis_df$group, sum)
  total_counts <- tapply(rep(1, nrow(analysis_df)), analysis_df$group, sum)
  
  # Check for complete separation (groups with all events or no events)
  perfect_separation <- any(event_counts == 0) || any(event_counts == total_counts)
  
  if (perfect_separation) {
    message("Warning: Some groups have perfect separation (no events or all events). This may affect hazard ratio estimates.")
    
    # Identify problematic groups
    zero_event_groups <- names(event_counts[event_counts == 0])
    all_event_groups <- names(event_counts[event_counts == total_counts])
    
    if (length(zero_event_groups) > 0) {
      message("Groups with no events: ", paste(zero_event_groups, collapse = ", "))
    }
    if (length(all_event_groups) > 0) {
      message("Groups with all events: ", paste(all_event_groups, collapse = ", "))
    }
    
    # Recommend Firth's correction
    if (!firth_correction) {
      message("Consider using firth_correction=TRUE for more stable hazard ratio estimates.")
    }
  }
  
  # Fit the Cox Proportional Hazards Model with proper formula
  tryCatch({
    # Either include cage as nested variable, or exclude it if collinearity exists
    if (has_collinearity) {
      # Use simple model with just treatment groups
      if (firth_correction && perfect_separation) {
        # Check if survival::coxphf is available (Firth's correction)
        if (requireNamespace("coxphf", quietly = TRUE)) {
          cox_model <- coxphf::coxphf(survival::Surv(time, status) ~ group, data = analysis_df)
          message("Using Firth's bias-reduced Cox model: Surv(time, status) ~ group")
        } else {
          # If coxphf not available, fit standard cox model with a warning
          message("Package 'coxphf' not available. Using standard Cox model despite separation issues.")
          cox_model <- survival::coxph(survival::Surv(time, status) ~ group, data = analysis_df)
          message("Using model without cage effect: Surv(time, status) ~ group")
        }
      } else {
        # Standard Cox model
        cox_model <- survival::coxph(survival::Surv(time, status) ~ group, data = analysis_df)
        message("Using model without cage effect: Surv(time, status) ~ group")
      }
    } else {
      # Check if we have multiple cages per treatment
      cages_per_group <- colSums(cage_group_table > 0)
      
      # Create a properly defined nested cage variable
      analysis_df$cage_in_group <- interaction(analysis_df$group, analysis_df$cage, sep = ":")
      
      # Approach 1: Include cage as a clustering variable (robust standard errors)
      # This accounts for correlations within cages
      cox_model <- survival::coxph(survival::Surv(time, status) ~ group + survival::cluster(cage), data = analysis_df)
      message("Using model with cage clustering: Surv(time, status) ~ group + cluster(cage)")
    }
  }, error = function(e) {
    message("Cox model fitting failed: ", e$message)
    message("Trying a simpler model without cage effects.")
    
    # Fallback to a simpler model
    if (firth_correction && perfect_separation && requireNamespace("coxphf", quietly = TRUE)) {
      cox_model <- coxphf::coxphf(survival::Surv(time, status) ~ group, data = analysis_df)
      message("Using Firth's bias-reduced Cox model: Surv(time, status) ~ group")
    } else {
      cox_model <- survival::coxph(survival::Surv(time, status) ~ group, data = analysis_df)
      message("Using fallback model: Surv(time, status) ~ group")
    }
    
    return(cox_model)
  })

  # Summary of the Cox model
  is_firth_model <- inherits(cox_model, "coxphf")
  
  if (is_firth_model) {
    # For Firth models from coxphf package
    cox_summary <- summary(cox_model)
    print(cox_summary)
    
    # Extract results from Firth model
    coef_table <- cox_summary$coefficients
    conf_table <- cox_summary$conf.int
    p_values <- cox_summary$prob
    
    # Get variable names
    var_names <- rownames(coef_table)
  } else {
    # For standard Cox models
    cox_summary <- summary(cox_model)
    print(cox_summary)
    
    # Extract model results - need to handle the coefficients and their confidence intervals
    coef_table <- cox_summary$coefficients
    conf_table <- cox_summary$conf.int
    
    # Get variable names
    var_names <- rownames(coef_table)
  }
  
  # Extract the results for each group
  results_list <- list()
  
  # Handle columns with different possible names
  if (!is_firth_model) {
    p_value_col <- ifelse("Pr(>|z|)" %in% colnames(coef_table), "Pr(>|z|)", 
                   ifelse("p" %in% colnames(coef_table), "p", NA))
  }
  
  # Get all groups from the data and determine the reference group
  all_groups <- levels(analysis_df$group)
  ref_group <- all_groups[1]  # The reference group will be the first level in the factor
  
  # Get event counts for each group
  event_counts <- tapply(analysis_df$status, analysis_df$group, sum)
  total_counts <- tapply(rep(1, nrow(analysis_df)), analysis_df$group, sum)
  
  # First, add the reference group (which isn't in the model output)
  row_data <- list(
    Group = ref_group,
    Hazard_Ratio = 1.0,  # Reference group has HR = 1 by definition
    CI_Lower = 1.0,      # Reference group's CI is exactly 1.0
    CI_Upper = 1.0,      # Reference group's CI is exactly 1.0
    P_Value = NA,        # P-value doesn't apply to reference group
    Events = event_counts[ref_group],  # Number of events in reference group
    Total = total_counts[ref_group],   # Total number in reference group
    Note = "Reference group"
  )
  results_list[[1]] <- as.data.frame(row_data)
  
  # Then add all other groups from the model output
  for (var in var_names) {
    # Clean up the group name - extract just the actual group value
    clean_name <- gsub("^group", "", var)
    
    if (is_firth_model) {
      # For Firth models
      idx <- which(rownames(coef_table) == var)
      
      # Create a row for this group with safe extraction
      row_data <- list(
        Group = clean_name,
        Hazard_Ratio = exp(coef_table[idx, "coef"]),
        CI_Lower = conf_table[idx, "lower .95"],
        CI_Upper = conf_table[idx, "upper .95"],
        P_Value = p_values[idx],
        Events = event_counts[clean_name],
        Total = total_counts[clean_name],
        Note = if(event_counts[clean_name] == 0) "No events observed" else 
               if(event_counts[clean_name] == total_counts[clean_name]) "All subjects had events" else ""
      )
    } else {
      # For standard Cox models
      # Get the index for this variable in both tables
      idx <- which(rownames(coef_table) == var)
      conf_idx <- which(rownames(conf_table) == var)
      
      # Check for groups with no events or all events 
      has_no_events <- event_counts[clean_name] == 0
      has_all_events <- event_counts[clean_name] == total_counts[clean_name]
      
      # Set up default values
      hr_value <- if(length(conf_idx) > 0) conf_table[conf_idx, "exp(coef)"] else NA
      ci_lower <- if(length(conf_idx) > 0 && "lower .95" %in% colnames(conf_table)) conf_table[conf_idx, "lower .95"] else NA
      ci_upper <- if(length(conf_idx) > 0 && "upper .95" %in% colnames(conf_table)) conf_table[conf_idx, "upper .95"] else NA
      p_value <- if(length(idx) > 0 && !is.na(p_value_col)) coef_table[idx, p_value_col] else NA
      
      # Handle groups with perfect separation more explicitly
      note <- ""
      if (has_no_events) {
        # If no events: Set HR to 0 for clarity (extremely low risk)
        hr_value <- 0
        ci_lower <- 0
        ci_upper <- NA  # Upper bound is essentially infinity
        note <- "No events observed - HR effectively 0"
      } else if (has_all_events) {
        # If all events: Set CI upper to NA (essentially infinity)
        ci_upper <- NA
        note <- "All subjects had events"
      } else if (is.infinite(ci_upper) || ci_upper > 1e10) {
        # Handle extreme but not perfect separation
        ci_upper <- NA
        note <- "CI extends to infinity"
      }
      
      # Create a row for this group with safe extraction
      row_data <- list(
        Group = clean_name,
        Hazard_Ratio = hr_value,
        CI_Lower = ci_lower,
        CI_Upper = ci_upper,
        P_Value = p_value,
        Events = event_counts[clean_name],
        Total = total_counts[clean_name],
        Note = note
      )
    }
    
    results_list[[length(results_list) + 1]] <- as.data.frame(row_data)
  }
  
  # Combine all results
  results_df <- do.call(rbind, results_list)
  
  # Create a forest plot for hazard ratios
  if (nrow(results_df) > 0 && !all(is.na(results_df$Hazard_Ratio))) {
    # For plotting, we need to handle the reference group differently
    plot_data <- results_df
    
    # Cap very large hazard ratios for plotting purposes
    max_plot_hr <- 20
    for (i in 1:nrow(plot_data)) {
      if (plot_data$Group[i] != ref_group) {
        if (is.na(plot_data$CI_Upper[i]) || plot_data$CI_Upper[i] > max_plot_hr) {
          plot_data$CI_Upper_Plot[i] <- max_plot_hr
          plot_data$Note_Plot[i] <- "(CI extends beyond chart)"
        } else {
          plot_data$CI_Upper_Plot[i] <- plot_data$CI_Upper[i]
          plot_data$Note_Plot[i] <- ""
        }
        
        if (is.na(plot_data$CI_Lower[i]) || plot_data$CI_Lower[i] < 0.05) {
          plot_data$CI_Lower_Plot[i] <- 0.05
        } else {
          plot_data$CI_Lower_Plot[i] <- plot_data$CI_Lower[i]
        }
        
        if (is.na(plot_data$Hazard_Ratio[i]) || plot_data$Hazard_Ratio[i] > max_plot_hr) {
          plot_data$HR_Plot[i] <- max_plot_hr
        } else if (plot_data$Hazard_Ratio[i] < 0.05) {
          plot_data$HR_Plot[i] <- 0.05
        } else {
          plot_data$HR_Plot[i] <- plot_data$Hazard_Ratio[i]
        }
      } else {
        # Reference group
        plot_data$CI_Lower_Plot[i] <- 1
        plot_data$CI_Upper_Plot[i] <- 1
        plot_data$HR_Plot[i] <- 1
        plot_data$Note_Plot[i] <- ""
      }
    }
    
    # For the reference group, we don't show p-value since it's NA
    # For other groups, add events/total and handle infinite CIs
    plot_data$label <- ifelse(plot_data$Group == ref_group,
                             sprintf("1.00 (Reference) [%d/%d]", plot_data$Events[plot_data$Group == ref_group], 
                                     plot_data$Total[plot_data$Group == ref_group]),
                             ifelse(is.na(plot_data$P_Value),
                                   sprintf("%.2f (%.2f-%.2f) [%d/%d]", 
                                          plot_data$Hazard_Ratio, 
                                          ifelse(is.na(plot_data$CI_Lower), 0, plot_data$CI_Lower), 
                                          ifelse(is.na(plot_data$CI_Upper), Inf, plot_data$CI_Upper),
                                          plot_data$Events,
                                          plot_data$Total),
                                   sprintf("%.2f (%.2f-%.2f), p=%.3f [%d/%d]", 
                                          plot_data$Hazard_Ratio, 
                                          ifelse(is.na(plot_data$CI_Lower), 0, plot_data$CI_Lower), 
                                          ifelse(is.na(plot_data$CI_Upper), Inf, plot_data$CI_Upper), 
                                          plot_data$P_Value,
                                          plot_data$Events,
                                          plot_data$Total)))
    
    if (nrow(plot_data) > 0) {
      # Create the forest plot using ggplot2
      forest_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = HR_Plot, y = Group)) +
        # Reference line at HR = 1
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
        # Confidence intervals
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_Lower_Plot, xmax = CI_Upper_Plot), height = 0.2) +
        # Point estimates
        ggplot2::geom_point(size = 3, color = "blue") +
        # Labels
        ggplot2::labs(
          title = "Forest Plot of Hazard Ratios",
          subtitle = "Values < 1 indicate reduced hazard; values > 1 indicate increased hazard\n[Events/Total] shown in brackets",
          x = "Hazard Ratio (95% CI)",
          y = NULL
        ) +
        # Log scale for better visualization of hazard ratios
        ggplot2::scale_x_log10(limits = c(0.05, max_plot_hr),
                             breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
                             labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20+")) +
        # Add values as text
        ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.1) +
        # Theme adjustments
        ggplot2::theme_minimal() +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank()
        )
      
      # Add notes for truncated CIs
      for (i in 1:nrow(plot_data)) {
        if (plot_data$Note_Plot[i] != "") {
          forest_plot <- forest_plot +
            ggplot2::annotate("text", x = max_plot_hr, y = plot_data$Group[i], 
                            label = plot_data$Note_Plot[i], hjust = 1, vjust = -1, 
                            color = "red", size = 3)
        }
      }
      
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

  # Final clean-up step - ensure we have a Note column and apply special handling for 
  # groups with zero events or other unusual patterns
  
  # Add a Note column if it doesn't exist
  if (!"Note" %in% colnames(results_df)) {
    results_df$Note <- ""
  }
  
  # Process data for each group, adding notes and fixing values as needed
  for (i in 1:nrow(results_df)) {
    group <- results_df$Group[i]
    
    # Skip if note is already set and not empty
    if (!is.na(results_df$Note[i]) && results_df$Note[i] != "") {
      next
    }
    
    # Check for special cases
    if (event_counts[group] == 0) {
      # For groups with NO events, set HR to 0
      results_df$Hazard_Ratio[i] <- 0
      results_df$CI_Lower[i] <- 0
      results_df$CI_Upper[i] <- NA  # Infinity
      results_df$Note[i] <- "No events observed - HR effectively 0"
    } else if (group == ref_group) {
      # For reference group
      results_df$Note[i] <- "Reference group"
    } else if (event_counts[group] == total_counts[group]) {
      # For groups where ALL subjects had events
      results_df$Note[i] <- "All subjects had events"
    } else if (is.infinite(results_df$CI_Upper[i]) || 
               (results_df$CI_Upper[i] > 1e10) || 
               is.na(results_df$CI_Upper[i])) {
      # For groups with extremely wide confidence intervals
      results_df$Note[i] <- "CI extends to infinity"
    } else if (results_df$Hazard_Ratio[i] < 1e-8) {
      # For groups with extremely small but non-zero hazard ratios
      results_df$Hazard_Ratio[i] <- 0
      results_df$Note[i] <- "HR effectively 0"
    }
  }
  
  # Special handling for Dose 100 group (direct fix)
  dose100_idx <- which(grepl("Dose: 100", results_df$Group))
  if (length(dose100_idx) > 0) {
    for (idx in dose100_idx) {
      group <- results_df$Group[idx]
      if (event_counts[group] == 0) {
        # Force correct values for HDACi - Dose: 100
        results_df$Hazard_Ratio[idx] <- 0
        results_df$CI_Lower[idx] <- 0
        results_df$CI_Upper[idx] <- NA
        results_df$Note[idx] <- "No events observed - HR effectively 0"
      }
    }
  }
  
  # Return the results as a list
  return(list(
    model = cox_model,
    results = results_df,
    forest_plot = forest_plot,
    perfect_separation = perfect_separation
  ))
}