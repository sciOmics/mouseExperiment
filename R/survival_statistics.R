# Copyright (c) 2025 Insight BioAnalytics. All rights reserved.
# Proprietary and confidential.

#' Perform Survival Analysis for Mouse Tumor Experiments
#'
#' @description
#' Performs comprehensive survival analysis using appropriate statistical methods
#' based on data characteristics. Automatically selects between Cox Proportional 
#' Hazards Model, Firth's bias-reduced estimation, or Log-Rank Test depending on
#' the presence of complete or quasi-complete separation in the data.
#'
#' @param df Data frame containing survival data
#' @param time_column Name of column containing time-to-event data. Default: "Day"
#' @param censor_column Name of column containing censoring indicator (1=event, 0=censored). Default: "Survival_Censor"
#' @param treatment_column Name of column containing treatment groups. Default: "Treatment"
#' @param cage_column Name of column containing cage identifiers. Default: "Cage"
#' @param id_column Name of column containing individual mouse identifiers. Default: "ID"
#' @param dose_column Optional name of column containing dose information. Default: NULL
#' @param reference_group Treatment group to use as reference. Default: NULL (uses first alphabetically)
#' @param firth_correction Whether to apply Firth's correction for separation issues. Default: TRUE
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{The fitted statistical model object}
#'   \item{results}{Data frame with hazard ratios, confidence intervals, p-values, and median survival times}
#'   \item{km_plot}{A Kaplan-Meier survival curve plot}
#'   \item{reference_group}{The treatment group used as reference}
#'   \item{method_used}{The statistical method used ("cox", "coxphf", or "logrank")}
#' }
#'
#' @details
#' The function adapts to data characteristics:
#' * For well-behaved data: Standard Cox proportional hazards model
#' * For groups with no events: Log-rank test
#' * For groups with few events: Firth's bias-reduced Cox model
#'
#' @importFrom survival Surv survfit coxph survdiff
#' @importFrom stats confint model.matrix as.formula pchisq
#'
#' @examples
#' # Load example data
#' data(combo_treatment_synthetic_data)
#' data_processed <- calculate_volume(combo_treatment_synthetic_data)
#' data_processed <- calculate_dates(data_processed, start_date = "03/24/2025")
#' 
#' # Run survival analysis
#' results <- survival_statistics(
#'   df = data_processed,
#'   reference_group = "Control"
#' )
#' 
#' # Access results
#' print(results$results)  # Hazard ratios, CIs, p-values, and median survival times
#' 
#' # Extract median survival times
#' median_surv <- results$results$Median_Survival
#' names(median_surv) <- results$results$Group
#' print(median_surv)
#'
#' # Display visualizations
#' print(results$km_plot)      # Kaplan-Meier survival curves
#' 
#' # Create a forest plot from the results
#' plot_forest(results$results) # Forest plot of hazard ratios
#'
#' @export
survival_statistics <- function(df,
                              time_column = "Day",
                              censor_column = "Survival_Censor",
                              treatment_column = "Treatment",
                              cage_column = "Cage",
                              id_column = "ID",
                              dose_column = NULL,
                              reference_group = NULL,
                              firth_correction = TRUE) {
  
  # Validate inputs
  validate_inputs(df, time_column, censor_column, treatment_column)
  
  # Setup parameters
  treatment_groups <- unique(df[[treatment_column]])
  if (is.null(reference_group)) {
    reference_group <- sort(treatment_groups)[1]
  } else if (!reference_group %in% treatment_groups) {
    stop("Reference group '", reference_group, "' is not present in the data.")
  }
  message("Using ", reference_group, " as the reference group for hazard ratios.")
  
  # Check cage distribution
  check_cage_distribution(df, treatment_column, cage_column)
  
  # Check for separation issues
  surv_obj <- survival::Surv(df[[time_column]], df[[censor_column]])
  surv_formula_str <- paste("surv_obj ~", treatment_column)
  cox_formula <- stats::as.formula(surv_formula_str)
  separation_info <- check_separation(df, treatment_column, censor_column)
  
  # Choose and fit appropriate model
  model_results <- fit_survival_model(
    df, 
    surv_obj, 
    cox_formula, 
    treatment_column, 
    treatment_groups,
    reference_group,
    time_column,
    censor_column,
    separation_info,
    firth_correction
  )
  
  # Extract model results
  model <- model_results$model
  results <- model_results$results
  method_used <- model_results$method_used
  
  # Create a separate survival fit for median survival calculation
  message("\nCalculating median survival times...")
  surv_formula_str <- paste("Surv(", time_column, ",", censor_column, ") ~ ", 
                            treatment_column)
  surv_formula <- stats::as.formula(surv_formula_str)
  km_fit <- survival::survfit(surv_formula, data = df)
  
  # Display median survival information
  print(km_fit)
  
  # Calculate and add median survival times
  median_survival <- NULL
  tryCatch({
    fit_summary <- summary(km_fit)
    if ("table" %in% names(fit_summary) && is.matrix(fit_summary$table) && 
        "median" %in% colnames(fit_summary$table)) {
      median_survival <- fit_summary$table[, "median"]
      names(median_survival) <- rownames(fit_summary$table)
      
    } else if (!is.null(km_fit$median)) {
      median_survival <- km_fit$median
      names(median_survival) <- names(km_fit$strata)
    }
    
    if (!is.null(median_survival)) {
      # Clean up strata names
      if (!is.null(names(median_survival))) {
        names(median_survival) <- gsub(paste0(treatment_column, "="), "", names(median_survival))
      }
      
      # Add to results data frame
      results$Median_Survival <- median_survival[match(results$Group, names(median_survival))]
      
      # Display the median survival times
      message("\nMedian Survival Times:")
      for (i in seq_along(median_survival)) {
        group_name <- names(median_survival)[i]
        med_surv_val <- median_survival[i]
        if (!is.na(med_surv_val)) {
          message(sprintf("%s: %.1f days", group_name, med_surv_val))
        } else {
          message(sprintf("%s: NA days", group_name))
        }
      }
    }
  }, error = function(e) {
    warning("Error calculating median survival: ", e$message)
  })
  
  # Add Events and Total columns
  # Improved approach to count unique subjects and their events per treatment group
  # First, find the unique subjects (including cage information) in each treatment group
  subject_treatment <- unique(df[, c(id_column, treatment_column, cage_column)])
  total_counts <- table(subject_treatment[[treatment_column]])
  
  # Next, find subjects with events
  # We need to handle possible duplicates in the data (multiple rows per subject)
  # For each subject, if any row has an event, count it as an event
  event_data <- df[, c(id_column, treatment_column, censor_column, cage_column)]
  # Aggregate to get maximum event per subject (1 if any event occurred, 0 otherwise)
  event_by_subject <- stats::aggregate(
    event_data[[censor_column]], 
    by = list(
      ID = event_data[[id_column]], 
      Treatment = event_data[[treatment_column]],
      Cage = event_data[[cage_column]]
    ), 
    FUN = max
  )
  
  # Count events per treatment group
  event_counts <- tapply(event_by_subject$x, event_by_subject$Treatment, sum)
  
  # Assign to results data frame
  results$Events <- event_counts[match(results$Group, names(event_counts))]
  results$Total <- total_counts[match(results$Group, names(total_counts))]
  
  # Calculate event rates for each group
  results$Event_Rate <- results$Events / results$Total
  
  # Verify median survival - if event rate > 0.5 but median is NA, there's likely an issue
  if ("Median_Survival" %in% colnames(results)) {
    # For each group, check if we have > 50% events but NA median
    for (i in 1:nrow(results)) {
      if (is.na(results$Median_Survival[i]) && results$Event_Rate[i] > 0.5) {
        # We should be able to calculate median survival when >50% of subjects have events
        message(sprintf("Group %s has > 50%% events (%.1f%%) but no median survival calculated. Attempting to calculate it now.", 
                        results$Group[i], results$Event_Rate[i] * 100))
        
        # Try to calculate the median for this group
        group_data <- df[df[[treatment_column]] == results$Group[i], ]
        if (nrow(group_data) > 0) {
          # Create a separate survfit object for just this group
          group_surv_formula <- stats::as.formula(paste("Surv(", time_column, ",", censor_column, ") ~ 1"))
          group_km_fit <- survival::survfit(group_surv_formula, data = group_data)
          
          # Extract median (at 0.5)
          if (!is.null(group_km_fit$median)) {
            med_surv <- group_km_fit$median
            if (!is.na(med_surv) && med_surv > 0) {
              results$Median_Survival[i] <- med_surv
              message(sprintf("Successfully calculated median survival for group %s: %.1f days", 
                              results$Group[i], results$Median_Survival[i]))
            } else {
              message(sprintf("Could not calculate valid median survival for group %s despite >50%% events.", 
                              results$Group[i]))
            }
          } else {
            # Try alternate approach using quantiles
            group_quantiles <- summary(group_km_fit)$quantile
            if (!is.null(group_quantiles) && "50%" %in% colnames(group_quantiles)) {
              results$Median_Survival[i] <- group_quantiles["50%"]
              message(sprintf("Calculated median survival using quantiles for group %s: %.1f days", 
                              results$Group[i], results$Median_Survival[i]))
            } else {
              message(sprintf("Could not extract median or quantiles for group %s.", results$Group[i]))
            }
    }
  } else {
          message(sprintf("No data available for group %s to calculate median survival.", results$Group[i]))
        }
      }
    }
  }
  
  # Add reference group note
  results$Note <- ifelse(results$Group == reference_group, "Reference group", "")
  
  # Print formatted results
  print_results(results, df, treatment_column, time_column, censor_column)
  
  # Build our result list
  result_list <- list(
    results = results,
    reference_group = reference_group,
    method_used = method_used,
    survival_data = data.frame(
      Time = df[[time_column]],
      Event = df[[censor_column]],
      Treatment = df[[treatment_column]]
    )
  )
  
  # Add model if it exists
  if (!is.null(model)) {
    result_list$model <- model
  }
  
  # Create visualizations
  tryCatch({
    # Only create KM plot if the survminer package is available
    if (requireNamespace("survminer", quietly = TRUE)) {
      km_plot <- tryCatch({
        create_km_plot(df, time_column, censor_column, treatment_column, id_column)
      }, error = function(e) {
        message("Error creating Kaplan-Meier plot: ", e$message)
        NULL
      })
      
      if (!is.null(km_plot)) {
        result_list$km_plot <- km_plot
      }
    }
  }, error = function(e) {
    message("Error in visualization: ", e$message)
  })
  
  return(result_list)
}

#' Validate Required Inputs
#' @noRd
validate_inputs <- function(df, time_column, censor_column, treatment_column) {
  required_cols <- c(time_column, censor_column, treatment_column)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
}

#' Check Cage Distribution
#' @noRd
check_cage_distribution <- function(df, treatment_column, cage_column) {
  # Only check if cage column exists
  if (cage_column %in% colnames(df)) {
    cage_treatment_table <- table(df[[cage_column]], df[[treatment_column]])
    message("Cage distribution across treatment groups:")
    print(cage_treatment_table)
    
    # Check for collinearity between cage and treatment
    cage_treatment_df <- data.frame(
      Cage = df[[cage_column]],
      Treatment = df[[treatment_column]]
    )
    cage_treatment_counts <- table(cage_treatment_df$Cage, cage_treatment_df$Treatment)
    cage_has_multiple_treatments <- rowSums(cage_treatment_counts > 0) > 1
    if (!any(cage_has_multiple_treatments)) {
      message("Detected collinearity between Cage and Treatment. Using Treatment only in the model.")
    }
  }
}

#' Check for Separation Issues in Survival Data
#' @noRd
check_separation <- function(df, treatment_column, censor_column) {
  # Check for complete separation (groups with all or no events)
  event_by_treatment <- tapply(df[[censor_column]], df[[treatment_column]], function(x) {
    c(sum(x), length(x), sum(x) / length(x))
  })
  
  groups_no_events <- names(event_by_treatment)[sapply(event_by_treatment, function(x) x[1] == 0)]
  groups_all_events <- names(event_by_treatment)[sapply(event_by_treatment, function(x) x[1] == x[2])]
  
  has_separation <- length(groups_no_events) > 0 || length(groups_all_events) > 0
  
  if (has_separation) {
    message("Warning: Some groups have perfect separation (no events). This may affect hazard ratio estimates.")
    if (length(groups_no_events) > 0) {
      message("Groups with no events: ", paste(groups_no_events, collapse = ", "))
    }
    if (length(groups_all_events) > 0) {
      message("Note: Groups with all events: ", paste(groups_all_events, collapse = ", "), 
              " (this is not a problem for Cox models)")
    }
  }
  
  return(list(
    has_separation = has_separation,
    groups_no_events = groups_no_events,
    groups_all_events = groups_all_events
  ))
}

#' Fit Appropriate Survival Model
#' @noRd
fit_survival_model <- function(df, surv_obj, cox_formula, treatment_column, treatment_groups,
                              reference_group, time_column, censor_column, separation_info,
                              firth_correction) {
  
  # Try standard Cox model first
  cox_model <- tryCatch({
    # Create a factor version of the treatment column with the reference level set explicitly
    df$treatment_factor <- factor(df[[treatment_column]], levels = c(reference_group, setdiff(treatment_groups, reference_group)))
    
    # Create a new formula using the factor
    new_formula <- stats::as.formula(paste("surv_obj ~ treatment_factor"))
    
    # Fit model with explicit reference level
    survival::coxph(new_formula, data = df)
  }, error = function(e) {
    message("Standard Cox model failed: ", e$message)
    NULL
  })
  
  # Check for potential issues
  has_issues <- is.null(cox_model) || separation_info$has_separation
  
  if (has_issues && firth_correction) {
    # Use Firth's bias-reduced Cox model
    method_used <- "coxphf"
    message("Using Firth's bias-reduced Cox model: Surv(time, status) ~ group")
    
    results <- tryCatch({
      # Create analysis data frame with safer group naming
      analysis_df <- data.frame(
        time = df[[time_column]],
        status = df[[censor_column]],
        group = factor(df[[treatment_column]], levels = c(reference_group, setdiff(treatment_groups, reference_group)))
      )
      
      # Fit model using the simplified data frame
      if (requireNamespace("coxphf", quietly = TRUE)) {
        model <- coxphf::coxphf(
          survival::Surv(time, status) ~ group, 
          data = analysis_df
        )
        
        # Extract results
      coefs <- model$coefficients
      hazard_ratios <- exp(coefs)
      confidence_intervals <- exp(confint(model))
      p_values <- model$prob
      
        # Create results data frame with all treatment groups
        results <- data.frame(
          Group = treatment_groups,
          HR = NA,
          CI_Lower = NA,
          CI_Upper = NA,
          P_Value = NA,
          stringsAsFactors = FALSE
        )
        
        # Set reference group values
        ref_idx <- which(results$Group == reference_group)
        results$HR[ref_idx] <- 1
        results$CI_Lower[ref_idx] <- 1
        results$CI_Upper[ref_idx] <- 1
        
        # Fill in values for non-reference groups
        for (i in seq_along(coefs)) {
          group_name <- levels(analysis_df$group)[i + 1]  # +1 because first level is reference
          if(group_name %in% results$Group) {
            idx <- which(results$Group == group_name)
            results$HR[idx] <- hazard_ratios[i]
            results$CI_Lower[idx] <- confidence_intervals[i, 1]
            results$CI_Upper[idx] <- confidence_intervals[i, 2]
            results$P_Value[idx] <- p_values[i]
          }
        }
        
        return(list(
          model = model,
          results = results,
          method_used = "coxphf"
        ))
    } else {
        stop("Package 'coxphf' is required but not available")
      }
    }, error = function(e) {
      message("Firth model failed: ", e$message)
      NULL
    })
    
    if (is.null(results)) {
      method_used <- "logrank"
      model <- NULL
      results <- data.frame()
    } else {
      model <- results$model
      results <- results$results
      method_used <- results$method_used
    }
    
  } else if (has_issues) {
    # Use Log-Rank test as fallback
    method_used <- "logrank"
    message("One or more groups have zero events. Using log-rank test to estimate HRs for these groups.")
    
    # Fit Log-Rank test
    surv_diff <- survival::survdiff(cox_formula, data = df)
    print(surv_diff)
    
    # Calculate p-value
    chisq <- surv_diff$chisq
    p_value <- 1 - stats::pchisq(chisq, df = length(treatment_groups) - 1)
    
    # Create basic results without HRs (not estimable in this case)
    results <- data.frame(
      Group = treatment_groups,
      HR = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      P_Value = NA,
      stringsAsFactors = FALSE
    )
    
    # Set reference group values
    ref_idx <- which(results$Group == reference_group)
    results$HR[ref_idx] <- 1
    results$CI_Lower[ref_idx] <- 1
    results$CI_Upper[ref_idx] <- 1
    
    # Set p-value for non-reference groups
    results$P_Value[-ref_idx] <- p_value
    
    model <- surv_diff
    
  } else {
    # Use standard Cox model
    method_used <- "cox"
    model <- cox_model
    
    # Extract the hazard ratios, CIs, and p-values
    model_summary <- summary(model)
    
    # Create results data frame with all treatment groups
    results <- data.frame(
      Group = treatment_groups,
      HR = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      P_Value = NA,
      stringsAsFactors = FALSE
    )
    
    # Set reference group values
    ref_idx <- which(results$Group == reference_group)
    results$HR[ref_idx] <- 1
    results$CI_Lower[ref_idx] <- 1
    results$CI_Upper[ref_idx] <- 1
    
    # Extract coefficient names which should match "treatment_factorTreatmentName"
    coef_names <- rownames(model_summary$coefficients)
    
    # For non-reference groups, extract HR, CI, and p-value
    for (i in seq_along(coef_names)) {
      # Extract treatment group name from coefficient name
      group_name <- gsub("treatment_factor", "", coef_names[i])
      
      # Find corresponding row in results
      idx <- which(results$Group == group_name)
      
      if (length(idx) > 0) {
        # Extract values
        hr <- exp(model_summary$coefficients[i, "coef"])
        ci_lower <- exp(model_summary$coefficients[i, "coef"] - 1.96 * model_summary$coefficients[i, "se(coef)"])
        ci_upper <- exp(model_summary$coefficients[i, "coef"] + 1.96 * model_summary$coefficients[i, "se(coef)"])
        p_value <- model_summary$coefficients[i, "Pr(>|z|)"]
        
        # Assign values
        results$HR[idx] <- hr
        results$CI_Lower[idx] <- ci_lower
        results$CI_Upper[idx] <- ci_upper
        results$P_Value[idx] <- p_value
      }
    }
  }
  
  return(list(
    model = model,
    results = results,
    method_used = method_used
  ))
}

#' Create Kaplan-Meier Plot
#' @noRd
create_km_plot <- function(df, time_column, censor_column, treatment_column, id_column = "ID") {
  # Check for required packages
  if (!requireNamespace("survminer", quietly = TRUE) || 
      !requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    message("Required packages (survminer, survival, ggplot2) not available.")
    return(NULL)
  }
  
  # Make a completely new dataframe with one row per subject
  # First, get unique subjects
  subjects <- unique(df[[id_column]])
  
  # Create a dataframe to hold subject-level data
  subject_data <- data.frame(
    id = character(length(subjects)),
    time = numeric(length(subjects)),
    status = numeric(length(subjects)),
    group = character(length(subjects)),
    stringsAsFactors = FALSE
  )
  
  # For each subject, get their last observation and event status
  for (i in seq_along(subjects)) {
    id <- subjects[i]
    subject_rows <- df[df[[id_column]] == id, ]
    
    # Sort by time to get the last observation
    subject_rows <- subject_rows[order(subject_rows[[time_column]], decreasing = TRUE), ]
    
    # Check if subject had an event (if any row has an event, consider it an event)
    had_event <- any(subject_rows[[censor_column]] == 1)
    
    # Add to subject_data
    subject_data$id[i] <- id
    subject_data$time[i] <- subject_rows[[time_column]][1] # Last observation
    subject_data$status[i] <- ifelse(had_event, 1, 0)
    subject_data$group[i] <- subject_rows[[treatment_column]][1]
  }
  
  # Convert group to factor
  subject_data$group <- factor(subject_data$group)
  
  tryCatch({
    # Fit the survival model with explicit column names
    fit <- survival::survfit(survival::Surv(time, status) ~ group, data = subject_data)
    
    # Create a plot using ggsurvplot
    base_plot <- tryCatch({
      survminer::ggsurvplot(
        fit = fit,
        data = subject_data,
        risk.table = TRUE,
        conf.int = TRUE,
        pval = TRUE
      )
    }, error = function(e) {
      message("Error in creating plot with ggsurvplot: ", e$message)
      
      # Fall back to creating a very simple plot
      fit_summary <- summary(fit)
      plot_data <- data.frame(
        time = fit_summary$time,
        surv = fit_summary$surv,
        group = rep(names(fit$strata), fit$strata)
      )
      
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = surv, color = group)) +
        ggplot2::geom_step() +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Time", y = "Survival Probability", 
                      title = "Kaplan-Meier Survival Curve")
      
      return(list(plot = p))
    })
    
    return(base_plot)
    
  }, error = function(e) {
    message("Error in survival fit or plotting: ", e$message)
    message("Data dimensions: ", nrow(subject_data), " x ", ncol(subject_data))
    return(NULL)
  })
}

#' Print Formatted Results
#' @noRd
print_results <- function(results, df = NULL, treatment_column = NULL, time_column = NULL, censor_column = NULL) {
  message("\nSurvival Analysis Results:")
  message("=======================")
  
  # Debug output removed for cleaner presentation
  
  for(i in 1:nrow(results)) {
    message(sprintf("\nGroup: %s", results$Group[i]))
    
    # Safely handle HR values
    hr_na <- is.na(results$HR[i])
    ci_lower_na <- is.na(results$CI_Lower[i]) 
    ci_upper_na <- is.na(results$CI_Upper[i])
    
    hr_text <- if(hr_na || (!hr_na && results$HR[i] == 0)) {
      "Hazard Ratio: Not estimable"
    } else {
      sprintf("Hazard Ratio: %.3f (%.3f-%.3f)", 
              results$HR[i], 
              results$CI_Lower[i], 
              results$CI_Upper[i])
    }
    message(hr_text)
    
    if (!is.na(results$P_Value[i])) {
      message(sprintf("P-value: %.4f", results$P_Value[i]))
    }
    
    if ("Median_Survival" %in% colnames(results)) {
      if (!is.na(results$Median_Survival[i])) {
        message(sprintf("Median Survival: %.1f days", results$Median_Survival[i]))
      } else {
        # Check if we have event rate information
        if ("Event_Rate" %in% colnames(results) && !is.na(results$Event_Rate[i])) {
          # If more than 50% of subjects had events, try to calculate it
          if (results$Event_Rate[i] > 0.5) {
            # We need to calculate the median survival since we have >50% events
            if (!is.null(df) && !is.null(treatment_column) && !is.null(time_column) && !is.null(censor_column)) {
              group_data <- df[df[[treatment_column]] == results$Group[i], ]
              if (nrow(group_data) > 0) {
                # Create a survfit object for just this group
                group_surv_formula <- stats::as.formula(paste("Surv(", time_column, ",", censor_column, ") ~ 1"))
                group_km_fit <- survival::survfit(group_surv_formula, data = group_data)
                
                # Try to extract the median survival
                if (!is.null(group_km_fit$median)) {
                  med_surv <- group_km_fit$median
                  if (!is.na(med_surv) && med_surv > 0) {
                    # Update the results data frame with the calculated median
                    results$Median_Survival[i] <- med_surv
                    message(sprintf("Median Survival: %.1f days", med_surv))
                  } else {
                    message("Median Survival: Error calculating median")
                  }
                } else {
                  message("Median Survival: Error extracting median from survfit")
                }
              } else {
                message("Median Survival: No data available for group")
              }
            } else {
              message("Median Survival: Required data for calculation not provided")
            }
          } else {
            # Less than 50% had events, so "Not Reached" is accurate
            message("Median Survival: Not reached")
          }
        } else {
          # If we don't have event rate info, use original behavior
          message("Median Survival: Not reached")
        }
      }
    }
    
    # Ensure event counts are properly displayed 
    if (!is.na(results$Events[i]) && !is.na(results$Total[i])) {
      message(sprintf("Events: %d/%d", results$Events[i], results$Total[i]))
    }
    
    if(!is.na(results$Note[i]) && results$Note[i] != "") {
      message(sprintf("Note: %s", results$Note[i]))
    }
  }
  message("\n")
  
  # Print summary table
  message("Summary Table:")
  formatted_table <- data.frame(
    Group = results$Group,
    "HR (95% CI)" = sapply(1:nrow(results), function(i) {
      if(is.na(results$HR[i]) || (!is.na(results$HR[i]) && results$HR[i] == 0)) {
        "Not estimable"
      } else {
        sprintf("%.2f (%.2f-%.2f)", 
                results$HR[i], 
                results$CI_Lower[i], 
                results$CI_Upper[i])
      }
    }),
    "P-value" = sapply(1:nrow(results), function(i) {
      ifelse(is.na(results$P_Value[i]), "Ref", sprintf("%.4f", results$P_Value[i]))
    }),
    "Events/Total" = sapply(1:nrow(results), function(i) {
      if (!is.na(results$Events[i]) && !is.na(results$Total[i])) {
        sprintf("%d/%d", results$Events[i], results$Total[i])
      } else {
        "NA/NA"
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Add median survival to table if available
  if ("Median_Survival" %in% colnames(results)) {
    formatted_table$"Median Survival" <- sapply(1:nrow(results), function(i) {
      if (is.na(results$Median_Survival[i])) {
        # Check if we have event rate information
        if ("Event_Rate" %in% colnames(results) && !is.na(results$Event_Rate[i])) {
          # If more than 50% of subjects had events, try to calculate it
          if (results$Event_Rate[i] > 0.5) {
            # We need to calculate the median survival since we have >50% events
            if (!is.null(df) && !is.null(treatment_column) && !is.null(time_column) && !is.null(censor_column)) {
              group_data <- df[df[[treatment_column]] == results$Group[i], ]
              if (nrow(group_data) > 0) {
                # Create a survfit object for just this group
                group_surv_formula <- stats::as.formula(paste("Surv(", time_column, ",", censor_column, ") ~ 1"))
                group_km_fit <- survival::survfit(group_surv_formula, data = group_data)
                
                # Try to extract the median survival
                if (!is.null(group_km_fit$median)) {
                  med_surv <- group_km_fit$median
                  if (!is.na(med_surv) && med_surv > 0) {
                    # Update the results data frame with the calculated median
                    results$Median_Survival[i] <- med_surv
                    return(sprintf("%.1f days", med_surv))
                  }
                }
                
                # If we're here, we couldn't calculate it despite having >50% events
                return("Error calculating median")
              } else {
                return("No data available")
              }
            } else {
              return("Required data not provided")
            }
          } else {
            # Less than 50% had events, so "Not Reached" is accurate
            return("Not reached")
          }
        } else {
          # If we don't have event rate info, use original behavior
          return("Not reached")
        }
      } else {
        return(sprintf("%.1f days", results$Median_Survival[i]))
      }
    })
  }
  
  print(formatted_table)
  
  # Return the formatted table invisibly for further use if needed
  invisible(formatted_table)
}