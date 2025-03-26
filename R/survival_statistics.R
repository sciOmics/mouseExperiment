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
#'   \item{forest_plot}{A forest plot visualizing hazard ratios with confidence intervals}
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
#' print(results$forest_plot)  # Forest plot of hazard ratios
#' print(results$km_plot)      # Kaplan-Meier survival curves
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
  cox_formula <- as.formula(paste("surv_obj ~", treatment_column))
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
  surv_formula <- as.formula(paste("Surv(", time_column, ",", censor_column, ") ~ ", 
                                  treatment_column))
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
        message(sprintf("%s: %.1f days", names(median_survival)[i], median_survival[i]))
      }
    }
  }, error = function(e) {
    warning("Error calculating median survival: ", e$message)
  })
  
  # Add Events and Total columns
  event_counts <- tapply(df[[censor_column]], df[[treatment_column]], sum)
  total_counts <- table(df[[treatment_column]])
  results$Events <- event_counts[match(results$Group, names(event_counts))]
  results$Total <- total_counts[match(results$Group, names(total_counts))]
  
  # Add reference group note
  results$Note <- ifelse(results$Group == reference_group, "Reference group", "")
  
  # Print formatted results
  print_results(results)
  
  # Create visualization plots
  km_plot <- create_km_plot(df, surv_obj, cox_formula, treatment_column)
  forest_plot <- create_forest_plot(results, title = "Hazard Ratios with 95% CIs")
  
  # Return all results
  return(list(
    model = model,
    results = results,
    forest_plot = forest_plot,
    km_plot = km_plot,
    reference_group = reference_group,
    method_used = method_used
  ))
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
    survival::coxph(cox_formula, data = df)
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
      fit_firth_model(df, treatment_column, time_column, censor_column, 
                     treatment_groups, reference_group)
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
      HR = ifelse(treatment_groups == reference_group, 1, NA),
      CI_Lower = ifelse(treatment_groups == reference_group, 1, NA),
      CI_Upper = ifelse(treatment_groups == reference_group, 1, NA),
      P_Value = ifelse(treatment_groups == reference_group, NA, p_value),
      stringsAsFactors = FALSE
    )
    
    model <- surv_diff
    
  } else {
    # Use standard Cox model
    method_used <- "cox"
    model <- cox_model
    
    # Extract the hazard ratios, CIs, and p-values
    model_summary <- summary(model)
    results <- extract_cox_results(model_summary, treatment_groups, reference_group)
  }
  
  return(list(
    model = model,
    results = results,
    method_used = method_used
  ))
}

#' Fit Firth's Bias-Reduced Cox Model
#' @noRd
fit_firth_model <- function(df, treatment_column, time_column, censor_column, 
                           treatment_groups, reference_group) {
  if (!requireNamespace("coxphf", quietly = TRUE)) {
    stop("Package 'coxphf' is required for Firth correction. Please install it.")
  }
  
  # Create a model matrix
  mm <- stats::model.matrix(~ df[[treatment_column]] - 1)
  colnames(mm) <- gsub("df\\[\\[treatment_column\\]\\]", "", colnames(mm))
  
  # Prepare data and formula
  df_coxphf <- cbind(df[, c(time_column, censor_column)], mm)
  var_names <- colnames(mm)
  formula_str <- paste("survival::Surv(", time_column, ", ", censor_column, ") ~ ", 
                      paste(var_names, collapse = " + "))
  coxphf_formula <- as.formula(formula_str)
  
  # Fit the model
  model <- coxphf::coxphf(coxphf_formula, data = df_coxphf)
  print(model)
  
  # Extract results
  coefs <- model$coefficients
  hazard_ratios <- exp(coefs)
  confidence_intervals <- exp(confint(model))
  p_values <- model$prob
  
  # Fix any naming inconsistencies
  groups <- var_names
  groups <- gsub("\\bPD1\\b", "aPD1", groups)
  
  # Create results data frame
  results <- data.frame(
    Group = treatment_groups,
    HR = ifelse(treatment_groups == reference_group, 1, hazard_ratios),
    CI_Lower = ifelse(treatment_groups == reference_group, 1, confidence_intervals[, 1]),
    CI_Upper = ifelse(treatment_groups == reference_group, 1, confidence_intervals[, 2]),
    P_Value = ifelse(treatment_groups == reference_group, NA, p_values),
    stringsAsFactors = FALSE
  )
  
  # Ensure proper ordering
  results <- results[match(treatment_groups, results$Group), ]
  rownames(results) <- NULL
  
  return(list(
    model = model,
    results = results
  ))
}

#' Extract Results from Standard Cox Model
#' @noRd
extract_cox_results <- function(model_summary, treatment_groups, reference_group) {
  # Extract statistics
  hazard_ratios <- exp(model_summary$coefficients[, "coef"])
  ci_lower <- exp(model_summary$coefficients[, "coef"] - 1.96 * model_summary$coefficients[, "se(coef)"])
  ci_upper <- exp(model_summary$coefficients[, "coef"] + 1.96 * model_summary$coefficients[, "se(coef)"])
  p_values <- model_summary$coefficients[, "Pr(>|z|)"]
  
  # Handle group names
  groups <- names(hazard_ratios)
  groups <- gsub("\\bPD1\\b", "aPD1", groups)
  
  # Create results data frame
  results <- data.frame(
    Group = treatment_groups,
    HR = ifelse(treatment_groups == reference_group, 1, hazard_ratios),
    CI_Lower = ifelse(treatment_groups == reference_group, 1, ci_lower),
    CI_Upper = ifelse(treatment_groups == reference_group, 1, ci_upper),
    P_Value = ifelse(treatment_groups == reference_group, NA, p_values),
    stringsAsFactors = FALSE
  )
  
  # Ensure proper ordering
  results <- results[match(treatment_groups, results$Group), ]
  rownames(results) <- NULL
  
  return(results)
}

#' Create Kaplan-Meier Plot
#' @noRd
create_km_plot <- function(df, surv_obj, cox_formula, treatment_column) {
  if (!requireNamespace("survminer", quietly = TRUE)) {
    message("Package 'survminer' not available. Kaplan-Meier plot not created.")
    return(NULL)
  }
  
  # Create a fit object for the K-M curve
  km_fit <- survival::survfit(cox_formula, data = df)
  
  # Create the K-M plot
  survminer::ggsurvplot(
    km_fit,
    data = df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    legend.labs = levels(as.factor(df[[treatment_column]])),
    xlab = "Time",
    risk.table.height = 0.25,
    ggtheme = ggplot2::theme_minimal()
  )
}

#' Print Formatted Results
#' @noRd
print_results <- function(results) {
  message("\nSurvival Analysis Results:")
  message("=======================")
  
  # Debug median survival data
  if ("Median_Survival" %in% colnames(results)) {
    message("Median survival data found in results")
    message(paste("Median survival values:", paste(results$Median_Survival, collapse=", ")))
  } else {
    message("No median survival column in results")
  }
  
  for(i in 1:nrow(results)) {
    message(sprintf("\nGroup: %s", results$Group[i]))
    
    hr_text <- if(is.na(results$HR[i]) || results$HR[i] == 0) {
      "Hazard Ratio: Not estimable"
    } else {
      sprintf("Hazard Ratio: %.3f (%.3f-%.3f)", 
             results$HR[i], results$CI_Lower[i], results$CI_Upper[i])
    }
    message(hr_text)
    
    if (!is.na(results$P_Value[i])) {
      message(sprintf("P-value: %.4f", results$P_Value[i]))
    }
    
    if ("Median_Survival" %in% colnames(results)) {
      if (!is.na(results$Median_Survival[i])) {
        message(sprintf("Median Survival: %.1f days", results$Median_Survival[i]))
      } else {
        message("Median Survival: Not reached")
      }
    }
    
    message(sprintf("Events: %d/%d", results$Events[i], results$Total[i]))
    
    if(results$Note[i] != "") {
      message(sprintf("Note: %s", results$Note[i]))
    }
  }
  message("\n")
  
  # Print summary table
  message("Summary Table:")
  formatted_table <- data.frame(
    Group = results$Group,
    "HR (95% CI)" = sapply(1:nrow(results), function(i) {
      if(is.na(results$HR[i]) || results$HR[i] == 0) {
        "Not estimable"
      } else {
        sprintf("%.2f (%.2f-%.2f)", 
               results$HR[i], results$CI_Lower[i], results$CI_Upper[i])
      }
    }),
    "P-value" = sapply(1:nrow(results), function(i) {
      ifelse(is.na(results$P_Value[i]), "Ref", sprintf("%.4f", results$P_Value[i]))
    }),
    "Events/Total" = sprintf("%d/%d", results$Events, results$Total),
    stringsAsFactors = FALSE
  )
  
  # Add median survival to table if available
  if ("Median_Survival" %in% colnames(results)) {
    formatted_table$"Median Survival" <- sapply(1:nrow(results), function(i) {
      ifelse(is.na(results$Median_Survival[i]), 
            "Not reached", 
            sprintf("%.1f days", results$Median_Survival[i]))
    })
  }
  
  print(formatted_table)
}

#' Create Forest Plot for Hazard Ratios
#' @noRd
create_forest_plot <- function(results, title = "Forest Plot") {
  # Check if we have valid results for a forest plot
  if (is.null(results) || !all(c("Group", "HR", "CI_Lower", "CI_Upper") %in% colnames(results))) {
    return(NULL)
  }
  
  # Remove any rows with NA values for plotting
  plot_data <- results[!is.na(results$HR), ]
  
  # If no valid rows remain, return NULL
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Format hazard ratios and CIs
  plot_data$HR_CI <- sprintf("%.2f (%.2f-%.2f)", 
                            plot_data$HR, plot_data$CI_Lower, plot_data$CI_Upper)
  
  # Create the forest plot
  ggplot2::ggplot(plot_data, ggplot2::aes(x = HR, y = Group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_x_continuous(trans = "log10", breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    ggplot2::labs(
      title = title,
      x = "Hazard Ratio (log scale)",
      y = "Treatment Group"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = max(HR) * 1.5, label = HR_CI),
      hjust = 0
    )
}