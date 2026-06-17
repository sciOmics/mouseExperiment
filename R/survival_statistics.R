# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

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
#' @param verbose Whether to print analysis details to the console. Default: TRUE
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{The fitted statistical model object}
#'   \item{results}{Data frame with hazard ratios, confidence intervals, p-values, and median survival times}
#'   \item{reference_group}{The treatment group used as reference}
#'   \item{method_used}{The statistical method used ("cox", "coxphf", or "logrank")}
#'   \item{ph_test}{A \code{cox.zph} object (Schoenfeld-residual proportional-hazards
#'     test) when the standard Cox path is used; \code{NULL} for the Firth and
#'     log-rank fallbacks. A small global p-value indicates the PH assumption is
#'     violated and hazard ratios should be interpreted with caution.}
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
#' @export
survival_statistics <- function(df,
                              time_column = "Day",
                              censor_column = "Survival_Censor",
                              treatment_column = "Treatment",
                              cage_column = "Cage",
                              id_column = "ID",
                              dose_column = NULL,
                              reference_group = NULL,
                              firth_correction = TRUE,
                              verbose = TRUE) {
  
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
    firth_correction,
    verbose = verbose
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
  if (isTRUE(verbose)) message(paste(utils::capture.output(print(km_fit)), collapse = "\n"))
  
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
  if (isTRUE(verbose)) print_results(results, df, treatment_column, time_column, censor_column)
  
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

  # Add cox.zph proportional-hazards check when available (cox path only)
  if (!is.null(model_results$ph_test)) {
    result_list$ph_test <- model_results$ph_test
  }

  # Add concordance / C-index when available (cox path only)
  if (!is.null(model_results$c_index)) {
    result_list$c_index <- model_results$c_index
  }

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
    message(paste(utils::capture.output(print(cage_treatment_table)), collapse = "\n"))
    
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
                              firth_correction, verbose = TRUE) {
  
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

        # Fallback: for non-reference groups where coxphf failed to converge
        # (p_value is NA), compute a pairwise log-rank p-value instead.
        non_ref_na <- which(results$Group != reference_group & is.na(results$P_Value))
        for (idx in non_ref_na) {
          grp <- results$Group[idx]
          pair_df <- df[df[[treatment_column]] %in% c(reference_group, grp), , drop = FALSE]
          pair_df[[treatment_column]] <- factor(pair_df[[treatment_column]])
          tryCatch({
            pair_formula <- stats::as.formula(
              paste0("survival::Surv(", time_column, ", ", censor_column, ") ~ ", treatment_column)
            )
            lr <- survival::survdiff(pair_formula, data = pair_df)
            results$P_Value[idx] <- round(1 - stats::pchisq(lr$chisq, df = 1), 4)
            message(sprintf("  Pairwise log-rank p-value used for %s (coxphf did not converge): %.4f",
                            grp, results$P_Value[idx]))
          }, error = function(e) {
            message(sprintf("  Could not compute fallback p-value for %s: %s", grp, e$message))
          })
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
      method_used <- results$method_used
      model <- results$model
      results <- results$results
    }
    
  } else if (has_issues) {
    # Use pairwise log-rank tests as fallback (one per non-reference group).
    # An omnibus logrank p-value must not be assigned to all groups — it is a
    # single test for any difference and is not a valid per-comparison p-value.
    method_used <- "logrank"
    message("One or more groups have zero events. Using pairwise log-rank tests.")

    # Omnibus test (stored separately for reference; not used as per-group p-values)
    surv_diff <- survival::survdiff(cox_formula, data = df)
    if (isTRUE(verbose)) message(paste(utils::capture.output(print(surv_diff)), collapse = "\n"))
    omnibus_p <- 1 - stats::pchisq(surv_diff$chisq, df = length(treatment_groups) - 1L)

    # Create basic results without HRs (not estimable when a group has zero events)
    results <- data.frame(
      Group    = treatment_groups,
      HR       = NA,
      CI_Lower = NA,
      CI_Upper = NA,
      P_Value  = NA,
      stringsAsFactors = FALSE
    )

    ref_idx <- which(results$Group == reference_group)
    results$HR[ref_idx]       <- 1
    results$CI_Lower[ref_idx] <- 1
    results$CI_Upper[ref_idx] <- 1

    # Pairwise log-rank: compare each treatment group vs. reference individually
    for (grp in treatment_groups[treatment_groups != reference_group]) {
      pair_data <- df[df[[treatment_column]] %in% c(reference_group, grp), ]
      pair_p <- tryCatch({
        pair_diff <- survival::survdiff(cox_formula, data = pair_data)
        1 - stats::pchisq(pair_diff$chisq, df = 1L)
      }, error = function(e) omnibus_p)
      results$P_Value[results$Group == grp] <- pair_p
    }

    model <- surv_diff
    
  } else {
    # Use standard Cox model
    method_used <- "cox"
    model <- cox_model

    # Proportional-hazards check via Schoenfeld residuals (cox.zph).
    # Surface global + per-covariate test for downstream display.
    ph_test <- tryCatch(survival::cox.zph(model), error = function(e) NULL)
    if (!is.null(ph_test) && isTRUE(verbose)) {
      global_p <- tryCatch(ph_test$table["GLOBAL", "p"], error = function(e) NA_real_)
      if (is.finite(global_p) && global_p < 0.05) {
        message(sprintf("Proportional-hazards check: global p = %.4f (< 0.05) — PH assumption may be violated.", global_p))
      }
    }

    # Concordance (C-index) — survival analogue of AUC-ROC; values around
    # 0.5 indicate no discrimination, 1.0 perfect discrimination.
    c_index <- tryCatch(
      survival::concordance(model),
      error = function(e) NULL
    )

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
        # Use summary(coxph)$conf.int directly — already contains exp(coef),
        # lower .95, and upper .95 at the proper qnorm(0.975) ≈ 1.959964.
        ci_row   <- model_summary$conf.int[i, , drop = TRUE]
        hr       <- unname(ci_row["exp(coef)"])
        ci_lower <- unname(ci_row["lower .95"])
        ci_upper <- unname(ci_row["upper .95"])
        p_value  <- model_summary$coefficients[i, "Pr(>|z|)"]
        
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
    method_used = method_used,
    ph_test = if (exists("ph_test", inherits = FALSE)) ph_test else NULL,
    c_index = if (exists("c_index", inherits = FALSE)) c_index else NULL
  ))
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
        message("Median Survival: Not reached")
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
      is_ref <- !is.na(results$Note[i]) && results$Note[i] == "Reference group"
      if (is_ref) {
        "Ref"
      } else if (is.na(results$P_Value[i])) {
        "NC"  # Not Converged
      } else {
        sprintf("%.4f", results$P_Value[i])
      }
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
    formatted_table$"Median Survival" <- sapply(seq_len(nrow(results)), function(i) {
      if (is.na(results$Median_Survival[i])) "Not reached"
      else sprintf("%.1f days", results$Median_Survival[i])
    })
  }
  
  message(paste(utils::capture.output(print(formatted_table)), collapse = "\n"))
  
  # Return the formatted table invisibly for further use if needed
  invisible(formatted_table)
}