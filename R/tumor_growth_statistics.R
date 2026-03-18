# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Analyze Tumor Growth Using Various Statistical Methods
#'
#' @importFrom utils head tail
#' @importFrom rlang sym !!
#' @importFrom dplyr group_by summarize mutate arrange filter
#'
#' @description
#' This function provides a comprehensive statistical analysis for tumor growth data using
#' different statistical approaches. It supports multiple methods including linear mixed-effects
#' models, generalized additive mixed models, and area under the curve analysis.
#'
#' @param df A data frame containing tumor growth data. Must include columns for time, volume,
#'        treatment group, cage identifier, and individual subject ID. Additional columns can be included.
#' @param time_column A character string specifying the column name for time points (e.g., "Day"). 
#'        This column should contain numeric values representing time since treatment started. Default is "Day".
#' @param volume_column A character string specifying the column name for tumor volume measurements (e.g., "Volume").
#'        This column should contain numeric tumor volume measurements. Default is "Volume".
#' @param treatment_column A character string specifying the column name for treatment groups (e.g., "Treatment").
#'        This column should contain categorical treatment identifiers. Default is "Treatment".
#' @param cage_column A character string specifying the column name for the cage identifier (e.g., "Cage").
#'        This column should contain cage identifiers. The function automatically tests for collinearity with treatment.
#'        Default is "Cage".
#' @param id_column A character string specifying the column name for individual subject identifiers (e.g., "ID").
#'        This column should contain unique identifiers for each animal/subject. Default is "ID".
#' @param dose_column Optional. A character string specifying the column name for dose levels, if available.
#'        This allows for dose-response analysis. Default is NULL.
#' @param transform A character string specifying the transformation to apply to volume data. 
#'        Options are "log", "sqrt", "none". Default is "log", which is recommended for exponential growth.
#' @param polynomial_degree Integer specifying the polynomial degree for time (day) effects in the model. 
#'        Default is 1 (linear). Increase to 2 or 3 to model non-linear growth patterns.
#' @param model_type A character string specifying the type of model to fit. Options are: 
#'        "lme4" (standard linear mixed effects model using lme4 package),
#'        "auc" (area under the curve analysis).
#'        Default is "lme4".
#' @param random_effects_specification A character string specifying the random effects structure.
#'        "intercept_only" (default): (1|ID) - random intercepts by subject
#'        "slope": (Day|ID) - random intercepts and slopes by subject
#'        "none": No random effects
#' @param handle_cage_effects Method for handling cage effects: "include_if_not_collinear", "always_include", 
#'        "never_include", or "as_random_effect". Default is "include_if_not_collinear".
#' @param auc_method Method for AUC calculation: "trapezoidal" or "last_observation". Default is "trapezoidal".
#' @param reference_group Optional. A character string specifying which treatment group should be used as the reference
#'        for statistical comparisons. If NULL, the first treatment group alphabetically will be used.
#' @param return_model Boolean. Should the full fitted model be returned? Default is TRUE.
#' @param include_diagnostics Boolean. Should model diagnostic information be included? Default is TRUE.
#' @param plots Boolean. Should standard plots be generated and returned? Default is TRUE.
#' @param verbose Boolean. Should detailed information be printed during analysis? Default is FALSE.
#' @param extrapolation_points Number of points to extrapolate for each subject (default: 0)
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{model}{The fitted statistical model (lme4 or auc object)}
#'   \item{anova}{ANOVA table for fixed effects with Type III tests}
#'   \item{summary}{Detailed model summary with parameter estimates}
#'   \item{pairwise_comparisons}{Results of pairwise comparisons between treatment groups, showing treatment differences, p-values, and confidence intervals}
#'   \item{treatment_effects}{Estimated treatment effects on tumor growth, showing adjusted means for each treatment group}
#'   \item{growth_rates}{Data frame containing growth rates for each subject, calculated as the slope of log-volume over time. Higher values indicate faster tumor growth.}
#'   \item{cage_analysis}{Analysis of cage effects, including:
#'     \itemize{
#'       \item{collinearity_test}{Chi-squared test result for collinearity between cage and treatment}
#'       \item{effects}{Data frame of cage-level statistics (mean and SD of volume by cage and treatment)}
#'     }
#'   }
#'   \item{model_selection}{Results of model selection process, including:
#'     \itemize{
#'       \item{aic}{AIC values for each model specification}
#'       \item{bic}{BIC values for each model specification}
#'       \item{selected_model}{Name of the model with lowest BIC (most parsimonious)}
#'     }
#'   }
#'   \item{diagnostics}{Model diagnostic information including residuals, random effects, and variance components}
#'   \item{auc_data}{Area under the curve data and statistics (when model_type="auc")}
#'   \item{plots}{List of standard plots for visualizing results}
#'   \item{data_summary}{Descriptive statistics of the processed data by treatment and time point, including mean, SD, and SEM of volumes}
#' }
#'
#' @details
#' This function offers several statistical approaches for analyzing tumor growth data. The choice of model
#' depends on the experimental design, growth patterns, and research questions:
#'
#' 1. Linear mixed-effects models (model_type="lme4") are suitable for most tumor growth experiments and
#'    account for the correlation structure of repeated measurements on the same subjects. 
#'    The default formula is log(Volume) ~ Day * Treatment + (1|ID).
#'
#' 2. Area under the curve analysis (model_type="auc") reduces each subject's growth curve to a single
#'    summary metric (AUC), which is then compared between treatment groups.
#'
#' The function automatically handles:
#' - Missing data patterns using appropriate methods
#' - Checking for and addressing cage effects
#' - Transformations to normalize growth curves
#' - Polynomial terms for non-linear growth patterns
#' - Diagnostic checks of model assumptions
#' - Growth rate analysis for each treatment group
#' - Cage effect analysis and collinearity testing
#' - Model selection based on AIC/BIC criteria
#' - Treatment effect estimation with proper reference group handling
#'
#' @examples
#' # Load example data
#' data(combo_treatment_synthetic_data)
#' tumor_data <- calculate_volume(combo_treatment_synthetic_data)
#' tumor_data <- calculate_dates(tumor_data, start_date = "03/24/2025")
#' 
#' # Basic analysis with default settings
#' results <- tumor_growth_statistics(tumor_data)
#' 
#' # View ANOVA table to assess treatment effects
#' print(results$anova)
#' 
#' # View pairwise comparisons
#' print(results$pairwise_comparisons)
#' 
#' # Plot adjusted means for each treatment group
#' print(results$plots$adjusted_means)
#' 
#' # Analyze with different model specification
#' results_gam <- tumor_growth_statistics(
#'   tumor_data,
#'   model_type = "gam",
#'   transform = "log"
#' )
#'
#' @export
tumor_growth_statistics <- function(df,
                                  time_column = "Day",
                                  volume_column = "Volume",
                                  treatment_column = "Treatment",
                                  cage_column = "Cage",
                                  id_column = "ID",
                                  dose_column = NULL,
                                  transform = c("log", "sqrt", "none"),
                                  polynomial_degree = 1,
                                  model_type = c("lme4", "auc"),
                                  random_effects_specification = c("intercept_only", "slope", "none"),
                                  handle_cage_effects = c("include_if_not_collinear", "always_include", 
                                                        "never_include", "as_random_effect"),
                                  auc_method = c("trapezoidal", "last_observation"),
                                  reference_group = NULL,
                                  return_model = TRUE,
                                  include_diagnostics = TRUE,
                                  plots = TRUE,
                                  verbose = FALSE,
                                  extrapolation_points = 0) {
  # Check for required packages
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Please install the lme4 package: install.packages('lme4')")
  }
  
  # Match arguments
  transform <- match.arg(transform)
  model_type <- match.arg(model_type, c("lme4", "auc"))
  random_effects_specification <- match.arg(random_effects_specification)
  handle_cage_effects <- match.arg(handle_cage_effects)
  auc_method <- match.arg(auc_method)
  
  if (verbose) {
    cat("Analyzing tumor growth data...\n")
    cat("Model type:", model_type, "\n")
    cat("Transform:", transform, "\n")
  }
  
  # Check if required columns exist
  required_cols <- c(time_column, volume_column, treatment_column, id_column, cage_column)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Set reference group if not specified
  treatment_groups <- unique(df[[treatment_column]])
  if (is.null(reference_group)) {
    reference_group <- sort(treatment_groups)[1]
  } else if (!reference_group %in% treatment_groups) {
    stop("Reference group '", reference_group, "' is not present in the data.")
  }
  if (verbose) {
    cat("Using", reference_group, "as reference group for statistical comparisons\n")
  }
  
  # Extrapolate data points if requested
  if (extrapolation_points > 0) {
    if (verbose) cat("Extrapolating points for subjects with missing data at the last timepoint\n")
    
    # Find the true maximum day across all subjects (global maximum day of the study)
    true_max_day <- max(df[[time_column]], na.rm = TRUE)
    if (verbose) cat("True maximum day of the study:", true_max_day, "\n")
    
    # Make sure all data has the Extrapolated column
    if (!"Extrapolated" %in% colnames(df)) {
      df$Extrapolated <- FALSE
    }
    
    # Create a list to store results for each subject
    subjects_with_extrapolation <- list()
    
    # Process each unique subject
    unique_subjects <- unique(paste(df[[id_column]], df[[treatment_column]], df[[cage_column]], sep="__"))
    
    for (subject_id in unique_subjects) {
      # Parse the composite ID
      id_parts <- strsplit(subject_id, "__")[[1]]
      id <- id_parts[1]
      treatment <- id_parts[2]
      cage <- id_parts[3]
      
      # Get data for this subject
      subject_data <- df[df[[id_column]] == id & 
                        df[[treatment_column]] == treatment & 
                        df[[cage_column]] == cage, ]
      
      # Get the max day for this subject
      max_subj_day <- max(subject_data[[time_column]])
      
      # Only extrapolate if the subject doesn't have data on the true max day
      if (max_subj_day < true_max_day) {
        # Need at least 2 points for extrapolation
        if (nrow(subject_data) >= 2) {
          # Use the last 3 points (or all if less than 3) to fit a linear model
          n_points <- min(3, nrow(subject_data))
          subject_data <- subject_data[order(subject_data[[time_column]]), ]
          last_points <- tail(subject_data, n_points)
          
          # Try to fit model
          tryCatch({
            lm_fit <- stats::lm(paste(volume_column, "~", time_column), data = last_points)
            
            # Predict at true_max_day
            new_data <- data.frame(time = true_max_day)
            names(new_data) <- time_column
            predicted_volume <- max(0, as.numeric(predict(lm_fit, newdata = new_data)))
            
            # Create a new row for the extrapolated point
            new_row <- subject_data[1, ]
            new_row[[time_column]] <- true_max_day
            new_row[[volume_column]] <- predicted_volume
            new_row$Extrapolated <- TRUE
            
            # Add the new extrapolated point to the subject data
            subject_data <- rbind(subject_data, new_row)
            
            if (verbose) {
              cat("Extrapolated subject", id, "from day", max_subj_day, "to day", true_max_day, "\n")
            }
          }, error = function(e) {
            if (verbose) {
              cat("Failed to extrapolate subject", id, ":", conditionMessage(e), "\n")
            }
          })
        }
      }
      
      # Store the processed subject data
      subjects_with_extrapolation[[subject_id]] <- subject_data
    }
    
    # Combine all subject data
    df <- do.call(rbind, subjects_with_extrapolation)
    
    # Count extrapolated subjects for verbose output
    if (verbose) {
      extrapolated_subjects <- unique(df$Extrapolated[df$Extrapolated])
      n_extrapolated <- length(extrapolated_subjects)
      if (n_extrapolated > 0) {
        cat("Successfully extrapolated", n_extrapolated, "subjects to day", true_max_day, "\n")
      } else {
        cat("No subjects needed or qualified for extrapolation to day", true_max_day, "\n")
      }
    }
  }
  
  # Create a copy of the data for analysis
  analysis_df <- df
  
  # Create a copy for AUC calculation before transformation
  auc_df <- df
  
  # Apply transformations if needed
  if (transform == "log") {
    # Add a small constant to avoid log(0) issues
    analysis_df[[volume_column]] <- log(analysis_df[[volume_column]] + 1)
    if (verbose) cat("Applied log(x+1) transformation to volume data\n")
  } else if (transform == "sqrt") {
    analysis_df[[volume_column]] <- sqrt(analysis_df[[volume_column]])
    if (verbose) cat("Applied square root transformation to volume data\n")
  }
  
  # Growth rate analysis
  growth_rates <- tryCatch({
    # Split data by treatment, ID, and cage to ensure unique subjects
    split_data <- split(analysis_df, list(analysis_df[[treatment_column]], 
                                        analysis_df[[id_column]], 
                                        analysis_df[[cage_column]]))
    
    # Calculate growth rates for each subject
    growth_rates_list <- lapply(split_data, function(subject_data) {
      if (nrow(subject_data) >= 3) {
        # Sort by time
        subject_data <- subject_data[order(subject_data[[time_column]]), ]
        
        # Calculate log volume
        log_volume <- log1p(subject_data[[volume_column]])
        
        # Fit linear model
        model <- stats::lm(log_volume ~ subject_data[[time_column]])
        growth_rate <- stats::coef(model)[2]
        
        # Return data frame with results including cage information
        data.frame(
          Treatment = unique(subject_data[[treatment_column]]),
          ID = unique(subject_data[[id_column]]),
          Cage = unique(subject_data[[cage_column]]),
          growth_rate = growth_rate
        )
      } else {
        NULL
      }
    })
    
    # Combine results
    do.call(rbind, growth_rates_list)
  }, error = function(e) {
    warning("Error calculating growth rates: ", e$message)
    NULL
  })

  # Cage effect analysis
  cage_analysis <- list()
  
  # Test for collinearity between cage and treatment
  cage_treatment_table <- table(analysis_df[[cage_column]], analysis_df[[treatment_column]])
  cage_analysis$collinearity_test <- tryCatch({
    stats::chisq.test(cage_treatment_table)
  }, error = function(e) {
    warning("Error in chi-square test: ", e$message)
    NULL
  })
  
  # Calculate cage-level effects
  cage_effects <- tryCatch({
    # Split data by cage and treatment
    split_data <- split(analysis_df, list(analysis_df[[cage_column]], analysis_df[[treatment_column]]))
    
    # Calculate statistics for each group
    cage_effects_list <- lapply(split_data, function(group_data) {
      if (nrow(group_data) > 0) {
        # Only create entries for non-empty groups
        data.frame(
          Cage = unique(group_data[[cage_column]]),
          Treatment = unique(group_data[[treatment_column]]),
          mean_volume = mean(group_data[[volume_column]], na.rm = TRUE),
          sd_volume = stats::sd(group_data[[volume_column]], na.rm = TRUE),
          n = nrow(group_data)
        )
      } else {
        # Skip empty groups
        NULL
      }
    })
    
    # Filter out NULL results and combine
    cage_effects_list <- cage_effects_list[!sapply(cage_effects_list, is.null)]
    if (length(cage_effects_list) > 0) {
      do.call(rbind, cage_effects_list)
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Error calculating cage effects: ", e$message)
    NULL
  })
  
  cage_analysis$effects <- cage_effects

  # Model selection
  model_selection <- list()
  
  # Fit different random effects specifications
  models <- list()
  
  # Base model (intercept only)
  models$intercept_only <- tryCatch({
    lme4::lmer(
      stats::as.formula(paste(volume_column, "~", time_column, "*", treatment_column, "+ (1|", id_column, ")")),
      data = analysis_df,
      control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                                check.nobs.vs.nRE = "ignore")
    )
  }, warning = function(w) {
    if (grepl("boundary", w$message)) {
      warning("Boundary (singular) fit detected in intercept-only model")
    }
    return(NULL)
  })
  
  # Random slope model
  models$random_slope <- tryCatch({
    lme4::lmer(
      stats::as.formula(paste(volume_column, "~", time_column, "*", treatment_column, "+ (", time_column, "|", id_column, ")")),
      data = analysis_df,
      control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                                check.nobs.vs.nRE = "ignore")
    )
  }, warning = function(w) {
    if (grepl("boundary", w$message)) {
      warning("Boundary (singular) fit detected in random slope model")
    }
    return(NULL)
  })
  
  # Remove any NULL models
  models <- models[!sapply(models, is.null)]
  
  if (length(models) > 0) {
    # Compare models using AIC and BIC
    model_selection$aic <- sapply(models, stats::AIC)
    model_selection$bic <- sapply(models, stats::BIC)
    
    # Select best model based on BIC (more conservative)
    best_model <- names(which.min(model_selection$bic))
    model <- models[[best_model]]
    
    if (verbose) {
      cat("Model selection results:\n")
      cat("AIC:", paste(names(model_selection$aic), "=", round(model_selection$aic, 2), collapse = ", "), "\n")
      cat("BIC:", paste(names(model_selection$bic), "=", round(model_selection$bic, 2), collapse = ", "), "\n")
      cat("Selected model:", best_model, "\n")
    }
  } else {
    warning("No valid models could be fitted. Using intercept-only model as fallback.")
    model <- lme4::lmer(
      stats::as.formula(paste(volume_column, "~", time_column, "*", treatment_column, "+ (1|", id_column, ")")),
      data = analysis_df,
      control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                                check.nobs.vs.nRE = "ignore")
    )
    model_selection$aic <- stats::AIC(model)
    model_selection$bic <- stats::BIC(model)
    model_selection$selected_model <- "intercept_only"
  }

  # Diagnostic plots
  diagnostics <- list()
  
  if (include_diagnostics) {
    # Residual plots
    diagnostics$residuals <- list(
      fitted = stats::fitted(model),
      residuals = stats::residuals(model),
      qq_plot = stats::qqnorm(stats::residuals(model))
    )
    
    # Random effects plots
    diagnostics$random_effects <- list(
      intercepts = lme4::ranef(model)[[id_column]],
      slopes = if (best_model == "random_slope") {
        lme4::ranef(model)[[id_column]][, 2]
      } else NULL
    )
    
    # Variance components
    diagnostics$variance_components <- lme4::VarCorr(model)
  }

  # Create a basic summary of the data
  data_summary <- tryCatch({
    # Split data by treatment and time
    split_data <- split(analysis_df, list(analysis_df[[treatment_column]], analysis_df[[time_column]]))
    
    # Calculate summary statistics for each group
    summary_list <- lapply(split_data, function(group_data) {
      if (nrow(group_data) > 0) {
        # Only create entries for non-empty groups
        data.frame(
          Treatment = unique(group_data[[treatment_column]]),
          Day = unique(group_data[[time_column]]),
          mean_volume = mean(group_data[[volume_column]], na.rm = TRUE),
          sd_volume = stats::sd(group_data[[volume_column]], na.rm = TRUE),
          n = nrow(group_data)
        )
      } else {
        # Skip empty groups
        NULL
      }
    })
    
    # Filter out NULL results and combine
    summary_list <- summary_list[!sapply(summary_list, is.null)]
    if (length(summary_list) > 0) {
      # Combine results and calculate standard error
      summary_df <- do.call(rbind, summary_list)
      summary_df$sem_volume <- summary_df$sd_volume / sqrt(summary_df$n)
      summary_df
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Error calculating data summary: ", e$message)
    NULL
  })

  # Only for auc model type or when additional AUC analysis is requested
  if (model_type == "auc" || (model_type == "lme4" && include_diagnostics)) {
    # Calculate AUC for each subject
    if (verbose) cat("Calculating AUC for each subject\n")
    
    # Calculate AUC using the trapezoidal rule
    calculate_auc <- function(time_values, volume_values) {
      # Sort by time
      sorted_indices <- order(time_values)
      time_values <- time_values[sorted_indices]
      volume_values <- volume_values[sorted_indices]
      
      # Need at least 2 points to calculate AUC
      if (length(time_values) < 2) {
        return(NA)
      }
      
      # Calculate AUC using the trapezoidal rule
      auc <- 0
      for (i in 2:length(time_values)) {
        time_diff <- time_values[i] - time_values[i-1]
        avg_height <- (volume_values[i] + volume_values[i-1]) / 2
        auc <- auc + (time_diff * avg_height)
      }
      
      return(auc)
    }
    
    # For each unique ID-Treatment-Cage combination, create a unique identifier
    # This ensures proper distinction of mice even when they share the same ID but are in different cages
    # First find all unique combinations
    unique_combinations <- unique(auc_df[, c(id_column, treatment_column, cage_column)])
    # Create a mapping of these combinations to sequential numbers
    unique_combinations$unique_id <- 1:nrow(unique_combinations)
    # Merge back with the original data to assign the correct unique ID to each row
    auc_df_with_id <- merge(auc_df, unique_combinations, by=c(id_column, treatment_column, cage_column))
    # Use this unique_id for processing
    composite_id <- paste(auc_df_with_id[[id_column]], auc_df_with_id[[treatment_column]], auc_df_with_id[[cage_column]], sep = "_")
    auc_data <- data.frame()
    
    # Get max experiment time to determine if extrapolation is needed
    max_experiment_time <- max(auc_df[[time_column]])
    
    for (unique_id in unique(composite_id)) {
      # Extract data for this unique ID
      id_parts <- strsplit(unique_id, "_")[[1]]
      actual_id <- id_parts[1]
      treatment <- id_parts[2]
      cage <- id_parts[3]
      
      if (length(id_parts) > 3) {
        # Handle the case where treatment has underscores (e.g., "Drug_A")
        treatment <- paste(id_parts[2:(length(id_parts)-1)], collapse = "_")
        cage <- id_parts[length(id_parts)]
      }
      
      subject_data <- auc_df_with_id[composite_id == unique_id, ]
      subject_data <- subject_data[order(subject_data[[time_column]]), ]
      
      # Calculate AUC using trapezoidal method
      auc_value <- calculate_auc(subject_data[[time_column]], subject_data[[volume_column]])
      
      # Check if this subject's data contains any extrapolated points
      has_extrapolated <- FALSE
      if ("Extrapolated" %in% colnames(subject_data)) {
        has_extrapolated <- any(subject_data$Extrapolated, na.rm = TRUE)
      }
      
      # Get the true last observation time (excluding extrapolated points)
      if (has_extrapolated && any(subject_data$Extrapolated, na.rm = TRUE)) {
        # If there are extrapolated points, get the max day from non-extrapolated points
        non_extrapolated_data <- subject_data[!subject_data$Extrapolated, ]
        true_last_observation <- max(non_extrapolated_data[[time_column]])
      } else {
        # If no extrapolated points, just use the max day
        true_last_observation <- max(subject_data[[time_column]])
      }
      
      # Make sure we have NumPoints data
      if (!"Extrapolated" %in% colnames(subject_data)) {
        n_points <- nrow(subject_data)
      } else {
        n_points <- nrow(subject_data[!subject_data$Extrapolated, ])
      }
      
      # Add to results
      auc_data <- rbind(auc_data, data.frame(
        ID = actual_id,
        Treatment = treatment,
        Cage = cage,
        Group = treatment, # Added Group column for compatibility with plot_auc
        AUC = auc_value,
        Last_Day = true_last_observation,
        First_Day = min(subject_data[[time_column]]),
        Extrapolated = has_extrapolated,
        NumPoints = n_points # Count only non-extrapolated points
      ))
    }
    
    # Calculate summary statistics
    auc_summary <- stats::aggregate(AUC ~ Treatment, data = auc_data, 
                                  FUN = function(x) c(Mean = mean(x), 
                                                    SD = stats::sd(x), 
                                                    N = length(x),
                                                    SEM = stats::sd(x)/sqrt(length(x))))
    auc_summary <- do.call(data.frame, auc_summary)
    
    # Create the AUC analysis list
    auc_analysis <- list(
      individual = auc_data,
      summary = auc_summary
    )
  }

  # Return the results for AUC model
  if (model_type == "auc") {
    # Create ANOVA model for AUC
    auc_model <- stats::aov(AUC ~ Treatment, data = auc_analysis$individual)
    anova_table <- stats::anova(auc_model)
    
    # Create pairwise comparisons for AUC using Welch's t-tests
    # This is more appropriate for AUC analysis where variances between groups may differ
    treatments <- unique(auc_analysis$individual$Treatment)
    pairwise_results <- list()
    pairwise_data <- list()
    
    # Generate all pairwise combinations
    pairs <- utils::combn(treatments, 2, simplify = FALSE)
    
    for(pair in pairs) {
      # Extract data for each treatment in the pair
      group1_data <- auc_analysis$individual$AUC[auc_analysis$individual$Treatment == pair[1]]
      group2_data <- auc_analysis$individual$AUC[auc_analysis$individual$Treatment == pair[2]]
      
      # Check if we have enough data for a t-test
      if(length(group1_data) < 2 || length(group2_data) < 2) {
        # Not enough data, create a placeholder result
        pairwise_results[[paste(pair[1], "-", pair[2])]] <- list(
          comparison = paste(pair[1], "-", pair[2]),
          mean_diff = ifelse(length(group1_data) > 0 && length(group2_data) > 0,
                            mean(group1_data) - mean(group2_data), NA),
          t_value = NA,
          df = NA,
          p_value = NA,
          ci_lower = NA,
          ci_upper = NA
        )
        
        # Store placeholder data
        pairwise_data[[paste(pair[1], "-", pair[2])]] <- list(
          group1 = pair[1],
          group2 = pair[2],
          data1 = group1_data,
          data2 = group2_data,
          result = list(
            statistic = NA,
            parameter = NA,
            p.value = NA,
            conf.int = c(NA, NA),
            estimate = NA
          )
        )
      } else {
        # We have enough data, perform Welch's t-test
        t_test_result <- stats::t.test(group1_data, group2_data, var.equal = FALSE)
        
        # Store results
        pairwise_results[[paste(pair[1], "-", pair[2])]] <- list(
          comparison = paste(pair[1], "-", pair[2]),
          mean_diff = mean(group1_data) - mean(group2_data),
          t_value = t_test_result$statistic,
          df = t_test_result$parameter,
          p_value = t_test_result$p.value,
          ci_lower = t_test_result$conf.int[1],
          ci_upper = t_test_result$conf.int[2]
        )
        
        # Store data for the posthoc object
        pairwise_data[[paste(pair[1], "-", pair[2])]] <- list(
          group1 = pair[1],
          group2 = pair[2],
          data1 = group1_data,
          data2 = group2_data,
          result = t_test_result
        )
      }
    }
    
    # Create data frame from pairwise results
    pairwise_df <- do.call(rbind, lapply(names(pairwise_results), function(comp) {
      res <- pairwise_results[[comp]]
      data.frame(
        comparison = res$comparison,
        estimate = res$mean_diff,
        t_value = res$t_value,
        df = res$df,
        p_value = res$p_value,
        ci_lower = res$ci_lower,
        ci_upper = res$ci_upper,
        stringsAsFactors = FALSE
      )
    }))
    
    # Apply Bonferroni correction for multiple comparisons
    pairwise_df$p_adjusted <- stats::p.adjust(pairwise_df$p_value, method = "bonferroni")
    
    # Handle reference group - ensure it exists in treatments
    if (!is.null(reference_group) && reference_group %in% treatments) {
      # Move reference group comparisons to the top
      ref_comparisons <- grep(paste0("^", reference_group, " -|^[^-]+ - ", reference_group, "$"), 
                           pairwise_df$comparison)
      if (length(ref_comparisons) > 0) {
        # Reorder rows to put reference group comparisons first
        pairwise_df <- rbind(
          pairwise_df[ref_comparisons, ],
          pairwise_df[-ref_comparisons, ]
        )
      }
    }
    
    # Extract treatment effects from AUC
    treatments <- unique(auc_analysis$individual$Treatment)
    treatment_effects <- data.frame(
      Treatment = treatments,
      Mean_AUC = sapply(treatments, function(t) {
        mean(auc_analysis$individual$AUC[auc_analysis$individual$Treatment == t])
      }),
      SD = sapply(treatments, function(t) {
        stats::sd(auc_analysis$individual$AUC[auc_analysis$individual$Treatment == t])
      }),
      N = sapply(treatments, function(t) {
        sum(auc_analysis$individual$Treatment == t)
      }),
      stringsAsFactors = FALSE
    )
    
    # Add reference indicator
    treatment_effects$Reference <- rep(FALSE, nrow(treatment_effects))
    if (!is.null(reference_group) && reference_group %in% treatments) {
      treatment_effects$Reference[treatment_effects$Treatment == reference_group] <- TRUE
    }
    
    # Create diagnostic plots
    if (include_diagnostics) {
      # Residual plots for AUC model
      diagnostics <- list(
        residuals = list(
          fitted = stats::fitted(auc_model),
          residuals = stats::residuals(auc_model),
          qq_plot = stats::qqnorm(stats::residuals(auc_model))
        )
      )
    } else {
      diagnostics <- NULL
    }
    
    # Create a descriptive summary
    analysis_summary <- list(
      analysis_type = "Area Under the Curve (AUC) Analysis",
      data_description = list(
        subjects = length(unique(paste(auc_df[[id_column]], auc_df[[treatment_column]], auc_df[[cage_column]], sep="_"))),
        treatment_groups = length(unique(auc_df[[treatment_column]])),
        time_points = length(unique(auc_df[[time_column]])),
        reference_group = reference_group
      ),
      methods = list(
        volume_transformation = transform,
        auc_calculation_method = auc_method,
        statistical_test = "One-way ANOVA on AUC values",
        posthoc_method = "Welch's t-tests with Bonferroni adjustment for multiple comparisons",
        individual_calculation = paste("AUC calculated using", auc_method, "method for each subject"),
        growth_rate_calculation = paste0(
          "Growth rates are calculated by fitting a linear regression model to log1p-transformed volume data over time for each subject. ",
          "The slope coefficient from this model represents the exponential growth rate. ",
          "A value of 0.1 indicates approximately 10% tumor volume increase per day. ",
          "Only subjects with 3 or more time points are included in growth rate calculations."
        )
      ),
      notes = c(
        if(transform != "none") paste("Volume data was", transform, "transformed prior to analysis") else "No transformation applied to volume data",
        "Composite IDs were created by combining subject ID, treatment group, and cage information to ensure correct AUC values",
        "Welch's t-tests are used for pairwise comparisons to account for potentially unequal variances between treatment groups"
      )
    )
    
    # Create posthoc object for compatibility with existing code
    posthoc <- list(
      method = "Welch's t-tests with Bonferroni adjustment",
      pairwise = pairwise_df,
      data = pairwise_data
    )
    
    # Return results for AUC model
    results <- list(
      model = if (return_model) auc_model else NULL,
      anova = anova_table,
      summary = analysis_summary,
      posthoc = posthoc,
      treatment_effects = treatment_effects,
      growth_rates = growth_rates,
      cage_analysis = cage_analysis,
      auc_analysis = auc_analysis,
      data_summary = data_summary,
      diagnostics = diagnostics
    )
    
    return(results)
  } else {
    # Create ANOVA table
    if (requireNamespace("car", quietly = TRUE)) {
      anova_table <- car::Anova(model, type = "III")
      anova_type <- "Type III ANOVA (car package)"
    } else {
      anova_table <- stats::anova(model)
      anova_type <- "ANOVA (stats package)"
      warning("Package 'car' not available. Using stats::anova instead of Type III tests.")
    }
    
    # Create pairwise comparisons
    if (requireNamespace("emmeans", quietly = TRUE)) {
      # Set up emmeans with reference group, averaging over time points
      lsmeans_obj <- emmeans::emmeans(model, specs = treatment_column, at = list(Day = mean(analysis_df$Day)))
      
      # Extract treatment effects
      emm_summary <- summary(lsmeans_obj)
      treatment_effects <- data.frame(
        Group = emm_summary[[treatment_column]],
        Adjusted_Mean = emm_summary$emmean,
        SE = emm_summary$SE,
        DF = emm_summary$df,
        Lower_CL = emm_summary$lower.CL,
        Upper_CL = emm_summary$upper.CL
      )

      # Add reference group note
      treatment_effects$Note <- ifelse(treatment_effects$Group == reference_group, "Reference group", "")

      # Reorder treatment effects to put reference group first
      ref_idx <- which(treatment_effects$Group == reference_group)
      if (ref_idx > 1) {
        treatment_effects <- rbind(
          treatment_effects[ref_idx, ],
          treatment_effects[-ref_idx, ]
        )
      }

      # Format numeric columns
      treatment_effects$Adjusted_Mean <- round(treatment_effects$Adjusted_Mean, 3)
      treatment_effects$SE <- round(treatment_effects$SE, 3)
      treatment_effects$Lower_CL <- round(treatment_effects$Lower_CL, 3)
      treatment_effects$Upper_CL <- round(treatment_effects$Upper_CL, 3)

      # Create contrasts with reference group
      contrasts <- list()
      other_groups <- setdiff(treatment_groups, reference_group)
      for (group in other_groups) {
        contrast_coef <- numeric(length(treatment_groups))
        ref_idx <- which(treatment_groups == reference_group)
        group_idx <- which(treatment_groups == group)
        contrast_coef[ref_idx] <- -1
        contrast_coef[group_idx] <- 1
        contrasts[[paste(group, "-", reference_group)]] <- contrast_coef
      }

      # Calculate pairwise comparisons
      pairwise_comp <- emmeans::contrast(lsmeans_obj, method = contrasts)
      
      posthoc_method <- "Estimated marginal means with pairwise contrasts"

    } else {
      pairwise_comp <- NULL
      treatment_effects <- NULL
      posthoc_method <- NA
      warning("Package 'emmeans' not available. Pairwise comparisons and treatment effects not calculated.")
    }
    
    # Create plots
    if (plots) {
      plot_data <- data_summary
      plots_list <- list(
        data_summary = plot_data,
        growth_rates = growth_rates,
        cage_effects = cage_effects,
        diagnostics = diagnostics
      )
      
      # Add adjusted means plot if available
      if (!is.null(treatment_effects)) {
        plots_list$adjusted_means <- treatment_effects
      }
    } else {
      plots_list <- NULL
    }
    
    # Create a descriptive summary
    analysis_summary <- list(
      analysis_type = "Linear Mixed Effects Model Analysis",
      data_description = list(
        subjects = length(unique(paste(analysis_df[[id_column]], analysis_df[[treatment_column]], analysis_df[[cage_column]], sep="_"))),
        treatment_groups = length(unique(analysis_df[[treatment_column]])),
        time_points = length(unique(analysis_df[[time_column]])),
        reference_group = reference_group
      ),
      model_specification = list(
        fixed_effects = paste(volume_column, "~", time_column, "*", treatment_column),
        random_effects = switch(random_effects_specification,
                               "intercept_only" = paste("(1|", id_column, ")"),
                               "slope" = paste("(", time_column, "|", id_column, ")"),
                               "none" = "None"),
        polynomial_degree = polynomial_degree,
        cage_effects = handle_cage_effects
      ),
      model_selection = list(
        criteria = "Model selected based on minimum BIC",
        selected_model = if(!is.null(model_selection$selected_model)) model_selection$selected_model else "Default model",
        models_compared = if(!is.null(model_selection$bic)) names(model_selection$bic) else "None"
      ),
      methods = list(
        volume_transformation = transform,
        anova_method = anova_type,
        posthoc_method = posthoc_method,
        growth_rate_calculation = paste0(
          "Growth rates are calculated by fitting a linear regression model to log1p-transformed volume data over time for each subject. ",
          "The slope coefficient from this model represents the exponential growth rate. ",
          "A value of 0.1 indicates approximately 10% tumor volume increase per day. ",
          "Only subjects with 3 or more time points are included in growth rate calculations."
        )
      ),
      notes = c(
        if(transform != "none") paste("Volume data was", transform, "transformed prior to analysis") else "No transformation applied to volume data",
        if(!is.null(cage_analysis$collinearity_test)) {
          paste("Cage-treatment collinearity test p-value:", format(cage_analysis$collinearity_test$p.value, digits=4))
        } else "No cage collinearity test performed"
      )
    )
    
    # Return the results
    results <- list(
      model = if (return_model) model else NULL,
      anova = anova_table,
      summary = analysis_summary,
      pairwise_comparisons = pairwise_comp,
      treatment_effects = treatment_effects,
      growth_rates = growth_rates,
      cage_analysis = cage_analysis,
      model_selection = model_selection,
      diagnostics = if (include_diagnostics) diagnostics else NULL,
      data_summary = data_summary,
      plots = plots_list
    )
    
    return(results)
  }
}