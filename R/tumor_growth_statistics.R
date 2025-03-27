# Copyright (c) 2025 Insight BioAnalytics. All rights reserved.
# Proprietary and confidential.

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
    if (verbose) cat("Extrapolating", extrapolation_points, "points for each subject\n")
    
    # Split data by subject
    df_split <- split(df, list(df[[id_column]], df[[treatment_column]]))
    
    # Function to extrapolate for one subject
    extrapolate_subject <- function(subject_data) {
      if (nrow(subject_data) < 2) return(subject_data)  # Need at least 2 points for extrapolation
      
      # Fit linear model to last 3 points (or all points if less than 3)
      n_points <- min(3, nrow(subject_data))
      last_points <- tail(subject_data, n_points)
      lm_fit <- stats::lm(paste(volume_column, "~", time_column), data = last_points)
      
      # Create new time points for extrapolation
      last_time <- max(subject_data[[time_column]])
      new_times <- seq(from = last_time + diff(range(subject_data[[time_column]])) / nrow(subject_data),
                      length.out = extrapolation_points,
                      by = diff(range(subject_data[[time_column]])) / nrow(subject_data))
      
      # Create data frame for prediction
      new_data <- data.frame(time = new_times)
      names(new_data) <- time_column
      
      # Predict new volumes and ensure they're non-negative
      new_volumes <- pmax(0, stats::predict(lm_fit, newdata = new_data))
      
      # Create new rows
      new_rows <- subject_data[1:extrapolation_points, ]
      new_rows[[time_column]] <- new_times
      new_rows[[volume_column]] <- new_volumes
      
      # Combine original and extrapolated data
      rbind(subject_data, new_rows)
    }
    
    # Apply extrapolation to each subject
    df_extrapolated <- do.call(rbind, lapply(df_split, extrapolate_subject))
    df <- df_extrapolated[order(df_extrapolated[[time_column]]), ]
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
    # Split data by treatment and ID
    split_data <- split(analysis_df, list(analysis_df[[treatment_column]], analysis_df[[id_column]]))
    
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
        
        # Return data frame with results
        data.frame(
          Treatment = unique(subject_data[[treatment_column]]),
          ID = unique(subject_data[[id_column]]),
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
    
    # Calculate AUC for each subject function
    calculate_auc <- function(subject_data) {
      # Sort by time
      subject_data <- subject_data[order(subject_data[[time_column]]), ]
      
      if (nrow(subject_data) < 2) {
        return(NA) # Need at least 2 points for AUC
      }
      
      # Calculate AUC using trapezoidal rule
      if (auc_method == "trapezoidal") {
        times <- subject_data[[time_column]]
        volumes <- subject_data[[volume_column]]
        
        # Calculate AUC using trapezoidal rule
        auc <- 0
        for (i in 2:length(times)) {
          dt <- times[i] - times[i-1]
          auc <- auc + dt * (volumes[i] + volumes[i-1]) / 2
        }
        return(auc)
      } else if (auc_method == "last_observation") {
        # Last observation carried forward
        # Simply take the latest time point and its volume
        latest <- subject_data[which.max(subject_data[[time_column]]), ]
        return(latest[[volume_column]])
      }
    }
    
    # Calculate AUC for each subject and collect metadata
    auc_list <- list()
    subjects <- unique(auc_df[[id_column]])
    
    # Calculate AUC for each subject with metadata
    for (i in seq_along(subjects)) {
      s <- subjects[i]
      subject_data <- auc_df[auc_df[[id_column]] == s, ]
      
      # Get the AUC value
      auc_result <- calculate_auc(subject_data)
      
      if (!is.na(auc_result)) {
        # Store AUC and metadata in list
        auc_list[[i]] <- data.frame(
          ID = s,
          Treatment = unique(subject_data[[treatment_column]]),
          Group = unique(subject_data[[treatment_column]]), # Add Group column for plot_auc compatibility
          AUC = as.numeric(auc_result),
          Last_Day = max(subject_data[[time_column]]),
          First_Day = min(subject_data[[time_column]]),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Combine into a single data frame
    individual_auc <- do.call(rbind, auc_list)
    
    # Calculate summary statistics
    auc_summary <- stats::aggregate(AUC ~ Treatment, data = individual_auc, 
                                  FUN = function(x) c(Mean = mean(x), 
                                                    SD = stats::sd(x), 
                                                    N = length(x),
                                                    SEM = stats::sd(x)/sqrt(length(x))))
    auc_summary <- do.call(data.frame, auc_summary)
    
    # Create the AUC analysis list
    auc_analysis <- list(
      individual = individual_auc,
      summary = auc_summary
    )
  }

  # For the AUC model type
  if (model_type == "auc") {
    # Create ANOVA model for AUC
    auc_model <- stats::aov(AUC ~ Treatment, data = auc_analysis$individual)
    anova_table <- stats::anova(auc_model)
    
    # Create pairwise comparisons for AUC
    if (requireNamespace("emmeans", quietly = TRUE)) {
      # Create "posthoc" for backward compatibility
      posthoc <- tryCatch({
        pairwise <- emmeans::emmeans(auc_model, pairwise ~ Treatment, adjust = "bonferroni")
        pairwise
      }, error = function(e) {
        warning("Error creating posthoc comparisons: ", e$message)
        NULL
      })
      
      # Extract treatment effects from AUC
      treatment_effects <- data.frame(
        Group = auc_summary$Treatment,
        Adjusted_Mean = auc_summary$AUC.Mean,
        SE = auc_summary$AUC.SEM,
        Lower_CL = auc_summary$AUC.Mean - 1.96 * auc_summary$AUC.SEM,
        Upper_CL = auc_summary$AUC.Mean + 1.96 * auc_summary$AUC.SEM,
        stringsAsFactors = FALSE
      )
      
      # Add reference group note
      treatment_effects$Note <- ifelse(treatment_effects$Group == reference_group, "Reference group", "")
      
      # Reorder to put reference group first
      ref_idx <- which(treatment_effects$Group == reference_group)
      if (ref_idx > 1) {
        treatment_effects <- rbind(
          treatment_effects[ref_idx, ],
          treatment_effects[-ref_idx, ]
        )
      }
    } else {
      posthoc <- NULL
      treatment_effects <- NULL
      warning("Package 'emmeans' not available. Pairwise comparisons not calculated.")
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
    
    # Return results for AUC model
    results <- list(
      model = if (return_model) auc_model else NULL,
      anova = anova_table,
      summary = summary(auc_model),
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
    } else {
      anova_table <- stats::anova(model)
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

    } else {
      pairwise_comp <- NULL
      treatment_effects <- NULL
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

    # Return the results
    results <- list(
      model = if (return_model) model else NULL,
      anova = anova_table,
      summary = summary(model),
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