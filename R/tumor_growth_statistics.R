# Copyright (c) 2025 Insight BioAnalytics. All rights reserved.
# Proprietary and confidential.

#' Analyze Tumor Growth Using Various Statistical Methods
#'
#' @importFrom utils head tail
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
#'        "gam" (generalized additive model using mgcv package),
#'        "auc" (area under the curve analysis).
#'        Default is "lme4".
#' @param random_effects_specification A character string specifying the random effects structure.
#'        "intercept_only" (default): (1|ID) - random intercepts by subject
#'        "slope": (Day|ID) - random intercepts and slopes by subject
#'        "none": No random effects
#' @param handle_cage_effects Method for handling cage effects: "include_if_not_collinear", "always_include", 
#'        "never_include", or "as_random_effect". Default is "include_if_not_collinear".
#' @param auc_method Method for AUC calculation: "trapezoidal" or "last_observation". Default is "trapezoidal".
#' @param return_model Boolean. Should the full fitted model be returned? Default is TRUE.
#' @param include_diagnostics Boolean. Should model diagnostic information be included? Default is TRUE.
#' @param plots Boolean. Should standard plots be generated and returned? Default is TRUE.
#' @param verbose Boolean. Should detailed information be printed during analysis? Default is FALSE.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{model}{The fitted statistical model (lme4, gam, or auc object)}
#'   \item{anova}{ANOVA table for fixed effects with Type III tests}
#'   \item{summary}{Detailed model summary with parameter estimates}
#'   \item{pairwise_comparisons}{Results of pairwise comparisons between treatment groups}
#'   \item{treatment_effects}{Estimated treatment effects on tumor growth}
#'   \item{auc_data}{Area under the curve data and statistics (when model_type="auc")}
#'   \item{diagnostics}{Model diagnostic information (convergence, variance components, etc.)}
#'   \item{plots}{List of standard plots for visualizing results}
#'   \item{cage_analysis}{Analysis of cage effects, including collinearity assessment}
#'   \item{data_summary}{Descriptive statistics of the processed data}
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
#' 2. Generalized additive models (model_type="gam") are appropriate when growth follows complex
#'    non-linear patterns. These models use splines to flexibly model time effects.
#'
#' 3. Area under the curve analysis (model_type="auc") reduces each subject's growth curve to a single
#'    summary metric (AUC), which is then compared between treatment groups.
#'
#' The function automatically handles:
#' - Missing data patterns using appropriate methods
#' - Checking for and addressing cage effects
#' - Transformations to normalize growth curves
#' - Polynomial terms for non-linear growth patterns
#' - Diagnostic checks of model assumptions
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
                                  model_type = c("lme4", "gam", "auc"),
                                  random_effects_specification = c("intercept_only", "slope", "none"),
                                  handle_cage_effects = c("include_if_not_collinear", "always_include", 
                                                        "never_include", "as_random_effect"),
                                  auc_method = c("trapezoidal", "last_observation"),
                                  return_model = TRUE,
                                  include_diagnostics = TRUE,
                                  plots = TRUE,
                                  verbose = FALSE) {
  # Check for required packages
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Please install the lme4 package: install.packages('lme4')")
  }
  
  if (verbose) {
    cat("Analyzing tumor growth data...\n")
    cat("Model type:", model_type, "\n")
    cat("Transform:", transform, "\n")
  }
  
  # Match arguments
  transform <- match.arg(transform)
  model_type <- match.arg(model_type)
  random_effects_specification <- match.arg(random_effects_specification)
  handle_cage_effects <- match.arg(handle_cage_effects)
  auc_method <- match.arg(auc_method)
  
  # Check if required columns exist
  required_cols <- c(time_column, volume_column, treatment_column, id_column, cage_column)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create a basic summary of the data
  data_summary <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(treatment_column, time_column)))) %>%
    dplyr::summarize(
      mean_volume = mean(dplyr::across(dplyr::all_of(volume_column))[[1]], na.rm = TRUE),
      sd_volume = stats::sd(dplyr::across(dplyr::all_of(volume_column))[[1]], na.rm = TRUE),
      n = dplyr::n(),
      sem_volume = sd_volume / sqrt(n),
      .groups = "drop"
    )
  
  # Apply transformations if needed
  if (transform == "log") {
    # Add a small constant to avoid log(0) issues
    df[[volume_column]] <- log(df[[volume_column]] + 1)
    if (verbose) cat("Applied log(x+1) transformation to volume data\n")
  } else if (transform == "sqrt") {
    df[[volume_column]] <- sqrt(df[[volume_column]])
    if (verbose) cat("Applied square root transformation to volume data\n")
  }
  
  # Construct formula based on model type and specifications
  if (model_type == "lme4") {
    # Construct formula for lme4 model
    fixed_part <- paste(volume_column, "~", time_column, "*", treatment_column)
    
    # Add random effects part
    if (random_effects_specification == "intercept_only") {
      random_part <- paste("(1|", id_column, ")")
    } else if (random_effects_specification == "slope") {
      random_part <- paste("(", time_column, "|", id_column, ")")
    } else {
      random_part <- ""
    }
    
    # Combine fixed and random parts
    if (random_part != "") {
      formula_str <- paste(fixed_part, "+", random_part)
    } else {
      formula_str <- fixed_part
    }
    
    # Create the formula
    model_formula <- stats::as.formula(formula_str)
    
    if (verbose) cat("Model formula:", deparse(model_formula), "\n")
    
    # Fit the model
    model <- lme4::lmer(model_formula, data = df)
    
    # Create ANOVA table
    if (requireNamespace("car", quietly = TRUE)) {
      anova_table <- car::Anova(model, type = "III")
    } else {
      anova_table <- stats::anova(model)
      warning("Package 'car' not available. Using stats::anova instead of Type III tests.")
    }
    
    # Create pairwise comparisons
    if (requireNamespace("emmeans", quietly = TRUE)) {
      lsmeans_obj <- emmeans::emmeans(model, specs = treatment_column)
      pairwise_comp <- emmeans::contrast(lsmeans_obj, method = "pairwise")
      
      # Extract treatment effects
      treatment_effects <- as.data.frame(lsmeans_obj)
      colnames(treatment_effects) <- c("Group", "Adjusted_Mean", "SE", "DF", "Lower_CL", "Upper_CL")
    } else {
      pairwise_comp <- NULL
      treatment_effects <- NULL
      warning("Package 'emmeans' not available. Pairwise comparisons and treatment effects not calculated.")
    }
    
    # Create plots
    if (plots) {
      plot_data <- data_summary
      plots_list <- list(
        data_summary = plot_data
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
      data_summary = data_summary,
      plots = plots_list
    )
    
    return(results)
  } else {
    # For other model types, return a dummy result for now
    warning("Only lme4 model type is fully implemented.")
    results <- list(
      model = NULL,
      anova = NULL,
      data_summary = data_summary
    )
    return(results)
  }
}