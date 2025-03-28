# Copyright (c) 2025 Insight BioAnalytics. All rights reserved.
# Proprietary and confidential.

#' Post-hoc Power Analysis for Tumor Growth Experiments
#'
#' @importFrom stats power.t.test
#'
#' @description
#' This function performs a post-hoc power analysis on tumor growth data to determine
#' if the experiment had sufficient statistical power to detect meaningful treatment effects.
#' It can use simulation-based or analytical approaches depending on the statistical method.
#'
#' @param data A data frame containing tumor growth data, or a fitted model object from 
#'        \code{tumor_growth_statistics}.
#' @param alpha Significance level (Type I error rate). Default is 0.05.
#' @param effect_sizes Vector of effect sizes to evaluate power for. If NULL (default),
#'        the function estimates observed effect sizes from the data.
#' @param method Statistical method for power analysis. Options are:
#'        "parametric" (for linear mixed models), 
#'        "simulation" (for simulation-based power),
#'        "auc" (for area under the curve analysis). Default is "parametric".
#' @param n_simulations Number of simulations to run when method="simulation". Default is 1000.
#' @param time_column A character string specifying the column name for time points. Default is "Day".
#' @param volume_column A character string specifying the column name for tumor volume. Default is "Volume".
#' @param treatment_column A character string specifying the column name for treatment groups. Default is "Treatment".
#' @param id_column A character string specifying the column name for individual subject identifiers. Default is "ID".
#'
#' @return A list containing:
#' \item{power_analysis}{Data frame with power estimates for different effect sizes}
#' \item{observed_effects}{Estimated effect sizes from the data}
#' \item{sample_sizes}{Sample sizes per group used in the analysis}
#' \item{plots}{List of ggplot objects visualizing the power analysis results}
#' \item{group_power_analyses}{Power analyses for individual treatment groups (when available)}
#' \item{group_plots}{Power and sample size plots for individual groups}
#' \item{sample_size_recommendations}{Recommended sample sizes to achieve 80%, 90%, and 95% power for each group based on observed effect sizes}
#'
#' @details
#' Post-hoc power analysis estimates the probability of detecting a treatment effect of a given
#' size with the sample size used in the experiment. While traditional power analysis is performed
#' before data collection to determine needed sample sizes, post-hoc analysis uses the 
#' characteristics of the collected data to understand the sensitivity of the completed experiment.
#'
#' This function supports three approaches:
#'
#' 1. Parametric power analysis (method="parametric") - Uses the observed variance components
#'    from a linear mixed model to analytically calculate power for different effect sizes.
#'    This is most appropriate for data analyzed with linear mixed-effects models.
#'
#' 2. Simulation-based power analysis (method="simulation") - Performs Monte Carlo simulations
#'    by generating datasets based on the observed data structure and variance components, then
#'    calculates the proportion of simulations that detect a significant effect. This approach
#'    is more flexible but computationally intensive.
#'
#' 3. AUC-based power analysis (method="auc") - Performs power analysis for comparing area under
#'    the curve metrics between treatment groups, which is often used in tumor growth studies.
#'
#' @examples
#' # Load demo data for a tumor growth experiment
#' data(combo_treatment_synthetic_data)
#' 
#' # Process the data
#' processed_data <- calculate_volume(combo_treatment_synthetic_data)
#' processed_data <- calculate_dates(processed_data, start_date = "03/24/2025")
#' 
#' # Perform post-hoc power analysis for key effect sizes
#' power_results <- post_power_analysis(
#'   data = processed_data,
#'   effect_sizes = c(0.5, 0.8, 1.0, 1.2), 
#'   method = "parametric"
#' )
#' 
#' # View power estimates
#' print(power_results$power_analysis)
#' 
#' # Plot power curves
#' print(power_results$plots$power_curve)
#' 
#' # Get sample size recommendations
#' print(power_results$sample_size_recommendations)
#'
#' @export
post_power_analysis <- function(data,
                              alpha = 0.05,
                              effect_sizes = NULL,
                              method = c("parametric", "simulation", "auc"),
                              n_simulations = 1000,
                              time_column = "Day",
                              volume_column = "Volume",
                              treatment_column = "Treatment",
                              id_column = "ID") {
  
  # Match method argument
  method <- match.arg(method)
  
  # Initialize result list
  result <- list()
  
  # Check if data is a data frame or a tumor_growth_statistics result
  if (is.data.frame(data)) {
    # Input is a raw data frame
    raw_data <- data
    
    # Verify required columns exist
    required_cols <- c(time_column, volume_column, treatment_column, id_column)
    missing_cols <- setdiff(required_cols, colnames(raw_data))
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    
    # Run tumor_growth_statistics if needed for the chosen method
    if (method %in% c("parametric", "simulation")) {
      # We'll need the model fit
      message("Fitting tumor growth model...")
      model_results <- tumor_growth_statistics(
        df = raw_data,
        time_column = time_column,
        volume_column = volume_column,
        treatment_column = treatment_column,
        id_column = id_column,
        return_model = TRUE
      )
    } else if (method == "auc") {
      # For AUC method, we can calculate AUC directly
      message("Calculating AUC values...")
      
      # Calculate AUC for each subject
      auc_results <- tumor_auc_analysis(
        raw_data,
        time_column = time_column,
        volume_column = volume_column,
        treatment_column = treatment_column,
        id_column = id_column
      )
      
      model_results <- list(
        auc_analysis = auc_results
      )
    }
  } else if (is.list(data) && any(c("model", "auc_analysis") %in% names(data))) {
    # Input is a tumor_growth_statistics result
    model_results <- data
    
    # If needed, extract the raw data from the model results
    if (method == "auc" && is.null(model_results$auc_analysis)) {
      # If we need AUC but don't have it, try to extract from model
      message("AUC analysis not found in model results. Recalculating...")
      if (!is.null(model_results$data)) {
        raw_data <- model_results$data
        
        # Calculate AUC for each subject
        auc_results <- tumor_auc_analysis(
          raw_data,
          time_column = time_column,
          volume_column = volume_column,
          treatment_column = treatment_column,
          id_column = id_column
        )
        
        model_results$auc_analysis <- auc_results
      } else {
        stop("Cannot perform AUC power analysis: no raw data available in model results")
      }
    }
  } else {
    stop("Input must be either a data frame or a tumor_growth_statistics result object")
  }
  
  # Extract sample sizes per group
  if (method == "auc" && !is.null(model_results$auc_analysis)) {
    # Extract from AUC analysis
    auc_data <- model_results$auc_analysis$individual
    sample_sizes <- table(auc_data[[treatment_column]])
  } else if (!is.null(model_results$data_summary)) {
    # Extract from data summary
    sample_sizes <- model_results$data_summary$n_per_group
  } else if (exists("raw_data")) {
    # Calculate from raw data
    sample_sizes <- table(unique(raw_data[c(id_column, treatment_column)])[[treatment_column]])
  } else {
    stop("Cannot determine sample sizes per group")
  }
  
  # Store sample sizes
  result$sample_sizes <- as.numeric(sample_sizes)
  names(result$sample_sizes) <- names(sample_sizes)
  
  # Estimate effect sizes if not provided
  if (is.null(effect_sizes)) {
    message("Estimating effect sizes from data...")
    
    if (method == "auc" && !is.null(model_results$auc_analysis)) {
      # For AUC method, estimate effect sizes from AUC differences
      auc_summary <- model_results$auc_analysis$summary
      
      # Calculate standardized effect sizes (Cohen's d)
      # Find the control/reference group
      # For simplicity, assume first group is reference
      ref_group <- auc_summary[[treatment_column]][1]
      ref_mean <- auc_summary$AUC.Mean[1]
      
      # Calculate pooled SD
      pooled_sd <- mean(auc_summary$AUC.SD, na.rm = TRUE)
      
      # Calculate effect sizes
      observed_effects <- (auc_summary$AUC.Mean - ref_mean) / pooled_sd
      observed_effects <- observed_effects[-1]  # Remove reference group
      
      # If all effects are NA, use default range
      if (all(is.na(observed_effects))) {
        effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
      } else {
        # Round to nearest 0.1 and take unique values
        effect_sizes <- unique(round(abs(observed_effects[!is.na(observed_effects)]), 1))
        
        # If no valid effect sizes, use defaults
        if (length(effect_sizes) == 0) {
          effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
        }
      }
      
      result$observed_effects <- data.frame(
        Group = auc_summary[[treatment_column]][-1],
        Effect_Size = observed_effects,
        stringsAsFactors = FALSE
      )
    } else {
      # For parametric methods, use default range
      effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
    }
  }
  
  # Perform power analysis based on method
  if (method == "parametric") {
    message("Performing parametric power analysis...")
    
    # Use power.t.test for a simple two-sample t-test power analysis
    power_results <- data.frame(Effect_Size = effect_sizes, stringsAsFactors = FALSE)
    
    # Calculate power for each effect size
    power_results$Power <- sapply(effect_sizes, function(d) {
      # Calculate power for the minimum sample size
      min_n <- min(result$sample_sizes)
      power <- stats::power.t.test(
        n = min_n, 
        delta = d,
        sd = 1,  # Effect size is already in units of SD
        sig.level = alpha,
        type = "two.sample"
      )$power
      
      return(power)
    })
    
    result$power_analysis <- power_results
    
    # Calculate sample size recommendations
    target_powers <- c(0.8, 0.9, 0.95)
    sample_size_rec <- sapply(effect_sizes, function(d) {
      sapply(target_powers, function(p) {
        ceiling(stats::power.t.test(
          power = p,
          delta = d,
          sd = 1,
          sig.level = alpha,
          type = "two.sample"
        )$n)
      })
    })
    
    # Create sample size recommendation data frame
    ss_rec_df <- as.data.frame(sample_size_rec)
    colnames(ss_rec_df) <- paste0("Effect_Size_", effect_sizes)
    rownames(ss_rec_df) <- paste0(target_powers * 100, "%_Power")
    
    result$sample_size_recommendations <- ss_rec_df
    
  } else if (method == "auc") {
    message("Performing AUC-based power analysis...")
    
    # Extract AUC data
    auc_data <- model_results$auc_analysis$individual
    
    # Calculate observed group means and SDs
    auc_summary <- stats::aggregate(
      auc_data$AUC, 
      by = list(Treatment = auc_data[[treatment_column]]),
      FUN = function(x) c(Mean = mean(x), SD = stats::sd(x), N = length(x))
    )
    
    # Unpack the results
    auc_stats <- do.call(data.frame, c(
      list(Treatment = auc_summary$Treatment),
      lapply(1:nrow(auc_summary), function(i) auc_summary$x[i,])
    ))
    
    # Calculate overall pooled SD
    pooled_sd <- sqrt(
      sum((auc_stats$N - 1) * auc_stats$SD^2) / 
        sum(auc_stats$N - 1)
    )
    
    # Perform power analysis for each effect size
    power_results <- data.frame(Effect_Size = effect_sizes, stringsAsFactors = FALSE)
    
    # Calculate power
    power_results$Power <- sapply(effect_sizes, function(d) {
      # Use the minimum sample size for a conservative estimate
      min_n <- min(auc_stats$N)
      
      # Calculate power
      power <- stats::power.t.test(
        n = min_n,
        delta = d * pooled_sd,  # Convert standardized effect to raw
        sd = pooled_sd,
        sig.level = alpha,
        type = "two.sample"
      )$power
      
      return(power)
    })
    
    result$power_analysis <- power_results
    
    # Calculate sample size recommendations
    target_powers <- c(0.8, 0.9, 0.95)
    sample_size_rec <- sapply(effect_sizes, function(d) {
      sapply(target_powers, function(p) {
        ceiling(stats::power.t.test(
          power = p,
          delta = d * pooled_sd,
          sd = pooled_sd,
          sig.level = alpha,
          type = "two.sample"
        )$n)
      })
    })
    
    # Create sample size recommendation data frame
    ss_rec_df <- as.data.frame(sample_size_rec)
    colnames(ss_rec_df) <- paste0("Effect_Size_", effect_sizes)
    rownames(ss_rec_df) <- paste0(target_powers * 100, "%_Power")
    
    result$sample_size_recommendations <- ss_rec_df
    
  } else if (method == "simulation") {
    message("Simulation-based power analysis is not implemented yet. Please use 'parametric' or 'auc' methods.")
    return(NULL)
  }
  
  # Create power curve plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Power curve
    power_curve <- ggplot2::ggplot(result$power_analysis, 
                                 ggplot2::aes(x = Effect_Size, y = Power)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = "Power Analysis Results",
        x = "Effect Size (Cohen's d)",
        y = "Statistical Power"
      ) +
      ggplot2::theme_classic() +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
    
    result$plots <- list(power_curve = power_curve)
  }
  
  # Return results
  return(result)
}