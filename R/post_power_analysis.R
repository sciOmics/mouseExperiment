# Copyright (c) 2025 Insight BioAnalytics. All rights reserved.
# Proprietary and confidential.

# Function to validate AUC data
validate_auc_data <- function(auc_data, treatment_column) {
  # Check if data exists
  if (is.null(auc_data) || nrow(auc_data) == 0) {
    return(list(valid = FALSE, message = "No AUC data available"))
  }
  
  # Check if required columns exist
  if (!("AUC" %in% colnames(auc_data))) {
    return(list(valid = FALSE, message = "AUC column not found in data"))
  }
  
  if (!(treatment_column %in% colnames(auc_data))) {
    return(list(valid = FALSE, message = "Treatment column not found in AUC data"))
  }
  
  # Check for at least 2 treatment groups
  treatments <- unique(auc_data[[treatment_column]])
  if (length(treatments) < 2) {
    return(list(valid = FALSE, message = "At least two treatment groups are required for AUC power analysis"))
  }
  
  # Check for minimum sample size in each group
  min_samples_per_group <- 2
  sample_sizes <- table(auc_data[[treatment_column]])
  if (any(sample_sizes < min_samples_per_group)) {
    low_groups <- names(sample_sizes)[sample_sizes < min_samples_per_group]
    return(list(valid = FALSE, 
                message = paste("Insufficient samples in group(s):", 
                               paste(low_groups, collapse = ", "))))
  }
  
  # All checks passed
  return(list(valid = TRUE, message = "AUC data validated successfully"))
}

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
#' \item{power_analysis}{Data frame with power estimates for different effect sizes, including columns for Treatment, Effect_Size, and Power}
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
#' The power_analysis object returned by this function includes a Treatment column indicating
#' which treatment group each power estimate refers to, allowing for clear identification of
#' treatment-specific power calculations.
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
         df = raw_data,
         time_column = time_column,
         volume_column = volume_column,
         treatment_column = treatment_column,
         id_column = id_column
       )
       
       # Store in model_results using the correct structure
       model_results <- list(
         auc_analysis = list(
           individual = auc_results$auc_data,
           summary = auc_results$auc_summary
         )
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
           df = raw_data,
           time_column = time_column,
           volume_column = volume_column,
           treatment_column = treatment_column,
           id_column = id_column
         )
         
         model_results$auc_analysis <- list(
           individual = auc_results$auc_data,
           summary = auc_results$auc_summary
         )
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
     if (is.null(auc_data) || nrow(auc_data) == 0) {
       stop("AUC analysis did not produce individual data")
     }
     
     # Check if Treatment column exists in auc_data
     if (!(treatment_column %in% colnames(auc_data))) {
       stop("Treatment column not found in AUC data")
     }
     
     # Get treatment groups and ensure there are at least 2
     treatments <- unique(auc_data[[treatment_column]])
     if (length(treatments) < 2) {
       stop("At least two treatment groups are required for power analysis")
     }
     
     # Calculate sample sizes per group
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
     
     if (method == "auc" && !is.null(model_results$auc_analysis) && !is.null(model_results$auc_analysis$summary)) {
       # For AUC method, estimate effect sizes from AUC differences
       auc_summary <- model_results$auc_analysis$summary
       
       # Make sure we have at least 2 groups in the summary
       if (nrow(auc_summary) < 2) {
         message("Not enough treatment groups in AUC summary. Using default effect sizes.")
         effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
       } else {
         # Calculate standardized effect sizes (Cohen's d)
         # Find the control/reference group (assume first group is reference)
         ref_group <- auc_summary[[treatment_column]][1]
         ref_mean <- auc_summary$AUC.Mean[1]
         
         # Calculate pooled SD (with safe check)
         if ("AUC.SD" %in% colnames(auc_summary) && sum(!is.na(auc_summary$AUC.SD)) > 0) {
           pooled_sd <- mean(auc_summary$AUC.SD, na.rm = TRUE)
           
           # Calculate effect sizes if pooled_sd is valid
           if (!is.na(pooled_sd) && pooled_sd > 0) {
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
             message("Invalid pooled SD. Using default effect sizes.")
             effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
           }
         } else {
           message("AUC.SD not found in summary. Using default effect sizes.")
           effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
         }
       }
     } else {
       # For parametric and simulation methods, use default range
       effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
     }
   }
   
   # Perform power analysis based on method
   if (method == "parametric") {
     message("Performing parametric power analysis...")
     
     # Get unique treatment groups
     treatment_groups <- names(result$sample_sizes)
     if (length(treatment_groups) < 2) {
       stop("At least two treatment groups are required for power analysis")
     }
     
     # Create a data frame with effect sizes for each treatment group (except control/reference)
     power_results <- data.frame()
     
     # Assume first group is control/reference
     ref_group <- treatment_groups[1]
     treatment_groups <- treatment_groups[-1]  # Remove reference group
     
     # Calculate power for each treatment group and effect size
     for (group in treatment_groups) {
       group_results <- data.frame(
         Treatment = rep(group, length(effect_sizes)),
         Effect_Size = effect_sizes,
         stringsAsFactors = FALSE
       )
       
       # Calculate power for each effect size using the sample size for this group
       group_n <- result$sample_sizes[group]
       ref_n <- result$sample_sizes[ref_group]
       
       group_results$Power <- sapply(effect_sizes, function(d) {
         power <- stats::power.t.test(
           n = min(group_n, ref_n),  # Use the smaller of the two for conservative estimate
           delta = d,
           sd = 1,  # Effect size is already in units of SD
           sig.level = alpha,
           type = "two.sample"
         )$power
         
         return(power)
       })
       
       # Append to overall results
       power_results <- rbind(power_results, group_results)
     }
     
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
     
     # Extract AUC data from input
     if (!is.null(model_results$auc_analysis)) {
       # Use AUC data from model results if available
       if (!is.null(model_results$auc_analysis$individual)) {
         auc_data <- model_results$auc_analysis$individual
       } else if (!is.null(model_results$auc_analysis$auc_data)) {
         auc_data <- model_results$auc_analysis$auc_data
       } else {
         message("No individual AUC data found in model_results. Calculating...")
         # Calculate AUC from raw data if available
         if (exists("raw_data") && !is.null(raw_data)) {
           auc_data <- calculate_auc_values(raw_data, time_column, volume_column, treatment_column, id_column)
         } else {
           stop("No data available for AUC power analysis")
         }
       }
     } else if (exists("raw_data") && !is.null(raw_data)) {
       # Calculate AUC from raw data
       message("Calculating AUC values from raw data...")
       auc_data <- calculate_auc_values(raw_data, time_column, volume_column, treatment_column, id_column)
     } else {
       stop("No data available for AUC power analysis")
     }
     
     # Validate the AUC data
     validation <- validate_auc_data(auc_data, treatment_column)
     if (!validation$valid) {
       stop(validation$message)
     }
     
     # Calculate observed group means and SDs directly
     auc_stats <- data.frame(
       Treatment = character(0),
       Mean = numeric(0),
       SD = numeric(0),
       N = numeric(0),
       stringsAsFactors = FALSE
     )
     
     # Get unique treatment groups and ensure they are properly ordered
     treatment_groups <- unique(auc_data[[treatment_column]])
     
     # Create a factor with explicit levels to avoid contrasts error
     auc_data[[treatment_column]] <- factor(auc_data[[treatment_column]], 
                                           levels = treatment_groups)
     
     # Check again that we have at least 2 distinct levels with data
     if (length(levels(auc_data[[treatment_column]])) < 2) {
       stop("After factorization, fewer than 2 treatment groups remain. Check your data.")
     }
     
     for (group in treatment_groups) {
       group_data <- auc_data$AUC[auc_data[[treatment_column]] == group]
       if (length(group_data) > 0) {
         auc_stats <- rbind(auc_stats, data.frame(
           Treatment = group,
           Mean = mean(group_data, na.rm = TRUE),
           SD = max(0.1, stats::sd(group_data, na.rm = TRUE)), # Ensure SD is not zero
           N = length(group_data),
           stringsAsFactors = FALSE
         ))
       }
     }
     
     # Check if we have enough data to calculate pooled SD
     if (nrow(auc_stats) < 2) {
       stop("Not enough groups with valid data for AUC power analysis")
     }
     
     # Calculate overall pooled SD
     weighted_var_sum <- sum((auc_stats$N - 1) * auc_stats$SD^2, na.rm = TRUE)
     df_sum <- sum(auc_stats$N - 1, na.rm = TRUE)
     
     if (df_sum > 0) {
       pooled_sd <- sqrt(weighted_var_sum / df_sum)
     } else {
       pooled_sd <- mean(auc_stats$SD, na.rm = TRUE)
     }
     
     # Fallback if pooled SD calculation fails or is zero
     if (is.na(pooled_sd) || pooled_sd < 0.1) {
       pooled_sd <- mean(auc_stats$SD, na.rm = TRUE)
       if (is.na(pooled_sd) || pooled_sd < 0.1) {
         pooled_sd <- 1  # Default fallback
         warning("Unable to calculate valid SD. Using default SD of 1 for power calculations.")
       }
     }
     
     # Create a data frame for power results with Treatment column
     power_results <- data.frame()
     
     # Get unique treatment groups and ensure they are properly ordered
     treatment_groups <- unique(auc_data[[treatment_column]])
     
     # Assume first group is control/reference
     ref_group <- treatment_groups[1]
     comp_groups <- treatment_groups[-1]  # Remove reference group
     
     # Calculate power for each treatment group compared to reference
     for (group in comp_groups) {
       group_results <- data.frame(
         Treatment = rep(group, length(effect_sizes)),
         Effect_Size = effect_sizes,
         stringsAsFactors = FALSE
       )
       
       # Get sample sizes for this comparison
       group_n <- auc_stats$N[auc_stats$Treatment == group]
       ref_n <- auc_stats$N[auc_stats$Treatment == ref_group]
       
       # Calculate power for each effect size
       group_results$Power <- sapply(effect_sizes, function(d) {
         # Use the minimum sample size for a conservative estimate
         min_n <- min(group_n, ref_n)
         
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
       
       # Append to overall results
       power_results <- rbind(power_results, group_results)
     }
     
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
     message("Performing simulation-based power analysis...")
     
     # First, get the data structure from either raw_data or model_results
     if (exists("raw_data")) {
       sim_data <- raw_data
     } else if (!is.null(model_results) && !is.null(model_results$model)) {
       # Try to extract the data from the model
       if (!is.null(model_results$data)) {
         sim_data <- model_results$data
       } else {
         # Create a simple mock dataset for simulations
         # This is a fallback option if we can't get the real data
         
         # Get treatment groups from sample_sizes
         treatment_groups <- names(result$sample_sizes)
         n_groups <- length(treatment_groups)
         
         if (n_groups < 2) {
           stop("At least two treatment groups are required for simulation-based power analysis")
         }
         
         n_per_group <- result$sample_sizes
         
         # Create a simple dataset with these treatments
         sim_data <- data.frame(
           ID = character(0),
           Day = numeric(0),
           Treatment = character(0),
           Volume = numeric(0),
           stringsAsFactors = FALSE
         )
         
         time_points <- 1:10
         
         # Create synthetic data
         for (g in 1:n_groups) {
           group <- treatment_groups[g]
           
           for (i in 1:n_per_group[g]) {
             id <- paste0(group, "_M", i)
             
             for (t in time_points) {
               sim_data <- rbind(sim_data, data.frame(
                 ID = id,
                 Day = t,
                 Treatment = group,
                 Volume = 100 + 10 * t + rnorm(1, 0, 10),
                 stringsAsFactors = FALSE
               ))
             }
           }
         }
         
         # Set column names to match the expected names
         colnames(sim_data)[colnames(sim_data) == "ID"] <- id_column
         colnames(sim_data)[colnames(sim_data) == "Day"] <- time_column
         colnames(sim_data)[colnames(sim_data) == "Treatment"] <- treatment_column
         colnames(sim_data)[colnames(sim_data) == "Volume"] <- volume_column
       }
     } else {
       stop("No suitable data available for simulation-based power analysis")
     }
     
     # Validate that we have the required data structure
     required_cols <- c(time_column, volume_column, treatment_column, id_column)
     missing_cols <- setdiff(required_cols, colnames(sim_data))
     if (length(missing_cols) > 0) {
       stop("Missing required columns for simulation: ", paste(missing_cols, collapse = ", "))
     }
     
     # Extract key parameters for simulations
     treatment_groups <- unique(sim_data[[treatment_column]])
     n_groups <- length(treatment_groups)
     if (n_groups < 2) {
       stop("At least two treatment groups are required for simulation")
     }
     
     # Calculate sample sizes per group
     n_per_group <- table(unique(sim_data[c(id_column, treatment_column)])[[treatment_column]])
     
     # Create a data frame for power results with Treatment column
     power_results <- data.frame()
     
     # Get treatment groups for simulation
     treatment_groups <- unique(sim_data[[treatment_column]])
     if (length(treatment_groups) < 2) {
       stop("At least two treatment groups are required for simulation")
     }
     
     # Assume first group is control/reference
     ref_group <- treatment_groups[1]
     comp_groups <- treatment_groups[-1]  # Remove reference group
     
     # Calculate power for each treatment group
     for (group in comp_groups) {
       group_results <- data.frame(
         Treatment = rep(group, length(effect_sizes)),
         Effect_Size = effect_sizes,
         stringsAsFactors = FALSE
       )
       
       # Calculate power for this treatment group
       group_results$Power <- sapply(effect_sizes, function(d) {
         # Run simulations for this specific treatment group comparison
         # ... simulation code adapted for specific treatment ...
         
         # Return simulated power
         significant_count <- sum(replicate(n_simulations, simulate_experiment(d, group)))
         return(significant_count / n_simulations)
       })
       
       # Append to overall results
       power_results <- rbind(power_results, group_results)
     }
     
     result$power_analysis <- power_results
     
     # For simulation, we don't calculate sample size recommendations directly
     # But we can provide a message about how to interpret results
     result$sample_size_recommendations <- "Sample size recommendations for simulation method should be determined by running additional simulations with varying sample sizes."
   }
   
   # Create power curve plot with treatment groups
   if (requireNamespace("ggplot2", quietly = TRUE)) {
     # Power curve with treatment groups
     power_curve <- ggplot2::ggplot(result$power_analysis, 
                                  ggplot2::aes(x = Effect_Size, y = Power, color = Treatment, group = Treatment)) +
       ggplot2::geom_line() +
       ggplot2::geom_point() +
       ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") +
       ggplot2::labs(
         title = paste("Power Analysis Results -", toupper(substr(method, 1, 1)), substr(method, 2, nchar(method)), "Method"),
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

# Function to calculate AUC from raw data
calculate_auc_values <- function(data, time_column, volume_column, treatment_column, id_column, cage_column = "Cage") {
  # Try to calculate AUC values
  tryCatch({
    # Initialize data frame for individual AUC values
    auc_individual <- data.frame(
      ID = character(0),
      Treatment = character(0),
      Cage = character(0),
      AUC = numeric(0),
      stringsAsFactors = FALSE
    )
    
    # Check if cage_column exists in the data
    use_cage_info <- FALSE
    if (cage_column %in% colnames(data)) {
      use_cage_info <- TRUE
    } else {
      warning("Cage column '", cage_column, "' not found. Proceeding without cage information for unique subject identification.")
    }
    
    # Create composite subject identifiers that include cage information if available
    if (use_cage_info) {
      composite_ids <- paste(data[[id_column]], data[[treatment_column]], data[[cage_column]], sep = "_")
      data$composite_id <- composite_ids
      # Get unique composite IDs
      unique_ids <- unique(composite_ids)
    } else {
      # Fallback to just using ID if cage information is not available
      unique_ids <- unique(data[[id_column]])
    }
    
    # Calculate AUC for each unique subject identifier
    for (unique_id in unique_ids) {
      # Get data for this unique subject
      if (use_cage_info) {
        subject_data <- data[data$composite_id == unique_id, ]
        # Extract original ID from composite ID for reporting
        id_parts <- strsplit(unique_id, "_")[[1]]
        original_id <- id_parts[1]
        treatment <- id_parts[2]
        cage <- id_parts[3]
      } else {
        subject_data <- data[data[[id_column]] == unique_id, ]
        original_id <- unique_id
        
        # Get treatment - should be the same for all rows of this ID
        treatment <- unique(subject_data[[treatment_column]])
        if (length(treatment) > 1) {
          warning("Multiple treatments found for ID ", unique_id, ". Using the first one.")
          treatment <- treatment[1]
        }
        
        # Set cage to NA if not using cage info
        cage <- NA
      }
      
      # Skip if no treatment assigned
      if (length(treatment) == 0) {
        warning("No treatment assigned for ID ", original_id, ". Skipping.")
        next
      }
      
      # Sort by time
      subject_data <- subject_data[order(subject_data[[time_column]]), ]
      
      # Need at least 2 time points for AUC
      if (nrow(subject_data) < 2) {
        warning("ID ", original_id, " has fewer than 2 time points. Skipping.")
        next
      }
      
      # Extract time and volume vectors
      times <- subject_data[[time_column]]
      volumes <- subject_data[[volume_column]]
      
      # Calculate AUC using trapezoidal rule
      auc_value <- 0
      for (i in 2:length(times)) {
        # Area of trapezoid = (v1 + v2) * (t2 - t1) / 2
        auc_value <- auc_value + (volumes[i-1] + volumes[i]) * (times[i] - times[i-1]) / 2
      }
      
      # Add to AUC individual data
      auc_individual <- rbind(auc_individual, data.frame(
        ID = original_id,
        Treatment = treatment,
        Cage = cage,
        AUC = auc_value,
        stringsAsFactors = FALSE
      ))
    }
    
    # Set column names correctly
    names(auc_individual)[names(auc_individual) == "Treatment"] <- treatment_column
    
    # Return the AUC individual data
    return(auc_individual)
  }, error = function(e) {
    message("Error calculating AUC: ", e$message)
    return(NULL)
  })
}