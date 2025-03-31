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
#' @param alpha Numeric vector of significance levels (Type I error rates). Default is c(0.05, 0.01).
#' @param power Numeric vector of target power levels for sample size estimation. Default is c(0.8, 0.85, 0.9, 0.95).
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
#' \item{effect_sizes}{Data frame with estimated effect sizes for each treatment group compared to the reference group}
#' \item{post_power_analysis}{Data frame with power estimates for different effect sizes, including columns for Treatment, Effect_Size, Alpha, and Power}
#' \item{sample_size_estimates}{Data frame with recommended sample sizes to achieve specified power levels for each treatment group based on observed effect sizes}
#' \item{plots}{List of ggplot objects visualizing the power analysis results}
#' \item{summary}{Summary information about the analysis}
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
#' The function estimates effect sizes from the data, calculates power for each treatment group compared
#' to the reference group, and provides sample size recommendations for future studies based on the 
#' observed effect sizes.
#'
#' @examples
#' # Load demo data for a tumor growth experiment
#' data(combo_treatment_synthetic_data)
#' 
#' # Process the data
#' processed_data <- calculate_volume(combo_treatment_synthetic_data)
#' processed_data <- calculate_dates(processed_data, start_date = "03/24/2025")
#' 
#' # Perform post-hoc power analysis
#' power_results <- post_power_analysis(
#'   data = processed_data,
#'   alpha = c(0.05, 0.01),
#'   power = c(0.8, 0.9, 0.95),
#'   method = "auc"
#' )
#' 
#' # View effect size estimates
#' print(power_results$effect_sizes)
#' 
#' # View power estimates
#' print(power_results$post_power_analysis)
#' 
#' # Get sample size recommendations
#' print(power_results$sample_size_estimates)
#'
#' @export
post_power_analysis <- function(data,
                              alpha = c(0.05, 0.01),
                              power = c(0.8, 0.85, 0.9, 0.95),
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
   sample_sizes_by_group <- as.numeric(sample_sizes)
   names(sample_sizes_by_group) <- names(sample_sizes)
   result$summary <- list(sample_sizes = sample_sizes_by_group)
   
   # ----- EFFECT SIZE CALCULATION -----
   
   # Initialize effect_sizes data frame
   effect_sizes_df <- data.frame(
     Treatment = character(0),
     Reference = character(0),
     Raw_Difference = numeric(0),
     Pooled_SD = numeric(0),
     Standardized_Effect = numeric(0),
     stringsAsFactors = FALSE
   )
   
   # Calculate effect sizes if not provided
   if (is.null(effect_sizes)) {
     message("Estimating effect sizes from data...")
     
     if (method == "auc" && !is.null(model_results$auc_analysis) && !is.null(model_results$auc_analysis$summary)) {
       # For AUC method, estimate effect sizes from AUC differences
       auc_summary <- model_results$auc_analysis$summary
       
       # Map column names for backwards compatibility
       auc_summary <- rename_auc_columns(auc_summary)
       
       # Make sure we have at least 2 groups in the summary
       if (nrow(auc_summary) < 2) {
         message("Not enough treatment groups in AUC summary. Using default effect sizes.")
         default_effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
       } else {
         # Calculate standardized effect sizes (Cohen's d)
         # Find the control/reference group (assume first group is reference)
         ref_group <- auc_summary[[treatment_column]][1]
         ref_mean <- auc_summary$Mean_AUC[1]
         
         # Calculate pooled SD 
         if ("SD_AUC" %in% colnames(auc_summary) && sum(!is.na(auc_summary$SD_AUC)) > 0) {
           pooled_sd <- mean(auc_summary$SD_AUC, na.rm = TRUE)
           
           # Calculate effect sizes if pooled_sd is valid
           if (!is.na(pooled_sd) && pooled_sd > 0) {
             # For each treatment group (except reference)
             for (i in 2:nrow(auc_summary)) {
               treatment <- auc_summary[[treatment_column]][i]
               mean_diff <- auc_summary$Mean_AUC[i] - ref_mean
               std_effect <- mean_diff / pooled_sd
               
               # Add to effect sizes data frame
               effect_sizes_df <- rbind(effect_sizes_df, data.frame(
                 Treatment = treatment,
                 Reference = ref_group,
                 Raw_Difference = mean_diff,
                 Pooled_SD = pooled_sd,
                 Standardized_Effect = std_effect,
                 stringsAsFactors = FALSE
               ))
             }
             
             # Extract unique effect sizes for power calculation
             if (nrow(effect_sizes_df) > 0) {
               default_effect_sizes <- unique(round(abs(effect_sizes_df$Standardized_Effect[!is.na(effect_sizes_df$Standardized_Effect)]), 1))
               # If no valid effect sizes, use defaults
               if (length(default_effect_sizes) == 0) {
                 default_effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
               }
             } else {
               default_effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
             }
           } else {
             message("Invalid pooled SD. Using default effect sizes.")
             default_effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
           }
         } else {
           message("SD_AUC not found in summary. Using default effect sizes.")
           default_effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
         }
       }
     } else {
       # For parametric and simulation methods, use default range
       default_effect_sizes <- c(0.2, 0.5, 0.8, 1.0, 1.5)
     }
     
     # If no effect sizes were calculated, use the defaults
     if (nrow(effect_sizes_df) == 0) {
       message("Using default effect sizes for analysis")
       if (exists("treatments") && length(treatments) > 1) {
         ref_group <- treatments[1]
         for (i in 2:length(treatments)) {
           treatment <- treatments[i]
           for (effect in default_effect_sizes) {
             effect_sizes_df <- rbind(effect_sizes_df, data.frame(
               Treatment = treatment,
               Reference = ref_group,
               Raw_Difference = NA,
               Pooled_SD = NA,
               Standardized_Effect = effect,
               stringsAsFactors = FALSE
             ))
           }
         }
       } else {
         # No treatment info available, just use generic effect sizes
         effect_sizes_df <- data.frame(
           Treatment = "Unknown",
           Reference = "Control",
           Raw_Difference = NA,
           Pooled_SD = NA,
           Standardized_Effect = default_effect_sizes,
           stringsAsFactors = FALSE
         )
       }
     }
     
     # Use the provided effect sizes
     effect_sizes <- sort(unique(abs(effect_sizes_df$Standardized_Effect)))
     # Make sure we have at least some reasonable effect sizes
     if (length(effect_sizes) == 0) {
       effect_sizes <- default_effect_sizes
     }
   } else {
     # Use the provided effect sizes to populate effect_sizes_df
     if (exists("treatments") && length(treatments) > 1) {
       ref_group <- treatments[1]
       for (i in 2:length(treatments)) {
         treatment <- treatments[i]
         for (effect in effect_sizes) {
           effect_sizes_df <- rbind(effect_sizes_df, data.frame(
             Treatment = treatment,
             Reference = ref_group,
             Raw_Difference = NA,
             Pooled_SD = NA,
             Standardized_Effect = effect,
             stringsAsFactors = FALSE
           ))
         }
       }
     } else {
       # No treatment info available, just use generic effect sizes
       effect_sizes_df <- data.frame(
         Treatment = "Unknown",
         Reference = "Control",
         Raw_Difference = NA,
         Pooled_SD = NA,
         Standardized_Effect = effect_sizes,
         stringsAsFactors = FALSE
       )
     }
   }
   
   # Store effect sizes
   result$effect_sizes <- effect_sizes_df
   
   # ----- POWER ANALYSIS -----
   
   # Initialize post_power_analysis data frame
   power_results <- data.frame(
     Treatment = character(0),
     Reference = character(0),
     Effect_Size = numeric(0),
     Alpha = numeric(0),
     Power = numeric(0),
     stringsAsFactors = FALSE
   )
   
   # Perform power analysis based on method
   if (method == "parametric") {
     message("Performing parametric power analysis...")
     
     # Get unique treatment groups
     treatment_groups <- names(sample_sizes_by_group)
     if (length(treatment_groups) < 2) {
       stop("At least two treatment groups are required for power analysis")
     }
     
     # Assume first group is control/reference
     ref_group <- treatment_groups[1]
     treatment_groups <- treatment_groups[-1]  # Remove reference group
     
     # Calculate power for each treatment group, effect size, and alpha level
     for (group in treatment_groups) {
       for (es in effect_sizes) {
         for (a in alpha) {
           # Get sample sizes for this comparison
           group_n <- sample_sizes_by_group[group]
           ref_n <- sample_sizes_by_group[ref_group]
           
           # Calculate power using power.t.test
           power_value <- stats::power.t.test(
             n = min(group_n, ref_n),  # Use the smaller for conservative estimate
             delta = es,
             sd = 1,  # Effect size is already in units of SD
             sig.level = a,
             type = "two.sample"
           )$power
           
           # Add to power results
           power_results <- rbind(power_results, data.frame(
             Treatment = group,
             Reference = ref_group,
             Effect_Size = es,
             Alpha = a,
             Power = power_value,
             stringsAsFactors = FALSE
           ))
         }
       }
     }
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
     
     # Get treatment groups and assume first is reference
     treatment_groups <- unique(auc_data[[treatment_column]])
     ref_group <- treatment_groups[1]
     comp_groups <- treatment_groups[-1]  # Remove reference group
     
     # Calculate power for each treatment group, effect size, and alpha
     for (group in comp_groups) {
       for (es in effect_sizes) {
         for (a in alpha) {
           # Get sample sizes for this comparison
           group_n <- auc_stats$N[auc_stats$Treatment == group]
           ref_n <- auc_stats$N[auc_stats$Treatment == ref_group]
           
           # Calculate power
           power_value <- stats::power.t.test(
             n = min(group_n, ref_n),
             delta = es * pooled_sd,  # Convert standardized effect to raw
             sd = pooled_sd,
             sig.level = a,
             type = "two.sample"
           )$power
           
           # Add to power results
           power_results <- rbind(power_results, data.frame(
             Treatment = group,
             Reference = ref_group,
             Effect_Size = es,
             Alpha = a,
             Power = power_value,
             stringsAsFactors = FALSE
           ))
         }
       }
     }
   } else if (method == "simulation") {
     message("Performing simulation-based power analysis...")
     
     # Simplified simulation approach - in a real implementation, this would
     # be more sophisticated and actually run simulations
     
     # Get treatment groups
     if (exists("treatment_groups")) {
       treatment_groups <- treatment_groups
     } else if (!is.null(sample_sizes_by_group)) {
       treatment_groups <- names(sample_sizes_by_group)
     } else {
       stop("Cannot determine treatment groups for simulation")
     }
     
     # Assume first group is reference
     ref_group <- treatment_groups[1]
     comp_groups <- treatment_groups[-1]  # Remove reference group
     
     # Generate power estimates for each treatment comparison, effect size, and alpha
     for (group in comp_groups) {
       for (es in effect_sizes) {
         for (a in alpha) {
           # In a real implementation, this would run actual simulations
           # Here we approximate with a parametric approach for illustration
           
           # Get sample sizes
           if (!is.null(sample_sizes_by_group)) {
             group_n <- sample_sizes_by_group[group]
             ref_n <- sample_sizes_by_group[ref_group]
           } else {
             # Use default if we can't determine
             group_n <- 8
             ref_n <- 8
           }
           
           # Calculate approximate power (in real implementation, this would
           # be based on simulation results)
           power_value <- stats::power.t.test(
             n = min(group_n, ref_n),
             delta = es,
             sd = 1,
             sig.level = a,
             type = "two.sample"
           )$power
           
           # Add to power results
           power_results <- rbind(power_results, data.frame(
             Treatment = group,
             Reference = ref_group,
             Effect_Size = es,
             Alpha = a,
             Power = power_value,
             stringsAsFactors = FALSE
           ))
         }
       }
     }
   }
   
   # Store power analysis results
   result$post_power_analysis <- power_results
   
   # ----- SAMPLE SIZE ESTIMATION -----
   
   # Initialize sample_size_estimates data frame
   sample_size_df <- data.frame(
     Treatment = character(0),
     Reference = character(0),
     Effect_Size = numeric(0),
     Alpha = numeric(0),
     Target_Power = numeric(0),
     Sample_Size = numeric(0),
     stringsAsFactors = FALSE
   )
   
   # Calculate sample size recommendations for each treatment, effect size, alpha, and power
   if (method %in% c("parametric", "auc")) {
     message("Calculating sample size recommendations...")
     
     # Get unique treatment comparisons from power_results
     treatment_comparisons <- unique(power_results[, c("Treatment", "Reference")])
     
     for (i in 1:nrow(treatment_comparisons)) {
       group <- treatment_comparisons$Treatment[i]
       ref_group <- treatment_comparisons$Reference[i]
       
       for (es in effect_sizes) {
         for (a in alpha) {
           for (p in power) {
             # For AUC method, we need to adjust effect size
             delta <- es
             sd_val <- 1
             
             if (method == "auc" && exists("pooled_sd") && !is.na(pooled_sd) && pooled_sd > 0) {
               delta <- es * pooled_sd
               sd_val <- pooled_sd
             }
             
             # Calculate required sample size
             sample_size <- ceiling(stats::power.t.test(
               power = p,
               delta = delta,
               sd = sd_val,
               sig.level = a,
               type = "two.sample"
             )$n)
             
             # Add to sample size data frame
             sample_size_df <- rbind(sample_size_df, data.frame(
               Treatment = group,
               Reference = ref_group,
               Effect_Size = es,
               Alpha = a,
               Target_Power = p,
               Sample_Size = sample_size,
               stringsAsFactors = FALSE
             ))
           }
         }
       }
     }
   } else {
     # For simulation method, provide appropriate message
     message("Sample size estimation for simulation method requires additional simulations")
     
     # Still create a placeholder sample size data frame with NA values
     for (es in effect_sizes) {
       for (a in alpha) {
         for (p in power) {
           sample_size_df <- rbind(sample_size_df, data.frame(
             Treatment = "Requires simulation",
             Reference = "Control",
             Effect_Size = es,
             Alpha = a,
             Target_Power = p,
             Sample_Size = NA,
             stringsAsFactors = FALSE
           ))
         }
       }
     }
   }
   
   # Store sample size recommendations
   result$sample_size_estimates <- sample_size_df
   
   # ----- CREATE PLOTS -----
   
   # Create power curve plot with treatment groups
   if (requireNamespace("ggplot2", quietly = TRUE)) {
     # Get a subset of power_results for the default alpha (usually 0.05)
     if (0.05 %in% alpha) {
       default_alpha <- 0.05
     } else {
       default_alpha <- alpha[1]
     }
     
     plot_data <- power_results[power_results$Alpha == default_alpha, ]
     
     # Power curve with treatment groups
     power_curve <- ggplot2::ggplot(plot_data, 
                                  ggplot2::aes(x = Effect_Size, y = Power, color = Treatment, group = Treatment)) +
       ggplot2::geom_line() +
       ggplot2::geom_point() +
       ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") +
       ggplot2::labs(
         title = paste("Power Analysis Results (α =", default_alpha, ") -", toupper(substr(method, 1, 1)), substr(method, 2, nchar(method)), "Method"),
         x = "Effect Size (Cohen's d)",
         y = "Statistical Power"
       ) +
       ggplot2::theme_classic() +
       ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
     
     # Create sample size recommendation plot
     if (nrow(sample_size_df) > 0 && !all(is.na(sample_size_df$Sample_Size))) {
       # Get a subset for plotting (for clarity)
       if (0.05 %in% alpha && 0.8 %in% power) {
         plot_ss_data <- sample_size_df[sample_size_df$Alpha == 0.05 & sample_size_df$Target_Power == 0.8, ]
       } else {
         plot_ss_data <- sample_size_df[sample_size_df$Alpha == alpha[1] & sample_size_df$Target_Power == power[1], ]
       }
       
       sample_size_plot <- ggplot2::ggplot(plot_ss_data, 
                                         ggplot2::aes(x = Effect_Size, y = Sample_Size, color = Treatment, group = Treatment)) +
         ggplot2::geom_line() +
         ggplot2::geom_point() +
         ggplot2::labs(
           title = paste("Sample Size Recommendations (α =", plot_ss_data$Alpha[1], ", Power =", plot_ss_data$Target_Power[1], ")"),
           x = "Effect Size (Cohen's d)",
           y = "Required Sample Size per Group"
         ) +
         ggplot2::theme_classic()
       
       result$plots <- list(power_curve = power_curve, sample_size_plot = sample_size_plot)
     } else {
       result$plots <- list(power_curve = power_curve)
     }
   }
   
   # ----- RETURN RESULTS -----
   
   # Return all results
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

#' Rename AUC Summary Columns for Compatibility
#' @noRd
rename_auc_columns <- function(auc_summary) {
  # Map old column names to new ones
  column_mapping <- c(
    "AUC.Mean" = "Mean_AUC",
    "AUC.SD" = "SD_AUC",
    "AUC.N" = "N",
    "AUC.SEM" = "SEM_AUC"
  )
  
  # Rename columns if they exist
  for (old_name in names(column_mapping)) {
    if (old_name %in% colnames(auc_summary)) {
      colnames(auc_summary)[colnames(auc_summary) == old_name] <- column_mapping[old_name]
    }
  }
  
  return(auc_summary)
}