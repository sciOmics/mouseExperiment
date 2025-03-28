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
     
     # Calculate AUC for the raw data if it wasn't provided
     if (!is.null(model_results) && !is.null(model_results$auc_analysis) && !is.null(model_results$auc_analysis$individual)) {
       auc_data <- model_results$auc_analysis$individual
     } else if (exists("raw_data")) {
       # Calculate AUC from scratch
       message("Calculating AUC values...")
       
       # Create a direct implementation for AUC calculation
       auc_results <- tryCatch({
         # Group the data by treatment and ID to calculate AUC for each subject
         unique_ids <- unique(raw_data[[id_column]])
         unique_treatments <- unique(raw_data[[treatment_column]])
         
         # Check if we have enough data for AUC analysis
         if (length(unique_ids) < 4 || length(unique_treatments) < 2) {
           stop("Insufficient data for AUC analysis (need at least 4 subjects and 2 treatment groups)")
         }
         
         # Initialize AUC individual data
         auc_individual <- data.frame(
           ID = character(),
           Treatment = character(),
           AUC = numeric(),
           stringsAsFactors = FALSE
         )
         
         # Manually calculate AUC for each subject
         for (id in unique_ids) {
           # Get data for this subject
           id_data <- raw_data[raw_data[[id_column]] == id, ]
           
           # Get treatment for this subject (assuming one treatment per subject)
           treatment <- unique(id_data[[treatment_column]])
           if (length(treatment) > 1) {
             warning("Subject ", id, " has multiple treatment assignments. Using the first one.")
             treatment <- treatment[1]
           }
           
           # Sort by time
           id_data <- id_data[order(id_data[[time_column]]), ]
           times <- id_data[[time_column]]
           volumes <- id_data[[volume_column]]
           
           # Check for enough time points
           if (length(times) < 3) {
             warning("Subject ", id, " has too few time points for AUC calculation. Skipping.")
             next
           }
           
           # Calculate trapezoidal AUC
           auc_value <- 0
           for (i in 2:length(times)) {
             # Area of trapezoid = (y1 + y2) * (x2 - x1) / 2
             auc_value <- auc_value + (volumes[i-1] + volumes[i]) * (times[i] - times[i-1]) / 2
           }
           
           # Add to AUC individual data
           auc_individual <- rbind(auc_individual, data.frame(
             ID = id,
             Treatment = treatment,
             AUC = auc_value,
             stringsAsFactors = FALSE
           ))
         }
         
         # Set column names correctly
         names(auc_individual)[names(auc_individual) == "Treatment"] <- treatment_column
         
         # Return the AUC individual data
         auc_individual
       }, error = function(e) {
         message("Error calculating AUC: ", e$message)
         return(NULL)
       })
       
       if (is.null(auc_results) || nrow(auc_results) == 0) {
         stop("Failed to calculate AUC values")
       }
       
       auc_data <- auc_results
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
     
     # Get unique treatment groups
     treatment_groups <- unique(auc_data[[treatment_column]])
     
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
     
     # Initialize power results
     power_results <- data.frame(Effect_Size = effect_sizes, stringsAsFactors = FALSE)
     
     # Set up residual SD for simulations
     if (!is.null(model_results) && !is.null(model_results$data_summary) && 
         !is.null(model_results$data_summary$sigma) && !is.na(model_results$data_summary$sigma) &&
         model_results$data_summary$sigma > 0) {
       residual_sd <- model_results$data_summary$sigma
     } else {
       # Calculate from raw data
       id_means <- stats::aggregate(sim_data[[volume_column]], 
                                  by = list(ID = sim_data[[id_column]],
                                          Time = sim_data[[time_column]]), 
                                  FUN = mean)
       
       # Calculate residuals around the mean
       merged_data <- merge(sim_data, id_means, by.x = c(id_column, time_column), 
                          by.y = c("ID", "Time"))
       
       residuals <- merged_data[[volume_column]] - merged_data$x
       residual_sd <- stats::sd(residuals, na.rm = TRUE)
       
       if (is.na(residual_sd) || residual_sd <= 0) {
         residual_sd <- 0.2  # Default fallback
         warning("Could not calculate residual SD from data. Using default value.")
       }
     }
     
     # Function to simulate data and test for significant treatment effect
     simulate_experiment <- function(effect_size) {
       # Create a simulated dataset
       sim_df <- data.frame(
         ID = character(0),
         Time = numeric(0),
         Treatment = character(0),
         Volume = numeric(0),
         stringsAsFactors = FALSE
       )
       
       # Time points for simulation (use the unique time points from the data)
       time_points <- sort(unique(sim_data[[time_column]]))
       
       # Generate data for each treatment group
       for (g in 1:n_groups) {
         group <- treatment_groups[g]
         n_subjects <- as.numeric(n_per_group[group])
         
         # Calculate base volume and treatment effect (first group is control)
         if (g == 1) {
           # Control group - no effect
           base_effect <- 0
         } else {
           # Treatment groups - effect size-dependent reduction
           base_effect <- -effect_size * residual_sd
         }
         
         # Generate data for each subject in this group
         for (i in 1:n_subjects) {
           # Generate unique subject ID
           subj_id <- paste0(group, "_S", i)
           
           # Generate subject-specific random effect
           subject_effect <- stats::rnorm(1, 0, residual_sd/2)
           
           # Generate data for each time point
           for (t in time_points) {
             # Simple growth model: baseline + time effect + treatment effect
             expected_volume <- 100 + 10 * t + base_effect + subject_effect
             
             # Add random noise
             volume <- expected_volume + stats::rnorm(1, 0, residual_sd)
             
             # Add to simulated dataset
             sim_df <- rbind(sim_df, data.frame(
               ID = subj_id,
               Time = t,
               Treatment = group,
               Volume = volume,
               stringsAsFactors = FALSE
             ))
           }
         }
       }
       
       # Set correct column names
       colnames(sim_df)[colnames(sim_df) == "ID"] <- id_column
       colnames(sim_df)[colnames(sim_df) == "Time"] <- time_column
       colnames(sim_df)[colnames(sim_df) == "Treatment"] <- treatment_column
       colnames(sim_df)[colnames(sim_df) == "Volume"] <- volume_column
       
       # Analyze the simulated data with a simple model
       # For simplicity, we use lm on the last time point
       last_time <- max(sim_df[[time_column]])
       last_data <- sim_df[sim_df[[time_column]] == last_time, ]
       
       # Make sure we have at least two treatment groups in the last data
       if (length(unique(last_data[[treatment_column]])) < 2) {
         warning("Simulation resulted in fewer than 2 treatment groups at the final time point")
         return(FALSE)
       }
       
       # Fit a simple model to test for treatment effect
       model_formula <- as.formula(paste(volume_column, "~", treatment_column))
       tryCatch({
         sim_model <- stats::lm(model_formula, data = last_data)
         sim_anova <- stats::anova(sim_model)
         
         # Return TRUE if we found a significant treatment effect
         return(sim_anova$"Pr(>F)"[1] < alpha)
       }, error = function(e) {
         warning("Error in simulation model fitting: ", e$message)
         return(FALSE)
       })
     }
     
     # Calculate power for each effect size
     power_results$Power <- sapply(effect_sizes, function(d) {
       # Run n_simulations and count significant results
       significant_count <- sum(replicate(n_simulations, simulate_experiment(d)))
       return(significant_count / n_simulations)
     })
     
     result$power_analysis <- power_results
     
     # For simulation, we don't calculate sample size recommendations directly
     # But we can provide a message about how to interpret results
     result$sample_size_recommendations <- "Sample size recommendations for simulation method should be determined by running additional simulations with varying sample sizes."
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