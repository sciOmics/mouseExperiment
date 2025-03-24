#' Post-hoc Power Analysis for Tumor Growth Experiments
#'
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
#'    the tumor growth curves between groups using t-tests or ANOVA, which is appropriate when
#'    the AUC method was used for the primary analysis.
#'
#' The function returns both numerical estimates of power and visualization plots to help
#' interpret the results.
#'
#' @note
#' Post-hoc power analysis should be interpreted with caution, as it is conditional on the
#' observed data and can sometimes lead to misleading conclusions if the effect size is
#' estimated from the same data used for the original hypothesis test.
#'
#' @examples
#' # Load example data
#' data(synthetic_data)
#' 
#' # Add volume calculation
#' df <- calculate_volume(synthetic_data)
#' 
#' # Run post-hoc power analysis on raw data
#' power_results <- post_power_analysis(df)
#' 
#' # View power estimates
#' power_results$power_analysis
#' 
#' # Visualize power curves
#' power_results$plots$power_curve
#' 
#' # Alternatively, run analysis on a fitted model
#' model_results <- tumor_growth_statistics(df)
#' power_from_model <- post_power_analysis(model_results)
#'
#' @seealso \code{\link{tumor_growth_statistics}} for analyzing tumor growth data,
#'          \code{\link[pwr]{pwr.t.test}} for basic power calculations
#'
#' @import ggplot2
#' @importFrom stats var sd power.t.test power.anova.test simulate model.matrix vcov
#' @importFrom lme4 VarCorr fixef getME
#' @importFrom dplyr group_by summarize n
#' @importFrom MASS mvrnorm
#' @export
post_power_analysis <- function(data, 
                               alpha = 0.05,
                               effect_sizes = NULL,
                               method = "parametric",
                               n_simulations = 1000,
                               time_column = "Day",
                               volume_column = "Volume",
                               treatment_column = "Treatment",
                               id_column = "ID") {
  
  # Check if input is a model object from tumor_growth_statistics
  is_model_object <- FALSE
  if(is.list(data) && any(names(data) == "model")) {
    is_model_object <- TRUE
    model <- data$model
    method_from_results <- data$method
    
    # If method wasn't specified, use the one from the results
    if(method == "parametric" && !is.null(method_from_results)) {
      if(method_from_results == "lme") {
        method <- "parametric"
      } else if(method_from_results == "gamm") {
        method <- "simulation"
      } else if(method_from_results == "auc") {
        method <- "auc"
      }
    }
  }
  
  # Validate method
  valid_methods <- c("parametric", "simulation", "auc")
  if(!method %in% valid_methods) {
    stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
  }
  
  # If data is raw data (not a model object), do some preprocessing
  if(!is_model_object) {
    # Check that necessary columns exist
    required_columns <- c(time_column, volume_column, treatment_column, id_column)
    missing_cols <- required_columns[!required_columns %in% colnames(data)]
    
    if(length(missing_cols) > 0) {
      stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
    }
    
    # Create a Group variable if needed
    if(!"Group" %in% colnames(data)) {
      data$Group <- data[[treatment_column]]
    }
  } else {
    # Extract data from model object
    if(method_from_results == "lme" || method_from_results == "gamm") {
      data <- model@frame
    } else if(method_from_results == "auc") {
      data <- data$auc_analysis$individual
    }
  }
  
  # Calculate sample sizes per group
  if("Group" %in% colnames(data)) {
    sample_sizes <- as.data.frame(table(data$Group))
    names(sample_sizes) <- c("Group", "N")
  } else if(treatment_column %in% colnames(data)) {
    sample_sizes <- as.data.frame(table(data[[treatment_column]]))
    names(sample_sizes) <- c("Group", "N")
  } else {
    # If we can't determine groups, just count unique IDs
    n_subjects <- length(unique(data[[id_column]]))
    sample_sizes <- data.frame(Group = "All", N = n_subjects)
  }
  
  # Calculate observed effect sizes if not provided
  if(is.null(effect_sizes)) {
    if(is_model_object) {
      observed_effects <- calculate_observed_effects(data, model, method_from_results)
    } else {
      observed_effects <- calculate_observed_effects_from_raw(data, time_column, volume_column, treatment_column, id_column)
    }
    
    # Generate a sequence of effect sizes around the observed effect
    if(length(observed_effects) > 0) {
      main_effect <- observed_effects$standardized_effect[1]  # Use the first (usually main treatment effect)
      effect_sizes <- seq(main_effect * 0.25, main_effect * 2, length.out = 10)
    } else {
      # If we couldn't estimate effects, use reasonable defaults
      effect_sizes <- seq(0.2, 1.2, by = 0.2)  # Cohen's d values from small to large
    }
  }
  
  # Initialize results list
  results <- list(
    power_analysis = NULL,
    observed_effects = NULL,
    sample_sizes = sample_sizes,
    plots = list()
  )
  
  # Store observed effects
  if(exists("observed_effects")) {
    results$observed_effects <- observed_effects
  }
  
  # Perform power analysis based on chosen method
  if(method == "parametric") {
    power_results <- parametric_power_analysis(data, effect_sizes, alpha, is_model_object, model)
  } else if(method == "simulation") {
    power_results <- simulation_power_analysis(data, effect_sizes, alpha, n_simulations, is_model_object, model)
  } else if(method == "auc") {
    power_results <- auc_power_analysis(data, effect_sizes, alpha, is_model_object, model)
  }
  
  # Store power analysis results
  results$power_analysis <- power_results
  
  # Create visualizations
  results$plots <- create_power_plots(power_results, observed_effects, sample_sizes)
  
  # Print summary
  cat("\n=== Post-hoc Power Analysis Results ===\n")
  cat("Method:", method, "\n")
  cat("Sample sizes per group:\n")
  print(sample_sizes)
  cat("\nPower Analysis Results:\n")
  print(power_results)
  
  # Return results
  return(results)
}

# Helper function to calculate observed effect sizes from model object
calculate_observed_effects <- function(data, model, method) {
  if(method == "lme") {
    # For linear mixed model
    fixed_effects <- lme4::fixef(model)
    
    # Get the variance-covariance matrix
    vcov_matrix <- vcov(model)
    
    # Extract treatment effects (those with "Group" in the name)
    group_effects <- fixed_effects[grep("Group", names(fixed_effects))]
    
    # Get standard errors for the effects
    group_se <- sqrt(diag(vcov_matrix)[grep("Group", names(fixed_effects))])
    
    # Calculate residual variance
    residual_var <- attr(lme4::VarCorr(model), "sc")^2
    
    # Calculate random effects variance
    random_var <- sum(sapply(lme4::VarCorr(model), function(x) x[1]))
    
    # Calculate total variance
    total_var <- residual_var + random_var
    
    # Calculate standardized effect sizes (Cohen's d)
    standardized_effects <- group_effects / sqrt(total_var)
    
    # Create data frame of effects
    observed_effects <- data.frame(
      effect = names(group_effects),
      estimate = group_effects,
      std_error = group_se,
      standardized_effect = standardized_effects
    )
    
  } else if(method == "gamm") {
    # For generalized additive mixed model
    # This is more complex and requires approximation
    summary_model <- summary(model)
    parametric_effects <- summary_model$p.table
    
    # Extract group effects
    group_effects <- parametric_effects[grep("Group", rownames(parametric_effects)), "Estimate"]
    group_se <- parametric_effects[grep("Group", rownames(parametric_effects)), "Std. Error"]
    
    # Approximate standardized effect sizes
    # For GAMs this is more complex, using t-values as approximate standardized effects
    standardized_effects <- parametric_effects[grep("Group", rownames(parametric_effects)), "t value"]
    
    # Create data frame of effects
    observed_effects <- data.frame(
      effect = rownames(parametric_effects)[grep("Group", rownames(parametric_effects))],
      estimate = group_effects,
      std_error = group_se,
      standardized_effect = standardized_effects / 2  # Dividing by 2 for more conservative estimate
    )
    
  } else if(method == "auc") {
    # For AUC analysis
    auc_data <- data
    
    # Calculate mean and SD for each group
    group_summary <- aggregate(AUC ~ Group, data = auc_data, 
                             FUN = function(x) c(Mean = mean(x), SD = sd(x), N = length(x)))
    group_summary <- do.call(data.frame, group_summary)
    
    # Calculate pooled SD
    pooled_sd <- sqrt(sum((group_summary$AUC.N - 1) * group_summary$AUC.SD^2) / 
                      sum(group_summary$AUC.N - 1))
    
    # Calculate effect sizes between all pairs of groups
    n_groups <- nrow(group_summary)
    effect_pairs <- list()
    
    if(n_groups > 1) {
      for(i in 1:(n_groups-1)) {
        for(j in (i+1):n_groups) {
          group1 <- as.character(group_summary$Group[i])
          group2 <- as.character(group_summary$Group[j])
          
          # Cohen's d
          d <- abs(group_summary$AUC.Mean[i] - group_summary$AUC.Mean[j]) / pooled_sd
          
          effect_pairs[[length(effect_pairs) + 1]] <- data.frame(
            effect = paste(group1, "vs", group2),
            estimate = group_summary$AUC.Mean[i] - group_summary$AUC.Mean[j],
            std_error = pooled_sd * sqrt(1/group_summary$AUC.N[i] + 1/group_summary$AUC.N[j]),
            standardized_effect = d
          )
        }
      }
      observed_effects <- do.call(rbind, effect_pairs)
    } else {
      observed_effects <- data.frame(
        effect = character(),
        estimate = numeric(),
        std_error = numeric(),
        standardized_effect = numeric()
      )
    }
  } else {
    # If method is unknown
    observed_effects <- data.frame(
      effect = character(),
      estimate = numeric(),
      std_error = numeric(),
      standardized_effect = numeric()
    )
  }
  
  return(observed_effects)
}

# Helper function to calculate effect sizes from raw data
calculate_observed_effects_from_raw <- function(data, time_column, volume_column, treatment_column, id_column) {
  # Ensure we have a Group variable
  if(!"Group" %in% colnames(data)) {
    data$Group <- data[[treatment_column]]
  }
  
  # If we have volume data over time, calculate AUC for each subject
  if(all(c(time_column, volume_column, id_column) %in% colnames(data))) {
    # Calculate AUC for each subject
    subjects <- unique(data[[id_column]])
    
    calculate_auc <- function(subject_data) {
      # Sort by time
      subject_data <- subject_data[order(subject_data[[time_column]]), ]
      
      # Calculate AUC using trapezoidal rule
      times <- subject_data[[time_column]]
      volumes <- subject_data[[volume_column]]
      
      auc <- 0
      for (i in 2:length(times)) {
        dt <- times[i] - times[i-1]
        auc <- auc + dt * (volumes[i] + volumes[i-1]) / 2
      }
      
      return(auc)
    }
    
    auc_data <- data.frame(
      ID = subjects,
      Group = data$Group[match(subjects, data[[id_column]])],
      AUC = sapply(subjects, function(s) {
        subject_data <- data[data[[id_column]] == s, ]
        calculate_auc(subject_data)
      }),
      stringsAsFactors = FALSE
    )
    
    # Calculate mean and SD for each group
    group_summary <- aggregate(AUC ~ Group, data = auc_data, 
                             FUN = function(x) c(Mean = mean(x), SD = sd(x), N = length(x)))
    group_summary <- do.call(data.frame, group_summary)
    
    # Calculate pooled SD
    pooled_sd <- sqrt(sum((group_summary$AUC.N - 1) * group_summary$AUC.SD^2) / 
                      sum(group_summary$AUC.N - 1))
    
    # Calculate effect sizes between all pairs of groups
    n_groups <- nrow(group_summary)
    effect_pairs <- list()
    
    if(n_groups > 1) {
      for(i in 1:(n_groups-1)) {
        for(j in (i+1):n_groups) {
          group1 <- as.character(group_summary$Group[i])
          group2 <- as.character(group_summary$Group[j])
          
          # Cohen's d
          d <- abs(group_summary$AUC.Mean[i] - group_summary$AUC.Mean[j]) / pooled_sd
          
          effect_pairs[[length(effect_pairs) + 1]] <- data.frame(
            effect = paste(group1, "vs", group2),
            estimate = group_summary$AUC.Mean[i] - group_summary$AUC.Mean[j],
            std_error = pooled_sd * sqrt(1/group_summary$AUC.N[i] + 1/group_summary$AUC.N[j]),
            standardized_effect = d
          )
        }
      }
      observed_effects <- do.call(rbind, effect_pairs)
    } else {
      observed_effects <- data.frame(
        effect = character(),
        estimate = numeric(),
        std_error = numeric(),
        standardized_effect = numeric()
      )
    }
  } else {
    # If we don't have time and volume data, just compare final volumes between groups
    observed_effects <- data.frame(
      effect = character(),
      estimate = numeric(),
      std_error = numeric(),
      standardized_effect = numeric()
    )
  }
  
  return(observed_effects)
}

# Function to perform parametric power analysis
parametric_power_analysis <- function(data, effect_sizes, alpha, is_model_object, model) {
  if(is_model_object) {
    # Extract needed information from model
    n_per_group <- table(model@frame$Group)
    
    # Get variance components for power calculation
    residual_var <- attr(lme4::VarCorr(model), "sc")^2
    random_var <- sum(sapply(lme4::VarCorr(model), function(x) x[1]))
    total_var <- residual_var + random_var
    
    # Mean group size
    avg_n <- mean(n_per_group)
    n_groups <- length(n_per_group)
  } else {
    # Calculate from raw data
    n_per_group <- table(data$Group)
    avg_n <- mean(n_per_group)
    n_groups <- length(n_per_group)
    
    # Estimate variance components (this is a rough approximation)
    # Calculate total variance
    if("Volume" %in% colnames(data)) {
      total_var <- var(data$Volume, na.rm = TRUE)
    } else {
      # Use a reasonable default if we don't have volume data
      total_var <- 1
    }
  }
  
  # Calculate power for each effect size
  power_results <- data.frame(
    effect_size = effect_sizes,
    power = numeric(length(effect_sizes))
  )
  
  for(i in 1:length(effect_sizes)) {
    if(n_groups == 2) {
      # For 2 groups, use t-test power
      power_obj <- stats::power.t.test(
        n = avg_n,
        delta = effect_sizes[i] * sqrt(total_var),
        sd = sqrt(total_var),
        sig.level = alpha,
        type = "two.sample",
        alternative = "two.sided"
      )
      power_results$power[i] <- power_obj$power
    } else {
      # For >2 groups, use ANOVA power
      power_obj <- stats::power.anova.test(
        groups = n_groups,
        n = avg_n,
        between.var = total_var * effect_sizes[i]^2,
        within.var = total_var,
        sig.level = alpha
      )
      power_results$power[i] <- power_obj$power
    }
  }
  
  return(power_results)
}

# Function to perform simulation-based power analysis
simulation_power_analysis <- function(data, effect_sizes, alpha, n_simulations, is_model_object, model) {
  # This is a simplified simulation approach
  if(is_model_object) {
    # Extract parameters from model
    if(class(model)[1] == "lmerMod") {
      # For linear mixed models
      fixed_effects <- lme4::fixef(model)
      random_effects_var <- lme4::VarCorr(model)
      residual_var <- attr(random_effects_var, "sc")^2
      
      # Design matrix components
      X <- lme4::getME(model, "X")
      Z <- lme4::getME(model, "Z")
      
      # Get groups and time points
      groups <- model@frame$Group
      times <- model@frame[[names(model@frame)[1]]]  # Assuming time is first column
      
    } else if(class(model)[1] == "gam") {
      # For GAM models, we need a different approach
      # This is complex and usually requires specialized implementations
      # Converting to approximation with parametric power analysis
      return(parametric_power_analysis(data, effect_sizes, alpha, is_model_object, model))
    }
  } else {
    # For raw data, build a simple model first
    if(!all(c("Group", "ID", "Volume") %in% colnames(data))) {
      stop("Raw data must contain Group, ID, and Volume columns for simulation power analysis")
    }
    
    # Fit a simple linear model on final volume
    final_data <- aggregate(Volume ~ Group + ID, data = data, FUN = max)
    lm_model <- lm(Volume ~ Group, data = final_data)
    
    # Extract parameters
    fixed_effects <- coef(lm_model)
    residual_var <- summary(lm_model)$sigma^2
    
    # Design matrix
    X <- model.matrix(~ Group, data = final_data)
    
    # Groups
    groups <- final_data$Group
  }
  
  # Initialize power results
  power_results <- data.frame(
    effect_size = effect_sizes,
    power = numeric(length(effect_sizes))
  )
  
  # Run simulations for each effect size
  for(i in 1:length(effect_sizes)) {
    # Scale the treatment effect by the desired effect size
    if(is_model_object) {
      # Find treatment effects (those with "Group" in name)
      treatment_idx <- grep("Group", names(fixed_effects))
      if(length(treatment_idx) > 0) {
        # Scale the effects
        scaled_effects <- fixed_effects
        scaled_effects[treatment_idx] <- fixed_effects[treatment_idx] * effect_sizes[i]
      } else {
        scaled_effects <- fixed_effects
      }
    } else {
      # For simple model
      scaled_effects <- fixed_effects
      if(length(fixed_effects) > 1) {
        scaled_effects[-1] <- fixed_effects[-1] * effect_sizes[i]
      }
    }
    
    # Count significant results
    significant_count <- 0
    
    for(sim in 1:n_simulations) {
      # Generate simulated data
      if(is_model_object && class(model)[1] == "lmerMod") {
        # Simulate random effects
        random_effects <- list()
        for(re in names(random_effects_var)) {
          var_matrix <- as.matrix(random_effects_var[[re]])
          n_levels <- dim(var_matrix)[1]
          random_effects[[re]] <- MASS::mvrnorm(n_levels, rep(0, nrow(var_matrix)), var_matrix)
        }
        
        # Simulate residuals
        epsilon <- rnorm(nrow(X), 0, sqrt(residual_var))
        
        # Combine fixed and random effects
        y_sim <- X %*% scaled_effects + epsilon
        
        # Refit model
        sim_data <- data.frame(
          y = y_sim,
          Group = groups,
          Time = times
        )
        
        sim_model <- lm(y ~ Group, data = sim_data)
        
      } else {
        # For simple model
        epsilon <- rnorm(nrow(X), 0, sqrt(residual_var))
        y_sim <- X %*% scaled_effects + epsilon
        
        sim_data <- data.frame(
          y = y_sim,
          Group = groups
        )
        
        sim_model <- lm(y ~ Group, data = sim_data)
      }
      
      # Test for significance (do any group coefficients have p < alpha?)
      p_values <- summary(sim_model)$coefficients[-1, "Pr(>|t|)"]
      if(any(p_values < alpha)) {
        significant_count <- significant_count + 1
      }
    }
    
    # Calculate power
    power_results$power[i] <- significant_count / n_simulations
  }
  
  return(power_results)
}

# Function to perform AUC-based power analysis
auc_power_analysis <- function(data, effect_sizes, alpha, is_model_object, model) {
  if(is_model_object) {
    # For results from tumor_growth_statistics with AUC method
    auc_data <- data$auc_analysis$individual
    groups <- auc_data$Group
    auc_values <- auc_data$AUC
  } else {
    # Calculate AUC from raw data
    if(!all(c("Group", "ID", "Day", "Volume") %in% colnames(data))) {
      # Try with common alternative names
      possible_time_cols <- c("Day", "Time", "Days", "TimePoint")
      possible_volume_cols <- c("Volume", "TumorVolume", "Tumor_Volume", "Size")
      
      time_col <- NULL
      for(col in possible_time_cols) {
        if(col %in% colnames(data)) {
          time_col <- col
          break
        }
      }
      
      volume_col <- NULL
      for(col in possible_volume_cols) {
        if(col %in% colnames(data)) {
          volume_col <- col
          break
        }
      }
      
      if(is.null(time_col) || is.null(volume_col)) {
        stop("Raw data must contain time and volume columns for AUC power analysis")
      }
    } else {
      time_col <- "Day"
      volume_col <- "Volume"
    }
    
    # Calculate AUC for each subject
    subjects <- unique(data$ID)
    
    calculate_auc <- function(subject_data) {
      # Sort by time
      subject_data <- subject_data[order(subject_data[[time_col]]), ]
      
      # Calculate AUC using trapezoidal rule
      times <- subject_data[[time_col]]
      volumes <- subject_data[[volume_col]]
      
      auc <- 0
      for (i in 2:length(times)) {
        dt <- times[i] - times[i-1]
        auc <- auc + dt * (volumes[i] + volumes[i-1]) / 2
      }
      
      return(auc)
    }
    
    auc_data <- data.frame(
      ID = subjects,
      Group = data$Group[match(subjects, data$ID)],
      AUC = sapply(subjects, function(s) {
        subject_data <- data[data$ID == s, ]
        calculate_auc(subject_data)
      }),
      stringsAsFactors = FALSE
    )
    
    groups <- auc_data$Group
    auc_values <- auc_data$AUC
  }
  
  # Calculate mean and SD for each group
  group_summary <- aggregate(AUC ~ Group, data = auc_data, 
                           FUN = function(x) c(Mean = mean(x), SD = sd(x), N = length(x)))
  group_summary <- do.call(data.frame, group_summary)
  
  # Calculate pooled SD
  pooled_sd <- sqrt(sum((group_summary$AUC.N - 1) * group_summary$AUC.SD^2) / 
                    sum(group_summary$AUC.N - 1))
  
  # Initialize power results
  power_results <- data.frame(
    effect_size = effect_sizes,
    power = numeric(length(effect_sizes))
  )
  
  # Count unique groups
  n_groups <- length(unique(groups))
  
  # Calculate power for each effect size
  for(i in 1:length(effect_sizes)) {
    if(n_groups == 2) {
      # For 2 groups, use t-test power
      power_obj <- stats::power.t.test(
        n = min(group_summary$AUC.N),
        delta = effect_sizes[i] * pooled_sd,
        sd = pooled_sd,
        sig.level = alpha,
        type = "two.sample",
        alternative = "two.sided"
      )
      power_results$power[i] <- power_obj$power
    } else {
      # For >2 groups, use ANOVA power
      power_obj <- stats::power.anova.test(
        groups = n_groups,
        n = min(group_summary$AUC.N),
        between.var = pooled_sd^2 * effect_sizes[i]^2,
        within.var = pooled_sd^2,
        sig.level = alpha
      )
      power_results$power[i] <- power_obj$power
    }
  }
  
  return(power_results)
}

# Function to create visualization plots
create_power_plots <- function(power_results, observed_effects, sample_sizes) {
  # Power curve plot
  power_curve <- ggplot2::ggplot(power_results, ggplot2::aes(x = effect_size, y = power)) +
    ggplot2::geom_line(size = 1.2, color = "blue") +
    ggplot2::geom_point(size = 3, color = "blue") +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    ggplot2::annotate("text", x = min(power_results$effect_size), y = 0.82, 
                    label = "80% Power", hjust = 0, color = "red") +
    ggplot2::labs(
      title = "Power Analysis Curve",
      x = "Standardized Effect Size (Cohen's d)",
      y = "Statistical Power (1 - β)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
  
  # If we have observed effects, add them to the plot
  if(!is.null(observed_effects) && nrow(observed_effects) > 0) {
    # Find where the observed effects are on the power curve
    for(i in 1:nrow(observed_effects)) {
      effect_size <- observed_effects$standardized_effect[i]
      
      # Only add if effect size is within the range of our power curve
      if(effect_size >= min(power_results$effect_size) && 
         effect_size <= max(power_results$effect_size)) {
        
        # Interpolate to find the power at this effect size
        power_at_effect <- approx(power_results$effect_size, 
                                power_results$power, 
                                xout = effect_size)$y
        
        # Add point and text
        power_curve <- power_curve +
          ggplot2::geom_point(data = data.frame(x = effect_size, y = power_at_effect),
                           ggplot2::aes(x = x, y = y), 
                           color = "green4", size = 4, shape = 17) +
          ggplot2::annotate("text", x = effect_size, y = power_at_effect + 0.05,
                          label = paste("Observed:", observed_effects$effect[i]),
                          color = "green4", angle = 0, hjust = 0, size = 3)
      }
    }
  }
  
  # Sample size vs power plot
  # Generate data for a range of sample sizes
  n_range <- seq(5, max(30, max(sample_sizes$N) * 1.5), by = 5)
  
  # Use one representative effect size (median of those tested)
  median_effect <- median(power_results$effect_size)
  
  # Calculate power for different sample sizes
  n_power_data <- data.frame(
    n = n_range,
    power = sapply(n_range, function(n) {
      if(nrow(sample_sizes) == 2) {
        # Two groups - use t-test
        power_obj <- stats::power.t.test(
          n = n,
          delta = median_effect,
          sd = 1,  # Standardized effect uses SD=1
          sig.level = 0.05,
          type = "two.sample",
          alternative = "two.sided"
        )
        return(power_obj$power)
      } else {
        # Multiple groups - use ANOVA
        power_obj <- stats::power.anova.test(
          groups = nrow(sample_sizes),
          n = n,
          between.var = median_effect^2,
          within.var = 1,
          sig.level = 0.05
        )
        return(power_obj$power)
      }
    })
  )
  
  # Create the plot
  n_power_plot <- ggplot2::ggplot(n_power_data, ggplot2::aes(x = n, y = power)) +
    ggplot2::geom_line(size = 1.2, color = "purple") +
    ggplot2::geom_point(size = 3, color = "purple") +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = paste("Power vs. Sample Size (Effect Size =", round(median_effect, 2), ")"),
      x = "Sample Size per Group",
      y = "Statistical Power (1 - β)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
  
  # Add points for current sample sizes
  for(i in 1:nrow(sample_sizes)) {
    current_n <- sample_sizes$N[i]
    # Interpolate to find power at this sample size
    power_at_n <- approx(n_power_data$n, n_power_data$power, xout = current_n)$y
    
    if(!is.na(power_at_n)) {
      n_power_plot <- n_power_plot +
        ggplot2::geom_point(data = data.frame(x = current_n, y = power_at_n),
                         ggplot2::aes(x = x, y = y), 
                         color = "orange", size = 4, shape = 18) +
        ggplot2::annotate("text", x = current_n, y = power_at_n + 0.05,
                        label = paste("Current n =", current_n, 
                                      "\nGroup:", sample_sizes$Group[i]),
                        color = "orange", hjust = 0, size = 3)
    }
  }
  
  # Combine plots
  plots <- list(
    power_curve = power_curve,
    sample_size_curve = n_power_plot,
    combined = ggpubr::ggarrange(power_curve, n_power_plot, 
                              ncol = 2, 
                              labels = c("A", "B"))
  )
  
  return(plots)
}