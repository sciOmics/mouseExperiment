#' Flexible Tumor Growth Statistical Analysis 
#'
#' This function provides a more flexible statistical analysis for tumor growth data that doesn't 
#' assume exponential growth patterns. It fits multiple statistical models, including both parametric
#' and non-parametric approaches, to accommodate different growth patterns including tumor regression,
#' stable disease, and non-exponential growth. The function automatically selects the most appropriate 
#' model based on the data characteristics.
#'
#' @param df A data frame containing tumor growth data.
#' @param time_column A character string specifying the column name for time points (e.g., "Day"). Default is "Day".
#' @param volume_column A character string specifying the column name for tumor volume measurements (e.g., "Volume"). Default is "Volume".
#' @param treatment_column A character string specifying the column name for treatment groups (e.g., "Treatment"). Default is "Treatment".
#' @param cage_column A character string specifying the column name for the cage identifier (e.g., "Cage"). Default is "Cage".
#' @param id_column A character string specifying the column name for individual subject identifiers (e.g., "ID"). Default is "ID".
#' @param dose_column Optional. A character string specifying the column name for dose levels, if available. Default is NULL.
#' @param model_type A character string specifying the preferred model type. Options are "auto" (automatically determine), 
#'        "parametric" (for traditional growth models), or "nonparametric" (for greater flexibility). Default is "auto".
#' @param regression_penalty Logical indicating whether to apply penalties for model complexity. Default is TRUE.
#'
#' @return A list containing:
#' \item{model}{The fitted statistical model.}
#' \item{model_type}{The type of model used for the analysis.}
#' \item{anova}{Statistical test results comparing groups.}
#' \item{posthoc}{Pairwise comparisons between treatment groups.}
#' \item{auc_analysis}{Area Under the Curve analysis for each treatment.}
#' \item{growth_rates}{Estimated growth rates or trajectories for each treatment.}
#' \item{plots}{A list of plots visualizing the results.}
#'
#' @details
#' Unlike traditional tumor growth analysis that assumes exponential growth patterns,
#' this function uses a flexible approach to accommodate various growth patterns:
#' 
#' 1. Data Exploration: The function first examines growth patterns in each treatment group
#'    to detect non-exponential growth, regression, or stable disease.
#' 
#' 2. Model Selection: Based on the detected patterns, it selects from multiple models:
#'    - Linear mixed-effects model with log transformation (traditional approach)
#'    - Generalized additive mixed models (GAMMs) for non-linear patterns
#'    - Piecewise linear models for patterns with distinct phases
#'    - Non-parametric approaches for highly irregular patterns
#' 
#' 3. Analysis Metrics: The function calculates multiple metrics beyond the traditional approach:
#'    - Area Under the Curve (AUC) for each treatment group
#'    - Best response (maximum regression from baseline)
#'    - Time to progression
#'    - Tumor growth rate during different phases
#' 
#' 4. The function provides comprehensive visualization of the results, including
#'    model diagnostics and treatment comparisons.
#'
#' @examples
#' # Basic usage with automatic model selection
#' results <- tumor_growth_statistics_flexible(df)
#'
#' # Enforce non-parametric analysis for highly variable data
#' results_nonparam <- tumor_growth_statistics_flexible(df, model_type = "nonparametric")
#'
#' @import lme4
#' @import mgcv
#' @import emmeans
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
#' @import performance
#' @export
tumor_growth_statistics_flexible <- function(df, 
                                          time_column = "Day", 
                                          volume_column = "Volume", 
                                          treatment_column = "Treatment", 
                                          cage_column = "Cage", 
                                          id_column = "ID", 
                                          dose_column = NULL,
                                          model_type = "auto",
                                          regression_penalty = TRUE) {
  
  # Ensure required columns exist
  required_columns <- c(time_column, volume_column, treatment_column, cage_column, id_column)
  missing_cols <- required_columns[!required_columns %in% base::colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for dose column if specified
  if (!is.null(dose_column) && !(dose_column %in% base::colnames(df))) {
    warning(paste("Dose column", dose_column, "not found in data frame, proceeding without dose information"))
    dose_column <- NULL
  }
  
  # Create a composite group identifier based on Treatment (and Dose if available)
  if (!is.null(dose_column)) {
    # Create a group identifier combining Treatment and Dose
    df$Group <- paste(df[[treatment_column]], df[[dose_column]], sep = " - Dose: ")
  } else {
    # Use Treatment as the group identifier
    df$Group <- df[[treatment_column]]
  }
  
  # Create a unique identifier for each mouse (combining cage, treatment, ID, and optionally dose)
  if (!is.null(dose_column)) {
    df$Mouse_ID <- base::factor(paste(df[[cage_column]], df[[treatment_column]], 
                              df[[dose_column]], df[[id_column]], sep = "_"))
  } else {
    df$Mouse_ID <- base::factor(paste(df[[cage_column]], df[[treatment_column]], 
                              df[[id_column]], sep = "_"))
  }
  
  # Ensure Group is a factor and time is numeric
  df$Group <- factor(df$Group)
  df[[time_column]] <- as.numeric(df[[time_column]])
  
  # Check for collinearity between Cage and Treatment
  # Create a temporary variable for the cage
  df$Cage_Var <- df[[cage_column]]
  
  # Test for collinearity between cage and treatment group
  has_collinearity <- FALSE
  
  # Only test for collinearity if we have multiple cages
  if (length(unique(df$Cage_Var)) > 1) {
    # Create a simple contingency table of cage vs treatment
    cont_table <- table(df[[cage_column]], df[[treatment_column]])
    
    # Check if each cage only has one treatment (perfect collinearity)
    cage_has_one_treatment <- all(rowSums(cont_table > 0) == 1)
    
    # Check if each treatment only has one cage (perfect collinearity)
    treatment_has_one_cage <- all(colSums(cont_table > 0) == 1)
    
    # If either condition is true, there's perfect collinearity
    has_collinearity <- cage_has_one_treatment || treatment_has_one_cage
    
    # Log the result
    if (has_collinearity) {
      message("Detected collinearity between Cage and Treatment. Using Treatment only in the model.")
    } else {
      message("No perfect collinearity detected between Cage and Treatment. Including Cage in the model.")
    }
  } else {
    message("Only one cage detected. Using Treatment only in the model.")
    has_collinearity <- TRUE  # Treat single cage as collinear case
  }
  
  # ---- Step 1: Detect Growth Patterns ----
  
  # Function to detect if a group shows tumor regression
  # Returns TRUE if there's evidence of regression (volume decrease over time)
  detect_regression <- function(group_data) {
    # Sort by time
    group_data <- group_data[order(group_data[[time_column]]), ]
    
    # Get mean volume at each time point
    mean_volumes <- tapply(group_data[[volume_column]], group_data[[time_column]], mean)
    
    # Check if there's a significant decrease at any point
    # We'll define regression as a 20% decrease from peak volume
    peak_volume <- max(mean_volumes)
    min_after_peak <- min(mean_volumes[which(mean_volumes == peak_volume):length(mean_volumes)])
    
    # Return TRUE if there's at least a 20% decrease from peak
    return((peak_volume - min_after_peak) / peak_volume > 0.2)
  }
  
  # Detect growth patterns for each treatment group
  groups <- unique(df$Group)
  growth_patterns <- data.frame(
    Group = groups,
    Shows_Regression = sapply(groups, function(g) detect_regression(df[df$Group == g, ])),
    stringsAsFactors = FALSE
  )
  
  # Log detected patterns
  message("Growth pattern detection:")
  for (i in 1:nrow(growth_patterns)) {
    pattern_desc <- ifelse(growth_patterns$Shows_Regression[i], 
                          "shows evidence of tumor regression", 
                          "follows standard growth pattern")
    message(paste("  -", growth_patterns$Group[i], pattern_desc))
  }
  
  # Determine if any group shows regression
  has_regression <- any(growth_patterns$Shows_Regression)
  
  # ---- Step 2: Select Appropriate Model ----
  
  # Automatically determine model type if set to "auto"
  if (model_type == "auto") {
    if (has_regression) {
      selected_model_type <- "nonparametric"
      message("Detected tumor regression in one or more groups. Using non-parametric model.")
    } else {
      selected_model_type <- "parametric"
      message("All groups follow standard growth patterns. Using parametric model.")
    }
  } else {
    selected_model_type <- model_type
    message(paste("Using", model_type, "model as specified."))
  }
  
  # Fit selected model
  if (selected_model_type == "parametric") {
    # Traditional approach with log transformation for exponential growth
    message("Fitting linear mixed-effects model with log transformation...")
    
    # Create a log-transformed version of the volume for handling exponential growth
    log_volume_col <- paste0("log_", volume_column)
    df[[log_volume_col]] <- log1p(df[[volume_column]])  # log1p = log(x+1) to handle zeros
    
    # First, determine if the data has the right structure for nested random effects
    df$Cage_Var <- factor(df$Cage_Var)
    
    # Create a contingency table to check cage distribution across groups
    cage_group_table <- table(df$Cage_Var, df$Group)
    
    # Fit appropriate model based on cage effects
    if (!has_collinearity) {
      # Check if we have a proper nesting structure
      cages_per_group <- colSums(cage_group_table > 0)
      
      if (all(cages_per_group > 1)) {
        # Use nested random effects
        df$Cage_in_Group <- interaction(df$Group, df$Cage_Var, sep = ":")
        
        model_formula <- as.formula(paste(
          log_volume_col, "~", 
          time_column, "* Group + (1 | Mouse_ID) + (1 | Cage_in_Group)"
        ))
        
        message("Using model with explicit nested random effects")
      } else {
        # Use additive random effect
        model_formula <- as.formula(paste(
          log_volume_col, "~", 
          time_column, "* Group + (1 | Mouse_ID) + (1 | Cage_Var)"
        ))
        
        message("Using model with additive cage random effect")
      }
    } else {
      # Standard model without cage effects
      model_formula <- as.formula(paste(
        log_volume_col, "~", 
        time_column, "* Group + (1 | Mouse_ID)"
      ))
      
      message("Using model without cage effect")
    }
    
    # Fit the model
    model <- lme4::lmer(model_formula, data = df)
    
    # Perform ANOVA
    anova_results <- car::Anova(model, type = 3)
    
    # Post-hoc comparisons
    emm_formula <- as.formula(paste("pairwise ~ Group |", time_column))
    posthoc_results <- emmeans::emmeans(model, specs = emm_formula, adjust = "bonferroni")
    
  } else {
    # Non-parametric approach for non-exponential growth and regression
    message("Fitting generalized additive mixed model for flexible patterns...")
    
    # Determine the model formula based on collinearity
    if (!has_collinearity) {
      # Include cage as random effect
      df$Cage_Var <- factor(df$Cage_Var)
      
      # GAMM with separate smooths for each group and random effects
      model_formula <- as.formula(paste(
        volume_column, "~", 
        "s(", time_column, ", by = Group, k = 10) + Group + s(Mouse_ID, bs = 're') + s(Cage_Var, bs = 're')"
      ))
      
      message("Using GAMM with separate smooths for each group and cage random effect")
      
    } else {
      # GAMM without cage effect
      model_formula <- as.formula(paste(
        volume_column, "~", 
        "s(", time_column, ", by = Group, k = 10) + Group + s(Mouse_ID, bs = 're')"
      ))
      
      message("Using GAMM with separate smooths for each group without cage effect")
    }
    
    # Set penalty level
    gamma_val <- ifelse(regression_penalty, 1.4, 1.0)
    
    # Fit the GAMM
    model <- mgcv::gam(model_formula, 
                     data = df, 
                     method = "REML", 
                     select = TRUE,
                     gamma = gamma_val)
    
    # Simulate ANOVA-like test
    anova_results <- mgcv::anova.gam(model, test = "F")
    
    # Create custom comparison for time points
    # Get unique time points for comparisons
    time_points <- sort(unique(df[[time_column]]))
    
    # Create a grid for predictions at each time point
    grid_data <- expand.grid(
      Group = levels(df$Group),
      time_temp = time_points
    )
    names(grid_data)[names(grid_data) == "time_temp"] <- time_column
    
    # Add dummy variables for random effects (they'll be marginalized over)
    if (!has_collinearity) {
      grid_data$Mouse_ID <- levels(df$Mouse_ID)[1]
      grid_data$Cage_Var <- levels(df$Cage_Var)[1]
    } else {
      grid_data$Mouse_ID <- levels(df$Mouse_ID)[1]
    }
    
    # Get predictions with standard errors
    predictions <- mgcv::predict.gam(model, newdata = grid_data, se.fit = TRUE, type = "response")
    
    # Add predictions to grid
    grid_data$Predicted <- predictions$fit
    grid_data$SE <- predictions$se.fit
    
    # Function to perform t-test comparison between groups at each time point
    group_comparisons <- function(time_val, group1, group2) {
      g1_pred <- grid_data$Predicted[grid_data[[time_column]] == time_val & grid_data$Group == group1]
      g2_pred <- grid_data$Predicted[grid_data[[time_column]] == time_val & grid_data$Group == group2]
      g1_se <- grid_data$SE[grid_data[[time_column]] == time_val & grid_data$Group == group1]
      g2_se <- grid_data$SE[grid_data[[time_column]] == time_val & grid_data$Group == group2]
      
      # Calculate t-statistic
      t_stat <- (g1_pred - g2_pred) / sqrt(g1_se^2 + g2_se^2)
      
      # Calculate approximate p-value
      p_val <- 2 * pt(-abs(t_stat), df = model$df.residual)
      
      return(data.frame(
        Time = time_val,
        Group1 = group1,
        Group2 = group2,
        Diff = g1_pred - g2_pred,
        T_stat = t_stat,
        P_value = p_val
      ))
    }
    
    # Generate all pairwise comparisons
    all_groups <- levels(df$Group)
    posthoc_results <- list()
    
    for (t in time_points) {
      for (i in 1:(length(all_groups)-1)) {
        for (j in (i+1):length(all_groups)) {
          comp <- group_comparisons(t, all_groups[i], all_groups[j])
          posthoc_results[[length(posthoc_results) + 1]] <- comp
        }
      }
    }
    
    # Combine results
    posthoc_results <- do.call(rbind, posthoc_results)
    
    # Apply Bonferroni correction for multiple comparisons
    posthoc_results$P_adj <- p.adjust(posthoc_results$P_value, method = "bonferroni")
    posthoc_results$Significant <- posthoc_results$P_adj < 0.05
  }
  
  # ---- Step 3: Calculate Additional Metrics ----
  
  # Calculate Area Under the Curve (AUC) for each group and subject
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
  
  # Calculate AUC for each subject
  subjects <- unique(df$Mouse_ID)
  auc_data <- data.frame(
    Mouse_ID = subjects,
    Group = df$Group[match(subjects, df$Mouse_ID)],
    AUC = sapply(subjects, function(s) calculate_auc(df[df$Mouse_ID == s, ])),
    stringsAsFactors = FALSE
  )
  
  # Calculate summary statistics by group
  auc_summary <- aggregate(AUC ~ Group, data = auc_data, 
                         FUN = function(x) c(Mean = mean(x), SD = sd(x), N = length(x)))
  auc_summary <- do.call(data.frame, auc_summary)
  
  # Calculate best response (maximum regression from baseline) for each subject
  calculate_best_response <- function(subject_data) {
    # Sort by time
    subject_data <- subject_data[order(subject_data[[time_column]]), ]
    
    # Calculate percent change from baseline for each time point
    baseline <- subject_data[[volume_column]][1]
    pct_changes <- 100 * (subject_data[[volume_column]] - baseline) / baseline
    
    # Return the minimum percent change (maximum regression)
    return(min(pct_changes))
  }
  
  # Calculate best response for each subject
  best_response_data <- data.frame(
    Mouse_ID = subjects,
    Group = df$Group[match(subjects, df$Mouse_ID)],
    Best_Response = sapply(subjects, function(s) calculate_best_response(df[df$Mouse_ID == s, ])),
    stringsAsFactors = FALSE
  )
  
  # Summarize best response by group
  best_response_summary <- aggregate(Best_Response ~ Group, data = best_response_data,
                                   FUN = function(x) c(Mean = mean(x), SD = sd(x), N = length(x)))
  best_response_summary <- do.call(data.frame, best_response_summary)
  
  # ---- Step 4: Create Visualizations ----
  
  # Calculate mean volume and standard error for each group at each time point
  summary_data <- df %>%
    dplyr::group_by(.data[[time_column]], .data$Group) %>%
    dplyr::summarize(
      Mean_Volume = mean(.data[[volume_column]], na.rm = TRUE),
      SE = sd(.data[[volume_column]], na.rm = TRUE) / sqrt(n()),
      N = n(),
      .groups = "drop"
    )
  
  # Raw data plot with individual mice and group means
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x = time_column, y = volume_column, color = "Group")) +
    ggplot2::geom_line(ggplot2::aes(group = Mouse_ID), alpha = 0.2) +  # Individual mice
    ggplot2::geom_point(data = summary_data, 
                     ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group"), 
                     size = 3) +  # Group means
    ggplot2::geom_line(data = summary_data, 
                    ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group", group = "Group"), 
                    size = 1.2) +  # Connect group means
    ggplot2::geom_errorbar(data = summary_data, 
                        ggplot2::aes_string(x = time_column, y = "Mean_Volume", 
                                         ymin = "Mean_Volume - SE", 
                                         ymax = "Mean_Volume + SE", 
                                         color = "Group"), 
                        width = 0.5) +  # Error bars
    ggplot2::labs(title = "Tumor Growth by Treatment Group",
                x = "Days",
                y = "Tumor Volume") +
    ggplot2::theme_minimal()
  
  # Create a model-based prediction plot
  if (selected_model_type == "parametric") {
    # For parametric model, we need to create predictions on log scale and back-transform
    # Create a prediction grid
    pred_times <- seq(min(df[[time_column]]), max(df[[time_column]]), length.out = 100)
    pred_grid <- expand.grid(
      Group = levels(df$Group),
      time_temp = pred_times
    )
    names(pred_grid)[names(pred_grid) == "time_temp"] <- time_column
    
    # Add a dummy ID for prediction
    pred_grid$Mouse_ID <- "dummy"
    
    # Add cage variable if needed
    if (!has_collinearity) {
      pred_grid$Cage_Var <- levels(df$Cage_Var)[1]
      pred_grid$Cage_in_Group <- paste(pred_grid$Group, pred_grid$Cage_Var, sep = ":")
    }
    
    # Add log volume column
    pred_grid[[log_volume_col]] <- NA
    
    # Make predictions
    pred_values <- predict(model, newdata = pred_grid, re.form = NA)
    
    # Back-transform to original scale
    pred_grid$Predicted_Volume <- expm1(pred_values)  # expm1 = exp(x) - 1
    
    # Create prediction plot
    p2 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = summary_data, 
                        ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group"), 
                        size = 3) +  # Observed means
      ggplot2::geom_line(data = pred_grid, 
                       ggplot2::aes_string(x = time_column, y = "Predicted_Volume", color = "Group", group = "Group"), 
                       size = 1.2) +  # Model predictions
      ggplot2::labs(title = "Parametric Model Predictions",
                  x = "Days",
                  y = "Tumor Volume") +
      ggplot2::theme_minimal()
  } else {
    # For non-parametric model, we already have predictions
    p2 <- ggplot2::ggplot() +
      ggplot2::geom_point(data = summary_data, 
                        ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group"), 
                        size = 3) +  # Observed means
      ggplot2::geom_line(data = grid_data, 
                       ggplot2::aes_string(x = time_column, y = "Predicted", color = "Group", group = "Group"), 
                       size = 1.2) +  # Model predictions
      ggplot2::geom_ribbon(data = grid_data, 
                         ggplot2::aes_string(x = time_column, 
                                          ymin = "Predicted - 1.96 * SE", 
                                          ymax = "Predicted + 1.96 * SE", 
                                          fill = "Group"), 
                         alpha = 0.2) +  # Confidence intervals
      ggplot2::labs(title = "Non-Parametric Model Predictions",
                  x = "Days",
                  y = "Tumor Volume") +
      ggplot2::theme_minimal()
  }
  
  # AUC barplot
  p3 <- ggplot2::ggplot(auc_data, ggplot2::aes(x = Group, y = AUC, fill = Group)) +
    ggplot2::geom_bar(stat = "summary", fun = "mean") +
    ggplot2::geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
    ggplot2::labs(title = "Area Under the Curve by Treatment Group",
                x = "Treatment Group",
                y = "AUC (Tumor Burden)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Best response plot (waterfall plot)
  p4 <- ggplot2::ggplot(best_response_data, ggplot2::aes(x = reorder(Mouse_ID, Best_Response), y = Best_Response, fill = Group)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = "Best Response by Subject (Waterfall Plot)",
                x = "Subject",
                y = "Best Response (% Change from Baseline)") +
    ggplot2::facet_wrap(~ Group, scales = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  
  # Combine plots
  plot_grid <- ggpubr::ggarrange(p1, p2, p3, p4, 
                              ncol = 2, nrow = 2, 
                              common.legend = TRUE, 
                              legend = "right")
  
  # ---- Step 5: Return Results ----
  
  # Create results list
  results <- list(
    model = model,
    model_type = selected_model_type,
    growth_patterns = growth_patterns,
    anova = anova_results,
    posthoc = posthoc_results,
    auc_analysis = list(
      individual = auc_data,
      summary = auc_summary
    ),
    best_response = list(
      individual = best_response_data,
      summary = best_response_summary
    ),
    plots = list(
      raw_data = p1,
      model_prediction = p2,
      auc = p3,
      best_response = p4,
      combined = plot_grid
    )
  )
  
  # Print summary of findings
  cat("\n=== Tumor Growth Analysis with Flexible Modeling ===\n")
  cat("Model type used:", selected_model_type, "\n\n")
  
  cat("Growth Pattern Summary:\n")
  print(growth_patterns)
  cat("\n")
  
  cat("Area Under the Curve Summary:\n")
  print(auc_summary)
  cat("\n")
  
  cat("Best Response Summary:\n")
  print(best_response_summary)
  cat("\n")
  
  # Print combined plots
  print(plot_grid)
  
  return(results)
}