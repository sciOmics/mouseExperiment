#' Analyze Tumor Growth Using Various Statistical Methods
#'
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
#'        This allows for dose-response analysis within treatment groups. Default is NULL.
#' @param method A character string specifying which statistical model to use. Options are:
#'        "lme" (linear mixed-effects with log transformation),
#'        "gamm" (generalized additive mixed model), or
#'        "auc" (area under the curve analysis). Default is "lme".
#' @param regression_penalty Logical indicating whether to apply penalties for model complexity in the GAMM model.
#'        If TRUE (default), applies a gamma value of 1.4 to prevent overfitting. If FALSE, uses the default gamma of 1.0.
#'        Only used when method="gamm".
#'
#' @return A list containing:
#' \item{model}{The fitted statistical model (lmer object for "lme", gam object for "gamm", NULL for "auc" method).}
#' \item{method}{The analysis method used (string indicating which method was used).}
#' \item{anova}{Statistical test results comparing groups (Anova object for "lme", anova.gam results for "gamm", 
#'              aov results for "auc").}
#' \item{posthoc}{Pairwise comparisons between treatment groups (emmeans object for "lme", data frame of comparisons 
#'               for "gamm", pairwise.t.test results for "auc").}
#' \item{auc_analysis}{Area Under the Curve analysis for each treatment. Contains two elements:
#'                   \code{individual} (data frame with AUC values for each subject) and
#'                   \code{summary} (data frame with mean, SD, and N for each treatment group).}
#' \item{plots}{A list of plots visualizing the results, containing the following ggplot objects:
#'             \code{raw_data} (tumor growth curves),
#'             method-specific plots like \code{transformed_data}, \code{model_prediction}, \code{auc_barplot}, etc., and 
#'             \code{combined} (a composite plot with the main visualizations).}
#'
#' @details
#' This function supports different analytical approaches depending on the nature of your tumor growth data:
#' 
#' 1. Linear Mixed-Effects Model ("lme"): 
#'    - Applies a log transformation to tumor volume to accommodate exponential growth
#'    - Uses Mouse_ID as a random effect to account for repeated measures
#'    - Intelligently handles cage effects based on collinearity testing
#'    - Performs Type III ANOVA and Bonferroni-adjusted pairwise comparisons
#'    - Best for data following expected exponential growth patterns
#'
#' 2. Generalized Additive Mixed Model ("gamm"):
#'    - Uses flexible non-linear smoothing terms for each treatment group
#'    - Provides better fit for non-exponential growth patterns, including regression
#'    - Still accounts for random effects of Mouse_ID and Cage
#'    - Calculates pairwise comparisons at each time point
#'    - Best for data showing complex growth patterns or tumor regression
#'
#' 3. Area Under the Curve Analysis ("auc"):
#'    - Calculates the AUC for each individual subject using the trapezoidal method
#'    - Compares AUCs between treatment groups
#'    - Treats the AUC as a summary metric of tumor burden over time
#'    - Performs statistical comparisons using ANOVA and t-tests
#'    - Best for simplifying analysis to a single metric of overall tumor burden
#'
#' The function automatically handles collinearity between cage and treatment variables
#' and provides comprehensive visualization of the results. It creates a unique identifier for each 
#' mouse by combining cage, treatment, and ID information, and uses this to track individual animals 
#' across time points.
#'
#' For cage effects, the function checks if there is perfect collinearity between cage and treatment 
#' (e.g., if each cage only contains animals from one treatment group). If collinearity exists, 
#' only treatment is included in the model. If no collinearity exists, cage effects are included 
#' either as nested random effects (if multiple cages per treatment) or as additive random effects.
#'
#' All methods produce publication-quality plots that visualize both the raw data and model-specific 
#' results. The plots can be accessed from the returned list and further customized as needed.
#'
#' @examples
#' # Load example data
#' data(synthetic_data)
#' 
#' # Using linear mixed-effects model (default)
#' results_lme <- tumor_growth_statistics(synthetic_data)
#' 
#' # View summary of LME results
#' summary(results_lme$model)
#' 
#' # Access the ANOVA results
#' results_lme$anova
#' 
#' # Using generalized additive mixed model for complex growth patterns
#' results_gamm <- tumor_growth_statistics(synthetic_data, method = "gamm")
#' 
#' # Using area under the curve analysis 
#' results_auc <- tumor_growth_statistics(synthetic_data, method = "auc")
#' 
#' # Compare AUC values between treatment groups
#' results_auc$anova
#' results_auc$posthoc
#' 
#' # Custom column names
#' results <- tumor_growth_statistics(df, 
#'                                 time_column = "TimePoint", 
#'                                 volume_column = "TumorSize",
#'                                 treatment_column = "Group",
#'                                 method = "lme")
#'
#' @seealso \code{\link[lme4]{lmer}} for details on linear mixed-effects models, 
#'          \code{\link[mgcv]{gam}} for details on generalized additive models, 
#'          \code{\link{plot_tumor_growth}} for visualizing tumor growth curves
#'
#' @import lme4
#' @import mgcv
#' @import emmeans
#' @import ggplot2
#' @import ggpubr
#' @import dplyr
#' @import performance
#' @importFrom stats pt p.adjust aov pairwise.t.test
#' @importFrom car Anova
#' @export
tumor_growth_statistics <- function(df, 
                                  time_column = "Day", 
                                  volume_column = "Volume", 
                                  treatment_column = "Treatment", 
                                  cage_column = "Cage", 
                                  id_column = "ID", 
                                  dose_column = NULL,
                                  method = "lme",
                                  regression_penalty = TRUE) {
  
  # Ensure required columns exist
  required_columns <- c(time_column, volume_column, treatment_column, cage_column, id_column)
  missing_cols <- required_columns[!required_columns %in% base::colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check method parameter
  valid_methods <- c("lme", "gamm", "auc")
  if (!method %in% valid_methods) {
    stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
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
  
  # Create a contingency table to check cage distribution across groups
  df$Cage_Var <- factor(df$Cage_Var)
  cage_group_table <- table(df$Cage_Var, df$Group)
  print("Cage distribution across treatment groups:")
  print(cage_group_table)
  
  # Initialize results list
  results <- list(
    model = NULL,
    method = method,
    anova = NULL,
    posthoc = NULL,
    auc_analysis = NULL,
    plots = list()
  )
  
  # Calculate AUC for each subject (used by all methods)
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
  
  # Store basic AUC analysis in results
  results$auc_analysis <- list(
    individual = auc_data,
    summary = auc_summary
  )
  
  # Calculate mean volume and standard error for each group at each time point (for plotting)
  summary_data <- df %>%
    dplyr::group_by(.data[[time_column]], .data$Group) %>%
    dplyr::summarize(
      Mean_Volume = mean(.data[[volume_column]], na.rm = TRUE),
      SE = sd(.data[[volume_column]], na.rm = TRUE) / sqrt(n()),
      N = n(),
      .groups = "drop"
    )
  
  # Method-specific analysis
  if (method == "lme") {
    # Linear Mixed-Effects Model with log transformation
    message("Fitting linear mixed-effects model with log transformation...")
    
    # Create a log-transformed version of the volume for handling exponential growth
    log_volume_col <- paste0("log_", volume_column)
    df[[log_volume_col]] <- log1p(df[[volume_column]])  # log1p = log(x+1) to handle zeros
    
    # Determine cage effect structure
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
    lme_model <- lme4::lmer(model_formula, data = df)
    print(summary(lme_model))
    
    # Type III ANOVA
    anova_results <- car::Anova(lme_model, type = 3)
    print(anova_results)
    
    # Post-hoc comparisons
    emm_formula <- as.formula(paste("pairwise ~ Group |", time_column))
    posthoc_results <- emmeans::emmeans(lme_model, specs = emm_formula, adjust = "bonferroni")
    print(posthoc_results)
    
    # Store results
    results$model <- lme_model
    results$anova <- anova_results
    results$posthoc <- posthoc_results
    
    # Model diagnostic plots
    residuals_df <- base::data.frame(Fitted = stats::fitted(lme_model), 
                                    Residuals = stats::residuals(lme_model))
    
    p_resid <- ggplot2::ggplot(residuals_df, ggplot2::aes(x = Fitted, y = Residuals)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_smooth(method = "loess", color = "blue") +
      ggplot2::labs(title = "Residuals vs. Fitted", x = "Fitted Values", y = "Residuals") +
      ggplot2::theme_minimal()
    
    p_qq <- ggplot2::ggplot(residuals_df, ggplot2::aes(sample = Residuals)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() +
      ggplot2::labs(title = "Q-Q Plot of Residuals") +
      ggplot2::theme_minimal()
    
    # Original scale plot
    p_raw <- ggplot2::ggplot(df, ggplot2::aes_string(x = time_column, y = volume_column, color = "Group")) +
      ggplot2::geom_line(ggplot2::aes(group = Mouse_ID), alpha = 0.3) +
      ggplot2::stat_summary(fun = base::mean, geom = "line", size = 1.2) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 2) +
      ggplot2::labs(title = "Tumor Growth Over Time (Original Scale)", 
                   x = "Days Since Injection", y = "Tumor Volume") +
      ggplot2::theme_minimal()
    
    # Log scale plot
    p_log <- ggplot2::ggplot(df, ggplot2::aes_string(x = time_column, y = log_volume_col, color = "Group")) +
      ggplot2::geom_line(ggplot2::aes(group = Mouse_ID), alpha = 0.3) +
      ggplot2::stat_summary(fun = base::mean, geom = "line", size = 1.2) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 2) +
      ggplot2::labs(title = "Log-Transformed Tumor Growth Over Time", 
                   x = "Days Since Injection", y = "Log(Tumor Volume + 1)") +
      ggplot2::theme_minimal()
    
    # Store plots
    results$plots <- list(
      raw_data = p_raw,
      transformed_data = p_log,
      residuals = p_resid,
      qq_plot = p_qq,
      combined = ggpubr::ggarrange(p_raw, p_log, p_resid, p_qq, 
                                 ncol = 2, nrow = 2, 
                                 labels = c("A", "B", "C", "D"))
    )
    
    # Check model diagnostics
    model_diagnostics <- performance::model_performance(lme_model)
    print(model_diagnostics)
    
  } else if (method == "gamm") {
    # Generalized Additive Mixed Model for flexible patterns
    message("Fitting generalized additive mixed model for flexible patterns...")
    
    # Determine the model formula based on collinearity
    if (!has_collinearity) {
      # Include cage as random effect
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
    gamm_model <- mgcv::gam(model_formula, 
                         data = df, 
                         method = "REML", 
                         select = TRUE,
                         gamma = gamma_val)
    
    print(summary(gamm_model))
    
    # Simulate ANOVA-like test
    anova_results <- mgcv::anova.gam(gamm_model, test = "F")
    print(anova_results)
    
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
    predictions <- mgcv::predict.gam(gamm_model, newdata = grid_data, se.fit = TRUE, type = "response")
    
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
      p_val <- 2 * pt(-abs(t_stat), df = gamm_model$df.residual)
      
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
    
    # Store results
    results$model <- gamm_model
    results$anova <- anova_results
    results$posthoc <- posthoc_results
    
    # Raw data plot
    p_raw <- ggplot2::ggplot(df, ggplot2::aes_string(x = time_column, y = volume_column, color = "Group")) +
      ggplot2::geom_line(ggplot2::aes(group = Mouse_ID), alpha = 0.2) +
      ggplot2::geom_point(data = summary_data, 
                       ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group"), 
                       size = 3) +
      ggplot2::geom_line(data = summary_data, 
                      ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group", group = "Group"), 
                      size = 1.2) +
      ggplot2::geom_errorbar(data = summary_data, 
                          ggplot2::aes_string(x = time_column, y = "Mean_Volume", 
                                           ymin = "Mean_Volume - SE", 
                                           ymax = "Mean_Volume + SE", 
                                           color = "Group"), 
                          width = 0.5) +
      ggplot2::labs(title = "Tumor Growth by Treatment Group",
                  x = "Days",
                  y = "Tumor Volume") +
      ggplot2::theme_minimal()
    
    # Model prediction plot
    p_pred <- ggplot2::ggplot() +
      ggplot2::geom_point(data = summary_data, 
                       ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group"), 
                       size = 3) +
      ggplot2::geom_line(data = grid_data, 
                      ggplot2::aes_string(x = time_column, y = "Predicted", color = "Group", group = "Group"), 
                      size = 1.2) +
      ggplot2::geom_ribbon(data = grid_data, 
                        ggplot2::aes_string(x = time_column, 
                                         ymin = "Predicted - 1.96 * SE", 
                                         ymax = "Predicted + 1.96 * SE", 
                                         fill = "Group"), 
                        alpha = 0.2) +
      ggplot2::labs(title = "GAMM Model Predictions",
                  x = "Days",
                  y = "Tumor Volume") +
      ggplot2::theme_minimal()
    
    # Residuals plot
    residuals_df <- data.frame(
      Fitted = fitted(gamm_model),
      Residuals = residuals(gamm_model)
    )
    
    p_resid <- ggplot2::ggplot(residuals_df, ggplot2::aes(x = Fitted, y = Residuals)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_smooth(method = "loess", color = "blue") +
      ggplot2::labs(title = "Residuals vs. Fitted", x = "Fitted Values", y = "Residuals") +
      ggplot2::theme_minimal()
    
    # QQ plot
    p_qq <- ggplot2::ggplot(residuals_df, ggplot2::aes(sample = Residuals)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() +
      ggplot2::labs(title = "Q-Q Plot of Residuals") +
      ggplot2::theme_minimal()
    
    # Store plots
    results$plots <- list(
      raw_data = p_raw,
      model_prediction = p_pred,
      residuals = p_resid,
      qq_plot = p_qq,
      combined = ggpubr::ggarrange(p_raw, p_pred, p_resid, p_qq, 
                                 ncol = 2, nrow = 2, 
                                 labels = c("A", "B", "C", "D"))
    )
    
  } else if (method == "auc") {
    # Area Under the Curve Analysis
    message("Performing Area Under the Curve (AUC) analysis...")
    
    # AUC already calculated above, now we perform statistical tests
    
    # ANOVA on AUC values
    auc_anova <- stats::aov(AUC ~ Group, data = auc_data)
    anova_results <- summary(auc_anova)
    print(anova_results)
    
    # Pairwise t-tests
    posthoc_results <- stats::pairwise.t.test(auc_data$AUC, auc_data$Group, 
                                            p.adjust.method = "bonferroni")
    print(posthoc_results)
    
    # Store results
    results$anova <- anova_results
    results$posthoc <- posthoc_results
    
    # Raw data plot
    p_raw <- ggplot2::ggplot(df, ggplot2::aes_string(x = time_column, y = volume_column, color = "Group")) +
      ggplot2::geom_line(ggplot2::aes(group = Mouse_ID), alpha = 0.2) +
      ggplot2::geom_point(data = summary_data, 
                       ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group"), 
                       size = 3) +
      ggplot2::geom_line(data = summary_data, 
                      ggplot2::aes_string(x = time_column, y = "Mean_Volume", color = "Group", group = "Group"), 
                      size = 1.2) +
      ggplot2::geom_errorbar(data = summary_data, 
                          ggplot2::aes_string(x = time_column, y = "Mean_Volume", 
                                           ymin = "Mean_Volume - SE", 
                                           ymax = "Mean_Volume + SE", 
                                           color = "Group"), 
                          width = 0.5) +
      ggplot2::labs(title = "Tumor Growth by Treatment Group",
                  x = "Days",
                  y = "Tumor Volume") +
      ggplot2::theme_minimal()
    
    # AUC barplot
    p_auc <- ggplot2::ggplot(auc_data, ggplot2::aes(x = Group, y = AUC, fill = Group)) +
      ggplot2::geom_bar(stat = "summary", fun = "mean") +
      ggplot2::geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
      ggplot2::labs(title = "Area Under the Curve by Treatment Group",
                  x = "Treatment Group",
                  y = "AUC (Tumor Burden)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    # Individual AUC values
    p_indiv <- ggplot2::ggplot(auc_data, ggplot2::aes(x = Group, y = AUC, color = Group)) +
      ggplot2::geom_boxplot(alpha = 0.3) +
      ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
      ggplot2::labs(title = "Individual AUC Values by Treatment Group",
                  x = "Treatment Group",
                  y = "AUC (Tumor Burden)") +
      ggplot2::theme_minimal()
    
    # AUC summary plot
    p_summary <- ggplot2::ggplot(auc_summary, ggplot2::aes(x = Group, y = AUC.Mean, fill = Group)) +
      ggplot2::geom_col() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = AUC.Mean - AUC.SD/sqrt(AUC.N), 
                                       ymax = AUC.Mean + AUC.SD/sqrt(AUC.N)), 
                          width = 0.2) +
      ggplot2::labs(title = "Summary of AUC by Treatment Group",
                  x = "Treatment Group",
                  y = "Mean AUC (± SEM)") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    # Store plots
    results$plots <- list(
      raw_data = p_raw,
      auc_barplot = p_auc,
      individual_auc = p_indiv,
      auc_summary = p_summary,
      combined = ggpubr::ggarrange(p_raw, p_auc, p_indiv, p_summary, 
                                 ncol = 2, nrow = 2, 
                                 labels = c("A", "B", "C", "D"))
    )
  }
  
  # Print combined plots before returning results
  print(results$plots$combined)
  
  # Return complete results list
  return(results)
}