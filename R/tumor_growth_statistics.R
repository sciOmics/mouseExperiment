#' Analyze Tumor Growth Using Linear Mixed Effects Model with Log Transformation
#'
#' This function fits a linear mixed-effects model (LME) to assess statistical significance
#' in tumor growth data over time, accounting for repeated measures within subjects.
#' It applies a log transformation to tumor volume data to handle the exponential nature
#' of tumor growth. The function also generates diagnostic plots and post-hoc comparisons.
#'
#' @param df A data frame containing tumor growth data.
#' @param time_column A character string specifying the column name for time points (e.g., "Day"). Default is "Day".
#' @param volume_column A character string specifying the column name for tumor volume measurements (e.g., "Volume"). Default is "Volume".
#' @param treatment_column A character string specifying the column name for treatment groups (e.g., "Treatment"). Default is "Treatment".
#' @param cage_column A character string specifying the column name for the cage identifier (e.g., "Cage"). Default is "Cage". The function will test for collinearity between Cage and Treatment variables. If there's no collinearity, Cage will be included as a nested random effect within Treatment in the model.
#' @param id_column A character string specifying the column name for individual subject identifiers (e.g., "ID"). Default is "ID".
#' @param dose_column Optional. A character string specifying the column name for dose levels, if available. Default is NULL.
#'
#' @return A list containing:
#' \item{model}{The fitted linear mixed-effects model (using log-transformed volume).}
#' \item{anova}{Type III ANOVA table for fixed effects.}
#' \item{posthoc}{Pairwise comparisons of estimated marginal means (adjusted for multiple comparisons).}
#' \item{plots}{A set of diagnostic and tumor growth visualization plots.}
#'
#' @details
#' - The function applies a log transformation (log(x+1)) to tumor volume data to account for exponential growth patterns.
#' - The function creates a unique identifier by combining `Group` and `ID` to account for repeated measures.
#' - A linear mixed-effects model (`lmer`) is used, with `Mouse_ID` as a random effect.
#' - The function handles cage effects in a sophisticated manner:
#'   - It prints a contingency table showing the distribution of cages across treatment groups.
#'   - If collinearity exists (e.g., each cage contains only one treatment group), only Treatment effect is used.
#'   - If no collinearity is detected, the function checks if there are multiple cages per treatment group:
#'     - With multiple cages per treatment: Uses explicit nested random effects with `(1 | Cage_in_Group)`,
#'       where `Cage_in_Group` is created using the `interaction()` function.
#'     - Otherwise: Uses an additive random effect for Cage with `(1 | Cage_Var)`.
#'   - This approach properly accounts for the hierarchical structure of the data.
#'   - The exact model formula used is reported as a message when the function runs.
#' - The function performs a Type III ANOVA to assess the significance of time, treatment group, and their interaction.
#' - Post-hoc pairwise comparisons are performed using `emmeans`, adjusted with the Bonferroni correction.
#' - The function generates diagnostic plots to check model assumptions, including:
#'   - Tumor growth over time with individual trajectories and group means (on original scale).
#'   - Log-transformed tumor growth over time (to visualize the transformed data).
#'   - Residuals vs. fitted values for checking homoscedasticity.
#'   - Q-Q plot for assessing normality of residuals.
#'   - Model assumption checks.
#'
#' @examples
#' # Example usage:
#' # analyze_tumor_growth(df)
#'
#' @import lme4
#' @import emmeans
#' @import ggplot2
#' @import ggpubr
#' @import performance
#' @export

tumor_growth_statistics <- function(df, time_column = "Day", volume_column = "Volume", 
                             treatment_column = "Treatment", cage_column = "Cage", id_column = "ID", 
                             dose_column = NULL) {
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

  # Create a log-transformed version of the volume for handling exponential growth
  log_volume_col <- paste0("log_", volume_column)
  df[[log_volume_col]] <- log1p(df[[volume_column]])  # log1p = log(x+1) to handle zeros
  
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
      message("No perfect collinearity detected between Cage and Treatment. Including Cage as a nested variable.")
    }
  } else {
    message("Only one cage detected. Using Treatment only in the model.")
    has_collinearity <- TRUE  # Treat single cage as collinear case
  }
  
  # First, determine if the data has the right structure for nested random effects
  # We need to verify there are multiple cages per treatment group
  df$Cage_Var <- factor(df$Cage_Var)
  df$Group <- factor(df$Group)
  
  # Create a contingency table to check cage distribution across groups
  cage_group_table <- table(df$Cage_Var, df$Group)
  print("Cage distribution across treatment groups:")
  print(cage_group_table)
  
  # Create nested cage variable if needed
  if (!has_collinearity) {
    # Check if we have a proper nesting structure
    # We need multiple cages per treatment group for proper nesting
    cages_per_group <- colSums(cage_group_table > 0)
    
    if (all(cages_per_group > 1)) {
      # We have proper nesting - multiple cages per treatment
      # Use the hierarchical structure with nested random effects
      
      # Create the nested cage identifier
      df$Cage_in_Group <- interaction(df$Group, df$Cage_Var, sep = ":")
      
      # Fit a Linear Mixed Effects Model with explicit nested random effects
      lme_formula <- as.formula(paste(
        log_volume_col, "~", 
        time_column, "* Group + (1 | Mouse_ID) + (1 | Cage_in_Group)"
      ))
      
      message("Using model with explicit nested random effects: ", deparse(lme_formula))
      lme_model <- lme4::lmer(lme_formula, data = df)
      
    } else {
      # We don't have proper nesting structure, but no collinearity
      # Include cage as an additive random effect
      
      lme_formula <- as.formula(paste(
        log_volume_col, "~", 
        time_column, "* Group + (1 | Mouse_ID) + (1 | Cage_Var)"
      ))
      
      message("Using model with additive cage random effect: ", deparse(lme_formula))
      lme_model <- lme4::lmer(lme_formula, data = df)
    }
  } else {
    # Fit a standard Linear Mixed Effects Model without cage effect
    lme_formula <- as.formula(paste(
      log_volume_col, "~", 
      time_column, "* Group + (1 | Mouse_ID)"
    ))
    
    message("Using model without cage effect: ", deparse(lme_formula))
    lme_model <- lme4::lmer(lme_formula, data = df)
  }

  # Model summary
  print(summary(lme_model))

  # Type III ANOVA to assess significance of fixed effects
  anova_results <- car::Anova(lme_model, type = 3)
  print(anova_results)

  # Post-hoc comparisons (estimated marginal means)
  # Create formula for emmeans that uses the constructed Group variable
  emm_formula <- as.formula(paste("pairwise ~ Group |", time_column))
  emmeans_results <- emmeans::emmeans(lme_model, specs = emm_formula, adjust = "bonferroni")
  print(emmeans_results)

  # ----- Visualization -----

  # 1. Tumor Growth Over Time (Original Scale)
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x = time_column, y = volume_column, color = "Group")) +
    ggplot2::geom_line(ggplot2::aes(group = Mouse_ID), alpha = 0.3) +  # Individual mice
    ggplot2::stat_summary(fun = base::mean, geom = "line", size = 1.2) +  # Group mean
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 2) +  # Error bars
    ggplot2::labs(title = "Tumor Growth Over Time (Original Scale)", 
                 x = "Days Since Injection", y = "Tumor Volume") +
    ggplot2::theme_minimal()

  # 2. Log-Transformed Tumor Growth Over Time
  p2 <- ggplot2::ggplot(df, ggplot2::aes_string(x = time_column, y = log_volume_col, color = "Group")) +
    ggplot2::geom_line(ggplot2::aes(group = Mouse_ID), alpha = 0.3) +  # Individual mice
    ggplot2::stat_summary(fun = base::mean, geom = "line", size = 1.2) +  # Group mean
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 2) +  # Error bars
    ggplot2::labs(title = "Log-Transformed Tumor Growth Over Time", 
                 x = "Days Since Injection", y = "Log(Tumor Volume + 1)") +
    ggplot2::theme_minimal()

  # 3. Residual Diagnostics for Model Fit
  residuals_df <- base::data.frame(Fitted = stats::fitted(lme_model), Residuals = stats::residuals(lme_model))
  p3 <- ggplot2::ggplot(residuals_df, ggplot2::aes(x = Fitted, y = Residuals)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(method = "loess", color = "blue") +
    ggplot2::labs(title = "Residuals vs. Fitted", x = "Fitted Values", y = "Residuals") +
    ggplot2::theme_minimal()

  # 4. Q-Q Plot for Normality of Residuals
  p4 <- ggplot2::ggplot(residuals_df, ggplot2::aes(sample = Residuals)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line() +
    ggplot2::labs(title = "Q-Q Plot of Residuals") +
    ggplot2::theme_minimal()

  # Arrange plots in a grid
  plot_grid <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))

  # Print plots
  print(plot_grid)

  # Check model diagnostics - basic normality and residual info without requiring the 'see' package
  # Instead of using check_model which requires the 'see' package
  model_diagnostics <- performance::model_performance(lme_model)
  print(model_diagnostics)

  return(list(model = lme_model, anova = anova_results, posthoc = emmeans_results, plots = plot_grid))
}

