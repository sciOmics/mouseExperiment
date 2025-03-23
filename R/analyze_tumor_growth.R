#' Analyze Tumor Growth Using Linear Mixed Effects Model
#'
#' This function fits a linear mixed-effects model (LME) to assess statistical significance
#' in tumor growth data over time, accounting for repeated measures within subjects.
#' It also generates diagnostic plots and post-hoc comparisons.
#'
#' @param df A data frame containing tumor growth data.
#' @param time_column A character string specifying the column name for time points (e.g., "Day"). Default is "Day".
#' @param volume_column A character string specifying the column name for tumor volume measurements (e.g., "Volume"). Default is "Volume".
#' @param group_column A character string specifying the column name for treatment groups (e.g., "Group"). Default is "Group".
#' @param id_column A character string specifying the column name for individual subject identifiers (e.g., "ID"). Default is "ID".
#'
#' @return A list containing:
#' \item{model}{The fitted linear mixed-effects model.}
#' \item{anova}{Type III ANOVA table for fixed effects.}
#' \item{posthoc}{Pairwise comparisons of estimated marginal means (adjusted for multiple comparisons).}
#' \item{plots}{A set of diagnostic and tumor growth visualization plots.}
#'
#' @details
#' - The function creates a unique identifier by combining `Group` and `ID` to account for repeated measures.
#' - A linear mixed-effects model (`lmer`) is used, with `Mouse_ID` as a random effect.
#' - The function performs a Type III ANOVA to assess the significance of time, treatment group, and their interaction.
#' - Post-hoc pairwise comparisons are performed using `emmeans`, adjusted with the Bonferroni correction.
#' - The function generates diagnostic plots to check model assumptions, including:
#'   - Tumor growth over time with individual trajectories and group means.
#'   - Residuals vs. fitted values for checking homoscedasticity.
#'   - Q-Q plot for assessing normality of residuals.
#'   - Model assumption checks using `check_model()`.
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

analyze_tumor_growth <- function(df, time_column = "Day", volume_column = "Volume", group_column = "Group", id_column = "ID") {
  # Ensure required columns exist
  if (!all(c(time_column, volume_column, group_column, id_column) %in% base::colnames(df))) {
    stop("One or more specified columns are not found in the dataframe.")
  }

  # Create a unique identifier for each mouse
  df$Mouse_ID <- base::factor(paste(df[[group_column]], df[[id_column]], sep = "_"))

  # Fit a Linear Mixed Effects Model
  lme_model <- lme4::lmer(as.formula(base::paste(volume_column, "~", time_column, "*", group_column, "+ (1 | Mouse_ID)")), data = df)

  # Model summary
  print(summary(lme_model))

  # Type III ANOVA to assess significance of fixed effects
  anova_results <- car::Anova(lme_model, type = 3)
  print(anova_results)

  # Post-hoc comparisons (estimated marginal means)
  emmeans_results <- emmeans::emmeans(lme_model, pairwise ~ df[[group_column]] | df[[time_column]], adjust = "bonferroni")
  print(emmeans_results)

  # ----- Visualization -----

  # 1. Tumor Growth Over Time
  p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x = time_column, y = volume_column, color = group_column)) +
    ggplot2::geom_line(ggplot2::aes(group = Mouse_ID), alpha = 0.3) +  # Individual mice
    ggplot2::stat_summary(fun = base::mean, geom = "line", size = 1.2) +  # Group mean
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 2) +  # Error bars
    ggplot2::labs(title = "Tumor Growth Over Time", x = "Days Since Injection", y = "Tumor Volume") +
    ggplot2::theme_minimal()

  # 2. Residual Diagnostics for Model Fit
  residuals_df <- base::data.frame(Fitted = stats::fitted(lme_model), Residuals = stats::residuals(lme_model))
  p2 <- ggplot2::ggplot(residuals_df, ggplot2::aes(x = Fitted, y = Residuals)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(method = "loess", color = "blue") +
    ggplot2::labs(title = "Residuals vs. Fitted", x = "Fitted Values", y = "Residuals") +
    ggplot2::theme_minimal()

  # 3. Q-Q Plot for Normality of Residuals
  p3 <- ggplot2::ggplot(residuals_df, ggplot2::aes(sample = Residuals)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line() +
    ggplot2::labs(title = "Q-Q Plot of Residuals") +
    ggplot2::theme_minimal()

  # Arrange plots in a grid
  plot_grid <- cowplot::ggarrange(p1, p2, p3, ncol = 2, nrow = 2, labels = c("A", "B", "C"))

  # Print plots
  print(plot_grid)

  # Check model assumptions
  print(performance::check_model(lme_model))

  return(list(model = lme_model, anova = anova_results, posthoc = emmeans_results, plots = plot_grid))
}

