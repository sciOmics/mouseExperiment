#' Create Caterpillar Plot for Frequentist Models
#'
#' This function creates a caterpillar plot (coefficient plot) for linear mixed-effects models,
#' displaying fixed effects coefficients with confidence intervals.
#'
#' @param model A linear mixed-effects model object from lme4::lmer() or a similar function
#' @param title A title for the plot (default: "Fixed Effects with 95% Confidence Intervals")
#' @param colors A vector of colors for different coefficient groups (default: blues and reds)
#' @param show_intercept Logical, whether to include the intercept in the plot (default: TRUE)
#' @param ci_level Confidence level for intervals (default: 0.95 for 95% CI)
#'
#' @return A ggplot object representing the caterpillar plot
#' 
#' @details
#' The function extracts fixed effects coefficients and their confidence intervals from
#' linear mixed-effects models. The output is a caterpillar plot with the following features:
#' - Point estimates for each coefficient
#' - Error bars showing confidence intervals
#' - Vertical line at zero for reference
#' - Optional grouping by coefficient type
#'
#' @examples
#' \dontrun{
#' # After running tumor growth statistics
#' results <- tumor_growth_statistics(df)
#' 
#' # Create caterpillar plot from the model
#' cat_plot <- plot_caterpillar(results$model)
#' print(cat_plot)
#' 
#' # With customizations
#' cat_plot2 <- plot_caterpillar(results$model, 
#'                              title = "Treatment Effects on Tumor Growth",
#'                              show_intercept = FALSE,
#'                              ci_level = 0.9)
#' }
#'
#' @import ggplot2
#' @export
plot_caterpillar <- function(model, title = "Fixed Effects with 95% Confidence Intervals", 
                            colors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"), 
                            show_intercept = TRUE, ci_level = 0.95) {
  
  # Check if the input is a valid model
  if (!inherits(model, c("lmerMod", "lm", "glm", "glmerMod"))) {
    stop("Model must be a linear or generalized linear mixed-effects model (lmerMod, lm, glm, glmerMod)")
  }
  
  # Get model summary
  mod_summary <- summary(model)
  
  # Extract fixed effects coefficients and standard errors
  coefs <- stats::coef(mod_summary)
  
  if ("Estimate" %in% colnames(coefs)) {
    # For lmer models
    coef_names <- rownames(coefs)
    estimates <- coefs[, "Estimate"]
    std_errors <- coefs[, "Std. Error"]
  } else {
    # For other model types
    coef_names <- names(stats::coef(model))
    estimates <- stats::coef(model)
    std_errors <- sqrt(diag(stats::vcov(model)))
  }
  
  # Calculate confidence intervals using t-distribution
  # Extract degrees of freedom based on model type
  if (inherits(model, "lm") && !inherits(model, c("lmerMod", "glmerMod"))) {
    df_resid <- stats::df.residual(model)
  } else {
    # For mixed models, approximate df as n - p
    df_resid <- stats::nobs(model) - length(estimates)
  }
  t_value <- stats::qt(1 - (1 - ci_level) / 2, df = df_resid)
  ci_lower <- estimates - t_value * std_errors
  ci_upper <- estimates + t_value * std_errors
  
  # Remove intercept if requested
  if (!show_intercept) {
    intercept_idx <- grep("^\\(Intercept\\)$", coef_names)
    if (length(intercept_idx) > 0) {
      coef_names <- coef_names[-intercept_idx]
      estimates <- estimates[-intercept_idx]
      ci_lower <- ci_lower[-intercept_idx]
      ci_upper <- ci_upper[-intercept_idx]
    }
  }
  
  # Create effect type groups
  effect_type <- rep("Other", length(coef_names))
  
  # Main effects (no interaction terms)
  main_effect_pattern <- "^[^:]*$"
  effect_type[grepl(main_effect_pattern, coef_names) & !grepl("^\\(Intercept\\)$", coef_names)] <- "Main Effect"
  
  # Interaction effects
  interaction_pattern <- ":"
  effect_type[grep(interaction_pattern, coef_names)] <- "Interaction"
  
  # Intercept
  effect_type[grep("^\\(Intercept\\)$", coef_names)] <- "Intercept"
  
  # Prepare data for plotting
  plot_data <- data.frame(
    term = factor(coef_names, levels = rev(coef_names)),  # Reverse for better display
    estimate = estimates,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    effect_type = effect_type
  )
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = estimate, y = term, color = effect_type)) +
    # Add vertical reference line at zero
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Add confidence intervals
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
    # Add point estimates
    ggplot2::geom_point(size = 3) +
    # Labels and title
    ggplot2::labs(
      title = title,
      x = "Coefficient Estimate",
      y = NULL,
      color = "Effect Type"
    ) +
    # Set colors for different effect types
    ggplot2::scale_color_manual(values = colors) +
    # Theme customizations
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      axis.text.y = ggplot2::element_text(hjust = 0),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}

