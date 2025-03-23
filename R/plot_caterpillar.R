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
  
  # Calculate confidence intervals
  z_value <- stats::qnorm(1 - (1 - ci_level) / 2)
  ci_lower <- estimates - z_value * std_errors
  ci_upper <- estimates + z_value * std_errors
  
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
  effect_type[grep(main_effect_pattern, coef_names) & !grepl("^\\(Intercept\\)$", coef_names)] <- "Main Effect"
  
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

#' Create Caterpillar Plot for Bayesian Models
#'
#' This function creates a caterpillar plot (coefficient plot) for Bayesian models 
#' from the brms package, displaying posterior distributions of model parameters.
#'
#' @param model A brmsfit object from brms::brm()
#' @param title A title for the plot (default: "Posterior Distributions with 95% Credible Intervals")
#' @param colors A vector of colors for different parameter groups (default: blues and reds)
#' @param show_intercept Logical, whether to include the intercept in the plot (default: TRUE)
#' @param ci_level Credible interval level (default: 0.95 for 95% CI)
#' @param show_random Logical, whether to show group-level (random) effects (default: FALSE)
#'
#' @return A ggplot object representing the caterpillar plot
#' 
#' @details
#' The function extracts posterior distributions from Bayesian models fit with brms.
#' The output is a caterpillar plot with the following features:
#' - Posterior median for each parameter
#' - Error bars showing credible intervals
#' - Vertical line at zero for reference
#' - Optional grouping by parameter type
#'
#' @examples
#' \dontrun{
#' # After running Bayesian tumor growth statistics
#' bayes_results <- tumor_growth_statistics_bayes(df)
#' 
#' # Create caterpillar plot from the Bayesian model
#' bayes_cat_plot <- plot_caterpillar_bayes(bayes_results$model)
#' print(bayes_cat_plot)
#' 
#' # With customizations
#' bayes_cat_plot2 <- plot_caterpillar_bayes(bayes_results$model, 
#'                                         title = "Treatment Effects on Tumor Growth",
#'                                         show_intercept = FALSE,
#'                                         ci_level = 0.89)
#' }
#'
#' @import ggplot2
#' @importFrom brms fixef 
#' @export
plot_caterpillar_bayes <- function(model, title = "Posterior Distributions with 95% Credible Intervals", 
                                  colors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"), 
                                  show_intercept = TRUE, ci_level = 0.95, show_random = FALSE) {
  
  # Check if the input is a valid brmsfit model
  if (!inherits(model, "brmsfit")) {
    stop("Model must be a brmsfit object from brms package")
  }
  
  # Calculate probabilities for credible intervals
  lower_prob <- (1 - ci_level) / 2
  upper_prob <- 1 - lower_prob
  
  # Extract fixed effects (population-level effects in brms terminology)
  fixed_effects <- brms::fixef(model, probs = c(lower_prob, 0.5, upper_prob))
  
  param_names <- rownames(fixed_effects)
  estimates <- fixed_effects[, "Estimate"]
  ci_lower <- fixed_effects[, paste0(lower_prob * 100, "%")]
  ci_upper <- fixed_effects[, paste0(upper_prob * 100, "%")]
  
  # Remove intercept if requested
  if (!show_intercept) {
    intercept_idx <- grep("^Intercept$", param_names)
    if (length(intercept_idx) > 0) {
      param_names <- param_names[-intercept_idx]
      estimates <- estimates[-intercept_idx]
      ci_lower <- ci_lower[-intercept_idx]
      ci_upper <- ci_upper[-intercept_idx]
    }
  }
  
  # Create parameter type groups
  param_type <- rep("Other", length(param_names))
  
  # Main effects (no interaction terms)
  param_type[grep("^Intercept$", param_names)] <- "Intercept"
  
  # Main effects (no interaction terms, not intercept)
  main_effect_pattern <- "^[^:]*$"
  param_type[grep(main_effect_pattern, param_names) & !grepl("^Intercept$", param_names)] <- "Main Effect"
  
  # Interaction effects
  interaction_pattern <- ":"
  param_type[grep(interaction_pattern, param_names)] <- "Interaction"
  
  # Prepare data for plotting
  plot_data <- data.frame(
    term = factor(param_names, levels = rev(param_names)),  # Reverse for better display
    estimate = estimates,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    param_type = param_type
  )
  
  # Add group-level effects (random effects) if requested
  if (show_random) {
    # Check if there are random effects in the model
    if (!is.null(model$ranef)) {
      # Get random effects summary
      ranef_summary <- brms::ranef(model, probs = c(lower_prob, 0.5, upper_prob))
      
      # Process each random effects group
      for (group in names(ranef_summary)) {
        for (effect in dimnames(ranef_summary[[group]])[[2]]) {
          for (level in dimnames(ranef_summary[[group]])[[1]]) {
            # Extract values
            est <- ranef_summary[[group]][level, effect, "Estimate"]
            low <- ranef_summary[[group]][level, effect, paste0(lower_prob * 100, "%")]
            upp <- ranef_summary[[group]][level, effect, paste0(upper_prob * 100, "%")]
            
            # Add to plot data
            term_name <- paste0(group, "_", effect, "_", level)
            plot_data <- rbind(plot_data, data.frame(
              term = term_name,
              estimate = est,
              ci_lower = low,
              ci_upper = upp,
              param_type = "Random Effect"
            ))
          }
        }
      }
      
      # Re-factor terms to maintain order
      plot_data$term <- factor(plot_data$term, levels = rev(as.character(plot_data$term)))
    }
  }
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = estimate, y = term, color = param_type)) +
    # Add vertical reference line at zero
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Add credible intervals
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
    # Add point estimates (posterior medians)
    ggplot2::geom_point(size = 3) +
    # Labels and title
    ggplot2::labs(
      title = title,
      x = "Parameter Estimate",
      y = NULL,
      color = "Parameter Type"
    ) +
    # Set colors for different parameter types
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