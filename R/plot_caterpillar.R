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

#' Create Forest Plot for Hazard Ratios or Effect Sizes
#'
#' @description
#' Creates a forest plot showing hazard ratios or effect sizes with confidence intervals.
#' Points represent estimates, and horizontal lines represent 95% confidence intervals.
#' This function works with both survival analysis results (hazard ratios) and
#' tumor growth statistics (effect sizes).
#'
#' @param results_df Data frame containing results (hazard ratios or effect sizes).
#' @param group_order Optional vector specifying the order in which groups should appear (from top to bottom).
#' @param title Optional title for the plot. Use "Forest Plot of Effect Sizes" for tumor growth data
#'        or "Forest Plot of Hazard Ratios" (default) for survival data.
#' @param show_text Logical, whether to display values with confidence intervals and p-values below group names. Default: TRUE.
#'
#' @return A ggplot2 object with the forest plot visualization.
#' @import ggplot2
#' @export
#'
#' @examples
#' # First run survival analysis
#' results <- survival_statistics(df = tumor_data)
#' 
#' # Then create forest plot with default ordering
#' forest_plot(results$results)
#' 
#' # Custom group ordering
#' forest_plot(results$results, 
#'            group_order = c("Control", "Drug A", "Drug B"))
forest_plot <- function(results_df, group_order = NULL, title = "Forest Plot of Hazard Ratios", show_text = TRUE) {
  if (nrow(results_df) == 0 || all(is.na(results_df$Hazard_Ratio))) {
    warning("Cannot create forest plot: No hazard ratio data available")
    return(NULL)
  }
  
  # Prepare plot data
  plot_data <- results_df
  ref_group <- results_df$Group[results_df$Note == "Reference group"]
  max_plot_hr <- 20
  
  # Reorder groups if specified
  if (!is.null(group_order)) {
    # Handle groups with dose information
    data_groups <- unique(plot_data$Group)
    
    # Check if we have groups with "Dose:" in them
    has_dose_info <- any(grepl("Dose:", data_groups))
    
    if (has_dose_info) {
      # Extract base treatment names and handle dose groups
      expanded_order <- c()
      for (base_group in group_order) {
        # Add the base group if it exists exactly
        if (base_group %in% data_groups) {
          expanded_order <- c(expanded_order, base_group)
        }
        
        # Find and add any groups starting with this base name followed by " - Dose:"
        dose_groups <- grep(paste0("^", base_group, " - Dose:"), data_groups, value = TRUE)
        if (length(dose_groups) > 0) {
          # Sort dose groups by dose level if possible
          dose_values <- as.numeric(gsub(paste0("^", base_group, " - Dose: "), "", dose_groups))
          if (!any(is.na(dose_values))) {
            # Sort by dose value if all could be converted to numbers
            dose_groups <- dose_groups[order(dose_values)]
          }
          expanded_order <- c(expanded_order, dose_groups)
        }
      }
      
      # Find any remaining groups
      missing_groups <- setdiff(data_groups, expanded_order)
      if (length(missing_groups) > 0) {
        expanded_order <- c(expanded_order, missing_groups)
      }
      
      # Use the expanded order
      group_order <- expanded_order
    } else {
      # For data without dose info, use original approach
      # Check if all specified groups exist in the data
      missing_groups <- setdiff(group_order, data_groups)
      if (length(missing_groups) > 0) {
        warning("Some specified groups not found in data: ", 
                paste(missing_groups, collapse = ", "))
      }
      
      # Ensure all data groups are included in the ordering
      missing_in_order <- setdiff(data_groups, group_order)
      if (length(missing_in_order) > 0) {
        group_order <- c(group_order, missing_in_order)
      }
    }
    
    # Convert Group to factor with specified levels for ordering
    # We reverse the order for plotting (top to bottom)
    plot_data$Group <- factor(plot_data$Group, levels = rev(group_order))
    # Store the original order (without reversal) for text labels
    original_group_order <- group_order
  }
  
  # Add plotting columns
  plot_data$CI_Upper_Plot <- NA
  plot_data$CI_Lower_Plot <- NA
  plot_data$HR_Plot <- NA
  plot_data$Note_Plot <- ""
  
  # Process each row for plotting
  for (i in 1:nrow(plot_data)) {
    if (plot_data$Group[i] != ref_group) {
      # Non-reference groups
      if (is.na(plot_data$CI_Upper[i]) || plot_data$CI_Upper[i] > max_plot_hr) {
        plot_data$CI_Upper_Plot[i] <- max_plot_hr
        plot_data$Note_Plot[i] <- "(CI extends beyond chart)"
      } else {
        plot_data$CI_Upper_Plot[i] <- plot_data$CI_Upper[i]
      }
      
      if (is.na(plot_data$CI_Lower[i]) || plot_data$CI_Lower[i] < 0.05) {
        plot_data$CI_Lower_Plot[i] <- 0.05
      } else {
        plot_data$CI_Lower_Plot[i] <- plot_data$CI_Lower[i]
      }
      
      if (is.na(plot_data$Hazard_Ratio[i]) || plot_data$Hazard_Ratio[i] > max_plot_hr) {
        plot_data$HR_Plot[i] <- max_plot_hr
      } else if (plot_data$Hazard_Ratio[i] < 0.05) {
        plot_data$HR_Plot[i] <- 0.05
      } else {
        plot_data$HR_Plot[i] <- plot_data$Hazard_Ratio[i]
      }
    } else {
      # Reference group
      plot_data$CI_Lower_Plot[i] <- 1
      plot_data$CI_Upper_Plot[i] <- 1
      plot_data$HR_Plot[i] <- 1
    }
  }
  
  # Create labels without [Events/Total]
  # Detect if this is an effect size plot (from tumor LME) or hazard ratio plot (from survival)
  using_effect_size <- any(grepl("effect", tolower(title)))
  
  # Adjust label prefix based on the type of plot
  ratio_prefix <- ifelse(using_effect_size, "ES: ", "HR: ")
  
  plot_data$label <- ifelse(plot_data$Group == ref_group,
                          sprintf("%s1.00 (Reference)", ratio_prefix), 
                          ifelse(is.na(plot_data$P_Value),
                                sprintf("%s%.2f (%.2f-%.2f)", 
                                       ratio_prefix,
                                       plot_data$Hazard_Ratio, 
                                       ifelse(is.na(plot_data$CI_Lower), 0, plot_data$CI_Lower), 
                                       ifelse(is.na(plot_data$CI_Upper), Inf, plot_data$CI_Upper)),
                                sprintf("%s%.2f (%.2f-%.2f), p=%.3f", 
                                       ratio_prefix,
                                       plot_data$Hazard_Ratio, 
                                       ifelse(is.na(plot_data$CI_Lower), 0, plot_data$CI_Lower), 
                                       ifelse(is.na(plot_data$CI_Upper), Inf, plot_data$CI_Upper), 
                                       plot_data$P_Value)))
  
  # Create the forest plot
  forest_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = HR_Plot, y = Group)) +
    # Reference line at HR = 1
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    # Confidence intervals
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_Lower_Plot, xmax = CI_Upper_Plot), height = 0.2) +
    # Point estimates (black instead of blue)
    ggplot2::geom_point(size = 3, color = "black") +
    # Labels
    ggplot2::labs(
      title = title,
      x = ifelse(any(grepl("effect", tolower(title))), 
                "Effect Size (95% CI)", 
                "Hazard Ratio (95% CI)"),
      y = NULL
    ) +
    # Log scale
    ggplot2::scale_x_log10(limits = c(0.05, max_plot_hr),
                         breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
                         labels = c("0.1", "0.2", "0.5", "1", "2", "5", "10", "20+")) +
    # Theme settings 
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = if(show_text) 30 else 10)
    )
  
  # Add HR, CI, and p-value text labels if requested
  if (show_text) {
    # Use the correct group levels for positioning
    # group_levels should be in reverse order to match the plot
    group_levels <- levels(plot_data$Group)
    
    # Check if we have any group levels before creating the data frame
    if (length(group_levels) > 0) {
      # Use a separate data frame for the labels
      # The y_pos should match the position on the plot (from top to bottom)
      label_data <- data.frame(
        y_pos = seq_along(group_levels),  # Positions from 1 to n
        group_name = group_levels,        # These are already in reversed order
        stringsAsFactors = FALSE
      )
    
      # Match the labels with the corresponding groups
      for (i in 1:nrow(label_data)) {
        group_name <- label_data$group_name[i]
        matching_row <- which(plot_data$Group == group_name)
        if (length(matching_row) > 0) {
          label_data$label[i] <- plot_data$label[matching_row[1]]
        } else {
          label_data$label[i] <- ""
        }
      }
      
      # Modify the plot with a better approach - put labels on y-axis
      # Create a named vector for new y-axis labels
      new_labels <- setNames(
        paste0(group_levels, "\n", label_data$label), 
        group_levels
      )
      
      # Apply the new labels to the y-axis
      forest_plot <- forest_plot + 
        ggplot2::scale_y_discrete(labels = new_labels) +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(lineheight = 0.9),
          plot.margin = ggplot2::margin(l = 5, r = 10, t = 10, b = 10)
        ) +
        ggplot2::coord_cartesian(clip = "off")
    }
  }
  
  # Add notes for truncated CIs
  # Get the numeric positions for each group
  group_levels <- levels(plot_data$Group)
  if (length(group_levels) > 0) {
    group_pos <- data.frame(
      group_name = group_levels,
      pos = seq_along(group_levels),  # Positions should match plot (1 to n from top to bottom)
      stringsAsFactors = FALSE
    )
  } else {
    # Create an empty data frame if no group levels
    group_pos <- data.frame(
      group_name = character(0),
      pos = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  for (i in 1:nrow(plot_data)) {
    if (plot_data$Note_Plot[i] != "") {
      # Find the numeric position for this group
      group_name <- as.character(plot_data$Group[i])
      pos <- group_pos$pos[group_pos$group_name == group_name]
      
      forest_plot <- forest_plot +
        ggplot2::annotate("text", x = max_plot_hr, y = pos, 
                        label = plot_data$Note_Plot[i], hjust = 1, vjust = -1, 
                        color = "red", size = 3)
    }
  }
  
  return(forest_plot)
}

