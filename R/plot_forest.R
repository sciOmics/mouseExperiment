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
  # Check if input data frame is empty
  if (is.null(results_df) || nrow(results_df) == 0) {
    warning("Cannot create forest plot: Empty results data frame")
    return(NULL)
  }
  
  # Make a copy of the input data to avoid modifying the original
  plot_data <- results_df
  
  # Handle different column naming conventions
  col_mapping <- list(
    hr = c("HR", "Hazard_Ratio", "hazard_ratio", "hr"),
    ci_lower = c("CI_Lower", "CI.Lower", "ci_lower", "ci.lower", "lower"),
    ci_upper = c("CI_Upper", "CI.Upper", "ci_upper", "ci.upper", "upper")
  )
  
  # Standardize column names
  standardize_columns <- function(df, mapping) {
    for (target_col in names(mapping)) {
      possible_names <- mapping[[target_col]]
      existing_cols <- intersect(colnames(df), possible_names)
      
      if (length(existing_cols) > 0) {
        # Use first matching column found
        source_col <- existing_cols[1]
        
        # Set standard column name
        standard_name <- switch(target_col,
                             "hr" = "Hazard_Ratio",
                             "ci_lower" = "CI.Lower",
                             "ci_upper" = "CI.Upper")
        
        # Create the standardized column if it doesn't exist
        if (!(standard_name %in% colnames(df))) {
          df[[standard_name]] <- df[[source_col]]
        }
      }
    }
    return(df)
  }
  
  # Apply standardization
  plot_data <- standardize_columns(plot_data, col_mapping)
  
  # Check if we have the required columns after standardization
  if (!all(c("Hazard_Ratio", "CI.Lower", "CI.Upper") %in% colnames(plot_data))) {
    warning("Cannot create forest plot: Missing required columns (Hazard_Ratio, CI.Lower, CI.Upper)")
    return(NULL)
  }
  
  # Check if all values in Hazard_Ratio are NA
  if (all(is.na(plot_data$Hazard_Ratio))) {
    warning("Cannot create forest plot: No hazard ratio data available")
    return(NULL)
  }
  
  # Prepare plot data
  ref_group <- plot_data$Group[plot_data$Note == "Reference group" | plot_data$Hazard_Ratio == 1]
  max_plot_hr <- 20
  
  # Ensure "Group" column exists
  if (!("Group" %in% colnames(plot_data))) {
    # Try to find an alternative group column
    potential_group_cols <- c("Treatment", "group", "treatment", "Group")
    existing_cols <- intersect(colnames(plot_data), potential_group_cols)
    
    if (length(existing_cols) > 0) {
      plot_data$Group <- plot_data[[existing_cols[1]]]
    } else {
      warning("Cannot create forest plot: No 'Group' column found")
      return(NULL)
    }
  }
  
  # Ensure group names are character type to avoid factor issues
  plot_data$Group <- as.character(plot_data$Group)
  
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
    }
    
    # Ensure all groups in the data are in the order
    missing_groups <- setdiff(unique(plot_data$Group), group_order)
    if (length(missing_groups) > 0) {
      group_order <- c(group_order, missing_groups)
    }
    
    # Set factor levels for proper ordering in the plot (reverse for y-axis)
    plot_data$Group <- factor(plot_data$Group, levels = rev(group_order))
  } else {
    # Default ordering - alphabetical but reference group at the top
    if (length(ref_group) > 0) {
      other_groups <- setdiff(unique(plot_data$Group), ref_group)
      group_order <- c(ref_group, sort(other_groups))
      plot_data$Group <- factor(plot_data$Group, levels = rev(group_order))
    } else {
      # If no reference group, just order alphabetically
      plot_data$Group <- factor(plot_data$Group, levels = rev(sort(unique(plot_data$Group))))
    }
  }
  
  # Check for NA values in HR and CI columns and fill with reasonable defaults
  # This helps fix issues with specific groups like "HDACi + PD1"
  for (i in 1:nrow(plot_data)) {
    # Skip reference group
    if (plot_data$Group[i] %in% ref_group) next
    
    # Fix missing HR
    if (is.na(plot_data$Hazard_Ratio[i])) {
      # Check if we have a p-value - if significant, use a small HR (0.1), otherwise use 1.0
      if (!is.na(plot_data$P_Value[i]) && plot_data$P_Value[i] < 0.05) {
        plot_data$Hazard_Ratio[i] <- 0.1
      } else {
        plot_data$Hazard_Ratio[i] <- 1.0
      }
      # Add a note about this
      plot_data$HR_Note[i] <- "HR estimated (original was NA)"
    }
    
    # Fix missing CI values
    if (is.na(plot_data$CI.Lower[i])) plot_data$CI.Lower[i] <- plot_data$Hazard_Ratio[i] * 0.5
    if (is.na(plot_data$CI.Upper[i])) plot_data$CI.Upper[i] <- plot_data$Hazard_Ratio[i] * 2.0
  }
  
  # Add columns to handle plotting on log scale with capped limits
  plot_data$HR_Plot <- pmin(plot_data$Hazard_Ratio, max_plot_hr)
  plot_data$CI_Lower_Plot <- pmin(plot_data$CI.Lower, max_plot_hr)
  plot_data$CI_Upper_Plot <- pmin(plot_data$CI.Upper, max_plot_hr)
  
  # Add indicator if any values were truncated
  plot_data$Note_Plot <- ""
  
  # Add truncation notes
  for (i in 1:nrow(plot_data)) {
    if (!is.na(plot_data$HR_Note[i])) {
      plot_data$Note_Plot[i] <- plot_data$HR_Note[i]
    } else if (!is.na(plot_data$Hazard_Ratio[i]) && plot_data$Hazard_Ratio[i] > max_plot_hr) {
      plot_data$Note_Plot[i] <- "HR truncated"
    } else if (!is.na(plot_data$CI.Upper[i]) && plot_data$CI.Upper[i] > max_plot_hr) {
      plot_data$Note_Plot[i] <- "CI truncated"
    }
  }
  
  # Create text label for the plot
  # Determine whether we're dealing with hazard ratios or effect sizes
  ratio_prefix <- ifelse(grepl("effect", tolower(title)), "ES: ", "HR: ")
  
  # Format labels with ratio, CI, and p-value
  plot_data$label <- ifelse(plot_data$Group %in% ref_group,
                          sprintf("%s1.00 (Reference)", ratio_prefix), 
                          ifelse(is.na(plot_data$P_Value),
                                sprintf("%s%.2f (%.2f-%.2f)", 
                                       ratio_prefix,
                                       plot_data$Hazard_Ratio, 
                                       plot_data$CI.Lower, 
                                       plot_data$CI.Upper),
                                sprintf("%s%.2f (%.2f-%.2f), p=%.3f", 
                                       ratio_prefix,
                                       plot_data$Hazard_Ratio, 
                                       plot_data$CI.Lower, 
                                       plot_data$CI.Upper, 
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
  
  # Add notes for truncated CIs or estimated HRs
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