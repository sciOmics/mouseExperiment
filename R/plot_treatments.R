#' Plot Treatment Schedule with Day Markers
#'
#' @param treatment_data A data frame containing treatment schedule information
#' @param tumor_growth_data A data frame containing tumor measurements (used to align x-axis)
#' @param day_column The name of the column with number of days for treatment_data
#' @param treatment_column The name of the column with the treatment identifier
#' @param color_palette Optional vector of colors to use for different treatments
#' @param treatment_order Optional character vector specifying the order in which treatments should appear (top to bottom)
#'
#' @return A ggplot of treatment schedule with filled triangles for treatment days
#' @export
#'
#' @examples
#' # Load included datasets
#' data(combo_treatment_synthetic_data)
#' data(combo_treatment_schedule)
#' 
#' # Process the tumor data
#' tumor_data <- calculate_volume(combo_treatment_synthetic_data)
#' tumor_data <- calculate_dates(tumor_data, start_date = "03/24/2025")
#' 
#' # Process treatment schedule data (if needed)
#' treatment_schedule <- calculate_dates(combo_treatment_schedule, 
#'                                     start_date = "03/24/2025")
#' 
#' # Plot treatment schedule
#' plot_treatments(treatment_schedule, tumor_data)
#' 
#' # With custom colors
#' custom_colors <- c("Control" = "gray", "Drug A" = "blue", 
#'                   "Drug B" = "red", "Combo" = "purple")
#' plot_treatments(treatment_schedule, tumor_data, 
#'                color_palette = custom_colors)
#'
#' # With custom treatment order
#' treatment_order <- c("Combo", "Drug B", "Drug A", "Control")
#' plot_treatments(treatment_schedule, tumor_data, 
#'                treatment_order = treatment_order)
plot_treatments <- function(treatment_data, 
                           tumor_growth_data, 
                           day_column = "Day", 
                           treatment_column = "Treatment",
                           color_palette = NULL,
                           treatment_order = NULL) {
  
  # Input validation
  if (!all(c(day_column, treatment_column) %in% colnames(treatment_data))) {
    stop("Missing required columns in treatment_data: ", 
         paste(c(day_column, treatment_column)[!c(day_column, treatment_column) %in% colnames(treatment_data)], 
               collapse = ", "))
  }
  
  # Verify day_column exists in tumor_growth_data for x-axis alignment
  if (!(day_column %in% colnames(tumor_growth_data))) {
    warning("Day column not found in tumor_growth_data. X-axis may not align with tumor growth plots.")
  }
  
  # Create a data frame for plotting
  plot_df <- treatment_data[, c(day_column, treatment_column)]
  
  # Remove any duplicate treatment-day combinations
  plot_df <- unique(plot_df)
  
  # Apply custom treatment ordering if provided
  if (!is.null(treatment_order)) {
    # Validate treatment order - all values must exist in the data
    missing_treatments <- treatment_order[!treatment_order %in% unique(plot_df[[treatment_column]])]
    if (length(missing_treatments) > 0) {
      warning("The following treatments in treatment_order are not found in the data: ", 
             paste(missing_treatments, collapse = ", "))
    }
    
    # Create a factor with specified levels to control the order (in reverse for y-axis top-down)
    plot_df[[treatment_column]] <- factor(plot_df[[treatment_column]], 
                                        levels = rev(treatment_order))
  } else {
    # If no custom order, just create a regular factor (default ordering - alphabetical)
    plot_df[[treatment_column]] <- factor(plot_df[[treatment_column]])
  }
  
  # Get x-axis range from tumor growth data if available
  if (day_column %in% colnames(tumor_growth_data)) {
    x_range <- range(tumor_growth_data[[day_column]], na.rm = TRUE)
  } else {
    x_range <- range(treatment_data[[day_column]], na.rm = TRUE)
  }
  
  # Create the plot
  plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[day_column]], y = .data[[treatment_column]])) +
    ggplot2::geom_point(
      shape = 25,  # Filled triangle pointing down
      size = 4,
      ggplot2::aes(fill = .data[[treatment_column]])
    ) +
    ggplot2::scale_x_continuous(
      "Day", 
      limits = c(0, x_range[2] * 1.05),  # Start at day 0, add some padding at the end
      breaks = seq(0, ceiling(x_range[2]), by = max(1, round(x_range[2]/10)))  # Sensible breaks
    ) +
    ggplot2::ylab("Treatment") +
    ggplot2::ggtitle("Treatment Schedule") +
    ggplot2::theme_classic() +
    # Moderate spacing between groups on y-axis (prevent bottom group from dipping below x-axis)
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0.08, 0.08))) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90"),
      legend.position = "none",  # Hide legend since y-axis already shows treatment names
      # Standard margin for y-axis text
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 3)),
      # Slightly reduced line height for y-axis text
      axis.text = ggplot2::element_text(lineheight = 0.9),
      # Add horizontal lines to make it easier to follow treatments across days
      panel.grid.major.x = ggplot2::element_line(color = "gray95", linetype = "dotted")
    )
  
  # Apply custom color palette if provided
  if (!is.null(color_palette)) {
    plot <- plot + ggplot2::scale_fill_manual(values = color_palette)
  }
  
  return(plot)
}