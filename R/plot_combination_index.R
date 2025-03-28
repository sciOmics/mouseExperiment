#' Plot Combination Index Over Time
#'
#' Creates a line plot visualizing the Combination Index (CI) over time,
#' with clear indication of synergy (CI < 0.85), additivity (0.85 ≤ CI ≤ 1.15), 
#' and antagonism (CI > 1.15) regions.
#'
#' @param synergy_summary A data frame from the synergy_stats function containing 
#'        Time_Point and Combination_Index values over time.
#' @param line_color Color of the Combination Index line. Default is "purple".
#' @param title Optional custom title for the plot. If NULL, a default title is used.
#' @param x_label Label for x-axis. Default is "Days".
#' @param y_label Label for y-axis. Default is "Combination Index".
#' @param y_limits Optional numeric vector of length 2 specifying y-axis limits.
#'        Default is c(0, 2) to match the common CI scale.
#'
#' @return A ggplot2 object representing the Combination Index plot.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_rect annotate scale_y_continuous
#' @importFrom ggplot2 labs theme_classic theme element_text element_blank
#'
#' @examples
#' \dontrun{
#' # First run the synergy analysis
#' synergy_results <- synergy_stats(
#'   df = tumor_data,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Combination"
#' )
#' 
#' # Create the Combination Index plot
#' plot_combination_index(synergy_results$synergy_summary)
#' 
#' # With custom title and y-axis limits
#' plot_combination_index(
#'   synergy_results$synergy_summary,
#'   title = "Drug Combination Synergy Analysis",
#'   y_limits = c(0, 1.5)
#' )
#' }
#'
#' @export
plot_combination_index <- function(synergy_summary, 
                                 line_color = "purple",
                                 title = "Combination Index Over Time",
                                 x_label = "Days",
                                 y_label = "Combination Index",
                                 y_limits = c(0, 2)) {
  
  # Input validation
  if (!is.data.frame(synergy_summary)) {
    stop("synergy_summary must be a data frame")
  }
  
  # Check required columns
  required_columns <- c("Time_Point", "Combination_Index")
  missing_columns <- setdiff(required_columns, colnames(synergy_summary))
  
  if (length(missing_columns) > 0) {
    stop("synergy_summary is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # Define region boundaries
  synergy_threshold <- 0.85    # Below this is synergistic
  antagonism_threshold <- 1.15 # Above this is antagonistic
  
  # Create the base plot with background regions
  p <- ggplot2::ggplot() +
    # Add region for antagonism (CI > 1.15) - Red background
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = Inf, 
                ymin = antagonism_threshold, ymax = Inf),
      fill = "pink", alpha = 0.2
    ) +
    # Add region for additivity (0.85 <= CI <= 1.15) - Grey background
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = Inf, 
                ymin = synergy_threshold, ymax = antagonism_threshold),
      fill = "grey", alpha = 0.4
    ) +
    # Add region for synergy (CI < 0.85) - Blue background
    ggplot2::geom_rect(
      ggplot2::aes(xmin = -Inf, xmax = Inf, 
                ymin = -Inf, ymax = synergy_threshold),
      fill = "lightblue", alpha = 0.2
    ) +
    
    # Add the CI line
    ggplot2::geom_line(
      data = synergy_summary,
      ggplot2::aes(x = Time_Point, y = Combination_Index),
      color = line_color, size = 1.5
    ) +
    
    # Add labels for regions in the background
    ggplot2::annotate(
      "text", 
      x = mean(range(synergy_summary$Time_Point)), 
      y = (y_limits[2] + antagonism_threshold) / 2,
      label = "Antagonistic",
      color = "red",
      size = 6,
      alpha = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = mean(range(synergy_summary$Time_Point)), 
      y = (antagonism_threshold + synergy_threshold) / 2,
      label = "Additive",
      color = "black",
      size = 6
    ) +
    ggplot2::annotate(
      "text",
      x = mean(range(synergy_summary$Time_Point)), 
      y = synergy_threshold / 2,
      label = "Synergistic",
      color = "blue",
      size = 6,
      alpha = 0.8
    ) +
    
    # Set y-axis limits and formatting
    ggplot2::scale_y_continuous(limits = y_limits) +
    
    # Add labels and styling
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  return(p)
} 