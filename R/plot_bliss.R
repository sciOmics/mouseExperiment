#' Plot Bliss Synergy Analysis Over Time
#'
#' Creates a line plot showing tumor growth inhibition percentages over time
#' for individual drugs and their combination, along with the calculated Bliss value.
#'
#' @param synergy_summary A data frame from the synergy_stats function containing TGI values
#'        over time for treatment groups and Bliss expected values.
#' @param colors Optional named vector of colors for the lines. If NULL, default colors are used.
#' @param title Optional custom title for the plot. If NULL, a default title is used.
#' @param add_legend Logical indicating whether to display the legend. Default is TRUE.
#' @param x_label Label for x-axis. Default is "Days".
#' @param y_label Label for y-axis. Default is "Percent Tumor Growth Inhibition".
#' @param y_limits Optional numeric vector of length 2 specifying y-axis limits.
#'
#' @return A ggplot2 object representing the Bliss synergy plot.
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual labs theme_minimal theme
#' @importFrom ggplot2 element_text element_line margin
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
#' # Create the Bliss plot
#' plot_bliss(synergy_results$synergy_summary)
#' 
#' # With custom colors
#' custom_colors <- c("Drug A" = "blue", "Drug B" = "red", 
#'                  "Combination" = "purple", "Bliss Value" = "grey")
#' plot_bliss(synergy_results$synergy_summary, colors = custom_colors)
#' }
#'
#' @export
plot_bliss <- function(synergy_summary, 
                      colors = NULL,
                      title = "Tumor Growth Inhibition and Bliss Synergy",
                      add_legend = TRUE,
                      x_label = "Days",
                      y_label = "Percent Tumor Growth Inhibition",
                      y_limits = NULL,
                      drug_a_label = "Drug A",
                      drug_b_label = "Drug B",
                      combo_label = "Combination",
                      expected_label = "Bliss Value") {
  
  # Input validation
  if (!is.data.frame(synergy_summary)) {
    stop("synergy_summary must be a data frame")
  }
  
  # Check required columns: we need at least Time_Point, TGI values, and Bliss expected
  required_columns <- c("Time_Point", "TGI_Drug_A", "TGI_Drug_B", "TGI_Combo", "Bliss_Expected_TGI")
  missing_columns <- setdiff(required_columns, colnames(synergy_summary))
  
  if (length(missing_columns) > 0) {
    stop("synergy_summary is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # Create a tidy format data frame for plotting
  group_levels <- c(drug_a_label, drug_b_label, combo_label, expected_label)
  plot_data <- data.frame(
    Time = rep(synergy_summary$Time_Point, 4),
    TGI = c(synergy_summary$TGI_Drug_A,
            synergy_summary$TGI_Drug_B,
            synergy_summary$TGI_Combo,
            synergy_summary$Bliss_Expected_TGI),
    Group = factor(c(rep(drug_a_label, nrow(synergy_summary)),
                     rep(drug_b_label, nrow(synergy_summary)),
                     rep(combo_label, nrow(synergy_summary)),
                     rep(expected_label, nrow(synergy_summary))),
                   levels = group_levels)
  )
  
  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- stats::setNames(
      c("blue", "red", "purple", "grey"),
      group_levels
    )
  }
  
  # Create line style values using group labels
  lw_vals <- stats::setNames(c(1, 1, 1, 2), group_levels)
  alpha_vals <- stats::setNames(c(1, 1, 1, 0.5), group_levels)
  
  # Create the line plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = TGI, color = Group, group = Group)) +
    # Add lines with appropriate styling
    ggplot2::geom_line(ggplot2::aes(linewidth = Group, alpha = Group)) +
    ggplot2::scale_linewidth_manual(values = lw_vals) +
    ggplot2::scale_alpha_manual(values = alpha_vals) +
    ggplot2::scale_color_manual(values = colors) +
    
    # Add labels and styling
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = if(add_legend) "right" else "none",
      legend.title = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(1, "cm"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.5),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
  
  # Apply y-axis limits if provided
  if (!is.null(y_limits) && length(y_limits) == 2) {
    p <- p + ggplot2::ylim(y_limits)
  }
  
  return(p)
} 