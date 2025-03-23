#' Plot Tumor Growth Curves by Treatment Group
#'
#' @param df A data frame containing tumor measurements
#' @param volume_column The name of the column storing tumor volume measurements
#' @param day_column The name of the column with number of days since the beginning of the experiment for each observation
#' @param group_column The name of the column with the group indicator
#' @param ID_column The name of the column with the individual mouse identifier. For example, ear-tag indicator.
#' @param group_summary_line Boolean. Should a line for the group average be plotted? Default is TRUE.
#'
#' @return A ggplot of tumor growth curves colored by group
#' @export
#'
#' @examples
#' \dontrun{
#' data(synthetic_data)
#' df <- calculate_volume(synthetic_data)
#' df <- calculate_dates(df, start_date = "2022-02-24")
#' plot_tumor_growth(df)
#' }
plot_tumor_growth <- function(df, volume_column = "Volume", day_column = "Day", group_column = "Group", ID_column = "ID", group_summary_line = TRUE) {
  
  # Input validation
  req_cols <- c(volume_column, day_column, group_column, ID_column)
  if (!all(req_cols %in% colnames(df))) {
    stop("Missing required columns in data frame: ", 
         paste(req_cols[!req_cols %in% colnames(df)], collapse = ", "))
  }
  
  # Create an interaction term as a string
  interaction_columns <- paste(as.character(group_column), as.character(ID_column), sep = ",")
  inter <- paste0('interaction(', paste0(interaction_columns),')')
  
  # Base plot with individual growth curves
  plot <- ggplot2::ggplot(df, ggplot2::aes_string(x = day_column, y = volume_column, group = inter)) +
    ggplot2::geom_line(ggplot2::aes_string(color = group_column), alpha = 0.5, size = 0.5) +
    ggplot2::geom_point(ggplot2::aes_string(color = group_column), alpha = 0.5, shape = "square")
  
  # Add group summary line if requested
  if (group_summary_line) {
    plot <- plot + ggplot2::stat_summary(
      ggplot2::aes_string(group = group_column, color = group_column),
      fun = mean, geom = "line", size = 1.4
    )
  }
  
  # Add styling and labels
  plot <- plot +
    ggplot2::ylab(bquote("Tumor Volume"(mm^3))) +
    ggplot2::xlab("Day") +
    ggplot2::ggtitle("Tumor Growth by Treatment Group") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_classic()
  
  return(plot)
}
