#' Title
#'
#' @param df A data frame containing tumor measurements
#' @param volume_column The name of the column storing tumor volumne measurements
#' @param day_column The name of the column with number of days since the beginning of the experiment for each observation
#' @param group_column The name of the column with the group indicator
#' @param ID_column The name of the column with the individual mouse indentifier. For example, ear-tage indicator.
#' @param group_summary_line Boolean. Should a line for the group average be plotted?
#'
#' @return A ggplot of tumor growth curves color by group
#' @export
#'
#' @examples plot_tumor_growth(synthetic_data)

plot_tumor_growth = function(df, volume_column = "Volume", day_column = "Day", group_column = "Group", ID_column = "ID", group_summary_line = TRUE) {

  #create an interaction term as a string
  interaction_columns = paste(as.character(group_column), as.character(ID_column), sep = ",")
  inter = paste0('interaction(', paste0(interaction_columns),')')

  #plot
  ggplot2::ggplot(df, ggplot2::aes_string(x = day_column, y = volume_column, group = inter)) +
    ggplot2::geom_line(ggplot2::aes_string(color = group_column), alpha = 0.5, size = 0.5) +
    ggplot2::geom_point(ggplot2::aes_string(color = group_column), alpha = 0.5, shape = "square") +
    ggplot2::stat_summary(aes_string(group=group_column, color = group_column), fun = mean, geom = "line", size = 1.4) +
    ggplot2::ylab(bquote("Tumor Volume"(mm^3))) +
    ggplot2::scale_x_continuous(expand=c(0,0)) +
    ggplot2::scale_y_continuous(expand=c(0,0)) +
    ggplot2::theme_classic()
}
