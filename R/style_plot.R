#' Style a plot as the dashboard does
#'
#' The dashboard's plot options as one function: a theme, colours and shapes
#' for the groups, titles and axis labels, the text size, the legend's place,
#' a log scale and the x-axis range. The R script in the dashboard's Download
#' Results calls it with the options you chose, so its plots look like the
#' dashboard's (R22 proposal F9). Each argument left `NULL` leaves that part of
#' the plot as it is.
#'
#' @param p A ggplot.
#' @param theme `"classic"`, `"minimal"` or `"bw"`; NULL keeps the plot's own.
#' @param colors Named colours for the groups, applied to colour and fill.
#' @param shapes Named point shapes (integers) for the groups.
#' @param title,subtitle,x_label,y_label,legend_title Text; NULL or "" keeps the
#'   plot's own.
#' @param title_align `"center"`, `"left"` or `"right"`; NULL keeps the plot's own.
#' @param text_size Base text size, in points.
#' @param legend_position `"right"`, `"left"`, `"top"`, `"bottom"`, `"none"`, or
#'   inside the plot: `"inplot_tl"`, `"inplot_tr"`, `"inplot_bl"`, `"inplot_br"`.
#' @param y_scale `"linear"` or `"log10"`.
#' @param x_limits Two numbers, the x-axis range shown (the data outside it are
#'   kept, as with [ggplot2::coord_cartesian()]).
#' @return The styled ggplot.
#' @export
#' @examples
#' p <- plot_tumor_growth(master_synthetic_data, cage_column = "Cage")
#' style_plot(p, theme = "classic", title = "Tumour growth", text_size = 14,
#'            legend_position = "bottom", y_scale = "log10")
style_plot <- function(p, theme = NULL, colors = NULL, shapes = NULL, title = NULL,
                       subtitle = NULL, title_align = NULL,
                       text_size = NULL, x_label = NULL, y_label = NULL, legend_title = NULL,
                       legend_position = NULL, y_scale = c("linear", "log10"), x_limits = NULL) {
  y_scale <- match.arg(y_scale)
  given <- function(x) length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  if (!is.null(theme)) {
    p <- p + switch(match.arg(theme, c("classic", "minimal", "bw")),
                    classic = ggplot2::theme_classic(), minimal = ggplot2::theme_minimal(),
                    bw = ggplot2::theme_bw())
  }
  if (length(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors) + ggplot2::scale_fill_manual(values = colors)
  }
  if (length(shapes)) p <- p + ggplot2::scale_shape_manual(values = shapes)
  if (given(title))    p <- p + ggplot2::labs(title = title)
  if (given(subtitle)) p <- p + ggplot2::labs(subtitle = subtitle)
  if (given(x_label))  p <- p + ggplot2::labs(x = x_label)
  if (given(y_label))  p <- p + ggplot2::labs(y = y_label)
  if (given(legend_title)) p <- p + ggplot2::labs(fill = legend_title, color = legend_title,
                                                  shape = legend_title)
  if (!is.null(title_align)) {
    hjust <- c(left = 0, center = 0.5, right = 1)[[match.arg(title_align, c("center", "left", "right"))]]
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = hjust),
                            plot.subtitle = ggplot2::element_text(hjust = hjust))
  }
  if (!is.null(text_size)) p <- p + ggplot2::theme(text = ggplot2::element_text(size = text_size))
  if (!is.null(legend_position)) {
    inside <- c(inplot_tl = "left top", inplot_tr = "right top",
                inplot_bl = "left bottom", inplot_br = "right bottom")
    if (legend_position %in% names(inside)) {
      just <- strsplit(inside[[legend_position]], " ")[[1L]]
      p <- p + ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = c(if (just[1] == "left") 0.02 else 0.98,
                                   if (just[2] == "top") 0.98 else 0.02),
        legend.justification = just,
        legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", 0.8),
                                                  colour = NA))
    } else {
      p <- p + ggplot2::theme(legend.position = legend_position)
    }
  }
  if (y_scale == "log10") {
    p <- p + ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.05)))
  }
  if (length(x_limits) == 2L && all(is.finite(x_limits))) {
    p <- suppressMessages(p + ggplot2::coord_cartesian(xlim = x_limits))
  }
  p
}
