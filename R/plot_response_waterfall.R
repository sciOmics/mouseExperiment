#' Best-response waterfall plot
#'
#' One bar per animal: its best response (% change in volume from baseline,
#' see [response_categories()]), sorted from the largest increase to the
#' largest decrease. Dashed lines mark the mRECIST thresholds for best
#' response: +35 % (progressive disease above), -50 % (partial response
#' below) and -95 % (complete response below). Bars above `cap` are drawn at
#' `cap`, and a caption says how many.
#'
#' @param responses The result of [response_categories()].
#' @param fill Colour the bars by `"treatment"` (default) or `"category"`.
#' @param colors Optional named colours for the treatments or categories.
#' @param cap Upper limit of the y axis, in %; `NULL` for none.
#' @param title Plot title.
#' @return A ggplot.
#' @export
#' @examples
#' rc <- response_categories(master_synthetic_data, cage_column = "Cage")
#' plot_response_waterfall(rc)
plot_response_waterfall <- function(responses, fill = c("treatment", "category"),
                                    colors = NULL, cap = 200,
                                    title = "Best response by animal") {
  fill <- match.arg(fill)
  d <- responses$per_animal
  d <- d[is.finite(d$Best_Response), , drop = FALSE]
  if (!nrow(d)) stop("No animal has an evaluable response.", call. = FALSE)
  d <- d[order(-d$Best_Response), , drop = FALSE]
  d$Rank <- seq_len(nrow(d))
  d$Shown <- if (is.null(cap)) d$Best_Response else pmin(d$Best_Response, cap)
  d$Fill <- if (fill == "treatment") d$Treatment else
    factor(d$Category, levels = c("mPD", "mSD", "mPR", "mCR"))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[["Rank"]], y = .data[["Shown"]],
                                       fill = .data[["Fill"]])) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
    ggplot2::geom_hline(yintercept = c(35, -50, -95), linetype = "dashed",
                        colour = "grey55") +
    ggplot2::labs(title = title, x = "Animals, by best response",
                  y = "Best response (% change from baseline)",
                  fill = if (fill == "treatment") "Treatment" else "mRECIST") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank())
  n_over <- if (is.null(cap)) 0L else sum(d$Best_Response > cap)
  if (n_over > 0L) {
    p <- p + ggplot2::labs(caption = sprintf(
      "%d animal%s above %s%% drawn at %s%% (largest %.0f%%).", n_over,
      if (n_over == 1L) "" else "s", format(cap), format(cap), max(d$Best_Response)))
  }
  if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
  p
}
