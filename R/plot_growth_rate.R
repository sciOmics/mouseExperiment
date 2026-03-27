#' Plot Growth Rates from Tumor Growth Experiments
#'
#' This function creates a visualization for tumor growth rates from the output of the 
#' tumor_growth_statistics function. It produces a single plot showing growth rates by 
#' treatment group as individual points with optional mean lines and error bars.
#'
#' @param growth_data A data frame containing growth rate values for individual subjects.
#'        Must include columns: 'Treatment' (treatment group) and 'growth_rate' (growth rate values).
#'        May also include 'ID' for more detailed visualization.
#' @param colors An optional vector of colors for treatment groups. If NULL, default ggplot colors are used.
#' @param title An optional title for the plot. If NULL, a default title is used.
#' @param show_mean Logical indicating whether to display a mean line for each group. Default is TRUE.
#' @param error_bar_type Character string specifying the type of error bars to display. Options are:
#'        "none" (no error bars),
#'        "SEM" (standard error of the mean),
#'        "SD" (standard deviation),
#'        "CI" (95% confidence interval). Default is "SEM".
#' @param group_order Optional vector specifying the order of groups to display. 
#'        If NULL, groups are ordered alphabetically. Default is NULL.
#' @param point_size Numeric value specifying the size of the individual data points. Default is 3.
#' @param jitter_width Numeric value specifying the width of the jitter for individual points. Default is 0.2.
#'
#' @return A ggplot object representing the growth rates plot.
#'
#' @importFrom ggplot2 ggplot aes geom_jitter stat_summary geom_errorbar theme_minimal
#' @importFrom ggplot2 labs scale_fill_manual scale_color_manual ggtitle
#' @importFrom stats aggregate qt
#'
#' @examples
#' \dontrun{
#' # Calculate tumor growth statistics
#' result <- tumor_growth_statistics(
#'   data, 
#'   time_column = "Day", 
#'   volume_column = "Volume",
#'   id_column = "Mouse_ID", 
#'   treatment_column = "Treatment"
#' )
#' 
#' # Extract growth rates data
#' growth_data <- result$growth_rates
#' 
#' # Basic usage
#' plot_growth_rate(growth_data)
#' 
#' # With custom error bars and no mean line
#' plot_growth_rate(growth_data, show_mean = FALSE, error_bar_type = "SD")
#' 
#' # Create a custom color palette for treatment groups
#' colors <- c("Control" = "gray", "Treatment A" = "blue", 
#'            "Treatment B" = "green", "Treatment C" = "red")
#' plot_growth_rate(growth_data, colors = colors)
#'
#' # Specify a custom order for the groups
#' plot_growth_rate(growth_data, group_order = c("Control", "Treatment B", "Treatment A", "Treatment C"))
#' }
#'
#' @export
plot_growth_rate <- function(growth_data, 
                     colors = NULL,
                     title = "Tumor Growth Rates by Treatment Group",
                     show_mean = TRUE,
                     error_bar_type = c("SEM", "SD", "CI", "none"),
                     group_order = NULL,
                     point_size = 3,
                     jitter_width = 0.2) {
  
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function. Please install it.")
  }
  
  # Match arguments
  error_bar_type <- match.arg(error_bar_type)
  
  # Validate inputs
  if (!is.data.frame(growth_data)) {
    stop("growth_data must be a data frame")
  }
  
  # Check required columns
  required_cols <- c("Treatment", "growth_rate") 
  if (!all(required_cols %in% colnames(growth_data))) {
    # Try auto-detection with alternative column names
    if ("Group" %in% colnames(growth_data) && !("Treatment" %in% colnames(growth_data))) {
      growth_data$Treatment <- growth_data$Group
    }
    
    if (!all(required_cols %in% colnames(growth_data))) {
      stop("growth_data must contain columns: ", paste(required_cols, collapse = ", "))
    }
  }
  
  # Check if we need to create or use composite IDs that include cage information
  id_column_to_use <- "ID"
  if ("ID" %in% colnames(growth_data)) {
    # Case 1: ID column already contains composite IDs with underscores (e.g., "01_Control_1")
    if (any(grepl("_", growth_data$ID))) {
      # Keep using the existing ID column as it already contains composite IDs
      id_column_to_use <- "ID"
    } 
    # Case 2: Need to create composite IDs using ID, Treatment, and Cage
    else if ("Cage" %in% colnames(growth_data)) {
      # Create composite ID
      growth_data$composite_id <- paste(growth_data$ID, growth_data$Treatment, growth_data$Cage, sep="_")
      id_column_to_use <- "composite_id"
    }
  }
  
  # Check if the source data might have duplicate mice (same ID in different cages)
  if (id_column_to_use == "composite_id") {
    # Count unique IDs vs unique composite IDs
    n_ids <- length(unique(growth_data$ID))
    n_composite_ids <- length(unique(growth_data$composite_id))
    
    if (n_composite_ids > n_ids) {
      cat("Note: Detected", n_composite_ids, "unique mice when accounting for cage information vs", 
          n_ids, "when using ID alone.\n")
    }
  }
  
  # Order groups if specified, otherwise sort by group name for consistency
  if (is.null(group_order)) {
    group_levels <- unique(growth_data$Treatment)
  } else {
    # Validate group_order
    if (!all(group_order %in% growth_data$Treatment)) {
      missing_groups <- setdiff(group_order, growth_data$Treatment)
      warning("Some groups in group_order not found in data: ", 
              paste(missing_groups, collapse = ", "))
      group_order <- intersect(group_order, growth_data$Treatment)
    }
    if (length(group_order) < length(unique(growth_data$Treatment))) {
      missing_groups <- setdiff(unique(growth_data$Treatment), group_order)
      warning("Not all groups specified in group_order. Adding missing groups at the end: ", 
              paste(missing_groups, collapse = ", "))
      group_order <- c(group_order, missing_groups)
    }
    group_levels <- group_order
  }
  
  # Set up colors if provided
  if (!is.null(colors)) {
    has_names <- !is.null(names(colors)) && any(names(colors) != "")
    if (is.vector(colors) && !has_names) {
      # If just a vector of colors without names, match to groups
      groups <- unique(growth_data$Treatment)
      if (length(colors) < length(groups)) {
        warning("Not enough colors provided. Recycling colors.")
        colors <- rep(colors, length.out = length(groups))
      }
      names(colors) <- groups
    }
    
    # Check if all groups have a color
    missing_groups <- setdiff(unique(growth_data$Treatment), names(colors))
    if (length(missing_groups) > 0) {
      warning("Missing colors for groups: ", paste(missing_groups, collapse = ", "), 
              ". Using default colors.")
    }
  }
  
  # Set factor levels for consistent ordering
  growth_data$Treatment <- factor(growth_data$Treatment, levels = group_levels)
  
  # Make a copy of growth_data that has a Subject column for clarity
  plot_data <- growth_data
  plot_data$Subject <- if (id_column_to_use == "composite_id") plot_data$composite_id else plot_data$ID
  
  # Calculate summary statistics for error bars
  summary_stats <- stats::aggregate(growth_rate ~ Treatment, data = plot_data, 
                               FUN = function(x) {
                                 n <- length(x)
                                 mean_x <- mean(x, na.rm = TRUE)
                                 sd_x <- stats::sd(x, na.rm = TRUE)
                                 sem_x <- sd_x / sqrt(n)
                                 ci_x <- sem_x * stats::qt(0.975, df = n-1)
                                 
                                 c(
                                   Mean = mean_x,
                                   SD = sd_x,
                                   SEM = sem_x,
                                   CI_lower = mean_x - ci_x,
                                   CI_upper = mean_x + ci_x,
                                   N = n
                                 )
                               })
  summary_stats <- do.call(data.frame, summary_stats)
  
  # Print summary of mice per treatment group
  cat("Mice per treatment group in plot_growth_rate:\n")
  for (treatment in levels(plot_data$Treatment)) {
    n_mice <- length(unique(plot_data$Subject[plot_data$Treatment == treatment]))
    cat("  ", treatment, ": ", n_mice, " mice\n", sep="")
  }
  
  # Create the base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Treatment, y = growth_rate, color = Treatment))
  
  # Add points
  p <- p + ggplot2::geom_jitter(ggplot2::aes(shape = Treatment),
                             width = jitter_width, height = 0,
                             size = point_size, alpha = 0.7)
  
  # Add mean line if requested
  if (show_mean) {
    p <- p + ggplot2::stat_summary(fun = mean, geom = "crossbar", 
                                 width = 0.5, size = 0.5, color = "black")
  }
  
  # Add error bars if requested
  if (error_bar_type != "none") {
    if (error_bar_type == "SEM") {
      p <- p + ggplot2::geom_errorbar(
        data = summary_stats,
        ggplot2::aes(x = Treatment, y = growth_rate.Mean, 
                   ymin = growth_rate.Mean - growth_rate.SEM, 
                   ymax = growth_rate.Mean + growth_rate.SEM),
        width = 0.3, color = "black", inherit.aes = FALSE
      )
    } else if (error_bar_type == "SD") {
      p <- p + ggplot2::geom_errorbar(
        data = summary_stats,
        ggplot2::aes(x = Treatment, y = growth_rate.Mean, 
                   ymin = growth_rate.Mean - growth_rate.SD, 
                   ymax = growth_rate.Mean + growth_rate.SD),
        width = 0.3, color = "black", inherit.aes = FALSE
      )
    } else if (error_bar_type == "CI") {
      p <- p + ggplot2::geom_errorbar(
        data = summary_stats,
        ggplot2::aes(x = Treatment, y = growth_rate.Mean, 
                   ymin = growth_rate.CI_lower, 
                   ymax = growth_rate.CI_upper),
        width = 0.3, color = "black", inherit.aes = FALSE
      )
    }
  }
  
  # Add a horizontal line at y = 0 to indicate the boundary between tumor growth and shrinkage
  p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5)
  
  # Add labels and theme
  p <- p + ggplot2::labs(
    title = title,
    x = "Treatment Group",
    y = "Growth Rate (per day)"
  ) + ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }
  
  return(p)
} 