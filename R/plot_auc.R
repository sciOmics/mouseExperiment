#' Plot AUC (Area Under the Curve) Data
#'
#' This function creates a visualization for Area Under the Curve (AUC) data from tumor growth experiments.
#' It produces a single plot showing AUC values by treatment group as individual points with optional
#' mean lines and error bars.
#'
#' @param auc_data A data frame containing AUC values for individual subjects.
#'        Must include columns: 'Group' (treatment group) and 'AUC' (AUC values).
#'        May also include 'ID', 'Cage', or 'Mouse_ID' columns for more detailed visualization.
#' @param colors An optional vector of colors for treatment groups. If NULL, default ggplot colors are used.
#' @param title An optional title for the plot. If NULL, a default title is used.
#' @param show_mean Logical indicating whether to display a mean line for each group. Default is TRUE.
#' @param error_bar_type Character string specifying the type of error bars to display. Options are:
#'        "none" (no error bars),
#'        "SEM" (standard error of the mean),
#'        "SD" (standard deviation),
#'        "CI" (95% confidence interval). Default is "SEM".
#' @param extrapolated_column Character string specifying the column name that indicates whether a data point 
#'        was extrapolated. Values should be TRUE/FALSE or 1/0. If NULL, all points are treated as non-extrapolated. 
#'        Default is NULL.
#' @param group_order Optional vector specifying the order of groups to display. 
#'        If NULL, groups are ordered alphabetically. Default is NULL.
#' @param point_size Numeric value specifying the size of the individual data points. Default is 3.
#' @param jitter_width Numeric value specifying the width of the jitter for individual points. Default is 0.2.
#'
#' @return A ggplot object representing the AUC plot.
#'
#' @importFrom ggplot2 ggplot aes geom_jitter stat_summary geom_errorbar theme_minimal
#' @importFrom ggplot2 labs scale_fill_manual scale_color_manual ggtitle
#' @importFrom stats aggregate qt
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' auc_data <- data.frame(
#'   ID = paste0("Mouse", 1:20),
#'   Group = rep(c("Control", "Treatment A", "Treatment B", "Treatment C"), each=5),
#'   AUC = c(runif(5, 10, 15), runif(5, 8, 12), runif(5, 5, 10), runif(5, 3, 8))
#' )
#' 
#' # Basic usage
#' plot_auc(auc_data)
#' 
#' # With custom error bars and no mean line
#' plot_auc(auc_data, show_mean = FALSE, error_bar_type = "SD")
#' 
#' # With extrapolated data points marked differently
#' auc_data$Extrapolated <- rep(c(TRUE, FALSE), 10)
#' plot_auc(auc_data, extrapolated_column = "Extrapolated")
#' 
#' # Create a custom color palette for treatment groups
#' colors <- c("Control" = "gray", "Treatment A" = "blue", 
#'            "Treatment B" = "green", "Treatment C" = "red")
#' plot_auc(auc_data, colors = colors)
#'
#' # Specify a custom order for the groups
#' plot_auc(auc_data, group_order = c("Control", "Treatment B", "Treatment A", "Treatment C"))
#' }
#'
#' @export
plot_auc <- function(auc_data, 
                     colors = NULL,
                     title = "AUC Values by Treatment Group",
                     show_mean = TRUE,
                     error_bar_type = c("SEM", "SD", "CI", "none"),
                     extrapolated_column = NULL,
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
  if (!is.data.frame(auc_data)) {
    stop("auc_data must be a data frame")
  }
  
  # Check required columns
  required_cols <- c("Group", "AUC")
  if (!all(required_cols %in% colnames(auc_data))) {
    stop("auc_data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Check extrapolated column if provided
  if (!is.null(extrapolated_column)) {
    if (!extrapolated_column %in% colnames(auc_data)) {
      warning("Specified extrapolated_column '", extrapolated_column, "' not found in data. All points will be treated as non-extrapolated.")
      extrapolated_column <- NULL
    }
  }
  
  # Order groups if specified, otherwise sort by group name for consistency
  if (is.null(group_order)) {
    group_levels <- unique(auc_data$Group)
  } else {
    # Validate group_order
    if (!all(group_order %in% auc_data$Group)) {
      missing_groups <- setdiff(group_order, auc_data$Group)
      warning("Some groups in group_order not found in data: ", 
              paste(missing_groups, collapse = ", "))
      group_order <- intersect(group_order, auc_data$Group)
    }
    if (length(group_order) < length(unique(auc_data$Group))) {
      missing_groups <- setdiff(unique(auc_data$Group), group_order)
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
      groups <- unique(auc_data$Group)
      if (length(colors) < length(groups)) {
        warning("Not enough colors provided. Recycling colors.")
        colors <- rep(colors, length.out = length(groups))
      }
      names(colors) <- groups
    }
    
    # Check if all groups have a color
    missing_groups <- setdiff(unique(auc_data$Group), names(colors))
    if (length(missing_groups) > 0) {
      warning("Missing colors for groups: ", paste(missing_groups, collapse = ", "), 
              ". Using default colors.")
    }
  }
  
  # Set factor levels for consistent ordering
  auc_data$Group <- factor(auc_data$Group, levels = group_levels)
  
  # Calculate summary statistics for error bars
  summary_stats <- stats::aggregate(AUC ~ Group, data = auc_data, 
                               FUN = function(x) {
                                 n <- length(x)
                                 mean_x <- mean(x)
                                 sd_x <- stats::sd(x)
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
  
  # Create the base plot
  p <- ggplot2::ggplot(auc_data, ggplot2::aes(x = Group, y = AUC, color = Group))
  
  # Add extrapolated/non-extrapolated points if column is provided
  if (!is.null(extrapolated_column)) {
    # Convert to logical if needed
    if (is.numeric(auc_data[[extrapolated_column]])) {
      auc_data[[extrapolated_column]] <- as.logical(auc_data[[extrapolated_column]])
    }
    
    # Split data for different point styles
    auc_data_extrapolated <- auc_data[auc_data[[extrapolated_column]], ]
    auc_data_nonextrapolated <- auc_data[!auc_data[[extrapolated_column]], ]
    
    # Add non-extrapolated points (filled circles)
    if (nrow(auc_data_nonextrapolated) > 0) {
      p <- p + ggplot2::geom_jitter(data = auc_data_nonextrapolated,
                                   width = jitter_width, height = 0, 
                                   size = point_size, alpha = 0.7,
                                   shape = 16) # filled circle
    }
    
    # Add extrapolated points (empty circles)
    if (nrow(auc_data_extrapolated) > 0) {
      p <- p + ggplot2::geom_jitter(data = auc_data_extrapolated,
                                   width = jitter_width, height = 0, 
                                   size = point_size, alpha = 0.7,
                                   shape = 1) # empty circle
    }
  } else {
    # All points treated as non-extrapolated (filled circles)
    p <- p + ggplot2::geom_jitter(width = jitter_width, height = 0, 
                                 size = point_size, alpha = 0.7,
                                 shape = 16) # filled circle
  }
  
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
        ggplot2::aes(x = Group, y = AUC.Mean, 
                   ymin = AUC.Mean - AUC.SEM, 
                   ymax = AUC.Mean + AUC.SEM),
        width = 0.3, color = "black", inherit.aes = FALSE
      )
    } else if (error_bar_type == "SD") {
      p <- p + ggplot2::geom_errorbar(
        data = summary_stats,
        ggplot2::aes(x = Group, y = AUC.Mean, 
                   ymin = AUC.Mean - AUC.SD, 
                   ymax = AUC.Mean + AUC.SD),
        width = 0.3, color = "black", inherit.aes = FALSE
      )
    } else if (error_bar_type == "CI") {
      p <- p + ggplot2::geom_errorbar(
        data = summary_stats,
        ggplot2::aes(x = Group, y = AUC.Mean, 
                   ymin = AUC.CI_lower, 
                   ymax = AUC.CI_upper),
        width = 0.3, color = "black", inherit.aes = FALSE
      )
    }
  }
  
  # Add labels and theme
  p <- p + ggplot2::labs(
    title = title,
    x = "Treatment Group",
    y = "AUC (Tumor Burden)"
  ) + ggplot2::theme_classic() +
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