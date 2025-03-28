#' Plot Area Under the Curve (AUC) Values by Treatment Group
#'
#' @description
#' Creates a visualization of Area Under the Curve (AUC) values by treatment group,
#' with options for displaying mean values, error bars, and differentiation between
#' extrapolated and non-extrapolated data points.
#'
#' @param auc_data A data frame containing AUC values for each subject.
#'        Must have columns for AUC values and group/treatment information.
#' @param auc_column Name of the column containing AUC values. Default is "AUC".
#' @param group_column Name of the column containing group/treatment information. Default is "Group".
#' @param title Plot title. Default is "Area Under the Curve (AUC) by Treatment Group".
#' @param show_mean Logical indicating whether to show lines for group means. Default is FALSE.
#' @param error_bar_type Type of error bars to display. Options are "none", "SEM" (standard error of mean),
#'        "SD" (standard deviation), or "CI" (95% confidence interval). Default is "none".
#' @param extrapolated_column Name of column indicating whether a data point
#'        was extrapolated. Values should be TRUE/FALSE or 1/0. If NULL, all points are treated as non-extrapolated. 
#'        Default is NULL.
#' @param group_order Optional vector specifying the order of groups to display. 
#'        If NULL, groups are ordered alphabetically. Default is NULL.
#' @param caption Optional caption text for the plot. Default is NULL.
#'
#' @return A ggplot object that can be further customized or displayed.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_boxplot geom_jitter
#' @importFrom ggplot2 stat_summary theme_classic labs position_dodge
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' # Generate some example data
#' auc_data <- data.frame(
#'   Subject = paste0("S", 1:20),
#'   Group = rep(c("Control", "Treatment"), each = 10),
#'   AUC = c(rnorm(10, 100, 15), rnorm(10, 70, 10)),
#'   Extrapolated = c(rep(FALSE, 15), rep(TRUE, 5))
#' )
#'
#' # Basic plot
#' plot_auc(auc_data)
#'
#' # With mean, SEM error bars, and differentiation of extrapolated points
#' plot_auc(
#'   auc_data,
#'   title = "Tumor Growth AUC by Group",
#'   show_mean = TRUE,
#'   error_bar_type = "SEM",
#'   extrapolated_column = "Extrapolated"
#' )
#' }
#'
#' @export
plot_auc <- function(auc_data,
                    auc_column = "AUC",
                    group_column = "Group",
                    title = "Area Under the Curve (AUC) by Treatment Group",
                    show_mean = FALSE,
                    error_bar_type = c("none", "SEM", "SD", "CI"),
                    extrapolated_column = NULL,
                    group_order = NULL,
                    caption = NULL) {
  
  # Check inputs
  if (!is.data.frame(auc_data)) {
    stop("auc_data must be a data frame")
  }
  
  if (!auc_column %in% colnames(auc_data)) {
    stop(sprintf("Column '%s' not found in auc_data", auc_column))
  }
  
  if (!group_column %in% colnames(auc_data)) {
    stop(sprintf("Column '%s' not found in auc_data", group_column))
  }
  
  # Match argument for error_bar_type
  error_bar_type <- match.arg(error_bar_type)
  
  # Set group order if specified
  if (!is.null(group_order)) {
    if (!all(unique(auc_data[[group_column]]) %in% group_order)) {
      warning("Not all groups from data are in group_order. Using alphabetical order.")
      auc_data[[group_column]] <- factor(auc_data[[group_column]])
    } else {
      auc_data[[group_column]] <- factor(auc_data[[group_column]], levels = group_order)
    }
  } else {
    # Default to alphabetical order
    auc_data[[group_column]] <- factor(auc_data[[group_column]])
  }
  
  # Determine if we have extrapolation information
  has_extrapolation <- FALSE
  if (!is.null(extrapolated_column) && extrapolated_column %in% colnames(auc_data)) {
    has_extrapolation <- TRUE
    
    # Convert to logical if it's not already
    if (!is.logical(auc_data[[extrapolated_column]])) {
      auc_data[[extrapolated_column]] <- as.logical(auc_data[[extrapolated_column]])
    }
    
    # Add caption about extrapolated points if one was not provided
    if (is.null(caption)) {
      caption <- "Open circles represent extrapolated values"
    }
  }
  
  # Create the base plot
  p <- ggplot2::ggplot(auc_data, 
                      ggplot2::aes(x = .data[[group_column]], 
                                  y = .data[[auc_column]], 
                                  color = .data[[group_column]])) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA)
  
  # Add points with or without extrapolation indication
  if (has_extrapolation) {
    p <- p + ggplot2::geom_jitter(
      ggplot2::aes(shape = .data[[extrapolated_column]]),
      position = ggplot2::position_jitter(width = 0.2),
      size = 2.5,
      alpha = 0.8
    ) +
      ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1),
                                name = "Extrapolated")
  } else {
    p <- p + ggplot2::geom_jitter(
      position = ggplot2::position_jitter(width = 0.2),
      size = 2.5,
      alpha = 0.8
    )
  }
  
  # Add error bars if requested
  if (error_bar_type != "none") {
    # Function to calculate error bar values
    error_fun <- switch(error_bar_type,
                        "SEM" = function(x) {
                          m <- mean(x, na.rm = TRUE)
                          sem <- stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
                          return(c(y = m, ymin = m - sem, ymax = m + sem))
                        },
                        "SD" = function(x) {
                          m <- mean(x, na.rm = TRUE)
                          s <- stats::sd(x, na.rm = TRUE)
                          return(c(y = m, ymin = m - s, ymax = m + s))
                        },
                        "CI" = function(x) {
                          m <- mean(x, na.rm = TRUE)
                          sem <- stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
                          ci <- 1.96 * sem
                          return(c(y = m, ymin = m - ci, ymax = m + ci))
                        })
    
    # Add error bars
    p <- p + ggplot2::stat_summary(
      fun.data = error_fun,
      geom = "errorbar",
      width = 0.2,
      size = 1,
      color = "black"
    )
    
    # Add mean points/lines if requested
    if (show_mean) {
      p <- p + ggplot2::stat_summary(
        fun = mean,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black"
      )
      
      # Add horizontal line at each mean
      group_means <- stats::aggregate(
        formula(paste(auc_column, "~", group_column)),
        data = auc_data,
        FUN = mean
      )
      
      # Convert to long format for geom_segment
      segments_data <- do.call(rbind, lapply(1:nrow(group_means), function(i) {
        g <- as.character(group_means[[group_column]][i])
        y <- group_means[[auc_column]][i]
        data.frame(
          group = g,
          y = y,
          x = as.integer(factor(g)) - 0.25,
          xend = as.integer(factor(g)) + 0.25
        )
      }))
      
      p <- p + ggplot2::geom_segment(
        data = segments_data,
        ggplot2::aes(
          x = x,
          xend = xend,
          y = y,
          yend = y
        ),
        color = "black",
        size = 1
      )
    }
  }
  
  # Apply classic theme and labels
  p <- p + ggplot2::theme_classic() +
    ggplot2::labs(
      title = title,
      x = "Treatment Group",
      y = "Area Under the Curve",
      caption = caption
    ) +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.caption = ggplot2::element_text(hjust = 0, face = "italic")
    )
  
  return(p)
} 