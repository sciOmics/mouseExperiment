#' Create Kaplan-Meier Survival Plot
#'
#' @description Creates a Kaplan-Meier survival plot from longitudinal tumor growth data.
#' The function handles complex data with treatment groups, dose levels, and cage information.
#'
#' @param data Data frame containing survival data with measurements over time
#' @param time_column Name of column containing time information (default: "Day")
#' @param censor_column Name of column containing censoring information (1=event/death, 0=censored) (default: "Survival_Censor")
#' @param treatment_column Name of column containing treatment group information (default: "Treatment")
#' @param id_column Name of column containing individual subject identifiers (default: "ID")
#' @param cage_column Optional name of column containing cage identifiers (default: NULL)
#' @param dose_column Optional name of column containing dose information (default: NULL)
#' @param colors Optional named vector of colors for groups (default: NULL)
#' @param show_risk_table Logical indicating whether to display risk table (default: FALSE)
#' @param show_censoring Logical indicating whether to show censoring marks (default: TRUE)
#' @param title Plot title (default: NULL)
#' @param subtitle Plot subtitle (default: NULL)
#' @param xlab X-axis label (default: "Time (Days)")
#' @param ylab Y-axis label (default: "Survival Probability")
#' @param font_size Base font size for plot text (default: 12)
#' @param font_family Font family for plot text (default: "sans")
#' @param xlim Optional vector of x-axis limits (default: NULL)
#' @param ylim Optional vector of y-axis limits (default: NULL)
#' @param xbreaks Optional vector of x-axis breaks (default: NULL)
#' @param ybreaks Optional vector of y-axis breaks (default: NULL)
#'
#' @return A survminer::ggsurvplot object containing the survival plot and risk table
#'
#' @import survival
#' @import survminer
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' plot_survival(survival_data)
#'
#' # Specify custom column names and colors
#' plot_survival(survival_data, 
#'              time_column = "Time", 
#'              censor_column = "Death",
#'              treatment_column = "Group",
#'              colors = c("Control" = "blue", "Treatment" = "red"))
#'
#' # Include dose information and risk table with custom axis limits
#' plot_survival(dose_data,
#'              dose_column = "Dose",
#'              show_risk_table = TRUE,
#'              title = "Survival by Treatment and Dose",
#'              subtitle = "Kaplan-Meier Analysis",
#'              xlim = c(0, 100),
#'              ylim = c(0, 1),
#'              xbreaks = seq(0, 100, by = 20),
#'              ybreaks = seq(0, 1, by = 0.2),
#'              show_censoring = TRUE)
#' }
plot_survival <- function(data,
                        time_column = "Day",
                        censor_column = "Survival_Censor",
                        treatment_column = "Treatment",
                        id_column = "ID",
                        cage_column = NULL,
                        dose_column = NULL,
                        colors = NULL,
                        show_risk_table = FALSE,
                        show_censoring = TRUE,
                        title = NULL,
                        subtitle = NULL,
                        xlab = "Time (Days)",
                        ylab = "Survival Probability",
                        font_size = 12,
                        font_family = "sans",
                        xlim = NULL,
                        ylim = NULL,
                        xbreaks = NULL,
                        ybreaks = NULL) {
  
  # Input validation
  required_cols <- c(time_column, censor_column, treatment_column, id_column)
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check optional columns
  if (!is.null(cage_column) && !(cage_column %in% colnames(data))) {
    warning("Cage column '", cage_column, "' not found, proceeding without cage information")
    cage_column <- NULL
  }
  
  if (!is.null(dose_column) && !(dose_column %in% colnames(data))) {
    warning("Dose column '", dose_column, "' not found, proceeding without dose information")
    dose_column <- NULL
  }
  
  # Create unique subject identifier
  if (!is.null(cage_column)) {
    data$subject <- paste(data[[cage_column]], data[[id_column]], sep = ":")
  } else {
    data$subject <- data[[id_column]]
  }
  
  # Create group variable
  if (!is.null(dose_column)) {
    data$group <- paste(data[[treatment_column]], " - Dose:", data[[dose_column]])
  } else {
    data$group <- data[[treatment_column]]
  }
  data$group <- factor(data$group)
  
  # Get last observation for each subject
  data <- data[order(data$subject, data[[time_column]]), ]
  last_obs <- lapply(unique(data$subject), function(s) {
    subject_data <- data[data$subject == s, ]
    subject_data[nrow(subject_data), ]
  })
  survival_data <- do.call(rbind, last_obs)
  
  # Prepare survival data
  survival_data$time <- survival_data[[time_column]]
  survival_data$status <- survival_data[[censor_column]]
  
  # Fit survival model
  fit <- survival::survfit(Surv(time, status) ~ group, data = survival_data)
  
  # Prepare plot arguments
  plot_args <- list(
    fit,
    data = survival_data,
    pval = FALSE,  # Never show p-values
    conf.int = FALSE,
    risk.table = show_risk_table,
    tables.height = 0.25,
    legend.title = "Group",
    legend = "left",
    xlab = xlab,
    ylab = ylab,
    risk.table.title = "Number at risk",
    risk.table.col = "black",
    risk.table.y.text = TRUE,
    risk.table.height = 0.25,
    fontsize = font_size,
    font.family = font_family,
    censor = show_censoring,  # Control censoring marks
    # Risk table theme
    tables.theme = ggplot2::theme_classic() +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = font_size * 0.7),
        plot.title = ggplot2::element_text(size = font_size * 0.8),
        plot.margin = ggplot2::unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
        axis.title = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = font_size * 0.7)
      ),
    # Main plot theme
    ggtheme = ggplot2::theme_classic() + 
      ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.2),
        legend.justification = c(0, 0),
        legend.background = ggplot2::element_rect(fill = "white", color = NA),
        legend.key.size = ggplot2::unit(1, "lines"),
        plot.margin = ggplot2::unit(c(0.5, 0.5, 0.1, 0.5), "cm"),
        plot.title = ggplot2::element_text(size = font_size * 1.2),
        plot.subtitle = ggplot2::element_text(size = font_size * 0.9)
      )
  )
  
  # Add title and subtitle if provided
  if (!is.null(title)) {
    plot_args$title <- title
  }
  if (!is.null(subtitle)) {
    plot_args$subtitle <- subtitle
  }
  
  # Handle custom colors
  if (!is.null(colors)) {
    # Get group levels
    group_levels <- levels(survival_data$group)
    
    # If unnamed colors are provided
    if (is.null(names(colors))) {
      # Only use as many colors as there are groups
      colors <- colors[1:min(length(colors), length(group_levels))]
      # Assign names to the colors
      names(colors) <- group_levels[1:length(colors)]
    }
    
    # Check if any groups are missing colors
    missing_groups <- setdiff(group_levels, names(colors))
    if (length(missing_groups) > 0) {
      warning("Colors not specified for groups: ", 
              paste(missing_groups, collapse = ", "), 
              ". Default colors will be used for these.")
    }
    
    # Add the palette to the plot arguments
    plot_args$palette <- colors
  }
  
  # Create the plot
  surv_plot <- do.call(survminer::ggsurvplot, plot_args)
  
  # Fix risk table labels if shown
  if (show_risk_table && !is.null(surv_plot$table)) {
    surv_plot$table <- surv_plot$table + 
      ggplot2::scale_y_discrete(labels = function(x) gsub("^group=", "", x)) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank()
      )
  }
  
  # Fix legend labels
  surv_plot$plot <- surv_plot$plot + 
    ggplot2::scale_color_discrete(labels = function(x) gsub("^group=", "", x)) +
    ggplot2::scale_fill_discrete(labels = function(x) gsub("^group=", "", x))
  
  # Apply axis limits and breaks if provided
  if (!is.null(xlim)) {
    surv_plot$plot <- surv_plot$plot + ggplot2::coord_cartesian(xlim = xlim)
  }
  if (!is.null(ylim)) {
    surv_plot$plot <- surv_plot$plot + ggplot2::coord_cartesian(ylim = ylim)
  }
  if (!is.null(xbreaks)) {
    surv_plot$plot <- surv_plot$plot + ggplot2::scale_x_continuous(breaks = xbreaks)
  }
  if (!is.null(ybreaks)) {
    surv_plot$plot <- surv_plot$plot + ggplot2::scale_y_continuous(breaks = ybreaks)
  }
  
  # Ensure risk table colors match plot colors
  if (show_risk_table && !is.null(surv_plot$table)) {
    # Don't try to extract colors from the plot as it may cause errors
    # Instead, use the same palette for both
    if (!is.null(colors) && length(names(colors)) > 0) {
      # If we have a named color palette, use it
      surv_plot$table <- surv_plot$table + 
        ggplot2::scale_color_manual(values = colors) +
        ggplot2::scale_fill_manual(values = colors)
    }
  }
  
  return(surv_plot)
}