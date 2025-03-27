#' Plot AUC (Area Under the Curve) Data
#'
#' This function creates visualizations for Area Under the Curve (AUC) data from tumor growth experiments.
#' It produces multiple plots showing AUC values by treatment group, including bar plots, box plots, 
#' and individual subject plots.
#'
#' @param auc_data A data frame containing AUC values for individual subjects.
#'        Must include columns: 'Group' (treatment group) and 'AUC' (AUC values).
#'        May also include 'ID', 'Cage', or 'Mouse_ID' columns for more detailed plots.
#' @param plot_type Character string specifying which type of plot to return. Options are:
#'        "all" (returns a list of all plots),
#'        "summary" (group summary bar plot),
#'        "boxplot" (boxplot with jittered points),
#'        "individual" (individual subject bars),
#'        "combined" (a composite of all plots). Default is "combined".
#' @param colors An optional vector of colors for treatment groups. If NULL, default ggplot colors are used.
#' @param title An optional main title for the combined plot. If NULL, a default title is used.
#' @param comparison_lines Logical indicating whether to add lines showing significant comparisons.
#'        Default is FALSE.
#' @param comparison_data A data frame containing pairwise comparison results, required if 
#'        comparison_lines=TRUE. Should have columns: 'group1', 'group2', 'p.value', 'significance'.
#' @param id_column Character string specifying the column name for subject identifiers. 
#'        If NULL, tries to find 'ID', 'Mouse_ID', or creates a row number. Default is NULL.
#' @param cage_column Character string specifying the column name for cage identifiers.
#'        If NULL, tries to find 'Cage', or omits cage information. Default is NULL.
#' @param group_order Optional vector specifying the order of groups to display in plots.
#'        If NULL, groups are ordered alphabetically. Default is NULL.
#'
#' @return If plot_type is "all", returns a list of ggplot objects. Otherwise, returns a single ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar geom_boxplot geom_jitter theme_minimal theme
#' @importFrom ggplot2 element_text labs scale_fill_manual scale_color_manual ggtitle
#' @importFrom ggpubr ggarrange
#' @importFrom stats aggregate
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
#' # Display just the summary bar plot
#' plot_auc(auc_data, plot_type = "summary")
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
                     plot_type = "combined",
                     colors = NULL,
                     title = NULL,
                     comparison_lines = FALSE,
                     comparison_data = NULL,
                     id_column = NULL,
                     cage_column = NULL,
                     group_order = NULL) {
  
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function. Please install it.")
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is needed for this function. Please install it.")
  }
  
  # Validate inputs
  if (!is.data.frame(auc_data)) {
    stop("auc_data must be a data frame")
  }
  
  # Check required columns
  required_cols <- c("Group", "AUC")
  if (!all(required_cols %in% colnames(auc_data))) {
    stop("auc_data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Determine ID column if not provided
  if (is.null(id_column)) {
    possible_id_cols <- c("ID", "Mouse_ID", "Subject")
    existing_id_cols <- intersect(possible_id_cols, colnames(auc_data))
    
    if (length(existing_id_cols) > 0) {
      id_column <- existing_id_cols[1]
    } else {
      # If no ID column exists, create a row number
      auc_data$ID <- 1:nrow(auc_data)
      id_column <- "ID"
    }
  } else if (!id_column %in% colnames(auc_data)) {
    warning("Specified id_column '", id_column, "' not found in data. Using row numbers instead.")
    auc_data$ID <- 1:nrow(auc_data)
    id_column <- "ID"
  }
  
  # Determine cage column if not provided
  if (is.null(cage_column)) {
    possible_cage_cols <- c("Cage", "Cage_ID", "CageID")
    existing_cage_cols <- intersect(possible_cage_cols, colnames(auc_data))
    
    if (length(existing_cage_cols) > 0) {
      cage_column <- existing_cage_cols[1]
    } else {
      # If cage information is explicitly needed but not available, warn the user
      if (plot_type == "individual") {
        warning("No cage column found. Individual plot will only show IDs without cage information.")
      }
      cage_column <- NULL
    }
  } else if (!cage_column %in% colnames(auc_data)) {
    warning("Specified cage_column '", cage_column, "' not found in data.")
    cage_column <- NULL
  }
  
  # Create display ID combining cage and ID if both are available
  if (!is.null(cage_column)) {
    auc_data$Display_ID <- paste(auc_data[[cage_column]], auc_data[[id_column]], sep = "-")
  } else {
    auc_data$Display_ID <- auc_data[[id_column]]
  }
  
  # Calculate summary statistics
  auc_summary <- stats::aggregate(AUC ~ Group, data = auc_data, 
                               FUN = function(x) c(Mean = mean(x), 
                                                  SD = stats::sd(x), 
                                                  N = length(x),
                                                  SEM = stats::sd(x)/sqrt(length(x))))
  auc_summary <- do.call(data.frame, auc_summary)
  
  # Order groups if specified, otherwise sort by group name for consistency
  if (is.null(group_order)) {
    auc_summary <- auc_summary[order(auc_summary$Group), ]
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
    # Reorder summary based on specified order
    auc_summary$Group <- factor(auc_summary$Group, levels = group_order)
    auc_summary <- auc_summary[order(auc_summary$Group), ]
    group_levels <- group_order
  }
  
  # Set up colors if provided
  if (!is.null(colors)) {
    if (is.vector(colors) && !is.named(colors)) {
      # If just a vector of colors, match to groups
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
  
  # Sort data for individual bar plot
  auc_data <- auc_data[order(auc_data$Group, auc_data$AUC), ]
  
  # Set factor levels for consistent ordering
  auc_data$Group <- factor(auc_data$Group, levels = group_levels)
  auc_data$Display_ID <- factor(auc_data$Display_ID, levels = unique(auc_data$Display_ID[order(auc_data$Group)]))
  
  # Create the plots
  plots <- list()
  
  # Raw data scatterplot with means
  plots$raw_data <- ggplot2::ggplot(auc_data, ggplot2::aes(x = Group, y = AUC, color = Group)) +
    ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
    ggplot2::stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", linetype = "dashed") +
    ggplot2::labs(title = "AUC Values by Treatment Group",
                x = "Treatment Group",
                y = "AUC (Tumor Burden)") +
    ggplot2::theme_minimal()
  
  # Summary bar plot
  plots$summary <- ggplot2::ggplot(auc_summary, ggplot2::aes(x = Group, y = AUC.Mean, fill = Group)) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = AUC.Mean - AUC.SEM, 
                                     ymax = AUC.Mean + AUC.SEM), 
                        width = 0.2) +
    ggplot2::labs(title = "Mean AUC by Treatment Group",
                x = "Treatment Group",
                y = "Mean AUC (± SEM)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Boxplot with jittered points
  plots$boxplot <- ggplot2::ggplot(auc_data, ggplot2::aes(x = Group, y = AUC, color = Group)) +
    ggplot2::geom_boxplot(alpha = 0.3) +
    ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.7) +
    ggplot2::labs(title = "AUC Distribution by Treatment Group",
                x = "Treatment Group",
                y = "AUC (Tumor Burden)") +
    ggplot2::theme_minimal()
  
  # Individual subject bar plot
  plots$individual <- ggplot2::ggplot(auc_data, ggplot2::aes(x = Display_ID, y = AUC, fill = Group)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(title = "Individual Subject AUC Values",
                x = if (!is.null(cage_column)) "Subject ID (Cage-ID)" else "Subject ID",
                y = "AUC (Tumor Burden)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 7),
                 legend.position = "top") +
    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    plots$raw_data <- plots$raw_data + ggplot2::scale_color_manual(values = colors)
    plots$summary <- plots$summary + ggplot2::scale_fill_manual(values = colors)
    plots$boxplot <- plots$boxplot + ggplot2::scale_color_manual(values = colors)
    plots$individual <- plots$individual + ggplot2::scale_fill_manual(values = colors)
  }
  
  # Add comparison lines if requested
  if (comparison_lines && !is.null(comparison_data)) {
    if (!all(c("group1", "group2", "p.value") %in% colnames(comparison_data))) {
      warning("comparison_data must contain columns: group1, group2, p.value")
    } else {
      # Add significance stars or p-values to the summary plot
      # This implementation varies based on the exact format of your comparison data
      # A simple example for adding text labels:
      for (i in 1:nrow(comparison_data)) {
        sig_text <- if ("significance" %in% colnames(comparison_data)) {
          comparison_data$significance[i]
        } else {
          if (comparison_data$p.value[i] < 0.001) "***"
          else if (comparison_data$p.value[i] < 0.01) "**"
          else if (comparison_data$p.value[i] < 0.05) "*"
          else "ns"
        }
        
        # Add text to summary plot 
        # Note: This is a simplified approach and might need adjustment
        plots$summary <- plots$summary + 
          ggplot2::annotate("text", 
                          x = (which(levels(auc_data$Group) == comparison_data$group1[i]) + 
                              which(levels(auc_data$Group) == comparison_data$group2[i])) / 2, 
                          y = max(auc_summary$AUC.Mean) * 1.1,
                          label = sig_text)
      }
    }
  }
  
  # Create combined plot
  if (is.null(title)) {
    title <- "AUC Analysis of Tumor Growth"
  }
  
  # Create a single plot with all subplots arranged in a grid
  plots$combined <- ggpubr::ggarrange(
    plots$raw_data, plots$summary, 
    plots$boxplot, plots$individual,
    ncol = 2, nrow = 2,
    common.legend = TRUE,
    legend = "top",
    labels = c("A", "B", "C", "D")
  )
  
  # Add a main title for the combined plot
  main_title_plot <- ggplot2::ggplot() + 
    ggplot2::annotate(geom = "text", x = 0, y = 0, label = title, size = 6, fontface = "bold") + 
    ggplot2::theme_void()
  
  # Combine the title and plots
  plots$combined <- ggpubr::ggarrange(
    main_title_plot, 
    plots$combined, 
    ncol = 1, 
    heights = c(0.1, 0.9)
  )
  
  # Return the appropriate plot(s)
  if (plot_type == "all") {
    return(plots)
  } else if (plot_type %in% names(plots)) {
    return(plots[[plot_type]])
  } else {
    warning("Invalid plot_type '", plot_type, "'. Returning combined plot.")
    return(plots$combined)
  }
} 