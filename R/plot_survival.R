#' Plot Kaplan-Meier Survival Curves
#'
#' This function generates Kaplan-Meier survival curves for different treatment groups in the dataset.
#'
#' @param df A data frame containing survival data.
#' @param time_column A string specifying the name of the column representing time to event. Default is "Day".
#' @param censor_column A string specifying the name of the column indicating censoring (1 = event occurred, 0 = censored). Default is "Survival_Censor".
#' @param treatment_column A string specifying the name of the column representing treatment groups. Default is "Treatment".
#' @param cage_column A string specifying the name of the column with the cage identifier. Default is "Cage".
#' @param id_column A string specifying the name of the column with the individual mouse identifier. Default is "ID".
#' @param dose_column Optional. A string specifying the name of the column with dose information. Default is NULL.
#' @param colors Optional. A named vector of colors for each group, or a vector of colors to be assigned to groups in the order they appear. If NULL, default ggplot2 colors are used.
#' @param palette Optional. A color palette function or vector that can be used with scale_color_manual. If provided, overrides the colors parameter.
#' @param palette_indices Optional. If providing a palette, these indices specify which palette colors to use for which groups. If NULL, uses colors in order of groups.
#' @param ggtheme Optional. A ggplot2 theme to use for the plot. Default is theme_classic().
#' @param show_risk_table Logical. Whether to display the risk table below the plot. Default is TRUE.
#' @param show_pvalue Logical. Whether to display the p-value on the plot. Default is FALSE.
#' @param show_legend Logical. Whether to display the legend on the plot. Default is FALSE.
#' @param legend_title Character. Title for the legend. Default is NULL (no title).
#' @param risk_table_title Character. Title for the risk table. Default is NULL (no title).
#'
#' @return A Kaplan-Meier survival plot.
#' @import survival survminer
#' @export
#'
#' @examples
#' df <- data.frame(
#'   Day = c(10, 20, 30, 40, 50, 60, 70, 80),
#'   Survival_Censor = c(1, 0, 1, 1, 0, 1, 0, 1),
#'   Treatment = c("Control", "Control", "Drug", "Drug", "Control", "Drug", "Control", "Drug"),
#'   Cage = c(1, 1, 2, 2, 1, 2, 1, 2),
#'   ID = c("A", "B", "C", "D", "E", "F", "G", "H")
#' )
#'
#' # Basic usage
#' plot_survival(df)
#' 
#' # With custom colors
#' plot_survival(df, colors = c("Control" = "blue", "Drug" = "red"))
#' 
#' # With palette and specific indices
#' palette <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF")
#' plot_survival(df, palette = palette, palette_indices = c(2, 4))
#' 
#' # Show legend
#' plot_survival(df, show_legend = TRUE, legend_title = "Treatment Group")
#' 
#' # With dose information
#' df_with_dose <- data.frame(
#'   Day = c(10, 20, 30, 40, 50, 60, 70, 80),
#'   Survival_Censor = c(1, 0, 1, 1, 0, 1, 0, 1),
#'   Treatment = c("Control", "Control", "Drug", "Drug", "Control", "Drug", "Control", "Drug"),
#'   Cage = c(1, 1, 2, 2, 1, 2, 1, 2),
#'   ID = c("A", "B", "C", "D", "E", "F", "G", "H"),
#'   Dose = c(0, 0, 10, 50, 0, 100, 0, 10)
#' )
#' 
#' # Plot with dose information and custom colors
#' plot_survival(df_with_dose, dose_column = "Dose", 
#'               colors = c("Control - Dose: 0" = "blue", 
#'                          "Drug - Dose: 10" = "red", 
#'                          "Drug - Dose: 50" = "green", 
#'                          "Drug - Dose: 100" = "purple"))

plot_survival = function(df, time_column = "Day", censor_column = "Survival_Censor", 
                      treatment_column = "Treatment", cage_column = "Cage", id_column = "ID",
                      dose_column = NULL, colors = NULL, palette = NULL, palette_indices = NULL,
                      ggtheme = ggplot2::theme_classic(), show_risk_table = TRUE, 
                      show_pvalue = FALSE, show_legend = FALSE, legend_title = NULL,
                      risk_table_title = NULL, fix_risk_table_labels = TRUE) {

  # Input validation
  required_columns <- c(time_column, censor_column, treatment_column, cage_column, id_column)
  missing_cols <- required_columns[!required_columns %in% base::colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for dose column if specified
  if (!is.null(dose_column) && !(dose_column %in% base::colnames(df))) {
    warning(paste("Dose column", dose_column, "not found in data frame, proceeding without dose information"))
    dose_column <- NULL
  }
  
  # Create a subject identifier to ensure each mouse is only counted once
  df$subject_id <- paste(df[[cage_column]], df[[id_column]], sep = "_")
  
  # Create a composite group identifier based on Treatment (and Dose if available)
  if (!is.null(dose_column)) {
    # Create a group identifier combining Treatment and Dose
    df$group <- paste(df[[treatment_column]], df[[dose_column]], sep = " - Dose: ")
  } else {
    # Use Treatment as the group identifier
    df$group <- df[[treatment_column]]
  }
  
  # Aggregate data to subject level - keep only the last time point per subject
  # First sort by subject_id and time to ensure we get the latest record per subject
  df <- df[order(df$subject_id, df[[time_column]]), ]
  
  # Get the last entry for each subject
  subjects <- unique(df$subject_id)
  last_records <- list()
  
  for (subject in subjects) {
    subject_data <- df[df$subject_id == subject, ]
    last_records[[length(last_records) + 1]] <- subject_data[nrow(subject_data), ]
  }
  
  # Combine into a single dataframe with one row per subject
  df_aggregated <- do.call(rbind, last_records)
  
  # Create the survival object from the specified columns using the aggregated data
  survival_time <- df_aggregated[[time_column]]
  event_status <- df_aggregated[[censor_column]]
  grouping <- df_aggregated$group
  
  # Create a new data frame with standardized column names
  analysis_df <- data.frame(
    time = survival_time,
    status = event_status,
    group = grouping
  )
  
  # Fit the Kaplan-Meier survival curve with a fixed formula
  surv_fit <- survival::survfit(survival::Surv(time, status) ~ group, data = analysis_df)

  # Set default legend title if not provided
  if (is.null(legend_title)) {
    legend_title <- ifelse(!is.null(dose_column), "Treatment and Dose", "Treatment")
  }
  
  # Build the list of plot arguments
  plot_args <- list(
    surv_fit,
    data = analysis_df,
    pval = show_pvalue,
    conf.int = FALSE,
    risk.table = show_risk_table,
    ggtheme = ggtheme,
    legend = if(show_legend) "top" else "none",
    legend.title = legend_title,
    risk.table.title = risk_table_title,
    xlab = "Time",
    ylab = "Survival Probability",
    tables.height = 0.3,  # Set a reasonable height for the tables
    risk.table.col = "strata", # Use strata colors
    # The following options help with risk table label display
    risk.table.y.text = TRUE,  # Show the y-axis text
    # Customize the tables theme to properly format y-axis labels
    tables.theme = ggplot2::theme(
      # The axis.text.y formatter function doesn't work directly, but we'll fix it later
      axis.text.y = ggplot2::element_text(hjust = 1, size = 10)
    )
  )
  
  # Handle coloring based on the provided parameters
  unique_groups <- unique(analysis_df$group)
  
  if (!is.null(palette) && !is.null(palette_indices)) {
    # If palette and indices are provided, use them for specific colors
    if (length(palette_indices) < length(unique_groups)) {
      warning("Not enough palette indices provided. Using available indices and reverting to default colors for remaining groups.")
      # Pad with NULL to let ggsurvplot use default colors for remaining groups
      color_values <- palette[palette_indices[1:min(length(palette_indices), length(unique_groups))]]
    } else {
      # Use the specified indices from the palette
      color_values <- palette[palette_indices[1:length(unique_groups)]]
    }
    # Add color parameter to plot args
    plot_args$palette <- color_values
    
    # Print the colors being used for transparency
    message("Using the following palette colors:")
    for (i in 1:length(unique_groups)) {
      if (i <= length(color_values)) {
        message(paste("  -", unique_groups[i], ":", color_values[i]))
      } else {
        message(paste("  -", unique_groups[i], ": default color"))
      }
    }
    
  } else if (!is.null(palette)) {
    # If only palette is provided with no indices, use colors in order
    if (length(palette) < length(unique_groups)) {
      warning("Palette has fewer colors than groups. Using available colors and reverting to default colors for remaining groups.")
    }
    plot_args$palette <- palette
    
  } else if (!is.null(colors)) {
    # If custom colors are provided
    if (is.null(names(colors)) && length(colors) >= length(unique_groups)) {
      # If unnamed vector with enough colors, assign them in order
      names(colors) <- unique_groups[1:length(colors)]
    } else if (is.null(names(colors))) {
      # If unnamed vector with not enough colors, use default
      warning("Not enough unnamed colors provided. Reverting to default colors.")
      colors <- NULL
    } else {
      # If named vector, check if all groups are covered
      missing_groups <- setdiff(unique_groups, names(colors))
      if (length(missing_groups) > 0) {
        warning("Some groups don't have assigned colors: ", 
                paste(missing_groups, collapse = ", "), ". Default colors will be used for these groups.")
      }
    }
    
    if (!is.null(colors)) {
      plot_args$palette <- colors
    }
  }
  
  # Add custom formatting function for risk table
  if (show_risk_table && fix_risk_table_labels) {
    # Create a clean labels formatter that removes "group=" prefix
    format_labels <- function(x) gsub("^group=", "", x)
    
    # Add this formatter to the survminer arguments
    plot_args$risk.table.y.text.col = FALSE  # Don't color the text by group
    plot_args$ylab.tables = list(risk.table = "")  # No title on y-axis in risk table
    
    # Create a function to process the ggsurvplot output to fix group labels
    clean_labels <- function(p) {
      if (!is.null(p$table) && inherits(p$table, "ggplot")) {
        p$table <- p$table + 
          ggplot2::scale_y_discrete(labels = format_labels)
      }
      return(p)
    }
  } else {
    # If not fixing labels, just return the plot as is
    clean_labels <- function(p) p
  }
  
  # Create the plot with the constructed arguments
  survplot <- do.call(survminer::ggsurvplot, plot_args)
  
  # Apply the label cleaning function
  survplot <- clean_labels(survplot)
  
  # Return the survplot object
  return(survplot)
}