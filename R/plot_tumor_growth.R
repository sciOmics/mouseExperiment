#' Plot Tumor Growth Curves by Treatment Group
#'
#' @param df A data frame containing tumor measurements
#' @param volume_column The name of the column storing tumor volume measurements
#' @param day_column The name of the column with number of days since the beginning of the experiment for each observation
#' @param treatment_column The name of the column with the treatment indicator
#' @param cage_column The name of the column with the cage identifier
#' @param ID_column The name of the column with the individual mouse identifier
#' @param dose_column Optional. The name of the column with dose information (if present in the data)
#' @param survival_column The name of the column indicating survival status (1 for death, 0 for alive)
#' @param extrapolate_volumes Boolean. Should missing volumes for mice that have died be extrapolated? Default is FALSE.
#' @param group_summary_line Boolean. Should a line for the treatment group average be plotted? Default is TRUE.
#'
#' @return A ggplot of tumor growth curves colored by treatment group
#' @export
#'
#' @examples
#' \dontrun{
#' data(synthetic_data)
#' df <- calculate_volume(synthetic_data)
#' df <- calculate_dates(df, start_date = "2022-02-24")
#' 
#' # Standard plot using Treatment column
#' plot_tumor_growth(df, treatment_column = "Treatment")
#' 
#' # Plot with extrapolated volumes for deceased mice
#' plot_tumor_growth(df, treatment_column = "Treatment", extrapolate_volumes = TRUE)
#' 
#' # For data with dose information
#' dose_data <- calculate_volume(dose_data)
#' dose_data <- calculate_dates(dose_data, start_date = "24-Mar", date_format = "%d-%b", year = 2023)
#' plot_tumor_growth(dose_data, treatment_column = "Treatment", dose_column = "Dose")
#' }
plot_tumor_growth <- function(df, volume_column = "Volume", day_column = "Day", 
                             treatment_column = "Treatment", cage_column = "Cage", ID_column = "ID", 
                             dose_column = NULL,
                             survival_column = "Survival_Censor", 
                             extrapolate_volumes = FALSE,
                             group_summary_line = TRUE) {
  
  # Input validation
  req_cols <- c(volume_column, day_column, treatment_column, cage_column, ID_column)
  if (!all(req_cols %in% colnames(df))) {
    stop("Missing required columns in data frame: ", 
         paste(req_cols[!req_cols %in% colnames(df)], collapse = ", "))
  }
  
  # Check for survival column if extrapolation is requested
  if (extrapolate_volumes && !(survival_column %in% colnames(df))) {
    stop(paste("Column", survival_column, "not found in data frame, but is needed for extrapolation"))
  }
  
  # Check for dose column if specified
  if (!is.null(dose_column) && !(dose_column %in% colnames(df))) {
    warning(paste("Dose column", dose_column, "not found in data frame, proceeding without dose information"))
    dose_column <- NULL
  }
  
  # Create a copy of the dataframe for plotting
  plot_df <- df
  
  # Create a composite group identifier based on Treatment (and Dose if available)
  if (!is.null(dose_column)) {
    # Create a group identifier combining Treatment and Dose
    plot_df$Group <- paste(plot_df[[treatment_column]], plot_df[[dose_column]], sep = " - Dose: ")
  } else {
    # Use Treatment as the group identifier
    plot_df$Group <- plot_df[[treatment_column]]
  }
  
  # Create a composite mouse identifier based on Cage, Treatment, ID (and Dose if available)
  if (!is.null(dose_column)) {
    plot_df$Mouse_ID <- paste(plot_df[[cage_column]], plot_df[[treatment_column]], 
                             plot_df[[dose_column]], plot_df[[ID_column]], sep = "_")
  } else {
    plot_df$Mouse_ID <- paste(plot_df[[cage_column]], plot_df[[treatment_column]], 
                             plot_df[[ID_column]], sep = "_")
  }
  
  # Ensure the group column is treated as a factor
  plot_df$Group <- factor(plot_df$Group)
  
  # Add an Extrapolated flag column that will be used later
  plot_df$Extrapolated <- FALSE
  
  # Extrapolate missing volumes for mice that have died if requested
  if (extrapolate_volumes) {
    # Get all unique time points (days)
    all_days <- sort(unique(plot_df[[day_column]]))
    max_day <- max(all_days)
    
    # For each unique mouse (using the Mouse_ID composite identifier), check if they died
    # and extrapolate values if needed
    plot_with_extrap <- plot_df  # Start with the original data
    
    # Get unique mouse identifiers
    unique_mice <- unique(plot_df$Mouse_ID)
    
    # Loop through each mouse
    for (mouse_id in unique_mice) {
      # Get this mouse's data
      mouse_data <- plot_df[plot_df$Mouse_ID == mouse_id, ]
      
      # Check if this mouse died (has Survival_Censor = 1)
      if (any(mouse_data[[survival_column]] == 1)) {
        # Find when the mouse died
        death_rows <- which(mouse_data[[survival_column]] == 1)
        death_day <- min(mouse_data[[day_column]][death_rows])
        death_row <- mouse_data[mouse_data[[day_column]] == death_day, ]
        last_volume <- death_row[[volume_column]][1]
        
        # Calculate growth rate based on previous measurements
        growth_rate <- 0.1  # Default fallback growth rate
        
        # Sort mouse data by day
        mouse_data <- mouse_data[order(mouse_data[[day_column]]), ]
        
        # If there are at least 3 measurements, calculate growth rate
        if (nrow(mouse_data) >= 3) {
          death_index <- which(mouse_data[[day_column]] == death_day)
          if (death_index >= 3) {
            # Use the last 3 measurements 
            days <- mouse_data[[day_column]][(death_index-2):death_index]
            volumes <- mouse_data[[volume_column]][(death_index-2):death_index]
            
            # Only calculate if all volumes are positive
            if (all(volumes > 0)) {
              # Calculate exponential growth rate
              day_diffs <- diff(days)
              volume_ratios <- diff(log(volumes))
              growth_rates <- volume_ratios / day_diffs
              growth_rate <- mean(growth_rates, na.rm = TRUE)
            }
          }
        }
        
        # Find days after death for extrapolation
        future_days <- all_days[all_days > death_day]
        
        if (length(future_days) > 0) {
          # Create template row for extrapolation (using the death row as template)
          for (future_day in future_days) {
            # Copy the row from death day as template
            new_row <- death_row
            
            # Calculate days since death
            days_after_death <- future_day - death_day
            
            # Calculate extrapolated volume using exponential growth model
            extrapolated_volume <- last_volume * exp(growth_rate * days_after_death)
            
            # Update the new row with extrapolated data
            new_row[[day_column]] <- future_day
            new_row[[volume_column]] <- extrapolated_volume
            new_row$Extrapolated <- TRUE
            
            # Add to the plot dataframe
            plot_with_extrap <- rbind(plot_with_extrap, new_row)
          }
        }
      }
    }
    
    # Replace the plot dataframe with the one containing extrapolations
    plot_df <- plot_with_extrap
  }
  
  # Base plot with individual growth curves using the composite identifiers
  plot <- ggplot2::ggplot(plot_df, ggplot2::aes_string(x = day_column, y = volume_column, group = "Mouse_ID")) +
    ggplot2::geom_line(ggplot2::aes_string(color = "Group"), alpha = 0.5, size = 0.5) +
    ggplot2::geom_point(ggplot2::aes_string(color = "Group"), alpha = 0.5, shape = "square")
  
  # If extrapolation was done, mark extrapolated points differently
  if (extrapolate_volumes) {
    extrapolated_points <- plot_df[plot_df$Extrapolated == TRUE, ]
    if (nrow(extrapolated_points) > 0) {
      plot <- plot + 
        ggplot2::geom_point(
          data = extrapolated_points,
          ggplot2::aes_string(x = day_column, y = volume_column, color = "Group"),
          shape = 1,  # hollow circle
          alpha = 0.7,
          size = 1
        ) +
        ggplot2::geom_line(
          data = extrapolated_points,
          ggplot2::aes_string(x = day_column, y = volume_column, color = "Group", group = "Mouse_ID"),
          linetype = "dashed",
          alpha = 0.5,
          size = 0.5
        )
    }
  }
  
  # Add group summary line if requested
  if (group_summary_line) {
    plot <- plot + ggplot2::stat_summary(
      ggplot2::aes_string(group = "Group", color = "Group"),
      fun = mean, geom = "line", size = 1.4
    )
  }
  
  # Create title based on whether we're using dose information
  if (!is.null(dose_column)) {
    title <- paste("Tumor Growth by Treatment and Dose")
  } else {
    title <- "Tumor Growth by Treatment Group"
  }
  
  # Add styling and labels
  plot <- plot +
    ggplot2::ylab(bquote("Tumor Volume"(mm^3))) +
    ggplot2::xlab("Day") +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_classic()
  
  # Add a caption if extrapolation was used
  if (extrapolate_volumes) {
    plot <- plot + ggplot2::labs(
      caption = "Note: Dashed lines represent extrapolated tumor volumes after mouse death"
    )
  }
  
  return(plot)
}
