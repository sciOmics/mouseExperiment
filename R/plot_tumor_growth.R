#' Plot Tumor Growth Curves by Treatment Group
#' 
#' @importFrom utils tail
#'
#' @param df A data frame containing tumor measurements
#' @param volume_column The name of the column storing tumor volume measurements
#' @param day_column The name of the column with number of days since the beginning of the experiment for each observation
#' @param treatment_column The name of the column with the treatment indicator
#' @param cage_column The name of the column with the cage identifier
#' @param ID_column The name of the column with the individual mouse identifier
#' @param dose_column Optional. The name of the column with dose information (if present in the data)
#' @param survival_column The name of the column indicating survival status (1 for death, 0 for alive)
#' @param extrapolate_volumes Boolean. Should missing volumes for mice that have died be extrapolated? Extrapolated values are shown with dashed lines. Default is FALSE.
#' @param extrapolation_points Character or numeric. Number of data points to use for extrapolation: "all" uses all available data points, a numeric value uses that many recent points, with fallback options if needed. Default is "all".
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
#' 
#' # Combining multiple options: dose information and extrapolation
#' plot_tumor_growth(dose_data, 
#'                 treatment_column = "Treatment", 
#'                 dose_column = "Dose",
#'                 extrapolate_volumes = TRUE)
#' }
plot_tumor_growth <- function(df, volume_column = "Volume", day_column = "Day", 
                             treatment_column = "Treatment", cage_column = "Cage", ID_column = "ID", 
                             dose_column = NULL,
                             survival_column = "Survival_Censor", 
                             extrapolate_volumes = FALSE,
                             extrapolation_points = "all",
                             group_summary_line = TRUE,
                             point_size = 2) {
  
  # Input validation — cage_column is optional
  req_cols <- c(volume_column, day_column, treatment_column, ID_column)
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

  # When cage_column is NULL or absent, inject a synthetic placeholder so that
  # downstream df[[cage_column]] calls do not crash with get1index errors.
  if (is.null(cage_column) || !cage_column %in% colnames(plot_df)) {
    cage_column <- ".cage_placeholder"
    plot_df[[cage_column]] <- "1"
  }
  
  # Capture any user-specified factor ordering on the treatment column before
  # constructing the composite Group column so it can be preserved below.
  tx_levels <- if (is.factor(plot_df[[treatment_column]])) {
    levels(plot_df[[treatment_column]])
  } else {
    NULL
  }

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
  
  # Convert Group to a factor, preserving the user-specified treatment ordering.
  if (!is.null(dose_column)) {
    # For dose data the Group is a composite string:
    # Build ordered levels by iterating over treatment levels in order.
    if (!is.null(tx_levels)) {
      ordered_grp_levels <- unique(unlist(lapply(tx_levels, function(tx) {
        sort(unique(plot_df$Group[as.character(plot_df[[treatment_column]]) == tx]))
      })))
      plot_df$Group <- factor(plot_df$Group, levels = ordered_grp_levels)
    } else {
      plot_df$Group <- factor(plot_df$Group)
    }
  } else {
    # Non-dose path: apply treatment levels directly (Group == treatment column values).
    if (!is.null(tx_levels)) {
      plot_df$Group <- factor(plot_df$Group, levels = tx_levels)
    } else {
      plot_df$Group <- factor(plot_df$Group)
    }
  }
  
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
    
    # Handle extrapolation_points parameter
    using_all_points <- identical(extrapolation_points, "all")
    if (!using_all_points && !is.numeric(extrapolation_points)) {
      warning("extrapolation_points must be 'all' or a positive number. Defaulting to 'all'.")
      using_all_points <- TRUE
      extrapolation_method <- "all available data points"
    } else if (!using_all_points) {
      extrapolation_points <- as.integer(extrapolation_points)
      if (extrapolation_points < 1) {
        warning("extrapolation_points must be at least 1. Defaulting to 3.")
        extrapolation_points <- 3
      }
      extrapolation_method <- paste(extrapolation_points, "most recent data points")
    } else {
      extrapolation_method <- "all available data points"
    }
    
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
        
        # Calculate growth rate using log-linear regression
        growth_rate <- 0.1  # Default fallback growth rate
        
        # Sort mouse data by day
        mouse_data <- mouse_data[order(mouse_data[[day_column]]), ]
        
        # Extract time and volume data
        days <- mouse_data[[day_column]]
        volumes <- mouse_data[[volume_column]]
        
        # Use all available time points with positive volumes for growth rate calculation
        positive_indices <- which(volumes > 0)
        
        # Check if we have enough data for the requested extrapolation method
        if (length(positive_indices) >= 2) {
          # Determine which data points to use based on extrapolation_points parameter
          if (using_all_points) {
            # Use all data points with positive volumes
            selected_indices <- positive_indices
          } else {
            # Use the most recent n points where n = extrapolation_points
            if (length(positive_indices) >= extrapolation_points) {
              selected_indices <- tail(positive_indices, extrapolation_points)
            } else {
              # Not enough points, use what we have
              selected_indices <- positive_indices
            }
          }
          
          growth_times <- days[selected_indices]
          growth_volumes <- volumes[selected_indices]
          
          # If we have at least 2 points, proceed with calculation
          if (length(growth_times) >= 2) {
            if (using_all_points && length(growth_times) > 2) {
              # Use linear regression on log-transformed data for all points
              log_volumes <- log(growth_volumes)
              
              # Fit log-linear model
              growth_model <- tryCatch({
                lm(log_volumes ~ growth_times)
              }, error = function(e) {
                # If model fitting fails, fall back to simpler calculation
                NULL
              })
              
              if (!is.null(growth_model) && !is.na(coef(growth_model)[2])) {
                # Extract growth rate from model (slope of the log-linear model)
                growth_rate <- as.numeric(coef(growth_model)[2])
                
                # Check if growth rate is reasonable (between 0.01 and 0.5 per day)
                if (growth_rate < 0.01 || growth_rate > 0.5) {
                  # If growth rate is extreme, fall back to recent points method
                  if (length(positive_indices) >= 3) {
                    last_indices <- tail(positive_indices, 3)
                    day_diffs <- diff(days[last_indices])
                    volume_ratios <- diff(log(volumes[last_indices]))
                    growth_rates <- volume_ratios / day_diffs
                    growth_rate <- mean(growth_rates, na.rm = TRUE)
                  }
                }
              } else {
                # Fallback if regression model fails
                fallback_indices <- tail(positive_indices, min(3, length(positive_indices)))
                if (length(fallback_indices) >= 2) {
                  # Use simple growth rate between points
                  day_diffs <- diff(days[fallback_indices])
                  volume_ratios <- diff(log(volumes[fallback_indices]))
                  growth_rates <- volume_ratios / day_diffs
                  growth_rate <- mean(growth_rates, na.rm = TRUE)
                }
              }
            } else {
              # For specific number of points or fallback: use simpler calculation
              day_diffs <- diff(growth_times)
              volume_ratios <- diff(log(growth_volumes))
              growth_rates <- volume_ratios / day_diffs
              growth_rate <- mean(growth_rates, na.rm = TRUE)
            }
            
            # Final reasonableness check
            if (is.na(growth_rate) || growth_rate < 0.01 || growth_rate > 0.5) {
              growth_rate <- 0.1  # Use default if calculated rate is unreasonable
            }
          }
        }
        
        # Find days after death for extrapolation
        future_days <- all_days[all_days > death_day]
        
        if (length(future_days) > 0) {
          # Collect extrapolated rows in a list, then bind once
          extrap_rows <- vector("list", length(future_days))
          for (i in seq_along(future_days)) {
            future_day <- future_days[i]
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
            
            extrap_rows[[i]] <- new_row
          }
          # Bind all extrapolated rows at once
          plot_with_extrap <- rbind(plot_with_extrap, do.call(rbind, extrap_rows))
        }
      }
    }
    
    # Replace the plot dataframe with the one containing extrapolations
    plot_df <- plot_with_extrap
  }
  
  # Base plot with individual growth curves — non-extrapolated points only
  # fill = Group is included so that dual-fill shapes (e.g. pch 25, down-triangle)
  # render correctly when scale_fill_manual is added externally.
  plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[day_column]], y = .data[[volume_column]], group = Mouse_ID)) +
    ggplot2::geom_line(
      data = ~ subset(., !Extrapolated),
      ggplot2::aes(color = Group), alpha = 0.5, linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = ~ subset(., !Extrapolated),
      ggplot2::aes(color = Group, fill = Group, shape = Group), alpha = 0.5, size = point_size
    )
  
  # If extrapolation was done, mark extrapolated points differently
  if (extrapolate_volumes) {
    extrapolated_points <- plot_df[plot_df$Extrapolated == TRUE, ]
    if (nrow(extrapolated_points) > 0) {
      # Include each mouse's last real point so the dashed line connects to the observed data
      mice_with_extrap <- unique(extrapolated_points$Mouse_ID)
      last_real_points <- do.call(rbind, lapply(mice_with_extrap, function(mid) {
        real_pts <- plot_df[plot_df$Mouse_ID == mid & !plot_df$Extrapolated, ]
        if (nrow(real_pts) == 0) return(NULL)
        real_pts[which.max(real_pts[[day_column]]), ]
      }))
      extrap_line_data <- rbind(last_real_points, extrapolated_points)

      plot <- plot + 
        ggplot2::geom_point(
          data = extrapolated_points,
          ggplot2::aes(x = .data[[day_column]], y = .data[[volume_column]], color = Group),
          shape = 1,  # hollow circle
          alpha = 0.7,
          size = point_size
        ) +
        ggplot2::geom_line(
          data = extrap_line_data,
          ggplot2::aes(x = .data[[day_column]], y = .data[[volume_column]], color = Group, group = Mouse_ID),
          linetype = "dashed",
          alpha = 0.5,
          linewidth = 0.5
        )
    }
  }
  
  # Add group summary line if requested
  if (group_summary_line) {
    plot <- plot + ggplot2::stat_summary(
      ggplot2::aes(group = Group, color = Group),
      fun = function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else mean(x) },
      na.rm = TRUE, geom = "line", linewidth = 1.4
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
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)), limits = c(0, NA)) + # Start at 0, 5% right padding
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10)), limits = c(0, NA)) + # Start at 0, 10% top buffer
    ggplot2::theme_classic() +
    # Move legend to top-left with proper spacing from title
    ggplot2::theme(
      legend.position = c(0.1, 0.85),  # Nudged down from 0.9 to avoid overlap with title
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.title = ggplot2::element_blank(),  # Remove the legend title
      legend.margin = ggplot2::margin(t = 10)  # Add top margin to increase space from title
    )
  
  # Create captions based on what features are used
  captions <- c()
  
  # Add caption for extrapolation if used
  if (extrapolate_volumes) {
    captions <- c(captions, "Dashed lines represent extrapolated tumor volumes after mouse death")
  }
  
  # Apply all captions if any exist
  if (length(captions) > 0) {
    caption_text <- paste(captions, collapse = "\n")
    plot <- plot + ggplot2::labs(caption = caption_text)
  }
  
  return(plot)
}
