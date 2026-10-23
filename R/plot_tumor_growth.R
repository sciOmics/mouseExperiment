#' Plot Tumor Growth Curves by Treatment Group
#' 
#' @importFrom utils tail
#'
#' @param df A data frame containing tumor measurements
#' @param volume_column The name of the column storing tumor volume measurements
#' @param time_column The name of the column with number of days since the beginning of the experiment for each observation
#' @param treatment_column The name of the column with the treatment indicator
#' @param cage_column The name of the column with the cage identifier
#' @param id_column The name of the column with the individual mouse identifier
#' @param dose_column Optional. The name of the column with dose information (if present in the data)
#' @param event_column The name of the column indicating survival status (1 for death, 0 for alive)
#' @param extrapolate_volumes Boolean. Should missing volumes for mice that have died be extrapolated? Extrapolated values are shown with dashed lines. Default is FALSE.
#' @param extrapolation_points Character or numeric. Number of data points to use for extrapolation: "all" uses all available data points, a numeric value uses that many recent points, with fallback options if needed. Default is "all".
#' @param group_summary_line Boolean. Should a line for the treatment group average be plotted? Default is TRUE.
#' @param necrotic_column Optional. Name of a column in \code{df} containing
#'   necrotic flag values (Y/N, yes/no, 1/0, TRUE/FALSE, case-insensitive).
#'   When supplied, an × marker is overlaid on each necrotic observation and a
#'   caption note is added. Defaults to \code{NULL} (no necrotic marking).
#' @param point_size Size of the individual-observation points. Default 2.
#' @param volume_units Units for the y-axis label: "mm3" (default), "cm3", or
#'   \code{NULL} for no unit.
#'
#' @return A ggplot of tumor growth curves colored by treatment group
#' @export
#'
#' @examples
#' plot_tumor_growth(master_synthetic_data, necrotic_column = "Necrotic")
#'
#' # One agent's dose series, labelled by dose
#' plot_tumor_growth(
#'   subset(master_synthetic_data,
#'          Treatment %in% c("Vehicle", "Drug_A Low", "Drug_A Mid", "Drug_A High")),
#'   dose_column = "Dose")
plot_tumor_growth <- function(df, volume_column = "Volume", time_column = "Day",
                             treatment_column = "Treatment", cage_column = "Cage", id_column = "ID",
                             dose_column = NULL,
                             event_column = "Survival_Censor",
                             extrapolate_volumes = FALSE,
                             extrapolation_points = "all",
                             group_summary_line = TRUE,
                             necrotic_column = NULL,
                             point_size = 2,
                             volume_units = "mm3") {
  
  # Input validation — cage_column is optional
  req_cols <- c(volume_column, time_column, treatment_column, id_column)
  if (!all(req_cols %in% colnames(df))) {
    stop("Missing required columns in data frame: ", 
         paste(req_cols[!req_cols %in% colnames(df)], collapse = ", "))
  }
  
  # Check for survival column if extrapolation is requested
  if (extrapolate_volumes && !(event_column %in% colnames(df))) {
    stop(paste("Column", event_column, "not found in data frame, but is needed for extrapolation"))
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
                             plot_df[[dose_column]], plot_df[[id_column]], sep = "_")
  } else {
    plot_df$Mouse_ID <- paste(plot_df[[cage_column]], plot_df[[treatment_column]], 
                             plot_df[[id_column]], sep = "_")
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
    all_days <- sort(unique(plot_df[[time_column]]))
    
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
    } else if (!using_all_points) {
      extrapolation_points <- as.integer(extrapolation_points)
      if (extrapolation_points < 1) {
        warning("extrapolation_points must be at least 1. Defaulting to 3.")
        extrapolation_points <- 3
      }
    }
    
    # Loop through each mouse
    for (mouse_id in unique_mice) {
      # Get this mouse's data
      mouse_data <- plot_df[plot_df$Mouse_ID == mouse_id, ]
      
      # Check if this mouse died (has Survival_Censor = 1)
      if (any(mouse_data[[event_column]] == 1)) {
        # Find when the mouse died
        death_rows <- which(mouse_data[[event_column]] == 1)
        death_day <- min(mouse_data[[time_column]][death_rows])
        death_row <- mouse_data[mouse_data[[time_column]] == death_day, ]
        last_volume <- death_row[[volume_column]][1]
        
        # Calculate growth rate using log-linear regression
        growth_rate <- 0.1  # Default fallback growth rate
        
        # Sort mouse data by day
        mouse_data <- mouse_data[order(mouse_data[[time_column]]), ]
        
        # Extract time and volume data
        days <- mouse_data[[time_column]]
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
            new_row[[time_column]] <- future_day
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
  plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[time_column]], y = .data[[volume_column]], group = Mouse_ID)) +
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
        real_pts[which.max(real_pts[[time_column]]), ]
      }))
      extrap_line_data <- rbind(last_real_points, extrapolated_points)

      plot <- plot + 
        ggplot2::geom_point(
          data = extrapolated_points,
          ggplot2::aes(x = .data[[time_column]], y = .data[[volume_column]], color = Group),
          shape = 1,  # hollow circle
          alpha = 0.7,
          size = point_size
        ) +
        ggplot2::geom_line(
          data = extrap_line_data,
          ggplot2::aes(x = .data[[time_column]], y = .data[[volume_column]], color = Group, group = Mouse_ID),
          linetype = "dashed",
          alpha = 0.5,
          linewidth = 0.5
        )
    }
  }
  
  # Captions accumulate across optional features; necrotic and extrapolation both append here.
  captions <- c()

  # Overlay × markers on necrotic observations
  if (!is.null(necrotic_column) && necrotic_column %in% colnames(plot_df)) {
    necrotic_flag <- tolower(as.character(plot_df[[necrotic_column]])) %in%
      c("y", "yes", "1", "true")
    necrotic_pts <- plot_df[necrotic_flag & !plot_df$Extrapolated, ]
    if (nrow(necrotic_pts) > 0) {
      plot <- plot +
        ggplot2::geom_point(
          data = necrotic_pts,
          ggplot2::aes(x = .data[[time_column]], y = .data[[volume_column]]),
          shape  = 4,              # \u00d7 symbol
          color  = "black",
          size   = point_size * 1.8,
          stroke = 1.2,
          inherit.aes = FALSE
        )
      captions <- c(captions, "\u00d7 indicates necrotic observation")
    }
  }

  # Add group summary line if requested.
  # R20.61: the mean of the animals still measured falls as the largest
  # tumours leave the study (control: 2,479 to 1,696 mm3 as n went 8 to 1,
  # while the true mean rose to 4,949). Each arm's line now ends on the last
  # day at least half its animals, and at least 3, are measured -- the
  # evaluable-day rule of the analyses -- and the caption says so.
  if (group_summary_line) {
    meas <- plot_df[!plot_df$Extrapolated & is.finite(plot_df[[volume_column]]), , drop = FALSE]
    enrolled <- tapply(meas$Mouse_ID, meas$Group, function(x) length(unique(x)))
    mean_df <- do.call(rbind, lapply(split(meas, list(meas$Group, meas[[time_column]]), drop = TRUE),
      function(g) data.frame(Group = g$Group[1], Day = g[[time_column]][1],
                             Mean = mean(g[[volume_column]]),
                             N = length(unique(g$Mouse_ID)), stringsAsFactors = FALSE)))
    if (!is.null(mean_df) && nrow(mean_df)) {
      need <- pmax(pmin(3, enrolled[as.character(mean_df$Group)]),
                   ceiling(0.5 * enrolled[as.character(mean_df$Group)]))
      ok <- mean_df$N >= need
      # Only up to the first day an arm falls short: later days are survivors.
      o <- order(mean_df$Group, mean_df$Day)
      mean_df <- mean_df[o, ]
      ok <- ok[o]
      ok <- unlist(lapply(split(ok, mean_df$Group), cumprod), use.names = FALSE) == 1L
      if (any(!ok)) {
        captions <- c(captions, paste(
          "Mean lines end once fewer than half of an arm's animals (or fewer than 3)",
          "are still measured; after that the survivors' mean understates growth."))
      }
      mean_df <- mean_df[ok, , drop = FALSE]
      if (is.factor(plot_df$Group)) mean_df$Group <- factor(mean_df$Group, levels = levels(plot_df$Group))
      names(mean_df)[names(mean_df) == "Day"] <- time_column
      plot <- plot + ggplot2::geom_line(
        data = mean_df,
        ggplot2::aes(x = .data[[time_column]], y = .data[["Mean"]], group = Group, color = Group),
        linewidth = 1.4, inherit.aes = FALSE)
    }
  }

  # Create title based on whether we're using dose information
  if (!is.null(dose_column)) {
    title <- paste("Tumor Growth by Treatment and Dose")
  } else {
    title <- "Tumor Growth by Treatment Group"
  }
  
  # Add styling and labels
  plot <- plot +
    # R20.61: the unit as given (it was always mm3), and an axis that keeps
    # pre-treatment days (limits = c(0, NA) dropped negative days silently).
    ggplot2::ylab(switch(if (is.null(volume_units)) "" else volume_units,
                         mm3 = bquote("Tumor Volume"(mm^3)),
                         cm3 = bquote("Tumor Volume"(cm^3)),
                         "Tumor Volume")) +
    ggplot2::xlab("Day") +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05)),
                                limits = c(min(0, plot_df[[time_column]], na.rm = TRUE), NA)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.10)), limits = c(0, NA)) + # Start at 0, 10% top buffer
    ggplot2::theme_classic() +
    # Move legend to top-left with proper spacing from title
    ggplot2::theme(
      legend.position = c(0.1, 0.85),  # Nudged down from 0.9 to avoid overlap with title
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.title = ggplot2::element_blank(),  # Remove the legend title
      legend.margin = ggplot2::margin(t = 10)  # Add top margin to increase space from title
    )
  
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
