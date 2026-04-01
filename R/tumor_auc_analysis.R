# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Perform Area Under the Curve (AUC) Analysis for Tumor Growth Data
#'
#' @description
#' This function calculates the area under the tumor growth curve for each subject
#' and performs statistical comparison between treatment groups.
#'
#' @param df Data frame containing tumor growth data
#' @param time_column Name of column containing time information (default: "Day")
#' @param volume_column Name of column containing tumor volume information (default: "Volume")
#' @param treatment_column Name of column containing treatment group information (default: "Treatment")
#' @param id_column Name of column containing individual subject identifiers (default: "ID")
#' @param cage_column Name of column containing cage information (default: "Cage"). This is used to create
#'        unique subject identifiers by combining ID, Treatment, and Cage, ensuring mice with the same ID 
#'        in different cages are treated as distinct subjects.
#' @param auc_method Method for calculating AUC: "trapezoidal" (default) or "last_observation"
#' @param extrapolation_points Number of most recent data points to use when calculating extrapolation
#'        for subjects with incomplete data. The function will use the last N points to fit the extrapolation curve.
#'        Default is 3. This value must be at least 2 to perform extrapolation.
#' @param reference_group Reference group for statistical comparisons (default: first alphabetically)
#' @param colors Optional named vector of colors for treatment groups in the plot (default: NULL)
#' @param point_size Size of the points in the plot (default: 2.5)
#' @param jitter_width Width of the jitter for the points in the plot (default: 0.2)
#'
#' @return A list containing:
#' \describe{
#'   \item{auc_data}{Data frame with AUC values for each subject}
#'   \item{auc_summary}{Summary statistics of AUC by treatment group}
#'   \item{auc_model}{Statistical model comparing AUC between groups}
#'   \item{auc_comparisons}{Pairwise comparisons between treatment groups}
#'   \item{auc_plot}{Plot of AUC by treatment group}
#' }
#'
#' @details
#' The function creates a composite subject identifier using ID, Treatment, and Cage (if available)
#' to ensure proper identification of mice, even when they share the same ID across different cages.
#' This is particularly important for experiments where mice are housed in multiple cages per treatment group.
#'
#' @importFrom stats aov TukeyHSD t.test
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point theme_classic labs
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(combo_treatment_synthetic_data)
#' tumor_data <- calculate_volume(combo_treatment_synthetic_data)
#' tumor_data <- calculate_dates(tumor_data, start_date = "03/24/2025")
#'
#' # Run AUC analysis
#' auc_results <- tumor_auc_analysis(tumor_data)
#'
#' # View results
#' auc_results$auc_summary
#' auc_results$auc_comparisons
#' print(auc_results$auc_plot)
#'
#' # With extrapolation settings and custom colors
#' treatment_colors <- c("Vehicle" = "#1f77b4", "Drug A" = "#ff7f0e", 
#'                       "Drug B" = "#2ca02c", "Combination" = "#d62728")
#' auc_results <- tumor_auc_analysis(
#'   tumor_data, 
#'   extrapolation_points = 4, 
#'   colors = treatment_colors,
#'   point_size = 3
#' )
#' }
#'
#' @export
tumor_auc_analysis <- function(df,
                              time_column = "Day",
                              volume_column = "Volume",
                              treatment_column = "Treatment",
                              id_column = "ID",
                              cage_column = "Cage",
                              auc_method = c("trapezoidal", "last_observation"),
                              extrapolation_points = 3,
                              reference_group = NULL,
                              colors = NULL,
                              point_size = 2.5,
                              jitter_width = 0.2) {
  
  # Match arguments
  auc_method <- match.arg(auc_method)
  
  # Validate extrapolation_points
  if (!is.numeric(extrapolation_points) || extrapolation_points < 2) {
    warning("extrapolation_points must be a number >= 2. Using default value of 3.")
    extrapolation_points <- 3
  }
  
  # Check for required columns
  required_cols <- c(time_column, volume_column, treatment_column, id_column)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for cage column
  use_cage_info <- FALSE
  if (cage_column %in% colnames(df)) {
    use_cage_info <- TRUE
  } else {
    warning("Cage column '", cage_column, "' not found. Proceeding without cage information for unique subject identification.")
  }
  
  # Create composite subject identifiers that include cage information if available
  if (use_cage_info) {
    composite_ids <- paste(df[[id_column]], df[[treatment_column]], df[[cage_column]], sep = "_")
    df$composite_id <- composite_ids
    # Get unique composite IDs instead of just subject IDs
    unique_ids <- unique(composite_ids)
  } else {
    # Fallback to just using ID if cage information is not available
    unique_ids <- unique(df[[id_column]])
  }
  
  # Set reference group if not specified
  treatments <- unique(df[[treatment_column]])
  if (is.null(reference_group)) {
    reference_group <- sort(treatments)[1]
  } else if (!reference_group %in% treatments) {
    stop("Reference group '", reference_group, "' not found in the data.")
  }
  
  # Calculate max experiment time for extrapolation detection
  max_experiment_time <- max(df[[time_column]])
  
  # Function to calculate AUC for one subject
  calculate_subject_auc <- function(subject_data, method) {
    # Sort by time
    subject_data <- subject_data[order(subject_data[[time_column]]), ]
    
    # Get the number of data points for this subject
    n_points <- nrow(subject_data)
    subject_max_time <- max(subject_data[[time_column]], na.rm = TRUE)
    
    # Determine if extrapolation is needed
    # Extrapolation is needed if the subject's last measurement time is less than
    # the maximum time in the experiment
    needs_extrapolation <- subject_max_time < max_experiment_time
    
    # Only attempt extrapolation if we have enough data points
    can_extrapolate <- n_points >= 2 && needs_extrapolation
    
    # By default, assume no extrapolation
    extrapolated_value <- 0 
    is_extrapolated <- FALSE
    
    if (method == "trapezoidal") {
      # Trapezoidal method
      times <- subject_data[[time_column]]
      volumes <- subject_data[[volume_column]]
      
      if (length(times) < 2) {
        return(list(auc = NA, extrapolated = NA)) # Need at least 2 points for AUC
      }
      
      # Calculate AUC using consolidated trapezoidal utility
      auc <- calculate_auc(times, volumes)
      
      # If extrapolation is needed and possible
      if (can_extrapolate) {
        # Limit the number of points used for extrapolation based on the parameter
        # and the available data points
        n_for_extrapolation <- min(extrapolation_points, n_points)
        
        # Get the last n points for the extrapolation
        extrap_data <- tail(subject_data, n_for_extrapolation)
        
        # Perform simple linear extrapolation using the last n points
        if (n_for_extrapolation >= 2) {
          # Fit a linear model to the last n points
          lm_fit <- stats::lm(
            formula = paste0(volume_column, " ~ ", time_column),
            data = extrap_data
          )
          
          # Get the coefficients for extrapolation
          intercept <- stats::coef(lm_fit)[1]
          slope <- stats::coef(lm_fit)[2]
          
          # Extrapolate from the last observed time to the max experiment time
          last_time <- subject_max_time
          last_volume <- extrap_data[[volume_column]][n_for_extrapolation]
          
          # Calculate the extrapolated AUC (trapezoidal rule for the extrapolated part)
          # For the trapezoidal rule, we need to calculate additional volume at max_experiment_time
          predicted_volume <- intercept + slope * max_experiment_time
          dt_extrapolation <- max_experiment_time - last_time
          extrapolated_value <- dt_extrapolation * (last_volume + predicted_volume) / 2
          
          # Mark that we used extrapolation
          is_extrapolated <- TRUE
        }
      }
      
      # Total AUC is the measured AUC plus any extrapolated component
      return(list(auc = auc + extrapolated_value, extrapolated = is_extrapolated))
      
    } else if (method == "last_observation") {
      # Last observation carried forward method
      latest <- subject_data[which.max(subject_data[[time_column]]), ]
      last_volume <- latest[[volume_column]]
      
      # For last observation method, extrapolation means extending the last volume
      # to the maximum experiment time
      if (can_extrapolate) {
        # Calculate the additional AUC from last observation to max experiment time
        dt_extrapolation <- max_experiment_time - subject_max_time
        extrapolated_value <- dt_extrapolation * last_volume
        is_extrapolated <- TRUE
      }
      
      # For LOCF method, AUC is the last volume (for the observed period) plus any extrapolation
      return(list(auc = last_volume + extrapolated_value, extrapolated = is_extrapolated))
    }
  }
  
  # Calculate AUC for each unique subject identifier
  auc_results <- list()
  
  for (unique_id in unique_ids) {
    # Get data for this unique subject
    if (use_cage_info) {
      subject_data <- df[df$composite_id == unique_id, ]
      # Extract original ID from composite ID for reporting
      id_parts <- strsplit(unique_id, "_")[[1]]
      original_id <- id_parts[1]
      treatment <- id_parts[2]
      cage <- id_parts[3]
    } else {
      subject_data <- df[df[[id_column]] == unique_id, ]
      original_id <- unique_id
      
      # Get treatment group for this subject
      treatment <- unique(subject_data[[treatment_column]])
      if (length(treatment) > 1) {
        warning("Subject ", unique_id, " has multiple treatment assignments. Using the first one.")
        treatment <- treatment[1]
      }
      
      # Set cage to NA if not using cage info
      cage <- NA
    }
    
    # Skip subjects with no data
    if (nrow(subject_data) == 0) {
      next
    }
    
    # Calculate AUC and determine if extrapolation was used
    result <- calculate_subject_auc(subject_data, auc_method)
    
    # Store result with additional information
    auc_results[[unique_id]] <- list(
      subject = original_id,
      treatment = treatment,
      cage = cage,
      auc = result$auc,
      extrapolated = result$extrapolated,
      n_points = nrow(subject_data),
      last_time = max(subject_data[[time_column]], na.rm = TRUE)
    )
  }
  
  # Convert to data frame
  auc_df <- do.call(rbind, lapply(auc_results, function(x) {
    data.frame(
      Subject = x$subject,
      Treatment = x$treatment,
      Cage = x$cage,
      AUC = x$auc,
      Extrapolated = x$extrapolated,
      NumPoints = x$n_points,
      LastTime = x$last_time,
      stringsAsFactors = FALSE
    )
  }))
  
  # Add Group column for compatibility with plot_auc
  auc_df$Group <- auc_df$Treatment
  
  # Remove NAs
  auc_df <- auc_df[!is.na(auc_df$AUC), ]
  
  # Calculate summary statistics
  auc_summary <- stats::aggregate(AUC ~ Treatment, data = auc_df, 
                            FUN = function(x) {
                              c(mean = mean(x, na.rm = TRUE),
                                sd = stats::sd(x, na.rm = TRUE),
                                n = sum(!is.na(x)),
                                sem = stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
                            })
  
  # Convert to cleaner format
  auc_summary <- data.frame(
    Treatment = auc_summary$Treatment,
    Mean_AUC = auc_summary$AUC[, "mean"],
    SD_AUC = auc_summary$AUC[, "sd"],
    N = auc_summary$AUC[, "n"],
    SEM_AUC = auc_summary$AUC[, "sem"],
    stringsAsFactors = FALSE
  )
  
  # Calculate extrapolation statistics
  extrapolation_stats <- stats::aggregate(Extrapolated ~ Treatment, data = auc_df, 
                                    FUN = function(x) {
                                      num_extrapolated <- sum(x, na.rm = TRUE)
                                      total_subjects <- length(x)
                                      pct_extrapolated <- 100 * num_extrapolated / total_subjects
                                      c(n_extrapolated = num_extrapolated,
                                        n_total = total_subjects,
                                        pct_extrapolated = pct_extrapolated)
                                    })
  
  # Convert extrapolation stats to data frame
  extrapolation_summary <- data.frame(
    Treatment = extrapolation_stats$Treatment,
    N_Extrapolated = extrapolation_stats$Extrapolated[, "n_extrapolated"],
    N_Total = extrapolation_stats$Extrapolated[, "n_total"],
    Pct_Extrapolated = extrapolation_stats$Extrapolated[, "pct_extrapolated"],
    stringsAsFactors = FALSE
  )
  
  # Merge extrapolation info into summary
  auc_summary <- merge(auc_summary, extrapolation_summary, by = "Treatment", all = TRUE)
  
  # Create ANOVA model
  auc_model <- stats::aov(AUC ~ Treatment, data = auc_df)
  
  # Perform pairwise comparisons
  auc_comparisons <- stats::TukeyHSD(auc_model, "Treatment")
  
  # Create plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Check if plot_auc function exists
    if (exists("plot_auc", mode = "function")) {
      tryCatch({
        auc_plot <- plot_auc(
          auc_data = auc_df,
          title = paste("Area Under the Curve (AUC) by Treatment Group\nMethod:", auc_method),
          show_mean = TRUE,
          error_bar_type = "SEM",
          extrapolated_column = "Extrapolated",
          colors = colors,
          point_size = point_size,
          jitter_width = jitter_width
        )
      }, error = function(e) {
        message("Error using plot_auc function: ", e$message)
        # Create a basic plot as fallback
        auc_plot <- ggplot2::ggplot(auc_df, ggplot2::aes(x = Treatment, y = AUC, color = Treatment)) +
          ggplot2::geom_boxplot(alpha = 0.7) +
          ggplot2::geom_jitter(ggplot2::aes(shape = Extrapolated), 
                             position = ggplot2::position_jitter(width = jitter_width), 
                             size = point_size,
                             alpha = 0.5) +
          ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1)) +
          # Add custom colors if provided
          (if (!is.null(colors)) ggplot2::scale_color_manual(values = colors) else NULL) +
          ggplot2::theme_classic() +
          ggplot2::labs(
            title = "Area Under the Curve (AUC) by Treatment Group",
            subtitle = paste("Method:", auc_method, 
                            ", Min Points for Extrapolation:", extrapolation_points),
            x = "Treatment Group",
            y = "Area Under the Curve"
          )
      })
    } else {
      # Fallback if plot_auc doesn't exist
      auc_plot <- ggplot2::ggplot(auc_df, ggplot2::aes(x = Treatment, y = AUC, color = Treatment)) +
        ggplot2::geom_boxplot(alpha = 0.7) +
        ggplot2::geom_jitter(ggplot2::aes(shape = Extrapolated), 
                           position = ggplot2::position_jitter(width = jitter_width), 
                           size = point_size,
                           alpha = 0.5) +
        ggplot2::scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1)) +
        # Add custom colors if provided
        (if (!is.null(colors)) ggplot2::scale_color_manual(values = colors) else NULL) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          title = "Area Under the Curve (AUC) by Treatment Group",
          subtitle = paste("Method:", auc_method, 
                         ", Min Points for Extrapolation:", extrapolation_points),
          x = "Treatment Group",
          y = "Area Under the Curve"
        )
    }
  } else {
    auc_plot <- NULL
    warning("Package 'ggplot2' is required for plotting but is not available.")
  }
  
  # Return results
  return(list(
    auc_data = auc_df,
    auc_summary = auc_summary,
    auc_model = auc_model,
    auc_comparisons = auc_comparisons,
    auc_plot = auc_plot
  ))
} 