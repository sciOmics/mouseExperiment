# Copyright (c) 2025 Insight BioAnalytics. All rights reserved.
# Proprietary and confidential.

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
#' @param auc_method Method for calculating AUC: "trapezoidal" (default) or "last_observation"
#' @param reference_group Reference group for statistical comparisons (default: first alphabetically)
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
#' }
#'
#' @export
tumor_auc_analysis <- function(df,
                              time_column = "Day",
                              volume_column = "Volume",
                              treatment_column = "Treatment",
                              id_column = "ID",
                              auc_method = c("trapezoidal", "last_observation"),
                              reference_group = NULL) {
  
  # Match arguments
  auc_method <- match.arg(auc_method)
  
  # Check for required columns
  required_cols <- c(time_column, volume_column, treatment_column, id_column)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Get unique subjects and treatment groups
  subjects <- unique(df[[id_column]])
  treatments <- unique(df[[treatment_column]])
  
  # Set reference group if not specified
  if (is.null(reference_group)) {
    reference_group <- sort(treatments)[1]
  } else if (!reference_group %in% treatments) {
    stop("Reference group '", reference_group, "' not found in the data.")
  }
  
  # Function to calculate AUC for one subject
  calculate_subject_auc <- function(subject_data, method) {
    # Sort by time
    subject_data <- subject_data[order(subject_data[[time_column]]), ]
    
    if (method == "trapezoidal") {
      # Trapezoidal method
      times <- subject_data[[time_column]]
      volumes <- subject_data[[volume_column]]
      
      if (length(times) < 2) {
        return(NA) # Need at least 2 points for AUC
      }
      
      # Calculate AUC using trapezoidal rule
      auc <- 0
      for (i in 2:length(times)) {
        dt <- times[i] - times[i-1]
        auc <- auc + dt * (volumes[i] + volumes[i-1]) / 2
      }
      return(auc)
      
    } else if (method == "last_observation") {
      # Last observation carried forward
      # Simply take the latest time point and its volume
      latest <- subject_data[which.max(subject_data[[time_column]]), ]
      return(latest[[volume_column]])
    }
  }
  
  # Calculate AUC for each subject
  auc_results <- list()
  for (subject in subjects) {
    subject_data <- df[df[[id_column]] == subject, ]
    
    # Skip subjects with no data
    if (nrow(subject_data) == 0) {
      next
    }
    
    # Get treatment group for this subject
    treatment <- unique(subject_data[[treatment_column]])
    if (length(treatment) > 1) {
      warning("Subject ", subject, " has multiple treatment assignments. Using the first one.")
      treatment <- treatment[1]
    }
    
    # Calculate AUC
    auc <- calculate_subject_auc(subject_data, auc_method)
    
    # Store result
    auc_results[[subject]] <- list(
      subject = subject,
      treatment = treatment,
      auc = auc
    )
  }
  
  # Convert to data frame
  auc_df <- do.call(rbind, lapply(auc_results, function(x) {
    data.frame(
      Subject = x$subject,
      Treatment = x$treatment,
      AUC = x$auc,
      stringsAsFactors = FALSE
    )
  }))
  
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
  
  # Create ANOVA model
  auc_model <- stats::aov(AUC ~ Treatment, data = auc_df)
  
  # Perform pairwise comparisons
  auc_comparisons <- stats::TukeyHSD(auc_model, "Treatment")
  
  # Create plot
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    auc_plot <- ggplot2::ggplot(auc_df, ggplot2::aes(x = Treatment, y = AUC, fill = Treatment)) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.2), alpha = 0.5) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = "Area Under the Curve (AUC) by Treatment Group",
        subtitle = paste("Method:", auc_method),
        x = "Treatment Group",
        y = "Area Under the Curve"
      )
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