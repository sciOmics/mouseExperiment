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
#' # Plot with dose information
#' plot_survival(df_with_dose, dose_column = "Dose")

plot_survival = function(df, time_column = "Day", censor_column = "Survival_Censor", 
                      treatment_column = "Treatment", cage_column = "Cage", id_column = "ID",
                      dose_column = NULL) {

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

  # Plot the survival curve using ggsurvplot
  survminer::ggsurvplot(surv_fit,
                      data = analysis_df,
                      pval = TRUE,
                      conf.int = FALSE,
                      risk.table = TRUE,
                      ggtheme = ggplot2::theme_classic(),
                      legend.title = ifelse(!is.null(dose_column), "Treatment and Dose", "Treatment"),
                      xlab = "Time",
                      ylab = "Survival Probability")

}
