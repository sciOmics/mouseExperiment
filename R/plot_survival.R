#' Plot Kaplan-Meier Survival Curves
#'
#' This function generates Kaplan-Meier survival curves for different groups in the dataset.
#'
#' @param df A data frame containing survival data.
#' @param time_column A string specifying the name of the column representing time to event. Default is "Day".
#' @param censor_column A string specifying the name of the column indicating censoring (1 = event occurred, 0 = censored). Default is "Survival_Censor".
#' @param group_column A string specifying the name of the column representing different groups for comparison. Default is "Group".
#'
#' @return A Kaplan-Meier survival plot.
#' @import survival survminer
#' @export
#'
#' @examples
#' df <- data.frame(
#'   Day = c(10, 20, 30, 40, 50, 60, 70, 80),
#'   Survival_Censor = c(1, 0, 1, 1, 0, 1, 0, 1),
#'   Group = c("A", "A", "B", "B", "A", "B", "A", "B")
#' )
#'
#' plot_survival(df)

plot_survival = function(df, time_column = "Day", censor_column = "Survival_Censor", group_column = "Group"){

  # Check required columns exist
  if (!all(c(time_column, censor_column, group_column) %in% colnames(df))) {
    stop("One or more specified columns are not found in the dataframe.")
  }

  # Create a survival object
  surv_obj = Surv(time = df[[time_column]], event = df[[censor_column]])

  # Fit the Kaplan-Meier survival curve
  surv_fit = survfit(surv_obj ~ df[[group_column]], data = df)

  # Plot the survival curve using ggsurvplot
  ggsurvplot(surv_fit,
             data = df,
             pval = TRUE,
             conf.int = TRUE,
             risk.table = TRUE,
             ggtheme = theme_classic(),
             legend.title = group_column,
             xlab = "Time",
             ylab = "Survival Probability")


}
