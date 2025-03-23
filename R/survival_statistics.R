#' Perform Survival Analysis and Return Hazard Ratios, Confidence Intervals, and P-Values
#'
#' This function performs survival analysis using the Cox Proportional Hazards Model.
#' It calculates the hazard ratio for different groups, along with the 95% confidence intervals
#' and p-values to assess the statistical significance of the group variable on survival outcomes.
#' The function uses the `survival` package to fit the Cox model and returns a summary of the results.
#'
#' @param df A data frame containing the data to be analyzed.
#' @param time_column A character string specifying the name of the column representing the time to event (default: "Day").
#' @param censor_column A character string specifying the name of the column representing the censoring indicator (1 = event, 0 = censored) (default: "Survival_Censor").
#' @param group_column A character string specifying the name of the column representing the treatment or group variable (default: "Group").
#' @param id_column A character string specifying the name of the column representing the individual identifier for each mouse (default: "ID").
#'
#' @return A list containing:
#' \item{model}{The fitted Cox Proportional Hazards model object.}
#' \item{results}{A data frame containing the following columns:
#' \itemize{
#'   \item{Group}{The treatment or group category.}
#'   \item{Hazard_Ratio}{The hazard ratio (exp(coef)) for the group.}
#'   \item{CI_Lower}{The lower bound of the 95% confidence interval for the hazard ratio.}
#'   \item{CI_Upper}{The upper bound of the 95% confidence interval for the hazard ratio.}
#'   \item{P_Value}{The p-value associated with the group effect in the Cox model.}
#' }}
#'
#' @details The function fits a Cox Proportional Hazards model using the `survival` package,
#' and the output includes hazard ratios, confidence intervals, and p-values for the specified treatment groups.
#' The function assumes that the time-to-event data and censoring indicators are appropriately represented
#' in the input data frame, and the model accounts for individual mice as random effects.
#'
#' @examples
#' # Example usage
#' results <- survival_statistics(df = tumor_data, time_column = "Day", censor_column = "Survival_Censor", group_column = "Group", id_column = "ID")
#' print(results$results)
#'
#' @import survival
#'
#' @export
survival_statistics <- function(df, time_column = "Day", censor_column = "Survival_Censor", group_column = "Group", id_column = "ID") {

  # Ensure required columns exist in the dataframe
  if (!all(c(time_column, censor_column, group_column, id_column) %in% colnames(df))) {
    stop("One or more specified columns are not found in the dataframe.")
  }

  # Create a unique identifier for each mouse (based on ID and Group)
  df$Mouse_ID <- factor(paste(df[[group_column]], df[[id_column]], sep = "_"))

  # Create a survival object using the time and censor columns
  surv_obj = survival::Surv(time = df[[time_column]], event = df[[censor_column]])

  # Fit the Cox Proportional Hazards Model
  cox_model = survival::coxph(surv_obj ~ df[[group_column]] + (1 | df$Mouse_ID), data = df)

  # Summary of the Cox model
  cox_summary = summary(cox_model)
  print(cox_summary)

  # Extracting the hazard ratio (exp(coef)) and its 95% Confidence Interval (CI)
  hr = exp(cox_summary$coefficients[, "coef"])
  hr_ci_lower = exp(cox_summary$coefficients[, "lower .95"])
  hr_ci_upper = exp(cox_summary$coefficients[, "upper .95"])

  # Adjusted p-values
  p_value = cox_summary$coefficients[, "Pr(>|z|)"]

  # Combine the hazard ratio, CI, and p-values into a data frame for output
  results_df = data.frame(
    Group = rownames(cox_summary$coefficients),
    Hazard_Ratio = hr,
    CI_Lower = hr_ci_lower,
    CI_Upper = hr_ci_upper,
    P_Value = p_value
  )

  # Return the results as a list
  return(list(
    model = cox_model,
    results = results_df
  ))
}
