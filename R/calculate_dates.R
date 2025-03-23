#' Calculate the number of days for measurements
#'
#' @param df A data frame containing longitudinal tumor measurements
#' @param start_date Date the tumors were injected
#' @param date_column The name of the column with dates of measurements
#'
#' @return Returns the original data frame with an appended column containing the number of days from the beginning of an experiment for each row/measurement
#' @export
#'
#' @examples calculate_dates(synthetic_data, start_date = "2022-02-24")

calculate_dates = function(df, start_date, date_column = "Date") {
  start_date = anytime::anytime(start_date)
  date_column = anytime::anytime(df[,date_column])
  Day = as.numeric(round(difftime(date_column, start_date, units = "days")))
  cbind(df, Day)
}
