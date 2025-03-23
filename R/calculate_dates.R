#' Calculate the number of days for measurements
#'
#' @param df A data frame containing longitudinal tumor measurements
#' @param start_date Date the tumors were injected
#' @param date_column The name of the column with dates of measurements
#' @param in_place Logical, whether to modify the input data frame (TRUE) or return a new data frame (FALSE, default)
#'
#' @return Returns a data frame with a Day column appended. If in_place=FALSE (default), 
#'         returns a new data frame; if in_place=TRUE, modifies the input data frame and returns it.
#' @export
#'
#' @examples
#' \dontrun{
#' # Return a new data frame with days
#' df_with_days <- calculate_dates(synthetic_data, start_date = "2022-02-24")
#' 
#' # Modify the original data frame in place
#' synthetic_data <- calculate_dates(synthetic_data, start_date = "2022-02-24", in_place = TRUE)
#' }
calculate_dates <- function(df, start_date, date_column = "Date", in_place = FALSE) {
  # Input validation
  if(!date_column %in% colnames(df)) {
    stop(paste("Column", date_column, "not found in data frame"))
  }
  
  # Choose whether to modify in place or create a copy
  if (in_place) {
    df_result <- df  # Reference to the original dataframe
  } else {
    df_result <- df  # Create a copy of the data frame
  }
  
  # Convert dates
  start_date_parsed <- anytime::anytime(start_date)
  dates_parsed <- anytime::anytime(df[,date_column])
  
  # Calculate days difference
  Day <- as.numeric(round(difftime(dates_parsed, start_date_parsed, units = "days")))
  
  # Add Day column to the result dataframe
  result <- cbind(df_result, Day)
  
  # If in_place is TRUE and we're in a function, update the original df in the parent environment
  if (in_place) {
    parent_frame <- parent.frame()
    df_name <- deparse(substitute(df))
    if (exists(df_name, envir = parent_frame)) {
      assign(df_name, result, envir = parent_frame)
    }
  }
  
  # Return the result
  return(result)
}
