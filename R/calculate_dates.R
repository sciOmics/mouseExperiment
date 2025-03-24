#' Calculate the number of days for measurements
#'
#' @param df A data frame containing longitudinal tumor measurements
#' @param start_date Date the tumors were injected
#' @param date_column The name of the column with dates of measurements
#' @param date_format Optional format string for parsing dates (e.g., "%d-%b" for "24-Mar" format, "%m/%d/%Y" for "03/24/2025" format)
#' @param year Optional year to use when date format doesn't include year
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
#' 
#' # With "DD-Mon" format dates
#' dose_data <- calculate_dates(dose_data, start_date = "24-Mar", 
#'                              date_format = "%d-%b", year = 2023)
#'                              
#' # With "MM/DD/YYYY" format dates
#' combo_data <- calculate_dates(combo_data, start_date = "03/24/2025", 
#'                              date_format = "%m/%d/%Y")
#' }
calculate_dates <- function(df, start_date, date_column = "Date", 
                            date_format = NULL, year = NULL, in_place = FALSE) {
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
  
  # Attempt to auto-detect common date formats if none provided
  if (is.null(date_format)) {
    # Try to detect the date format based on the first date in the dataset
    first_date <- as.character(df[1, date_column])
    
    # Also check start_date format to ensure we handle it correctly
    start_date_str <- as.character(start_date)
    
    # Check for MM/DD/YYYY format (e.g., "03/24/2025")
    if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_date)) {
      date_format <- "%m/%d/%Y"
      message("Auto-detected date format: MM/DD/YYYY")
    }
    # Check for DD-Mon format (e.g., "24-Mar" or "17-Mar")
    else if (grepl("^\\d{1,2}-[A-Za-z]{3}$", first_date) || 
             grepl("^\\d{1,2}-[A-Za-z]{3}$", start_date_str)) {
      date_format <- "%d-%b"
      message("Auto-detected date format: DD-Mon (requires year parameter)")
      if (is.null(year)) {
        warning("No year provided for DD-Mon format. Will use current year for calculation.")
        # Use current year if no year provided
        year <- as.numeric(format(Sys.Date(), "%Y"))
      }
    }
    # Check for YYYY-MM-DD format (e.g., "2023-03-24")
    else if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", first_date)) {
      date_format <- "%Y-%m-%d"
      message("Auto-detected date format: YYYY-MM-DD")
    }
  }
  
  # Parse dates based on whether a specific format was provided or detected
  if (!is.null(date_format)) {
    # Custom date format parsing
    
    # Check if year needs to be added
    if (grepl("%Y", date_format) || !is.null(year)) {
      # If format includes year or year is provided separately
      
      if (!is.null(year) && !grepl("%Y", date_format)) {
        # Append year to date strings if needed
        dates <- paste(df[,date_column], year, sep = "-")
        date_format_with_year <- paste(date_format, "%Y", sep = "-")
        dates_parsed <- as.POSIXct(strptime(dates, format = date_format_with_year))
      } else {
        # Format already includes year
        dates_parsed <- as.POSIXct(strptime(df[,date_column], format = date_format))
      }
      
      # Parse start date with same format
      if (!is.null(year) && !grepl("%Y", date_format)) {
        start_date_with_year <- paste(start_date, year, sep = "-")
        start_date_parsed <- as.POSIXct(strptime(start_date_with_year, format = date_format_with_year))
      } else {
        start_date_parsed <- as.POSIXct(strptime(start_date, format = date_format))
      }
    } else {
      # Format without year and no year provided (not recommended)
      warning("No year information in dates. Results may be incorrect.")
      dates_parsed <- as.POSIXct(strptime(df[,date_column], format = date_format))
      start_date_parsed <- as.POSIXct(strptime(start_date, format = date_format))
    }
  } else {
    # Use anytime for automatic parsing (original behavior)
    message("Using automatic date parsing. For more reliable results, specify date_format.")
    start_date_parsed <- anytime::anytime(start_date)
    dates_parsed <- anytime::anytime(df[,date_column])
  }
  
  # Check if date parsing was successful
  if (any(is.na(dates_parsed))) {
    warning("Some dates could not be parsed. Check your date format.")
  }
  
  if (is.na(start_date_parsed)) {
    # Try a fallback approach for start date
    tryCatch({
      message("Trying alternative approach for parsing start date...")
      # Check if start_date might be in DD-Mon format
      if (grepl("^\\d{1,2}-[A-Za-z]{3}$", start_date)) {
        # Use current year if none provided
        if (is.null(year)) {
          year <- as.numeric(format(Sys.Date(), "%Y"))
          message("Using current year (", year, ") for calculations.")
        }
        # Format with explicit year
        start_date_with_year <- paste(start_date, year, sep = "-")
        date_format_with_year <- "%d-%b-%Y"
        start_date_parsed <- as.POSIXct(strptime(start_date_with_year, format = date_format_with_year))
        
        # Also reparse the data dates to match format
        if (grepl("^\\d{1,2}-[A-Za-z]{3}$", df[1, date_column])) {
          dates <- paste(df[,date_column], year, sep = "-")
          dates_parsed <- as.POSIXct(strptime(dates, format = date_format_with_year))
        }
      } else {
        # Fall back to anytime as last resort
        start_date_parsed <- anytime::anytime(start_date)
      }
      
      # Check if successful now
      if (is.na(start_date_parsed)) {
        stop("Start date could not be parsed. Please specify date_format parameter explicitly.")
      }
    }, error = function(e) {
      stop(paste("Start date could not be parsed:", e$message, 
                 "\nTry specifying date_format parameter (e.g., date_format = '%d-%b', year = 2023)."))
    })
  }
  
  # Calculate days difference
  Day <- as.numeric(round(difftime(dates_parsed, start_date_parsed, units = "days")))
  
  # Add Day column to the result dataframe
  # First check if Day column already exists, and if so, remove it
  if ("Day" %in% colnames(df_result)) {
    df_result <- df_result[, !colnames(df_result) %in% "Day", drop = FALSE]
  }
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
