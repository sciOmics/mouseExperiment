#' Helper Functions for mouseExperimentDashboard
#' 
#' Utility functions for data validation, formatting, and processing
#'

#' Pick the best default reference group
#'
#' Returns the first group whose name matches common control-group terms
#' (case-insensitive partial match). If none match, returns the first group.
#'
#' @param groups Character vector of treatment group names
#' @return Single character string
#' @keywords internal
pick_reference_group <- function(groups) {
  if (length(groups) == 0) return("")
  control_patterns <- c("vehicle", "control", "igg", "dmso")
  for (g in groups) {
    if (any(grepl(paste(control_patterns, collapse = "|"), g, ignore.case = TRUE))) {
      return(g)
    }
  }
  groups[1]
}

#' Validate uploaded data file
#'
#' @param file_path Path to uploaded file
#' @param file_ext File extension
#'
#' @return List with valid (logical), message (char), and data (df)
#' @keywords internal
validate_uploaded_file <- function(file_path, file_ext) {
  tryCatch({
    if (tolower(file_ext) == "csv") {
      data <- read.csv(file_path, stringsAsFactors = FALSE)
    } else if (tolower(file_ext) %in% c("rda", "rdata")) {
      load(file_path, envir = .GlobalEnv)
      # Get the loaded object (assumes single object)
      data <- get(ls()[1])
    } else {
      return(list(
        valid = FALSE,
        message = "Unsupported file format. Please use CSV or RDA files.",
        data = NULL
      ))
    }
    
    # Basic validation
    if (!is.data.frame(data)) {
      return(list(
        valid = FALSE,
        message = "Loaded file is not a data frame.",
        data = NULL
      ))
    }
    
    if (nrow(data) == 0) {
      return(list(
        valid = FALSE,
        message = "Data frame is empty.",
        data = NULL
      ))
    }
    
    list(
      valid = TRUE,
      message = paste("Successfully loaded", nrow(data), "rows and", 
                      ncol(data), "columns"),
      data = data
    )
  }, error = function(e) {
    list(
      valid = FALSE,
      message = paste("Error loading file:", conditionMessage(e)),
      data = NULL
    )
  })
}

#' Detect column names matching common patterns
#'
#' @param col_names Character vector of column names
#'
#' @return List with detected column mappings
#' @keywords internal
detect_columns <- function(col_names) {
  col_names_lower <- tolower(col_names)
  
  detected <- list(
    id = detect_column(col_names, col_names_lower, c("id", "mouse_id", "animal_id", "subject")),
    day = detect_column(col_names, col_names_lower, c("day", "time", "timepoint", "day_number")),
    volume = detect_column(col_names, col_names_lower, c("volume", "tumor_volume", "tv")),
    length = detect_column(col_names, col_names_lower, c("length", "tumor_length")),
    width = detect_column(col_names, col_names_lower, c("width", "tumor_width")),
    treatment = detect_column(col_names, col_names_lower, c("treatment", "group", "treatment_group")),
    cage = detect_column(col_names, col_names_lower, c("cage", "cage_id")),
    date = detect_column(col_names, col_names_lower, c("date", "measurement_date")),
    dose = detect_column(col_names, col_names_lower, c("dose", "dose_level", "concentration", "conc")),
    survival_time = detect_column(col_names, col_names_lower, c("survival_time", "time_to_event", "days_to_event")),
    survival_event = detect_column(col_names, col_names_lower, c("survival_event", "event", "survival_censor", "censor"))
  )
  
  detected
}

#' Helper to find column in names
#' @keywords internal
detect_column <- function(col_names, col_names_lower, patterns) {
  for (pattern in patterns) {
    match_idx <- which(col_names_lower == pattern)
    if (length(match_idx) > 0) return(col_names[match_idx[1]])  # Return actual column name
  }
  NULL
}

#' Format numeric values for display
#'
#' @param x Numeric value
#' @param digits Number of decimal places
#'
#' @return Formatted character string
#' @keywords internal
format_number <- function(x, digits = 3) {
  if (is.na(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

#' Safely execute analysis function with error handling
#'
#' @param expr Expression to execute
#'
#' @return Result or error message
#' @keywords internal
safe_analysis <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      list(error = TRUE, message = conditionMessage(e))
    }
  )
}

#' Format a p-value for display
#'
#' Very small p-values are shown in scientific notation rather than being
#' rounded to "0". Uses \code{formatC} with \code{format="g"} so the output
#' automatically switches to scientific notation when needed.
#'
#' @param x Numeric p-value(s)
#' @param digits Significant digits (default 3)
#' @return Character string(s)
#' @keywords internal
format_pval <- function(x, digits = 3) {
  formatC(x, digits = digits, format = "g")
}

#' Format a data frame for tabular display
#'
#' Rounds ordinary numeric columns to \code{digits} decimal places.
#' Columns whose names match common p-value patterns are formatted with
#' \code{\link{format_pval}} so very small values appear in scientific
#' notation rather than being displayed as "0".
#'
#' @param df Data frame
#' @param digits Decimal places for regular numeric columns (default 4)
#' @return Data frame with formatted columns
#' @keywords internal
format_df_for_display <- function(df, digits = 4) {
  pval_pattern <- "p[._]?value|p[._]?adj|padj|^p$|Pr\\("
  for (col in names(df)) {
    if (!is.numeric(df[[col]])) next
    if (grepl(pval_pattern, col, ignore.case = TRUE)) {
      df[[col]] <- format_pval(df[[col]])
    } else {
      df[[col]] <- round(df[[col]], digits)
    }
  }
  df
}

#' Create a summary table from analysis results
#'
#' @param result Analysis result object
#'
#' @return Data frame with summary
#' @keywords internal
extract_summary_table <- function(result) {
  if (is.null(result) || isTRUE(result$error)) {
    return(NULL)
  }
  
  # Extract relevant summary information
  # This will be customized based on analysis type
  if (!is.null(result$summary)) {
    return(as.data.frame(result$summary))
  }
  
  NULL
}
