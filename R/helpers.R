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
  ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = "g"))
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
  # Flatten any array/table/matrix columns to plain atomic vectors first.
  # Columns returned by survival_statistics (Events, Total, Event_Rate) can be
  # table- or array-class objects; round() and datatable() both have trouble with
  # these non-standard types.
  for (col in names(df)) {
    v <- df[[col]]
    if (!is.null(dim(v)) || inherits(v, "table")) {
      df[[col]] <- as.vector(v)
    }
  }
  pval_pattern <- "p[._]?value|p[._]?adj|padj|^p$|Pr\\("
  for (col in names(df)) {
    if (!is.numeric(df[[col]])) next
    tryCatch({
      if (grepl(pval_pattern, col, ignore.case = TRUE)) {
        df[[col]] <- format_pval(df[[col]])
      } else {
        df[[col]] <- round(df[[col]], digits)
      }
    }, error = function(e) {
      # Leave the column as-is if formatting fails
    })
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
#' Build a named colour vector from colourInput values
#'
#' Reads shiny `input[[paste0(prefix, "_", make.names(group))]]` for each group
#' and returns a named character vector suitable for ggplot2's
#' `scale_color_manual(values = ...)` / `scale_fill_manual(values = ...)`.
#' Entries whose input value is NULL or empty are dropped.
#'
#' @param input  The Shiny `input` list.
#' @param prefix Character prefix used when naming the colourInput widgets,
#'   e.g. `"tg_color"` → input ID `"tg_color_<make.names(group)>"`.
#' @param groups Character vector of group names.
#' @return Named character vector of hex colours, or NULL if all are missing.
#' @keywords internal
build_color_vector <- function(input, prefix, groups) {
  if (is.null(groups) || length(groups) == 0L) return(NULL)
  colors <- vapply(groups, function(g) {
    val <- input[[paste0(prefix, "_", make.names(g))]]
    if (is.null(val) || !nzchar(val)) NA_character_ else val
  }, character(1L))
  colors <- colors[!is.na(colors)]
  if (length(colors) == 0L) NULL else colors
}

#' Generate n evenly-spaced ggplot2-style default colours for n groups
#'
#' @param n Integer number of colours.
#' @return Character vector of hex colour strings.
#' @keywords internal
default_group_colors <- function(n) {
  if (n <= 0L) return(character(0L))
  scales::hue_pal()(n)
}

#' Render colourInput widgets for a set of treatment groups
#'
#' Returns a \code{tagList} of \code{colourpicker::colourInput} widgets,
#' one per group, with auto-assigned ggplot2-style default colours.
#'
#' @param prefix    Input-ID prefix (e.g. \code{"tg_color"}).
#' @param groups    Character vector of group names.
#' @param header    Optional header string. Default \code{"Group Colours"}.
#' @return A \code{tagList}, or \code{NULL} when \code{groups} is empty.
#' @keywords internal
group_color_pickers_ui <- function(prefix, groups, header = "Group Colours") {
  if (is.null(groups) || length(groups) == 0L) return(NULL)
  defaults <- default_group_colors(length(groups))
  shiny::tagList(
    shiny::tags$h6(shiny::icon("palette"), " ", header),
    lapply(seq_along(groups), function(i) {
      colourpicker::colourInput(
        inputId = paste0(prefix, "_", make.names(groups[i])),
        label   = groups[i],
        value   = defaults[i]
      )
    })
  )
}

#' Render a drag-and-drop rank list for reordering treatment groups
#'
#' Returns a \code{sortable::rank_list} widget.  The current order is
#' read back from \code{input[[paste0(prefix, "_order")]]}.
#'
#' @param prefix    Input-ID prefix shared with colour pickers (e.g. \code{"tg"}).
#' @param groups    Character vector of group names in their default order.
#' @return A \code{sortable::rank_list} tag, or \code{NULL} when fewer than 2 groups.
#' @keywords internal
group_order_ui <- function(prefix, groups) {
  if (is.null(groups) || length(groups) < 2L) return(NULL)
  sortable::rank_list(
    text     = shiny::tags$span(shiny::icon("sort"), " Group order"),
    labels   = as.list(groups),
    input_id = paste0(prefix, "_order"),
    class    = "default-sortable"
  )
}

#' Return groups in the user-specified drag-and-drop order
#'
#' Reads \code{input[[paste0(prefix, "_order")]]} and validates it.
#' Falls back to \code{groups} (original order) if the input is absent or
#' does not contain exactly the same set of labels.
#'
#' @param input  Shiny \code{input} object.
#' @param prefix Input-ID prefix used when calling \code{group_order_ui}.
#' @param groups Character vector of group names in their default order.
#' @return Character vector of group names in the current user-selected order.
#' @keywords internal
get_ordered_groups <- function(input, prefix, groups) {
  if (is.null(groups)) return(groups)
  ord <- input[[paste0(prefix, "_order")]]
  if (!is.null(ord)) ord <- as.character(unlist(ord))
  if (!is.null(ord) &&
      length(ord) == length(groups) &&
      isTRUE(setequal(ord, groups))) {
    ord
  } else {
    groups
  }
}

#' Apply a legend position value (including in-plot coordinates) to a ggplot
#'
#' Translates the dashboard's legend position selector values into proper
#' ggplot2 theme settings, including the coordinate-based in-plot options.
#'
#' @param p   A ggplot2 object.
#' @param pos Character; one of "auto", "none", "top", "bottom", "left",
#'   "right", "inplot_tl" (top-left), or "inplot_tr" (top-right).
#' @return The modified ggplot2 object.
#' @keywords internal
apply_legend_pos <- function(p, pos) {
  if (is.null(pos) || pos == "auto") return(p)
  if (pos == "inplot_tl") {
    p + ggplot2::theme(
      legend.position    = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.background  = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), colour = NA)
    )
  } else if (pos == "inplot_tr") {
    p + ggplot2::theme(
      legend.position    = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background  = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), colour = NA)
    )
  } else {
    p + ggplot2::theme(legend.position = pos)
  }
}