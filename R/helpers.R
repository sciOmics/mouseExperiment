#' Helper Functions for mouseExperimentDashboard
#' 
#' Utility functions for data validation, formatting, and processing
#'

# ---- Bayesian progress modal ------------------------------------------------

#' Show a modal dialog while a Bayesian fit is running
#'
#' Uses shiny::showModal() which flushes to the browser immediately, before
#' the future_promise() computation starts. Dismiss with shiny::removeModal().
#'
#' @param message Title line shown in the dialog.
#' @param detail  Subtitle / detail line.
#' @keywords internal
show_analysis_modal <- function(
    message = "Analysis running…",
    detail  = "Results will appear automatically when complete.") {
  shiny::showModal(shiny::modalDialog(
    title     = NULL,
    footer    = NULL,
    easyClose = FALSE,
    size      = "s",
    shiny::tags$div(
      style = "text-align:center; padding:24px 8px 16px;",
      shiny::tags$i(
        class = "fa fa-spinner fa-spin fa-2x",
        style = "color:var(--primary-color); margin-bottom:16px; display:block;"
      ),
      shiny::tags$h4(message, style = "font-weight:600; margin-bottom:16px;"),
      shiny::tags$div(
        class = "progress",
        style = "height:6px; border-radius:3px; margin-bottom:16px;",
        shiny::tags$div(
          class = "progress-bar progress-bar-striped active",
          role  = "progressbar",
          style = "width:100%;"
        )
      ),
      shiny::tags$p(detail, style = "color:#888; font-size:0.9em; margin-bottom:0;")
    )
  ))
}


# ---- Subset filter (CODE_REVIEW.md Round 2 F.5) -----------------------------

#' Apply a global subset filter to uploaded data
#'
#' Reads exclusion lists from \code{rv$exclude_ids},
#' \code{rv$exclude_treatments}, and \code{rv$exclude_cages} (set by the
#' Data Upload tab's Subset Filter panel) and returns the data frame with
#' those rows dropped. Quietly returns the input unchanged when no
#' exclusions are active or the relevant columns aren't mapped.
#'
#' Every analysis module calls this on its uploaded data as the first
#' step, so users can re-run any analysis with a different subset by
#' simply editing the Data Upload filters and clicking Run — no need to
#' re-upload the CSV. CODE_REVIEW.md Round 2 F.5.
#'
#' @param df Raw uploaded data frame.
#' @param rv Shared reactiveValues list.
#' @return Filtered data frame.
#' @keywords internal
apply_subset_filter <- function(df, rv) {
  if (is.null(df) || nrow(df) == 0L) return(df)

  ex_ids <- shiny::isolate(rv$exclude_ids)        %||% character(0)
  ex_tx  <- shiny::isolate(rv$exclude_treatments) %||% character(0)
  ex_cg  <- shiny::isolate(rv$exclude_cages)      %||% character(0)

  id_col <- shiny::isolate(rv$settings$col_id)
  tx_col <- shiny::isolate(rv$settings$col_treatment)
  cg_col <- shiny::isolate(rv$settings$col_cage)

  if (length(ex_ids) > 0L && !is.null(id_col) && nzchar(id_col) &&
      id_col %in% colnames(df)) {
    df <- df[!(as.character(df[[id_col]]) %in% ex_ids), , drop = FALSE]
  }
  if (length(ex_tx) > 0L && !is.null(tx_col) && nzchar(tx_col) &&
      tx_col %in% colnames(df)) {
    df <- df[!(as.character(df[[tx_col]]) %in% ex_tx), , drop = FALSE]
  }
  if (length(ex_cg) > 0L && !is.null(cg_col) && nzchar(cg_col) &&
      cg_col %in% colnames(df)) {
    df <- df[!(as.character(df[[cg_col]]) %in% ex_cg), , drop = FALSE]
  }
  df
}

# ---- Alert helpers (consolidated validate/need pattern) ----------------------

#' Create a Bootstrap alert div for module status messages
#'
#' @param msg  Character message to display.
#' @param type One of "danger", "warning", "info", "success".
#' @return A \code{shiny::div} tag.
#' @keywords internal
status_alert <- function(msg, type = "danger") {
  icon_name <- switch(type,
    danger  = "exclamation-circle",
    warning = "exclamation-triangle",
    info    = "spinner",
    success = "check",
    "info-circle"
  )
  icon_cls <- if (type == "info") paste(icon_name, "fa-spin") else icon_name
  shiny::div(class = paste0("alert alert-", type),
             shiny::icon(icon_cls), msg)
}

# ---- Date Preprocessing (consolidated — Issue 2.3) ---------------------------

#' Detect the date format of a column and normalise data for calculate_dates()
#'
#' Supports four patterns:
#'   MM/DD/YYYY, YYYY-MM-DD, DD-Mon, and "Month Day" (e.g. "January 5").
#' For the "Month Day" pattern the year is inferred from \code{reference_date}
#' and both the data column and \code{reference_date} are reformatted so that
#' \code{mouseExperiment::calculate_dates()} can parse them.
#'
#' @param data          Data frame containing the date column.
#' @param date_column   Name of the date column.
#' @param reference_date Character string of the reference / start date.
#' @return A list with components \code{data} (possibly modified data frame),
#'   \code{date_format} (format string or \code{NULL}), and
#'   \code{start_date} (possibly modified reference date string).
#' @keywords internal
detect_date_format <- function(data, date_column, reference_date) {
  first_val <- as.character(data[1, date_column])
  start_date <- reference_date

  # Parse reference_date from its canonical YYYY-MM-DD or common dashboard formats
  ref_parsed_base <- tryCatch(
    as.Date(start_date,
            tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d-%b-%Y",
                           "%B %d %Y", "%B %d, %Y")),
    error = function(e) NA_Date_
  )
  if (is.na(ref_parsed_base)) ref_parsed_base <- Sys.Date()

  if (grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", first_val)) {
    fmt <- "%m/%d/%Y"
    start_date <- format(ref_parsed_base, "%m/%d/%Y")
  } else if (grepl("^\\d{4}-\\d{1,2}-\\d{1,2}$", first_val)) {
    fmt <- "%Y-%m-%d"
    start_date <- format(ref_parsed_base, "%Y-%m-%d")
  } else if (grepl("^\\d{1,2}-[A-Za-z]{3}$", first_val)) {
    fmt <- "%d-%b"
    start_date <- format(ref_parsed_base, "%d-%b")
  } else if (grepl("^[A-Za-z]+\\s+\\d{1,2}$", first_val)) {
    ref_year <- as.integer(format(ref_parsed_base, "%Y"))
    data[[date_column]] <- paste(data[[date_column]], ref_year)
    start_date <- format(ref_parsed_base, "%B %d %Y")
    fmt <- "%B %d %Y"
  } else {
    fmt <- NULL
  }

  list(data = data, date_format = fmt, start_date = start_date)
}

# ---- End Date Preprocessing -------------------------------------------------

#' Ensure data has been preprocessed (dates and volumes computed)
#'
#' If the canonical preprocessing from the Tumor Growth module has already run
#' (i.e. \code{rv$processed_data} is non-NULL), returns it along with the
#' stored column names. Otherwise performs inline preprocessing using
#' \code{rv$settings} and stores the result in \code{rv} for future callers.
#'
#' @param rv The shared reactiveValues list.
#' @return A list with \code{data}, \code{time_col}, \code{vol_col}, \code{tx_col}.
#' @keywords internal
ensure_preprocessed <- function(rv) {
  # Reuse canonical result if available. CODE_REVIEW.md Round 2 F.5:
  # always re-apply the Subset Filter to the cached frame so changes to
  # the filter take effect on the next analysis run without needing a
  # full re-upload + re-preprocess.
  if (!is.null(rv$processed_data)) {
    return(list(
      data     = apply_subset_filter(rv$processed_data, rv),
      time_col = rv$time_col,
      vol_col  = rv$vol_col,
      tx_col   = rv$tx_col
    ))
  }

  # Otherwise preprocess from scratch. Apply the Subset Filter to the
  # raw upload before any column derivation.
  data <- apply_subset_filter(rv$uploaded_data, rv)
  settings <- rv$settings

  # Date → Day
  date_fmt_warning <- NULL
  if (!is.null(settings$time_method) && settings$time_method == "from_date") {
    date_info <- detect_date_format(data, settings$col_date_for_day,
                                    settings$reference_date)
    if (is.null(date_info$date_format)) {
      date_fmt_warning <- paste0(
        "Could not auto-detect date format in column '",
        settings$col_date_for_day,
        "'. Please verify your date column or pre-compute the Day column."
      )
    }
    data <- mouseExperiment::calculate_dates(
      date_info$data,
      start_date  = date_info$start_date,
      date_column = settings$col_date_for_day,
      date_format = date_info$date_format
    )
    time_col <- "Day"
  } else {
    time_col <- settings$col_day
  }

  # Volume calculation
  if (!is.null(settings$volume_method) && settings$volume_method == "calculate") {
    len_col <- settings$col_length %||% "Length"
    wid_col <- settings$col_width  %||% "Width"
    vol_formula <- settings$volume_formula %||% "modified_ellipsoid"
    ht_col <- if (!is.null(settings$col_height) && nzchar(settings$col_height)) settings$col_height else NULL
    if (has_col(data, len_col) && has_col(data, wid_col)) {
      data <- mouseExperiment::calculate_volume(
        data,
        length_column = len_col,
        width_column  = wid_col,
        height_column = ht_col,
        formula       = vol_formula
      )
    }
    vol_col <- "Volume"
  } else {
    vol_col <- settings$col_volume %||% "Volume"
  }

  tx_col <- settings$col_treatment

  # Resolve optional columns (id, cage, dose) so other modules don't

  # depend on TG analysis having been run first.
  rv$id_col   <- if (!is.null(settings$col_id)   && nzchar(settings$col_id))   settings$col_id   else NULL
  rv$cage_col <- if (!is.null(settings$col_cage)  && nzchar(settings$col_cage)) settings$col_cage  else NULL
  rv$dose_col   <- if (!is.null(settings$col_dose)   && nzchar(settings$col_dose))   settings$col_dose   else NULL
  rv$weight_col <- if (!is.null(settings$col_weight) && nzchar(settings$col_weight)) settings$col_weight else NULL
  rv$sex_col    <- if (!is.null(settings$col_sex)    && nzchar(settings$col_sex))    settings$col_sex    else NULL

  # Cache for future callers
  rv$processed_data <- data
  rv$time_col       <- time_col
  rv$vol_col        <- vol_col
  rv$tx_col         <- tx_col

  list(data = data, time_col = time_col, vol_col = vol_col, tx_col = tx_col,
       date_fmt_warning = date_fmt_warning)
}

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
    # "mouse tag" / "ear tag" / "tag" cover non-standard ID column names
    id = detect_column(col_names, col_names_lower,
                       c("id", "mouse_id", "animal_id", "subject",
                         "mouse tag", "ear tag", "tag")),
    day = detect_column(col_names, col_names_lower,
                        c("day", "time", "timepoint", "day_number")),
    # "tumor volume" (space) covers columns named without underscore
    volume = detect_column(col_names, col_names_lower,
                           c("volume", "tumor_volume", "tumor volume", "tv")),
    # "tumor length" prefix catches "Tumor Length (mm)" via startsWith fallback
    length = detect_column(col_names, col_names_lower,
                           c("length", "tumor_length", "tumor length")),
    # "tumor width" prefix catches "Tumor Width (mm)" via startsWith fallback
    width  = detect_column(col_names, col_names_lower,
                           c("width", "tumor_width", "tumor width")),
    treatment = detect_column(col_names, col_names_lower,
                              c("treatment", "group", "treatment_group")),
    cage = detect_column(col_names, col_names_lower,
                         c("cage", "cage_id")),
    date = detect_column(col_names, col_names_lower,
                         c("date", "measurement_date")),
    dose = detect_column(col_names, col_names_lower,
                         c("dose", "dose_level", "concentration", "conc")),
    survival_time  = detect_column(col_names, col_names_lower,
                                   c("survival_time", "time_to_event", "days_to_event")),
    survival_event = detect_column(col_names, col_names_lower,
                                    c("survival_event", "event", "survival_censor", "censor")),
    weight = detect_column(col_names, col_names_lower,
                           c("weight", "body_weight", "body weight", "mouse_weight",
                             "mouse weight", "mass", "body_mass", "body mass", "bw")),
    sex    = detect_column(col_names, col_names_lower,
                           c("sex", "gender", "m/f")),
    necrotic = detect_column(col_names, col_names_lower,
                             c("necrotic", "necrosis", "necr",
                               "is_necrotic", "necrotic_flag"))
  )

  detected
}

#' Helper to find column in names
#' @keywords internal
detect_column <- function(col_names, col_names_lower, patterns) {
  # 1. Exact match (case-insensitive)
  for (pattern in patterns) {
    match_idx <- which(col_names_lower == pattern)
    if (length(match_idx) > 0) return(col_names[match_idx[1]])
  }
  # 2. Contains match — handles columns with prefixes/suffixes/units, e.g.
  #    "Tumor Length (mm)" or "Mean Tumor Volume" match their respective patterns
  for (pattern in patterns) {
    match_idx <- which(grepl(pattern, col_names_lower, fixed = TRUE))
    if (length(match_idx) > 0) return(col_names[match_idx[1]])
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
#' @param module Optional module name for structured logging.
#'
#' @return Result or error message
#' @keywords internal
translate_error <- function(msg) {
  if (grepl("not found in|not found$|object .* not found", msg, perl = TRUE))
    return("A required column was not found. Check your column mappings in Data Upload.")
  if (grepl("singular|singular fit", msg, ignore.case = TRUE))
    return(paste0(
      "Model fitting failed (singular fit). ",
      "Try increasing sample size or simplifying the random effects structure."
    ))
  if (grepl("subscript out of bounds", msg, fixed = TRUE))
    return(paste0(
      "An internal indexing error occurred. ",
      "Check that your data is complete and all required columns are mapped."
    ))
  if (grepl("non-conformable", msg, fixed = TRUE))
    return("Matrix dimension mismatch — check that all treatment groups have sufficient observations.")
  if (grepl("NaN|NA|missing values", msg, ignore.case = TRUE))
    return(paste0(
      "The model encountered missing or invalid values. ",
      "Check your data for NA values or zero-variance groups."
    ))
  msg
}

toggle_downloads <- function(session, ids, enabled) {
  fn <- if (isTRUE(enabled)) shinyjs::enable else shinyjs::disable
  invisible(lapply(ids, function(id) fn(session$ns(id))))
}

safe_analysis <- function(expr, module = "unknown") {
  t0 <- proc.time()[["elapsed"]]
  tryCatch(
    {
      result <- expr
      elapsed <- proc.time()[["elapsed"]] - t0
      log_info(sprintf("analysis_complete | module=%s | elapsed_s=%.2f", module, round(elapsed, 2)))
      result
    },
    error = function(e) {
      elapsed <- proc.time()[["elapsed"]] - t0
      raw_msg <- conditionMessage(e)
      log_error(sprintf("analysis_failed | module=%s | elapsed_s=%.2f | error=%s",
                        module, round(elapsed, 2), raw_msg))
      list(error = TRUE, message = translate_error(raw_msg))
    }
  )
}

# ---------------------------------------------------------------------------
# Validation helpers (3.4 / 3.5)
# ---------------------------------------------------------------------------

#' Check whether a column exists in a data frame
#'
#' Convenience guard used throughout dashboards before accessing columns
#' by name. Returns \code{FALSE} when \code{col} is NULL, empty string,
#' or not present in \code{colnames(data)}.
#'
#' @param data A data frame (or anything with \code{colnames}).
#' @param col  Column name (character scalar).
#' @return Logical scalar.
#' @keywords internal
has_col <- function(data, col) {
  !is.null(col) && nzchar(col) && col %in% colnames(data)
}

js_string <- function(x) {
  jsonlite::toJSON(as.character(x %||% ""), auto_unbox = TRUE)
}

#' Validate that required columns exist in a data frame
#'
#' Returns a character vector of missing column names (empty if all present).
#' Designed for use with \code{validate(need(length(missing) == 0, ...))}.
#'
#' @param data A data frame.
#' @param cols Named list or character vector of column names.
#'   Names (if any) are used in the error message.
#' @return Character vector of missing column names.
#' @keywords internal
check_required_cols <- function(data, cols) {
  cn <- colnames(data)
  missing <- vapply(cols, function(col) {
    is.null(col) || !nzchar(col) || !(col %in% cn)
  }, logical(1))
  cols[missing]
}

#' Safely generate a plot, returning an error placeholder on failure
#'
#' Wraps a plot expression in tryCatch. On success returns the plot object.
#' On failure logs via structured logger and returns a simple ggplot with the
#' error text, so users see a message instead of a blank panel.
#'
#' @param expr    Expression that produces a ggplot / base plot grob.
#' @param label   Short label for the plot (used in the error message).
#' @return A plot object.
#' @keywords internal
safe_plot <- function(expr, label = "Plot") {
  tryCatch(
    expr,
    error = function(e) {
      msg <- paste0(label, " generation failed: ", conditionMessage(e))
      log_caught_error("plot", label, e)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                          size = 5, hjust = 0.5) +
        ggplot2::theme_void()
    }
  )
}

#' Plot export dimension controls
#'
#' Generates Width (px), Height (px), and Resolution (DPI) numeric inputs that
#' are placed at the bottom of each analysis sidebar. Download handlers read
#' these with \code{input[[paste0(prefix, "_plot_width")]]} etc.
#'
#' @param prefix Input ID prefix, e.g. \code{"tg"}, \code{"surv"}.
#' @keywords internal
plot_export_ui <- function(prefix) {
  tagList(
    hr(style = "margin-top: 14px; margin-bottom: 8px;"),
    p(style = "font-weight: 600; margin-bottom: 4px; font-size: 0.88em;",
      "\U0001F4BE  Plot Export Settings"),
    fluidRow(
      column(6,
        numericInput(paste0(prefix, "_plot_width"),  "Width (px)",
                     value = 800, min = 200, max = 5000, step = 50)
      ),
      column(6,
        numericInput(paste0(prefix, "_plot_height"), "Height (px)",
                     value = 600, min = 200, max = 5000, step = 50)
      )
    ),
    numericInput(paste0(prefix, "_plot_dpi"), "Resolution (DPI)",
                 value = 100, min = 72, max = 600, step = 10)
  )
}

#' Back-transform values from analysis scale to original scale
#'
#' Inverse of the transform applied in `tumor_growth_statistics()` and the
#' Bayesian / GAM tumor-growth fits, which all use `log(Volume)` (not
#' `log1p`). For `"sqrt"`, returns `x^2`. For `"none"`, returns `x` unchanged.
#'
#' @param x Numeric vector on the transformed scale
#' @param transform One of "log", "sqrt", or "none"
#' @return Numeric vector on the original scale
#' @keywords internal
back_transform <- function(x, transform) {
  x <- suppressWarnings(as.numeric(x))
  switch(transform,
    "log"  = exp(x),
    "sqrt" = x ^ 2,
    "none" = x,
    x
  )
}

#' Axis-label suffix for an analysis-scale plot
#'
#' Returns " (log scale)", " (sqrt scale)", or "" so plot y-axis labels can
#' indicate that values are on the transformed scale without overriding any
#' user-supplied label.
#'
#' @param transform One of "log", "sqrt", or "none"
#' @return Character scalar suitable for appending to an axis label
#' @keywords internal
scale_suffix <- function(transform) {
  switch(transform,
    "log"  = " (log scale)",
    "sqrt" = " (sqrt scale)",
    ""
  )
}

#' Apply custom title, subtitle, text size, and alignment to a ggplot
#'
#' When title/subtitle are empty strings the existing labs() are preserved
#' (i.e. the hardcoded plot title is not removed).
#'
#' @param p A ggplot2 object.
#' @param title Character; if non-empty overrides the plot title.
#' @param subtitle Character; if non-empty adds a subtitle.
#' @param align One of "left", "center", "right" — alignment for title/subtitle.
#' @param text_size Numeric; base font size in pt.
#' @param xlab Character; if non-empty overrides the X-axis label.
#' @param ylab Character; if non-empty overrides the Y-axis label.
#' @param legend_lab Character; if non-empty overrides the legend title.
#' @return A ggplot2 object with the theme layers applied.
#' @keywords internal
apply_plot_labels <- function(p, title = "", subtitle = "", align = "center",
                              text_size = 12, xlab = "", ylab = "", legend_lab = "") {
  hjust <- switch(align, left = 0, center = 0.5, right = 1, 0.5)
  if (nzchar(trimws(title)))    p <- p + ggplot2::labs(title = title)
  if (nzchar(trimws(subtitle))) p <- p + ggplot2::labs(subtitle = subtitle)
  if (nzchar(trimws(xlab)))     p <- p + ggplot2::labs(x = xlab)
  if (nzchar(trimws(ylab)))     p <- p + ggplot2::labs(y = ylab)
  if (nzchar(trimws(legend_lab))) p <- p + ggplot2::labs(fill = legend_lab, color = legend_lab, shape = legend_lab)
  p + ggplot2::theme(
    text          = ggplot2::element_text(size = text_size),
    plot.title    = ggplot2::element_text(hjust = hjust),
    plot.subtitle = ggplot2::element_text(hjust = hjust)
  )
}

#' Format a p-value for display
#'
#' Very small p-values are shown in scientific notation rather than being
#' rounded to "0". Values that underflow to exactly 0 are reported as
#' "< 2.2e-16" (the smallest representable positive double). Uses
#' \code{formatC} with \code{format="g"} which automatically switches to
#' scientific notation for small values.
#'
#' @param x Numeric p-value(s)
#' @param digits Significant digits (default 3)
#' @return Character string(s)
#' @keywords internal
format_pval <- function(x, digits = 3) {
  limit <- .Machine$double.eps  # ~2.2e-16
  ifelse(
    is.na(x), NA_character_,
    ifelse(
      x == 0 | (!is.na(x) & x < limit),
      paste0("< ", formatC(limit, digits = 2, format = "e")),
      formatC(x, digits = digits, format = "g")
    )
  )
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
      # Log and leave column as-is if formatting fails
      log_caught_error("format", paste0("format_df_for_display col '", col, "'"), e)
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
  # Return NULL unless ALL groups have values — partial coverage causes
  # ggplot2 scale_colour_manual to map missing groups as grey/invisible.
  if (any(is.na(colors))) return(NULL)
  colors
}

#' Remap colour vector names from dose-composite to treatment-group names
#'
#' When groups are stored as \code{"Treatment - Dose: X"}, colour names
#' include the \code{" - Dose: "} separator.  Plot data (growth_rates, AUC) uses
#' the \code{TreatmentGroup} column which pastes treatment and dose with a
#' single space (e.g. \code{"Control 0"}).  This function converts
#' \code{" - Dose: "} to \code{" "} and deduplicates (keeping first).
#'
#' @param colors Named character vector (or NULL).
#' @return Named character vector with treatment-group names, or NULL.
#' @keywords internal
remap_colors_to_treatment <- function(colors) {
  if (is.null(colors)) return(NULL)
  nms <- gsub(" - Dose: ", " ", names(colors))
  dups <- duplicated(nms)
  out <- colors[!dups]
  names(out) <- nms[!dups]
  out
}

#' Remap shape vector names from dose-composite to treatment-group names
#'
#' Same logic as \code{remap_colors_to_treatment} but for the integer shape vector.
#'
#' @param shapes Named integer vector (or NULL).
#' @return Named integer vector with treatment-level names, or NULL.
#' @keywords internal
remap_shapes_to_treatment <- function(shapes) {
  if (is.null(shapes)) return(NULL)
  nms <- gsub(" - Dose: ", " ", names(shapes))
  dups <- duplicated(nms)
  out <- shapes[!dups]
  names(out) <- nms[!dups]
  out
}

#' Build a named integer vector of per-group point shapes from Shiny inputs
#'
#' Reads \code{input[[paste0(prefix, "_", make.names(group))]]} for each group.
#' Expected values are pch integers as character strings ("16", "15", etc.).
#'
#' @param input  Shiny \code{input} object.
#' @param prefix Input ID prefix, e.g. \code{"tg_shape"}.
#' @param groups Character vector of group names.
#' @return Named integer vector of pch codes, or NULL if none set.
#' @keywords internal
build_shape_vector <- function(input, prefix, groups) {
  if (is.null(groups) || length(groups) == 0L) return(NULL)
  shapes <- vapply(groups, function(g) {
    val <- input[[paste0(prefix, "_", make.names(g))]]
    if (is.null(val) || !nzchar(val)) NA_integer_ else as.integer(val)
  }, integer(1L))
  # Return NULL unless ALL groups have values — partial coverage causes
  # ggplot2 scale_shape_manual to map missing groups as NA shape (invisible).
  if (any(is.na(shapes))) return(NULL)
  shapes
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

#' Named palette presets (colorblind-friendly options)
#'
#' Returns a named character vector of hex colours for the requested palette.
#' Falls back to the ggplot2 default hue palette when \code{name} is
#' \code{"default"} or unrecognised.
#'
#' @param name Character palette name: \code{"default"}, \code{"viridis"},
#'   \code{"okabe_ito"}, \code{"set2"}, or \code{"wes_anderson"}.
#' @param n    Number of colours required.
#' @return Character vector of hex colour strings of length \code{n}.
#' @keywords internal
palette_colors <- function(name, n) {
  if (n <= 0L) return(character(0L))
  switch(name,
    viridis  = substring(scales::viridis_pal(option = "D")(n), 1L, 7L),
    okabe_ito = {
      oi <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
              "#0072B2", "#D55E00", "#CC79A7", "#999999")
      rep_len(oi, n)
    },
    set2 = {
      s2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
              "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
      rep_len(s2, n)
    },
    wes_anderson = {
      wa <- c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
      rep_len(wa, n)
    },
    default_group_colors(n)
  )
}

#' Set up an observer that applies palette presets to group colour pickers
#'
#' Call this once inside a module server.
#' When the user selects a non-custom preset, every per-group \code{colourInput}
#' is updated to the preset colour.  Selecting \code{"custom"} is a no-op so
#' the user can freely pick colours after choosing a starting palette.
#'
#' @param input   Shiny \code{input} object (module-scoped).
#' @param session Shiny \code{session} object (module-scoped).
#' @param prefix  Same \code{prefix} passed to \code{group_order_and_color_ui}.
#' @param groups_rv A reactive expression returning the current group names.
#' @keywords internal
observe_palette_preset <- function(input, session, prefix, groups_rv) {
  pal_id <- paste0(prefix, "_palette")
  observeEvent(input[[pal_id]], {
    pal_name <- input[[pal_id]]
    if (is.null(pal_name) || pal_name == "custom") return()
    groups <- groups_rv()
    if (is.null(groups) || length(groups) == 0L) return()
    cols <- palette_colors(pal_name, length(groups))
    for (i in seq_along(groups)) {
      color_id <- paste0(prefix, "_color_", make.names(groups[i]))
      colourpicker::updateColourInput(session, color_id, value = cols[i])
    }
  }, ignoreInit = TRUE)
}

#' Render colourInput widgets for a set of treatment groups
#'
#' Returns a \code{tagList} of \code{colourpicker::colourInput} widgets,
#' one per group, with auto-assigned ggplot2-style default colours.
#'
#' @param prefix    Input-ID prefix (e.g. \code{"tg_color"}).
#' @param groups    Character vector of group names.
#' @param header    Optional header string. Default \code{"Group Colors"}.
#' @return A \code{tagList}, or \code{NULL} when \code{groups} is empty.
#' @keywords internal
group_color_pickers_ui <- function(prefix, groups, header = "Group Colors") {
  if (is.null(groups) || length(groups) == 0L) return(NULL)
  defaults <- default_group_colors(length(groups))
  shiny::tagList(
    shiny::tags$h6(style = "font-size: 14px; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                   shiny::icon("palette"), " ", header),
    lapply(seq_along(groups), function(i) {
      colourpicker::colourInput(
        inputId = paste0(prefix, "_", make.names(groups[i])),
        label   = groups[i],
        value   = defaults[i]
      )
    })
  )
}

#' Synthesize a survival/censor column from dropout patterns
#'
#' When no censor column is available, this function infers which mice died
#' before the study end by checking whether each mouse's last observed timepoint
#' is before the global maximum timepoint. For each such mouse, the row at its
#' last observed timepoint is flagged with \code{1} (dead); all other rows
#' receive \code{0} (alive/censored-at-study-end).
#'
#' @param data        A data frame of per-mouse longitudinal measurements.
#' @param id_col      Name of the column containing mouse/subject IDs.
#' @param time_col    Name of the column containing the timepoint (e.g. Day).
#' @param new_col     Name to give the synthesized column. Default \code{"_auto_censor"}.
#' @return A list with:
#'   \describe{
#'     \item{\code{data}}{The input data frame with the new column appended.}
#'     \item{\code{col}}{The name of the new column.}
#'     \item{\code{n_flagged}}{Number of mice flagged as early dropouts.}
#'   }
#' @keywords internal
synthesize_censor_column <- function(data, id_col, time_col,
                                     new_col = "_auto_censor") {
  data[[new_col]] <- 0L
  study_end <- max(data[[time_col]], na.rm = TRUE)

  ids <- unique(data[[id_col]])
  n_flagged <- 0L
  for (id in ids) {
    idx <- which(data[[id_col]] == id)
    last_t <- max(data[[time_col]][idx], na.rm = TRUE)
    if (last_t < study_end) {
      # Mark only the very last row for this mouse
      last_idx <- idx[which.max(data[[time_col]][idx])]
      data[[new_col]][last_idx] <- 1L
      n_flagged <- n_flagged + 1L
    }
  }
  list(data = data, col = new_col, n_flagged = n_flagged)
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
    text     = "Group order",
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
#' @param pos Character; one of "none", "top", "bottom", "left",
#'   "right", "inplot_tl" (top-left), "inplot_tr" (top-right),
#'   "inplot_bl" (bottom-left), or "inplot_br" (bottom-right).
#' @return The modified ggplot2 object.
#' @keywords internal
apply_legend_pos <- function(p, pos) {
  if (is.null(pos)) return(p)
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
  } else if (pos == "inplot_bl") {
    p + ggplot2::theme(
      legend.position    = c(0.02, 0.02),
      legend.justification = c("left", "bottom"),
      legend.background  = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), colour = NA)
    )
  } else if (pos == "inplot_br") {
    p + ggplot2::theme(
      legend.position    = c(0.98, 0.02),
      legend.justification = c("right", "bottom"),
      legend.background  = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), colour = NA)
    )
  } else {
    p + ggplot2::theme(legend.position = pos)
  }
}

#' Render a drag-and-drop rank list with embedded color pickers
#'
#' Combines group ordering (drag-and-drop) with color selection in a single
#' control. Groups can be reordered by dragging, and colors are adjusted via
#' inline colour pickers. The order is stored in \code{input[[paste0(prefix, "_order")]]}
#' and colors are read from \code{input[[paste0(prefix, "_color_", make.names(group))]]}
#' for each group.
#'
#' @param prefix    Input-ID prefix (e.g. \code{"tg"}).
#' @param groups    Character vector of group names in their default order.
#' @return A \code{tagList} containing HTML structure and colour inputs, or
#'   \code{NULL} when fewer than 2 groups.
#' @keywords internal
group_order_and_color_ui <- function(prefix, groups, ns = NS(NULL)) {
  if (is.null(groups) || length(groups) < 2L) return(NULL)

  defaults <- default_group_colors(length(groups))
  list_id  <- ns(paste0(prefix, "_drag_list"))
  input_id <- ns(paste0(prefix, "_order"))
  pal_id   <- ns(paste0(prefix, "_palette"))

  # Closed pch shapes available for per-group selection.
  # Open shapes (0, 1, 2, ...) are reserved for extrapolated points.
  shape_choices <- c(
    "\u25cf Circle"        = "16",
    "\u25a0 Square"        = "15",
    "\u25b2 Triangle \u2191" = "17",
    "\u25bc Triangle \u2193" = "25",
    "\u25c6 Diamond"       = "18"
  )
  default_shapes <- c(16L, 15L, 17L, 25L, 18L)

  # Each row: drag handle + group label + colour picker + shape selector
  # data-value holds the exact group name string — read directly by JS,
  # no innerText parsing, no custom Shiny input type.
  items <- lapply(seq_along(groups), function(i) {
    g        <- groups[i]
    color_id <- ns(paste0(prefix, "_color_", make.names(g)))
    shape_id <- ns(paste0(prefix, "_shape_", make.names(g)))
    def_shape <- as.character(default_shapes[((i - 1L) %% length(default_shapes)) + 1L])
    shiny::tags$li(
      `data-value` = g,
      style = paste0(
        "display:flex; align-items:center; gap:6px;",
        "padding:5px 8px; margin-bottom:4px;",
        "background:#f8f9fa; border:1px solid #dee2e6;",
        "border-radius:4px; cursor:default;"
      ),
      shiny::tags$span(
        class = "drag-handle",
        style = "cursor:grab; color:#adb5bd; font-size:14px; user-select:none; flex-shrink:0;",
        "\u2630"   # grip / hamburger icon
      ),
      shiny::tags$span(
        style = "flex:1; font-size:13px; min-width:55px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;",
        g
      ),
      colourpicker::colourInput(
        inputId          = color_id,
        label            = NULL,
        value            = defaults[i],
        showColour       = "background",
        palette          = "square",
        allowTransparent = FALSE
      ),
      shiny::tags$div(
        style = "flex-shrink:0;",
        shiny::selectInput(
          inputId  = shape_id,
          label    = NULL,
          choices  = shape_choices,
          selected = def_shape,
          width    = "125px"
        )
      )
    )
  })

  # JS helper: read data-value attributes from the sortable <li> elements
  # and push the current order to Shiny.  Uses getElementById for reliability
  # (avoids depending on SortableJS `this.el` inside async callbacks).
  send_order_js <- htmlwidgets::JS(sprintf(
    "function(evt) {
       var el = document.getElementById('%s');
       if (el && typeof Shiny !== 'undefined' && Shiny.setInputValue) {
         var order = Array.prototype.slice.call(el.children).map(
           function(li) { return li.getAttribute('data-value'); }
         );
         Shiny.setInputValue('%s', order, {priority: 'event'});
       }
     }",
    list_id, input_id
  ))

  # Fire once after the sortable is rendered so the initial order reaches Shiny
  # (onLoad is NOT a standard SortableJS event and is silently ignored).
  init_js <- sprintf(
    "$(function(){
       var el = document.getElementById('%s');
       if (el && typeof Shiny !== 'undefined' && Shiny.setInputValue) {
         var order = Array.prototype.slice.call(el.children).map(
           function(li) { return li.getAttribute('data-value'); }
         );
         Shiny.setInputValue('%s', order, {priority: 'event'});
       }
     });",
    list_id, input_id
  )

  shiny::tagList(
    shiny::selectInput(
      pal_id, shiny::tags$span(shiny::icon("eye"), " Color Palette Presets"),
      choices  = c("Default" = "custom",
                   "Viridis" = "viridis",
                   "Okabe-Ito" = "okabe_ito",
                   "Set2" = "set2",
                   "Wes Anderson" = "wes_anderson"),
      selected = "custom",
      width    = "100%"
    ),
    shiny::tags$p(
      style = "font-size:11px; color:#6c757d; margin-bottom:6px;",
      shiny::icon("grip-vertical"), " Drag the \u2630 handle to reorder"
    ),
    shiny::tags$ul(
      id    = list_id,
      style = "padding:0; margin:0; list-style:none;",
      items
    ),
    sortable::sortable_js(
      css_id  = list_id,
      options = sortable::sortable_options(
        handle    = ".drag-handle",
        animation = 150L,
        onSort    = send_order_js
      )
    ),
    shiny::tags$script(shiny::HTML(init_js))
  )
}

