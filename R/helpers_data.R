#' Data ingestion, detection, validation and preparation helpers
#'
#' Split out of the former `helpers.R` (1,227 LOC) under CODE_REVIEW.md D.1.
#'
#' The split is organised by *testability seam*, not by line count. These are the
#' functions that decide what the data means before any model sees it: which
#' column is the ID, what a date string is, which group is the reference, and --
#' most consequentially -- which animals count as events. A defect here is
#' invisible downstream because every number the user reads is already
#' conditioned on it.
#'
#' D.1 originally read as "these files are too long". Length was not the problem:
#' `mod_tumor_growth.R` had already been split once (2,659 -> 1,877 LOC) and D.7
#' still sat at line 7 of the refactored file. The problem is that 29 of the 39
#' functions in `helpers.R` had no test referencing them, and length made that
#' easy not to notice.

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
  # R20.D4: the handler returned `NA_Date_`, which does not exist, so any
  # reference date outside these formats raised "object 'NA_Date_' not found".
  # The intended fallback -- today's date -- would count every Day from today,
  # which is silently wrong, so an unusable reference date is reported instead.
  if (is.null(start_date) || !nzchar(trimws(start_date))) {
    stop("No reference date is set. Enter the study's reference date (Day 0) ",
         "as YYYY-MM-DD, for example 2025-03-24.", call. = FALSE)
  }
  ref_parsed_base <- tryCatch(
    as.Date(start_date,
            tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%d-%b-%Y",
                           "%B %d %Y", "%B %d, %Y")),
    error = function(e) as.Date(NA)
  )
  if (is.na(ref_parsed_base)) {
    stop("Reference date '", start_date, "' is not a recognised date. Enter ",
         "it as YYYY-MM-DD, for example 2025-03-24.", call. = FALSE)
  }

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
    # CSV only. R data files (.rda/.RData) were removed in v26.09.001
    # (CODE_REVIEW.md R20.D9): the path never worked, and it loaded uploaded
    # objects into the process-wide global environment.
    if (tolower(file_ext) == "csv") {
      data <- read.csv(file_path, stringsAsFactors = FALSE)
    } else {
      return(list(
        valid = FALSE,
        message = "Unsupported file format. Please upload a CSV file.",
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
