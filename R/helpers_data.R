# Data ingestion, detection, validation and preparation helpers
#
# Split out of the former `helpers.R` (1,227 LOC) under CODE_REVIEW.md D.1.
#
# The split is organised by *testability seam*, not by line count. These are the
# functions that decide what the data means before any model sees it: which
# column is the ID, what a date string is, which group is the reference, and --
# most consequentially -- which animals count as events. A defect here is
# invisible downstream because every number the user reads is already
# conditioned on it.
#
# D.1 originally read as "these files are too long". Length was not the problem:
# `mod_tumor_growth.R` had already been split once (2,659 -> 1,877 LOC) and D.7
# still sat at line 7 of the refactored file. The problem is that 29 of the 39
# functions in `helpers.R` had no test referencing them, and length made that
# easy not to notice.

# ---- Subset filter (CODE_REVIEW.md Round 2 F.5) -----------------------------

#' One key per animal: treatment + ID (+ cage)
#'
#' Ear tags are commonly reused in each arm or cage, so an ID alone does not
#' identify an animal (CODE_REVIEW.md R20.D2, R20.D57). This is the dashboard's
#' copy of the package's internal `make_mouse_key()` convention (package T1):
#' the parts joined by "|||", treatment first.
#'
#' @param treatment,id Vectors of equal length.
#' @param cage Optional cage vector; NULL when no cage column is mapped.
#' @return Character vector of keys.
#' @keywords internal
animal_key <- function(treatment, id, cage = NULL) {
  parts <- list(as.character(treatment), as.character(id))
  if (!is.null(cage)) parts <- c(parts, list(as.character(cage)))
  do.call(paste, c(parts, sep = "|||"))
}

#' Human-readable animal label, e.g. "3 (Vehicle, cage C1)"
#' @inheritParams animal_key
#' @keywords internal
animal_label <- function(treatment, id, cage = NULL) {
  paste0(as.character(id), " (", as.character(treatment),
         if (!is.null(cage)) paste0(", cage ", as.character(cage)), ")")
}

#' Apply a global subset filter to uploaded data
#'
#' Reads exclusion lists from \code{rv$exclude_ids},
#' \code{rv$exclude_treatments}, and \code{rv$exclude_cages} (set by the
#' Data Upload tab's Subset Filter panel) and returns the data frame with
#' those rows dropped. Quietly returns the input unchanged when no
#' exclusions are active or the relevant columns aren't mapped.
#'
#' Excluded animals are \code{animal_key()} values (treatment + ID + cage):
#' excluding ear tag "3" used to drop animal 3 from every arm.
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
  args <- subset_filter_args(df, rv)
  if (is.null(args)) return(df)
  do.call(apply_exclusions, c(list(df = df), args))
}

#' The subset filter as [apply_exclusions()] arguments, or NULL when none is set
#'
#' Excluded animals are stored as [animal_key()]s; they are given to
#' apply_exclusions() as a data frame of arm, id and (with cages) cage.
#' Columns that are not mapped, or not in `df`, are left out, as before.
#' @param df The data the filter applies to.
#' @param rv Shared reactiveValues.
#' @keywords internal
subset_filter_args <- function(df, rv) {
  ex_ids <- shiny::isolate(rv$exclude_ids)        %||% character(0)
  ex_tx  <- shiny::isolate(rv$exclude_treatments) %||% character(0)
  ex_cg  <- shiny::isolate(rv$exclude_cages)      %||% character(0)
  if (!length(ex_ids) && !length(ex_tx) && !length(ex_cg)) return(NULL)
  s <- shiny::isolate(rv$settings)
  col <- function(x) if (has_col(df, x)) x
  cage <- col(s$col_cage)
  animals <- if (length(ex_ids)) {
    parts <- strsplit(ex_ids, "|||", fixed = TRUE)
    a <- data.frame(arm = vapply(parts, `[`, "", 1L), id = vapply(parts, `[`, "", 2L),
                    stringsAsFactors = FALSE)
    if (!is.null(cage)) a$cage <- vapply(parts, function(p) if (length(p) >= 3L) p[3L] else NA_character_, "")
    a
  }
  list(animals = animals,
       treatments = if (length(ex_tx)) ex_tx,
       cages = if (length(ex_cg)) ex_cg,
       treatment_column = s$col_treatment %||% "",
       id_column = s$col_id %||% "",
       cage_column = cage,
       dose_column = col(s$col_dose))
}

#' Sorted distinct values of the treatment column
#' @keywords internal
treatment_values <- function(data, tx_col) {
  if (!has_col(data, tx_col)) return(character(0))
  sort(unique(as.character(stats::na.omit(data[[tx_col]]))))
}

#' Each row's arm: the treatment, or "Treatment Dose" when a dose column is mapped
#'
#' Tumour Growth and Survival analyse each treatment-dose combination as its
#' own arm, labelled this way.
#' @keywords internal
arm_labels <- function(data, tx_col, dose_col = NULL) {
  if (has_col(data, dose_col)) {
    paste(as.character(data[[tx_col]]), as.character(data[[dose_col]]), sep = " ")
  } else as.character(data[[tx_col]])
}

#' The arms an analysis will see, after the subset filter
#' @param rv Shared reactiveValues.
#' @return Sorted character vector (empty when the treatment column is unmapped).
#' @keywords internal
analysis_group_labels <- function(rv) {
  rv$filter_applied_at  # re-evaluate after Apply / Clear
  df <- apply_subset_filter(rv$uploaded_data, rv)
  tx <- rv$settings$col_treatment
  if (is.null(df) || !has_col(df, tx)) return(character(0))
  df <- df[!is.na(df[[tx]]), , drop = FALSE]
  sort(unique(arm_labels(df, tx, rv$settings$col_dose)))
}

#' Reduce longitudinal data to one row per animal for survival analysis
#'
#' CODE_REVIEW.md R20.D2. Animals were keyed by ID (plus cage), so ear tags
#' reused across arms merged animals from different arms into one record with
#' the latest day, the largest event flag and an arbitrary treatment. They are
#' now keyed by arm + ID (+ cage). The event time was each animal's last day,
#' even when measurements continued after a flagged event; it is now the first
#' flagged day. When every row of an animal is flagged, the flag is an
#' animal-level status rather than an event row, and the last day is used.
#'
#' @param data Longitudinal data.
#' @param time_col,tx_col,id_col Column names (tx_col is the arm column).
#' @param cage_col Optional cage column.
#' @param event_col Optional 0/1 event column. Without it, each animal's last
#'   observation is returned.
#' @return One row per animal, with the event column set to the animal's event
#'   status and an `.animal_key` column.
#' @keywords internal
survival_animal_rows <- function(data, time_col, tx_col, id_col,
                                 cage_col = NULL, event_col = NULL) {
  key <- animal_key(data[[tx_col]], data[[id_col]],
                    if (has_col(data, cage_col)) data[[cage_col]])
  tt  <- suppressWarnings(as.numeric(data[[time_col]]))
  ord <- order(key, tt, na.last = TRUE)
  d <- data[ord, , drop = FALSE]; key <- key[ord]; tt <- tt[ord]
  ev <- if (has_col(d, event_col)) suppressWarnings(as.numeric(d[[event_col]])) else NULL
  rows <- split(seq_len(nrow(d)), factor(key, levels = unique(key)))
  pick <- vapply(rows, function(i) {
    if (!is.null(ev)) {
      e <- ev[i]
      flagged <- which(!is.na(e) & e == 1)
      if (length(flagged) && length(flagged) < sum(!is.na(e))) return(i[flagged[1L]])
    }
    ok <- which(!is.na(tt[i]))
    if (length(ok)) i[ok[length(ok)]] else i[length(i)]
  }, integer(1L))
  out <- d[pick, , drop = FALSE]
  if (!is.null(ev)) {
    out[[event_col]] <- vapply(rows, function(i) {
      e <- ev[i]
      if (any(!is.na(e))) max(e, na.rm = TRUE) else NA_real_
    }, numeric(1L))
  }
  out$.animal_key <- names(rows)
  rownames(out) <- NULL
  out
}

#' Event times from a tumour-volume threshold
#'
#' The first day each animal's volume reaches `threshold` is its event time
#' (Event = 1); an animal that never reaches it is censored at its last
#' observation (Event = 0). Keyed like [survival_animal_rows()] (R20.D2).
#'
#' @inheritParams survival_animal_rows
#' @param vol_col Volume column.
#' @param threshold Volume limit, in the data's units.
#' @return Data frame with `.animal_key`, `Time` and `Event`.
#' @keywords internal
volume_threshold_events <- function(data, time_col, tx_col, id_col,
                                    cage_col = NULL, vol_col, threshold) {
  key  <- animal_key(data[[tx_col]], data[[id_col]],
                     if (has_col(data, cage_col)) data[[cage_col]])
  tt   <- suppressWarnings(as.numeric(data[[time_col]]))
  vol  <- suppressWarnings(as.numeric(data[[vol_col]]))
  rows <- split(seq_len(nrow(data)), factor(key, levels = unique(key)))
  res  <- lapply(rows, function(i) {
    i   <- i[order(tt[i])]
    hit <- which(!is.na(vol[i]) & vol[i] >= threshold)
    if (length(hit)) c(tt[i][hit[1L]], 1)
    else c(suppressWarnings(max(tt[i], na.rm = TRUE)), 0)
  })
  data.frame(.animal_key = names(rows),
             Time  = vapply(res, `[`, numeric(1L), 1L),
             Event = as.integer(vapply(res, `[`, numeric(1L), 2L)),
             stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
}

# ---- Dose series (R20.D47) -----------------------------------------------------

#' A suggested dose series: one agent's arms
#'
#' Package R20.16 analyses one agent at several doses. The suggestion is the
#' largest set of single-agent arms (no "+" in the name) that share the first
#' word of their name, have one dose each, and differ in dose. A single arm is
#' its own series (one treatment name at every dose). Otherwise nothing is
#' suggested, and the user chooses.
#'
#' @param data Data with the treatment and dose columns.
#' @param tx_col,dose_col Column names.
#' @param arms The arms to choose from (the control excluded).
#' @return Character vector of arms.
#' @keywords internal
guess_dose_series <- function(data, tx_col, dose_col, arms) {
  if (length(arms) <= 1L) return(arms)
  single <- arms[!grepl("+", arms, fixed = TRUE)]
  if (!length(single)) return(character(0))
  stem <- sub("\\s.*$", "", single)
  n    <- table(stem)
  best <- names(n)[n == max(n)][1L]
  pick <- single[stem == best]
  dose <- suppressWarnings(as.numeric(data[[dose_col]]))
  one_dose <- vapply(pick, function(a) {
    v <- unique(dose[!is.na(data[[tx_col]]) & data[[tx_col]] == a & is.finite(dose)])
    if (length(v) == 1L) v else NA_real_
  }, numeric(1L))
  if (length(pick) >= 2L && all(is.finite(one_dose)) && !anyDuplicated(one_dose)) pick
  else character(0)
}

# ---- Removal reasons (R20-D.H) ------------------------------------------------

#' Categories for an optional removal-reason column
#'
#' R20-D.H: when the column is mapped, each of its values is assigned one of
#' these on Data Upload. Survival can then count deaths from other causes as
#' events (R20.D13), and the weight-loss threshold separates weight-loss
#' removals, planned ends and competing removals (package R20.23).
#' @keywords internal
removal_reason_categories <- c(
  "Tumour volume limit"           = "limit",
  "Death or moribund"             = "death",
  "Weight loss or body condition" = "weight",
  "Planned end or data cut"       = "planned",
  "Other"                         = "other")

#' A first guess at a removal reason's category, from keywords
#'
#' Only a starting point: every value's category is shown and can be changed
#' on Data Upload.
#' @param x Character vector of reason values.
#' @return Character vector of category codes.
#' @keywords internal
guess_removal_category <- function(x) {
  v <- tolower(as.character(x))
  out <- rep("other", length(v))
  rules <- list(
    weight  = "weight|body condition|\\bbcs\\b|emaciat|cachex",
    death   = "dead|death|died|moribund|natural",
    planned = "end of study|study end|terminal|scheduled|planned|data cut|cut-?off|\\beos\\b",
    limit   = "tumou?r|volume|size|burden")
  for (cat in rev(names(rules))) out[grepl(rules[[cat]], v)] <- cat
  out
}

#' The distinct non-empty values of a removal-reason column
#' @keywords internal
removal_reason_values <- function(data, col) {
  if (!has_col(data, col)) return(character(0))
  v <- trimws(as.character(data[[col]]))
  sort(unique(v[!is.na(v) & nzchar(v)]))
}

#' The mapped removal reasons, or NULL when no reason column is mapped
#'
#' @param rv The shared reactiveValues.
#' @return `list(column, map)`, where `map` is a named character vector from
#'   each reason value to its category code.
#' @keywords internal
removal_reason_spec <- function(rv) {
  col <- rv$settings$col_removal_reason %||% ""
  map <- unlist(rv$settings$removal_reason_map)
  if (!nzchar(col) || !length(map)) return(NULL)
  list(column = col, map = map)
}

#' The reason values mapped to some categories
#' @param spec Output of [removal_reason_spec()].
#' @param categories Category codes.
#' @keywords internal
reasons_in <- function(spec, categories) {
  if (is.null(spec)) return(character(0))
  names(spec$map)[spec$map %in% categories]
}

#' Each animal's removal reason: its last non-empty value
#'
#' Keyed like [survival_animal_rows()].
#' @inheritParams survival_animal_rows
#' @param reason_col The removal-reason column.
#' @return Named character vector, animal key to reason (NA when none).
#' @keywords internal
animal_removal_reasons <- function(data, tx_col, id_col, cage_col = NULL,
                                   reason_col, time_col = NULL) {
  key <- animal_key(data[[tx_col]], data[[id_col]],
                    if (has_col(data, cage_col)) data[[cage_col]])
  r   <- trimws(as.character(data[[reason_col]]))
  ord <- if (has_col(data, time_col))
    order(key, suppressWarnings(as.numeric(data[[time_col]])), na.last = TRUE)
  else seq_along(key)
  key <- key[ord]; r <- r[ord]
  ok  <- !is.na(r) & nzchar(r)
  last <- tapply(r[ok], factor(key[ok], levels = unique(key)),
                 function(x) x[length(x)])
  out <- stats::setNames(rep(NA_character_, length(unique(key))), unique(key))
  out[names(last)[!is.na(last)]] <- as.character(last[!is.na(last)])
  out
}

# ---- Date Preprocessing (consolidated — Issue 2.3) ---------------------------

# Regular expressions for the supported layouts, tried in this order. Month
# names are matched with .month_number() rather than strptime's %b/%B, which
# follow the server's locale.
.date_patterns <- c(
  ymd     = "^(\\d{4})[-/.](\\d{1,2})[-/.](\\d{1,2})(?:[ T].*)?$",
  numeric = "^(\\d{1,2})([-/.])(\\d{1,2})\\2(\\d{4})$",
  d_mon_y = "^(\\d{1,2})[- ]([A-Za-z]{3,9})\\.?[- ,]+(\\d{4})$",
  mon_d_y = "^([A-Za-z]{3,9})\\.?\\s+(\\d{1,2}),?\\s+(\\d{4})$",
  d_mon   = "^(\\d{1,2})[- ]([A-Za-z]{3,9})\\.?$",
  mon_d   = "^([A-Za-z]{3,9})\\.?\\s+(\\d{1,2})$"
)

.date_layout_labels <- c(
  ymd = "YYYY-MM-DD", dmy = "DD/MM/YYYY", mdy = "MM/DD/YYYY",
  d_mon_y = "DD-Mon-YYYY", mon_d_y = "Month DD, YYYY",
  d_mon = "DD-Mon (no year)", mon_d = "Month DD (no year)"
)

.month_number <- function(x) {
  match(tolower(substr(x, 1L, 3L)), tolower(month.abb))
}

.ymd_date <- function(y, m, d) {
  out <- rep(as.Date(NA), length(y))
  ok <- !is.na(y) & !is.na(m) & !is.na(d)
  out[ok] <- as.Date(sprintf("%04d-%02d-%02d", as.integer(y[ok]),
                             as.integer(m[ok]), as.integer(d[ok])),
                     format = "%Y-%m-%d")
  out
}

# Split x by layout. Returns the layout name for each value (NA when none
# matches) and the regmatches() pieces.
.match_date_layouts <- function(x) {
  layout <- rep(NA_character_, length(x))
  parts  <- vector("list", length(x))
  for (nm in names(.date_patterns)) {
    todo <- which(is.na(layout) & !is.na(x))
    if (!length(todo)) break
    m <- regmatches(x[todo], regexec(.date_patterns[[nm]], x[todo], perl = TRUE))
    hit <- lengths(m) > 0L
    layout[todo[hit]] <- nm
    parts[todo[hit]]  <- m[hit]
  }
  list(layout = layout, parts = parts)
}

# Resolve the order of d/d/yyyy values. `a` and `b` are the first and second
# numbers. Returns list(order, reason), or stops when the column cannot be
# read in one order or the data contradict the declared order.
.resolve_numeric_order <- function(a, b, examples, date_order, what) {
  first_big  <- which(a > 12L)
  second_big <- which(b > 12L)
  if (length(first_big) && length(second_big)) {
    stop("The ", what, " contains '", examples[first_big[1L]], "', which can only be ",
         "DD/MM/YYYY, and '", examples[second_big[1L]], "', which can only be ",
         "MM/DD/YYYY. Correct the values so the column uses one order.",
         call. = FALSE)
  }
  forced <- if (length(first_big)) "dmy" else if (length(second_big)) "mdy" else NA_character_
  if (!is.na(forced)) {
    ex <- examples[c(first_big, second_big)[1L]]
    if (date_order != "auto" && date_order != forced) {
      stop("The date order is set to ", .date_layout_labels[[date_order]],
           ", but ", what, " contains '", ex, "', which would need month ",
           if (forced == "dmy") a[first_big[1L]] else b[second_big[1L]],
           " in that order. Set the date order to ",
           .date_layout_labels[[forced]], " or 'Detect from the data'.",
           call. = FALSE)
    }
    return(list(order = forced, reason = sprintf(
      "'%s' can only be %s", ex, .date_layout_labels[[forced]])))
  }
  if (date_order == "auto") {
    stop("The dates in ", what, " could be DD/MM/YYYY or MM/DD/YYYY: no value ",
         "has a day above 12 (for example '", examples[1L], "'). Choose the ",
         "date order under the date column on Data Upload.", call. = FALSE)
  }
  list(order = date_order, reason = "the date order set on Data Upload")
}

#' Resolve a date column's layout and read every value
#'
#' CODE_REVIEW.md R20.D1. The layout used to be chosen from the first row
#' alone, and any `d/d/yyyy` value was assumed to be month-first, so a
#' DD/MM/YYYY export was read as MM/DD/YYYY without a warning: days 1-12 became
#' wrong dates and days 13-31 became NA. The whole column is now read. A first
#' number above 12 means day-first and a second number above 12 means
#' month-first; a column with neither is ambiguous, so the declared order is
#' used, or the user is asked to declare one. Year-less dates ("24-Mar",
#' "March 24") take the year that puts them nearest the reference date, so a
#' study that crosses New Year no longer gets negative days.
#'
#' Supported layouts: YYYY-MM-DD (also with / or . and a trailing time);
#' DD/MM/YYYY and MM/DD/YYYY (with /, . or -); 24-Mar-2025; March 24, 2025;
#' and the year-less 24-Mar and March 24. Month names are English.
#'
#' @param data Data frame containing the date column.
#' @param date_column Name of the date column.
#' @param reference_date The study's Day 0, preferably as YYYY-MM-DD. A
#'   d/d/yyyy reference date is read in the column's order.
#' @param date_order "auto", "dmy" or "mdy": how to read d/d/yyyy values
#'   when the data do not decide it.
#' @return A list: `dates` (Date, NA where a value could not be read),
#'   `reference` (Date), `layout` (for example "DD/MM/YYYY"), `reason` (how a
#'   d/d/yyyy order was decided, or NULL) and `notes` (warnings for the user).
#' @keywords internal
detect_date_format <- function(data, date_column, reference_date,
                               date_order = c("auto", "dmy", "mdy")) {
  date_order <- match.arg(date_order)
  # R20.D4: an unusable reference date is reported, never replaced.
  if (is.null(reference_date) || !nzchar(trimws(reference_date))) {
    stop("No reference date is set. Enter the study's reference date (Day 0) ",
         "as YYYY-MM-DD, for example 2025-03-24.", call. = FALSE)
  }
  if (is.null(date_column) || !date_column %in% names(data)) {
    stop("Date column '", date_column %||% "", "' is not in the data.",
         call. = FALSE)
  }
  what <- sprintf("column '%s'", date_column)
  raw  <- trimws(as.character(data[[date_column]]))
  raw[!is.na(raw) & !nzchar(raw)] <- NA_character_
  m <- .match_date_layouts(raw)
  if (all(is.na(m$layout))) {
    ex <- utils::head(raw[!is.na(raw)], 1L)
    stop("None of the values in ", what, " is a recognised date",
         if (length(ex)) paste0(" (for example '", ex, "')"), ". Use ",
         "YYYY-MM-DD, DD/MM/YYYY or MM/DD/YYYY, or map a Day column instead.",
         call. = FALSE)
  }
  piece <- function(idx, k) vapply(m$parts[idx], `[`, character(1L), k)

  # d/d/yyyy: decide the order from the whole column.
  num_idx <- which(m$layout == "numeric")
  order <- NULL; reason <- NULL
  if (length(num_idx)) {
    res <- .resolve_numeric_order(as.integer(piece(num_idx, 2L)),
                                  as.integer(piece(num_idx, 4L)),
                                  raw[num_idx], date_order, what)
    order <- res$order; reason <- res$reason
  }

  # The reference date: YYYY-MM-DD, a named month with its year, or d/d/yyyy
  # in the column's order (or the declared one).
  ref_raw <- trimws(reference_date)
  rm <- .match_date_layouts(ref_raw)
  bad_ref <- function() stop(
    "Reference date '", ref_raw, "' is not a recognised date. Enter it as ",
    "YYYY-MM-DD, for example 2025-03-24.", call. = FALSE)
  ref_layout <- if (is.na(rm$layout)) "none" else rm$layout
  ref <- switch(ref_layout,
    ymd = .ymd_date(rm$parts[[1L]][2L], rm$parts[[1L]][3L], rm$parts[[1L]][4L]),
    d_mon_y = .ymd_date(rm$parts[[1L]][4L], .month_number(rm$parts[[1L]][3L]),
                        rm$parts[[1L]][2L]),
    mon_d_y = .ymd_date(rm$parts[[1L]][4L], .month_number(rm$parts[[1L]][2L]),
                        rm$parts[[1L]][3L]),
    numeric = {
      a <- as.integer(rm$parts[[1L]][2L]); b <- as.integer(rm$parts[[1L]][4L])
      ro <- if (a > 12L) "dmy"
            else if (b > 12L) "mdy"
            else if (!is.null(order)) order
            else if (date_order != "auto") date_order
            else NULL
      if (is.null(ro)) {
        stop("Reference date '", ref_raw, "' could be DD/MM/YYYY or ",
             "MM/DD/YYYY. Enter it as YYYY-MM-DD, for example 2025-03-24.",
             call. = FALSE)
      }
      if (ro == "dmy") .ymd_date(rm$parts[[1L]][5L], b, a)
      else .ymd_date(rm$parts[[1L]][5L], a, b)
    },
    as.Date(NA))
  if (is.na(ref)) bad_ref()
  ref_year <- as.integer(format(ref, "%Y"))

  dates <- rep(as.Date(NA), length(raw))
  idx <- which(m$layout == "ymd")
  if (length(idx)) dates[idx] <- .ymd_date(piece(idx, 2L), piece(idx, 3L), piece(idx, 4L))
  if (length(num_idx)) {
    a <- piece(num_idx, 2L); b <- piece(num_idx, 4L); y <- piece(num_idx, 5L)
    dates[num_idx] <- if (order == "dmy") .ymd_date(y, b, a) else .ymd_date(y, a, b)
  }
  idx <- which(m$layout == "d_mon_y")
  if (length(idx)) dates[idx] <- .ymd_date(piece(idx, 4L), .month_number(piece(idx, 3L)), piece(idx, 2L))
  idx <- which(m$layout == "mon_d_y")
  if (length(idx)) dates[idx] <- .ymd_date(piece(idx, 4L), .month_number(piece(idx, 2L)), piece(idx, 3L))

  # Year-less values take the year that puts them nearest the reference date.
  notes <- character(0)
  yl <- which(m$layout %in% c("d_mon", "mon_d"))
  if (length(yl)) {
    is_dm <- m$layout[yl] == "d_mon"
    mon <- ifelse(is_dm, .month_number(piece(yl, 3L)), .month_number(piece(yl, 2L)))
    day <- ifelse(is_dm, piece(yl, 2L), piece(yl, 3L))
    cand <- lapply(ref_year + (-1L):1L, function(yy) .ymd_date(rep(yy, length(yl)), mon, day))
    gap  <- vapply(cand, function(d) abs(as.numeric(d - ref)), numeric(length(yl)))
    gap  <- matrix(gap, nrow = length(yl))
    best <- apply(gap, 1L, function(g) if (all(is.na(g))) NA_integer_ else which.min(g))
    for (k in seq_along(yl)) {
      if (!is.na(best[k])) dates[yl[k]] <- cand[[best[k]]][k]
    }
    other_year <- sum(!is.na(best) & best != 2L)
    if (other_year) {
      notes <- c(notes, sprintf(paste0(
        "Dates without a year were given the year that puts them nearest the ",
        "reference date (%s); %d fall in a different year from it, so the ",
        "study crosses New Year."), format(ref), other_year))
    }
  }

  kinds <- unique(stats::na.omit(ifelse(m$layout == "numeric", order, m$layout)))
  layout <- paste(unname(.date_layout_labels[kinds]), collapse = " + ")
  if (length(kinds) > 1L) {
    notes <- c(notes, sprintf(
      "The %s mixes date layouts (%s); each value was read in its own layout.",
      what, layout))
  }
  unread <- which(!is.na(raw) & is.na(dates))
  blank  <- sum(is.na(raw))
  if (length(unread) || blank) {
    notes <- c(notes, paste0(
      length(unread) + blank, " row(s) have no usable date",
      if (length(unread)) sprintf(" (for example '%s')", raw[unread[1L]]),
      "; they have no Day and are left out of analyses."))
  }
  list(dates = dates, reference = ref, layout = layout, reason = reason,
       notes = notes)
}

#' Add a study Day column computed from a date column
#'
#' Reads the dates with [detect_date_format()] and passes them to
#' `mouseExperiment::calculate_dates()` as unambiguous YYYY-MM-DD values, so
#' the backend defines Day and the dashboard only decides what the dates are.
#' Days before the reference date are reported (R20.D1).
#'
#' @inheritParams detect_date_format
#' The exported [study_days()] returns the data and raises the notes as a
#' warning; the dashboard shows them as a notification instead.
#'
#' @return A list: `data` (with a `Day` column), `layout`, `reason` and
#'   `notes`.
#' @keywords internal
study_days_info <- function(data, date_column, reference_date, date_order = "auto") {
  info <- detect_date_format(data, date_column, reference_date, date_order)
  tmp <- data
  tmp[[".me_date_iso"]] <- format(info$dates, "%Y-%m-%d")
  out <- withCallingHandlers(
    mouseExperiment::calculate_dates(
      tmp, start_date = format(info$reference, "%Y-%m-%d"),
      date_column = ".me_date_iso", date_format = "%Y-%m-%d"),
    # Unreadable values are already described in info$notes.
    warning = function(w) {
      if (grepl("could not be parsed", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  out[[".me_date_iso"]] <- NULL
  notes <- info$notes
  neg <- which(!is.na(out$Day) & out$Day < 0)
  if (length(neg)) {
    notes <- c(notes, sprintf(paste0(
      "%d measurement(s) fall before the reference date %s (earliest Day %s). ",
      "Check that the reference date is the study's Day 0."),
      length(neg), format(info$reference), format(min(out$Day[neg]))))
  }
  list(data = out, layout = info$layout, reason = info$reason, notes = notes)
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
#' @return A list with \code{data}, \code{time_col}, \code{vol_col}, \code{tx_col}
#'   and \code{steps}, the recorded preparation (see [prepare_upload()]).
#' @keywords internal
ensure_preprocessed <- function(rv) {
  # Reuse canonical result if available. CODE_REVIEW.md Round 2 F.5:
  # always re-apply the Subset Filter to the cached frame so changes to
  # the filter take effect on the next analysis run without needing a
  # full re-upload + re-preprocess. (Applying or clearing the filter also
  # clears the cache, so the recorded steps hold the current filter.)
  if (!is.null(rv$processed_data)) {
    return(list(
      data     = apply_subset_filter(rv$processed_data, rv),
      time_col = rv$time_col,
      vol_col  = rv$vol_col,
      tx_col   = rv$tx_col,
      steps    = rv$processed_steps
    ))
  }

  pu <- prepare_upload(rv)
  settings <- rv$settings
  tx_col <- settings$col_treatment

  # Resolve optional columns (id, cage, dose) so other modules don't
  # depend on TG analysis having been run first.
  rv$id_col   <- if (!is.null(settings$col_id)   && nzchar(settings$col_id))   settings$col_id   else NULL
  rv$cage_col <- if (!is.null(settings$col_cage)  && nzchar(settings$col_cage)) settings$col_cage  else NULL
  rv$dose_col   <- if (!is.null(settings$col_dose)   && nzchar(settings$col_dose))   settings$col_dose   else NULL
  rv$weight_col <- if (!is.null(settings$col_weight) && nzchar(settings$col_weight)) settings$col_weight else NULL
  rv$sex_col    <- if (!is.null(settings$col_sex)    && nzchar(settings$col_sex))    settings$col_sex    else NULL

  # Cache for future callers
  rv$processed_data  <- pu$data
  rv$processed_steps <- pu$steps
  rv$time_col        <- pu$time_col
  rv$vol_col         <- pu$vol_col
  rv$tx_col          <- tx_col

  list(data = pu$data, time_col = pu$time_col, vol_col = pu$vol_col, tx_col = tx_col,
       steps = pu$steps,
       date_fmt_warning = if (length(pu$notes)) paste(pu$notes, collapse = " "))
}

#' Prepare the upload, recording each step: exclusions, study days, volume
#'
#' Tumor Growth and [ensure_preprocessed()] both prepare an upload this way.
#' Each step is an exported function (R20-AB) and is recorded with
#' [record_call()] (R20-AC), so the R script can replay it from the uploaded
#' file. The exclusions come first, so excluded rows play no part in deciding
#' how dates are read.
#' @param rv Shared reactiveValues.
#' @return list(data, time_col, vol_col, steps, notes): `notes` are what
#'   reading the dates found (see [study_days_info()]).
#' @keywords internal
prepare_upload <- function(rv) {
  s <- rv$settings
  data <- rv$uploaded_data
  steps <- list()
  notes <- character(0)

  ex <- subset_filter_args(data, rv)
  if (!is.null(ex) && !is.null(data) && nrow(data) > 0L) {
    args <- c(list(df = data), ex)
    data <- do.call(apply_exclusions, args)
    steps <- c(steps, list(record_call("apply_exclusions", args, keep_data = FALSE)))
  }

  # Date -> Day. R20.D1: the whole column decides the day/month order, and
  # an ambiguous or contradictory layout stops instead of being guessed.
  if (identical(s$time_method, "from_date")) {
    args <- list(df = data, date_column = s$col_date_for_day,
                 reference_date = s$reference_date, date_order = s$date_order %||% "auto")
    info <- study_days_info(data, args$date_column, args$reference_date, args$date_order)
    data <- info$data
    notes <- info$notes
    steps <- c(steps, list(record_call("study_days", args, keep_data = FALSE)))
    time_col <- "Day"
  } else {
    time_col <- s$col_day
  }

  # Volume from calipers. Without both caliper columns nothing is calculated
  # (Tumor Growth reports that before it gets here).
  if (identical(s$volume_method, "calculate")) {
    len <- s$col_length %||% "Length"
    wid <- s$col_width  %||% "Width"
    if (has_col(data, len) && has_col(data, wid)) {
      args <- list(df = data, length_column = len, width_column = wid,
                   height_column = if (has_col(data, s$col_height)) s$col_height,
                   formula = s$volume_formula %||% "ellipsoid")   # calculate_volume()'s default
      data <- do.call(mouseExperiment::calculate_volume, args)
      steps <- c(steps, list(record_call("calculate_volume", args, keep_data = FALSE)))
    }
    vol_col <- "Volume"
  } else {
    vol_col <- s$col_volume %||% "Volume"
  }
  list(data = data, time_col = time_col, vol_col = vol_col, steps = steps, notes = notes)
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
      data <- read_study_csv(file_path)   # R20.D23
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

# read_study_csv() moved to R/prepare_study.R and is exported (R20-AB).

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
    # (A survival_time role was detected here and never used: time to event is
    # each animal's last measurement day, or the day it crosses the volume
    # limit. Removed in 26.10.009 with its Help row.)
    survival_event = detect_column(col_names, col_names_lower,
                                    c("survival_event", "event", "survival_censor", "censor")),
    weight = detect_column(col_names, col_names_lower,
                           c("weight", "body_weight", "body weight", "mouse_weight",
                             "mouse weight", "mass", "body_mass", "body mass", "bw")),
    sex    = detect_column(col_names, col_names_lower,
                           c("sex", "gender", "m/f")),
    necrotic = detect_column(col_names, col_names_lower,
                             c("necrotic", "necrosis", "necr",
                               "is_necrotic", "necrotic_flag")),
    # R20-D.H: why each animal left the study.
    removal_reason = detect_column(col_names, col_names_lower,
                                   c("removal_reason", "removal reason", "reason",
                                     "cause_of_death", "cause of death", "fate"))
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

#' The rows of one animal in the Tumour Growth drill-down
#'
#' CODE_REVIEW.md R20.D57: one animal is treatment + ID + cage. Matching on ID
#' and treatment alone plotted cage-mates that share an ear tag as one mouse
#' (29 rows from two cages on the Combo demo).
#'
#' @param data The analysed data.
#' @param rate_row One row of the backend's growth-rate table (columns
#'   `Treatment`, `ID` and, when a cage column is mapped, `Cage`).
#' @param id_col,tx_col,cage_col Column names in `data`.
#' @return The animal's rows of `data`.
#' @keywords internal
drilldown_rows <- function(data, rate_row, id_col, tx_col, cage_col = NULL) {
  keep <- as.character(data[[id_col]]) == as.character(rate_row$ID[1L]) &
          as.character(data[[tx_col]]) == as.character(rate_row$Treatment[1L])
  if ("Cage" %in% names(rate_row) && has_col(data, cage_col)) {
    keep <- keep & as.character(data[[cage_col]]) == as.character(rate_row$Cage[1L])
  }
  data[which(keep), , drop = FALSE]
}
