#' Check a study's data for likely errors
#'
#' Looks through longitudinal tumour data for entries that are probably
#' mistakes, entries worth a second look, and features of the design that
#' change how the results read. Nothing is changed or removed: the result
#' lists what was found, with each row's line in the file when `row_column`
#' gives it.
#'
#' **Problems** (probably errors):
#' * *Missing ID, arm or day*: a row that cannot be placed.
#' * *Duplicate measurement*: one animal measured twice on one day, with the
#'   same volume (counted twice) or different ones (which is right?).
#' * *Negative volume* or *non-positive body weight*.
#'
#' **To check:**
#' * *No volume*: a row with an animal and a day but no volume.
#' * *Large jump*: between consecutive measurements the volume rose more than
#'   `jump_fold`-fold or fell below `1 / jump_fold` of its value, both volumes
#'   at least the study's 10th percentile (small tumours are noisy).
#' * *Width above length*: caliper width greater than length; the volume
#'   formulas swap them, so check the columns are not reversed.
#' * *Missed measurement*: an animal measured before and after a day has no row
#'   on it, though at least three quarters of its arm's animals on study were
#'   measured then. It has no `Row`.
#' * *Weight drop*: body weight fell more than `weight_drop` (a fraction)
#'   between consecutive weighings.
#' * *Measured after its event*: rows after the first day an animal's event
#'   column is 1.
#' * *Single measurement*: an animal measured once, which adds nothing to a
#'   growth model.
#'
#' **Notes** (design, not errors):
#' * *Staggered start*: animals whose first day differs from their arm's
#'   usual first day.
#' * *Small arm*: fewer than 3 animals.
#' * *IDs reused across arms*: the analyses keep them apart by arm.
#' * *Cage labels reused across arms*: see [label_cages()].
#'
#' @param df Study data, one row per measurement, with study days and volumes
#'   (see [study_days()] and [calculate_volume()]).
#' @param time_column,volume_column,treatment_column,id_column Column names.
#' @param cage_column,dose_column,weight_column,event_column Optional column
#'   names, or NULL. With a dose, an arm is "Treatment Dose".
#' @param length_column,width_column Optional caliper columns, or NULL. Give
#'   them as measured: [calculate_volume()] swaps a width above the length.
#' @param row_column Optional column giving each row's line in the file;
#'   without it, rows are numbered as in `df`.
#' @param jump_fold The fold change between consecutive measurements counted as
#'   a large jump.
#' @param weight_drop The fractional drop in body weight between consecutive
#'   weighings counted as a weight drop.
#' @return A list:
#'   * `summary`: one row per check that found something, with `Severity`
#'     ("Problem", "Check" or "Note"), `Check`, `Animals`, `Rows` and `What`
#'     (what it means); no rows when nothing was found;
#'   * `flagged`: one row per flagged row, with `Row`, `Arm`, `ID`, `Cage`,
#'     `Day`, `Severity`, `Check` and `Detail`.
#' @export
#' @examples
#' d <- master_synthetic_data
#' d$Volume[5] <- -d$Volume[5]          # a sign error
#' d <- rbind(d, d[10, ])               # a row entered twice
#' qc <- check_study_data(d, cage_column = "Cage", weight_column = "Weight")
#' qc$summary
check_study_data <- function(df, time_column = "Day", volume_column = "Volume",
                             treatment_column = "Treatment", id_column = "ID",
                             cage_column = NULL, dose_column = NULL,
                             weight_column = NULL, event_column = NULL,
                             length_column = NULL, width_column = NULL,
                             row_column = NULL, jump_fold = 4, weight_drop = 0.15) {
  for (col in c(time_column, volume_column, treatment_column, id_column)) {
    if (!has_col(df, col)) stop("Column '", col, "' is not in the data.", call. = FALSE)
  }
  num <- function(col) if (has_col(df, col)) suppressWarnings(as.numeric(df[[col]]))
  d <- data.frame(
    Row = if (has_col(df, row_column)) df[[row_column]] else seq_len(nrow(df)),
    Arm = arm_labels(df, treatment_column, if (has_col(df, dose_column)) dose_column),
    ID = as.character(df[[id_column]]),
    Cage = if (has_col(df, cage_column)) as.character(df[[cage_column]]) else NA_character_,
    Day = num(time_column), Volume = num(volume_column), stringsAsFactors = FALSE)
  d$Arm[is.na(df[[treatment_column]])] <- NA_character_
  if (has_col(df, weight_column)) d$Weight <- num(weight_column)
  if (has_col(df, event_column)) d$Event <- num(event_column)
  if (has_col(df, length_column) && has_col(df, width_column)) {
    d$Length <- num(length_column)
    d$Width <- num(width_column)
  }
  d$Key <- if (has_col(df, cage_column)) make_mouse_key(d$Arm, d$ID, d$Cage) else make_mouse_key(d$Arm, d$ID)

  found <- list()
  flag <- function(rows, severity, check, detail) {
    if (!length(rows) || !any(rows)) return(invisible())
    i <- which(rows)
    found[[length(found) + 1L]] <<- data.frame(
      d[i, c("Row", "Arm", "ID", "Cage", "Day")], Severity = severity, Check = check,
      Detail = rep_len(detail[if (length(detail) > 1L) i else 1L], length(i)),
      stringsAsFactors = FALSE)
  }

  # ---- Problems -------------------------------------------------------------
  placed <- !is.na(d$Arm) & !is.na(d$ID) & nzchar(d$ID) & is.finite(d$Day)
  flag(!placed, "Problem", "Missing ID, arm or day", "The row cannot be placed in an animal's record.")
  p <- d[placed, , drop = FALSE]
  key_day <- paste(p$Key, p$Day)
  dup <- duplicated(key_day) | duplicated(key_day, fromLast = TRUE)
  same <- dup & stats::ave(p$Volume, key_day, FUN = function(v) length(unique(v))) == 1
  flag(placed & seq_len(nrow(d)) %in% which(placed)[dup & same], "Problem",
       "Duplicate measurement", "The same volume twice on one day: it is counted twice.")
  flag(placed & seq_len(nrow(d)) %in% which(placed)[dup & !same], "Problem",
       "Duplicate measurement", "Different volumes on one day: which is right?")
  flag(is.finite(d$Volume) & d$Volume < 0, "Problem", "Negative volume", "A volume below zero.")
  if (!is.null(d$Weight)) {
    flag(is.finite(d$Weight) & d$Weight <= 0, "Problem", "Non-positive body weight",
         "A body weight of zero or below.")
  }

  # ---- To check -------------------------------------------------------------
  flag(placed & !is.finite(d$Volume), "Check", "No volume",
       "The animal and day are given, but no volume.")
  if (!is.null(d$Length)) {
    flag(is.finite(d$Length) & is.finite(d$Width) & d$Width > d$Length, "Check",
         "Width above length", "The formulas swap them; check the columns are not reversed.")
  }
  ok <- placed & is.finite(d$Volume) & !duplicated(paste(d$Key, d$Day))
  o <- d[ok, , drop = FALSE]
  o$.i <- which(ok)
  o <- o[order(o$Key, o$Day), , drop = FALSE]
  same_animal <- c(FALSE, o$Key[-1L] == o$Key[-nrow(o)])
  prev_v <- c(NA, o$Volume[-nrow(o)])
  floor_v <- if (any(o$Volume > 0)) stats::quantile(o$Volume[o$Volume > 0], 0.1, names = FALSE) else Inf
  big <- same_animal & is.finite(prev_v) & prev_v >= floor_v & o$Volume >= floor_v &
    (o$Volume > jump_fold * prev_v | o$Volume < prev_v / jump_fold)
  flag(seq_len(nrow(d)) %in% o$.i[big], "Check", "Large jump",
       sprintf("Volume %s from %s to %s since the previous measurement.",
               ifelse(o$Volume > prev_v, "rose", "fell"), format(signif(prev_v, 3)),
               format(signif(o$Volume, 3)))[match(seq_len(nrow(d)), o$.i)])
  if (!is.null(d$Weight)) {
    prev_w <- c(NA, o$Weight[-nrow(o)])
    drop <- same_animal & is.finite(prev_w) & is.finite(o$Weight) & prev_w > 0 &
      o$Weight < (1 - weight_drop) * prev_w
    flag(seq_len(nrow(d)) %in% o$.i[drop], "Check", "Weight drop",
         sprintf("Body weight fell %.0f %% since the previous weighing.",
                 100 * (1 - o$Weight / prev_w))[match(seq_len(nrow(d)), o$.i)])
  }
  if (!is.null(d$Event)) {
    first_event <- stats::ave(ifelse(o$Event %in% 1, o$Day, Inf), o$Key, FUN = min)
    flag(seq_len(nrow(d)) %in% o$.i[o$Day > first_event], "Check", "Measured after its event",
         "The event column was 1 on an earlier day.")
  }
  n_rows <- table(o$Key)
  single <- o$Key %in% names(n_rows)[n_rows == 1L]
  flag(seq_len(nrow(d)) %in% o$.i[single], "Check", "Single measurement",
       "One measurement adds nothing to a growth model.")

  # Missed measurements: no row at all on a day most of the arm's animals on
  # study were measured (a row without a volume is "No volume" instead).
  first_day <- stats::ave(o$Day, o$Key, FUN = min)
  missed <- do.call(rbind, lapply(split(p, p$Arm), function(a) {
    last <- tapply(a$Day, a$Key, max)
    do.call(rbind, lapply(sort(unique(a$Day)), function(t) {
      seen <- unique(a$Key[a$Day == t])
      on <- union(names(last)[last > t & tapply(a$Day, a$Key, min)[names(last)] < t], seen)
      if (length(on) < 2L || length(seen) / length(on) < 0.75) return(NULL)
      gone <- setdiff(on, seen)
      if (length(gone)) data.frame(Key = gone, Day = t, stringsAsFactors = FALSE)
    }))
  }))
  if (!is.null(missed) && nrow(missed)) {
    m <- p[match(missed$Key, p$Key), c("Arm", "ID", "Cage")]
    found[[length(found) + 1L]] <- data.frame(Row = NA_integer_, m, Day = missed$Day,
      Severity = "Check", Check = "Missed measurement",
      Detail = "No row on a day most of its arm was measured, though measured before and after.",
      stringsAsFactors = FALSE)
  }

  # ---- Notes ----------------------------------------------------------------
  usual_first <- tapply(first_day[!duplicated(o$Key)], o$Arm[!duplicated(o$Key)],
                        function(x) as.numeric(names(which.max(table(x)))))
  late <- !duplicated(o$Key) & first_day != usual_first[o$Arm]
  flag(seq_len(nrow(d)) %in% o$.i[late], "Note", "Staggered start",
       "The animal's first day differs from its arm's usual first day.")
  animals_per_arm <- tapply(o$Key, o$Arm, function(k) length(unique(k)))
  small <- names(animals_per_arm)[animals_per_arm < 3L]
  flag(seq_len(nrow(d)) %in% o$.i[!duplicated(o$Key) & o$Arm %in% small], "Note", "Small arm",
       "The arm has fewer than 3 animals.")
  ids <- unique(o[, c("Arm", "ID")])
  shared_ids <- unique(ids$ID[duplicated(ids$ID)])
  flag(seq_len(nrow(d)) %in% o$.i[!duplicated(o$Key) & o$ID %in% shared_ids], "Note",
       "IDs reused across arms", "The analyses keep them apart by arm.")
  if (has_col(df, cage_column) && length(lab <- reused_cage_labels(o$Cage, o$Arm))) {
    flag(seq_len(nrow(d)) %in% o$.i[!duplicated(o$Key) & o$Cage %in% lab], "Note",
         "Cage labels reused across arms",
         "Read as cages holding several arms unless numbered within each arm (label_cages()).")
  }

  flagged <- if (length(found)) do.call(rbind, found) else
    data.frame(Row = integer(0), Arm = character(0), ID = character(0), Cage = character(0),
               Day = numeric(0), Severity = character(0), Check = character(0),
               Detail = character(0), stringsAsFactors = FALSE)
  sev <- c("Problem", "Check", "Note")
  flagged <- flagged[order(match(flagged$Severity, sev), flagged$Check, flagged$Row), , drop = FALSE]
  rownames(flagged) <- NULL
  if (!has_col(df, cage_column)) flagged$Cage <- NULL

  what <- c(
    "Missing ID, arm or day" = "Rows that cannot be placed in an animal's record",
    "Duplicate measurement" = "An animal measured twice on one day",
    "Negative volume" = "Volumes below zero",
    "Non-positive body weight" = "Body weights of zero or below",
    "No volume" = "Rows with an animal and a day but no volume",
    "Width above length" = "Caliper width greater than length",
    "Large jump" = paste0("Volume changes of more than ", format(jump_fold), "-fold between measurements"),
    "Weight drop" = paste0("Body weight falling more than ", round(100 * weight_drop), " % between weighings"),
    "Measured after its event" = "Rows after an animal's event",
    "Single measurement" = "Animals measured once",
    "Missed measurement" = "Gaps in an animal's schedule while it was on study",
    "Staggered start" = "Animals starting later than their arm",
    "Small arm" = "Arms with fewer than 3 animals",
    "IDs reused across arms" = "The same ID in several arms",
    "Cage labels reused across arms" = "Every arm's cage labels from one set")
  checks <- unique(flagged[, c("Severity", "Check")])
  summary <- data.frame(
    Severity = checks$Severity, Check = checks$Check,
    Animals = vapply(checks$Check, function(k) {
      f <- flagged[flagged$Check == k, , drop = FALSE]
      length(unique(paste(f$Arm, f$ID, f$Cage)))
    }, integer(1)),
    Rows = vapply(checks$Check, function(k) sum(flagged$Check == k), integer(1)),
    What = unname(what[checks$Check]), stringsAsFactors = FALSE)
  rownames(summary) <- NULL
  list(summary = summary, flagged = flagged)
}
