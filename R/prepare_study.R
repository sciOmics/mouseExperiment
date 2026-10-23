# Preparing a study's data: the dashboard's steps, as package functions
# (CODE_REVIEW.md R20-AB).
#
# Before any analysis the dashboard reads the file, removes the animals, arms
# and cages the user excluded, turns dates into study days, calculates volume
# from calipers, labels arms by treatment and dose, and, for survival, reduces
# each animal to its time and event. Until 0.37.0 these steps existed only in
# the dashboard, so the same file analysed with the package directly could get
# different inputs. The dashboard now calls these functions, and its R script
# can replay them from the uploaded file.

#' Read a study CSV file, whatever its delimiter, decimal mark and encoding
#'
#' Reads what spreadsheet programs write: comma-, semicolon- or tab-separated
#' values; a decimal comma with semicolons or tabs (as Excel writes in much of
#' Europe); UTF-8, Latin-1 or UTF-16. A row whose number of fields differs from
#' the header's stops with a message naming the line, usually a value with an
#' unquoted separator in it. This is how the dashboard reads an upload.
#'
#' @param path Path to the file.
#' @return A data frame. Column names are made syntactic, as by
#'   [utils::read.csv()]; empty fields and "NA" are missing values.
#' @export
#' @examples
#' f <- system.file("sample_data", "combo_treatment_synthetic_data.csv",
#'                  package = "mouseExperiment")
#' head(read_study_csv(f))
read_study_csv <- function(path) {
  raw <- readBin(path, "raw", n = min(file.size(path), 262144L))
  enc <- if (length(raw) >= 2L && raw[1] == as.raw(0xFF) && raw[2] == as.raw(0xFE)) {
    "UTF-16LE"
  } else if (length(raw) >= 2L && raw[1] == as.raw(0xFE) && raw[2] == as.raw(0xFF)) {
    "UTF-16BE"
  } else if (length(raw) >= 3L && all(raw[1:3] == as.raw(c(0xEF, 0xBB, 0xBF)))) {
    "UTF-8-BOM"
  } else if (validUTF8(rawToChar(raw[raw != as.raw(0)]))) {
    "UTF-8"
  } else {
    "latin1"
  }
  con <- file(path, encoding = enc)
  head_lines <- tryCatch(readLines(con, n = 50L, warn = FALSE), finally = close(con))
  head_lines <- head_lines[nzchar(trimws(head_lines))]
  if (!length(head_lines)) stop("The file is empty.", call. = FALSE)
  # Fields outside quotes, per line, for each candidate delimiter.
  unquoted <- gsub('"[^"]*"', "", head_lines)
  counts <- vapply(c(",", ";", "\t"), function(d) {
    n <- lengths(regmatches(unquoted, gregexpr(d, unquoted, fixed = TRUE)))
    if (all(n == n[1]) && n[1] > 0) n[1] else 0L
  }, integer(1))
  sep <- if (any(counts > 0)) names(counts)[which.max(counts)] else ","
  dec <- if (sep != "," && any(grepl("[0-9],[0-9]", unquoted))) "," else "."
  # Every row must have the header's number of fields.
  con <- file(path, encoding = enc)
  nf <- tryCatch(utils::count.fields(con, sep = sep, quote = "\"", blank.lines.skip = TRUE,
                                     comment.char = ""), finally = close(con))
  if (length(nf) > 1L && any(nf != nf[1], na.rm = TRUE)) {
    bad <- which(nf != nf[1])[1]
    stop(sprintf(paste0(
      "Line %d has %d fields but the header has %d. A value probably contains ",
      "an unquoted '%s'; put that value in quotes, or remove the '%s'."),
      bad, nf[bad], nf[1], if (sep == "\t") "tab" else sep, if (sep == "\t") "tab" else sep),
      call. = FALSE)
  }
  utils::read.table(path, sep = sep, dec = dec, header = TRUE, quote = "\"",
                    comment.char = "", fill = FALSE, stringsAsFactors = FALSE,
                    check.names = TRUE, strip.white = TRUE, fileEncoding = enc,
                    na.strings = c("NA", ""))
}

#' Remove excluded animals, arms and cages
#'
#' An animal is identified by its arm, its ID and, when the study has cages,
#' its cage: ear tags are often reused in each arm or cage, so an ID alone
#' would remove an animal from every arm.
#'
#' @param df Study data, one row per measurement.
#' @param animals Animals to remove: a data frame with columns `arm` and `id`,
#'   and `cage` when `cage_column` is given. `arm` is the treatment, or
#'   "Treatment Dose" when `dose_column` is given (as [label_arms()] writes it).
#' @param treatments Values of `treatment_column` to remove.
#' @param arms Arms to remove, as [label_arms()] writes them ("Drug A 10"): one
#'   dose of a treatment given at several.
#' @param cages Cages to remove: values of `cage_column`, or a data frame with
#'   columns `arm` and `cage` to remove one arm's cage when cages are numbered
#'   within each arm ("C1" in Vehicle but not "C1" in Drug A).
#' @param treatment_column,id_column Column names.
#' @param cage_column,dose_column Column names, or NULL when the study has none.
#' @return `df` without the excluded rows.
#' @export
#' @examples
#' d <- master_synthetic_data
#' nrow(d)
#' d <- apply_exclusions(d, animals = data.frame(arm = "Vehicle", id = 1),
#'                       treatments = "Drug_B")
#' nrow(d)
apply_exclusions <- function(df, animals = NULL, treatments = NULL, cages = NULL,
                             treatment_column = "Treatment", id_column = "ID",
                             cage_column = NULL, dose_column = NULL, arms = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  has_cage <- has_col(df, cage_column)
  if (!is.null(animals) && NROW(animals) > 0L && has_col(df, id_column) &&
      has_col(df, treatment_column)) {
    animals <- as.data.frame(animals, stringsAsFactors = FALSE)
    if (!all(c("arm", "id") %in% names(animals)) || (has_cage && !"cage" %in% names(animals))) {
      stop("`animals` needs columns arm and id", if (has_cage) ", and cage", ".", call. = FALSE)
    }
    drop <- animal_key(animals$arm, animals$id, if (has_cage) animals$cage)
    key  <- animal_key(arm_labels(df, treatment_column, dose_column), df[[id_column]],
                       if (has_cage) df[[cage_column]])
    df <- df[!(key %in% drop), , drop = FALSE]
  }
  if (length(treatments) && has_col(df, treatment_column)) {
    df <- df[!(as.character(df[[treatment_column]]) %in% as.character(treatments)), , drop = FALSE]
  }
  if (length(arms) && has_col(df, treatment_column)) {   # R22.16
    df <- df[!(arm_labels(df, treatment_column, dose_column) %in% as.character(arms)), , drop = FALSE]
  }
  if (NROW(cages) && has_cage) {
    cage <- as.character(df[[cage_column]])
    if (is.data.frame(cages)) {
      # One arm's cage, for cages numbered within each arm (R22.3).
      if (!all(c("arm", "cage") %in% names(cages))) {
        stop("`cages` as a data frame needs columns arm and cage.", call. = FALSE)
      }
      drop <- animal_key(cages$arm, cages$cage)
      df <- df[!(animal_key(arm_labels(df, treatment_column, dose_column), cage) %in% drop), ,
               drop = FALSE]
    } else {
      df <- df[!(cage %in% as.character(cages)), , drop = FALSE]
    }
  }
  df
}

#' Make cage labels unique across the study
#'
#' Cages are often numbered within each arm: C1 and C2 in Vehicle, and C1 and
#' C2 again in every other arm. The analyses read a cage label as one cage, so
#' those would be two cages that each hold every arm, and the cage effect would
#' be estimated across arms. `label_cages()` puts each cage's arm in front of
#' its label ("Vehicle C1"), so each cage has its own label. Run it after
#' [apply_exclusions()] and before the analyses.
#'
#' From the labels alone, cages numbered within each arm cannot be told apart
#' from cages that really hold several arms, so the choice is the study's. The
#' analyses warn when every arm's cage labels come from one set.
#'
#' @param df Study data.
#' @param cage_column The cage column.
#' @param treatment_column,dose_column The treatment column, and the dose
#'   column or NULL. With a dose, the arm is "Treatment Dose", as
#'   [label_arms()] writes it.
#' @return `df` with the cage column relabelled, as character. A missing cage
#'   stays missing.
#' @export
#' @examples
#' d <- data.frame(ID = 1:8, Treatment = rep(c("Vehicle", "Drug A"), each = 4),
#'                 Cage = rep(c("C1", "C1", "C2", "C2"), 2))
#' label_cages(d)$Cage
label_cages <- function(df, cage_column = "Cage", treatment_column = "Treatment",
                        dose_column = NULL) {
  for (col in c(cage_column, treatment_column)) {
    if (!has_col(df, col)) stop("Column '", col, "' is not in the data.", call. = FALSE)
  }
  cage <- as.character(df[[cage_column]])
  df[[cage_column]] <- ifelse(is.na(cage), NA_character_,
                              paste(arm_labels(df, treatment_column, dose_column), cage))
  df
}

#' Add a study Day column computed from a date column
#'
#' The whole column decides how dates are read: a value such as 13/03/2025
#' settles day-first for every date, and a column with no such value stops
#' unless `date_order` says which order it uses. Dates are counted from
#' `reference_date`, the study's Day 0. Values that cannot be read get no Day,
#' and dates before the reference date are reported; both come as a warning.
#'
#' @param df Study data.
#' @param date_column The column of dates.
#' @param reference_date The study's Day 0, preferably as YYYY-MM-DD.
#' @param date_order "auto", "dmy" or "mdy": how to read d/d/yyyy dates when
#'   the column does not decide it.
#' @return `df` with a numeric `Day` column.
#' @seealso [calculate_dates()], which reads one date format given to it.
#' @export
#' @examples
#' d <- data.frame(ID = 1:3, Date = c("24/03/2025", "27/03/2025", "31/03/2025"))
#' study_days(d, "Date", "2025-03-24")
study_days <- function(df, date_column, reference_date,
                       date_order = c("auto", "dmy", "mdy")) {
  info <- study_days_info(df, date_column, reference_date, match.arg(date_order))
  if (length(info$notes)) warning(paste(info$notes, collapse = " "), call. = FALSE)
  info$data
}

#' Label each row's arm by treatment and dose
#'
#' When arms differ by dose, as in a dose series, an arm is its treatment and
#' dose, for example "Drug_A 10". Without a dose column the arm is the
#' treatment.
#'
#' @param df Study data.
#' @param treatment_column The treatment column.
#' @param dose_column The dose column, or NULL.
#' @param name The new column's name.
#' @return `df` with the column `name`.
#' @export
#' @examples
#' d <- label_arms(master_synthetic_data, "Treatment", "Dose", name = "Arm")
#' table(d$Arm)
label_arms <- function(df, treatment_column = "Treatment", dose_column = NULL,
                       name = "Arm") {
  if (!has_col(df, treatment_column)) {
    stop("Treatment column '", treatment_column, "' is not in the data.", call. = FALSE)
  }
  df[[name]] <- arm_labels(df, treatment_column, dose_column)
  df
}

#' Each animal's survival time and event
#'
#' Reduces longitudinal data to one row per animal, keyed by arm, ID and cage,
#' with `Time` and `Event` (1 = event, 0 = censored):
#' * with `event_column`, an animal's event is at its first flagged day; an
#'   animal flagged on every row is an event at its last day;
#' * with `volume_limit`, an animal's event is the first day its tumour
#'   reached the limit, and an animal that never did is censored at its last
#'   measurement.
#'
#' With `removal_reason_column` and `death_reasons`, an animal removed for one
#' of those reasons (death from another cause) and not already an event
#' becomes one, at its last measurement.
#'
#' @param df Study data, one row per measurement.
#' @param time_column,treatment_column,id_column Column names. Give the arm
#'   column (see [label_arms()]) as `treatment_column` when arms differ by dose.
#' @param cage_column Cage column, or NULL.
#' @param event_column A 0/1 event column, or NULL when `volume_limit` defines
#'   events.
#' @param volume_column,volume_limit The volume column and the limit, in the
#'   data's units, when the limit defines events.
#' @param removal_reason_column,death_reasons The column of removal reasons,
#'   and the reasons that count as deaths.
#' @return One row per animal (its event row, or its last), with `Time`,
#'   `Event` and `.animal_key`. The attribute `death_events_added` is the
#'   number of deaths added from the removal reasons (NA when not used).
#' @export
#' @examples
#' ev <- survival_events(master_synthetic_data, volume_column = "Volume",
#'                       volume_limit = 1500)
#' table(ev$Treatment, ev$Event)
survival_events <- function(df, time_column = "Day", treatment_column = "Treatment",
                            id_column = "ID", cage_column = NULL, event_column = NULL,
                            volume_column = NULL, volume_limit = NULL,
                            removal_reason_column = NULL, death_reasons = NULL) {
  for (col in c(time_column, treatment_column, id_column)) {
    if (!has_col(df, col)) stop("Column '", col, "' is not in the data.", call. = FALSE)
  }
  if (!has_col(df, cage_column)) cage_column <- NULL
  by_volume <- !is.null(volume_limit)
  if (by_volume && !has_col(df, volume_column)) {
    stop("Events from a volume limit need the volume column.", call. = FALSE)
  }
  if (!by_volume && !has_col(df, event_column)) {
    stop("Give an event column, or a volume column and a volume limit.", call. = FALSE)
  }
  out <- survival_animal_rows(df, time_column, treatment_column, id_column, cage_column,
                              event_col = if (!by_volume) event_column)
  if (by_volume) {
    ev <- volume_threshold_events(df, time_column, treatment_column, id_column, cage_column,
                                  vol_col = volume_column, threshold = volume_limit)
    m <- match(out$.animal_key, ev$.animal_key)
    out$Time  <- ev$Time[m]
    out$Event <- ev$Event[m]
  } else {
    out$Time  <- suppressWarnings(as.numeric(out[[time_column]]))
    out$Event <- suppressWarnings(as.numeric(out[[event_column]]))
  }
  added <- NA_integer_
  if (length(death_reasons) && has_col(df, removal_reason_column)) {
    reasons <- animal_removal_reasons(df, treatment_column, id_column, cage_column,
                                      removal_reason_column, time_column)
    died <- names(reasons)[!is.na(reasons) & reasons %in% death_reasons]
    hit  <- out$.animal_key %in% died & (is.na(out$Event) | out$Event == 0)
    out$Event <- ifelse(hit, 1, out$Event)
    added <- sum(hit)
  }
  attr(out, "death_events_added") <- added
  out
}
