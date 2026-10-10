# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Evaluable days (CODE_REVIEW.md R20-K, R20.1, R20.17).
#
# Endpoint TGI, synergy, the therapeutic window and the dose-response endpoint
# are reported only on days with enough animals still on study in every arm.
# Past that day an arm's mean is an extrapolation beyond most of its own
# animals: on the Combo demo the old default evaluated the control arm at day
# 32, where one animal remained, and read 28,544 mm3 off the model against
# about 3,000 observed. The rule was confirmed by the maintainer on 2026-09-30
# and is shared by every function so that the dashboard's tabs agree.

ME_EVAL_MIN_FRACTION <- 0.5
ME_EVAL_MIN_ANIMALS  <- 3L

#' The evaluable-day rule in words
#' @noRd
#' @keywords internal
me_evaluable_rule_text <- function() {
  sprintf(paste0("An arm is evaluable on a day when at least %s%% of its ",
                 "enrolled animals, and at least %d animals, are still on study ",
                 "(measured on that day or later)."),
          format(100 * ME_EVAL_MIN_FRACTION), ME_EVAL_MIN_ANIMALS)
}

#' Days on which every arm has enough animals still on study
#'
#' An arm is \emph{evaluable} on a day when at least 50 % of its enrolled
#' animals, and at least 3 animals, are still on study, meaning measured on
#' that day or later. \code{\link{analyze_drug_synergy}},
#' \code{\link{analyze_drug_synergy_over_time}},
#' \code{\link{therapeutic_window_metric}},
#' \code{\link{dose_response_statistics}} and the AUC path of
#' \code{\link{tumor_growth_statistics}} use only days on which every arm they
#' compare is evaluable. By default they evaluate the last such day, and a
#' requested day that is not evaluable is an error. Beyond that day an arm's
#' mean would rest on a minority of its animals, or be extrapolated past them
#' (CODE_REVIEW.md R20.1).
#'
#' Because an animal counts as on study up to its last measurement, each arm's
#' evaluable days run from the first study day to a last evaluable day, and so
#' do the days on which every arm is evaluable.
#'
#' @param df Long data frame, one row per measurement.
#' @param treatment_column,id_column,time_column Column names.
#' @param cage_column Optional cage column, part of the animal key
#'   (treatment + ID + cage).
#' @param volume_column Optional. When given, only rows with a finite volume
#'   count as measurements.
#' @param arms Arms that must be evaluable. Default: every arm in \code{df}.
#' @return A list:
#'   \describe{
#'     \item{rule}{The rule in words.}
#'     \item{min_fraction, min_animals}{Its thresholds (0.5 and 3).}
#'     \item{arms}{The arms considered.}
#'     \item{table}{One row per arm and measured day: \code{Treatment},
#'       \code{Day}, \code{N_Enrolled}, \code{N_On_Study}, \code{Evaluable}.}
#'     \item{days}{The days on which every arm is evaluable.}
#'     \item{last_day}{The last of them, or \code{NA} when there is none.}
#'     \item{last_day_by_arm}{Each arm's last evaluable day (\code{NA} when it
#'       has none).}
#'     \item{excluded}{The other measured days, each with the arms that failed
#'       the rule, e.g. "Vehicle 2/8".}
#'   }
#' @examples
#' data(master_synthetic_data)
#' ev <- evaluable_days(master_synthetic_data, volume_column = "Volume")
#' ev$last_day
#' ev$excluded
#' @export
evaluable_days <- function(df,
                           treatment_column = "Treatment",
                           id_column        = "ID",
                           time_column      = "Day",
                           cage_column      = NULL,
                           volume_column    = NULL,
                           arms             = NULL) {
  need <- c(treatment_column, id_column, time_column, volume_column)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  has_cage <- !is.null(cage_column) && cage_column %in% names(df)
  d <- data.frame(
    MouseKey  = if (has_cage) {
      make_mouse_key(as.character(df[[treatment_column]]),
                     as.character(df[[id_column]]),
                     as.character(df[[cage_column]]))
    } else {
      make_mouse_key(as.character(df[[treatment_column]]),
                     as.character(df[[id_column]]))
    },
    Treatment = as.character(df[[treatment_column]]),
    Day       = suppressWarnings(as.numeric(df[[time_column]])),
    stringsAsFactors = FALSE
  )
  ok <- is.finite(d$Day) & !is.na(d$Treatment)
  if (!is.null(volume_column)) {
    ok <- ok & is.finite(suppressWarnings(as.numeric(df[[volume_column]])))
  }
  me_evaluability(d[ok, , drop = FALSE], arms)
}

#' The evaluable-day rule on a standardised frame
#'
#' @param d Data frame with `MouseKey`, `Treatment` and `Day` (measurements
#'   only).
#' @param arms Arms that must be evaluable; `NULL` for all.
#' @return See [evaluable_days()].
#' @noRd
#' @keywords internal
me_evaluability <- function(d, arms = NULL) {
  if (is.null(arms)) arms <- sort(unique(as.character(d$Treatment)))
  arms <- unique(as.character(arms))
  absent <- setdiff(arms, d$Treatment)
  if (length(absent)) {
    stop("No measurements for arm(s): ", paste(absent, collapse = ", "),
         ".", call. = FALSE)
  }
  d <- d[d$Treatment %in% arms, , drop = FALSE]
  last <- stats::aggregate(Day ~ MouseKey + Treatment, data = d, FUN = max)
  days <- sort(unique(d$Day))

  tab <- do.call(rbind, lapply(arms, function(a) {
    la    <- last$Day[last$Treatment == a]
    n_enr <- length(la)
    n_on  <- vapply(days, function(t) sum(la >= t), integer(1L))
    data.frame(
      Treatment  = a,
      Day        = days,
      N_Enrolled = n_enr,
      N_On_Study = n_on,
      Evaluable  = n_on >= ME_EVAL_MIN_FRACTION * n_enr &
                   n_on >= ME_EVAL_MIN_ANIMALS,
      stringsAsFactors = FALSE)
  }))
  rownames(tab) <- NULL

  all_ok  <- vapply(days, function(t) all(tab$Evaluable[tab$Day == t]),
                    logical(1L))
  ev_days <- days[all_ok]
  last_by <- vapply(arms, function(a) {
    x <- tab$Day[tab$Treatment == a & tab$Evaluable]
    if (length(x)) max(x) else NA_real_
  }, numeric(1L))

  excl_days <- days[!all_ok]
  excluded <- data.frame(
    Day  = excl_days,
    Arms = vapply(excl_days, function(t) {
      bad <- tab[tab$Day == t & !tab$Evaluable, , drop = FALSE]
      paste(sprintf("%s %d/%d", bad$Treatment, bad$N_On_Study, bad$N_Enrolled),
            collapse = "; ")
    }, character(1L)),
    stringsAsFactors = FALSE)

  list(
    rule            = me_evaluable_rule_text(),
    min_fraction    = ME_EVAL_MIN_FRACTION,
    min_animals     = ME_EVAL_MIN_ANIMALS,
    arms            = arms,
    table           = tab,
    days            = ev_days,
    last_day        = if (length(ev_days)) max(ev_days) else NA_real_,
    last_day_by_arm = last_by,
    excluded        = excluded
  )
}

#' Choose, or check, the evaluation day
#'
#' @param ev Output of `me_evaluability()`.
#' @param requested A day, or `NULL` for the last evaluable day. A day that was
#'   not measured is moved to the closest measured day, with a warning.
#' @return The evaluation day. Stops when it is not evaluable.
#' @noRd
#' @keywords internal
me_resolve_eval_day <- function(ev, requested = NULL) {
  if (is.null(requested)) {
    if (!length(ev$days)) stop(me_no_evaluable_day_msg(ev), call. = FALSE)
    return(ev$last_day)
  }
  requested <- as.numeric(requested)[1L]
  measured  <- sort(unique(ev$table$Day))
  if (!requested %in% measured) {
    closest <- measured[which.min(abs(measured - requested))]
    warning("Day ", requested, " was not measured; using the closest measured ",
            "day, ", closest, ".", call. = FALSE)
    requested <- closest
  }
  if (!requested %in% ev$days) {
    bad <- ev$table[ev$table$Day == requested & !ev$table$Evaluable, , drop = FALSE]
    stop("Day ", requested, " is not evaluable: ",
         paste(sprintf("%s has %d of %d animals on study", bad$Treatment,
                       bad$N_On_Study, bad$N_Enrolled), collapse = "; "),
         ". ", ev$rule, " ",
         if (length(ev$days)) paste0("The last evaluable day is ", ev$last_day, ".")
         else "No day is evaluable.",
         call. = FALSE)
  }
  requested
}

#' @noRd
#' @keywords internal
me_no_evaluable_day_msg <- function(ev) {
  first <- min(ev$table$Day)
  bad <- ev$table[ev$table$Day == first & !ev$table$Evaluable, , drop = FALSE]
  paste0("No day is evaluable for ", paste(ev$arms, collapse = ", "), ". ",
         ev$rule,
         if (nrow(bad)) paste0(" On day ", first, ": ",
                               paste(sprintf("%s %d/%d", bad$Treatment,
                                             bad$N_On_Study, bad$N_Enrolled),
                                     collapse = "; "), ".")
         else "")
}
