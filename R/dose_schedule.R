# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Dosing-schedule annotation — CODE_REVIEW.md D.15.
#
# Replaces `plot_treatments()`, which drew a standalone panel (y = arm, x = day,
# down-triangles at each dosing day) intended to be stacked beneath a growth
# curve with manually aligned axes. Two things made that awkward enough to remove:
# it required a separate schedule data frame uploaded as its own file, and
# keeping two independently-built plots aligned is fragile.
#
# The constraint that shapes the replacement: **dosing timing is not recoverable
# from measurement data.** In every dataset shipped with this package `Dose` is
# constant within an animal — it encodes dose *level*, not dose *timing* — and
# the dosing days (1, 5, 9, 13) are not even on the measurement grid
# (0, 2, 4, 6, ...). So the user must supply the timing, and the design question
# is not what the annotation looks like but how to get a handful of integers out
# of them without a file-upload path. Hence a single text field.
#
# The marks are drawn in a strip below the panel rather than across it. Vertical
# lines through the data region compete with the curves, which is precisely why
# the original put the schedule in its own panel; a margin strip keeps that
# separation without needing two plots.

#' Parse a dosing schedule from free text
#'
#' Accepts either a bare list of days applying to every treated group, or one
#' line per group. Blank input returns `NULL`, which callers treat as "no
#' annotation requested" rather than an error.
#'
#' @param text Character scalar. Either `"1, 5, 9, 13"` (applies to all groups in
#'   `groups`) or per-group lines such as `"DrugA: 1,5,9,13"`, one per line or
#'   separated by `;`. Group names are matched case-insensitively.
#' @param groups Optional character vector of treatment-group names present in
#'   the data. Used to expand the bare form and to validate the per-group form.
#' @return A data frame with `Treatment` and `Day`, carrying a `scope` attribute
#'   of `"shared"` or `"per_group"`; or `NULL` when nothing was supplied.
#'   Unmatched group names raise a warning and are dropped — silently ignoring
#'   them is the failure mode this package has repeatedly found (§K.2).
#' @export
#' @examples
#' parse_dose_schedule("1, 5, 9, 13", groups = c("Control", "DrugA"))
#' parse_dose_schedule("DrugA: 1,5,9\nCombo: 1,5", groups = c("DrugA", "Combo"))
parse_dose_schedule <- function(text, groups = NULL) {
  if (is.null(text) || !length(text)) return(NULL)
  text <- paste(as.character(text), collapse = "\n")
  if (!nzchar(trimws(text))) return(NULL)

  lines <- unlist(strsplit(text, "[\n;]+"))
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(NULL)

  parse_days <- function(x) {
    v <- suppressWarnings(as.numeric(trimws(unlist(strsplit(x, "[,[:space:]]+")))))
    v <- v[!is.na(v)]
    sort(unique(v))
  }

  has_label <- grepl(":", lines, fixed = TRUE)

  if (!any(has_label)) {
    # Bare form: one schedule for every group.
    days <- parse_days(paste(lines, collapse = ","))
    if (!length(days)) {
      warning("No dosing days could be read from '", text,
              "'. Expected numbers such as '1, 5, 9, 13'.", call. = FALSE)
      return(NULL)
    }
    grp <- if (is.null(groups) || !length(groups)) NA_character_ else groups
    out <- expand.grid(Treatment = grp, Day = days,
                       stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    out <- out[order(out$Treatment, out$Day), c("Treatment", "Day")]
    rownames(out) <- NULL
    return(structure(out, scope = "shared"))
  }

  # Per-group form. A line without a colon in an otherwise labelled block is
  # ambiguous, so it is reported rather than guessed at.
  if (any(!has_label)) {
    warning("Ignoring line(s) without a group label: ",
            paste(lines[!has_label], collapse = " | "),
            ". Use 'Group: days' on every line, or give days alone for all groups.",
            call. = FALSE)
    lines <- lines[has_label]
  }

  rows <- list()
  unmatched <- character(0)
  for (ln in lines) {
    parts <- strsplit(ln, ":", fixed = TRUE)[[1]]
    nm   <- trimws(parts[1])
    days <- parse_days(paste(parts[-1], collapse = ":"))
    if (!length(days)) {
      warning("No dosing days read for group '", nm, "'.", call. = FALSE)
      next
    }
    if (!is.null(groups) && length(groups)) {
      hit <- groups[tolower(groups) == tolower(nm)]
      if (!length(hit)) { unmatched <- c(unmatched, nm); next }
      nm <- hit[1]
    }
    rows[[length(rows) + 1L]] <- data.frame(Treatment = nm, Day = days,
                                            stringsAsFactors = FALSE)
  }
  if (length(unmatched)) {
    warning("Dosing schedule names not found in the data and dropped: ",
            paste(unique(unmatched), collapse = ", "),
            ". Present groups: ", paste(groups, collapse = ", "), call. = FALSE)
  }
  if (!length(rows)) return(NULL)

  out <- do.call(rbind, rows)
  out <- out[order(out$Treatment, out$Day), , drop = FALSE]
  rownames(out) <- NULL
  structure(out, scope = "per_group")
}

#' Choose between discrete marks and a dosing window
#'
#' A twenty-one-day daily schedule drawn as twenty-one triangles is noise; the
#' same schedule drawn as a shaded window reads instantly. Conversely a band
#' spanning four intermittent doses hides the fact that dosing was intermittent.
#'
#' @param schedule A data frame from [parse_dose_schedule()].
#' @param max_marks Above this many doses for any one group, prefer a window.
#' @return `"rug"` or `"band"`.
#' @export
dose_schedule_style <- function(schedule, max_marks = 8L) {
  if (is.null(schedule) || !nrow(schedule)) return("rug")
  per <- split(schedule$Day, schedule$Treatment)
  n_max <- max(vapply(per, length, integer(1)))
  # Median spacing of 1 day or less means effectively continuous dosing.
  gap_min <- suppressWarnings(min(vapply(per, function(d) {
    if (length(d) < 2L) return(Inf)
    stats::median(diff(sort(d)))
  }, numeric(1))))
  if (n_max > max_marks || (is.finite(gap_min) && gap_min <= 1)) "band" else "rug"
}

# Recover the y-axis transform so marks can be placed in a strip below the panel
# on log, sqrt or identity scales alike. ggplot2 renamed this field in 3.5, so
# both are checked; an unrecognised transform falls back to identity, which is
# wrong only in placement, never in the data.
.me_y_trans <- function(built) {
  sc <- tryCatch(built$layout$panel_scales_y[[1]], error = function(e) NULL)
  tr <- NULL
  if (!is.null(sc)) {
    tr <- tryCatch(sc$transformation, error = function(e) NULL)
    if (is.null(tr)) tr <- tryCatch(sc$trans, error = function(e) NULL)
  }
  if (is.null(tr) || is.null(tr$inverse)) list(inverse = identity) else tr
}

#' Annotate a time-axis plot with a dosing schedule
#'
#' Draws dosing days in a strip below the plotting panel: down-triangles for an
#' intermittent schedule, a filled window for a continuous one. Groups with
#' different schedules are stacked, one row each, so the annotation carries the
#' same information the removed `plot_treatments()` panel did without a second
#' plot to keep aligned.
#'
#' @param p A ggplot whose x aesthetic is study day.
#' @param schedule A data frame from [parse_dose_schedule()], or `NULL` to return
#'   `p` unchanged.
#' @param style `"auto"` (default), `"rug"`, or `"band"`.
#' @param colors Optional named colour vector keyed by treatment, so marks match
#'   the curve colours. Groups absent from it fall back to grey.
#' @param max_marks Passed to [dose_schedule_style()].
#' @param label Axis-side label for the strip; `NULL` for none.
#' @return The plot with the annotation added.
#' @export
annotate_dose_schedule <- function(p, schedule, style = c("auto", "rug", "band"),
                                   colors = NULL, max_marks = 8L,
                                   label = "Dosing") {
  if (is.null(schedule) || !is.data.frame(schedule) || !nrow(schedule)) return(p)
  if (!inherits(p, "ggplot")) return(p)
  style <- match.arg(style)
  if (identical(style, "auto")) style <- dose_schedule_style(schedule, max_marks)

  built <- tryCatch(ggplot2::ggplot_build(p), error = function(e) NULL)
  if (is.null(built)) return(p)
  yr <- tryCatch(built$layout$panel_params[[1]]$y.range, error = function(e) NULL)
  if (is.null(yr) || !all(is.finite(yr))) return(p)

  inv  <- .me_y_trans(built)$inverse
  span <- diff(yr)
  if (!is.finite(span) || span <= 0) return(p)

  grps <- unique(schedule$Treatment)
  grps <- grps[!is.na(grps)]

  # Stack only when the schedules actually differ. The bare input form expands to
  # one entry per group with identical days, so counting groups would draw N
  # identical rows and imply a per-arm difference that does not exist. Check the
  # days themselves, not the number of groups, so a per-group spec whose arms
  # happen to coincide also collapses to a single row.
  per_days <- lapply(split(schedule$Day, schedule$Treatment),
                     function(d) sort(unique(d)))
  all_same <- length(per_days) <= 1L ||
    length(unique(lapply(per_days, paste, collapse = ","))) == 1L
  stacked <- length(grps) > 1L && !all_same

  rows <- if (stacked) grps else NA_character_
  n_rows <- length(rows)

  pad  <- 0.07 * span               # gap between the data and the strip
  step <- 0.055 * span              # row pitch within the strip
  half <- 0.020 * span              # half-height of a band row

  y_at <- function(i) inv(yr[1] - pad - (i - 1L) * step)

  fill_for <- function(g) {
    if (is.na(g) || is.null(colors) || !length(colors)) return("grey35")
    if (!is.null(names(colors)) && g %in% names(colors)) unname(colors[[g]]) else "grey35"
  }

  layers <- list()
  for (i in seq_len(n_rows)) {
    g <- rows[i]
    days <- if (stacked) schedule$Day[schedule$Treatment == g] else unique(schedule$Day)
    days <- sort(unique(days))
    if (!length(days)) next
    col <- fill_for(g)

    if (identical(style, "band")) {
      layers <- c(layers, list(ggplot2::annotate(
        "rect", xmin = min(days), xmax = max(days),
        ymin = inv(yr[1] - pad - (i - 1L) * step - half),
        ymax = inv(yr[1] - pad - (i - 1L) * step + half),
        fill = col, alpha = 0.35)))
    } else {
      layers <- c(layers, list(ggplot2::annotate(
        "point", x = days, y = y_at(i),
        shape = 25, size = 2.4, fill = col, colour = col)))
    }
  }
  if (!length(layers)) return(p)

  lowest <- inv(yr[1] - pad - (n_rows - 1L) * step - 2 * half)

  # Row labels sit at the left edge of the panel when schedules differ, so a
  # reader can tell which arm each row belongs to without a legend.
  if (stacked) {
    xr <- tryCatch(built$layout$panel_params[[1]]$x.range, error = function(e) NULL)
    if (!is.null(xr) && all(is.finite(xr))) {
      for (i in seq_len(n_rows)) {
        layers <- c(layers, list(ggplot2::annotate(
          "text", x = xr[1], y = y_at(i), label = rows[i],
          hjust = 1.05, vjust = 0.5, size = 2.6, colour = "grey30")))
      }
    }
  } else if (!is.null(label) && nzchar(label)) {
    xr <- tryCatch(built$layout$panel_params[[1]]$x.range, error = function(e) NULL)
    if (!is.null(xr) && all(is.finite(xr))) {
      layers <- c(layers, list(ggplot2::annotate(
        "text", x = xr[1], y = y_at(1L), label = label,
        hjust = 1.05, vjust = 0.5, size = 2.6, colour = "grey30")))
    }
  }

  p +
    layers +
    ggplot2::expand_limits(y = lowest) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 22))
}
