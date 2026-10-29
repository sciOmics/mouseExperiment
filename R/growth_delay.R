#' Tumour growth delay
#'
#' The classic delay endpoint: **T - C**, the median time for an arm's tumours
#' to reach `volume` minus the reference arm's, in days, and the delay as a
#' percentage of the reference's median, \eqn{100 (T - C) / C}.
#'
#' * **Each animal's time** is the day its tumour reached `volume`,
#'   interpolated on the log scale between the last measurement below it and
#'   the first at or above it. An animal that never reached it is censored at
#'   its last measurement.
#' * **Each arm's median** is the Kaplan-Meier median of those times, so
#'   censored animals count; it is `Inf` ("not reached") when fewer than half
#'   the arm's animals reached the volume.
#' * **The interval** for T - C is a percentile bootstrap over animals within
#'   each arm; an upper limit of `Inf` means the treated median was not reached
#'   in enough resamples.
#' * **The test** compares each arm's times with the reference arm's by a
#'   log-rank test, with the p-values corrected over the arms
#'   (`p_adjust_method`). It treats animals as independent: when whole cages
#'   were randomised, it overstates the evidence.
#' * **The specific growth delay** is T - C divided by the reference arm's
#'   tumour volume doubling time (the median of its animals' doubling times,
#'   from [tumor_doubling_time()]): the delay in doublings, which can be
#'   compared across models that grow at different rates. Its interval is the
#'   delay's, divided by the same doubling time.
#'
#' Animals are keyed by arm, ID and cage, so reused ear tags stay apart.
#'
#' @param df Study data, one row per measurement.
#' @param volume The volume each tumour is timed to, in the data's units.
#' @param time_column,volume_column,treatment_column,id_column Column names.
#' @param cage_column Cage column, or NULL.
#' @param reference_group The reference arm. NULL picks the control by name,
#'   as the other analyses do.
#' @param n_boot Bootstrap resamples for the interval; 0 for none.
#' @param seed Seed for the resamples. The caller's random state is left as
#'   it was.
#' @param p_adjust_method Correction of the log-rank p-values over the arms
#'   compared with the reference: `"holm"`, `"bonferroni"`, `"fdr"` or
#'   `"none"`.
#' @return A list:
#'   * `delay_table`: per arm, `N`, `N_Reached`, `Median_Days`, `Delay_Days`
#'     (T - C), `Delay_Lower`, `Delay_Upper`, `Delay_Pct`, the specific growth
#'     delay `SGD` with `SGD_Lower` and `SGD_Upper`, and the log-rank
#'     `P_Value` and `P_Adjusted` against the reference;
#'   * `per_animal`: each animal's `Time` and `Reached` (1 or 0);
#'   * `volume`, `reference_group`, `doubling_time_reference` (days),
#'     `p_adjust_method`.
#' @seealso [survival_events()], which times tumours to a volume limit for
#'   survival analysis.
#' @export
#' @examples
#' gd <- growth_delay(master_synthetic_data, volume = 1000, cage_column = "Cage",
#'                    reference_group = "Vehicle", n_boot = 200, seed = 1)
#' gd$delay_table
growth_delay <- function(df, volume, time_column = "Day", volume_column = "Volume",
                         treatment_column = "Treatment", id_column = "ID",
                         cage_column = NULL, reference_group = NULL,
                         n_boot = 2000L, seed = NULL,
                         p_adjust_method = c("holm", "bonferroni", "fdr", "none")) {
  p_adjust_method <- match.arg(p_adjust_method)
  for (col in c(time_column, volume_column, treatment_column, id_column)) {
    if (!has_col(df, col)) stop("Column '", col, "' is not in the data.", call. = FALSE)
  }
  if (!is.numeric(volume) || length(volume) != 1L || !(volume > 0)) {
    stop("`volume` must be one positive number.", call. = FALSE)
  }
  has_cage <- has_col(df, cage_column)
  d <- data.frame(Treatment = as.character(df[[treatment_column]]),
                  ID = as.character(df[[id_column]]),
                  Cage = if (has_cage) as.character(df[[cage_column]]) else NA_character_,
                  Day = suppressWarnings(as.numeric(df[[time_column]])),
                  Volume = suppressWarnings(as.numeric(df[[volume_column]])),
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$Treatment) & is.finite(d$Day) & is.finite(d$Volume), , drop = FALSE]
  d$Key <- if (has_cage) make_mouse_key(d$Treatment, d$ID, d$Cage) else make_mouse_key(d$Treatment, d$ID)
  d <- d[order(d$Key, d$Day), , drop = FALSE]
  reference_group <- me_resolve_reference(d$Treatment, reference_group, "growth_delay")

  per_animal <- do.call(rbind, lapply(split(d, d$Key), function(a) {
    hit <- which(a$Volume >= volume)[1]
    time <- if (is.na(hit)) max(a$Day)
            else if (hit == 1L) a$Day[1]
            else gd_interpolate(a$Day[hit - 1L], a$Day[hit], a$Volume[hit - 1L], a$Volume[hit], volume)
    data.frame(Treatment = a$Treatment[1], ID = a$ID[1], Cage = a$Cage[1],
               Time = time, Reached = as.integer(!is.na(hit)), stringsAsFactors = FALSE)
  }))
  rownames(per_animal) <- NULL
  if (!has_cage) per_animal$Cage <- NULL

  arms <- c(reference_group, setdiff(sort(unique(per_animal$Treatment)), reference_group))
  by_arm <- split(per_animal, factor(per_animal$Treatment, levels = arms))
  med <- vapply(by_arm, function(a) gd_km_median(a$Time, a$Reached), numeric(1))

  boot <- if (n_boot >= 2L) me_with_seed(seed, vapply(seq_len(n_boot), function(b) {
    vapply(by_arm, function(a) {
      i <- sample.int(nrow(a), nrow(a), replace = TRUE)
      gd_km_median(a$Time[i], a$Reached[i])
    }, numeric(1))
  }, numeric(length(arms))))
  if (!is.null(boot)) boot <- matrix(boot, nrow = length(arms), dimnames = list(arms, NULL))

  c_med <- med[[reference_group]]
  delay_table <- do.call(rbind, lapply(arms, function(a) {
    ref <- identical(a, reference_group)
    delay <- if (ref) 0 else med[[a]] - c_med
    ci <- if (ref || is.null(boot) || !is.finite(c_med)) c(NA_real_, NA_real_) else {
      dd <- boot[a, ] - boot[reference_group, ]
      dd <- dd[!is.nan(dd)]
      stats::quantile(dd, c(0.025, 0.975), names = FALSE, type = 1, na.rm = TRUE)
    }
    data.frame(Treatment = a, N = nrow(by_arm[[a]]), N_Reached = sum(by_arm[[a]]$Reached),
               Median_Days = med[[a]],
               Delay_Days = if (is.finite(c_med)) delay else NA_real_,
               Delay_Lower = if (ref) NA_real_ else ci[1],
               Delay_Upper = if (ref) NA_real_ else ci[2],
               Delay_Pct = if (is.finite(c_med) && c_med > 0) 100 * delay / c_med else NA_real_,
               stringsAsFactors = FALSE)
  }))
  rownames(delay_table) <- NULL

  # R23.2: a log-rank test of each arm's times against the reference's,
  # corrected over the arms.
  others <- setdiff(arms, reference_group)
  p_raw <- vapply(others, function(a) gd_logrank_p(per_animal, reference_group, a), numeric(1))
  p_adj <- rep(NA_real_, length(p_raw))
  ok <- !is.na(p_raw)
  if (any(ok)) p_adj[ok] <- stats::p.adjust(p_raw[ok], method = p_adjust_method)
  delay_table$P_Value <- c(NA_real_, unname(p_raw))[match(delay_table$Treatment, arms)]
  delay_table$P_Adjusted <- c(NA_real_, p_adj)[match(delay_table$Treatment, arms)]

  # R23.2: the specific growth delay, T - C in the reference arm's doubling
  # times. The per-animal fits' own warnings (an exact fit, too few points)
  # concern those fits, not the delay.
  dt <- tryCatch(suppressWarnings(tumor_doubling_time(
                   df, time_column = time_column, volume_column = volume_column,
                   treatment_column = treatment_column, id_column = id_column,
                   cage_column = if (has_cage) cage_column)),
                 error = function(e) NULL)
  dt_ref <- if (is.data.frame(dt)) stats::median(dt$doubling_time[as.character(dt$Treatment) == reference_group],
                                                 na.rm = TRUE) else NA_real_
  if (!is.finite(dt_ref) || dt_ref <= 0) dt_ref <- NA_real_
  delay_table$SGD <- delay_table$Delay_Days / dt_ref
  delay_table$SGD_Lower <- delay_table$Delay_Lower / dt_ref
  delay_table$SGD_Upper <- delay_table$Delay_Upper / dt_ref
  delay_table <- delay_table[, c("Treatment", "N", "N_Reached", "Median_Days", "Delay_Days",
                                 "Delay_Lower", "Delay_Upper", "Delay_Pct", "SGD", "SGD_Lower",
                                 "SGD_Upper", "P_Value", "P_Adjusted")]

  me_result(list(delay_table = delay_table, per_animal = per_animal, volume = volume,
                 reference_group = reference_group, n_boot = as.integer(n_boot),
                 doubling_time_reference = dt_ref, p_adjust_method = p_adjust_method),
            "me_growth_delay")
}

#' The log-rank p-value of one arm's times against the reference's
#' @noRd
gd_logrank_p <- function(per_animal, reference, arm) {
  s <- per_animal[per_animal$Treatment %in% c(reference, arm), , drop = FALSE]
  if (length(unique(s$Treatment)) < 2L || !any(s$Reached == 1L)) return(NA_real_)
  t <- tryCatch(survival::survdiff(survival::Surv(Time, Reached) ~ Treatment, data = s),
                error = function(e) NULL)
  if (is.null(t) || !is.finite(t$chisq)) return(NA_real_)
  stats::pchisq(t$chisq, df = 1, lower.tail = FALSE)
}

#' The day a tumour reached `v`, on the log scale between two measurements
#' @noRd
gd_interpolate <- function(t0, t1, v0, v1, v) {
  if (!(v0 > 0) || !(v1 > v0)) return(t1)
  t0 + (t1 - t0) * (log(v) - log(v0)) / (log(v1) - log(v0))
}

#' Kaplan-Meier median time; Inf when not reached
#' @noRd
gd_km_median <- function(time, event) {
  if (!length(time)) return(NA_real_)
  if (!any(event == 1L)) return(Inf)
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  m <- unname(summary(fit)$table["median"])
  if (is.na(m)) Inf else m
}
