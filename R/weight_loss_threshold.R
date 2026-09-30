#' Weight Loss Time-to-Threshold Analysis
#'
#' Performs Kaplan-Meier and optional Cox PH analysis for time to a specified
#' percentage body weight loss threshold.
#'
#' @param df Data frame with longitudinal data.
#' @param weight_column Name of the body weight column.
#' @param time_column Name of the time/day column.
#' @param treatment_column Name of the treatment group column.
#' @param id_column Name of the mouse/subject ID column.
#' @param cage_column Name of the cage column. NULL to omit. Included in the
#'   composite mouse key so reused IDs across cages are not collapsed.
#' @param volume_column Name of the tumor volume column. NULL to skip tumor adjustment.
#' @param adjust_tumor_weight Logical; subtract estimated tumor weight.
#' @param tumor_density Density in g/cm³ (default 1.0).
#' @param volume_units Units of the volume column, \code{"mm3"} or \code{"cm3"}.
#'   Required when \code{adjust_tumor_weight = TRUE} and a volume column is used,
#'   because volume is converted to mass and subtracted from body weight: the
#'   wrong unit scales that correction by 1000. Units were inferred from the data
#'   before v0.25.0, which read small-tumour mm3 studies as cm3 (CODE_REVIEW.md
#'   R20.83); the data are now only checked against the declared unit, with a
#'   warning on disagreement and an error when the implied tumour mass exceeds
#'   half the body weight.
#' @param threshold Fractional weight loss threshold (default 0.20 = 20%).
#' @param baseline_day Day to use as baseline for initial weight. NULL = first observation per mouse.
#' @param reference_group Name of the control/reference group.
#' @param removal_reason_column Optional column recording why each animal left
#'   the study (read from its last row with a non-empty value). Without it,
#'   an animal whose record ends before the threshold is censored at its last
#'   day, whatever the cause (see "Removals"). CODE_REVIEW.md R20.23.
#' @param weight_loss_reasons Values of \code{removal_reason_column} meaning
#'   the animal was removed for weight loss or body condition. Such an animal
#'   counts as a weight-loss event at its last day even if its recorded weight
#'   stayed above the threshold.
#' @param planned_end_reasons Values meaning a planned end: the end of the
#'   study, a scheduled sacrifice or a data cut. These are ordinary
#'   (administrative) censoring. Every other non-empty reason, such as tumour
#'   burden or death, is a competing removal.
#' @return A list with: event_data (per animal, including \code{Censor_Type} —
#'   "event", "administrative" or "competing_removal" — and \code{Status_CR}),
#'   km_fit, km_summary, log_rank, cox_model, cox_summary, cox_method, ph_test,
#'   \code{cuminc} (Aalen-Johansen cumulative incidence with competing
#'   removals, when a removal-reason column identifies any),
#'   \code{censoring_summary}, \code{n_competing_risk} and \code{assumption}
#'   (how early ends were treated, in words).
#'
#' @section Removals:
#' Without a removal-reason column an animal whose record ends before it
#' reaches the threshold is censored at its last day. That assumes its later
#' weight-loss risk was like that of the animals still on study. Before
#' v0.27.0 every such record was a competing "removal", so staggered enrolment
#' or a data cut halved the Aalen-Johansen incidence (R20.23). With a
#' removal-reason column, removals for weight loss are events, planned ends
#' are censoring, and other removals (tumour burden, death) are competing
#' risks: \code{cuminc} then gives the cumulative incidence of weight loss
#' allowing for them, which \code{1 - km_fit} overstates.
#' @export
weight_loss_threshold <- function(df,
                                  weight_column    = "Weight",
                                  time_column      = "Day",
                                  treatment_column = "Treatment",
                                  id_column        = "ID",
                                  cage_column      = NULL,
                                  volume_column    = NULL,
                                  adjust_tumor_weight = TRUE,
                                  tumor_density    = 1.0,
                                  volume_units     = NULL,
                                  threshold        = 0.20,
                                  baseline_day     = NULL,
                                  reference_group  = NULL,
                                  removal_reason_column = NULL,
                                  weight_loss_reasons   = character(0),
                                  planned_end_reasons   = character(0)) {

  # --- Validate ---
  required <- c(weight_column, time_column, treatment_column, id_column)
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # --- Build working data ---
  # CODE_REVIEW.md R3.25 — cage is part of the mouse identity everywhere else
  # in the package; without it, same-numeric-ID mice in different cages within
  # one treatment arm collapse into a single subject.
  has_cage <- !is.null(cage_column) && cage_column %in% names(df)
  wd <- data.frame(
    ID        = as.character(df[[id_column]]),
    Treatment = as.character(df[[treatment_column]]),
    Cage      = if (has_cage) as.character(df[[cage_column]]) else "1",
    Day       = as.numeric(df[[time_column]]),
    Weight    = as.numeric(df[[weight_column]]),
    stringsAsFactors = FALSE
  )
  has_reason <- !is.null(removal_reason_column) &&
    removal_reason_column %in% names(df)
  if (!is.null(removal_reason_column) && !has_reason) {
    stop("Removal-reason column '", removal_reason_column, "' not found.",
         call. = FALSE)
  }
  wd$Reason <- if (has_reason) trimws(as.character(df[[removal_reason_column]])) else NA_character_

  has_volume <- !is.null(volume_column) && volume_column %in% names(df)
  if (adjust_tumor_weight && has_volume) {
    # CODE_REVIEW.md R3.30 — resolve units explicitly rather than assuming mm³.
    # R20.22: volume filled in on weighing days without a calliper reading.
    vol <- me_fill_volume(make_mouse_key(wd$Treatment, wd$ID, wd$Cage), wd$Day,
                          as.numeric(df[[volume_column]]))
    volume_units <- resolve_volume_units(vol, volume_units)
    tumor_mass <- volume_to_mass(vol, tumor_density, volume_units)
    check_tumor_mass_plausible(tumor_mass, wd$Weight, volume_units)
    wd$Weight <- wd$Weight - tumor_mass
  }

  wd <- wd[!is.na(wd$Weight) & !is.na(wd$Day), ]

  # Composite key: ensures IDs shared across treatment groups are treated as
  # distinct mice (common when ear tags / cage labels are reused per group).
  wd$.MouseKey <- make_mouse_key(wd$Treatment, wd$ID, wd$Cage)
  wd <- wd[order(wd$.MouseKey, wd$Day), ]

  # --- Compute baseline weight per mouse ---
  if (!is.null(baseline_day)) {
    bl <- wd[wd$Day == baseline_day, ]
    # For mice without an observation on baseline_day, use their earliest
    missing_keys <- setdiff(unique(wd$.MouseKey), unique(bl$.MouseKey))
    if (length(missing_keys) > 0) {
      fallback <- do.call(rbind, lapply(missing_keys, function(key) {
        sub <- wd[wd$.MouseKey == key, ]
        sub[1, , drop = FALSE]
      }))
      bl <- rbind(bl, fallback)
    }
  } else {
    bl <- do.call(rbind, lapply(unique(wd$.MouseKey), function(key) {
      sub <- wd[wd$.MouseKey == key, ]
      sub[1, , drop = FALSE]
    }))
  }
  baseline_weights <- stats::setNames(bl$Weight, bl$.MouseKey)

  # --- Determine event time per mouse ---
  # CODE_REVIEW.md R20.21 -- `bw` is unnamed here: a named baseline turned the
  # data frame's row names into mouse keys, coxphf() failed on them, and the
  # failure left a non-converged coxph fit reported as "cox".
  event_list <- lapply(unique(wd$.MouseKey), function(key) {
    sub <- wd[wd$.MouseKey == key, ]
    bw <- unname(baseline_weights[key])
    threshold_weight <- bw * (1 - threshold)
    hit <- which(sub$Weight <= threshold_weight)
    reasons <- sub$Reason[!is.na(sub$Reason) & nzchar(sub$Reason)]
    reason  <- if (length(reasons)) reasons[length(reasons)] else NA_character_
    if (!length(hit) && !is.na(reason) && reason %in% weight_loss_reasons) {
      # R20.23: removed for weight loss or body condition -- an event at the
      # animal's last day, whatever its last recorded weight.
      hit <- nrow(sub)
    }
    if (length(hit) > 0) {
      # Event: first day at or below threshold
      data.frame(
        ID        = sub$ID[1],
        Treatment = sub$Treatment[1],
        Baseline_Weight = bw,
        Time      = sub$Day[hit[1]],
        Event     = 1L,
        Reason    = reason,
        stringsAsFactors = FALSE
      )
    } else {
      # Censored: last observation day
      data.frame(
        ID        = sub$ID[1],
        Treatment = sub$Treatment[1],
        Baseline_Weight = bw,
        Time      = max(sub$Day),
        Event     = 0L,
        Reason    = reason,
        stringsAsFactors = FALSE
      )
    }
  })
  event_df <- do.call(rbind, event_list)
  rownames(event_df) <- NULL

  # CODE_REVIEW.md R3.26 / R20.23 -- an animal euthanised for tumour burden is
  # a competing risk for weight loss, and 1 - KM, which keeps it at risk,
  # overstates the incidence. But which early ends are such removals is known
  # only from a removal-reason column. Without one an early end is ordinary
  # censoring: before v0.27.0 every record ending before the global last day
  # was a competing "removal", so staggered enrolment or a data cut halved the
  # Aalen-Johansen incidence. With reasons, planned ends are censoring and
  # other removals (tumour burden, death) compete with weight loss.
  competing <- has_reason & event_df$Event == 0L & !is.na(event_df$Reason) &
    !event_df$Reason %in% c(weight_loss_reasons, planned_end_reasons)
  event_df$Censor_Type <- ifelse(event_df$Event == 1L, "event",
                          ifelse(competing, "competing_removal", "administrative"))
  event_df$Status_CR <- ifelse(event_df$Event == 1L, 1L, ifelse(competing, 2L, 0L))

  censoring_summary <- as.data.frame(
    table(Treatment = event_df$Treatment, Censor_Type = event_df$Censor_Type)
  )

  n_competing <- sum(event_df$Status_CR == 2L)
  assumption <- if (!has_reason) paste(
    "No removal-reason column: an animal whose record ends before it reaches",
    "the threshold is censored at its last day, which assumes its later risk",
    "was like that of the animals still on study.")
  else paste(
    "Removal reasons: removals for weight loss are events, planned ends are",
    "censored, and", n_competing, "other removal(s) (e.g. tumour burden, death)",
    "are competing risks.")
  if (n_competing > 0L) {
    message(n_competing, " animal(s) were removed for other reasons before ",
            "reaching the weight-loss threshold. They are a competing risk in ",
            "`cuminc`; the Kaplan-Meier curve treats them as censored and ",
            "overstates weight-loss incidence.")
  }

  # Aalen-Johansen cumulative incidence via survfit() on a multi-state factor.
  cuminc <- tryCatch({
    cr_df <- event_df
    cr_df$Status_MS <- factor(cr_df$Status_CR, levels = c(0L, 1L, 2L),
                              labels = c("censored", "weight_loss", "removed"))
    survival::survfit(survival::Surv(Time, Status_MS) ~ Treatment, data = cr_df)
  }, error = function(e) NULL)

  # Set reference group
  event_df$Treatment <- as.factor(event_df$Treatment)
  if (!is.null(reference_group) && reference_group %in% levels(event_df$Treatment)) {
    event_df$Treatment <- stats::relevel(event_df$Treatment, ref = reference_group)
  }

  # --- Kaplan-Meier ---
  km_fit <- survival::survfit(
    survival::Surv(Time, Event) ~ Treatment,
    data = event_df
  )

  km_summary <- summary(km_fit)

  # --- Log-rank test ---
  log_rank <- tryCatch(
    survival::survdiff(
      survival::Surv(Time, Event) ~ Treatment,
      data = event_df
    ),
    error = function(e) NULL
  )

  # --- Cox PH (if ≥2 groups) ---
  # Now mirrors survival_statistics(): tries coxph; if rare events / complete
  # separation prevent convergence, falls back to coxphf (Firth). Also runs
  # survival::cox.zph() on a successful coxph fit so PH violations are
  # surfaced for time-to-weight-loss endpoints (where treatment-induced
  # acute loss followed by recovery routinely violates PH).
  cox_model <- NULL
  cox_summary <- NULL
  cox_method  <- NA_character_
  ph_test     <- NULL
  if (length(levels(event_df$Treatment)) >= 2) {
    # Detect complete separation (a group with no events). Firth provides
    # bias-reduced estimates when standard Cox is unstable.
    ev_by_grp <- tapply(event_df$Event, event_df$Treatment,
                        function(x) sum(x, na.rm = TRUE))
    # CODE_REVIEW.md R3.27 — only zero-event groups cause separation. A group
    # where every animal has an event is estimable by standard Cox; treating it
    # as separation sent ordinary data down the Firth path.
    has_separation <- any(ev_by_grp == 0L, na.rm = TRUE)
    fit_cox <- function() tryCatch(
      survival::coxph(survival::Surv(Time, Event) ~ Treatment, data = event_df),
      error = function(e) NULL)

    # Under separation standard Cox does not converge, so it is fitted only
    # when Firth is unavailable or fails, and then with its warning.
    if (!has_separation) cox_model <- fit_cox()
    if (is.null(cox_model) && requireNamespace("coxphf", quietly = TRUE)) {
      firth_fit <- tryCatch(
        coxphf::coxphf(
          survival::Surv(Time, Event) ~ Treatment,
          data = event_df
        ),
        error = function(e) NULL
      )
      if (!is.null(firth_fit)) {
        cox_model   <- firth_fit
        cox_summary <- firth_fit
        cox_method  <- "coxphf"
      }
    }
    if (is.null(cox_model) && has_separation) cox_model <- fit_cox()

    if (!is.null(cox_model) && is.na(cox_method)) {
      cox_summary <- summary(cox_model)
      cox_method  <- "cox"
      ph_test <- tryCatch(survival::cox.zph(cox_model),
                          error = function(e) NULL)
    }
  }

  list(
    event_data   = event_df,
    km_fit       = km_fit,
    cuminc            = cuminc,
    censoring_summary = censoring_summary,
    n_competing_risk  = n_competing,
    assumption        = assumption,
    km_summary   = km_summary,
    log_rank     = log_rank,
    cox_model    = cox_model,
    cox_summary  = cox_summary,
    cox_method   = cox_method,
    ph_test      = ph_test,
    threshold    = threshold,
    baseline_day = baseline_day
  )
}
