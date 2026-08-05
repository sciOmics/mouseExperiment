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
#' @param volume_units Units of \code{volume_column}: "mm3", "cm3", or NULL
#'   (default) to auto-detect from the data and report the inference.
#' @param threshold Fractional weight loss threshold (default 0.20 = 20%).
#' @param baseline_day Day to use as baseline for initial weight. NULL = first observation per mouse.
#' @param reference_group Name of the control/reference group.
#' @return A list with: event_data (per animal, including \code{Censor_Type} —
#'   "event" / "administrative" / "early_removal" — and \code{Status_CR}),
#'   km_fit, km_summary, log_rank, cox_model, cox_summary, cox_method, ph_test,
#'   \code{cuminc} (Aalen-Johansen cumulative incidence treating removal before
#'   study end as a competing risk), \code{censoring_summary}, and
#'   \code{n_competing_risk}.
#'
#' @section Competing risks:
#' Animals removed for tumour burden before reaching the weight-loss threshold
#' are not independently censored — removal and weight loss share a cause. The
#' Kaplan-Meier estimator assumes a censored animal remains at risk, so
#' \code{1 - km_fit} overstates the cumulative incidence of weight loss
#' whenever such removals occur. Prefer \code{cuminc} for incidence statements
#' and use \code{censoring_summary} to see how much of the censoring is
#' informative.
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
                                  reference_group  = NULL) {

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

  has_volume <- !is.null(volume_column) && volume_column %in% names(df)
  if (adjust_tumor_weight && has_volume) {
    # CODE_REVIEW.md R3.30 — resolve units explicitly rather than assuming mm³.
    vol <- as.numeric(df[[volume_column]])
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
  event_list <- lapply(unique(wd$.MouseKey), function(key) {
    sub <- wd[wd$.MouseKey == key, ]
    bw <- baseline_weights[key]
    threshold_weight <- bw * (1 - threshold)
    hit <- which(sub$Weight <= threshold_weight)
    if (length(hit) > 0) {
      # Event: first day at or below threshold
      data.frame(
        ID        = sub$ID[1],
        Treatment = sub$Treatment[1],
        Baseline_Weight = bw,
        Time      = sub$Day[hit[1]],
        Event     = 1L,
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
        stringsAsFactors = FALSE
      )
    }
  })
  event_df <- do.call(rbind, event_list)

  # CODE_REVIEW.md R3.26 — censoring here is informative. An animal whose record
  # ends early usually ended because it was euthanised for tumour burden, which
  # is not independent of the weight-loss hazard: it is a COMPETING RISK.
  # 1 - KM therefore overestimates the cumulative incidence of weight loss,
  # because KM implicitly assumes a censored animal remains at risk.
  #
  # Distinguish the two censoring mechanisms so the reader can see how much of
  # the censoring is informative, and report an Aalen-Johansen cumulative
  # incidence alongside the KM curve when a competing event is identifiable.
  study_end <- max(wd$Day, na.rm = TRUE)
  event_df$Censor_Type <- ifelse(
    event_df$Event == 1L, "event",
    ifelse(event_df$Time >= study_end, "administrative", "early_removal")
  )

  # Multi-state status: 0 = still at risk at end, 1 = weight-loss event,
  # 2 = competing removal before the study ended.
  event_df$Status_CR <- ifelse(event_df$Event == 1L, 1L,
                        ifelse(event_df$Censor_Type == "early_removal", 2L, 0L))

  censoring_summary <- as.data.frame(
    table(Treatment = event_df$Treatment, Censor_Type = event_df$Censor_Type)
  )

  n_competing <- sum(event_df$Status_CR == 2L)
  if (n_competing > 0L) {
    message(n_competing, " animal(s) left the study before day ", study_end,
            " without reaching the weight-loss threshold. These are treated as ",
            "a competing risk in `cuminc`; the Kaplan-Meier curve treats them ",
            "as ordinary censoring and will overstate weight-loss incidence.")
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
    cox_model <- tryCatch({
      survival::coxph(
        survival::Surv(Time, Event) ~ Treatment,
        data = event_df
      )
    }, error = function(e) NULL)

    # Detect complete separation (a group with 0 or all events). Firth
    # provides bias-reduced estimates when standard Cox is unstable.
    ev_by_grp <- tapply(event_df$Event, event_df$Treatment,
                        function(x) sum(x, na.rm = TRUE))
    # CODE_REVIEW.md R3.27 — only zero-event groups cause separation. A group
    # where every animal has an event is estimable by standard Cox; treating it
    # as separation sent ordinary data down the Firth path.
    has_separation <- any(ev_by_grp == 0L, na.rm = TRUE)

    needs_firth <- is.null(cox_model) || has_separation
    if (needs_firth && requireNamespace("coxphf", quietly = TRUE)) {
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
