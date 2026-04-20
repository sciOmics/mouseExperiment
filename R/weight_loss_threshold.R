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
#' @param volume_column Name of the tumor volume column. NULL to skip tumor adjustment.
#' @param adjust_tumor_weight Logical; subtract estimated tumor weight.
#' @param tumor_density Density in g/cm³ (default 1.0).
#' @param threshold Fractional weight loss threshold (default 0.20 = 20%).
#' @param baseline_day Day to use as baseline for initial weight. NULL = first observation per mouse.
#' @param reference_group Name of the control/reference group.
#' @return A list with: event_data, km_fit, km_summary, log_rank, cox_model, cox_summary.
#' @export
weight_loss_threshold <- function(df,
                                  weight_column    = "Weight",
                                  time_column      = "Day",
                                  treatment_column = "Treatment",
                                  id_column        = "ID",
                                  volume_column    = NULL,
                                  adjust_tumor_weight = TRUE,
                                  tumor_density    = 1.0,
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
  wd <- data.frame(
    ID        = as.character(df[[id_column]]),
    Treatment = as.character(df[[treatment_column]]),
    Day       = as.numeric(df[[time_column]]),
    Weight    = as.numeric(df[[weight_column]]),
    stringsAsFactors = FALSE
  )

  has_volume <- !is.null(volume_column) && volume_column %in% names(df)
  if (adjust_tumor_weight && has_volume) {
    vol <- as.numeric(df[[volume_column]])
    wd$Weight <- wd$Weight - (vol / 1000 * tumor_density)
  }

  wd <- wd[!is.na(wd$Weight) & !is.na(wd$Day), ]
  wd <- wd[order(wd$ID, wd$Day), ]

  # --- Compute baseline weight per mouse ---
  if (!is.null(baseline_day)) {
    bl <- wd[wd$Day == baseline_day, ]
    # For mice without an observation on baseline_day, use their earliest
    missing_ids <- setdiff(unique(wd$ID), unique(bl$ID))
    if (length(missing_ids) > 0) {
      fallback <- do.call(rbind, lapply(missing_ids, function(id) {
        sub <- wd[wd$ID == id, ]
        sub[1, , drop = FALSE]
      }))
      bl <- rbind(bl, fallback)
    }
  } else {
    bl <- do.call(rbind, lapply(unique(wd$ID), function(id) {
      sub <- wd[wd$ID == id, ]
      sub[1, , drop = FALSE]
    }))
  }
  baseline_weights <- stats::setNames(bl$Weight, bl$ID)

  # --- Determine event time per mouse ---
  event_list <- lapply(unique(wd$ID), function(id) {
    sub <- wd[wd$ID == id, ]
    bw <- baseline_weights[id]
    threshold_weight <- bw * (1 - threshold)
    hit <- which(sub$Weight <= threshold_weight)
    if (length(hit) > 0) {
      # Event: first day at or below threshold
      data.frame(
        ID        = id,
        Treatment = sub$Treatment[1],
        Baseline_Weight = bw,
        Time      = sub$Day[hit[1]],
        Event     = 1L,
        stringsAsFactors = FALSE
      )
    } else {
      # Censored: last observation day
      data.frame(
        ID        = id,
        Treatment = sub$Treatment[1],
        Baseline_Weight = bw,
        Time      = max(sub$Day),
        Event     = 0L,
        stringsAsFactors = FALSE
      )
    }
  })
  event_df <- do.call(rbind, event_list)

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
  cox_model <- NULL
  cox_summary <- NULL
  if (length(levels(event_df$Treatment)) >= 2) {
    cox_model <- tryCatch({
      survival::coxph(
        survival::Surv(Time, Event) ~ Treatment,
        data = event_df
      )
    }, error = function(e) NULL)
    if (!is.null(cox_model)) {
      cox_summary <- summary(cox_model)
    }
  }

  list(
    event_data  = event_df,
    km_fit      = km_fit,
    km_summary  = km_summary,
    log_rank    = log_rank,
    cox_model   = cox_model,
    cox_summary = cox_summary,
    threshold   = threshold,
    baseline_day = baseline_day
  )
}
