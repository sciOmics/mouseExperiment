#' Weight-Corrected Tumor Growth Inhibition
#'
#' Computes TGI using only mice that stayed within a body weight safety threshold.
#' Prevents drugs from appearing effective simply because toxicity stopped tumor growth.
#'
#' @param df Data frame with longitudinal data.
#' @param weight_column Name of the body weight column.
#' @param volume_column Name of the tumor volume column.
#' @param time_column Name of the time/day column.
#' @param treatment_column Name of the treatment group column.
#' @param id_column Name of the mouse/subject ID column.
#' @param adjust_tumor_weight Logical; subtract estimated tumor weight.
#' @param tumor_density Density in g/cm³ (default 1.0).
#' @param endpoint_day Day at which TGI is evaluated. NULL uses the maximum
#'   observed day.
#' @param endpoint_method How each arm's endpoint volume is obtained:
#'   "model" (default, log-scale LMM marginal means using every observation),
#'   "last_obs", or "survivors" (pre-0.8.0 behaviour; warns, and understates
#'   efficacy because it conditions on survival). See CODE_REVIEW.md R3.5 / G.3.
#' @param reference_group Name of the control/reference group.
#' @param safety_threshold Fractional weight loss threshold (default 0.20 = 20%).
#' @return A list with: corrected_tgi, uncorrected_tgi, excluded_mice, comparison.
#' @export
weight_corrected_tgi <- function(df,
                                 weight_column    = "Weight",
                                 volume_column    = "Volume",
                                 time_column      = "Day",
                                 treatment_column = "Treatment",
                                 id_column        = "ID",
                                 adjust_tumor_weight = TRUE,
                                 tumor_density    = 1.0,
                                 reference_group  = NULL,
                                endpoint_day     = NULL,
                                endpoint_method  = c("model", "last_obs", "survivors"),
                                 safety_threshold = 0.20) {

  endpoint_method <- match.arg(endpoint_method)

  # --- Validate ---
  required <- c(weight_column, volume_column, time_column, treatment_column, id_column)
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
    Volume    = as.numeric(df[[volume_column]]),
    stringsAsFactors = FALSE
  )

  wd <- wd[!is.na(wd$Weight) & !is.na(wd$Day) & !is.na(wd$Volume), ]
  wd <- wd[order(wd$ID, wd$Day), ]

  groups <- unique(wd$Treatment)
  if (is.null(reference_group)) {
    ctrl_patterns <- c("control", "vehicle", "dmso", "pbs", "saline", "placebo")
    ref_match <- groups[tolower(groups) %in% ctrl_patterns]
    reference_group <- if (length(ref_match) > 0) ref_match[1] else groups[1]
  }

  # --- Compute net weight and baseline ---
  if (adjust_tumor_weight) {
    wd$Net_Weight <- wd$Weight - (wd$Volume / 1000 * tumor_density)
  } else {
    wd$Net_Weight <- wd$Weight
  }

  # Baseline per mouse+treatment at the earliest study day.
  # Aggregate by both ID and Treatment to handle IDs shared across groups.
  min_day <- min(wd$Day, na.rm = TRUE)
  baseline <- stats::aggregate(Net_Weight ~ ID + Treatment,
                               data = wd[wd$Day == min_day, ],
                               FUN = mean, na.rm = TRUE)
  names(baseline)[3] <- "Baseline_Weight"
  wd <- merge(wd, baseline, by = c("ID", "Treatment"), all.x = TRUE)
  wd$Pct_Loss <- (wd$Baseline_Weight - wd$Net_Weight) / wd$Baseline_Weight

  # --- Identify mice exceeding threshold ---
  max_loss <- stats::aggregate(Pct_Loss ~ ID + Treatment, data = wd, FUN = max, na.rm = TRUE)
  max_loss$Exceeded <- max_loss$Pct_Loss >= safety_threshold
  # Use composite key so shared IDs in different groups are excluded independently.
  excluded_keys <- make_mouse_key(
    max_loss$Treatment[max_loss$Exceeded],
    max_loss$ID[max_loss$Exceeded]
  )

  # --- Compute TGI (uncorrected, all mice) ---
  # CODE_REVIEW.md R3.5 / G.3 — endpoint means come from a log-scale LMM over
  # every observation rather than the raw mean among animals still observed at
  # the last day, which conditions on survival and understates efficacy.
  ep_all <- endpoint_volumes(
    wd, id_column = "ID", treatment_column = "Treatment", time_column = "Day",
    volume_column = "Volume", endpoint_day = endpoint_day,
    endpoint_method = endpoint_method
  )
  max_day     <- ep_all$endpoint_day
  uncorrected <- endpoint_tgi(ep_all$group_means, reference_group)

  # --- Compute TGI (corrected, excluding unsafe mice) ---
  safe_data <- wd[
    !make_mouse_key(wd$Treatment, wd$ID) %in% excluded_keys, ]
  ep_safe <- tryCatch(
    endpoint_volumes(safe_data, id_column = "ID",
                     treatment_column = "Treatment", time_column = "Day",
                     volume_column = "Volume", endpoint_day = max_day,
                     endpoint_method = endpoint_method),
    error = function(e) NULL
  )

  if (is.null(ep_safe) || nrow(ep_safe$group_means) == 0) {
    corrected <- uncorrected
    corrected$TGI <- NA_real_
    corrected$N <- 0L
  } else {
    corrected <- tryCatch(endpoint_tgi(ep_safe$group_means, reference_group),
                          error = function(e) NULL)
    if (is.null(corrected)) {
      corrected <- uncorrected
      corrected$TGI <- NA_real_
      corrected$N <- 0L
    }
  }

  # --- Exclusion summary ---
  excluded_summary <- stats::aggregate(Exceeded ~ Treatment, data = max_loss,
                                       FUN = sum)
  names(excluded_summary)[2] <- "N_Excluded"
  total_per_group <- stats::aggregate(Exceeded ~ Treatment, data = max_loss,
                                      FUN = length)
  names(total_per_group)[2] <- "N_Total"
  excluded_summary <- merge(excluded_summary, total_per_group, by = "Treatment")
  excluded_summary$N_Retained <- excluded_summary$N_Total - excluded_summary$N_Excluded

  # --- Comparison table ---
  comparison <- merge(
    uncorrected[, c("Treatment", "TGI", "N")],
    corrected[, c("Treatment", "TGI", "N")],
    by = "Treatment",
    suffixes = c("_Uncorrected", "_Corrected")
  )
  comparison$TGI_Difference <- comparison$TGI_Corrected - comparison$TGI_Uncorrected

  list(
    corrected_tgi   = corrected,
    uncorrected_tgi = uncorrected,
    excluded_mice   = max_loss[max_loss$Exceeded, c("ID", "Treatment", "Pct_Loss")],
    excluded_summary = excluded_summary,
    comparison      = comparison,
    safety_threshold = safety_threshold,
    reference_group  = reference_group
  )
}
