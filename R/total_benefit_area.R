#' Total Benefit Area (Integrated Efficacy-Toxicity)
#'
#' Computes an integrated benefit score per treatment group:
#' B = AUC_efficacy - lambda * AUC_toxicity.
#'
#' @param df Data frame with longitudinal data.
#' @param weight_column Name of the body weight column.
#' @param volume_column Name of the tumor volume column.
#' @param time_column Name of the time/day column.
#' @param treatment_column Name of the treatment group column.
#' @param id_column Name of the mouse/subject ID column.
#' @param adjust_tumor_weight Logical; subtract estimated tumor weight.
#' @param tumor_density Density in g/cm³ (default 1.0).
#' @param reference_group Name of the control/reference group.
#' @param lambda Trade-off parameter (default 1.0). Higher = penalize toxicity more.
#' @return A list with: benefit_table, efficacy_auc, toxicity_auc.
#' @export
total_benefit_area <- function(df,
                               weight_column    = "Weight",
                               volume_column    = "Volume",
                               time_column      = "Day",
                               treatment_column = "Treatment",
                               id_column        = "ID",
                               adjust_tumor_weight = TRUE,
                               tumor_density    = 1.0,
                               reference_group  = NULL,
                               lambda           = 1.0) {

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

  trap_auc <- function(time, value) {
    if (length(time) < 2) return(NA_real_)
    ord <- order(time)
    t <- time[ord]; v <- value[ord]
    sum(diff(t) * (v[-length(v)] + v[-1]) / 2)
  }

  # --- Efficacy AUC: TGI over time per group ---
  # Compute mean volume per group per day
  vol_agg <- stats::aggregate(Volume ~ Treatment + Day, data = wd, FUN = mean, na.rm = TRUE)
  ctrl_vol <- vol_agg[vol_agg$Treatment == reference_group, c("Day", "Volume")]
  names(ctrl_vol)[2] <- "Ctrl_Volume"

  efficacy_auc <- lapply(setdiff(groups, reference_group), function(g) {
    gv <- vol_agg[vol_agg$Treatment == g, c("Day", "Volume")]
    merged <- merge(gv, ctrl_vol, by = "Day")
    merged$TGI <- (1 - merged$Volume / merged$Ctrl_Volume) * 100
    merged$TGI[is.nan(merged$TGI)] <- 0
    auc_val <- trap_auc(merged$Day, merged$TGI)
    data.frame(Treatment = g, Efficacy_AUC = auc_val, stringsAsFactors = FALSE)
  })
  efficacy_df <- do.call(rbind, efficacy_auc)

  # --- Toxicity AUC: % weight loss over time per group ---
  if (adjust_tumor_weight) {
    wd$Net_Weight <- wd$Weight - (wd$Volume / 1000 * tumor_density)
  } else {
    wd$Net_Weight <- wd$Weight
  }

  # Per-group mean baseline
  baseline_grp <- stats::aggregate(Net_Weight ~ Treatment, data = wd,
    FUN = function(x) x[1])
  names(baseline_grp)[2] <- "Baseline_Weight"

  wt_agg <- stats::aggregate(Net_Weight ~ Treatment + Day, data = wd, FUN = mean, na.rm = TRUE)
  wt_agg <- merge(wt_agg, baseline_grp, by = "Treatment")
  wt_agg$Pct_Loss <- pmax((wt_agg$Baseline_Weight - wt_agg$Net_Weight) /
                            wt_agg$Baseline_Weight * 100, 0)

  toxicity_auc <- lapply(setdiff(groups, reference_group), function(g) {
    gw <- wt_agg[wt_agg$Treatment == g, ]
    auc_val <- trap_auc(gw$Day, gw$Pct_Loss)
    data.frame(Treatment = g, Toxicity_AUC = auc_val, stringsAsFactors = FALSE)
  })
  toxicity_df <- do.call(rbind, toxicity_auc)

  # --- Benefit score ---
  benefit <- merge(efficacy_df, toxicity_df, by = "Treatment")
  benefit$Benefit_Score <- benefit$Efficacy_AUC - lambda * benefit$Toxicity_AUC
  benefit <- benefit[order(-benefit$Benefit_Score), ]
  benefit$Rank <- seq_len(nrow(benefit))

  list(
    benefit_table  = benefit,
    efficacy_auc   = efficacy_df,
    toxicity_auc   = toxicity_df,
    lambda         = lambda,
    reference_group = reference_group
  )
}
