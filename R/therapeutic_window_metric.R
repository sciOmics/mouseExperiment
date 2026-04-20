#' Therapeutic Window Metric (TWM)
#'
#' Computes TWM = TGI / MaxWeightLoss% per treatment group.
#' When weight loss is negligible (≤ noise floor), the safety score equals TGI.
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
#' @param noise_floor Minimum weight loss % below which TWM = TGI (default 1.0).
#' @return A list with: twm_table, tgi_data, weight_loss_data.
#' @export
therapeutic_window_metric <- function(df,
                                      weight_column    = "Weight",
                                      volume_column    = "Volume",
                                      time_column      = "Day",
                                      treatment_column = "Treatment",
                                      id_column        = "ID",
                                      adjust_tumor_weight = TRUE,
                                      tumor_density    = 1.0,
                                      reference_group  = NULL,
                                      noise_floor      = 1.0) {

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

  if (adjust_tumor_weight) {
    wd$Weight <- wd$Weight - (wd$Volume / 1000 * tumor_density)
  }

  wd <- wd[!is.na(wd$Weight) & !is.na(wd$Day) & !is.na(wd$Volume), ]
  wd <- wd[order(wd$ID, wd$Day), ]

  groups <- unique(wd$Treatment)
  if (is.null(reference_group)) {
    # Pick first alphabetically or common control names
    ctrl_patterns <- c("control", "vehicle", "dmso", "pbs", "saline", "placebo")
    ref_match <- groups[tolower(groups) %in% ctrl_patterns]
    reference_group <- if (length(ref_match) > 0) ref_match[1] else groups[1]
  }

  # --- TGI per treatment group ---
  # Final timepoint mean volume
  max_day <- max(wd$Day, na.rm = TRUE)
  final <- wd[wd$Day == max_day, ]
  ctrl_mean_vol <- mean(final$Volume[final$Treatment == reference_group], na.rm = TRUE)

  tgi_data <- stats::aggregate(Volume ~ Treatment, data = final, FUN = mean)
  names(tgi_data)[2] <- "Mean_Volume"
  tgi_data$TGI <- (1 - tgi_data$Mean_Volume / ctrl_mean_vol) * 100
  tgi_data$TGI[tgi_data$Treatment == reference_group] <- 0

  # --- Max % weight loss per group ---
  # Per mouse: baseline weight, nadir weight, max % loss
  baseline <- stats::aggregate(Weight ~ ID + Treatment, data = wd,
                               FUN = function(x) x[1])
  names(baseline)[3] <- "Baseline_Weight"

  nadir <- stats::aggregate(Weight ~ ID + Treatment, data = wd,
                            FUN = min, na.rm = TRUE)
  names(nadir)[3] <- "Nadir_Weight"

  mouse_wl <- merge(baseline, nadir, by = c("ID", "Treatment"))
  mouse_wl$Pct_Loss <- (mouse_wl$Baseline_Weight - mouse_wl$Nadir_Weight) /
                        mouse_wl$Baseline_Weight * 100
  mouse_wl$Pct_Loss <- pmax(mouse_wl$Pct_Loss, 0)  # clamp negative (weight gain)

  group_wl <- stats::aggregate(Pct_Loss ~ Treatment, data = mouse_wl,
                               FUN = max, na.rm = TRUE)
  names(group_wl)[2] <- "Max_Pct_Weight_Loss"

  # --- TWM ---
  twm <- merge(tgi_data, group_wl, by = "Treatment")
  twm$TWM <- ifelse(
    twm$Max_Pct_Weight_Loss <= noise_floor,
    abs(twm$TGI),  # Safety score = TGI when weight loss negligible
    abs(twm$TGI) / twm$Max_Pct_Weight_Loss
  )
  twm$Safety_Note <- ifelse(
    twm$Max_Pct_Weight_Loss <= noise_floor,
    "Negligible weight loss",
    ""
  )
  twm <- twm[order(-twm$TWM), ]

  list(
    twm_table        = twm,
    tgi_data         = tgi_data,
    weight_loss_data = mouse_wl
  )
}
