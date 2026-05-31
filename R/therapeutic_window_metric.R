#' Therapeutic Window Metric (TWM)
#'
#' Computes TWM = TGI / MeanWeightLoss% per treatment group.
#' When weight loss is negligible (≤ noise floor), the safety score equals TGI.
#' The denominator uses the group mean of per-mouse maximum weight loss, not the
#' single worst-case mouse, to give a representative measure of typical toxicity.
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
#' @param noise_floor Minimum group-mean per-mouse maximum weight loss (percent)
#'   below which the ratio TGI / weight_loss\% is numerically unstable. When
#'   below this floor, \code{TWM = abs(TGI)} (pure efficacy score) instead of
#'   the ratio. Default 1.0 — an experimentally pragmatic threshold chosen to
#'   avoid division by near-zero weight loss in well-tolerated treatments. No
#'   formal clinical basis; users with experiment-specific noise estimates
#'   (e.g. scale precision) should tune accordingly.
#' @return A list with: twm_table, tgi_data, weight_loss_data.
#' @export
therapeutic_window_metric <- function(df,
                                      weight_column    = "Weight",
                                      volume_column    = "Volume",
                                      time_column      = "Day",
                                      treatment_column = "Treatment",
                                      id_column        = "ID",
                                      cage_column      = NULL,
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
  cage_col_present <- !is.null(cage_column) && cage_column %in% names(df)
  cage_vec <- if (cage_col_present) as.character(df[[cage_column]]) else "1"
  wd <- data.frame(
    ID        = as.character(df[[id_column]]),
    Treatment = as.character(df[[treatment_column]]),
    Cage      = cage_vec,
    Day       = as.numeric(df[[time_column]]),
    Weight    = as.numeric(df[[weight_column]]),
    Volume    = as.numeric(df[[volume_column]]),
    stringsAsFactors = FALSE
  )
  # Composite mouse key prevents collapsing same-numeric-ID mice across cages
  wd$MouseKey <- make_mouse_key(wd$ID, wd$Treatment, wd$Cage)

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
  # Filter to the earliest study day before aggregating so x[1] is ordered.
  # Aggregate by MouseKey (ID|||Treatment|||Cage) so reused IDs across cages
  # don't collapse — same fix class as Round 1 1.8 for weight_corrected_tgi.
  min_day <- min(wd$Day, na.rm = TRUE)
  baseline <- stats::aggregate(Weight ~ MouseKey + Treatment,
                               data = wd[wd$Day == min_day, ],
                               FUN = mean, na.rm = TRUE)
  names(baseline)[3] <- "Baseline_Weight"

  nadir <- stats::aggregate(Weight ~ MouseKey + Treatment, data = wd,
                            FUN = min, na.rm = TRUE)
  names(nadir)[3] <- "Nadir_Weight"

  mouse_wl <- merge(baseline, nadir, by = c("MouseKey", "Treatment"))
  mouse_wl$Pct_Loss <- (mouse_wl$Baseline_Weight - mouse_wl$Nadir_Weight) /
                        mouse_wl$Baseline_Weight * 100
  mouse_wl$Pct_Loss <- pmax(mouse_wl$Pct_Loss, 0)  # clamp negative (weight gain)

  # Group mean of per-mouse max weight loss: more representative of typical
  # toxicity than the single worst-case mouse (previously used max here).
  group_wl <- stats::aggregate(Pct_Loss ~ Treatment, data = mouse_wl,
                               FUN = mean, na.rm = TRUE)
  names(group_wl)[2] <- "Mean_Pct_Weight_Loss"

  # --- TWM ---
  twm <- merge(tgi_data, group_wl, by = "Treatment")
  twm$TWM <- ifelse(
    twm$Mean_Pct_Weight_Loss <= noise_floor,
    abs(twm$TGI),  # Safety score = TGI when weight loss negligible
    abs(twm$TGI) / twm$Mean_Pct_Weight_Loss
  )
  twm$Safety_Note <- ifelse(
    twm$Mean_Pct_Weight_Loss <= noise_floor,
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
