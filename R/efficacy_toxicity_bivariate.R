#' Efficacy-Toxicity Bivariate Analysis
#'
#' Computes per-mouse and per-group efficacy and toxicity metrics for
#' safety-efficacy scatter plots. Supports multiple efficacy metrics:
#' Final TGI, AUC of Tumor Volume, and Log-Cell Kill.
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
#' @param efficacy_metric One of "tgi", "tumor_auc", or "log_cell_kill".
#' @return A list with: per_mouse, per_group, efficacy_metric, reference_group.
#' @export
efficacy_toxicity_bivariate <- function(df,
                                        weight_column    = "Weight",
                                        volume_column    = "Volume",
                                        time_column      = "Day",
                                        treatment_column = "Treatment",
                                        id_column        = "ID",
                                        adjust_tumor_weight = TRUE,
                                        tumor_density    = 1.0,
                                        reference_group  = NULL,
                                        efficacy_metric  = c("tgi", "tumor_auc", "log_cell_kill")) {

  efficacy_metric <- match.arg(efficacy_metric)

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

  # --- Toxicity: max % weight loss per mouse ---
  wt <- wd
  if (adjust_tumor_weight) {
    wt$Net_Weight <- wt$Weight - (wt$Volume / 1000 * tumor_density)
  } else {
    wt$Net_Weight <- wt$Weight
  }

  baseline_w <- stats::aggregate(Net_Weight ~ ID, data = wt, FUN = function(x) x[1])
  names(baseline_w)[2] <- "Baseline_Weight"
  wt <- merge(wt, baseline_w, by = "ID", all.x = TRUE)
  wt$Pct_Weight_Loss <- (wt$Baseline_Weight - wt$Net_Weight) / wt$Baseline_Weight * 100

  max_wl <- stats::aggregate(Pct_Weight_Loss ~ ID + Treatment, data = wt, FUN = max, na.rm = TRUE)
  names(max_wl)[3] <- "Max_Pct_Weight_Loss"
  max_wl$Max_Pct_Weight_Loss <- pmax(max_wl$Max_Pct_Weight_Loss, 0)

  # --- Efficacy: per mouse ---
  trap_auc <- function(time, value) {
    if (length(time) < 2) return(NA_real_)
    ord <- order(time)
    t <- time[ord]; v <- value[ord]
    sum(diff(t) * (v[-length(v)] + v[-1]) / 2)
  }

  max_day <- max(wd$Day, na.rm = TRUE)
  ctrl_data <- wd[wd$Treatment == reference_group, ]
  ctrl_mean_final_vol <- mean(ctrl_data$Volume[ctrl_data$Day == max_day], na.rm = TRUE)

  # Control group doubling time for LCK
  ctrl_doubling_time <- NA_real_
  if (efficacy_metric == "log_cell_kill") {
    ctrl_growth <- tryCatch({
      ctrl_agg <- stats::aggregate(Volume ~ Day, data = ctrl_data, FUN = mean)
      ctrl_agg <- ctrl_agg[order(ctrl_agg$Day), ]
      if (nrow(ctrl_agg) >= 2 && all(ctrl_agg$Volume > 0)) {
        lm_fit <- stats::lm(log(Volume) ~ Day, data = ctrl_agg)
        growth_rate <- stats::coef(lm_fit)["Day"]
        if (!is.na(growth_rate) && growth_rate > 0) {
          log(2) / growth_rate
        } else NA_real_
      } else NA_real_
    }, error = function(e) NA_real_)
    ctrl_doubling_time <- ctrl_growth
  }

  # Per-mouse efficacy
  ids <- unique(wd$ID)
  eff_list <- lapply(ids, function(id) {
    sub <- wd[wd$ID == id, ]
    tx <- sub$Treatment[1]
    final_vol <- sub$Volume[sub$Day == max_day]
    if (length(final_vol) == 0) final_vol <- sub$Volume[nrow(sub)]

    efficacy_val <- switch(efficacy_metric,
      tgi = {
        (1 - final_vol[1] / ctrl_mean_final_vol) * 100
      },
      tumor_auc = {
        # Lower AUC = better efficacy; we invert so higher = better
        ctrl_auc <- trap_auc(ctrl_data$Day, ctrl_data$Volume)
        mouse_auc <- trap_auc(sub$Day, sub$Volume)
        if (!is.na(ctrl_auc) && ctrl_auc > 0) {
          (1 - mouse_auc / ctrl_auc) * 100
        } else NA_real_
      },
      log_cell_kill = {
        if (is.na(ctrl_doubling_time) || ctrl_doubling_time <= 0) {
          NA_real_
        } else {
          # Growth delay version: LCK = (T - C) / (3.32 × Td)
          # Use time for treated to reach a target volume vs control
          target_vol <- ctrl_mean_final_vol
          ctrl_time <- tryCatch({
            ctrl_agg <- stats::aggregate(Volume ~ Day, data = ctrl_data, FUN = mean)
            ctrl_agg <- ctrl_agg[order(ctrl_agg$Day), ]
            hit <- ctrl_agg$Day[ctrl_agg$Volume >= target_vol]
            if (length(hit) > 0) hit[1] else max(ctrl_agg$Day)
          }, error = function(e) NA_real_)

          treated_time <- tryCatch({
            hit <- sub$Day[sub$Volume >= target_vol]
            if (length(hit) > 0) hit[1] else max(sub$Day)
          }, error = function(e) NA_real_)

          if (!is.na(ctrl_time) && !is.na(treated_time)) {
            (treated_time - ctrl_time) / (3.32 * ctrl_doubling_time)
          } else NA_real_
        }
      }
    )

    data.frame(ID = id, Treatment = tx, Efficacy = efficacy_val,
               stringsAsFactors = FALSE)
  })
  eff_df <- do.call(rbind, eff_list)

  # --- Merge per-mouse ---
  per_mouse <- merge(max_wl, eff_df, by = c("ID", "Treatment"))

  # --- Per-group means ---
  per_group <- stats::aggregate(
    cbind(Max_Pct_Weight_Loss, Efficacy) ~ Treatment,
    data = per_mouse,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE))
  )
  pg <- data.frame(
    Treatment        = per_group$Treatment,
    Toxicity_Mean    = per_group$Max_Pct_Weight_Loss[, "mean"],
    Toxicity_SD      = per_group$Max_Pct_Weight_Loss[, "sd"],
    Efficacy_Mean    = per_group$Efficacy[, "mean"],
    Efficacy_SD      = per_group$Efficacy[, "sd"],
    stringsAsFactors = FALSE
  )

  list(
    per_mouse       = per_mouse,
    per_group       = pg,
    efficacy_metric = efficacy_metric,
    reference_group = reference_group,
    ctrl_doubling_time = ctrl_doubling_time
  )
}
