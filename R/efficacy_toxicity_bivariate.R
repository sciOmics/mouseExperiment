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
#' @param endpoint_day Day at which efficacy is evaluated. NULL uses the
#'   maximum observed day.
#' @param endpoint_method How the TGI denominator is obtained: "model"
#'   (default), "last_obs", or "survivors". See CODE_REVIEW.md R3.5 / G.3.
#' @param reference_group Name of the control/reference group.
#' @param efficacy_metric One of "tgi" or "tumor_auc".
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
                                        endpoint_day     = NULL,
                                        endpoint_method  = c("model", "last_obs", "survivors"),
                                        efficacy_metric  = c("tgi", "tumor_auc")) {

  endpoint_method <- match.arg(endpoint_method)

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

  # Composite key prevents collapsing reused IDs across treatments.
  # The efficacy arm below already uses make_mouse_key — bring toxicity inline.
  # Also filter to the earliest study day so x[1] is unambiguous (Round 1 1.7).
  wt$.MouseKey <- make_mouse_key(wt$Treatment, wt$ID)
  min_day_w <- min(wt$Day, na.rm = TRUE)
  baseline_w <- stats::aggregate(
    Net_Weight ~ .MouseKey,
    data = wt[wt$Day == min_day_w, ],
    FUN  = mean, na.rm = TRUE
  )
  names(baseline_w)[2] <- "Baseline_Weight"
  wt <- merge(wt, baseline_w, by = ".MouseKey", all.x = TRUE)
  wt$Pct_Weight_Loss <- (wt$Baseline_Weight - wt$Net_Weight) / wt$Baseline_Weight * 100

  max_wl <- stats::aggregate(
    Pct_Weight_Loss ~ .MouseKey + Treatment,
    data = wt, FUN = max, na.rm = TRUE
  )
  names(max_wl)[c(1, 3)] <- c("MouseKey", "Max_Pct_Weight_Loss")
  # Recover ID for display from the composite key (Treatment|||ID format)
  max_wl$ID <- vapply(strsplit(max_wl$MouseKey, "\\|\\|\\|"),
                      function(p) p[2], character(1L))
  max_wl <- max_wl[, c("ID", "Treatment", "Max_Pct_Weight_Loss")]
  max_wl$Max_Pct_Weight_Loss <- pmax(max_wl$Max_Pct_Weight_Loss, 0)

  # --- Efficacy: per mouse ---
  # CODE_REVIEW.md R3.5 / G.3 — the TGI denominator was the raw control mean
  # among animals still observed at the global last day, which conditions on
  # survival and biases TGI downward. Take it from the shared endpoint helper
  # instead, which by default reads the arm's geometric mean at the endpoint day
  # off a log-scale LMM fitted to every observation.
  ep <- endpoint_volumes(
    wd, id_column = "ID", treatment_column = "Treatment", time_column = "Day",
    volume_column = "Volume", endpoint_day = endpoint_day,
    endpoint_method = endpoint_method
  )
  max_day <- ep$endpoint_day
  ctrl_data <- wd[wd$Treatment == reference_group, ]
  ctrl_mean_final_vol <- ep$group_means$Mean_Volume[
    ep$group_means$Treatment == reference_group]
  if (length(ctrl_mean_final_vol) != 1L || !is.finite(ctrl_mean_final_vol) ||
      ctrl_mean_final_vol <= 0) {
    stop("Reference group '", reference_group, "' has no usable endpoint volume.",
         call. = FALSE)
  }

  # Control AUC: mean of per-mouse AUCs (pooling all observations produces a
  # nonsensical integral when multiple mice share the same timepoints).
  ctrl_mouse_keys <- unique(make_mouse_key(ctrl_data$Treatment, ctrl_data$ID))
  ctrl_aucs <- vapply(ctrl_mouse_keys, function(k) {
    s <- ctrl_data[make_mouse_key(ctrl_data$Treatment, ctrl_data$ID) == k, ]
    calculate_auc(s$Day, s$Volume)
  }, numeric(1))
  ctrl_mean_auc <- mean(ctrl_aucs, na.rm = TRUE)

  # Per-mouse efficacy — use composite key so IDs shared across groups are
  # treated as distinct mice (matches the pattern in body_weight_auc.R)
  wd$.MouseKey <- make_mouse_key(wd$Treatment, wd$ID)
  mouse_keys <- unique(wd$.MouseKey)
  eff_list <- lapply(mouse_keys, function(key) {
    sub <- wd[wd$.MouseKey == key, ]
    tx <- sub$Treatment[1]
    id <- sub$ID[1]
    # R3.5 — an animal removed before max_day has no row at max_day, and the
    # old `sub$Volume[sub$Day == max_day]` yielded numeric(0) for it. Use its
    # last observation at or before the endpoint day so it still contributes.
    in_window <- sub[sub$Day <= max_day, , drop = FALSE]
    in_window <- in_window[order(in_window$Day), ]
    final_vol <- if (nrow(in_window) > 0) {
      in_window$Volume[nrow(in_window)]
    } else sub$Volume[nrow(sub)]

    efficacy_val <- switch(efficacy_metric,
      tgi = {
        (1 - final_vol[1] / ctrl_mean_final_vol) * 100
      },
      tumor_auc = {
        # Lower AUC = better efficacy; we invert so higher = better
        mouse_auc <- calculate_auc(sub$Day, sub$Volume)
        if (!is.na(ctrl_mean_auc) && ctrl_mean_auc > 0) {
          (1 - mouse_auc / ctrl_mean_auc) * 100
        } else NA_real_
      }
    )

    data.frame(ID = id, Treatment = tx, Efficacy = efficacy_val,
               stringsAsFactors = FALSE)
  })
  eff_df <- do.call(rbind, eff_list)
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
    reference_group = reference_group
  )
}
