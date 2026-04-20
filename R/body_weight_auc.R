#' Body Weight AUC Analysis
#'
#' Computes trapezoidal area under the curve for body weight change over time,
#' with per-mouse and per-group summaries.
#'
#' @param df Data frame with longitudinal data.
#' @param weight_column Name of the body weight column.
#' @param time_column Name of the time/day column.
#' @param treatment_column Name of the treatment group column.
#' @param id_column Name of the mouse/subject ID column.
#' @param volume_column Name of the tumor volume column. NULL to skip tumor adjustment.
#' @param adjust_tumor_weight Logical; subtract estimated tumor weight.
#' @param tumor_density Density in g/cm³ (default 1.0).
#' @param reference_group Name of the control/reference group.
#' @return A list with: auc_per_mouse, auc_summary, comparisons, nadir_data.
#' @export
body_weight_auc <- function(df,
                            weight_column    = "Weight",
                            time_column      = "Day",
                            treatment_column = "Treatment",
                            id_column        = "ID",
                            volume_column    = NULL,
                            adjust_tumor_weight = TRUE,
                            tumor_density    = 1.0,
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

  # --- Per-mouse baseline and % change ---
  baseline <- stats::aggregate(Weight ~ ID, data = wd,
                               FUN = function(x) x[1])
  names(baseline)[2] <- "Baseline_Weight"
  wd <- merge(wd, baseline, by = "ID", all.x = TRUE)
  wd$Pct_Change <- (wd$Weight - wd$Baseline_Weight) / wd$Baseline_Weight * 100

  # --- Trapezoidal AUC per mouse ---
  trap_auc <- function(time, value) {
    if (length(time) < 2) return(NA_real_)
    ord <- order(time)
    t <- time[ord]
    v <- value[ord]
    sum(diff(t) * (v[-length(v)] + v[-1]) / 2)
  }

  ids <- unique(wd$ID)
  auc_list <- lapply(ids, function(id) {
    sub <- wd[wd$ID == id, ]
    data.frame(
      ID        = id,
      Treatment = sub$Treatment[1],
      Baseline_Weight = sub$Baseline_Weight[1],
      AUC_Weight      = trap_auc(sub$Day, sub$Weight),
      AUC_Pct_Change  = trap_auc(sub$Day, sub$Pct_Change),
      Nadir_Weight    = min(sub$Weight, na.rm = TRUE),
      Nadir_Day       = sub$Day[which.min(sub$Weight)],
      Nadir_Pct_Change = min(sub$Pct_Change, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  auc_df <- do.call(rbind, auc_list)

  # --- Group summaries ---
  auc_summary <- stats::aggregate(
    cbind(AUC_Pct_Change, Nadir_Pct_Change) ~ Treatment,
    data = auc_df,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                        sd   = stats::sd(x, na.rm = TRUE),
                        n    = length(x))
  )
  # Flatten the matrix columns
  flat <- do.call(data.frame, auc_summary)
  names(flat) <- c("Treatment",
                    "AUC_Mean", "AUC_SD", "AUC_N",
                    "Nadir_Mean", "Nadir_SD", "Nadir_N")
  flat$AUC_SEM   <- flat$AUC_SD / sqrt(flat$AUC_N)
  flat$Nadir_SEM <- flat$Nadir_SD / sqrt(flat$Nadir_N)

  # --- Pairwise comparisons (Welch t-test on AUC) ---
  groups <- unique(auc_df$Treatment)
  comparisons <- NULL
  if (length(groups) >= 2) {
    pairs <- utils::combn(groups, 2, simplify = FALSE)
    comp_list <- lapply(pairs, function(pair) {
      a <- auc_df$AUC_Pct_Change[auc_df$Treatment == pair[1]]
      b <- auc_df$AUC_Pct_Change[auc_df$Treatment == pair[2]]
      tt <- tryCatch(stats::t.test(a, b), error = function(e) NULL)
      if (is.null(tt)) return(NULL)
      data.frame(
        Group1   = pair[1],
        Group2   = pair[2],
        Mean_Diff = tt$estimate[1] - tt$estimate[2],
        P_Value  = tt$p.value,
        CI_Lower = tt$conf.int[1],
        CI_Upper = tt$conf.int[2],
        stringsAsFactors = FALSE
      )
    })
    comparisons <- do.call(rbind, comp_list[!vapply(comp_list, is.null, logical(1))])
  }

  list(
    auc_per_mouse = auc_df,
    auc_summary   = flat,
    comparisons   = comparisons,
    nadir_data    = auc_df[, c("ID", "Treatment", "Nadir_Weight",
                               "Nadir_Day", "Nadir_Pct_Change")],
    weight_data   = wd
  )
}
