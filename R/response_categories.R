#' Response categories (mRECIST) for each animal
#'
#' Classifies each animal's tumour response by the modified RECIST criteria
#' for patient-derived xenografts (Gao et al. 2015, *Nature Medicine* 21:1318):
#' from the % change in volume from baseline, \eqn{\Delta V_t = 100 (V_t -
#' V_0) / V_0}, over the days at least `min_day` after baseline,
#'
#' * **Best response:** the smallest \eqn{\Delta V_t};
#' * **Best average response:** the smallest running mean of
#'   \eqn{\Delta V_0, \ldots, \Delta V_t};
#'
#' and then
#'
#' | Category | Best response | Best average response |
#' |---|---|---|
#' | mCR (complete response) | < -95 % | < -40 % |
#' | mPR (partial response) | < -50 % | < -20 % |
#' | mSD (stable disease) | < 35 % | < 30 % |
#' | mPD (progressive disease) | otherwise | |
#'
#' An animal with no positive baseline volume, or not measured `min_day` or
#' more after baseline, is "Not evaluable". Animals are keyed by arm, ID and
#' cage, so reused ear tags stay apart.
#'
#' @param df Study data, one row per measurement.
#' @param time_column,volume_column,treatment_column,id_column Column names.
#' @param cage_column Cage column, or NULL.
#' @param min_day Days after baseline from which responses count (10 in
#'   Gao et al.); earlier days still enter the running mean.
#' @return A list:
#'   * `per_animal`: one row per animal with `Treatment`, `ID` (and `Cage`),
#'     `Baseline_Day`, `Baseline_Volume`, `Best_Response`, `Best_Response_Day`,
#'     `Best_Avg_Response`, `Category`;
#'   * `summary`: per arm, the number and percentage of animals in each
#'     category;
#'   * `criteria`: the thresholds, as a data frame.
#' @references Gao H, et al. (2015). High-throughput screening using
#'   patient-derived tumor xenografts to predict clinical trial drug response.
#'   *Nature Medicine* 21, 1318-1325.
#' @seealso [plot_response_waterfall()]
#' @export
#' @examples
#' rc <- response_categories(master_synthetic_data, cage_column = "Cage")
#' rc$summary
response_categories <- function(df, time_column = "Day", volume_column = "Volume",
                                treatment_column = "Treatment", id_column = "ID",
                                cage_column = NULL, min_day = 10) {
  for (col in c(time_column, volume_column, treatment_column, id_column)) {
    if (!has_col(df, col)) stop("Column '", col, "' is not in the data.", call. = FALSE)
  }
  has_cage <- has_col(df, cage_column)
  d <- data.frame(Treatment = as.character(df[[treatment_column]]),
                  ID = as.character(df[[id_column]]),
                  Cage = if (has_cage) as.character(df[[cage_column]]) else NA_character_,
                  Day = suppressWarnings(as.numeric(df[[time_column]])),
                  Volume = suppressWarnings(as.numeric(df[[volume_column]])),
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$Treatment) & is.finite(d$Day) & is.finite(d$Volume), , drop = FALSE]
  d$Key <- if (has_cage) make_mouse_key(d$Treatment, d$ID, d$Cage) else make_mouse_key(d$Treatment, d$ID)
  d <- d[order(d$Key, d$Day), , drop = FALSE]

  per_animal <- do.call(rbind, lapply(split(d, d$Key), function(a) {
    v0 <- a$Volume[1]; t0 <- a$Day[1]
    out <- data.frame(Treatment = a$Treatment[1], ID = a$ID[1], Cage = a$Cage[1],
                      Baseline_Day = t0, Baseline_Volume = v0,
                      Best_Response = NA_real_, Best_Response_Day = NA_real_,
                      Best_Avg_Response = NA_real_, Category = "Not evaluable",
                      stringsAsFactors = FALSE)
    if (!(v0 > 0)) return(out)
    dv <- 100 * (a$Volume - v0) / v0
    run_mean <- cumsum(dv) / seq_along(dv)
    late <- a$Day - t0 >= min_day
    if (!any(late)) return(out)
    i <- which(late)[which.min(dv[late])]
    out$Best_Response     <- dv[i]
    out$Best_Response_Day <- a$Day[i]
    out$Best_Avg_Response <- min(run_mean[late])
    out$Category <- mrecist_category(out$Best_Response, out$Best_Avg_Response)
    out
  }))
  rownames(per_animal) <- NULL
  if (!has_cage) per_animal$Cage <- NULL
  per_animal <- per_animal[order(per_animal$Treatment, per_animal$ID), , drop = FALSE]
  rownames(per_animal) <- NULL

  cats <- c("mCR", "mPR", "mSD", "mPD", "Not evaluable")
  tab <- table(factor(per_animal$Treatment, levels = unique(per_animal$Treatment)),
               factor(per_animal$Category, levels = cats))
  summary <- data.frame(Treatment = rownames(tab), N = as.integer(rowSums(tab)),
                        stringsAsFactors = FALSE)
  for (k in cats) {
    summary[[k]] <- as.integer(tab[, k])
  }
  n_eval <- summary$N - summary[["Not evaluable"]]
  summary$Response_Rate <- ifelse(n_eval > 0, 100 * (summary$mCR + summary$mPR) / n_eval, NA_real_)

  me_result(list(per_animal = per_animal, summary = summary,
                 criteria = data.frame(Category = c("mCR", "mPR", "mSD"),
                                       Best_Response_Below = c(-95, -50, 35),
                                       Best_Avg_Response_Below = c(-40, -20, 30)),
                 min_day = min_day), "me_responses")
}

#' The mRECIST category for best and best-average responses
#' @noRd
mrecist_category <- function(best, best_avg) {
  ifelse(best < -95 & best_avg < -40, "mCR",
  ifelse(best < -50 & best_avg < -20, "mPR",
  ifelse(best < 35 & best_avg < 30, "mSD", "mPD")))
}
