#' Handle necrotic tumor observations
#'
#' Internal helper called by \code{tumor_growth_statistics()} and
#' \code{bayesian_tumor_growth()} when a necrotic flag column is present.
#' Normalises the flag, computes a per-group summary, and optionally removes
#' necrotic rows or adds a covariate column before model fitting.
#'
#' @param df Data frame (untransformed; after extrapolation).
#' @param necrotic_column Name of the column containing the necrotic flag.
#'   Recognised values: "Y"/"N", "yes"/"no", "1"/"0", "TRUE"/"FALSE"
#'   (case-insensitive).
#' @param necrotic_handling One of \code{"exclude"}, \code{"covariate"}, or
#'   \code{"none"}.
#' @param treatment_column Name of the treatment group column.
#' @param id_column Name of the animal ID column.
#' @param time_column Name of the time/day column.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{data}}{Processed data frame. When \code{necrotic_handling =
#'       "exclude"}, necrotic rows are removed. When \code{"covariate"}, a
#'       column \code{necrotic_cov_flag} (integer 0/1) is appended and retained
#'       for the formula builder.}
#'     \item{\code{necrosis_summary}}{Data frame with one row per treatment
#'       group: Treatment, N_Observations, N_Necrotic, Pct_Necrotic,
#'       N_Mice_Total, N_Mice_With_Necrosis, First_Necrotic_Day.}
#'   }
#' @keywords internal
tgs_handle_necrosis <- function(df, necrotic_column, necrotic_handling,
                                treatment_column, id_column, time_column) {
  raw        <- df[[necrotic_column]]
  is_necrotic <- tolower(as.character(raw)) %in% c("y", "yes", "1", "true")
  df$necrotic_cov_flag <- as.integer(is_necrotic)

  # Per-group summary computed on the full data (before any exclusion)
  groups <- split(seq_len(nrow(df)), df[[treatment_column]])
  necrosis_summary <- do.call(rbind, lapply(groups, function(idx) {
    sub_df  <- df[idx, , drop = FALSE]
    nec_sub <- is_necrotic[idx]
    data.frame(
      Treatment            = sub_df[[treatment_column]][1],
      N_Observations       = nrow(sub_df),
      N_Necrotic           = sum(nec_sub),
      Pct_Necrotic         = round(100 * mean(nec_sub), 1),
      N_Mice_Total         = length(unique(sub_df[[id_column]])),
      N_Mice_With_Necrosis = length(unique(sub_df[[id_column]][nec_sub])),
      First_Necrotic_Day   = if (any(nec_sub))
                               min(sub_df[[time_column]][nec_sub], na.rm = TRUE)
                             else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(necrosis_summary) <- NULL

  if (necrotic_handling == "exclude") {
    n_out <- sum(is_necrotic)
    df    <- df[!is_necrotic, , drop = FALSE]
    attr(df, ".n_necrotic_excluded") <- n_out
  }
  # "covariate": necrotic_cov_flag stays in df; formula builder appends it
  # "none":      necrotic_cov_flag is present but not used in the formula

  list(data = df, necrosis_summary = necrosis_summary)
}
