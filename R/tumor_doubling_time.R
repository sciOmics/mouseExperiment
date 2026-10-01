#' Calculate Tumor Doubling Time
#'
#' Estimates the time for tumor volume to double based on exponential growth
#' model fitted to individual subject data.
#'
#' @param df Data frame containing tumor measurements.
#' @param time_column Name of the time column. Default "Day".
#' @param volume_column Name of the volume column. Default "Volume".
#' @param treatment_column Name of the treatment column. Default "Treatment".
#' @param id_column Name of the subject ID column. Default "ID".
#' @param cage_column Optional name of the cage column. When supplied, animals are keyed
#'   by ID, treatment and cage, so IDs reused across cages or arms stay distinct.
#'
#' @return A data frame with columns: ID, Treatment, doubling_time, growth_rate, r_squared.
#'
#' @details
#' For each subject, fits \eqn{\log(V) = \beta_0 + \beta_1 \cdot t} via OLS.
#' Doubling time = \eqn{\log(2) / \beta_1}.
#' Subjects with fewer than 3 time points or non-positive growth rates are returned with NA.
#'
#' @examples
#' \dontrun{
#' dt <- tumor_doubling_time(master_synthetic_data)
#' }
#'
#' @export
tumor_doubling_time <- function(df,
                                time_column = "Day",
                                volume_column = "Volume",
                                treatment_column = "Treatment",
                                id_column = "ID",
                                cage_column = NULL) {
  # Validate inputs
  required <- c(time_column, volume_column, treatment_column, id_column)
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  # Composite mouse key so IDs reused across cages/treatments don't collapse.
  cage_col_present <- !is.null(cage_column) && cage_column %in% colnames(df)
  cage_vec <- if (cage_col_present) df[[cage_column]] else rep("1", nrow(df))
  mouse_keys_all <- make_mouse_key(
    df[[id_column]], df[[treatment_column]], cage_vec
  )
  subject_keys <- unique(mouse_keys_all)
  result_list <- vector("list", length(subject_keys))

  for (i in seq_along(subject_keys)) {
    key   <- subject_keys[i]
    sdata <- df[mouse_keys_all == key, , drop = FALSE]
    sdata <- sdata[order(sdata[[time_column]]), ]

    sid       <- sdata[[id_column]][1]
    treatment <- sdata[[treatment_column]][1]
    cage_val  <- if (cage_col_present) sdata[[cage_column]][1] else NA_character_

    # Need positive volumes for log transform, and at least 3 points
    positive <- sdata[[volume_column]] > 0
    sdata_pos <- sdata[positive, ]

    if (nrow(sdata_pos) < 3) {
      result_list[[i]] <- data.frame(
        ID = sid,
        Treatment = treatment,
        Cage = cage_val,
        doubling_time = NA_real_,
        growth_rate = NA_real_,
        r_squared = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    log_vol <- log(sdata_pos[[volume_column]])
    times <- sdata_pos[[time_column]]

    fit <- tryCatch(
      stats::lm(log_vol ~ times),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      growth_rate <- NA_real_
      r_sq <- NA_real_
    } else {
      growth_rate <- stats::coef(fit)[2]
      r_sq <- summary(fit)$r.squared
    }

    dt_val <- if (!is.na(growth_rate) && growth_rate > 0) log(2) / growth_rate else NA_real_

    result_list[[i]] <- data.frame(
      ID = sid,
      Treatment = treatment,
      Cage = cage_val,
      doubling_time = dt_val,
      growth_rate = growth_rate,
      r_squared = r_sq,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, result_list)
}
