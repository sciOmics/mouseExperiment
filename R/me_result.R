#' S3 Class for mouseExperiment Analysis Results
#'
#' @description
#' The \code{me_result} S3 class provides a consistent wrapper around analysis
#' outputs from mouseExperiment functions.
#'
#' @details
#' Every analysis function returns an \code{me_result} object with at least:
#' \itemize{
#'   \item \code{analysis_type}: character label (e.g. "tumor_growth", "survival", "synergy")
#'   \item \code{data}: the input data used
#'   \item \code{results}: a named list of analysis outputs
#'   \item \code{plots}: a named list of ggplot objects (may be NULL)
#'   \item \code{summary}: a \code{data.frame} of key findings
#'   \item \code{call}: the matched call
#'   \item \code{timestamp}: when analysis was run
#' }
#'
#' @param x An \code{me_result} object.
#' @param ... Additional arguments passed to methods.
#' @name me_result
NULL

#' Create an me_result object
#'
#' @param analysis_type Character label for the analysis type.
#' @param data The input data frame.
#' @param results Named list of analysis outputs.
#' @param plots Named list of ggplot objects (default NULL).
#' @param summary_df Data frame of key summary findings (default NULL).
#' @param call The matched call (default NULL).
#' @param ... Additional named elements to include.
#'
#' @return An S3 object of class \code{me_result}.
#' @export
new_me_result <- function(analysis_type,
                          data = NULL,
                          results = list(),
                          plots = NULL,
                          summary_df = NULL,
                          call = NULL,
                          ...) {
  structure(
    list(
      analysis_type = analysis_type,
      data = data,
      results = results,
      plots = plots,
      summary = summary_df,
      call = call,
      timestamp = Sys.time(),
      ...
    ),
    class = "me_result"
  )
}

#' @rdname me_result
#' @export
print.me_result <- function(x, ...) {
  cat("mouseExperiment Analysis Result\n")
  cat("  Type      :", x$analysis_type, "\n")
  cat("  Timestamp :", format(x$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
  if (!is.null(x$data)) {
    cat("  Data      :", nrow(x$data), "rows x", ncol(x$data), "cols\n")
  }
  if (!is.null(x$summary) && is.data.frame(x$summary)) {
    cat("  Summary   :", nrow(x$summary), "rows\n")
  }
  if (!is.null(x$plots)) {
    cat("  Plots     :", length(x$plots), "\n")
  }
  cat("  Elements  :", paste(names(x$results), collapse = ", "), "\n")
  invisible(x)
}

#' @rdname me_result
#' @export
summary.me_result <- function(object, ...) {
  if (!is.null(object$summary) && is.data.frame(object$summary)) {
    object$summary
  } else {
    cat("No summary data frame available.\n")
    invisible(NULL)
  }
}

#' Extract plots from an me_result
#'
#' @param x An \code{me_result} object.
#' @param which Character name of the plot to extract, or NULL for all plots.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object or named list of ggplot objects.
#' @export
plot.me_result <- function(x, which = NULL, ...) {
  if (is.null(x$plots) || length(x$plots) == 0) {
    message("No plots available in this result.")
    return(invisible(NULL))
  }
  if (!is.null(which)) {
    if (which %in% names(x$plots)) {
      return(x$plots[[which]])
    } else {
      stop("Plot '", which, "' not found. Available: ",
           paste(names(x$plots), collapse = ", "), call. = FALSE)
    }
  }
  x$plots
}

#' Export model diagnostics from an me_result
#'
#' @param x An \code{me_result} object.
#' @param file Path to write diagnostics CSV. If NULL, returns data frame.
#'
#' @return A data frame of diagnostic metrics (invisibly if written to file).
#' @export
export_diagnostics <- function(x, file = NULL) {
  if (!inherits(x, "me_result")) {
    stop("x must be an me_result object", call. = FALSE)
  }
  
  diag_df <- data.frame(
    metric = character(0),
    value = numeric(0),
    stringsAsFactors = FALSE
  )
  
  results <- x$results
  
  # Extract model-level diagnostics if present
  if (!is.null(results$model)) {
    model <- results$model
    if (inherits(model, c("lm", "lmerMod"))) {
      diag_df <- rbind(diag_df, data.frame(
        metric = c("AIC", "BIC", "logLik", "n_obs"),
        value = c(stats::AIC(model), stats::BIC(model),
                  as.numeric(stats::logLik(model)), stats::nobs(model)),
        stringsAsFactors = FALSE
      ))
    }
    if (inherits(model, "lm") && !inherits(model, "lmerMod")) {
      s <- summary(model)
      diag_df <- rbind(diag_df, data.frame(
        metric = c("R_squared", "Adj_R_squared", "F_statistic"),
        value = c(s$r.squared, s$adj.r.squared, s$fstatistic[1]),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Extract ANOVA p-values if present
  if (!is.null(results$anova_p)) {
    diag_df <- rbind(diag_df, data.frame(
      metric = "anova_p_value",
      value = results$anova_p,
      stringsAsFactors = FALSE
    ))
  }
  
  if (!is.null(file)) {
    utils::write.csv(diag_df, file = file, row.names = FALSE)
    message("Diagnostics written to: ", file)
    return(invisible(diag_df))
  }
  diag_df
}

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
#' dt <- tumor_doubling_time(my_data)
#' }
#'
#' @export
tumor_doubling_time <- function(df,
                                time_column = "Day",
                                volume_column = "Volume",
                                treatment_column = "Treatment",
                                id_column = "ID") {
  # Validate inputs
  required <- c(time_column, volume_column, treatment_column, id_column)
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  
  subjects <- unique(df[[id_column]])
  result_list <- vector("list", length(subjects))
  
  for (i in seq_along(subjects)) {
    sid <- subjects[i]
    sdata <- df[df[[id_column]] == sid, ]
    sdata <- sdata[order(sdata[[time_column]]), ]
    
    treatment <- sdata[[treatment_column]][1]
    
    # Need positive volumes for log transform, and at least 3 points
    positive <- sdata[[volume_column]] > 0
    sdata_pos <- sdata[positive, ]
    
    if (nrow(sdata_pos) < 3) {
      result_list[[i]] <- data.frame(
        ID = sid,
        Treatment = treatment,
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
      doubling_time = dt_val,
      growth_rate = growth_rate,
      r_squared = r_sq,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, result_list)
}

#' Repeated-Measures ANOVA for Tumor Growth
#'
#' Performs a repeated-measures ANOVA on longitudinal tumor volume data,
#' testing the effect of treatment, time, and their interaction.
#'
#' @param df Data frame with tumor measurements.
#' @param time_column Name of the time column. Default "Day".
#' @param volume_column Name of the volume column. Default "Volume".
#' @param treatment_column Name of the treatment column. Default "Treatment".
#' @param id_column Name of the subject ID column. Default "ID".
#' @param transform Transformation to apply: "log" (default), "sqrt", or "none".
#'
#' @return An \code{me_result} object containing the ANOVA table and model.
#'
#' @details
#' Uses \code{lme4::lmer()} to fit a linear mixed model with treatment, time,
#' and their interaction as fixed effects and subject as a random intercept.
#' Then applies \code{stats::anova()} (via lmerTest's Satterthwaite method)
#' to produce the ANOVA table.
#'
#' @examples
#' \dontrun{
#' result <- repeated_measures_anova(my_data)
#' print(result)
#' }
#'
#' @export
repeated_measures_anova <- function(df,
                                    time_column = "Day",
                                    volume_column = "Volume",
                                    treatment_column = "Treatment",
                                    id_column = "ID",
                                    transform = c("log", "sqrt", "none")) {
  transform <- match.arg(transform)
  
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop("Package 'lmerTest' is required for Satterthwaite ANOVA.", call. = FALSE)
  }
  
  # Validate columns
  required <- c(time_column, volume_column, treatment_column, id_column)
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  
  # Transform
  df$y <- switch(transform,
    "log"  = log(df[[volume_column]] + 1),
    "sqrt" = sqrt(df[[volume_column]]),
    "none" = df[[volume_column]]
  )
  
  # Ensure factors
  df[[treatment_column]] <- factor(df[[treatment_column]])
  df[[id_column]] <- factor(df[[id_column]])
  
  # Fit model using lmerTest for Satterthwaite ddf
  formula_str <- paste0("y ~ ", treatment_column, " * ", time_column, " + (1|", id_column, ")")
  model <- lmerTest::lmer(stats::as.formula(formula_str), data = df)
  
  anova_table <- stats::anova(model, type = "III")
  
  new_me_result(
    analysis_type = "repeated_measures_anova",
    data = df,
    results = list(
      model = model,
      anova_table = as.data.frame(anova_table),
      transform = transform
    ),
    summary_df = as.data.frame(anova_table),
    call = match.call()
  )
}
