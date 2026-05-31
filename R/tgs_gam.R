# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Fit a generalized additive mixed model via gamm4
#'
#' Internal helper used by \code{tumor_growth_statistics()} and
#' \code{analyze_body_weight()} when \code{model_type = "gam"}. Returns a
#' \code{gamm4} object plus model-selection metadata in the same shape as
#' \code{tgs_fit_lme4_models()} so downstream wrappers can branch minimally.
#'
#' The model has the form
#' \deqn{y ~ Treatment + s(Day, by = Treatment, k = k) + (cage)}
#' with random effects \code{(1 | ID)} or \code{(1 | Cage) + (1 | ID)}.
#' Basis dimension \code{k} is auto-chosen as \code{min(10, n_days - 1)} but
#' at least 3, so the smoother can't ask for more knots than data supports.
#'
#' @param analysis_df Data frame on the modelling scale (post-transform).
#' @param response_column Column name of the response variable on the
#'   modelling scale (e.g. transformed Volume or Net_Weight).
#' @param time_column,treatment_column,id_column,cage_column Column names.
#' @param handle_cage_effects One of \code{"include_if_not_collinear"},
#'   \code{"always_include"}, \code{"never_include"}, \code{"as_random_effect"}.
#' @param cage_collinear Logical from \code{tgs_compute_cage_effects()}.
#' @param verbose Logical.
#' @param include_necrotic_covariate Logical; adds \code{necrotic_cov_flag}.
#' @return List with \code{model} (gamm4 result), \code{model_selection}
#'   (AIC/BIC of the mer component + \code{selected_model = "gam"}), and
#'   \code{best_model = "gam"}. Returns \code{NULL} when gamm4 errors.
#' @keywords internal
#' @noRd
tgs_fit_gamm4_model <- function(analysis_df,
                                response_column,
                                time_column,
                                treatment_column,
                                id_column,
                                cage_column,
                                handle_cage_effects,
                                cage_collinear,
                                verbose,
                                include_necrotic_covariate = FALSE) {

  if (!requireNamespace("gamm4", quietly = TRUE)) {
    stop(
      "Package 'gamm4' is required for model_type = 'gam'.\n",
      "Install it with: install.packages('gamm4')"
    )
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop(
      "Package 'mgcv' is required for model_type = 'gam'.\n",
      "Install it with: install.packages('mgcv')"
    )
  }

  bt <- function(x) paste0("`", x, "`")

  # Choose basis dimension from the number of unique days, clamped to [3, 10].
  # A study with only 4 timepoints can't support k = 10 (mgcv will error).
  n_days <- length(unique(analysis_df[[time_column]]))
  k_val  <- max(3L, min(10L, n_days - 1L))

  # Cage in fixed vs random vs not modelled
  include_cage_fixed <- switch(handle_cage_effects,
    include_if_not_collinear = !cage_collinear,
    always_include           = TRUE,
    never_include            = FALSE,
    as_random_effect         = FALSE
  )
  include_cage_random <- handle_cage_effects == "as_random_effect"

  # Build fixed-part formula: Treatment main effect + factor-by smoother.
  # The Treatment main effect is needed because s(Day, by=Treatment) is
  # centred (sum-to-zero) — without the main effect each group's mean is
  # absorbed into the smoother intercept and emmeans contrasts become wrong.
  fixed_part <- paste0(
    bt(response_column), " ~ ",
    bt(treatment_column),
    " + s(", bt(time_column), ", by = ", bt(treatment_column),
    ", k = ", k_val, ", bs = \"tp\")"
  )
  if (include_cage_fixed) {
    fixed_part <- paste(fixed_part, "+", bt(cage_column))
  }
  if (isTRUE(include_necrotic_covariate)) {
    fixed_part <- paste(fixed_part, "+ necrotic_cov_flag")
  }

  # Random-effects formula for gamm4 (must be a one-sided formula).
  random_part <- if (include_cage_random) {
    stats::as.formula(paste0(
      "~ (1 | ", bt(cage_column), ") + (1 | ", bt(id_column), ")"
    ))
  } else {
    stats::as.formula(paste0("~ (1 | ", bt(id_column), ")"))
  }

  formula_obj <- stats::as.formula(fixed_part)

  if (isTRUE(verbose)) {
    message("GAMM4 fixed: ", deparse1(formula_obj))
    message("GAMM4 random: ", deparse1(random_part))
    message("Smoother basis dimension k = ", k_val,
            " (auto-chosen from ", n_days, " unique time points)")
  }

  # Fit. gamm4 returns list(mer = lmer model, gam = mgcv gam-like object).
  fit <- tryCatch(
    gamm4::gamm4(
      formula  = formula_obj,
      random   = random_part,
      data     = analysis_df,
      REML     = TRUE
    ),
    error = function(e) {
      warning("gamm4() failed: ", conditionMessage(e),
              "\nReturning NULL; caller should fall back.")
      NULL
    }
  )

  if (is.null(fit)) return(NULL)

  model_selection <- list(
    aic            = stats::AIC(fit$mer),
    bic            = stats::BIC(fit$mer),
    selected_model = "gam",
    k_basis        = k_val
  )

  list(
    model           = fit,
    model_selection = model_selection,
    best_model      = "gam",
    formula_string  = paste0(deparse1(formula_obj), "  random: ", deparse1(random_part))
  )
}


#' Treatment effects (EMMs) at the study-mean day, GAM version
#'
#' Returns a data frame with one row per treatment group and columns
#' \code{Group}, \code{Adjusted_Mean}, \code{SE}, \code{DF}, \code{Lower_CL},
#' \code{Upper_CL}, \code{Note} — the same shape as the LMM
#' \code{treatment_effects} so the dashboard's render path is shared.
#'
#' @param gam_obj The \code{$gam} element of a \code{gamm4()} result.
#' @param treatment_column,time_column Column names.
#' @param mean_day Numeric scalar — day at which to marginalise.
#' @param reference_group Character — used for the Reference Note column and
#'   to move that row to the top.
#' @keywords internal
#' @noRd
tgs_gam_treatment_effects <- function(gam_obj, treatment_column, time_column,
                                      mean_day, reference_group) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  emm <- tryCatch(
    emmeans::emmeans(
      gam_obj,
      specs = treatment_column,
      at    = stats::setNames(list(mean_day), time_column)
    ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)

  emm_summary <- as.data.frame(summary(emm))
  te <- data.frame(
    Group         = as.character(emm_summary[[treatment_column]]),
    Adjusted_Mean = round(emm_summary$emmean,   3),
    SE            = round(emm_summary$SE,       3),
    DF            = emm_summary$df,
    Lower_CL      = round(emm_summary$lower.CL, 3),
    Upper_CL      = round(emm_summary$upper.CL, 3),
    Note          = ifelse(
      as.character(emm_summary[[treatment_column]]) == reference_group,
      "Reference group", ""
    ),
    stringsAsFactors = FALSE
  )

  # Move the reference row to the top so summary tables read consistently.
  ref_idx <- which(te$Group == reference_group)
  if (length(ref_idx) > 0 && ref_idx[1] > 1L) {
    te <- rbind(te[ref_idx[1], ], te[-ref_idx[1], ])
    rownames(te) <- NULL
  }
  te
}


#' Treatment EMMs at the 5 study-day quantiles (time-resolved), GAM version
#'
#' One row per (Group, Day) combination at the min / Q1 / median / Q3 / max
#' study days. Mirrors the LMM \code{treatment_effects_over_time} field.
#'
#' @param gam_obj,treatment_column,time_column See above.
#' @param day_range Numeric vector of all observed study days.
#' @keywords internal
#' @noRd
tgs_gam_emm_time <- function(gam_obj, treatment_column, time_column, day_range) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  quant_days <- unique(round(stats::quantile(
    day_range, probs = c(0, 0.25, 0.5, 0.75, 1), type = 1
  )))

  emm_time <- tryCatch(
    emmeans::emmeans(
      gam_obj,
      specs = c(treatment_column, time_column),
      at    = stats::setNames(list(quant_days), time_column)
    ),
    error = function(e) NULL
  )
  if (is.null(emm_time)) return(NULL)

  df <- as.data.frame(summary(emm_time))
  names(df)[names(df) == "emmean"]   <- "Adjusted_Mean"
  names(df)[names(df) == "lower.CL"] <- "Lower_CL"
  names(df)[names(df) == "upper.CL"] <- "Upper_CL"
  df
}


#' Pairwise contrasts vs the reference group, by day, for a GAM
#'
#' Returns one row per (Day, contrast) — the requested "difference of smooths
#' at quantile days" output. Matches the LMM \code{pairwise_comparisons} shape
#' so the dashboard's table renderer is shared.
#'
#' @param gam_obj,treatment_column,time_column See above.
#' @param day_range,reference_group See above.
#' @keywords internal
#' @noRd
tgs_gam_pairwise <- function(gam_obj, treatment_column, time_column,
                             day_range, reference_group) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  quant_days <- unique(round(stats::quantile(
    day_range, probs = c(0, 0.25, 0.5, 0.75, 1), type = 1
  )))

  emm_time <- tryCatch(
    emmeans::emmeans(
      gam_obj,
      specs = c(treatment_column, time_column),
      at    = stats::setNames(list(quant_days), time_column)
    ),
    error = function(e) NULL
  )
  if (is.null(emm_time)) return(NULL)

  # Per-day pairwise comparisons against the reference group.
  pc <- tryCatch(
    emmeans::contrast(
      emm_time,
      method = "trt.vs.ctrl",
      ref    = reference_group,
      by     = time_column
    ),
    error = function(e) NULL
  )
  if (is.null(pc)) return(NULL)

  df <- as.data.frame(summary(pc))
  # Keep the column name in line with the LMM emmeans output
  if ("p.value" %in% names(df)) names(df)[names(df) == "p.value"] <- "p_value"
  df
}


#' Smooth-term significance table (mgcv) for the GAM ANOVA panel
#'
#' Wraps the \code{s.table} from \code{summary(gam_obj)} into a tidy
#' data frame for the dashboard's ANOVA tab.
#'
#' @keywords internal
#' @noRd
tgs_gam_anova_table <- function(gam_obj) {
  smry <- tryCatch(summary(gam_obj), error = function(e) NULL)
  if (is.null(smry) || is.null(smry$s.table)) return(NULL)

  st <- as.data.frame(smry$s.table)
  st$Term <- rownames(st)
  rownames(st) <- NULL
  st <- st[, c("Term", setdiff(names(st), "Term"))]
  # Standardise the p-value column name
  if ("p-value" %in% names(st)) names(st)[names(st) == "p-value"] <- "p_value"
  st
}


#' Diagnostics for a fitted GAMM4 result
#'
#' Returns a list with the same field set as the LMM \code{diagnostics}
#' element so the dashboard tab can render either path. Adds GAM-specific
#' fields: \code{k_check} (basis-dimension adequacy) and \code{deviance_explained}.
#'
#' @keywords internal
#' @noRd
tgs_gam_diagnostics <- function(fit, id_column) {
  gam_obj <- fit$gam
  mer_obj <- fit$mer

  resids <- stats::residuals(gam_obj)
  fitted <- stats::fitted(gam_obj)
  qq     <- stats::qqnorm(resids, plot.it = FALSE)

  # k-index check: low k.index (<<1) or p-value near 0 suggests the
  # basis dimension is too low.
  k_check <- tryCatch({
    utils::capture.output(kc <- mgcv::k.check(gam_obj))
    as.data.frame(kc)
  }, error = function(e) NULL)

  dev_expl <- tryCatch(
    summary(gam_obj)$dev.expl,
    error = function(e) NA_real_
  )

  random_effects <- tryCatch({
    re <- lme4::ranef(mer_obj)
    # mer_obj's random-effect grouping name is set internally by gamm4 and
    # may not match `id_column` verbatim, so we just return the first
    # available random-effect data frame as "intercepts".
    list(
      intercepts = if (length(re) > 0) re[[1]] else NULL,
      slopes     = NULL
    )
  }, error = function(e) list(intercepts = NULL, slopes = NULL))

  list(
    residuals = list(
      fitted    = fitted,
      residuals = resids,
      qq_plot   = qq
    ),
    random_effects      = random_effects,
    variance_components = tryCatch(lme4::VarCorr(mer_obj), error = function(e) NULL),
    k_check             = k_check,
    deviance_explained  = dev_expl
  )
}
