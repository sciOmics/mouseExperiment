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
#' @param cage_handling Output of \code{resolve_cage_handling()}: a list with
#'   \code{fixed} / \code{random} logicals and a \code{reason} string.
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
                                cage_handling,
                                verbose,
                                include_necrotic_covariate = FALSE) {


  bt <- function(x) paste0("`", x, "`")

  # Choose basis dimension from the number of unique days, clamped to [3, 10].
  # A study with only 4 timepoints can't support k = 10 (mgcv will error).
  n_days <- length(unique(analysis_df[[time_column]]))
  k_val  <- max(3L, min(10L, n_days - 1L))

  # CODE_REVIEW.md R3.17 — cage placement comes from the resolved design
  # structure, not a chi-square p-value.
  include_cage_fixed  <- isTRUE(cage_handling$fixed)
  include_cage_random <- isTRUE(cage_handling$random)

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

  # gamm4 requires the `by` variable in s(x, by = ...) to be a factor.
  # If treatment / cage are character vectors at this point, gamm4 fails
  # with the cryptic message "Can't find by variable". Coerce here so the
  # GAMM4 fit doesn't depend on the caller having done it. Same for
  # cage_column when it's used as a fixed-effect by-variable.
  if (!is.factor(analysis_df[[treatment_column]])) {
    analysis_df[[treatment_column]] <- factor(analysis_df[[treatment_column]])
  }
  if ((include_cage_fixed || include_cage_random) &&
      !is.null(cage_column) && cage_column %in% colnames(analysis_df) &&
      !is.factor(analysis_df[[cage_column]])) {
    analysis_df[[cage_column]] <- factor(analysis_df[[cage_column]])
  }

  # Fit. gamm4 returns list(mer = lmer model, gam = mgcv gam-like object).
  fit_err <- NULL
  fit <- tryCatch(
    gamm4::gamm4(
      formula  = formula_obj,
      random   = random_part,
      data     = analysis_df,
      REML     = TRUE
    ),
    error = function(e) {
      fit_err <<- conditionMessage(e)
      warning("gamm4() failed: ", fit_err,
              "\nReturning NULL; caller should fall back.")
      NULL
    }
  )

  if (is.null(fit)) {
    # Surface the underlying gamm4 message via an attribute so the caller's
    # stop() can include it. Without this the dashboard's safe_analysis
    # only ever sees the generic "failed to fit" stop message and the
    # actual cause is hidden in server logs.
    out <- NULL
    attr(out, "gamm4_error") <- if (!is.null(fit_err)) fit_err else "unknown gamm4 failure"
    return(out)
  }

  # gamm4's $gam component is a "stub" mgcv::gam object missing two things
  # that downstream emmeans dispatch needs:
  #   * class vector lacks "glm" and "lm" (emmeans dispatches by class)
  #   * `$call` is NULL (emmeans::recover_data.gam reads class(object$call)
  #     and rejects with "Can't handle an object of class 'NULL'")
  # Both downstream `tgs_gam_treatment_effects` and
  # `tgs_gam_treatment_effects_over_time` call emmeans on `fit$gam`; without
  # this patch they silently return NULL and the dashboard's TG GAM result
  # has empty Treatment Effects + Pairwise Comparisons tables. Patch here so
  # both consumers benefit.
  fit <- patch_gamm4_stub(fit)

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
                             day_range, reference_group,
                             comparison_spec = NULL,
                             custom_contrasts = NULL) {
  # CODE_REVIEW.md R3.1 — this used to hardcode trt.vs.ctrl with emmeans'
  # default by-day dunnettx, so `p_adjust_method` never reached it and the
  # family was silently "within one day" while five correlated days were
  # reported together.
  if (is.null(comparison_spec)) {
    comparison_spec <- resolve_comparison_spec("vs_reference", "bonferroni")
  }

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

  # For the p-value family methods, ask emmeans for unadjusted values and
  # adjust across every returned (contrast x day) cell afterwards — the five
  # quantile days are reported and read together, so that is the real family.
  # Tukey / Dunnett cannot be re-derived across strata without the full joint
  # covariance, so those stay within-day and are labelled as such.
  spec_for_emmeans <- comparison_spec
  if (!is.na(comparison_spec$padjust_method)) {
    spec_for_emmeans$emmeans_adjust <- "none"
  }

  pc <- build_requested_contrasts(
    emm_time, spec_for_emmeans,
    reference_group  = reference_group,
    custom_contrasts = custom_contrasts,
    by               = time_column
  )
  if (is.null(pc)) return(NULL)

  df <- as.data.frame(summary(pc))
  # Keep the column name in line with the LMM emmeans output
  if ("p.value" %in% names(df)) names(df)[names(df) == "p.value"] <- "p_value"

  me_adjust_across_by(df, comparison_spec)
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

  # Concurvity check — GAM analog of multicollinearity for smooths.
  # Values close to 1 indicate one smooth can be predicted from the
  # others. Returns one row per smooth × concurvity-component (worst,
  # observed, estimate) when computable.
  concurvity_df <- tryCatch({
    cc <- mgcv::concurvity(gam_obj, full = FALSE)
    # cc is a list of three matrices (worst, observed, estimate). Each
    # column is a term in the model. Reshape to long for display.
    if (!is.list(cc) || length(cc) == 0L) return(NULL)
    nms <- names(cc)
    do.call(rbind, lapply(nms, function(metric) {
      mat <- cc[[metric]]
      data.frame(
        Component = metric,
        Term      = colnames(mat),
        Value     = round(diag(mat), 4L),
        stringsAsFactors = FALSE
      )
    }))
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

  # Standard residual diagnostic ggplots so the dashboard's TG
  # Diagnostics tab renders for GAM the same way as for LME4 / AUC.
  rd <- build_residual_diagnostic_plots(gam_obj,
                                        title_prefix = "Tumour growth GAMM")

  list(
    residuals = list(
      fitted    = fitted,
      residuals = resids,
      qq_plot   = qq
    ),
    random_effects           = random_effects,
    variance_components      = tryCatch(lme4::VarCorr(mer_obj), error = function(e) NULL),
    k_check                  = k_check,
    concurvity               = concurvity_df,
    deviance_explained       = dev_expl,
    diag_qq_plot             = rd$diag_qq_plot,
    diag_resid_fitted_plot   = rd$diag_resid_fitted_plot,
    diag_scale_location_plot = rd$diag_scale_location_plot
  )
}


#' Repair a gamm4 `$gam` stub so emmeans can dispatch on it
#'
#' gamm4's `$gam` element is missing two things emmeans needs: `c("glm", "lm")`
#' in its class vector, and a non-NULL `$call` (emmeans::recover_data.gam reads
#' `class(object$call)` and rejects with "Can't handle an object of class
#' 'NULL'"). Extracted from `tgs_fit_gamm4_model()` in CODE_REVIEW.md R3.4 so
#' `analyze_body_weight()` gets the same treatment — it fits gamm4 inline and
#' silently returned empty marginal-means tables without this.
#'
#' @param fit A `gamm4` result (`list(mer =, gam =)`).
#' @return The same object with `$gam` patched.
#' @noRd
#' @keywords internal
patch_gamm4_stub <- function(fit) {
  if (is.null(fit) || is.null(fit$gam)) return(fit)
  if (!all(c("glm", "lm") %in% class(fit$gam))) {
    class(fit$gam) <- unique(c(class(fit$gam), "glm", "lm"))
  }
  if (is.null(fit$gam$call)) {
    fit$gam$call <- call("gam", formula = fit$gam$formula,
                         data = quote(analysis_df))
  }
  fit
}
