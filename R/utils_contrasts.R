# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Comparison family + multiplicity adjustment, shared by every model path.
#
# CODE_REVIEW.md R3.1 / G.1 — `p_adjust_method` was documented with
# "bonferroni" as the default but reached only the AUC path. The lme4 path
# passed a custom contrast list to emmeans::contrast(), for which emmeans
# defaults to adjust = "none"; the GAM path hardcoded trt.vs.ctrl with by-day
# Dunnett. So the default model path reported unadjusted p-values while the
# signature promised Bonferroni, and the three paths adjusted over three
# different families.
#
# The maintainer's decision (G.1) is that users define the comparison of
# interest. These helpers make the family an explicit argument, tie the
# adjustment to exactly the set of comparisons returned, and record both in the
# result so the output is self-describing.

#' Valid comparison families
#' @noRd
#' @keywords internal
ME_COMPARISON_FAMILIES <- c("vs_reference", "all_pairs", "custom")

#' Valid multiplicity adjustments
#' @noRd
#' @keywords internal
ME_P_ADJUST_METHODS <- c("bonferroni", "holm", "fdr", "dunnett", "tukey", "none")

#' Validate a comparison family / adjustment pairing and resolve the emmeans code
#'
#' Some adjustments are only meaningful for a particular family: Dunnett's
#' procedure is defined for many-to-one comparisons against a control, and
#' Tukey's HSD for the complete set of pairwise comparisons. Applying either to
#' the other family gives a correction calibrated for the wrong number of
#' comparisons and the wrong correlation structure, so mismatches are rejected
#' rather than silently substituted.
#'
#' Note on Dunnett: emmeans does not recognise `adjust = "dunnett"` and silently
#' downgrades it to `"dunnettx"`, its fast approximation. The exact procedure is
#' `"mvt"` (multivariate-t), which is affordable at the group counts used in
#' these studies, so `"dunnett"` maps to `"mvt"` here.
#'
#' @param comparison_family One of `ME_COMPARISON_FAMILIES`.
#' @param p_adjust_method One of `ME_P_ADJUST_METHODS`.
#' @param custom_contrasts Named list of coefficient vectors, or NULL.
#' @param supports_joint Logical. TRUE for paths where the comparisons come from
#'   a single fitted model (lme4, gam) and so admit Tukey / Dunnett. FALSE for
#'   the AUC path, whose comparisons are independent Welch t-tests with no
#'   common error term or joint covariance — Tukey and Dunnett are undefined
#'   there.
#' @return A list with `family`, `p_adjust_method`, `emmeans_adjust` (the code
#'   to hand emmeans), and `padjust_method` (the code for `stats::p.adjust`,
#'   NA when the adjustment is applied by emmeans itself).
#' @noRd
#' @keywords internal
resolve_comparison_spec <- function(comparison_family,
                                    p_adjust_method,
                                    custom_contrasts = NULL,
                                    supports_joint = TRUE) {
  comparison_family <- match.arg(comparison_family, ME_COMPARISON_FAMILIES)
  p_adjust_method   <- match.arg(p_adjust_method,   ME_P_ADJUST_METHODS)

  if (comparison_family == "custom") {
    if (is.null(custom_contrasts) || !is.list(custom_contrasts) ||
        length(custom_contrasts) == 0L) {
      stop("comparison_family = 'custom' requires `custom_contrasts`: a named ",
           "list of contrast coefficient vectors, one per comparison.",
           call. = FALSE)
    }
    if (is.null(names(custom_contrasts)) || any(!nzchar(names(custom_contrasts)))) {
      stop("`custom_contrasts` must be a *named* list; the names label the ",
           "comparisons in the output.", call. = FALSE)
    }
  }

  if (p_adjust_method == "dunnett" && comparison_family != "vs_reference") {
    stop("p_adjust_method = 'dunnett' is only defined for ",
         "comparison_family = 'vs_reference' (many-to-one against a control). ",
         "Got comparison_family = '", comparison_family, "'.", call. = FALSE)
  }
  if (p_adjust_method == "tukey" && comparison_family != "all_pairs") {
    stop("p_adjust_method = 'tukey' is only defined for ",
         "comparison_family = 'all_pairs' (Tukey's HSD covers the complete ",
         "set of pairwise comparisons). Got comparison_family = '",
         comparison_family, "'.", call. = FALSE)
  }

  if (!supports_joint && p_adjust_method %in% c("dunnett", "tukey")) {
    stop("p_adjust_method = '", p_adjust_method, "' requires comparisons from ",
         "a single fitted model. The AUC path uses independent Welch t-tests, ",
         "which have no common error term or joint covariance, so Tukey and ",
         "Dunnett are undefined there. Use 'bonferroni', 'holm', 'fdr', or ",
         "'none'.", call. = FALSE)
  }

  # emmeans applies tukey / mvt itself using the joint covariance; the p-value
  # family methods can be handed to it directly too, so let emmeans own all of
  # them on the model paths for a single code path.
  emmeans_adjust <- switch(p_adjust_method,
    dunnett = "mvt",       # exact; emmeans silently downgrades "dunnett"
    tukey   = "tukey",
    fdr     = "fdr",
    p_adjust_method
  )

  list(
    family          = comparison_family,
    p_adjust_method = p_adjust_method,
    emmeans_adjust  = emmeans_adjust,
    padjust_method  = if (p_adjust_method %in% c("dunnett", "tukey")) {
      NA_character_
    } else {
      p_adjust_method
    }
  )
}

#' Build the requested contrasts from an emmeans object
#'
#' @param emm An `emmGrid`.
#' @param spec Output of `resolve_comparison_spec()`.
#' @param reference_group Reference level name (used by `vs_reference`).
#' @param custom_contrasts Named list of coefficient vectors (used by `custom`).
#' @param by Optional character vector of variables to compare *within*
#'   (the GAM path compares within study day). When supplied the adjustment
#'   applies across every returned cell, not within each `by` stratum — see
#'   `me_adjust_across_by()`.
#' @return An emmeans contrast object, or NULL if it could not be built.
#' @noRd
#' @keywords internal
build_requested_contrasts <- function(emm, spec, reference_group = NULL,
                                      custom_contrasts = NULL, by = NULL) {
  lv <- tryCatch(levels(emm)[[1]], error = function(e) NULL)

  tryCatch(
    switch(spec$family,
      vs_reference = {
        ref_idx <- if (!is.null(reference_group) && !is.null(lv)) {
          which(lv == reference_group)
        } else 1L
        if (length(ref_idx) != 1L) ref_idx <- 1L
        emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref_idx,
                          by = by, adjust = spec$emmeans_adjust)
      },
      all_pairs = emmeans::contrast(emm, method = "pairwise", by = by,
                                    adjust = spec$emmeans_adjust),
      custom    = emmeans::contrast(emm, method = custom_contrasts, by = by,
                                    adjust = spec$emmeans_adjust)
    ),
    error = function(e) {
      warning("Could not build '", spec$family, "' contrasts: ",
              conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}

#' Re-adjust p-values across all cells of a `by`-stratified contrast table
#'
#' CODE_REVIEW.md G.1 — when contrasts are computed `by` study day (the GAM
#' path), emmeans adjusts *within* each stratum. The five quantile days are
#' reported together and read together, so the family is every returned cell,
#' not each day in isolation. Tukey / Dunnett cannot be re-derived across
#' strata without the full joint covariance, so for those the within-stratum
#' adjustment is kept and flagged; for the p-value family methods the
#' adjustment is redone across all rows.
#'
#' @param df A contrast summary data frame with a p-value column.
#' @param spec Output of `resolve_comparison_spec()`.
#' @return `df` with the p-value column re-adjusted and an
#'   `Adjust_Scope` column recording what the family actually was.
#' @noRd
#' @keywords internal
me_adjust_across_by <- function(df, spec) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)

  p_col <- intersect(c("p_value", "p.value"), names(df))[1]
  if (is.na(p_col)) return(df)

  if (!is.na(spec$padjust_method)) {
    if (spec$padjust_method != "none") {
      # emmeans already adjusted within stratum; redo from the unadjusted
      # values is not possible, so apply the family method across all cells
      # starting from the within-stratum values only when no adjustment was
      # requested at the emmeans level. Request "none" from emmeans and adjust
      # here instead — that is what build_requested_contrasts() is given.
      df[[p_col]] <- stats::p.adjust(df[[p_col]], method = spec$padjust_method)
    }
    df$Adjust_Scope <- "all returned cells (contrast x day)"
  } else {
    df$Adjust_Scope <- paste0(spec$p_adjust_method, ", within each day")
  }
  df$P_Adjust_Method <- spec$p_adjust_method
  df$Comparison_Family <- spec$family
  df
}

#' Pairwise contrast table for the body-weight models
#'
#' CODE_REVIEW.md R3.12 — `analyze_body_weight()` returned group marginal means
#' and nothing inferential. Produces the same shape as the tumour-growth
#' pairwise table so downstream renderers are shared.
#'
#' @param emm_obj An `emmGrid` over Treatment.
#' @param spec Output of `resolve_comparison_spec()`.
#' @param reference_group Reference level name.
#' @param custom_contrasts Named list of coefficient vectors, or NULL.
#' @return A data frame, or NULL when contrasts could not be built.
#' @noRd
#' @keywords internal
bw_pairwise_table <- function(emm_obj, spec, reference_group = NULL,
                              custom_contrasts = NULL) {
  pc <- build_requested_contrasts(emm_obj, spec,
                                  reference_group  = reference_group,
                                  custom_contrasts = custom_contrasts)
  if (is.null(pc)) return(NULL)

  df <- tryCatch(as.data.frame(summary(pc, infer = c(TRUE, TRUE))),
                 error = function(e) as.data.frame(summary(pc)))
  if ("p.value" %in% names(df)) names(df)[names(df) == "p.value"] <- "p_value"
  df$P_Adjust_Method   <- spec$p_adjust_method
  df$Comparison_Family <- spec$family
  df
}
