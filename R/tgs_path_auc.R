# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# AUC model path for tumor_growth_statistics()
#
# Extracted from R/tumor_growth_statistics.R as part of CODE_REVIEW.md Round 2
# D.1 (1000+ LOC files). The function body is large because the AUC path
# contains the full one-way-ANOVA + Welch pairwise + treatment effects +
# diagnostics + summary metadata blocks. Splitting it out leaves
# tumor_growth_statistics.R substantially smaller and groups all AUC-path
# logic in one file.
#
# Behaviour is bit-identical to the inline version — this is a pure
# code-organisation refactor with no semantic change.

#' AUC model path
#'
#' Internal helper consumed by \code{\link{tumor_growth_statistics}} when
#' \code{model_type = "auc"}.
#'
#' CODE_REVIEW.md R20.4 / R20-K. The AUC was each animal's trapezoid over its
#' own follow-up, compared by Welch t-tests. Animals removed at the volume
#' limit -- controls first -- integrated over the shortest windows, so in 53.5 %
#' of simulated studies an effective drug looked worse than control (R3.3's
#' fix was recorded but never applied here). The AUC is now the area under
#' each arm's fitted curve (the endpoint model: per-arm spline on the log
#' scale, random slopes) over one window common to every arm, from the first
#' study day to the last day on which every arm is evaluable. Arms are
#' compared by the ratio of their AUCs, with intervals and p-values from draws
#' of the model's fixed effects. The per-animal trapezoids stay in
#' \code{auc_analysis$individual} as a descriptive table only.
#'
#' @param auc_analysis Output of \code{tgs_compute_auc()} (descriptive).
#' @param growth_rates Output of \code{tgs_compute_growth_rates()}.
#' @param auc_df Untransformed working data frame.
#' @param transform,reference_group See \code{\link{tumor_growth_statistics}}.
#' @param comparison_spec Output of \code{resolve_comparison_spec()}.
#' @param return_model,include_diagnostics See
#'   \code{\link{tumor_growth_statistics}}.
#' @param auc_bootstrap_n,auc_permutations Ignored since v0.26.0, with a
#'   warning when positive: they belonged to the per-animal t-tests.
#' @param auc_bootstrap_seed Seed for the model draws; \code{NULL} uses a fixed
#'   seed so the reported interval is reproducible.
#' @param cage_analysis,data_summary,necrosis_summary Summaries assembled in
#'   the main function.
#' @param id_column,treatment_column,cage_column,time_column,volume_column
#'   Column names.
#' @noRd
tgs_path_auc <- function(auc_analysis,
                         growth_rates,
                         auc_df,
                         transform,
                         comparison_spec,
                         reference_group,
                         return_model,
                         include_diagnostics,
                         auc_bootstrap_n,
                         auc_bootstrap_seed,
                         auc_permutations,
                         cage_analysis,
                         data_summary,
                         necrosis_summary,
                         id_column,
                         treatment_column,
                         cage_column,
                         time_column,
                         volume_column) {

  if (isTRUE(auc_bootstrap_n > 0) || isTRUE(auc_permutations > 0)) {
    warning("auc_bootstrap_n and auc_permutations are ignored: since v0.26.0 ",
            "the AUC is model-based, with intervals and p-values from draws of ",
            "the model (CODE_REVIEW.md R20.4).", call. = FALSE)
  }

  d <- data.frame(
    MouseKey  = make_mouse_key(as.character(auc_df[[treatment_column]]),
                               as.character(auc_df[[id_column]]),
                               as.character(auc_df[[cage_column]])),
    Treatment = as.character(auc_df[[treatment_column]]),
    Day       = as.numeric(auc_df[[time_column]]),
    Volume    = as.numeric(auc_df[[volume_column]]),
    stringsAsFactors = FALSE)
  d <- d[is.finite(d$Day) & is.finite(d$Volume), , drop = FALSE]

  # The common window: first study day to the last day every arm is evaluable.
  evaluability <- me_evaluability(d)
  t_end   <- me_resolve_eval_day(evaluability, NULL)
  t_start <- min(d$Day)
  if (!(t_end > t_start)) {
    stop("The AUC window is empty: the last day on which every arm is ",
         "evaluable (", t_end, ") is the first study day. ", evaluability$rule,
         call. = FALSE)
  }
  em <- me_endpoint_model(d)
  if (is.null(em)) {
    stop("The AUC model could not be fitted to these data.", call. = FALSE)
  }
  arms <- evaluability$arms
  if (!is.null(reference_group) && !reference_group %in% arms) {
    stop("Reference group '", reference_group, "' is not in the data.", call. = FALSE)
  }

  grid <- seq(t_start, t_end, length.out = 201L)
  h <- diff(grid)
  trapz_rows <- function(M) as.numeric(M[, -1L, drop = FALSE] %*% (h / 2) +
                                       M[, -ncol(M), drop = FALSE] %*% (h / 2))
  seed  <- if (is.null(auc_bootstrap_seed)) 20260930L else auc_bootstrap_seed
  draws <- me_beta_draws(em, 4000L, seed)
  X_by  <- lapply(stats::setNames(arms, arms), function(a) me_endpoint_X(em, a, grid))
  auc_hat <- vapply(X_by, function(X) {
    trapz_rows(matrix(exp(as.numeric(X %*% em$beta)), nrow = 1L))
  }, numeric(1L))
  auc_draw <- vapply(X_by, function(X) trapz_rows(exp(draws %*% t(X))),
                     numeric(nrow(draws)))

  n_by <- table(factor(unique(d[, c("MouseKey", "Treatment")])$Treatment,
                       levels = arms))
  treatment_effects <- data.frame(
    Treatment = arms,
    AUC       = unname(auc_hat),
    Lower_CL  = apply(auc_draw, 2L, stats::quantile, 0.025, names = FALSE),
    Upper_CL  = apply(auc_draw, 2L, stats::quantile, 0.975, names = FALSE),
    N         = as.integer(n_by),
    Reference = if (is.null(reference_group)) FALSE else arms == reference_group,
    stringsAsFactors = FALSE)
  rownames(treatment_effects) <- NULL

  # Comparisons: the ratio of AUCs (log scale for the test), and the difference.
  pairs <- switch(comparison_spec$family,
    all_pairs = utils::combn(arms, 2, simplify = FALSE),
    vs_reference = {
      if (is.null(reference_group)) {
        warning("comparison_family = 'vs_reference' without a reference group; ",
                "using all pairs.", call. = FALSE)
        utils::combn(arms, 2, simplify = FALSE)
      } else {
        lapply(setdiff(arms, reference_group), function(g) c(g, reference_group))
      }
    },
    custom = stop("comparison_family = 'custom' is not supported for ",
                  "model_type = 'auc'. Use 'vs_reference' or 'all_pairs'.",
                  call. = FALSE)
  )
  pairwise_df <- do.call(rbind, lapply(pairs, function(pr) {
    lr_draw <- log(auc_draw[, pr[1]]) - log(auc_draw[, pr[2]])
    lr_hat  <- log(auc_hat[[pr[1]]]) - log(auc_hat[[pr[2]]])
    se      <- stats::sd(lr_draw)
    z       <- lr_hat / se
    diff_draw <- auc_draw[, pr[1]] - auc_draw[, pr[2]]
    data.frame(
      comparison  = paste(pr[1], "-", pr[2]),
      estimate    = auc_hat[[pr[1]]] - auc_hat[[pr[2]]],
      ci_lower    = stats::quantile(diff_draw, 0.025, names = FALSE),
      ci_upper    = stats::quantile(diff_draw, 0.975, names = FALSE),
      ratio       = exp(lr_hat),
      ratio_lower = exp(stats::quantile(lr_draw, 0.025, names = FALSE)),
      ratio_upper = exp(stats::quantile(lr_draw, 0.975, names = FALSE)),
      z_value     = z,
      p_value     = 2 * stats::pnorm(-abs(z)),
      stringsAsFactors = FALSE)
  }))
  pairwise_df$p_adjusted <- stats::p.adjust(pairwise_df$p_value,
                                            method = comparison_spec$padjust_method)
  pairwise_df$p_adjust_method   <- comparison_spec$p_adjust_method
  pairwise_df$Comparison_Family <- comparison_spec$family

  # Omnibus: Wald test that every arm has the same log AUC.
  ref_arm <- if (!is.null(reference_group)) reference_group else arms[1L]
  others  <- setdiff(arms, ref_arm)
  L_hat   <- log(auc_hat[others]) - log(auc_hat[[ref_arm]])
  L_draw  <- log(auc_draw[, others, drop = FALSE]) - log(auc_draw[, ref_arm])
  S       <- stats::cov(L_draw)
  wald    <- tryCatch(as.numeric(t(L_hat) %*% solve(S, L_hat)), error = function(e) NA_real_)
  anova_table <- data.frame(
    Term    = "Treatment",
    Chisq   = wald,
    Df      = length(others),
    p_value = stats::pchisq(wald, df = length(others), lower.tail = FALSE),
    Method  = "Wald test that every arm has the same model-based AUC (log scale)",
    stringsAsFactors = FALSE)

  window_note <- sprintf(paste0(
    "AUC is the area under each arm's fitted curve from day %s to day %s, ",
    "the last day on which every arm is evaluable. %s"),
    format(t_start), format(t_end), evaluability$rule)

  analysis_summary <- list(
    analysis_type = "Area Under the Curve (AUC) Analysis",
    data_description = list(
      subjects = length(unique(d$MouseKey)),
      treatment_groups = length(arms),
      time_points      = length(unique(d$Day)),
      reference_group  = reference_group
    ),
    methods = list(
      volume_transformation  = "log (the curve is fitted on the log scale and integrated on the volume scale)",
      transform_requested    = transform,
      auc_calculation_method = "area under each arm's fitted geometric-mean curve over a common window",
      auc_window             = c(start = t_start, end = t_end),
      model                  = me_endpoint_model_info(em)$formula,
      statistical_test       = anova_table$Method,
      posthoc_method         = paste0(
        "ratio of AUCs (", comparison_spec$family, "), Wald tests from draws of ",
        "the model, with ", comparison_spec$p_adjust_method, " adjustment"),
      individual_calculation = paste(
        "Per-animal trapezoidal AUCs are reported for description only: each",
        "covers the animal's own follow-up, so they are not comparable across",
        "arms with different attrition."),
      growth_rate_calculation = paste0(
        "Growth rates are calculated by fitting a linear regression model to log-transformed volume data over time for each subject ",
        "(non-positive volumes are replaced with half the smallest positive value before the log). ",
        "The slope coefficient from this model represents the exponential growth rate. ",
        "A value of 0.1 indicates approximately 10% tumor volume increase per day. ",
        "Only subjects with 3 or more time points are included in growth rate calculations."
      )
    ),
    notes = c(window_note,
              "Composite IDs combine subject ID, treatment group and cage.")
  )

  posthoc <- list(
    method   = analysis_summary$methods$posthoc_method,
    pairwise = pairwise_df
  )

  rd <- if (include_diagnostics)
    build_residual_diagnostic_plots(em$fit, title_prefix = "Tumour growth AUC model")
  else
    list(diag_qq_plot = NULL, diag_resid_fitted_plot = NULL,
         diag_scale_location_plot = NULL)

  auc_analysis$window <- c(start = t_start, end = t_end)
  auc_analysis$model_based <- treatment_effects

  list(
    model                = if (return_model) em$fit else NULL,
    model_type_used      = "auc",
    meta = me_result_meta(
      analysis_type     = "Area under each arm's fitted curve over the common evaluable window",
      model_type_used   = "auc",
      inference         = "frequentist",
      interval_type     = "confidence",
      # The curve is fitted on the log scale, but every reported number is on
      # the volume scale (volume x day); a consumer must not back-transform.
      transform_used    = "none",
      estimate_scale    = "AUC (volume x day, raw scale)",
      comparison_family = comparison_spec$family,
      p_adjust_method   = comparison_spec$p_adjust_method
    ),
    comparison_family    = comparison_spec$family,
    p_adjust_method_used = comparison_spec$p_adjust_method,
    transform_used       = "none",
    transform_requested  = transform,
    anova                = anova_table,
    summary              = analysis_summary,
    posthoc              = posthoc,
    pairwise_comparisons = me_pairwise_frame(
      posthoc$pairwise, comparison_spec,
      adjusted_col = "p_adjusted", raw_col = "p_value",
      adjust_scope = "all requested AUC comparisons"),
    treatment_effects    = treatment_effects,
    auc_window           = c(start = t_start, end = t_end),
    evaluability         = evaluability,
    endpoint_model       = me_endpoint_model_info(em),
    growth_rates         = growth_rates,
    cage_analysis        = cage_analysis,
    auc_analysis         = auc_analysis,
    data_summary         = data_summary,
    diagnostics          = NULL,
    variance_test        = NULL,
    diag_qq_plot             = rd$diag_qq_plot,
    diag_resid_fitted_plot   = rd$diag_resid_fitted_plot,
    diag_scale_location_plot = rd$diag_scale_location_plot,
    diag_re_qq_plot          = if (include_diagnostics)
      build_random_effects_qq_plot(em$fit, "MouseKey",
                                   title_prefix = "Tumour growth AUC model") else NULL,
    necrosis_summary     = necrosis_summary
  )
}
