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
#' \code{model_type = "auc"}. Returns the same result-list shape as the
#' inline code it replaces.
#'
#' @param auc_analysis Output of \code{tgs_compute_auc()}.
#' @param growth_rates Output of \code{tgs_compute_growth_rates()}.
#' @param auc_df Working data frame with composite IDs.
#' @param transform,reference_group See
#'   \code{\link{tumor_growth_statistics}}.
#' @param comparison_spec Output of \code{resolve_comparison_spec()} — the
#'   validated comparison family and multiplicity adjustment.
#' @param return_model,include_diagnostics See
#'   \code{\link{tumor_growth_statistics}}.
#' @param auc_permutations Number of label permutations for the permutation
#'   p-value (CODE_REVIEW.md H.2). 0 skips.
#' @param auc_bootstrap_n,auc_bootstrap_seed See
#'   \code{\link{tumor_growth_statistics}}.
#' @param cage_analysis,data_summary,necrosis_summary Already-computed
#'   summaries assembled in the main function.
#' @param id_column,treatment_column,cage_column,time_column Column names
#'   from the main function call.
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
                         time_column) {

  # CODE_REVIEW.md R3.22 — the omnibus test must make the same variance
  # assumption as the pairwise tests. `aov()` assumes homoscedasticity while the
  # pairwise tests deliberately use var.equal = FALSE because, per this
  # function's own note, "variances between treatment groups may differ". Use
  # Welch's ANOVA for the omnibus so the two agree, and report a
  # Brown-Forsythe variance-homogeneity check so the choice is evidence-based
  # rather than assumed. The `aov` fit is retained for residual diagnostics
  # only (the QQ / residual plots below need a model object).
  auc_model   <- stats::aov(AUC ~ Treatment, data = auc_analysis$individual)
  welch_anova <- tryCatch(
    stats::oneway.test(AUC ~ Treatment, data = auc_analysis$individual,
                       var.equal = FALSE),
    error = function(e) NULL
  )
  anova_table <- if (!is.null(welch_anova)) {
    data.frame(
      Term    = "Treatment",
      F_value = unname(welch_anova$statistic),
      df_num  = unname(welch_anova$parameter[1]),
      df_den  = unname(welch_anova$parameter[2]),
      p_value = welch_anova$p.value,
      Method  = "Welch's one-way ANOVA (unequal variances)",
      stringsAsFactors = FALSE
    )
  } else {
    stats::anova(auc_model)
  }

  # Brown-Forsythe test of variance homogeneity across groups (median-centred
  # Levene, robust to non-normality). Reported, not acted on: the analysis
  # already does not assume equal variances.
  variance_test <- tryCatch({
    lt <- car::leveneTest(AUC ~ factor(Treatment),
                          data = auc_analysis$individual, center = stats::median)
    data.frame(
      Test    = "Brown-Forsythe (median-centred Levene)",
      F_value = lt[1, "F value"],
      df_num  = lt[1, "Df"],
      df_den  = lt[2, "Df"],
      p_value = lt[1, "Pr(>F)"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)

  # Create pairwise comparisons for AUC using Welch's t-tests
  # This is more appropriate for AUC analysis where variances between groups may differ
  treatments       <- unique(auc_analysis$individual$Treatment)
  pairwise_results <- list()
  pairwise_data    <- list()

  # CODE_REVIEW.md R3.21 — the comparison family is now explicit. This path
  # previously always built all k(k-1)/2 pairs and adjusted over that family
  # even when a reference group was supplied and the table was then reordered
  # to put reference comparisons first, so the same `p_adjust_method` meant
  # materially different stringency here than on the lme4 path (10 tests vs 4
  # at k = 5).
  pairs <- switch(comparison_spec$family,
    all_pairs = utils::combn(treatments, 2, simplify = FALSE),
    vs_reference = {
      if (is.null(reference_group) || !reference_group %in% treatments) {
        warning("comparison_family = 'vs_reference' but the reference group is ",
                "not present in the AUC data; falling back to all pairs.",
                call. = FALSE)
        utils::combn(treatments, 2, simplify = FALSE)
      } else {
        lapply(setdiff(treatments, reference_group),
               function(g) c(g, reference_group))
      }
    },
    custom = stop("comparison_family = 'custom' is not supported for ",
                  "model_type = 'auc': custom contrasts are coefficient ",
                  "vectors over model terms, and the AUC path runs independent ",
                  "Welch t-tests rather than fitting a joint model. Use ",
                  "'vs_reference' or 'all_pairs'.", call. = FALSE)
  )

  for (pair in pairs) {
    group1_data <- auc_analysis$individual$AUC[auc_analysis$individual$Treatment == pair[1]]
    group2_data <- auc_analysis$individual$AUC[auc_analysis$individual$Treatment == pair[2]]

    if (length(group1_data) < 2 || length(group2_data) < 2) {
      # Not enough data, create a placeholder result
      pairwise_results[[paste(pair[1], "-", pair[2])]] <- list(
        comparison       = paste(pair[1], "-", pair[2]),
        mean_diff        = ifelse(length(group1_data) > 0 && length(group2_data) > 0,
                                  mean(group1_data) - mean(group2_data), NA),
        t_value          = NA,
        df               = NA,
        p_value          = NA,
        ci_lower         = NA,
        ci_upper         = NA,
        boot_ci_lower    = NA_real_,
        boot_ci_upper    = NA_real_,
        perm_p_value     = NA_real_
      )
      pairwise_data[[paste(pair[1], "-", pair[2])]] <- list(
        group1 = pair[1],
        group2 = pair[2],
        data1  = group1_data,
        data2  = group2_data,
        result = list(
          statistic = NA,  parameter = NA,  p.value = NA,
          conf.int  = c(NA, NA), estimate = NA
        )
      )
    } else {
      t_test_result <- stats::t.test(group1_data, group2_data, var.equal = FALSE)

      # Optional non-parametric bootstrap CI for the mean difference.
      # Welch's CI is approximate at small N when per-mouse AUC is skewed
      # (especially with extrapolation); the bootstrap percentile CI is
      # the recommended honest alternative.
      boot_ci <- if (auc_bootstrap_n > 0L) {
        tgs_boot_diff_ci(group1_data, group2_data,
                         n_boot = auc_bootstrap_n,
                         seed   = auc_bootstrap_seed)
      } else c(lower = NA_real_, upper = NA_real_)

      # CODE_REVIEW.md H.2 — permutation p-value alongside Welch's. Interval
      # from the bootstrap, p-value from the permutation.
      perm_p <- if (auc_permutations > 0L) {
        # Offset the seed per comparison so the pairs do not all reuse the
        # identical permutation draws (reproducible, but not degenerate).
        pair_seed <- if (is.null(auc_bootstrap_seed)) NULL else {
          as.integer(auc_bootstrap_seed) + length(pairwise_results)
        }
        tgs_perm_diff_p(group1_data, group2_data,
                        n_perm = auc_permutations,
                        seed   = pair_seed)
      } else NA_real_

      pairwise_results[[paste(pair[1], "-", pair[2])]] <- list(
        comparison       = paste(pair[1], "-", pair[2]),
        mean_diff        = mean(group1_data) - mean(group2_data),
        t_value          = t_test_result$statistic,
        df               = t_test_result$parameter,
        p_value          = t_test_result$p.value,
        ci_lower         = t_test_result$conf.int[1],
        ci_upper         = t_test_result$conf.int[2],
        boot_ci_lower    = unname(boot_ci["lower"]),
        boot_ci_upper    = unname(boot_ci["upper"]),
        perm_p_value     = perm_p
      )
      pairwise_data[[paste(pair[1], "-", pair[2])]] <- list(
        group1 = pair[1],
        group2 = pair[2],
        data1  = group1_data,
        data2  = group2_data,
        result = t_test_result
      )
    }
  }

  # Create data frame from pairwise results
  pairwise_df <- do.call(rbind, lapply(names(pairwise_results), function(comp) {
    res <- pairwise_results[[comp]]
    data.frame(
      comparison    = res$comparison,
      estimate      = res$mean_diff,
      t_value       = res$t_value,
      df            = res$df,
      p_value       = res$p_value,
      ci_lower      = res$ci_lower,
      ci_upper      = res$ci_upper,
      boot_ci_lower = res$boot_ci_lower,
      boot_ci_upper = res$boot_ci_upper,
      perm_p_value  = res$perm_p_value,
      stringsAsFactors = FALSE
    )
  }))

  # Apply the multiple-comparison correction over exactly the family built
  # above, and record what was applied so the result is self-describing.
  pairwise_df$p_adjusted        <- stats::p.adjust(
    pairwise_df$p_value, method = comparison_spec$padjust_method)
  if (any(is.finite(pairwise_df$perm_p_value))) {
    pairwise_df$perm_p_adjusted <- stats::p.adjust(
      pairwise_df$perm_p_value, method = comparison_spec$padjust_method)
  }
  pairwise_df$p_adjust_method   <- comparison_spec$p_adjust_method
  pairwise_df$Comparison_Family <- comparison_spec$family

  # Handle reference group - ensure it exists in treatments
  if (!is.null(reference_group) && reference_group %in% treatments) {
    ref_comparisons <- grep(
      paste0("^", reference_group, " -|^[^-]+ - ", reference_group, "$"),
      pairwise_df$comparison
    )
    if (length(ref_comparisons) > 0) {
      pairwise_df <- rbind(
        pairwise_df[ref_comparisons,  , drop = FALSE],
        pairwise_df[-ref_comparisons, , drop = FALSE]
      )
    }
  }

  # Extract treatment effects from AUC
  treatment_effects <- data.frame(
    Treatment = treatments,
    Mean_AUC = sapply(treatments, function(t) {
      mean(auc_analysis$individual$AUC[auc_analysis$individual$Treatment == t])
    }),
    SD = sapply(treatments, function(t) {
      stats::sd(auc_analysis$individual$AUC[auc_analysis$individual$Treatment == t])
    }),
    N = sapply(treatments, function(t) {
      sum(auc_analysis$individual$Treatment == t)
    }),
    stringsAsFactors = FALSE
  )

  # Add reference indicator
  treatment_effects$Reference <- rep(FALSE, nrow(treatment_effects))
  if (!is.null(reference_group) && reference_group %in% treatments) {
    treatment_effects$Reference[treatment_effects$Treatment == reference_group] <- TRUE
  }

  # Create diagnostic plots
  if (include_diagnostics) {
    diagnostics <- list(
      residuals = list(
        fitted    = stats::fitted(auc_model),
        residuals = stats::residuals(auc_model),
        qq_plot   = stats::qqnorm(stats::residuals(auc_model), plot.it = FALSE)
      )
    )
  } else {
    diagnostics <- NULL
  }

  # Create a descriptive summary
  analysis_summary <- list(
    analysis_type = "Area Under the Curve (AUC) Analysis",
    data_description = list(
      subjects = length(unique(make_mouse_key(
        auc_df[[id_column]], auc_df[[treatment_column]], auc_df[[cage_column]]
      ))),
      treatment_groups = length(unique(auc_df[[treatment_column]])),
      time_points      = length(unique(auc_df[[time_column]])),
      reference_group  = reference_group
    ),
    methods = list(
      # CODE_REVIEW.md R3.16 — the AUC path takes its working copy of the data
      # *before* the transform is applied, so AUC is always computed on the raw
      # volume scale regardless of `transform`. Reporting the requested
      # transform here put a false statement ("Volume data was log transformed
      # prior to analysis") into the dashboard methods panel and the HTML
      # report export. Report what actually happened.
      volume_transformation  = "none (AUC is computed on the raw volume scale)",
      transform_requested    = transform,
      auc_calculation_method = "trapezoidal",
      statistical_test       = "Welch's one-way ANOVA on per-animal AUC values (unequal variances)",
      posthoc_method         = paste0(
        "Welch's t-tests (", comparison_spec$family, ") with ",
        comparison_spec$p_adjust_method, " adjustment for multiple comparisons"
      ),
      individual_calculation = "AUC calculated using trapezoidal rule for each subject",
      growth_rate_calculation = paste0(
        "Growth rates are calculated by fitting a linear regression model to log-transformed volume data over time for each subject ",
        "(non-positive volumes are replaced with half the smallest positive value before the log). ",
        "The slope coefficient from this model represents the exponential growth rate. ",
        "A value of 0.1 indicates approximately 10% tumor volume increase per day. ",
        "Only subjects with 3 or more time points are included in growth rate calculations."
      )
    ),
    notes = c(
      # CODE_REVIEW.md R3.16 — AUC is integrated on the raw volume scale; the
      # `transform` argument does not reach this path.
      if (transform != "none") {
        paste0("AUC was computed on the RAW volume scale. The requested '",
               transform, "' transform applies to the lme4 / gam model paths ",
               "only and was not applied here.")
      } else {
        "No transformation applied to volume data"
      },
      "Composite IDs were created by combining subject ID, treatment group, and cage information to ensure correct AUC values",
      "Welch's t-tests are used for pairwise comparisons to account for potentially unequal variances between treatment groups"
    )
  )

  # Create posthoc object for compatibility with existing code
  posthoc <- list(
    method   = paste0("Welch's t-tests (", comparison_spec$family, ") with ",
                      comparison_spec$p_adjust_method, " adjustment"),
    pairwise = pairwise_df,
    data     = pairwise_data
  )

  # CODE_REVIEW.md DIAGNOSTICS bug (1) — also build ggplot diagnostics for
  # the AUC path so the dashboard's TG Diagnostics tab renders correctly
  # regardless of model_type.
  rd <- if (include_diagnostics)
    build_residual_diagnostic_plots(auc_model,
                                    title_prefix = "Tumour growth AUC")
  else
    list(diag_qq_plot = NULL,
         diag_resid_fitted_plot = NULL,
         diag_scale_location_plot = NULL)

  # Return results for AUC model
  list(
    model                = if (return_model) auc_model else NULL,
    # CODE_REVIEW.md R3.29 — every bayesian_* function returns these; no
    # frequentist path did, so a caller handed a result object had no
    # programmatic way to learn what scale the numbers are on. "none" here is
    # the truth for the AUC path (see R3.16).
    model_type_used      = "auc",
    meta = me_result_meta(
      analysis_type     = "Area under the curve (per-animal trapezoidal)",
      model_type_used   = "auc",
      inference         = "frequentist",
      interval_type     = "confidence",
      # R3.16 -- AUC is integrated on the RAW scale whatever `transform` says.
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
    pairwise_comparisons = posthoc$pairwise,
    treatment_effects    = treatment_effects,
    growth_rates         = growth_rates,
    cage_analysis        = cage_analysis,
    auc_analysis         = auc_analysis,
    data_summary         = data_summary,
    diagnostics          = diagnostics,
    variance_test        = variance_test,
    diag_qq_plot             = rd$diag_qq_plot,
    diag_resid_fitted_plot   = rd$diag_resid_fitted_plot,
    diag_scale_location_plot = rd$diag_scale_location_plot,
    diag_re_qq_plot          = NULL,  # AUC fit is OLS, no random effects
    necrosis_summary     = necrosis_summary
  )
}
