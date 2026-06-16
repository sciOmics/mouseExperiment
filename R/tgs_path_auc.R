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
#' @param transform,p_adjust_method,reference_group See
#'   \code{\link{tumor_growth_statistics}}.
#' @param return_model,include_diagnostics See
#'   \code{\link{tumor_growth_statistics}}.
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
                         p_adjust_method,
                         reference_group,
                         return_model,
                         include_diagnostics,
                         auc_bootstrap_n,
                         auc_bootstrap_seed,
                         cage_analysis,
                         data_summary,
                         necrosis_summary,
                         id_column,
                         treatment_column,
                         cage_column,
                         time_column) {

  # Create ANOVA model for AUC
  auc_model   <- stats::aov(AUC ~ Treatment, data = auc_analysis$individual)
  anova_table <- stats::anova(auc_model)

  # Create pairwise comparisons for AUC using Welch's t-tests
  # This is more appropriate for AUC analysis where variances between groups may differ
  treatments       <- unique(auc_analysis$individual$Treatment)
  pairwise_results <- list()
  pairwise_data    <- list()

  pairs <- utils::combn(treatments, 2, simplify = FALSE)

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
        boot_ci_upper    = NA_real_
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

      pairwise_results[[paste(pair[1], "-", pair[2])]] <- list(
        comparison       = paste(pair[1], "-", pair[2]),
        mean_diff        = mean(group1_data) - mean(group2_data),
        t_value          = t_test_result$statistic,
        df               = t_test_result$parameter,
        p_value          = t_test_result$p.value,
        ci_lower         = t_test_result$conf.int[1],
        ci_upper         = t_test_result$conf.int[2],
        boot_ci_lower    = unname(boot_ci["lower"]),
        boot_ci_upper    = unname(boot_ci["upper"])
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
      stringsAsFactors = FALSE
    )
  }))

  # Apply multiple comparison correction
  pairwise_df$p_adjusted      <- stats::p.adjust(pairwise_df$p_value, method = p_adjust_method)
  pairwise_df$p_adjust_method <- p_adjust_method

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
      volume_transformation  = transform,
      auc_calculation_method = "trapezoidal",
      statistical_test       = "One-way ANOVA on AUC values",
      posthoc_method         = paste0(
        "Welch's t-tests with ", p_adjust_method,
        " adjustment for multiple comparisons"
      ),
      individual_calculation = "AUC calculated using trapezoidal rule for each subject",
      growth_rate_calculation = paste0(
        "Growth rates are calculated by fitting a linear regression model to log1p-transformed volume data over time for each subject. ",
        "The slope coefficient from this model represents the exponential growth rate. ",
        "A value of 0.1 indicates approximately 10% tumor volume increase per day. ",
        "Only subjects with 3 or more time points are included in growth rate calculations."
      )
    ),
    notes = c(
      if (transform != "none") paste("Volume data was", transform, "transformed prior to analysis") else "No transformation applied to volume data",
      "Composite IDs were created by combining subject ID, treatment group, and cage information to ensure correct AUC values",
      "Welch's t-tests are used for pairwise comparisons to account for potentially unequal variances between treatment groups"
    )
  )

  # Create posthoc object for compatibility with existing code
  posthoc <- list(
    method   = paste0("Welch's t-tests with ", p_adjust_method, " adjustment"),
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
    diag_qq_plot             = rd$diag_qq_plot,
    diag_resid_fitted_plot   = rd$diag_resid_fitted_plot,
    diag_scale_location_plot = rd$diag_scale_location_plot,
    diag_re_qq_plot          = NULL,  # AUC fit is OLS, no random effects
    necrosis_summary     = necrosis_summary
  )
}
