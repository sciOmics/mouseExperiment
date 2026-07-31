#' Analyze Drug Combination Synergy in Tumor Growth
#'
#' This function tests for synergistic effects of drug combinations in tumor growth data.
#' It compares the observed combination effect against several expected interaction models,
#' including Bliss independence and Loewe additivity. The function calculates synergy scores
#' and performs statistical tests to determine if the combination shows synergistic,
#' additive, or antagonistic effects.
#'
#' @param df A data frame containing tumor growth data.
#' @param treatment_column A character string specifying the column name for treatment groups. Default is "Treatment".
#' @param volume_column A character string specifying the column name for tumor volume measurements. Default is "Volume".
#' @param time_column A character string specifying the column name for time points. Default is "Day".
#' @param drug_a_name A character string specifying the name of the first single agent treatment group.
#' @param drug_b_name A character string specifying the name of the second single agent treatment group.
#' @param combo_name A character string specifying the name of the combination treatment group.
#' @param control_name A character string specifying the name of the control/vehicle group. Default is "Control".
#' @param eval_time_point Optional. A numeric value specifying a specific time point to evaluate synergy.
#'        If NULL (default), the function will use the last time point in the data.
#' @param verbose Logical. If TRUE, prints detailed results to the console.
#'        Default is TRUE for interactive use; set to FALSE for programmatic/dashboard use.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{summary}{A data frame summarizing the tumor growth inhibition (TGI) for each treatment and synergy metrics.}
#'   \item{bliss_independence}{Results of the Bliss independence model, including expected vs. observed effects.}
#'   \item{loewe_additivity}{Results of the Loewe additivity model (linear dose-response assumption).}
#'   \item{combination_index}{The Loewe-based combination index (CI), where CI < 1 indicates synergy, CI = 1 indicates additivity, and CI > 1 indicates antagonism.}
#'   \item{statistical_test}{Results of statistical tests comparing observed vs. expected effects.}
#'   \item{plot_data}{Data prepared for plotting, to be used with plot_drug_synergy function.}
#' }
#'
#' @details
#' The function calculates tumor growth inhibition (TGI) for each treatment group relative to the control.
#' It then applies several models to test for synergy:
#'
#' 1. Bliss Independence Model: Assumes drugs act independently through different mechanisms.
#'    Expected effect = EA + EB - (EA * EB), where EA and EB are the effects of drug A and B alone.
#'
#' 2. Loewe Additivity Model: Under a linear dose-response assumption (appropriate when only
#'    single-dose TGI data is available), the expected combination effect equals the sum of
#'    individual fractional effects, capped at 1.0. The Combination Index is calculated as
#'    CI = min(FE_A + FE_B, 1) / FE_combo (Berenbaum, 1989).
#'
#' 3. Combination Index: CI < 1 indicates synergy, CI = 1 indicates additivity, CI > 1 indicates antagonism.
#'
#' The function performs statistical tests to determine if the observed combination effect
#' significantly differs from the expected effect under these models.
#'
#' @section Assumptions and Limitations:
#' \strong{Bliss Independence applied to TGI:} Bliss Independence was formulated for the
#' probability of cell death, not for proportional growth inhibition. Applying it to TGI is a
#' common pragmatic choice but carries a ceiling effect: when individual drug TGIs are large
#' (each > 50\%), the Bliss expected combined TGI approaches 100\%, making it nearly impossible
#' to demonstrate synergy by this criterion regardless of the true biological interaction.
#' Interpret Bliss results cautiously when individual-agent TGIs exceed 50\%.
#'
#' \strong{Point estimates and labels:} the \code{synergy_label} and
#' \code{synergy_interpretation} strings are derived from three group means via
#' fixed thresholds. They are descriptive summaries, not test results. Read them
#' alongside \code{synergy_ci} (mouse-level bootstrap 95\% intervals) and
#' \code{group_n}: a "Strong Synergy" label whose \code{Bliss_Excess_FE}
#' interval spans zero is not evidence of synergy. For a model-based posterior
#' treatment of the same question, use \code{\link{bayesian_synergy}}.
#'
#' \strong{Loewe Additivity single-dose approximation:} The CI formula
#' \code{min(FE_A + FE_B, 1) / FE_combo} assumes a linear dose-response relationship, which
#' is appropriate when only single-dose data are available but is an approximation. Without
#' dose-response curves the true IC50 values are unknown, so the result should be interpreted
#' as a qualitative synergy indicator rather than a precise mechanistic estimate.
#'
#' @examples
#' # Example with synthetic dataset
#' data(combo_treatment_synthetic_data)
#' data_processed <- calculate_volume(combo_treatment_synthetic_data)
#' data_processed <- calculate_dates(data_processed, start_date = "03/24/2025")
#' 
#' synergy_results <- analyze_drug_synergy(
#'   df = data_processed,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Combo",
#'   control_name = "Control"
#' )
#' 
#' # Print the summary
#' print(synergy_results$summary)
#' 
#' # Create and display the synergy visualization
#' synergy_plot <- plot_drug_synergy(synergy_results)
#' print(synergy_plot)
#'
#' @import dplyr
#' @import ggplot2
#' @param id_column Column identifying individual animals. Used to resample
#'   mice for the bootstrap; also used to report per-group n.
#' @param endpoint_method How each arm's volume at \code{eval_time_point} is
#'   obtained: "model" (default, log-scale LMM marginal means using every
#'   observation), "last_obs", or "survivors" (pre-0.8.0 behaviour; conditions
#'   on survival and understates TGI). See CODE_REVIEW.md R3.5 / G.3.
#' @param ci_thresholds Length-2 numeric giving the Loewe combination-index band
#'   treated as indistinguishable from additive. Default \code{c(0.85, 1.15)} —
#'   a convention, not a derived quantity (CODE_REVIEW.md R3.28).
#' @param strong_synergy_delta Numeric. Excess fractional effect (Bliss *and*
#'   Loewe) above which the label is "Strong Synergy". Default 0.1, also a
#'   convention.
#' @param n_boot Integer >= 0. Number of bootstrap resamples used to attach
#'   confidence intervals to the synergy metrics. Mice are resampled with
#'   replacement *within* treatment group, including the control arm, and the
#'   whole statistic is recomputed per resample — so the interval propagates
#'   both the sampling error of each arm and the sampling error of the control
#'   mean that forms the TGI denominator (CODE_REVIEW.md R3.6 / R3.7 / G.6).
#'   Default 2000. Set 0 to skip.
#' @param boot_seed Optional integer seed for reproducible resampling.
#' @export
analyze_drug_synergy <- function(df, 
                               treatment_column = "Treatment",
                               volume_column = "Volume",
                               time_column = "Day",
                               drug_a_name,
                               drug_b_name,
                               combo_name,
                               control_name = "Control",
                               eval_time_point = NULL,
                               id_column = "ID",
                               endpoint_method = c("model", "last_obs", "survivors"),
                               ci_thresholds = c(0.85, 1.15),
                               strong_synergy_delta = 0.1,
                               n_boot = 2000L,
                               boot_seed = NULL,
                               verbose = TRUE) {

  endpoint_method <- match.arg(endpoint_method)
  
  # Input validation
  required_columns <- c(treatment_column, volume_column, time_column)
  missing_cols <- required_columns[!required_columns %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check that the specified groups exist in the data
  all_groups <- c(drug_a_name, drug_b_name, combo_name, control_name)
  missing_groups <- all_groups[!all_groups %in% unique(df[[treatment_column]])]
  
  if (length(missing_groups) > 0) {
    stop("The following specified groups do not exist in the treatment column: ", 
         paste(missing_groups, collapse = ", "))
  }
  
  # If no specific time point is provided, use the last time point
  if (is.null(eval_time_point)) {
    eval_time_point <- max(df[[time_column]])
    message(paste("No specific evaluation time point provided. Using the last time point:", eval_time_point))
  } else {
    # Check that the specified time point exists
    if (!eval_time_point %in% df[[time_column]]) {
      closest_time <- df[[time_column]][which.min(abs(df[[time_column]] - eval_time_point))]
      warning(paste("Specified time point", eval_time_point, "not found in data.",
                   "Using closest available time point:", closest_time))
      eval_time_point <- closest_time
    }
  }
  
  # CODE_REVIEW.md R3.5 / G.3 — filtering to rows observed exactly at
  # eval_time_point conditions on survival to that day. Animals leave because
  # their tumours grew, so this selects the slowest growers, hardest in the
  # control arm, and biases every TGI downward. Route through the shared
  # endpoint helper: the default reads each arm's geometric mean at the
  # evaluation day off a log-scale LMM fitted to every observation, and the
  # per-mouse table it returns keeps all animals for the bootstrap.
  ep <- endpoint_volumes(
    df, id_column = id_column, treatment_column = treatment_column,
    time_column = time_column, volume_column = volume_column,
    endpoint_day = eval_time_point, endpoint_method = endpoint_method
  )
  # Per-mouse endpoint volumes (all animals) for the group means and bootstrap.
  analysis_data <- if (!is.null(ep$per_mouse)) {
    data.frame(Treatment = ep$per_mouse$Treatment,
               Volume    = ep$per_mouse$Volume,
               stringsAsFactors = FALSE)
  } else {
    pm <- endpoint_volumes(
      df, id_column = id_column, treatment_column = treatment_column,
      time_column = time_column, volume_column = volume_column,
      endpoint_day = eval_time_point, endpoint_method = "last_obs")$per_mouse
    data.frame(Treatment = pm$Treatment, Volume = pm$Volume,
               stringsAsFactors = FALSE)
  }
  names(analysis_data) <- c(treatment_column, volume_column)

  # CODE_REVIEW.md R3.7 — tapply()'s default na.rm = FALSE meant a single
  # missing volume at the evaluation day silently NA'd that group's mean and
  # every quantity derived from it.
  group_means <- stats::setNames(ep$group_means$Mean_Volume,
                                ep$group_means$Treatment)
  group_n <- stats::setNames(ep$group_means$N, ep$group_means$Treatment)

  # Fail loudly rather than propagating NA when a named arm is absent.
  for (nm in c(control_name, drug_a_name, drug_b_name, combo_name)) {
    if (!nm %in% names(group_means) || !is.finite(group_means[[nm]])) {
      stop("Treatment group '", nm, "' has no usable volume observations at ",
           "day ", eval_time_point, ".", call. = FALSE)
    }
  }
  
  # Extract mean volumes for each group
  control_mean <- group_means[control_name]
  drug_a_mean <- group_means[drug_a_name]
  drug_b_mean <- group_means[drug_b_name]
  combo_mean <- group_means[combo_name]
  
  # Calculate Tumor Growth Inhibition (TGI) for each treatment
  tgi_a <- 100 * (1 - (drug_a_mean / control_mean))
  tgi_b <- 100 * (1 - (drug_b_mean / control_mean))
  tgi_combo <- 100 * (1 - (combo_mean / control_mean))
  
  # Convert TGI to fractional effect (FE)
  fe_a <- tgi_a / 100
  fe_b <- tgi_b / 100
  fe_combo <- tgi_combo / 100
  
  # CODE_REVIEW.md R3.6 / R3.7 / G.6 — attach mouse-level bootstrap intervals to
  # every synergy metric. Resampling the control arm too propagates the TGI
  # denominator's uncertainty, which was previously treated as a known constant.
  per_mouse_vols <- list(
    control = as.numeric(analysis_data[[volume_column]][
      analysis_data[[treatment_column]] == control_name]),
    drug_a  = as.numeric(analysis_data[[volume_column]][
      analysis_data[[treatment_column]] == drug_a_name]),
    drug_b  = as.numeric(analysis_data[[volume_column]][
      analysis_data[[treatment_column]] == drug_b_name]),
    combo   = as.numeric(analysis_data[[volume_column]][
      analysis_data[[treatment_column]] == combo_name])
  )
  synergy_ci <- if (n_boot > 0L) {
    synergy_bootstrap(per_mouse_vols, n_boot = as.integer(n_boot),
                      seed = boot_seed)
  } else NULL

  # Calculate expected effect using Bliss Independence model.
  # Shared scalar/vector formula in R/utils_synergy.R.
  bliss_expected_fe <- synergy_bliss_expected(fe_a, fe_b)
  bliss_expected_tgi <- bliss_expected_fe * 100

  # Calculate expected effect using Loewe Additivity model
  # Under a linear dose-response assumption (Berenbaum, 1989), the Loewe expected
  # combination effect equals the sum of individual fractional effects.
  # Capped at 1.0 since fractional effects cannot exceed 100% inhibition.
  loewe_expected_fe <- min(fe_a + fe_b, 1.0)
  loewe_expected_tgi <- loewe_expected_fe * 100
  
  # R14.2: a single agent that INCREASES tumour burden relative to control gives a
  # negative fractional effect, and the ratio below then goes negative. Because
  # the verdict thresholds are one-sided ("CI < 0.85 => synergistic"), any
  # negative CI was classified as synergy. A worked case: DrugA accelerating
  # growth, DrugB inert and the combination roughly matching control produced
  # CI = -29.2 reported as "Synergistic (CI < 0.85)" — the verdict exactly
  # inverted, on the arm configuration a reviewer would scrutinise hardest.
  #
  # Both synergy models are defined for *inhibitory* agents: Bliss multiplies
  # surviving fractions, Loewe adds dose-equivalents. Neither has a meaning when
  # an arm's effect is negative, so the honest answer is that the analysis does
  # not apply, not a number pushed through a threshold. The Toxicity path already
  # takes this view — TWM uses max(TGI, 0) — so this also removes an
  # inconsistency between the two.
  agents_inhibitory <- is.finite(fe_a) && is.finite(fe_b) && fe_a > 0 && fe_b > 0
  if (!agents_inhibitory) {
    ci_value <- NA_real_
    warning("Combination Index not computed: a single-agent arm did not inhibit ",
            "growth relative to control (fractional effect <= 0). Bliss and Loewe ",
            "are both defined for inhibitory agents; a synergy verdict here would ",
            "be meaningless.", call. = FALSE)
  } else if (fe_combo > 0) {
    ci_value <- loewe_expected_fe / fe_combo
  } else {
    ci_value <- NA_real_
    warning("Cannot calculate Combination Index: Combination effect is zero or negative",
            call. = FALSE)
  }
  
  # Determine synergy interpretation based on CI.
  # CODE_REVIEW.md R3.28 — the 0.85 / 1.15 band was hardcoded and undocumented.
  # It is a convention, not a derived quantity: the interval is a pragmatic
  # "indistinguishable from additive" zone, and different labs use different
  # widths. It is now an argument so the choice is explicit and tunable.
  lo <- min(ci_thresholds); hi <- max(ci_thresholds)
  if (is.na(ci_value)) {
    synergy_interpretation <- if (!agents_inhibitory) {
      "Not evaluable - a single agent did not inhibit growth relative to control"
    } else {
      "Cannot determine (CI calculation error)"
    }
  } else if (ci_value < lo) {
    synergy_interpretation <- paste0("Synergistic (CI < ", lo, ")")
  } else if (ci_value <= hi) {
    synergy_interpretation <- paste0("Additive (", lo, " <= CI <= ", hi, ")")
  } else {
    synergy_interpretation <- paste0("Antagonistic (CI > ", hi, ")")
  }
  
  # Calculate the difference between observed and expected effects
  bliss_difference <- fe_combo - bliss_expected_fe
  loewe_difference <- fe_combo - loewe_expected_fe
  
  # Determine synergy label based on both models.
  # CODE_REVIEW.md R3.28 — `strong_synergy_delta` replaces the hardcoded 0.1.
  # NOTE these labels are computed from three group means with no uncertainty
  # attached; see the `synergy_ci` element and the Assumptions section for why
  # a label alone should not drive a decision.
  d <- strong_synergy_delta
  if (!agents_inhibitory) {
    # Same reasoning as the CI guard above: the Bliss/Loewe differences are large
    # and positive here precisely *because* a single agent did harm, so the label
    # logic would read the harm as evidence of synergy.
    synergy_label <- "Not evaluable (a single agent did not inhibit growth)"
  } else if (bliss_difference > d && loewe_difference > d) {
    synergy_label <- "Strong Synergy"
  } else if (bliss_difference > 0 && loewe_difference > 0) {
    synergy_label <- "Synergy"
  } else if (bliss_difference > -d && loewe_difference > -d) {
    synergy_label <- "Additivity"
  } else {
    synergy_label <- "Antagonism"
  }
  
  # Statistical test for synergy
  # Perform t-test to compare combo group vs. each single agent
  t_test_a_combo <- t.test(
    analysis_data[[volume_column]][analysis_data[[treatment_column]] == combo_name],
    analysis_data[[volume_column]][analysis_data[[treatment_column]] == drug_a_name]
  )
  
  t_test_b_combo <- t.test(
    analysis_data[[volume_column]][analysis_data[[treatment_column]] == combo_name],
    analysis_data[[volume_column]][analysis_data[[treatment_column]] == drug_b_name]
  )
  
  # Create a data frame for summary results
  summary_df <- data.frame(
    Treatment = c(drug_a_name, drug_b_name, combo_name, "Bliss Expected", "Loewe Expected"),
    Mean_Volume = c(drug_a_mean, drug_b_mean, combo_mean, 
                   control_mean * (1 - bliss_expected_fe), 
                   control_mean * (1 - loewe_expected_fe)),
    TGI_Percent = c(tgi_a, tgi_b, tgi_combo, bliss_expected_tgi, loewe_expected_tgi),
    Fractional_Effect = c(fe_a, fe_b, fe_combo, bliss_expected_fe, loewe_expected_fe)
  )
  
  # Add synergy metrics to a separate data frame
  synergy_metrics <- data.frame(
    Metric = c("Bliss Difference", "Loewe Difference", "Combination Index", "Interpretation"),
    Value = c(bliss_difference, loewe_difference, ci_value, synergy_interpretation)
  )
  
  # Create statistical test summary
  stat_tests <- data.frame(
    Comparison = c(paste("Combo vs", drug_a_name), paste("Combo vs", drug_b_name)),
    P_Value = c(t_test_a_combo$p.value, t_test_b_combo$p.value),
    Significant = c(t_test_a_combo$p.value < 0.05, t_test_b_combo$p.value < 0.05)
  )
  
  # Prepare data for plotting (will be used by plot_drug_synergy function)
  plot_data <- data.frame(
    Treatment = factor(c(drug_a_name, drug_b_name, combo_name, "Bliss Expected", "Loewe Expected"),
                     levels = c(drug_a_name, drug_b_name, "Bliss Expected", "Loewe Expected", combo_name)),
    TGI = c(tgi_a, tgi_b, bliss_expected_tgi, loewe_expected_tgi, tgi_combo),
    Type = c("Observed", "Observed", "Expected", "Expected", "Observed")
  )
  
  # Print results (only when verbose)
  if (isTRUE(verbose)) {
    message("\n=== Drug Combination Synergy Analysis ===")
    message("Evaluation time point: ", eval_time_point, "\n")
    
    message("Treatment Mean Volumes:")
    message("Control (", control_name, "): ", round(control_mean, 2))
    message(drug_a_name, ": ", round(drug_a_mean, 2), " (TGI: ", round(tgi_a, 1), "%)")
    message(drug_b_name, ": ", round(drug_b_mean, 2), " (TGI: ", round(tgi_b, 1), "%)")
    message(combo_name, ": ", round(combo_mean, 2), " (TGI: ", round(tgi_combo, 1), "%)\n")
    
    message("Expected Effects:")
    message("Bliss Independence: TGI = ", round(bliss_expected_tgi, 1), "%")
    message("Loewe Additivity: TGI = ", round(loewe_expected_tgi, 1), "%\n")
    
    message("Synergy Assessment:")
    message("Bliss Difference: ", round(bliss_difference * 100, 1), "% (", 
               ifelse(bliss_difference > 0, "Synergy", "No Synergy"), ")")
    message("Loewe Difference: ", round(loewe_difference * 100, 1), "% (", 
               ifelse(loewe_difference > 0, "Synergy", "No Synergy"), ")")
    message("Combination Index: ", round(ci_value, 2), " (", synergy_interpretation, ")\n")
    
    message("Statistical Tests:")
    message("Combo vs ", drug_a_name, ": p = ", round(t_test_a_combo$p.value, 4), 
               ifelse(t_test_a_combo$p.value < 0.05, " (Significant)", " (Not Significant)"))
    message("Combo vs ", drug_b_name, ": p = ", round(t_test_b_combo$p.value, 4), 
               ifelse(t_test_b_combo$p.value < 0.05, " (Significant)", " (Not Significant)"), "\n")
    
    message("Overall Assessment: ", synergy_label, "\n")
  }
  
  # CODE_REVIEW.md DIAGNOSTICS gap (6) — frequentist synergy had no
  # checkable model-assumption diagnostics for the t-test path. Surface a
  # per-group Q-Q plot of per-mouse volumes (Welch t robustness check)
  # and a boxplot of per-mouse volumes (outlier check).
  diag_group_qq_plot <- NULL
  diag_group_boxplot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    diag_group_qq_plot <- tryCatch({
      ad <- analysis_data
      # Build a Q-Q per group via stratified qqnorm computations.
      grp_col <- treatment_column
      val_col <- volume_column
      groups <- unique(ad[[grp_col]])
      long <- do.call(rbind, lapply(groups, function(g) {
        v <- ad[[val_col]][ad[[grp_col]] == g]
        v <- v[is.finite(v)]
        if (length(v) < 3L) return(NULL)
        qq <- stats::qqnorm(v, plot.it = FALSE)
        data.frame(group = g, theoretical = qq$x, sample = qq$y,
                   stringsAsFactors = FALSE)
      }))
      if (is.null(long) || nrow(long) == 0L) NULL
      else ggplot2::ggplot(long, ggplot2::aes(x = .data[["theoretical"]],
                                              y = .data[["sample"]])) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::geom_smooth(method = "lm", se = FALSE,
                             colour = "red", linetype = "dashed",
                             formula = y ~ x) +
        ggplot2::facet_wrap(~ group, scales = "free") +
        ggplot2::theme_classic() +
        ggplot2::labs(
          title = "Per-group Q-Q of tumor volumes",
          subtitle = paste("At evaluation day", eval_time_point,
                           "— heavy tails or curvature warn t-test assumptions"),
          x = "Theoretical Quantiles", y = "Sample Quantiles")
    }, error = function(e) NULL)

    diag_group_boxplot <- tryCatch({
      ggplot2::ggplot(analysis_data,
                      ggplot2::aes(x = .data[[treatment_column]],
                                   y = .data[[volume_column]])) +
        ggplot2::geom_boxplot(outlier.colour = "red", outlier.alpha = 0.8) +
        ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.4) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          title = "Per-group tumor volumes at evaluation day",
          subtitle = paste("Day", eval_time_point,
                           "— outliers in red. Inspect before reporting synergy."),
          x = "Treatment", y = "Tumor volume")
    }, error = function(e) NULL)
  }

  # Return a list with all results
  return(list(
    summary = summary_df,
    synergy_metrics = synergy_metrics,
    # CODE_REVIEW.md R3.6 / R3.7 — bootstrap percentile CIs for every metric,
    # plus per-group n so a reader can weigh the point estimates.
    synergy_ci  = synergy_ci,
    group_n     = group_n,
    attrition   = ep$attrition,
    endpoint_method = ep$method,
    eval_time_point = eval_time_point,
    thresholds  = list(ci_thresholds = ci_thresholds,
                       strong_synergy_delta = strong_synergy_delta),
    bliss_independence = list(
      expected_effect = bliss_expected_fe,
      observed_effect = fe_combo,
      difference = bliss_difference,
      synergy = bliss_difference > 0
    ),
    loewe_additivity = list(
      expected_effect = loewe_expected_fe,
      observed_effect = fe_combo,
      difference = loewe_difference,
      synergy = loewe_difference > 0
    ),
    combination_index = list(
      ci = ci_value,
      interpretation = synergy_interpretation
    ),
    statistical_tests = stat_tests,
    overall_assessment = synergy_label,
    evaluation_time_point = eval_time_point,
    plot_data = plot_data, # Keep the plot data for later plotting
    # Per-group diagnostics (Welch t robustness check) — CODE_REVIEW.md
    # DIAGNOSTICS gap (6).
    diag_group_qq_plot = diag_group_qq_plot,
    diag_group_boxplot = diag_group_boxplot,
    # Additional data needed for plotting
    drug_a_name = drug_a_name,
    drug_b_name = drug_b_name,
    combo_name = combo_name,
    control_name = control_name
  ))
}

#' Plot Drug Combination Synergy Analysis
#'
#' Creates a bar plot visualizing tumor growth inhibition (TGI) for different treatment groups
#' and synergy metrics from a drug combination analysis.
#'
#' @param synergy_results Results object from analyze_drug_synergy function
#' @param custom_title Optional custom title for the plot
#' @param custom_colors Optional named vector of custom colors for plot elements
#'
#' @return A ggplot2 object visualizing the synergy analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # First run the analysis
#' results <- analyze_drug_synergy(
#'   df = tumor_data,
#'   drug_a_name = "Drug A",
#'   drug_b_name = "Drug B", 
#'   combo_name = "Drug A + Drug B",
#'   control_name = "Vehicle"
#' )
#' 
#' # Then create the plot
#' plot_drug_synergy(results)
#' 
#' # With custom title
#' plot_drug_synergy(results, custom_title = "Custom Analysis Title")
#' }
plot_drug_synergy <- function(synergy_results, custom_title = NULL, custom_colors = NULL) {
  # Validate input
  if (!is.list(synergy_results) || is.null(synergy_results$plot_data)) {
    stop("Input must be a valid result object from analyze_drug_synergy()")
  }
  
  # Extract the plot data
  plot_data <- synergy_results$plot_data
  
  # Set title
  if (is.null(custom_title)) {
    title <- paste("Drug Combination Analysis at Day", synergy_results$evaluation_time_point)
  } else {
    title <- custom_title
  }
  
  # Set colors
  if (is.null(custom_colors)) {
    fill_colors <- c("Expected" = "lightblue", "Observed" = "darkblue")
  } else {
    fill_colors <- custom_colors
  }
  
  # Create the plot
  synergy_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Treatment, y = TGI, fill = Type)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "black") +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::labs(
      title = title,
      subtitle = paste("Synergy Assessment:", synergy_results$overall_assessment, 
                     "(CI =", round(synergy_results$combination_index$ci, 2), ")"),
      x = "Treatment",
      y = "Tumor Growth Inhibition (%)",
      fill = "Data Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(ggplot2::aes(label = round(TGI, 1)), vjust = -0.5)
  
  return(synergy_plot)
}