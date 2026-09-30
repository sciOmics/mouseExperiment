#' Analyze Drug Combination Synergy in Tumor Growth
#'
#' This function tests for synergistic effects of drug combinations in tumor growth data.
#' It compares the observed combination effect against several expected interaction models,
#' using Bliss independence. The function calculates synergy scores
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
#' @param eval_time_point Day at which to evaluate synergy. \code{NULL}
#'   (default): the last day on which all four arms are evaluable, i.e. at
#'   least 50 % of each arm's enrolled animals, and at least 3, are still on
#'   study (see \code{\link{evaluable_days}}). A requested day must be
#'   evaluable; one that was not measured moves to the closest measured day.
#'   Before v0.26.0 the default was the last day in the data, where the control
#'   arm had usually left the study and its mean was an extrapolation
#'   (CODE_REVIEW.md R20.1).
#' @param verbose Logical. If TRUE, prints detailed results to the console.
#'        Default is TRUE for interactive use; set to FALSE for programmatic/dashboard use.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{summary}{A data frame summarizing the tumor growth inhibition (TGI) for each treatment and synergy metrics.}
#'   \item{bliss_independence}{Results of the Bliss independence model, including expected vs. observed effects.
#'     When \code{bliss_applies} is FALSE, \code{expected_effect}, \code{difference} and \code{synergy}
#'     are NA.}
#'   \item{bliss_applies}{FALSE when a single agent did not inhibit growth relative to control.
#'     (Named \code{evaluable} before v0.26.0; renamed so that "evaluable" refers only to
#'     evaluable days.)
#'     Bliss independence is defined for inhibitory agents, so every Bliss quantity is then NA:
#'     the expectation and difference, the \code{synergy} flag, the "Bliss Expected" row of
#'     \code{summary} and the \code{Bliss_Excess_FE} interval in \code{synergy_ci}
#'     (CODE_REVIEW.md R14.2, R20.8).}
#'   \item{statistical_tests}{The combination against each single agent at the
#'     evaluation day, on the same estimand as the point estimates: under
#'     "model", a Wald test of the log volume ratio from the endpoint model;
#'     otherwise Welch t-tests on the per-animal volumes.}
#'   \item{synergy_ci, interval_method}{95 % intervals for the TGIs and the
#'     Bliss excess, and how they were obtained.}
#'   \item{evaluability}{The evaluable-day record: the rule, the day used, and
#'     the days and arms it excluded (see \code{\link{evaluable_days}}).}
#'   \item{endpoint_model}{How the endpoint model was fitted: time basis,
#'     random effects (and any fallback), excluded pre-palpable zeros.}
#'   \item{diag_group_qq_plot, diag_group_boxplot, diag_fit_plot}{Diagnostics
#'     for the estimand used: under "model", the model's residual Q-Q by arm
#'     and each arm's fitted curve over its data; otherwise per-animal Q-Q and
#'     box plots.}
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
#' \strong{Why there is no Combination Index.} A Loewe-style CI was removed in
#' v0.21.0. What it computed was not Loewe additivity: true Loewe is
#' dose-equivalence, \eqn{CI = d_A/D_A + d_B/D_B}, which needs a dose-response
#' curve per agent that single-dose designs do not provide. The stand-in,
#' \code{min(FE_A + FE_B, 1) / FE_combo}, is \emph{response additivity} -- a
#' different null that fails the sham-combination test: combining a drug with
#' itself it predicts twice the fractional effect, so an agent at FE = 0.5 should
#' reach 100 percent inhibition, and when it does not the method calls the drug
#' antagonistic with itself. Across the (FE_A, FE_B) grid with the combination
#' set to exactly Bliss-additive, 42 percent of cells were labelled antagonistic
#' and none synergistic. For a real Loewe analysis collect per-agent
#' dose-response curves and use \code{drc::isobole()}.
#'
#' The function performs statistical tests to determine if the observed combination effect
#' significantly differs from the expected effect under these models.
#'
#' @section Assumptions and Limitations:
#' \strong{Bliss Independence applied to TGI:} Bliss Independence was formulated for the
#' probability of cell death, not for proportional growth inhibition. Applying it to TGI is a
#' common pragmatic choice but carries a ceiling effect: when individual drug TGIs are large
#' (each > 50%), the Bliss expected combined TGI approaches 100%, making it nearly impossible
#' to demonstrate synergy by this criterion regardless of the true biological interaction.
#' Interpret Bliss results cautiously when individual-agent TGIs exceed 50%.
#'
#' \strong{Point estimates and labels:} \code{synergy_label} is derived from
#' fixed thresholds. They are descriptive summaries, not test results. Read them
#' alongside \code{synergy_ci} (95% intervals for the reported estimates) and
#' \code{group_n}: a "Strong Synergy" label whose \code{Bliss_Excess_FE}
#' interval spans zero is not evidence of synergy.
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
#' @param cage_column Optional cage column. Part of the animal key
#'   (treatment + ID + cage), so an ID reused in different cages of one arm
#'   is counted as different animals (CODE_REVIEW.md T1). NULL or absent means
#'   no cage information.
#' @param endpoint_method How each arm's volume at \code{eval_time_point} is
#'   obtained: "model" (default: each arm's geometric mean from a mixed model
#'   of log volume fitted to every observation, with a natural spline in time
#'   per arm and per-animal random slopes), "last_obs", or "survivors"
#'   (pre-0.8.0 behaviour; conditions on survival and understates TGI). See
#'   CODE_REVIEW.md R3.5 / G.3 and R20.1.
#' @param strong_synergy_delta Numeric. Bliss excess fractional effect above
#'   which the label is "Strong Synergy". Default 0.1, a convention rather than
#'   a derived quantity.
#' @param ci_thresholds Deprecated and ignored; retained so existing calls do
#'   not error. It configured the removed Combination Index band.
#' @param n_boot Integer >= 0. Number of draws behind the 95 % intervals in
#'   \code{synergy_ci}. Under \code{endpoint_method = "model"} they are draws
#'   of the endpoint model's fixed effects, so the intervals describe the
#'   reported model-based estimates; under the per-animal estimands they are
#'   bootstrap resamples of animals within arm, control included
#'   (CODE_REVIEW.md R3.6 / R3.7 / G.6, R20.2). Either way the point estimate
#'   is reported, not the median of the draws. Default 2000. Set 0 to skip.
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
                               cage_column = NULL,
                               endpoint_method = c("model", "last_obs", "survivors"),
                               ci_thresholds = c(0.85, 1.15),
                               strong_synergy_delta = 0.1,
                               n_boot = 2000L,
                               boot_seed = NULL,
                               verbose = TRUE) {

  endpoint_method <- match.arg(endpoint_method)

  # Input validation. The ID column is required (R20.15): without it every
  # animal in an arm shared one key and each arm became a single "mouse".
  required_columns <- c(treatment_column, volume_column, time_column, id_column)
  missing_cols <- required_columns[!required_columns %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the data frame: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.null(cage_column) && !cage_column %in% colnames(df)) cage_column <- NULL

  # Check that the specified groups exist in the data
  all_groups <- c(drug_a_name, drug_b_name, combo_name, control_name)
  missing_groups <- all_groups[!all_groups %in% unique(df[[treatment_column]])]
  if (length(missing_groups) > 0) {
    stop("The following specified groups do not exist in the treatment column: ",
         paste(missing_groups, collapse = ", "))
  }

  # CODE_REVIEW.md R20.1 / R20-K -- the default evaluation day was the last day
  # in the data, where the control arm (removed first, at the volume limit) was
  # an extrapolation past its own animals. The default is now the last day on
  # which all four arms are evaluable (>= 50 % and >= 3 animals on study), and
  # a requested day must be evaluable. R3.5 / G.3: the default estimand uses
  # every observation of every animal, not the survivors on that day.
  ep <- endpoint_volumes(
    df, id_column = id_column, treatment_column = treatment_column,
    time_column = time_column, volume_column = volume_column,
    cage_column = cage_column,
    endpoint_day = eval_time_point, endpoint_method = endpoint_method,
    arms = c(control_name, drug_a_name, drug_b_name, combo_name)
  )
  if (is.null(eval_time_point) && isTRUE(verbose)) {
    message("Evaluating at day ", ep$endpoint_day,
            ", the last day on which all four arms are evaluable.")
  }

  me_synergy_from_endpoint(
    ep, control_name = control_name, drug_a_name = drug_a_name,
    drug_b_name = drug_b_name, combo_name = combo_name,
    strong_synergy_delta = strong_synergy_delta,
    n_boot = n_boot, boot_seed = boot_seed, verbose = verbose)
}

#' Synergy quantities from one endpoint evaluation
#'
#' The core of [analyze_drug_synergy()], shared with
#' [analyze_drug_synergy_over_time()], which reuses one endpoint model across
#' days.
#'
#' CODE_REVIEW.md R20.2: the intervals, the tests and the diagnostics describe
#' the estimate that is reported. Under the model estimand they come from the
#' endpoint model (draws of its fixed effects; log-ratio contrasts; its
#' residuals); under the per-animal estimands from those animals (bootstrap;
#' Welch t-tests; per-arm Q-Q). Before, the model-based point estimates were
#' paired with a bootstrap and t-tests of each animal's last observation, so a
#' point estimate could sit outside its own interval.
#'
#' @param ep Output of `endpoint_volumes()` for the four arms.
#' @return The result list of [analyze_drug_synergy()].
#' @noRd
#' @keywords internal
me_synergy_from_endpoint <- function(ep, control_name, drug_a_name, drug_b_name,
                                     combo_name, strong_synergy_delta = 0.1,
                                     n_boot = 2000L, boot_seed = NULL,
                                     verbose = FALSE, warn_not_evaluable = TRUE) {
  eval_time_point <- ep$endpoint_day
  arms <- c(control_name, drug_a_name, drug_b_name, combo_name)
  model_based <- identical(ep$method, "model")

  # CODE_REVIEW.md R3.7 -- tapply()'s default na.rm = FALSE meant a single
  # missing volume at the evaluation day silently NA'd that group's mean and
  # every quantity derived from it.
  group_means <- stats::setNames(ep$group_means$Mean_Volume,
                                ep$group_means$Treatment)
  group_n <- stats::setNames(ep$group_means$N, ep$group_means$Treatment)

  # Fail loudly rather than propagating NA when a named arm is absent.
  for (nm in arms) {
    if (!nm %in% names(group_means) || !is.finite(group_means[[nm]])) {
      stop("Treatment group '", nm, "' has no usable volume observations at ",
           "day ", eval_time_point, ".", call. = FALSE)
    }
  }

  control_mean <- group_means[[control_name]]
  drug_a_mean  <- group_means[[drug_a_name]]
  drug_b_mean  <- group_means[[drug_b_name]]
  combo_mean   <- group_means[[combo_name]]

  # Tumor Growth Inhibition (TGI) and fractional effect (FE) per arm
  tgi_a     <- 100 * (1 - drug_a_mean / control_mean)
  tgi_b     <- 100 * (1 - drug_b_mean / control_mean)
  tgi_combo <- 100 * (1 - combo_mean  / control_mean)
  fe_a <- tgi_a / 100
  fe_b <- tgi_b / 100
  fe_combo <- tgi_combo / 100

  # Bliss independence (R/utils_synergy.R).
  bliss_expected_fe  <- synergy_bliss_expected(fe_a, fe_b)
  bliss_expected_tgi <- bliss_expected_fe * 100
  bliss_difference   <- fe_combo - bliss_expected_fe
  point <- c(TGI_A_pct = tgi_a, TGI_B_pct = tgi_b, TGI_Combo_pct = tgi_combo,
             Bliss_Excess_FE = bliss_difference)

  # Intervals for the reported estimate (R20.2). The point estimate is
  # reported, not the median of the draws.
  per_mouse_vols <- NULL
  if (!model_based) {
    pm <- ep$per_mouse
    per_mouse_vols <- list(
      control = pm$Volume[pm$Treatment == control_name],
      drug_a  = pm$Volume[pm$Treatment == drug_a_name],
      drug_b  = pm$Volume[pm$Treatment == drug_b_name],
      combo   = pm$Volume[pm$Treatment == combo_name])
  }
  synergy_ci <- if (n_boot > 0L) {
    if (model_based) {
      synergy_model_ci(ep$model, arms, eval_time_point, point,
                       n_draws = as.integer(n_boot), seed = boot_seed)
    } else {
      ci <- synergy_bootstrap(per_mouse_vols, n_boot = as.integer(n_boot),
                              seed = boot_seed)
      if (!is.null(ci)) ci$Estimate <- unname(point[ci$Metric])
      ci
    }
  } else NULL

  # R17.2: the Loewe / Combination Index path was removed. What it computed was
  # response additivity, which fails the sham-combination test (a drug combined
  # with itself looks antagonistic). Bliss independence is retained; its excess
  # carries an interval, so the verdict can be read against uncertainty.

  # R14.2: an agent that ACCELERATES growth has a negative fractional effect.
  # Bliss multiplies surviving fractions, so it is defined for inhibitory
  # agents; at negative effect the excess goes large and positive precisely
  # *because* an agent did harm. The honest output is that it does not apply.
  agents_inhibitory <- is.finite(fe_a) && is.finite(fe_b) && fe_a > 0 && fe_b > 0
  if (!agents_inhibitory && warn_not_evaluable) {
    warning("Synergy not evaluated: a single-agent arm did not inhibit growth ",
            "relative to control (fractional effect <= 0). Bliss independence is ",
            "defined for inhibitory agents; a synergy verdict here would be ",
            "meaningless.", call. = FALSE)
  }

  # CODE_REVIEW.md R20.8 -- when Bliss does not apply, none of its quantities
  # are reported.
  if (!agents_inhibitory) {
    bliss_expected_fe  <- NA_real_
    bliss_expected_tgi <- NA_real_
    bliss_difference   <- NA_real_
    if (!is.null(synergy_ci)) {
      bliss_rows <- grepl("^Bliss", synergy_ci$Metric)
      synergy_ci[bliss_rows, c("Estimate", "CI_Lower", "CI_Upper")] <- NA_real_
    }
  }

  # Synergy label from Bliss alone. CODE_REVIEW.md R3.28: the threshold is an
  # argument. The label is descriptive; the interval on Bliss_Excess_FE is the
  # quantity to read.
  d <- strong_synergy_delta
  synergy_label <- if (!agents_inhibitory) {
    "Not evaluable (a single agent did not inhibit growth)"
  } else if (bliss_difference > d) {
    "Strong Synergy"
  } else if (bliss_difference > 0) {
    "Synergy"
  } else if (bliss_difference > -d) {
    "Additivity"
  } else {
    "Antagonism"
  }

  # Combination vs each single agent, on the same estimand (R20.2).
  stat_tests <- if (model_based) {
    synergy_model_tests(ep$model, eval_time_point, combo_name,
                        c(drug_a_name, drug_b_name))
  } else {
    do.call(rbind, lapply(c(drug_a_name, drug_b_name), function(g) {
      tt <- stats::t.test(per_mouse_vols$combo,
                          ep$per_mouse$Volume[ep$per_mouse$Treatment == g])
      data.frame(Comparison = paste("Combo vs", g), P_Value = tt$p.value,
                 Significant = tt$p.value < 0.05,
                 Method = "Welch t-test on per-animal volumes",
                 stringsAsFactors = FALSE)
    }))
  }

  summary_df <- data.frame(
    Treatment = c(drug_a_name, drug_b_name, combo_name, "Bliss Expected"),
    Mean_Volume = c(drug_a_mean, drug_b_mean, combo_mean,
                    control_mean * (1 - bliss_expected_fe)),
    TGI_Percent = c(tgi_a, tgi_b, tgi_combo, bliss_expected_tgi),
    Fractional_Effect = c(fe_a, fe_b, fe_combo, bliss_expected_fe)
  )

  synergy_metrics <- data.frame(
    Metric = c("Bliss Difference", "Interpretation"),
    Value = c(bliss_difference, synergy_label)
  )

  plot_data <- data.frame(
    Treatment = factor(c(drug_a_name, drug_b_name, combo_name, "Bliss Expected"),
                       levels = c(drug_a_name, drug_b_name, "Bliss Expected", combo_name)),
    TGI = c(tgi_a, tgi_b, tgi_combo, bliss_expected_tgi),
    Type = c("Observed", "Observed", "Observed", "Expected")
  )

  if (isTRUE(verbose)) {
    message("\n=== Drug Combination Synergy Analysis ===")
    message("Evaluation day: ", eval_time_point, " (", ep$method, " estimand)\n")
    message("Control (", control_name, "): ", round(control_mean, 2))
    message(drug_a_name, ": ", round(drug_a_mean, 2), " (TGI: ", round(tgi_a, 1), "%)")
    message(drug_b_name, ": ", round(drug_b_mean, 2), " (TGI: ", round(tgi_b, 1), "%)")
    message(combo_name, ": ", round(combo_mean, 2), " (TGI: ", round(tgi_combo, 1), "%)\n")
    message("Bliss Independence: TGI = ",
            if (agents_inhibitory) paste0(round(bliss_expected_tgi, 1), "%") else "not evaluable")
    message("Bliss Difference: ",
            if (agents_inhibitory) paste0(round(bliss_difference * 100, 1), "%") else "not evaluable")
    message("Overall: ", synergy_label, "\n")
    for (i in seq_len(nrow(stat_tests))) {
      message(stat_tests$Comparison[i], ": p = ", signif(stat_tests$P_Value[i], 3))
    }
  }

  diag <- if (model_based) {
    synergy_model_diagnostics(ep$model, arms, eval_time_point)
  } else {
    synergy_permouse_diagnostics(ep$per_mouse, eval_time_point)
  }

  list(
    summary = summary_df,
    synergy_metrics = synergy_metrics,
    # Intervals for the reported estimate: model draws or an animal bootstrap,
    # per `interval_method` (R20.2).
    synergy_ci  = synergy_ci,
    interval_method = if (model_based) "draws from the endpoint model's fixed effects"
                      else "bootstrap of animals within arm",
    group_n     = group_n,
    attrition   = ep$attrition,
    endpoint_method = ep$method,
    eval_time_point = eval_time_point,
    # The evaluable-day record: the rule, the day used, and the days and arms
    # it excluded (R20-K).
    evaluability = ep$evaluability,
    endpoint_model = me_endpoint_model_info(ep$model),
    thresholds  = list(strong_synergy_delta = strong_synergy_delta),
    bliss_independence = list(
      expected_effect = bliss_expected_fe,
      observed_effect = fe_combo,
      difference = bliss_difference,
      # NA, not FALSE, when Bliss does not apply (R20.8).
      synergy = if (agents_inhibitory) bliss_difference > 0 else NA
    ),
    bliss_applies = agents_inhibitory,
    statistical_tests = stat_tests,
    overall_assessment = synergy_label,
    evaluation_time_point = eval_time_point,
    plot_data = plot_data,
    diag_group_qq_plot = diag$qq,
    diag_group_boxplot = diag$box,
    diag_fit_plot      = diag$fit,
    drug_a_name = drug_a_name,
    drug_b_name = drug_b_name,
    combo_name = combo_name,
    control_name = control_name
  )
}

#' Intervals for the synergy metrics from the endpoint model
#' @noRd
#' @keywords internal
synergy_model_ci <- function(em, arms, t, point, n_draws = 2000L, seed = NULL) {
  if (is.null(em) || n_draws < 2L) return(NULL)
  lm_ <- me_draw_logmeans(em, me_beta_draws(em, n_draws, seed), arms, t)
  ctrl <- lm_[, 1L]
  fe <- 1 - exp(lm_[, 2:4, drop = FALSE] - ctrl)
  bliss <- fe[, 3L] - synergy_bliss_expected(fe[, 1L], fe[, 2L])
  draws <- cbind(fe * 100, bliss)
  metrics <- c("TGI_A_pct", "TGI_B_pct", "TGI_Combo_pct", "Bliss_Excess_FE")
  do.call(rbind, lapply(seq_along(metrics), function(i) {
    q <- stats::quantile(draws[, i], c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    data.frame(Metric = metrics[i], Estimate = unname(point[[metrics[i]]]),
               CI_Lower = q[1], CI_Upper = q[2], n_boot = nrow(draws),
               stringsAsFactors = FALSE)
  }))
}

#' Combination vs single agents on the endpoint model's log scale
#' @noRd
#' @keywords internal
synergy_model_tests <- function(em, t, combo, others) {
  xc <- me_endpoint_X(em, combo, t)
  z <- stats::qnorm(0.975)
  do.call(rbind, lapply(others, function(g) {
    cv  <- xc - me_endpoint_X(em, g, t)
    est <- as.numeric(cv %*% em$beta)
    se  <- sqrt(as.numeric(cv %*% em$V %*% t(cv)))
    p   <- 2 * stats::pnorm(-abs(est / se))
    data.frame(Comparison = paste("Combo vs", g), P_Value = p,
               Significant = p < 0.05,
               Volume_Ratio = exp(est), CI_Lower = exp(est - z * se),
               CI_Upper = exp(est + z * se),
               Method = "Wald test of the log volume ratio, endpoint model",
               stringsAsFactors = FALSE)
  }))
}

#' Endpoint-model diagnostics: residual Q-Q by arm and the fitted curves
#' @noRd
#' @keywords internal
synergy_model_diagnostics <- function(em, arms, t) {
  if (is.null(em) || !requireNamespace("ggplot2", quietly = TRUE)) {
    return(list(qq = NULL, box = NULL, fit = NULL))
  }
  fr <- em$fit@frame
  res <- stats::residuals(em$fit)
  dd <- data.frame(Treatment = as.character(fr$Treatment), resid = res,
                   stringsAsFactors = FALSE)
  dd <- dd[dd$Treatment %in% arms, , drop = FALSE]
  qq <- tryCatch({
    long <- do.call(rbind, lapply(split(dd, dd$Treatment), function(g) {
      if (nrow(g) < 3L) return(NULL)
      q <- stats::qqnorm(g$resid / stats::sd(res), plot.it = FALSE)
      data.frame(group = g$Treatment[1], theoretical = q$x, sample = q$y)
    }))
    ggplot2::ggplot(long, ggplot2::aes(.data[["theoretical"]], .data[["sample"]])) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
      ggplot2::facet_wrap(~ group) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Endpoint model residuals by arm (standardised)",
                    subtitle = "Points near the line: the log-normal model fits",
                    x = "Theoretical quantiles", y = "Standardised residuals")
  }, error = function(e) NULL)
  fit <- tryCatch({
    grid <- seq(em$day_range[1], t, length.out = 60L)
    curves <- do.call(rbind, lapply(arms, function(a) data.frame(
      Treatment = a, Day = grid,
      logv = as.numeric(me_endpoint_X(em, a, grid) %*% em$beta))))
    obs <- data.frame(Treatment = as.character(fr$Treatment),
                      Day = fr$.day_c + em$day_mean, logv = fr$.logv,
                      stringsAsFactors = FALSE)
    obs <- obs[obs$Treatment %in% arms, , drop = FALSE]
    ggplot2::ggplot(obs, ggplot2::aes(.data[["Day"]], .data[["logv"]])) +
      ggplot2::geom_point(alpha = 0.35) +
      ggplot2::geom_line(data = curves, colour = "steelblue", linewidth = 1) +
      ggplot2::geom_vline(xintercept = t, linetype = "dashed") +
      ggplot2::facet_wrap(~ Treatment) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Endpoint model: each arm's fitted curve over its data",
                    subtitle = paste0("Log volume. The dashed line is the evaluation day (",
                                      t, "); the curve should follow the points up to it."),
                    x = "Day", y = "log(volume)")
  }, error = function(e) NULL)
  list(qq = qq, box = NULL, fit = fit)
}

#' Per-animal diagnostics for the per-animal estimands
#' @noRd
#' @keywords internal
synergy_permouse_diagnostics <- function(pm, t) {
  if (is.null(pm) || !requireNamespace("ggplot2", quietly = TRUE)) {
    return(list(qq = NULL, box = NULL, fit = NULL))
  }
  qq <- tryCatch({
    long <- do.call(rbind, lapply(split(pm, pm$Treatment), function(g) {
      v <- g$Volume[is.finite(g$Volume)]
      if (length(v) < 3L) return(NULL)
      q <- stats::qqnorm(v, plot.it = FALSE)
      data.frame(group = g$Treatment[1], theoretical = q$x, sample = q$y)
    }))
    ggplot2::ggplot(long, ggplot2::aes(.data[["theoretical"]], .data[["sample"]])) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "red",
                           linetype = "dashed", formula = y ~ x) +
      ggplot2::facet_wrap(~ group, scales = "free") +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Per-group Q-Q of per-animal volumes",
                    subtitle = paste("Evaluation day", t,
                                     "- heavy tails or curvature weaken the t-tests"),
                    x = "Theoretical Quantiles", y = "Sample Quantiles")
  }, error = function(e) NULL)
  box <- tryCatch({
    ggplot2::ggplot(pm, ggplot2::aes(.data[["Treatment"]], .data[["Volume"]])) +
      ggplot2::geom_boxplot(outlier.colour = "red", outlier.alpha = 0.8) +
      ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.4) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = "Per-animal volumes used at the evaluation day",
                    subtitle = paste("Day", t, "- outliers in red"),
                    x = "Treatment", y = "Tumor volume")
  }, error = function(e) NULL)
  list(qq = qq, box = box, fit = NULL)
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
      subtitle = paste("Synergy Assessment:", synergy_results$overall_assessment),
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