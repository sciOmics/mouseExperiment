# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Bayesian Linear Mixed-Effects Model for Body Weight
#'
#' Fits a Bayesian linear mixed-effects model to longitudinal body weight data
#' using \pkg{brms} (Bürkner, 2017). Provides full posterior distributions,
#' 95 % highest-posterior-density (HPD) credible intervals, a weight-loss
#' summary, and posterior weight trajectory plots. Complements the frequentist
#' \code{\link{analyze_body_weight}} function.
#'
#' @param df Data frame containing longitudinal body weight data.
#' @param weight_column Column name for body weight (grams). Default
#'   \code{"Weight"}.
#' @param time_column Column name for study day (numeric). Default \code{"Day"}.
#' @param treatment_column Column name for treatment group. Default
#'   \code{"Treatment"}.
#' @param id_column Column name for individual animal ID. Default \code{"ID"}.
#' @param cage_column Column name for cage ID, or \code{NULL} (default).
#'   When supplied and \code{include_cage_effect = TRUE}, a cage-level random
#'   intercept is added: \code{(1|cage/ID)} for intercept-only random effects,
#'   \code{(Day|ID) + (1|cage)} for random slopes.
#' @param volume_column Column name for tumor volume (mm³), used only when
#'   \code{adjust_tumor_weight = TRUE}. \code{NULL} to skip adjustment.
#' @param adjust_tumor_weight Logical. If \code{TRUE} and
#'   \code{volume_column} is supplied, subtract estimated tumor weight from
#'   gross body weight before modelling. Default \code{FALSE}.
#' @param tumor_density Numeric g/cm³ used to convert volume to mass (default
#'   \code{1.0}).
#' @param transform Weight transformation applied before modelling:
#'   \code{"none"} (default, appropriate for approximately normal weight
#'   values), \code{"log"}, or \code{"sqrt"}.
#' @param random_effects_specification Random effects structure:
#'   \code{"intercept_only"} (default, \code{(1 | ID)}) or \code{"slope"}
#'   (\code{(Day | ID)}).
#' @param reference_group Treatment group used as reference. Auto-detected
#'   (first alphabetically) if \code{NULL}.
#' @param prior_strength Prior preset:
#'   \describe{
#'     \item{\code{"skeptical"}}{(default) \eqn{b \sim N(0, 0.25)};
#'       \eqn{\text{sd}, \sigma \sim \text{Exponential}(2)}.}
#'     \item{\code{"weakly_informative"}}{\eqn{b \sim N(0, 1)};
#'       \eqn{\text{sd}, \sigma \sim \text{Exponential}(1)}.}
#'     \item{\code{"informative"}}{\eqn{b \sim N(0, 0.5)};
#'       \eqn{\text{sd}, \sigma \sim \text{Exponential}(2)}.}
#'     \item{\code{"diffuse"}}{\eqn{b \sim N(0, 2.5)};
#'       \eqn{\text{sd}, \sigma \sim \text{Exponential}(0.5)}.}
#'     \item{\code{"manual"}}{Use the \code{prior_*} arguments directly.}
#'   }
#' @param prior_b,prior_intercept,prior_sd,prior_sigma brms prior strings
#'   (used only when \code{prior_strength = "manual"}).
#' @param n_chains Number of MCMC chains. Default \code{4}.
#' @param n_warmup Warm-up (burn-in) iterations per chain. Default \code{1000}.
#' @param n_iter Post-warmup draws per chain. Default \code{500}.
#' @param seed Integer random seed. Default \code{42}.
#' @param include_cage_effect Logical. Include cage random intercept when
#'   \code{cage_column} is supplied? Default \code{TRUE}.
#' @param return_model Logical. Return the \code{brmsfit} object? Default
#'   \code{TRUE}.
#' @param plots Logical. Pre-compute and return \pkg{ggplot2} plot objects?
#'   Default \code{TRUE}.
#' @param verbose Logical. Show Stan messages? Default \code{FALSE}.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{model}}{\code{brmsfit} object, or \code{NULL} when
#'     \code{return_model = FALSE}.}
#'   \item{\code{model_type_used}}{Character \code{"bayes_bw"}.}
#'   \item{\code{transform_used}}{Transformation applied.}
#'   \item{\code{summary}}{Named list of analysis metadata.}
#'   \item{\code{posterior_summary}}{Data frame of fixed-effect posterior
#'     summaries (median, 95 % CrI, Rhat, Bulk_ESS, Tail_ESS).}
#'   \item{\code{treatment_effects}}{Posterior EMMs per group with columns
#'     \code{Group}, \code{Adjusted_Mean}, \code{SE}, \code{DF},
#'     \code{Lower_CrI}, \code{Upper_CrI}, \code{Note}.}
#'   \item{\code{pairwise_comparisons}}{Posterior contrasts vs reference.}
#'   \item{\code{mcmc_diagnostics}}{Data frame with Rhat and convergence flag
#'     per parameter.}
#'   \item{\code{weight_loss_summary}}{Data frame: per-group posterior
#'     percentage weight change from first to last study day (median, 95 %
#'     CrI). Weight-loss percentages are computed from the earliest to the
#'     latest study day present in the model data. This is a population-level
#'     prediction; individual dropout patterns do not affect the endpoint
#'     days used.}
#'   \item{\code{pp_check_plot}}{Posterior predictive density overlay.}
#'   \item{\code{posterior_dist_plot}}{Treatment-parameter posterior areas.}
#'   \item{\code{prior_posterior_plot}}{Prior vs posterior density overlay.}
#'   \item{\code{credible_intervals_plot}}{Forest plot of group EMMs.}
#'   \item{\code{mcmc_trace_plot}}{MCMC trace for treatment parameters.}
#'   \item{\code{residuals_plot}}{Residuals vs day, faceted by group.}
#'   \item{\code{weight_trajectory_plot}}{Posterior weight trajectory per
#'     group: median + 95 % CrI ribbon over time with observed points.}
#'   \item{\code{weight_data}}{The prepared analysis data frame.}
#' }
#'
#' @references
#' Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models
#' using Stan. \emph{Journal of Statistical Software}, 80(1), 1–28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @importFrom stats as.formula relevel quantile median
#' @importFrom ggplot2 ggplot aes geom_density geom_line geom_point
#'   geom_pointrange geom_ribbon geom_smooth geom_vline geom_hline
#'   scale_colour_manual scale_fill_manual facet_wrap labs theme theme_classic
#' @export
bayesian_body_weight <- function(
  df,
  weight_column                = "Weight",
  time_column                  = "Day",
  treatment_column             = "Treatment",
  id_column                    = "ID",
  cage_column                  = NULL,
  volume_column                = NULL,
  adjust_tumor_weight          = FALSE,
  tumor_density                = 1.0,
  transform                    = c("none", "log", "sqrt"),
  random_effects_specification = c("intercept_only", "slope"),
  reference_group              = NULL,
  prior_strength               = c("skeptical", "weakly_informative",
                                   "informative", "diffuse", "manual"),
  prior_b                      = NULL,
  prior_intercept              = NULL,
  prior_sd                     = NULL,
  prior_sigma                  = NULL,
  n_chains                     = 4L,
  n_warmup                     = 1000L,
  n_iter                       = 500L,
  seed                         = 42L,
  include_cage_effect          = TRUE,
  return_model                 = TRUE,
  plots                        = TRUE,
  verbose                      = FALSE
) {

  # ── Dependency check ───────────────────────────────────────────────────────
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for Bayesian analysis.\n",
      "Install it with: install.packages('brms')"
    )
  }

  transform                    <- match.arg(transform)
  random_effects_specification <- match.arg(random_effects_specification)
  prior_strength               <- match.arg(prior_strength)

  # ── Column validation ──────────────────────────────────────────────────────
  required_cols <- c(weight_column, time_column, treatment_column, id_column)
  missing_cols  <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ── Cage placeholder ───────────────────────────────────────────────────────
  cage_setup   <- setup_cage_column(df, cage_column)
  df           <- cage_setup$df
  cage_column  <- cage_setup$cage_column
  no_cage_mode <- cage_setup$no_cage_mode

  # ── Reference group ────────────────────────────────────────────────────────
  treatment_groups <- unique(as.character(df[[treatment_column]]))
  if (is.null(reference_group)) {
    reference_group <- sort(treatment_groups)[1]
  } else if (!reference_group %in% treatment_groups) {
    stop(
      "Reference group '", reference_group,
      "' not found in treatment column."
    )
  }
  df[[treatment_column]] <- stats::relevel(
    factor(df[[treatment_column]]), ref = reference_group
  )

  # ── Data preparation ───────────────────────────────────────────────────────
  analysis_df <- df

  # Optional tumor weight subtraction
  has_volume <- !is.null(volume_column) && volume_column %in% colnames(df)
  if (isTRUE(adjust_tumor_weight) && has_volume) {
    analysis_df$Net_Weight <- (
      as.numeric(df[[weight_column]]) -
        as.numeric(df[[volume_column]]) / 1000 * tumor_density
    )
  } else {
    analysis_df$Net_Weight <- as.numeric(df[[weight_column]])
  }

  # Remove missing response or time
  analysis_df <- analysis_df[
    !is.na(analysis_df$Net_Weight) & !is.na(analysis_df[[time_column]]),
  ]
  if (nrow(analysis_df) == 0L) {
    stop("No valid weight observations after removing missing values.")
  }

  # Apply transform
  if (transform == "log") {
    wt            <- analysis_df$Net_Weight
    positive_vals <- wt[is.finite(wt) & wt > 0]
    if (length(positive_vals) == 0L) {
      stop(
        "No positive values in '", weight_column,
        "' after filtering. Cannot apply log transform."
      )
    }
    wt[wt <= 0] <- min(positive_vals, na.rm = TRUE) / 2
    analysis_df$Net_Weight <- log(wt)
  } else if (transform == "sqrt") {
    analysis_df$Net_Weight <- sqrt(analysis_df$Net_Weight)
  }

  # ── Formula ────────────────────────────────────────────────────────────────
  use_cage_re <- isTRUE(include_cage_effect) && !no_cage_mode

  re_term <- if (use_cage_re) {
    if (random_effects_specification == "slope") {
      paste0(
        "(", time_column, " | ", id_column, ")",
        " + (1 | ", cage_column, ")"
      )
    } else {
      paste0("(1 | ", cage_column, "/", id_column, ")")
    }
  } else {
    if (random_effects_specification == "slope") {
      paste0("(", time_column, " | ", id_column, ")")
    } else {
      paste0("(1 | ", id_column, ")")
    }
  }

  brms_formula <- stats::as.formula(
    paste(
      "Net_Weight ~",
      treatment_column, "*", time_column, "+", re_term
    )
  )

  # ── Priors ─────────────────────────────────────────────────────────────────
  if (prior_strength == "manual") {
    if (any(is.null(c(prior_b, prior_intercept, prior_sd, prior_sigma)))) {
      stop(
        "When prior_strength = 'manual', ",
        "all four prior_* arguments must be supplied."
      )
    }
    selected_priors <- c(
      brms::prior_string(prior_b,         class = "b"),
      brms::prior_string(prior_intercept, class = "Intercept"),
      brms::prior_string(prior_sd,        class = "sd"),
      brms::prior_string(prior_sigma,     class = "sigma")
    )
  } else {
    pp       <- bayes_prior_params(prior_strength)
    b_sd     <- pp$b_sd
    exp_rate <- pp$exp_rate
    selected_priors <- c(
      brms::prior_string(
        paste0("normal(0, ", b_sd, ")"),        class = "b"
      ),
      brms::prior_string(
        paste0("normal(0, ", b_sd * 2.5, ")"),  class = "Intercept"
      ),
      brms::prior_string(
        paste0("exponential(", exp_rate, ")"),  class = "sd"
      ),
      brms::prior_string(
        paste0("exponential(", exp_rate, ")"),  class = "sigma"
      )
    )
  }

  # ── Fit ────────────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    message(
      "Fitting Bayesian LMM via brms (",
      n_chains, " chains × ", n_iter, " post-warmup draws, ",
      n_warmup, " warmup)..."
    )
  }

  model <- brms::brm(
    formula      = brms_formula,
    data         = analysis_df,
    prior        = selected_priors,
    sample_prior = "yes",
    chains       = as.integer(n_chains),
    cores        = as.integer(n_chains),
    iter         = as.integer(n_warmup + n_iter),
    warmup       = as.integer(n_warmup),
    seed         = as.integer(seed),
    silent       = if (isTRUE(verbose)) 0L else 2L,
    refresh      = if (isTRUE(verbose)) 100L else 0L
  )

  # ── Posterior summary ──────────────────────────────────────────────────────
  brms_smry     <- summary(model)
  fixed_df      <- as.data.frame(brms_smry$fixed)
  fixed_df      <- cbind(
    Parameter = rownames(fixed_df), fixed_df, stringsAsFactors = FALSE
  )
  rownames(fixed_df) <- NULL
  names(fixed_df)[names(fixed_df) == "l-95% CI"] <- "Lower_95_CrI"
  names(fixed_df)[names(fixed_df) == "u-95% CI"] <- "Upper_95_CrI"
  posterior_summary <- fixed_df

  # ── MCMC diagnostics ───────────────────────────────────────────────────────
  mcmc_diagnostics <- make_mcmc_diagnostics(posterior_summary,
                                            total_draws = n_chains * n_iter)
  nuts_diagnostics <- make_nuts_diagnostics(model)
  loo_diagnostics  <- bayes_loo(model)
  bayes_r2         <- bayes_r2_summary(model)

  # ── Treatment effects (emmeans) ────────────────────────────────────────────
  treatment_effects    <- NULL
  pairwise_comparisons <- NULL

  if (requireNamespace("emmeans", quietly = TRUE)) {
    mean_day <- mean(analysis_df[[time_column]])

    emm <- tryCatch(
      emmeans::emmeans(
        model, specs = treatment_column,
        at = stats::setNames(list(mean_day), time_column)
      ),
      error = function(e) {
        warning(
          "emmeans on brmsfit failed: ", conditionMessage(e),
          "\nTreatment effects not computed."
        )
        NULL
      }
    )

    if (!is.null(emm)) {
      emm_df    <- as.data.frame(summary(emm, point.est = stats::median))
      lower_col <- if ("lower.HPD" %in% names(emm_df)) {
        "lower.HPD"
      } else {
        "lower.CL"
      }
      upper_col <- if ("upper.HPD" %in% names(emm_df)) {
        "upper.HPD"
      } else {
        "upper.CL"
      }

      treatment_effects <- data.frame(
        Group         = as.character(emm_df[[treatment_column]]),
        Adjusted_Mean = round(emm_df$emmean,         3),
        SE            = round(emm_df$SE,              3),
        DF            = NA_real_,
        Lower_CrI     = round(emm_df[[lower_col]],   3),
        Upper_CrI     = round(emm_df[[upper_col]],   3),
        Note          = ifelse(
          as.character(emm_df[[treatment_column]]) == reference_group,
          "Reference group", ""
        ),
        stringsAsFactors = FALSE
      )

      ref_idx <- which(treatment_effects$Group == reference_group)
      if (length(ref_idx) > 0 && ref_idx[1] > 1L) {
        treatment_effects <- rbind(
          treatment_effects[ref_idx[1], ],
          treatment_effects[-ref_idx[1], ]
        )
        rownames(treatment_effects) <- NULL
      }

      pc <- tryCatch(
        emmeans::contrast(emm, method = "trt.vs.ctrl"),
        error = function(e) NULL
      )
      if (!is.null(pc)) {
        pc_df    <- as.data.frame(summary(pc, point.est = stats::median))
        lower_pc <- if ("lower.HPD" %in% names(pc_df)) {
          "lower.HPD"
        } else {
          "lower.CL"
        }
        upper_pc <- if ("upper.HPD" %in% names(pc_df)) {
          "upper.HPD"
        } else {
          "upper.CL"
        }
        se_col <- if ("SE" %in% names(pc_df)) round(pc_df$SE, 4) else NA_real_
        pairwise_comparisons <- data.frame(
          contrast  = as.character(pc_df$contrast),
          estimate  = round(pc_df$estimate,    4),
          SE        = se_col,
          Lower_CrI = round(pc_df[[lower_pc]], 4),
          Upper_CrI = round(pc_df[[upper_pc]], 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ── Weight loss summary (population-level posterior predictions) ───────────
  weight_loss_summary <- NULL
  study_days  <- sort(unique(analysis_df[[time_column]]))
  tx_groups   <- levels(analysis_df[[treatment_column]])
  first_day   <- study_days[1L]
  last_day    <- study_days[length(study_days)]

  make_pred_nd <- function(day) {
    nd <- data.frame(
      tx = tx_groups, dy = rep(day, length(tx_groups)),
      stringsAsFactors = FALSE
    )
    names(nd) <- c(treatment_column, time_column)
    nd
  }

  # Combine first-day and last-day prediction into a single posterior_epred
  # call (B5.1): halves one brms round-trip.
  n_groups     <- length(tx_groups)
  nd_endpoints <- rbind(make_pred_nd(first_day), make_pred_nd(last_day))
  epred_endpts <- tryCatch(
    brms::posterior_epred(
      model, newdata = nd_endpoints,
      re_formula = NA, allow_new_levels = TRUE
    ),
    error = function(e) NULL
  )
  epred_base <- if (!is.null(epred_endpts)) {
    epred_endpts[, seq_len(n_groups), drop = FALSE]
  } else {
    NULL
  }
  epred_end  <- if (!is.null(epred_endpts)) {
    epred_endpts[, seq_len(n_groups) + n_groups, drop = FALSE]
  } else {
    NULL
  }

  if (!is.null(epred_base) && !is.null(epred_end)) {
    eb  <- bayes_backtransform(epred_base, transform)
    ee  <- bayes_backtransform(epred_end,  transform)
    pct <- (ee - eb) / eb * 100

    weight_loss_summary <- do.call(rbind, lapply(
      seq_along(tx_groups),
      function(i) {
        q <- stats::quantile(pct[, i], c(0.025, 0.5, 0.975))
        data.frame(
          Group        = tx_groups[i],
          Pct_Change   = round(q["50%"],   1),
          Lower_CrI    = round(q["2.5%"],  1),
          Upper_CrI    = round(q["97.5%"], 1),
          Day_Baseline = first_day,
          Day_End      = last_day,
          stringsAsFactors = FALSE
        )
      }
    ))
    rownames(weight_loss_summary) <- NULL
  }

  # ── Plots ──────────────────────────────────────────────────────────────────
  pp_check_plot           <- NULL
  posterior_dist_plot     <- NULL
  prior_posterior_plot    <- NULL
  credible_intervals_plot <- NULL
  mcmc_trace_plot         <- NULL
  residuals_plot          <- NULL
  weight_trajectory_plot  <- NULL

  if (isTRUE(plots)) {

    pp_check_plot <- tryCatch(
      brms::pp_check(model, type = "dens_overlay", ndraws = 50),
      error = function(e) NULL
    )

    if (requireNamespace("bayesplot", quietly = TRUE)) {
      draws_arr <- tryCatch(brms::as.array(model), error = function(e) NULL)
      if (!is.null(draws_arr)) {
        all_pars <- dimnames(draws_arr)$variable
        safe_tx  <- gsub(
          "([.^$*+?()\\[\\]{}|])", "\\\\\\1", treatment_column
        )
        tx_pars  <- grep(
          paste0("^b_", safe_tx), all_pars, value = TRUE
        )
        if (length(tx_pars) > 0) {
          posterior_dist_plot <- tryCatch(
            bayesplot::mcmc_areas(draws_arr, pars = tx_pars, prob = 0.95),
            error = function(e) NULL
          )
          mcmc_trace_plot <- tryCatch(
            bayesplot::mcmc_trace(draws_arr, pars = tx_pars),
            error = function(e) NULL
          )
        }
      }
    }

    prior_posterior_plot <- tryCatch(
      bayes_prior_posterior_plot(model, treatment_column),
      error = function(e) NULL
    )

    if (!is.null(treatment_effects) && nrow(treatment_effects) > 0) {
      te        <- treatment_effects
      te$Group  <- factor(te$Group, levels = rev(te$Group))
      ref_val   <- te$Adjusted_Mean[
        as.character(te$Group) == reference_group
      ]
      is_ref    <- as.character(te$Group) == reference_group

      credible_intervals_plot <- ggplot2::ggplot(
        te,
        ggplot2::aes(
          x      = .data[["Adjusted_Mean"]],
          y      = .data[["Group"]],
          xmin   = .data[["Lower_CrI"]],
          xmax   = .data[["Upper_CrI"]],
          colour = is_ref
        )
      ) +
        ggplot2::geom_pointrange(size = 0.8, linewidth = 0.8) +
        ggplot2::geom_vline(
          xintercept = if (length(ref_val) > 0) ref_val[1] else NA_real_,
          linetype = "dashed", colour = "grey50", linewidth = 0.5
        ) +
        ggplot2::scale_colour_manual(
          values = c("TRUE" = "grey40", "FALSE" = "steelblue"),
          guide  = "none"
        ) +
        ggplot2::labs(
          title    = "Treatment Effects — 95 % Credible Intervals",
          subtitle = paste0(
            "Posterior medians at mean day (",
            round(mean(analysis_df[[time_column]]), 1),
            "); ", transform, " scale"
          ),
          x = paste0("Estimated marginal mean (", transform, " weight)"),
          y = NULL
        ) +
        ggplot2::theme_classic(base_size = 14)
    }

    residuals_plot <- tryCatch({
      resids   <- residuals(model, type = "ordinary")[, "Estimate"]
      resid_df <- data.frame(
        study_day = analysis_df[[time_column]],
        Residual  = resids,
        Treatment = analysis_df[[treatment_column]]
      )
      ggplot2::ggplot(
        resid_df,
        ggplot2::aes(x = .data[["study_day"]], y = .data[["Residual"]])
      ) +
        ggplot2::geom_point(alpha = 0.35, size = 1.5) +
        ggplot2::geom_smooth(
          method  = "loess", se = TRUE, formula = y ~ x,
          colour  = "steelblue", linewidth = 0.8,
          fill    = "steelblue", alpha = 0.15
        ) +
        ggplot2::geom_hline(
          yintercept = 0, linetype = "dashed",
          colour = "grey50", linewidth = 0.5
        ) +
        ggplot2::facet_wrap(~ Treatment) +
        ggplot2::labs(
          title    = "Residuals vs. Study Day",
          subtitle = paste0("Response: ", transform, " weight"),
          x        = time_column,
          y        = paste0("Residual (", transform, " scale)")
        ) +
        ggplot2::theme_classic(base_size = 14)
    }, error = function(e) NULL)

    # Weight trajectory: posterior median + 95 % CrI ribbon + observed points
    weight_trajectory_plot <- tryCatch({
      traj_nd  <- do.call(rbind, lapply(study_days, make_pred_nd))
      ep_traj  <- brms::posterior_epred(
        model, newdata = traj_nd,
        re_formula = NA, allow_new_levels = TRUE
      )
      ep_traj  <- bayes_backtransform(ep_traj, transform)
      traj_nd$Median <- apply(ep_traj, 2L, stats::median)
      traj_nd$Lower  <- apply(
        ep_traj, 2L, function(x) stats::quantile(x, 0.025)
      )
      traj_nd$Upper  <- apply(
        ep_traj, 2L, function(x) stats::quantile(x, 0.975)
      )
      traj_nd$Trt_label <- as.character(traj_nd[[treatment_column]])

      obs_y  <- bayes_backtransform(analysis_df$Net_Weight, transform)
      obs_df <- data.frame(
        Day_obs   = analysis_df[[time_column]],
        Trt_label = as.character(analysis_df[[treatment_column]]),
        Wt_obs    = obs_y,
        stringsAsFactors = FALSE
      )

      net_adj  <- isTRUE(adjust_tumor_weight) && has_volume
      wt_label <- if (net_adj) "Net weight" else "Weight"
      y_label  <- if (transform == "none") {
        paste0(wt_label, " (g)")
      } else {
        paste0("Weight (g) [", transform, " scale]")
      }

      ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = traj_nd,
          ggplot2::aes(
            x    = .data[[time_column]],
            ymin = .data[["Lower"]],
            ymax = .data[["Upper"]],
            fill = .data[["Trt_label"]]
          ),
          alpha = 0.2
        ) +
        ggplot2::geom_line(
          data = traj_nd,
          ggplot2::aes(
            x      = .data[[time_column]],
            y      = .data[["Median"]],
            colour = .data[["Trt_label"]]
          ),
          linewidth = 1
        ) +
        ggplot2::geom_point(
          data = obs_df,
          ggplot2::aes(
            x      = .data[["Day_obs"]],
            y      = .data[["Wt_obs"]],
            colour = .data[["Trt_label"]]
          ),
          alpha = 0.35, size = 1.5
        ) +
        ggplot2::labs(
          title    = "Body Weight Trajectory",
          subtitle = "Posterior median ± 95 % CrI; points = observed",
          x        = time_column,
          y        = y_label,
          colour   = treatment_column,
          fill     = treatment_column
        ) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::theme(legend.position = "top")
    }, error = function(e) NULL)
  }

  # ── Analysis summary metadata ──────────────────────────────────────────────
  analysis_summary <- list(
    analysis_type = "Bayesian Linear Mixed-Effects Model — Body Weight (brms)",
    data_description = list(
      subjects         = length(unique(analysis_df[[id_column]])),
      treatment_groups = length(unique(analysis_df[[treatment_column]])),
      time_points      = length(unique(analysis_df[[time_column]])),
      reference_group  = reference_group,
      adjust_tumor     = isTRUE(adjust_tumor_weight) && has_volume
    ),
    model_specification = list(
      fixed_effects        = paste(
        "Net_Weight ~", treatment_column, "*", time_column
      ),
      random_effects       = re_term,
      cage_effect_modelled = use_cage_re,
      prior_strength       = prior_strength,
      transform            = transform
    ),
    methods = list(
      engine = paste0(
        "brms (", n_chains, " chains × ",
        n_iter, " draws + ", n_warmup, " warmup, seed = ", seed, ")"
      ),
      prior_b = if (prior_strength == "manual") {
        prior_b
      } else {
        paste0("normal(0, ", b_sd, ")")
      },
      prior_intercept = if (prior_strength == "manual") {
        prior_intercept
      } else {
        paste0("normal(0, ", round(b_sd * 2.5, 2), ")")
      },
      prior_sd = if (prior_strength == "manual") {
        prior_sd
      } else {
        paste0("exponential(", exp_rate, ")")
      },
      prior_sigma = if (prior_strength == "manual") {
        prior_sigma
      } else {
        paste0("exponential(", exp_rate, ")")
      }
    )
  )

  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    model                   = if (isTRUE(return_model)) model else NULL,
    model_type_used         = "bayes_bw",
    transform_used          = transform,
    summary                 = analysis_summary,
    posterior_summary       = posterior_summary,
    treatment_effects       = treatment_effects,
    pairwise_comparisons    = pairwise_comparisons,
    mcmc_diagnostics        = mcmc_diagnostics,
    nuts_diagnostics        = nuts_diagnostics,
    loo_diagnostics         = loo_diagnostics,
    bayes_R2                = bayes_r2,
    weight_loss_summary     = weight_loss_summary,
    pp_check_plot           = pp_check_plot,
    posterior_dist_plot     = posterior_dist_plot,
    prior_posterior_plot    = prior_posterior_plot,
    credible_intervals_plot = credible_intervals_plot,
    mcmc_trace_plot         = mcmc_trace_plot,
    residuals_plot          = residuals_plot,
    weight_trajectory_plot  = weight_trajectory_plot,
    weight_data             = analysis_df
  )
}
