# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Bayesian Linear Mixed-Effects Model for Tumor Growth
#'
#' Fits a Bayesian linear mixed-effects model to longitudinal tumor volume data
#' using \pkg{brms} (Bürkner, 2017). Provides full posterior distributions,
#' 95 % highest-posterior-density (HPD) credible intervals, and posterior
#' predictive checks as complements to the frequentist
#' \code{\link{tumor_growth_statistics}} function.
#'
#' @param df Data frame containing tumor growth data.
#' @param time_column Column name for study day (numeric). Default \code{"Day"}.
#' @param volume_column Column name for tumor volume. Default \code{"Volume"}.
#' @param treatment_column Column name for treatment group. Default \code{"Treatment"}.
#' @param id_column Column name for individual animal ID. Default \code{"ID"}.
#' @param cage_column Column name for cage ID, or \code{NULL} (default).
#'   When supplied, cage effects are not modelled but the column is used for
#'   composite mouse-key construction (consistent with the \code{lme4} path).
#' @param dose_column Reserved for future use; currently ignored.
#' @param transform Volume transformation applied before modelling:
#'   \code{"log"} (default, recommended for exponential growth),
#'   \code{"sqrt"}, or \code{"none"}.
#' @param random_effects_specification Random effects structure:
#'   \code{"intercept_only"} (default, \code{(1 | ID)}) or
#'   \code{"slope"} (\code{(Day | ID)}, adds per-animal random slopes).
#' @param reference_group Treatment group used as the reference (Intercept)
#'   level. Auto-detected from common control names if \code{NULL}.
#' @param prior_strength Prior preset applied to all fixed effects:
#'   \describe{
#'     \item{\code{"skeptical"}}{(default) \eqn{b \sim N(0, 0.25)};
#'       \eqn{\text{sd}, \sigma \sim \text{Exponential}(2)}.
#'       Expresses prior belief that treatment effects are small; requires
#'       stronger data to support large estimated differences.}
#'     \item{\code{"weakly_informative"}}{\eqn{b \sim N(0, 1)};
#'       \eqn{\text{sd}, \sigma \sim \text{Exponential}(1)}.}
#'     \item{\code{"informative"}}{\eqn{b \sim N(0, 0.5)};
#'       \eqn{\text{sd}, \sigma \sim \text{Exponential}(2)}.}
#'     \item{\code{"diffuse"}}{\eqn{b \sim N(0, 2.5)};
#'       \eqn{\text{sd}, \sigma \sim \text{Exponential}(0.5)}.}
#'     \item{\code{"manual"}}{Use the \code{prior_b}, \code{prior_intercept},
#'       \code{prior_sd}, and \code{prior_sigma} arguments to specify brms
#'       prior strings directly.}
#'   }
#'   All presets use an intercept prior with SD = \eqn{2.5 \times \sigma_b}.
#' @param prior_b brms prior string for fixed-effect coefficients (class
#'   \code{"b"}), e.g. \code{"normal(0, 0.25)"} or \code{"student_t(3, 0, 0.5)"}.
#'   Only used when \code{prior_strength = "manual"}.
#' @param prior_intercept brms prior string for the model intercept, e.g.
#'   \code{"normal(0, 2.5)"}. Only used when \code{prior_strength = "manual"}.
#' @param prior_sd brms prior string for random-effect standard deviations
#'   (class \code{"sd"}), e.g. \code{"exponential(2)"}.
#'   Only used when \code{prior_strength = "manual"}.
#' @param prior_sigma brms prior string for the residual standard deviation
#'   (class \code{"sigma"}), e.g. \code{"exponential(2)"}.
#'   Only used when \code{prior_strength = "manual"}.
#' @param n_chains Number of MCMC chains. Default \code{4}.
#' @param n_iter Total iterations per chain (including warmup). Default \code{2000}.
#' @param seed Integer random seed for reproducibility. Default \code{42}.
#' @param return_model Logical. Return the \code{brmsfit} object? Default
#'   \code{TRUE}. Set \code{FALSE} to reduce memory when only summaries are
#'   needed.
#' @param plots Logical. Pre-compute and return \pkg{ggplot2} / \pkg{bayesplot}
#'   objects? Default \code{TRUE}. Set \code{FALSE} when plots are generated
#'   reactively from the stored model object (e.g. in Shiny).
#' @param verbose Logical. Show Stan compilation and sampling messages? Default
#'   \code{FALSE}.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{model}}{\code{brmsfit} object, or \code{NULL} when
#'     \code{return_model = FALSE}.}
#'   \item{\code{model_type_used}}{Character \code{"bayes"}.}
#'   \item{\code{transform_used}}{Transformation applied.}
#'   \item{\code{summary}}{Named list of analysis metadata (analysis_type,
#'     data_description, model_specification, methods).}
#'   \item{\code{posterior_summary}}{Data frame of fixed-effect posterior
#'     medians, 2.5 %–97.5 % CrI, Rhat, Bulk_ESS, Tail_ESS for all
#'     parameters.}
#'   \item{\code{treatment_effects}}{Data frame with columns
#'     \code{Group}, \code{Adjusted_Mean}, \code{SE}, \code{DF},
#'     \code{Lower_CL}, \code{Upper_CL}, \code{Note} — compatible with the
#'     \code{lme4} path. Values are posterior medians and 95 % HPD CrI of
#'     group-level EMMs at the mean study day.}
#'   \item{\code{pairwise_comparisons}}{Data frame of posterior contrasts vs
#'     the reference group (columns: \code{contrast}, \code{estimate},
#'     \code{Lower_CL}, \code{Upper_CL}).}
#'   \item{\code{mcmc_diagnostics}}{Data frame with \code{Parameter},
#'     \code{Rhat}, \code{Bulk_ESS}, \code{Tail_ESS}, \code{Converged}
#'     per fixed-effect parameter. Rhat > 1.01 is flagged as not converged.}
#'   \item{\code{pp_check_plot}}{Posterior predictive density overlay ggplot,
#'     or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{posterior_dist_plot}}{Posterior density areas for treatment
#'     parameters (bayesplot), or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{credible_intervals_plot}}{Forest plot of group EMMs with
#'     95 % CrI, or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{mcmc_trace_plot}}{MCMC trace plot for treatment parameters,
#'     or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{growth_rates}}{Per-animal exponential growth rates (same
#'     structure as the \code{lme4} path; estimated by log-linear OLS, not
#'     MCMC).}
#'   \item{\code{data_summary}}{Descriptive statistics by treatment group and
#'     day.}
#' }
#'
#' @section Model:
#' The model fitted (intercept-only random effects, log transform, linear time):
#' \deqn{
#'   \log V_{ij} = \mu + \beta_T T_i + \beta_t t_{ij} + \beta_{Tt} T_i t_{ij}
#'                 + u_{0i} + \varepsilon_{ij}
#' }
#' where \eqn{u_{0i} \sim N(0, \sigma_u^2)} and
#' \eqn{\varepsilon_{ij} \sim N(0, \sigma^2)}.
#'
#' @section Convergence:
#' Rhat values > 1.01 in \code{mcmc_diagnostics} indicate potential
#' non-convergence. Increasing \code{n_iter} or inspecting trace plots usually
#' resolves this.
#'
#' @section Assumptions and Limitations:
#' \itemize{
#'   \item Computation typically takes 3–12 minutes on modern hardware. The
#'     Stan model is compiled on the first call with a given formula and
#'     cached; subsequent calls with the same formula start faster.
#'   \item Credible intervals have a direct probability interpretation
#'     ("there is 95 % posterior probability the parameter lies in this
#'     interval"), unlike frequentist confidence intervals.
#'   \item Treatment effects and pairwise comparisons use estimated marginal
#'     means from \pkg{emmeans}, which draws from the posterior to marginalise.
#'     Results may differ slightly across calls due to MCMC variability.
#' }
#'
#' @references
#' Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models
#' using Stan. \emph{Journal of Statistical Software}, 80(1), 1–28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @importFrom stats as.formula relevel
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_vline scale_colour_manual
#'   labs theme_classic
#' @export
bayesian_tumor_growth <- function(
  df,
  time_column                  = "Day",
  volume_column                = "Volume",
  treatment_column             = "Treatment",
  id_column                    = "ID",
  cage_column                  = NULL,
  dose_column                  = NULL,
  transform                    = c("log", "sqrt", "none"),
  random_effects_specification = c("intercept_only", "slope"),
  reference_group              = NULL,
  prior_strength               = c("skeptical", "weakly_informative", "informative", "diffuse", "manual"),
  prior_b                      = NULL,
  prior_intercept              = NULL,
  prior_sd                     = NULL,
  prior_sigma                  = NULL,
  n_chains                     = 4L,
  n_iter                       = 2000L,
  seed                         = 42L,
  return_model                 = TRUE,
  plots                        = TRUE,
  verbose                      = FALSE,
  necrotic_column              = NULL,
  necrotic_handling            = c("exclude", "covariate", "none")
) {

  # ── Dependency checks ──────────────────────────────────────────────────────
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for Bayesian analysis.\n",
      "Install it with: install.packages('brms')"
    )
  }

  transform                    <- match.arg(transform)
  random_effects_specification <- match.arg(random_effects_specification)
  prior_strength               <- match.arg(prior_strength)
  necrotic_handling            <- match.arg(necrotic_handling)

  # ── Column validation ──────────────────────────────────────────────────────
  required_cols <- c(time_column, volume_column, treatment_column, id_column)
  missing_cols  <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Cage placeholder — mirrors lme4 path so make_mouse_key() always works
  no_cage_mode <- is.null(cage_column) || !cage_column %in% colnames(df)
  if (no_cage_mode) {
    cage_column       <- ".cage_placeholder"
    df[[cage_column]] <- "1"
  }

  # ── Reference group resolution ─────────────────────────────────────────────
  treatment_groups <- unique(as.character(df[[treatment_column]]))
  if (is.null(reference_group)) {
    reference_group <- sort(treatment_groups)[1]
  } else if (!reference_group %in% treatment_groups) {
    stop("Reference group '", reference_group, "' not found in treatment column.")
  }

  # Relevel so reference group is first (controls contrast coding in brms)
  df[[treatment_column]] <- stats::relevel(
    factor(df[[treatment_column]]), ref = reference_group
  )

  # ── Data preparation ───────────────────────────────────────────────────────
  auc_df      <- df          # untransformed copy for growth rates
  analysis_df <- df

  # Handle necrotic observations (before transformation)
  necrosis_summary <- NULL
  include_necrotic_covariate <- FALSE
  if (!is.null(necrotic_column) && necrotic_column %in% colnames(analysis_df)) {
    nec <- tgs_handle_necrosis(analysis_df, necrotic_column, necrotic_handling,
                               treatment_column, id_column, time_column)
    analysis_df      <- nec$data
    necrosis_summary <- nec$necrosis_summary
    include_necrotic_covariate <- necrotic_handling == "covariate"
  }

  if (transform == "log") {
    vol <- analysis_df[[volume_column]]
    vol[vol <= 0] <- min(vol[vol > 0], na.rm = TRUE) / 2
    analysis_df[[volume_column]] <- log(vol)
  } else if (transform == "sqrt") {
    analysis_df[[volume_column]] <- sqrt(analysis_df[[volume_column]])
  }

  # ── Formula ────────────────────────────────────────────────────────────────
  re_term <- if (random_effects_specification == "slope") {
    paste0("(", time_column, " | ", id_column, ")")
  } else {
    paste0("(1 | ", id_column, ")")
  }

  fixed_part <- paste(volume_column, "~", treatment_column, "*", time_column)
  if (isTRUE(include_necrotic_covariate)) {
    fixed_part <- paste(fixed_part, "+ necrotic_cov_flag")
  }
  brms_formula <- stats::as.formula(paste(fixed_part, "+", re_term))

  # ── Prior specification ────────────────────────────────────────────────────
  if (prior_strength == "manual") {
    if (any(is.null(c(prior_b, prior_intercept, prior_sd, prior_sigma)))) {
      stop("When prior_strength = 'manual', all four prior_* arguments must be supplied.")
    }
    selected_priors <- c(
      brms::prior_string(prior_b,         class = "b"),
      brms::prior_string(prior_intercept, class = "Intercept"),
      brms::prior_string(prior_sd,        class = "sd"),
      brms::prior_string(prior_sigma,     class = "sigma")
    )
  } else {
    b_sd     <- switch(prior_strength,
      skeptical          = 0.25,
      weakly_informative = 1,
      informative        = 0.5,
      diffuse            = 2.5
    )
    exp_rate <- switch(prior_strength,
      skeptical          = 2,
      weakly_informative = 1,
      informative        = 2,
      diffuse            = 0.5
    )
    selected_priors <- c(
      brms::prior_string(paste0("normal(0, ", b_sd,        ")"), class = "b"),
      brms::prior_string(paste0("normal(0, ", b_sd * 2.5,  ")"), class = "Intercept"),
      brms::prior_string(paste0("exponential(", exp_rate,  ")"), class = "sd"),
      brms::prior_string(paste0("exponential(", exp_rate,  ")"), class = "sigma")
    )
  }

  # ── Fit model ──────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    message("Fitting Bayesian LMM via brms (",
            n_chains, " chains × ", n_iter, " iter)...")
  }

  model <- brms::brm(
    formula = brms_formula,
    data    = analysis_df,
    prior   = selected_priors,
    chains  = as.integer(n_chains),
    iter    = as.integer(n_iter),
    seed    = as.integer(seed),
    silent  = if (isTRUE(verbose)) 0L else 2L,
    refresh = if (isTRUE(verbose)) 100L else 0L
  )

  # ── Posterior summary (fixed effects) ──────────────────────────────────────
  brms_smry     <- summary(model)
  fixed_df      <- as.data.frame(brms_smry$fixed)
  fixed_df      <- cbind(Parameter = rownames(fixed_df), fixed_df,
                         stringsAsFactors = FALSE)
  rownames(fixed_df) <- NULL
  # Rename CrI columns to consistent names
  names(fixed_df)[names(fixed_df) == "l-95% CI"] <- "Lower_95_CrI"
  names(fixed_df)[names(fixed_df) == "u-95% CI"] <- "Upper_95_CrI"
  posterior_summary <- fixed_df

  # ── MCMC diagnostics ───────────────────────────────────────────────────────
  mcmc_diagnostics <- data.frame(
    Parameter = posterior_summary$Parameter,
    Rhat      = round(posterior_summary$Rhat,     4),
    Bulk_ESS  = round(posterior_summary$Bulk_ESS, 0),
    Tail_ESS  = round(posterior_summary$Tail_ESS, 0),
    Converged = posterior_summary$Rhat <= 1.01,
    stringsAsFactors = FALSE
  )

  # ── Treatment effects and pairwise comparisons via emmeans ─────────────────
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
        warning("emmeans on brmsfit failed: ", conditionMessage(e),
                "\nTreatment effects not computed.")
        NULL
      }
    )

    if (!is.null(emm)) {
      emm_df <- as.data.frame(summary(emm, point.est = stats::median))

      # Column names differ between frequentist (lower.CL) and Bayesian (lower.HPD)
      lower_col <- if ("lower.HPD" %in% names(emm_df)) "lower.HPD" else "lower.CL"
      upper_col <- if ("upper.HPD" %in% names(emm_df)) "upper.HPD" else "upper.CL"

      treatment_effects <- data.frame(
        Group         = as.character(emm_df[[treatment_column]]),
        Adjusted_Mean = round(emm_df$emmean,          3),
        SE            = round(emm_df$SE,               3),
        DF            = NA_real_,
        Lower_CL      = round(emm_df[[lower_col]],     3),
        Upper_CL      = round(emm_df[[upper_col]],     3),
        Note          = ifelse(
          as.character(emm_df[[treatment_column]]) == reference_group,
          "Reference group", ""
        ),
        stringsAsFactors = FALSE
      )

      # Reference group row first
      ref_idx <- which(treatment_effects$Group == reference_group)
      if (length(ref_idx) > 0 && ref_idx[1] > 1L) {
        treatment_effects <- rbind(
          treatment_effects[ref_idx[1], ],
          treatment_effects[-ref_idx[1], ]
        )
        rownames(treatment_effects) <- NULL
      }

      # Pairwise contrasts vs reference
      pc <- tryCatch(
        emmeans::contrast(emm, method = "trt.vs.ctrl"),
        error = function(e) NULL
      )
      if (!is.null(pc)) {
        pc_df <- as.data.frame(summary(pc, point.est = stats::median))
        lower_pc <- if ("lower.HPD" %in% names(pc_df)) "lower.HPD" else "lower.CL"
        upper_pc <- if ("upper.HPD" %in% names(pc_df)) "upper.HPD" else "upper.CL"
        pairwise_comparisons <- data.frame(
          contrast  = as.character(pc_df$contrast),
          estimate  = round(pc_df$estimate,    4),
          SE        = if ("SE" %in% names(pc_df)) round(pc_df$SE, 4) else NA_real_,
          Lower_CL  = round(pc_df[[lower_pc]], 4),
          Upper_CL  = round(pc_df[[upper_pc]], 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ── Growth rates (per-animal OLS log-linear, same as lme4 path) ───────────
  growth_rates <- tgs_compute_growth_rates(
    auc_df, treatment_column, id_column,
    cage_column, time_column, volume_column
  )

  # ── Data summary ───────────────────────────────────────────────────────────
  data_summary <- tgs_compute_summary(
    analysis_df, treatment_column, time_column, volume_column
  )

  # ── Plots ──────────────────────────────────────────────────────────────────
  pp_check_plot           <- NULL
  posterior_dist_plot     <- NULL
  credible_intervals_plot <- NULL
  mcmc_trace_plot         <- NULL

  if (isTRUE(plots)) {

    # Posterior predictive check
    pp_check_plot <- tryCatch(
      brms::pp_check(model, type = "dens_overlay", ndraws = 50),
      error = function(e) NULL
    )

    # Treatment-parameter posterior densities and trace plots
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      draws_arr <- tryCatch(brms::as.array(model), error = function(e) NULL)

      if (!is.null(draws_arr)) {
        all_pars <- dimnames(draws_arr)$variable
        tx_pars  <- grep(
          paste0("^b_", gsub("([.^$*+?()\\[\\]{}|])", "\\\\\\1", treatment_column)),
          all_pars, value = TRUE
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

    # Credible intervals forest plot from treatment_effects
    if (!is.null(treatment_effects) && nrow(treatment_effects) > 0) {
      te      <- treatment_effects
      te$Group <- factor(te$Group, levels = rev(te$Group))
      ref_val  <- te$Adjusted_Mean[as.character(te$Group) == reference_group]
      is_ref   <- as.character(te$Group) == reference_group

      credible_intervals_plot <- ggplot2::ggplot(
          te,
          ggplot2::aes(
            x    = .data[["Adjusted_Mean"]],
            y    = .data[["Group"]],
            xmin = .data[["Lower_CL"]],
            xmax = .data[["Upper_CL"]],
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
          subtitle = paste0("Posterior medians at mean study day (", transform, " scale)"),
          x        = paste0("Estimated marginal mean (", transform, " volume)"),
          y        = NULL
        ) +
        ggplot2::theme_classic(base_size = 14)
    }
  }

  # ── Analysis summary metadata ──────────────────────────────────────────────
  re_label <- if (random_effects_specification == "slope") {
    paste0("(", time_column, " | ", id_column, ")")
  } else {
    paste0("(1 | ", id_column, ")")
  }

  analysis_summary <- list(
    analysis_type = "Bayesian Linear Mixed-Effects Model (brms)",
    data_description = list(
      subjects         = length(unique(make_mouse_key(
        auc_df[[id_column]],
        auc_df[[treatment_column]],
        auc_df[[cage_column]]
      ))),
      treatment_groups = length(unique(auc_df[[treatment_column]])),
      time_points      = length(unique(auc_df[[time_column]])),
      reference_group  = reference_group
    ),
    model_specification = list(
      fixed_effects  = paste(volume_column, "~",
                             treatment_column, "*", time_column),
      random_effects = re_label,
      prior_strength = prior_strength,
      transform      = transform
    ),
    methods = list(
      engine          = paste0("brms (", n_chains, " chains × ",
                               n_iter, " iterations, seed = ", seed, ")"),
      prior_b         = if (prior_strength == "manual") prior_b         else paste0("normal(0, ", b_sd, ")"),
      prior_intercept = if (prior_strength == "manual") prior_intercept else paste0("normal(0, ", round(b_sd * 2.5, 2), ")"),
      prior_sd        = if (prior_strength == "manual") prior_sd        else paste0("exponential(", exp_rate, ")"),
      prior_sigma     = if (prior_strength == "manual") prior_sigma     else paste0("exponential(", exp_rate, ")"),
      treatment_effects_note = paste0(
        "Estimated marginal means and 95 % HPD credible intervals ",
        "at mean study day (day ", round(mean(analysis_df[[time_column]]), 1),
        ") via emmeans."
      )
    )
  )

  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    model                   = if (isTRUE(return_model)) model else NULL,
    model_type_used         = "bayes",
    transform_used          = transform,
    summary                 = analysis_summary,
    posterior_summary       = posterior_summary,
    treatment_effects       = treatment_effects,
    pairwise_comparisons    = pairwise_comparisons,
    mcmc_diagnostics        = mcmc_diagnostics,
    pp_check_plot           = pp_check_plot,
    posterior_dist_plot     = posterior_dist_plot,
    credible_intervals_plot = credible_intervals_plot,
    mcmc_trace_plot         = mcmc_trace_plot,
    growth_rates            = growth_rates,
    data_summary            = data_summary,
    necrosis_summary        = necrosis_summary
  )
}
