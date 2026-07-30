# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Bayesian Parametric Survival Analysis for Mouse Tumor Experiments
#'
#' Fits a Bayesian parametric survival model via \pkg{brms} / Stan.
#' Complements the frequentist \code{\link{survival_statistics}} with full
#' posterior distributions, 95 \% HPD credible intervals, optional cage-level
#' frailty, and survival-curve plots. Four parametric families cover the most
#' common preclinical endpoint distributions.
#'
#' @param df Data frame containing one row per animal (or one row per event
#'   time if the animal experienced the event; censored animals contribute one
#'   row with their last observed day).
#' @param time_column Column name for time-to-event (days). Default
#'   \code{"Day"}.
#' @param event_column Column name for event indicator: \code{1} = event
#'   observed, \code{0} = right-censored. Default \code{"Survival_Censor"}.
#'   Matches the convention used by \code{\link{survival_statistics}}.
#' @param treatment_column Column name for treatment group. Default
#'   \code{"Treatment"}.
#' @param id_column Column name for animal ID. Default \code{"ID"}.
#' @param cage_column Column name for cage ID, or \code{NULL} (default).
#'   When supplied and \code{include_cage_effect = TRUE}, a cage-level
#'   frailty (\code{(1 | cage)}) is added to account for between-cage
#'   heterogeneity.
#' @param dose_column Reserved for future use; currently ignored.
#' @param family Parametric distribution for event times:
#'   \describe{
#'     \item{\code{"weibull"}}{(default) Flexible — handles increasing,
#'       decreasing, and constant hazards. Interpretable as both AFT (time
#'       ratio) and PH (hazard ratio).}
#'     \item{\code{"lognormal"}}{Log-normal AFT. Hazard rises then falls —
#'       biologically plausible for treated tumours.}
#'     \item{\code{"exponential"}}{Constant hazard. Special case of Weibull
#'       with shape = 1; use when a memoryless process is plausible.}
#'     \item{\code{"gamma"}}{Gamma AFT. More flexible tail than Weibull; good
#'       alternative when Weibull does not fit well.}
#'   }
#' @param include_cage_effect Logical. When \code{TRUE} (default) and
#'   \code{cage_column} is supplied, a cage-level random intercept (frailty)
#'   is added: \code{(1 | cage)}.
#' @param reference_group Treatment group used as the reference level.
#'   Auto-detected alphabetically if \code{NULL}.
#' @param prior_strength Prior preset for fixed-effect coefficients (log time
#'   ratios):
#'   \describe{
#'     \item{\code{"skeptical"}}{(default) \eqn{b \sim N(0, 0.25)}.
#'       Requires strong data to support large time-ratio differences.}
#'     \item{\code{"weakly_informative"}}{\eqn{b \sim N(0, 1)}.}
#'     \item{\code{"informative"}}{\eqn{b \sim N(0, 0.5)}.}
#'     \item{\code{"diffuse"}}{\eqn{b \sim N(0, 2.5)}.}
#'     \item{\code{"manual"}}{Specify all four \code{prior_*} arguments
#'       directly.}
#'   }
#' @param prior_b brms prior string for fixed-effect coefficients (class
#'   \code{"b"}). Only used when \code{prior_strength = "manual"}.
#' @param prior_intercept brms prior string for the model intercept. Only used
#'   when \code{prior_strength = "manual"}.
#' @param prior_sd brms prior string for the frailty standard deviation (class
#'   \code{"sd"}). Only used when \code{prior_strength = "manual"}.
#' @param prior_aux brms prior string for the auxiliary parameter: \code{shape}
#'   for Weibull and Gamma; \code{sigma} for log-normal; ignored for
#'   exponential. Only used when \code{prior_strength = "manual"}.
#' @param n_chains Number of MCMC chains. Default \code{4}.
#' @param n_warmup Warm-up (burn-in) iterations per chain. Default \code{1000}.
#' @param n_iter Post-warmup draws per chain. Default \code{500}.
#' @param seed Integer random seed for reproducibility. Default \code{42}.
#' @param return_model Logical. Return the \code{brmsfit} object? Default
#'   \code{TRUE}.
#' @param plots Logical. Pre-compute diagnostic and survival-curve plots?
#'   Default \code{TRUE}.
#' @param verbose Logical. Show Stan compilation and sampling messages? Default
#'   \code{FALSE}.
#' @param backend brms backend: \code{"rstan"} (default) or \code{"cmdstanr"}.
#'   See \code{\link{bayesian_tumor_growth}} for details.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{model}}{\code{brmsfit} object, or \code{NULL} when
#'     \code{return_model = FALSE}.}
#'   \item{\code{model_type_used}}{Character \code{"bayes_survival"}.}
#'   \item{\code{family_used}}{The parametric family fitted.}
#'   \item{\code{frailty_used}}{Logical — whether cage frailty was modelled.}
#'   \item{\code{summary}}{Named list of analysis metadata.}
#'   \item{\code{treatment_effects}}{Data frame with columns \code{Group},
#'     \code{Time_Ratio}, \code{Lower_CrI}, \code{Upper_CrI}, \code{HR}
#'     (hazard ratio; Weibull and exponential only, \code{NA} otherwise),
#'     \code{Median_Survival}, \code{Events}, \code{Total}, \code{Event_Rate},
#'     \code{Note}. Output schema mirrors \code{\link{survival_statistics}}.}
#'   \item{\code{posterior_summary}}{Data frame of fixed-effect posterior
#'     medians, 2.5 \%–97.5 \% CrI, Rhat, Bulk_ESS, Tail_ESS.}
#'   \item{\code{mcmc_diagnostics}}{Per-parameter Rhat, ESS, and convergence
#'     flags (Rhat > 1.01 flagged as not converged).}
#'   \item{\code{survival_data}}{Data frame with \code{Time}, \code{Event},
#'     \code{Treatment} — mirrors \code{\link{survival_statistics}}.}
#'   \item{\code{pp_check_plot}}{Posterior predictive density overlay, or
#'     \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{posterior_dist_plot}}{Treatment-parameter posterior density
#'     areas (bayesplot), or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{prior_posterior_plot}}{Density overlay of prior (grey) and
#'     posterior (blue) for each treatment-effect coefficient. A vertical
#'     dashed line marks zero (no effect). \code{NULL} when
#'     \code{plots = FALSE}.}
#'   \item{\code{mcmc_trace_plot}}{MCMC trace plot, or \code{NULL} when
#'     \code{plots = FALSE}.}
#'   \item{\code{survival_curve_plot}}{Parametric survival curves with 95 \%
#'     posterior credible bands overlaid on Kaplan-Meier step functions, or
#'     \code{NULL} when \code{plots = FALSE}.}
#' }
#'
#' @section Effect interpretation:
#' All families are fitted in the AFT (accelerated failure time)
#' parameterisation via \pkg{brms}. Treatment coefficients are on the
#' log(mean survival time)
#' scale; exponentiating gives \strong{time ratios} (TR): TR > 1 means the
#' treated group survives longer. For Weibull and exponential families, an
#' approximate \strong{hazard ratio} (HR = TR\eqn{^{-\hat{\kappa}}}, where
#' \eqn{\hat{\kappa}} is the posterior-median Weibull shape) is also reported;
#' HR < 1 indicates a protective treatment effect. Lognormal and gamma are
#' not proportional-hazards models; their \code{HR} column is \code{NA}.
#'
#' @section Frailty:
#' Cage-level frailty \code{(1 | cage)} models unobserved between-cage
#' heterogeneity (e.g., microbiome, handling, cohort). The frailty variance is
#' reported in the \code{posterior_summary} random-effects block.
#'
#' @references
#' Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models
#' using Stan. \emph{Journal of Statistical Software}, 80(1), 1–28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @seealso \code{\link{survival_statistics}},
#'   \code{\link{bayesian_tumor_growth}}
#'
#' @importFrom stats as.formula relevel median quantile qgamma pgamma pnorm
#'   setNames
#' @importFrom ggplot2 ggplot aes geom_density geom_line geom_ribbon geom_step
#'   scale_fill_manual scale_colour_manual scale_y_continuous
#'   facet_wrap labs theme theme_classic
#' @importFrom survival Surv survfit
#' @export
bayesian_survival <- function(
  df,
  time_column      = "Day",
  event_column     = "Survival_Censor",
  treatment_column = "Treatment",
  id_column        = "ID",
  cage_column      = NULL,
  dose_column      = NULL,
  family              = c("weibull", "lognormal", "exponential", "gamma"),
  include_cage_effect = TRUE,
  reference_group  = NULL,
  prior_strength   = c("skeptical", "weakly_informative", "informative",
                       "diffuse", "manual"),
  prior_b          = NULL,
  prior_intercept  = NULL,
  prior_sd         = NULL,
  prior_aux        = NULL,
  n_chains         = 4L,
  n_warmup         = 1000L,
  n_iter           = 500L,
  seed             = 42L,
  return_model     = TRUE,
  plots            = TRUE,
  verbose          = FALSE,
  backend          = c("rstan", "cmdstanr"),
  priors           = NULL,
  mcmc             = NULL
) {

  # ── Dependency check ───────────────────────────────────────────────────────
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package 'brms' is required.\nInstall with: install.packages('brms')")
  }

  family         <- match.arg(family)
  prior_strength <- match.arg(prior_strength)
  backend        <- resolve_brms_backend(backend)

  # CODE_REVIEW.md Round 2 D.2 — accept `priors = tg_priors()` /
  # `mcmc = tg_mcmc()` override objects. Note: `priors$sigma` is mapped
  # to `prior_aux` for survival families that use an auxiliary parameter
  # (shape for weibull, etc.); the user can also still set `prior_aux`
  # directly.
  .p <- .resolve_priors(priors, prior_strength, prior_b, prior_intercept,
                        prior_sd, prior_aux)
  prior_strength  <- .p$strength
  prior_b         <- .p$b
  prior_intercept <- .p$intercept
  prior_sd        <- .p$sd
  prior_aux       <- .p$sigma
  .m <- .resolve_mcmc(mcmc, n_chains, n_warmup, n_iter, seed, backend)
  n_chains <- .m$chains
  n_warmup <- .m$warmup
  n_iter   <- .m$iter
  seed     <- .m$seed
  backend  <- .m$backend

  # ── Column validation ──────────────────────────────────────────────────────
  required_cols <- c(time_column, event_column, treatment_column, id_column)
  missing_cols  <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!all(df[[event_column]] %in% c(0L, 1L, 0, 1, NA))) {
    stop("'", event_column,
         "' must contain only 0 (censored) and 1 (event).")
  }

  # ── Cage setup ─────────────────────────────────────────────────────────────
  cage_setup   <- setup_cage_column(df, cage_column)
  df           <- cage_setup$df
  cage_column  <- cage_setup$cage_column
  no_cage_mode <- cage_setup$no_cage_mode
  use_cage_re  <- isTRUE(include_cage_effect) && !no_cage_mode

  # ── Reference group ────────────────────────────────────────────────────────
  treatment_groups <- unique(as.character(df[[treatment_column]]))
  if (is.null(reference_group)) {
    reference_group <- sort(treatment_groups)[1]
  } else if (!reference_group %in% treatment_groups) {
    stop("Reference group '", reference_group,
         "' not found in treatment column.")
  }

  df[[treatment_column]] <- stats::relevel(
    factor(df[[treatment_column]]), ref = reference_group
  )
  treatment_levels <- levels(df[[treatment_column]])

  # ── Analysis data frame ────────────────────────────────────────────────────
  # brms cens() convention: 0 = observed event, 1 = right-censored
  # Input convention (matching survival package): 1 = event, 0 = censored
  analysis_df <- df
  analysis_df$.brms_cens <- 1L - as.integer(analysis_df[[event_column]])

  # ── Formula ────────────────────────────────────────────────────────────────
  re_term  <- if (use_cage_re) paste0(" + (1 | ", cage_column, ")") else ""
  brms_lhs <- paste0(time_column, " | cens(.brms_cens)")
  brms_formula <- stats::as.formula(
    paste(brms_lhs, "~", treatment_column, re_term)
  )

  # ── brms family object ─────────────────────────────────────────────────────
  brms_family <- switch(family,
    weibull     = brms::weibull(),
    lognormal   = brms::lognormal(),
    exponential = brms::exponential(),
    gamma       = brms::Gamma(link = "log")
  )

  # ── Prior specification ────────────────────────────────────────────────────
  if (prior_strength == "manual") {
    if (any(is.null(c(prior_b, prior_intercept, prior_sd)))) {
      stop("When prior_strength = 'manual', prior_b, prior_intercept, and ",
           "prior_sd must be supplied.")
    }
    # CODE_REVIEW.md R3.32 — only declare a class = "sd" prior when the model
    # actually HAS a random effect. brms rejects a prior that matches no model
    # parameter ("The following priors do not correspond to any model
    # parameter: sd ~ ..."), so with include_cage_effect = FALSE every fit
    # errored. Masked until now because these tests skipped whenever brms was
    # absent from the dev environment.
    selected_priors <- c(
      brms::prior_string(prior_b,         class = "b"),
      brms::prior_string(prior_intercept, class = "Intercept")
    )
    if (use_cage_re) {
      selected_priors <- c(selected_priors,
                           brms::prior_string(prior_sd, class = "sd"))
    }
    if (!is.null(prior_aux)) {
      aux_class <- if (family == "lognormal") "sigma" else "shape"
      selected_priors <- c(
        selected_priors,
        brms::prior_string(prior_aux, class = aux_class)
      )
    }
  } else {
    pp       <- bayes_prior_params(prior_strength)
    b_sd     <- pp$b_sd
    exp_rate <- pp$exp_rate
    # CODE_REVIEW.md R3.8 / G.5 — the Intercept here is on the log-time scale
    # (an AFT model), so for a study running to day 35 it centres near
    # log(35) ~ 3.6. A fixed normal(0, b_sd * 2.5) prior put almost no mass
    # there. Scale it to the observed event times. The `b` coefficients are
    # log-time ratios, which are unit-free, so the ladder applies directly and
    # there is no time covariate to divide by.
    t_obs <- as.numeric(analysis_df[[time_column]])
    t_obs <- t_obs[is.finite(t_obs) & t_obs > 0]
    log_t_med <- if (length(t_obs)) stats::median(log(t_obs)) else 0
    log_t_mad <- if (length(t_obs)) stats::mad(log(t_obs)) else 1
    if (!is.finite(log_t_mad) || log_t_mad <= 0) log_t_mad <- 1

    selected_priors <- c(
      brms::prior_string(paste0("normal(0, ", b_sd, ")"),
                         class = "b"),
      brms::prior_string(
        paste0("normal(", round(log_t_med, 4), ", ",
               round(2.5 * log_t_mad, 4), ")"),
        class = "Intercept")
    )
    # See R3.32 above — a "sd" prior with no random effect makes brms error.
    if (use_cage_re) {
      selected_priors <- c(selected_priors,
                           brms::prior_string(paste0("exponential(", exp_rate, ")"),
                                              class = "sd"))
    }
    if (family == "lognormal") {
      aux_class <- "sigma"
    } else if (family != "exponential") {
      aux_class <- "shape"
    } else {
      aux_class <- NULL
    }
    if (!is.null(aux_class)) {
      selected_priors <- c(
        selected_priors,
        brms::prior_string(paste0("exponential(", exp_rate, ")"),
                           class = aux_class)
      )
    }
  }

  # ── Fit model ──────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    message("Fitting Bayesian parametric survival model (",
            family, ") via brms (",
            n_chains, " chains × ", n_iter, " post-warmup draws, ",
            n_warmup, " warmup)...")
  }

  model <- brms::brm(
    formula      = brms_formula,
    data         = analysis_df,
    family       = brms_family,
    prior        = selected_priors,
    sample_prior = "yes",
    chains       = as.integer(n_chains),
    cores        = as.integer(n_chains),
    iter         = as.integer(n_warmup + n_iter),
    warmup       = as.integer(n_warmup),
    seed         = as.integer(seed),
    backend      = backend,
    silent       = if (isTRUE(verbose)) 0L else 2L,
    refresh      = if (isTRUE(verbose)) 100L else 0L
  )

  # ── Posterior summary (fixed effects) ──────────────────────────────────────
  brms_smry <- summary(model)
  fixed_df  <- as.data.frame(brms_smry$fixed)
  fixed_df  <- cbind(Parameter = rownames(fixed_df), fixed_df,
                     stringsAsFactors = FALSE)
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
  ppc_coverage     <- bayes_ppc_coverage(model)

  # ── Treatment effects table ────────────────────────────────────────────────
  treatment_effects <- bs_build_treatment_table(
    model, analysis_df, treatment_column, treatment_levels, reference_group,
    event_column, time_column, id_column, cage_column, family
  )

  # ── Plots ──────────────────────────────────────────────────────────────────
  pp_check_plot        <- NULL
  posterior_dist_plot  <- NULL
  prior_posterior_plot <- NULL
  mcmc_trace_plot      <- NULL
  survival_curve_plot  <- NULL

  if (isTRUE(plots)) {

    pp_check_plot <- tryCatch(
      brms::pp_check(model, type = "dens_overlay", ndraws = 50),
      error = function(e) NULL
    )

    if (requireNamespace("bayesplot", quietly = TRUE)) {
      draws_arr <- tryCatch(posterior::as_draws_array(model),
                            error = function(e) NULL)
      if (!is.null(draws_arr)) {
        all_pars <- dimnames(draws_arr)$variable
        safe_tx  <- gsub("([.^$*+?()\\[\\]{}|])", "\\\\\\1",
                         treatment_column, perl = TRUE)
        tx_pars  <- grep(paste0("^b_", safe_tx), all_pars, value = TRUE)
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

    # Prior vs posterior overlay
    prior_posterior_plot <- tryCatch(
      bayes_prior_posterior_plot(model, treatment_column),
      error = function(e) NULL
    )

    survival_curve_plot <- tryCatch(
      bs_survival_curves_plot(
        model, analysis_df, treatment_column, treatment_levels, reference_group,
        time_column, event_column, family
      ),
      error = function(e) NULL
    )
  }

  # ── Analysis summary metadata ──────────────────────────────────────────────
  frailty_term_str <- if (use_cage_re) paste0("(1 | ", cage_column, ")") else
    "none"
  prior_b_label <- if (prior_strength == "manual") prior_b else
    paste0("normal(0, ", b_sd, ")")
  prior_int_label <- if (prior_strength == "manual") prior_intercept else
    paste0("normal(0, ", round(b_sd * 2.5, 2), ")")
  analysis_summary <- list(
    analysis_type = paste0("Bayesian Parametric Survival (", family, ", brms)"),
    data_description = list(
      subjects         = length(unique(df[[id_column]])),
      treatment_groups = length(treatment_levels),
      total_events     = sum(df[[event_column]] == 1L, na.rm = TRUE),
      reference_group  = reference_group
    ),
    model_specification = list(
      formula        = deparse(brms_formula),
      family         = family,
      frailty        = use_cage_re,
      frailty_term   = frailty_term_str,
      prior_strength = prior_strength
    ),
    methods = list(
      engine          = paste0("brms (", n_chains, " chains × ",
                               n_iter, " draws + ", n_warmup,
                               " warmup, seed = ", seed, ")"),
      effect_type     = "Time_Ratio (AFT parameterisation)",
      prior_b         = prior_b_label,
      prior_intercept = prior_int_label
    )
  )

  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    model               = if (isTRUE(return_model)) model else NULL,
    model_type_used     = "bayes_survival",
    family_used         = family,
    frailty_used        = use_cage_re,
    summary             = analysis_summary,
    treatment_effects   = treatment_effects,
    posterior_summary   = posterior_summary,
    mcmc_diagnostics    = mcmc_diagnostics,
    nuts_diagnostics    = nuts_diagnostics,
    loo_diagnostics     = loo_diagnostics,
    bayes_R2            = bayes_r2,
    ppc_coverage        = ppc_coverage,
    survival_data       = data.frame(
      Time      = df[[time_column]],
      Event     = df[[event_column]],
      Treatment = df[[treatment_column]]
    ),
    pp_check_plot        = pp_check_plot,
    posterior_dist_plot  = posterior_dist_plot,
    prior_posterior_plot = prior_posterior_plot,
    mcmc_trace_plot      = mcmc_trace_plot,
    survival_curve_plot  = survival_curve_plot
  )
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Build treatment effects table from brms posterior
#' @noRd
bs_build_treatment_table <- function(
  model, analysis_df, treatment_column, treatment_levels, reference_group,
  event_column, time_column, id_column, cage_column, family
) {
  post <- tryCatch(brms::as_draws_df(model), error = function(e) NULL)
  if (is.null(post)) return(NULL)

  intercept_draws <- if ("b_Intercept" %in% names(post))
    post$b_Intercept else NULL
  shape_draws <- if ("shape" %in% names(post)) post$shape else NULL
  sigma_draws <- if ("sigma" %in% names(post)) post$sigma else NULL

  # Map non-reference treatment levels → posterior draw column names
  # Uses positional matching: factor level order matches coefficient order
  fixef_names  <- rownames(brms::fixef(model))
  tx_coef_names <- fixef_names[startsWith(fixef_names, treatment_column)]
  non_ref_levels <- treatment_levels[treatment_levels != reference_group]

  # Positional map (safe because both share the same factor level ordering)
  if (length(tx_coef_names) == length(non_ref_levels)) {
    level_to_coef <- stats::setNames(tx_coef_names, non_ref_levels)
  } else {
    # Fallback: name-based matching via make.names sanitisation
    level_to_coef <- stats::setNames(
      vapply(non_ref_levels, function(lvl) {
        expected <- paste0(treatment_column, stats::make.names(lvl))
        if (expected %in% tx_coef_names) expected else NA_character_
      }, character(1)),
      non_ref_levels
    )
  }

  # Initialise results table
  results <- data.frame(
    Group           = treatment_levels,
    Time_Ratio      = NA_real_,
    Lower_CrI       = NA_real_,
    Upper_CrI       = NA_real_,
    HR              = NA_real_,
    Median_Survival = NA_real_,
    Events          = NA_integer_,
    Total           = NA_integer_,
    Event_Rate      = NA_real_,
    Note            = "",
    stringsAsFactors = FALSE
  )

  # Reference group
  ref_idx <- which(results$Group == reference_group)
  results$Time_Ratio[ref_idx] <- 1
  results$Lower_CrI[ref_idx]  <- 1
  results$Upper_CrI[ref_idx]  <- 1
  results$HR[ref_idx]          <- 1
  results$Note[ref_idx]        <- "Reference group"
  if (!is.null(intercept_draws)) {
    results$Median_Survival[ref_idx] <- bs_median_from_draws(
      intercept_draws, shape_draws, sigma_draws, family
    )
  }

  # Non-reference groups
  for (grp in non_ref_levels) {
    idx      <- which(results$Group == grp)
    coef_col <- paste0("b_", level_to_coef[[grp]])

    if (!is.na(level_to_coef[[grp]]) && coef_col %in% names(post)) {
      coef_draws <- post[[coef_col]]
      tr_draws   <- exp(coef_draws)

      results$Time_Ratio[idx] <- round(stats::median(tr_draws), 4)
      results$Lower_CrI[idx]  <- round(stats::quantile(tr_draws, 0.025), 4)
      results$Upper_CrI[idx]  <- round(stats::quantile(tr_draws, 0.975), 4)

      # Hazard ratio: only meaningful for PH-compatible families
      if (family == "weibull" && !is.null(shape_draws)) {
        hr_draws <- exp(-shape_draws * coef_draws)
        results$HR[idx] <- round(stats::median(hr_draws), 4)
      } else if (family == "exponential") {
        # Exponential is Weibull(shape=1), so HR = 1/TR
        results$HR[idx] <- round(1 / results$Time_Ratio[idx], 4)
      }

      if (!is.null(intercept_draws)) {
        results$Median_Survival[idx] <- bs_median_from_draws(
          intercept_draws + coef_draws, shape_draws, sigma_draws, family
        )
      }
    }
  }

  # Event counts (mirrors survival_statistics approach)
  for (grp in treatment_levels) {
    grp_mask <- as.character(analysis_df[[treatment_column]]) == grp
    grp_data <- analysis_df[grp_mask, , drop = FALSE]
    idx      <- which(results$Group == grp)
    subjects <- unique(grp_data[, c(id_column, cage_column)])
    events   <- sum(grp_data[[event_column]] == 1L, na.rm = TRUE)
    total    <- nrow(subjects)
    results$Events[idx]     <- events
    results$Total[idx]      <- total
    results$Event_Rate[idx] <- if (total > 0)
      round(events / total, 3) else NA_real_
  }

  # Reference group row first
  ref_row <- which(results$Group == reference_group)
  if (length(ref_row) > 0 && ref_row[1] > 1L) {
    results <- rbind(results[ref_row[1], ], results[-ref_row[1], ])
    rownames(results) <- NULL
  }

  results
}


#' Compute posterior-median survival time from linear-predictor draws
#' @noRd
bs_median_from_draws <- function(linpred_draws, shape_draws,
                                 sigma_draws, family) {
  tryCatch({
    med_draws <- switch(family,
      lognormal = {
        # Identity link: linpred = mu (log-scale mean); median = exp(mu)
        exp(linpred_draws)
      },
      exponential = {
        # Log link: linpred = log(mean); median = mean * log(2)
        exp(linpred_draws) * log(2)
      },
      weibull = {
        # log link: linpred = log(mean); convert via Weibull scale parameter
        if (is.null(shape_draws)) return(NA_real_)
        mu    <- exp(linpred_draws)
        scale <- mu / gamma(1 + 1 / shape_draws)
        scale * log(2)^(1 / shape_draws)
      },
      gamma = {
        # log link: linpred = log(mean); rate = shape / mean
        if (is.null(shape_draws)) return(NA_real_)
        mu   <- exp(linpred_draws)
        rate <- shape_draws / mu
        mapply(function(sh, ra) stats::qgamma(0.5, shape = sh, rate = ra),
               shape_draws, rate)
      }
    )
    round(stats::median(med_draws, na.rm = TRUE), 2)
  }, error = function(e) NA_real_)
}


#' Build survival-curve plot with posterior credible bands + KM overlay
#' @noRd
bs_survival_curves_plot <- function(
  model, analysis_df, treatment_column, treatment_levels, reference_group,
  time_column, event_column, family
) {
  post <- tryCatch(brms::as_draws_df(model), error = function(e) NULL)
  if (is.null(post)) return(NULL)

  t_max  <- max(analysis_df[[time_column]], na.rm = TRUE) * 1.05
  t_grid <- seq(0.01, t_max, length.out = 200)

  n_draws  <- min(200L, nrow(post))
  draw_idx <- sample.int(nrow(post), n_draws)
  post_sub <- post[draw_idx, ]

  intercept_d <- if ("b_Intercept" %in% names(post_sub))
    post_sub$b_Intercept else rep(0, n_draws)
  shape_d <- if ("shape" %in% names(post_sub)) post_sub$shape else NULL
  sigma_d <- if ("sigma" %in% names(post_sub)) post_sub$sigma else NULL

  # Map non-reference levels → draw columns (positional)
  fixef_names   <- rownames(brms::fixef(model))
  tx_coef_names <- fixef_names[startsWith(fixef_names, treatment_column)]
  non_ref_levels <- treatment_levels[treatment_levels != reference_group]

  if (length(tx_coef_names) == length(non_ref_levels)) {
    level_to_coef <- stats::setNames(tx_coef_names, non_ref_levels)
  } else {
    level_to_coef <- stats::setNames(
      vapply(non_ref_levels, function(lvl) {
        expected <- paste0(treatment_column, stats::make.names(lvl))
        if (expected %in% tx_coef_names) expected else NA_character_
      }, character(1)),
      non_ref_levels
    )
  }

  # Compute survival curves per group
  curve_list <- lapply(treatment_levels, function(grp) {
    if (grp == reference_group) {
      lp_d <- intercept_d
    } else {
      coef_col <- paste0("b_", level_to_coef[[grp]])
      if (!is.na(level_to_coef[[grp]]) && coef_col %in% names(post_sub)) {
        lp_d <- intercept_d + post_sub[[coef_col]]
      } else {
        lp_d <- intercept_d
      }
    }

    surv_mat <- bs_survival_matrix(t_grid, lp_d, shape_d, sigma_d, family)

    data.frame(
      t         = t_grid,
      median    = apply(surv_mat, 2, stats::median, na.rm = TRUE),
      lower     = apply(surv_mat, 2, stats::quantile,
                        probs = 0.025, na.rm = TRUE),
      upper     = apply(surv_mat, 2, stats::quantile,
                        probs = 0.975, na.rm = TRUE),
      Treatment = grp,
      stringsAsFactors = FALSE
    )
  })
  curves_df <- do.call(rbind, curve_list)

  # Kaplan-Meier data for overlay
  km_df <- tryCatch({
    km_formula <- stats::as.formula(
      paste0("survival::Surv(", time_column, ", ", event_column, ") ~ ",
             treatment_column)
    )
    km_fit  <- survival::survfit(km_formula, data = analysis_df)
    km_smry <- summary(km_fit)
    data.frame(
      t         = km_smry$time,
      surv      = km_smry$surv,
      Treatment = gsub(paste0(treatment_column, "="), "",
                       as.character(km_smry$strata)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)

  # Build plot
  p <- ggplot2::ggplot(
    curves_df,
    ggplot2::aes(colour = .data[["Treatment"]], fill = .data[["Treatment"]])
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = .data[["t"]], ymin = .data[["lower"]],
                   ymax = .data[["upper"]]),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data[["t"]], y = .data[["median"]]),
      linewidth = 0.9
    )

  if (!is.null(km_df)) {
    p <- p + ggplot2::geom_step(
      data    = km_df,
      mapping = ggplot2::aes(x = .data[["t"]], y = .data[["surv"]],
                             colour = .data[["Treatment"]]),
      linetype = "dashed", linewidth = 0.5
    )
  }

  p +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(
      title    = "Bayesian Parametric Survival Curves",
      subtitle = paste0(
        tools::toTitleCase(family),
        " model — solid = posterior median, band = 95 % CrI,",
        " dashed = Kaplan–Meier"
      ),
      x      = time_column,
      y      = "Survival Probability",
      colour = treatment_column,
      fill   = treatment_column
    ) +
    ggplot2::theme_classic(base_size = 14)
}


#' Compute survival probability matrix from posterior draws
#'
#' Returns a matrix with rows = draws and columns = time points.
#' @noRd
bs_survival_matrix <- function(t_grid, linpred_draws, shape_draws,
                               sigma_draws, family) {
  n_draws <- length(linpred_draws)

  switch(family,
    lognormal = {
      # identity link: linpred = mu (log-scale mean)
      if (is.null(sigma_draws)) sigma_draws <- rep(1, n_draws)
      t(mapply(function(mu, sig) {
        stats::pnorm((log(t_grid) - mu) / sig, lower.tail = FALSE)
      }, linpred_draws, sigma_draws))
    },
    exponential = {
      # Log link: linpred = log(mean); S(t) = exp(-t / mean)
      t(sapply(linpred_draws, function(lp) {
        exp(-t_grid / exp(lp))
      }))
    },
    weibull = {
      if (is.null(shape_draws)) shape_draws <- rep(1, n_draws)
      t(mapply(function(lp, sh) {
        mu    <- exp(lp)
        scale <- mu / gamma(1 + 1 / sh)
        exp(-(t_grid / scale)^sh)
      }, linpred_draws, shape_draws))
    },
    gamma = {
      if (is.null(shape_draws)) shape_draws <- rep(1, n_draws)
      t(mapply(function(lp, sh) {
        mu   <- exp(lp)
        rate <- sh / mu
        stats::pgamma(t_grid, shape = sh, rate = rate, lower.tail = FALSE)
      }, linpred_draws, shape_draws))
    }
  )
}
