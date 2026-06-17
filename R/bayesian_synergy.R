# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file


#' Shared brms LMM fit for the two synergy entry points
#'
#' Encapsulates the steps that bayesian_synergy() and
#' bayesian_synergy_over_time() did identically (CODE_REVIEW.md G.6):
#' apply the volume transform, set up cage / id placeholders, build the
#' random-effects term and full brms formula, construct priors, fit the
#' model, and pull the standard suite of diagnostics
#' (posterior_summary, mcmc_diagnostics, nuts_diagnostics, loo_diagnostics,
#' bayes_R2, ppc_coverage).
#'
#' @param analysis_df Data frame already subset to the four analysis groups
#'   and with the Treatment factor levelled.
#' @param volume_column,treatment_column,time_column,id_column,cage_column
#'   Column names; \code{cage_column} may be NULL.
#' @param transform One of "log", "sqrt", "none".
#' @param include_cage_effect Logical.
#' @param re_spec "intercept_only" or "slope".
#' @param prior_str Prior preset name, or "manual".
#' @param manual_priors When \code{prior_str = "manual"}, a list with
#'   \code{prior_b}, \code{prior_intercept}, \code{prior_sd},
#'   \code{prior_sigma}.
#' @param n_chains,n_warmup,n_iter,seed,verbose Standard brms fit knobs.
#' @noRd
bs_fit_synergy_model <- function(analysis_df,
                                  volume_column, treatment_column,
                                  time_column, id_column, cage_column,
                                  transform = c("log", "sqrt", "none"),
                                  include_cage_effect = TRUE,
                                  re_spec = c("intercept_only", "slope"),
                                  prior_str,
                                  manual_priors = NULL,
                                  n_chains, n_warmup, n_iter,
                                  seed, verbose = FALSE,
                                  backend = "rstan") {
  transform <- match.arg(transform)
  re_spec   <- match.arg(re_spec)

  # ── Apply transform ──────────────────────────────────────────────────────
  if (transform == "log") {
    vol           <- analysis_df[[volume_column]]
    positive_vals <- vol[is.finite(vol) & vol > 0]
    if (length(positive_vals) == 0L) {
      stop(
        "No positive values in '", volume_column,
        "' after filtering. Cannot apply log transform."
      )
    }
    vol[vol <= 0] <- min(positive_vals, na.rm = TRUE) / 2
    analysis_df$Response <- log(vol)
  } else if (transform == "sqrt") {
    analysis_df$Response <- sqrt(analysis_df[[volume_column]])
  } else {
    analysis_df$Response <- analysis_df[[volume_column]]
  }

  # ── Cage / id setup ──────────────────────────────────────────────────────
  has_cage <- !is.null(cage_column) && cage_column %in% colnames(analysis_df)
  if (has_cage) {
    analysis_df$cage_var <- analysis_df[[cage_column]]
  } else {
    analysis_df$cage_var <- analysis_df[[id_column]]
  }
  analysis_df$id_var <- paste0(
    analysis_df$cage_var, "__", analysis_df[[id_column]]
  )

  # ── Random-effects term + full formula ───────────────────────────────────
  re_term <- if (re_spec == "slope") {
    if (has_cage && isTRUE(include_cage_effect)) {
      paste0("(", time_column, " | id_var) + (1 | cage_var)")
    } else {
      paste0("(", time_column, " | id_var)")
    }
  } else {
    if (has_cage && isTRUE(include_cage_effect)) {
      "(1 | cage_var/id_var)"
    } else {
      "(1 | id_var)"
    }
  }
  brms_formula <- stats::as.formula(
    paste0(
      "Response ~ ", treatment_column, " * ", time_column, " + ", re_term
    )
  )

  # ── Priors ───────────────────────────────────────────────────────────────
  if (prior_str == "manual") {
    if (is.null(manual_priors$prior_b) ||
          is.null(manual_priors$prior_intercept) ||
          is.null(manual_priors$prior_sd) ||
          is.null(manual_priors$prior_sigma)) {
      stop(
        "prior_strength = 'manual' requires prior_b, prior_intercept, ",
        "prior_sd, and prior_sigma."
      )
    }
    priors <- c(
      brms::prior_string(manual_priors$prior_b,         class = "b"),
      brms::prior_string(manual_priors$prior_intercept, class = "Intercept"),
      brms::prior_string(manual_priors$prior_sd,        class = "sd"),
      brms::prior_string(manual_priors$prior_sigma,     class = "sigma")
    )
  } else {
    pp       <- bayes_prior_params(prior_str)
    b_sd     <- pp$b_sd
    exp_rate <- pp$exp_rate
    priors   <- c(
      brms::prior_string(paste0("normal(0, ", b_sd,       ")"),
                         class = "b"),
      brms::prior_string(paste0("normal(0, ", b_sd * 2.5, ")"),
                         class = "Intercept"),
      brms::prior_string(paste0("exponential(", exp_rate, ")"),
                         class = "sd"),
      brms::prior_string(paste0("exponential(", exp_rate, ")"),
                         class = "sigma")
    )
  }

  # ── Fit ──────────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) message("Fitting brms model …")

  # Stan warnings (divergent transitions, max_treedepth, low ESS, Rhat) are
  # exactly the signals we want to surface — they should NOT be suppressed.
  # See nuts_diagnostics for the programmatic equivalents.
  model <- brms::brm(
    formula      = brms_formula,
    data         = analysis_df,
    prior        = priors,
    chains       = n_chains,
    cores        = n_chains,
    iter         = n_warmup + n_iter,
    warmup       = n_warmup,
    seed         = seed,
    backend      = backend,
    sample_prior = "yes",
    silent       = if (isTRUE(verbose)) 0L else 2L,
    refresh      = if (isTRUE(verbose)) 100L else 0L
  )

  posterior_summary <- build_posterior_summary(model)
  mcmc_diagnostics  <- make_mcmc_diagnostics(
    posterior_summary, total_draws = n_chains * n_iter
  )
  nuts_diagnostics  <- make_nuts_diagnostics(model)
  loo_diagnostics   <- bayes_loo(model)
  bayes_r2          <- bayes_r2_summary(model)
  ppc_coverage      <- bayes_ppc_coverage(model)

  list(
    model             = model,
    analysis_df       = analysis_df,
    has_cage          = has_cage,
    re_term           = re_term,
    brms_formula      = brms_formula,
    posterior_summary = posterior_summary,
    mcmc_diagnostics  = mcmc_diagnostics,
    nuts_diagnostics  = nuts_diagnostics,
    loo_diagnostics   = loo_diagnostics,
    bayes_R2          = bayes_r2,
    ppc_coverage      = ppc_coverage
  )
}

#' Bayesian Drug Combination Synergy Analysis
#'
#' Fits a Bayesian linear mixed-effects model to longitudinal tumor volume data
#' from a four-group drug combination experiment and propagates full posterior
#' uncertainty through the Bliss Independence and Loewe Additivity synergy
#' metrics via draw-wise arithmetic. Provides posterior probability of synergy
#' for both metrics.
#'
#' @param df Data frame with longitudinal tumor volume data containing at
#'   least the four groups named by \code{drug_a_name}, \code{drug_b_name},
#'   \code{combo_name}, and \code{control_name}.
#' @param volume_column Column name for tumor volume (mm³). Default
#'   \code{"Volume"}.
#' @param time_column Column name for study day (numeric). Default
#'   \code{"Day"}.
#' @param treatment_column Column name for treatment group. Default
#'   \code{"Treatment"}.
#' @param id_column Column name for individual animal ID. Default
#'   \code{"ID"}.
#' @param cage_column Column name for cage ID, or \code{NULL} (default).
#' @param drug_a_name Name of the first single-agent treatment group.
#' @param drug_b_name Name of the second single-agent treatment group.
#' @param combo_name Name of the combination treatment group.
#' @param control_name Name of the vehicle/control group. Default
#'   \code{"Control"}.
#' @param endpoint_day Study day at which synergy is evaluated. \code{NULL}
#'   (default) uses the last observed day in \code{df}.
#' @param transform Volume transformation before modelling: \code{"log"}
#'   (default), \code{"sqrt"}, or \code{"none"}.
#' @param random_effects_specification Random-effects structure for the brms
#'   model: \code{"intercept_only"} (default, \code{(1|ID)}) or
#'   \code{"slope"} (\code{(Day|ID)}).
#' @param prior_strength Prior preset for all fixed-effect coefficients:
#'   \describe{
#'     \item{\code{"skeptical"}}{(default)
#'       \eqn{b \sim N(0, 0.25)};
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
#' @param plots Logical. Return \pkg{ggplot2} plot objects? Default
#'   \code{TRUE}.
#' @param verbose Logical. Print progress messages? Default \code{FALSE}.
#' @param backend brms backend: \code{"rstan"} (default) or \code{"cmdstanr"}.
#'   See \code{\link{bayesian_tumor_growth}} for details. Applies to both
#'   \code{bayesian_synergy()} and \code{bayesian_synergy_over_time()}.
#'
#' @details
#' The model formula is \code{Volume_transformed ~ Treatment * Day + re},
#' where \code{re} is the random-effects term for animal-level intercepts (or
#' slopes). A single model is fitted on all four groups simultaneously so
#' group-specific uncertainty is fully propagated.
#'
#' @section Assumptions and Limitations:
#' \itemize{
#'   \item \strong{Bliss Independence on TGI:} Bliss Independence was
#'     formulated for probability of cell death, not for proportional
#'     inhibition (TGI). When individual drug effects are large (TGI > 0.5
#'     each), the expected combined effect can approach 1, making it nearly
#'     impossible to detect synergy by the Bliss criterion regardless of
#'     actual biological interaction.
#'   \item \strong{Loewe single-dose approximation:} The Loewe CI formula
#'     \eqn{\min(FE_A + FE_B, 1) / FE_{\text{combo}}} assumes a linear
#'     dose-response relationship. Without EC50 estimates for each drug, this
#'     is a single-dose approximation only. Results should be interpreted
#'     directionally rather than as a precise pharmacological CI.
#'   \item \strong{Posterior independence assumption:} Draws from the TG
#'     model are paired draw-by-draw to propagate uncertainty through the
#'     synergy metrics. This is exact (the model is fitted on all four groups
#'     jointly), but the derived FE values are computed from population-level
#'     predictions that marginalise over animal-level random effects.
#' }
#'
#' At \code{endpoint_day}, \code{brms::posterior_epred(..., re_formula = NA)}
#' yields an \eqn{S \times 4} matrix of population-level predictive draws
#' (S = number of posterior draws). These are back-transformed to the original
#' volume scale and used to compute per-draw TGI and synergy metrics.
#'
#' **Bliss Independence excess:**
#' \deqn{
#'   \Delta_{\text{Bliss}}^{(s)} =
#'   \text{FE}_{\text{combo}}^{(s)} -
#'   \bigl(\text{FE}_A^{(s)} + \text{FE}_B^{(s)} -
#'   \text{FE}_A^{(s)}\,\text{FE}_B^{(s)}\bigr)
#' }
#' Positive values indicate synergy (observed combo effect exceeds
#' Bliss-independence expectation). \eqn{P(\Delta_{\text{Bliss}} > 0)} is
#' reported as the posterior probability of Bliss synergy.
#'
#' **Loewe Combination Index (single-dose approximation):**
#' \deqn{
#'   \text{CI}^{(s)} =
#'   \frac{\min(\text{FE}_A^{(s)} + \text{FE}_B^{(s)},\; 1)}
#'        {\max(\text{FE}_{\text{combo}}^{(s)},\; \epsilon)}
#' }
#' CI < 1 indicates synergy; CI > 1 indicates antagonism.
#' \eqn{P(\text{CI} < 1)} is reported.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{model_type_used}}{Character \code{"bayes_synergy"}.}
#'   \item{\code{model}}{\code{brmsfit} object, or \code{NULL} when
#'     \code{return_model = FALSE}.}
#'   \item{\code{transform_used}}{The transform applied
#'     (\code{"log"}, \code{"sqrt"}, or \code{"none"}).}
#'   \item{\code{tgi_summary}}{Data frame: per-group TGI at the endpoint day —
#'     \code{Group}, \code{TGI_Median}, \code{TGI_Lower}, \code{TGI_Upper}.}
#'   \item{\code{bliss_summary}}{Named list: \code{Excess_Median},
#'     \code{Excess_Lower}, \code{Excess_Upper}, \code{P_Synergy}
#'     (\eqn{P(\Delta>0)}), \code{Expected_FE_Median},
#'     \code{Observed_FE_Median}.}
#'   \item{\code{loewe_summary}}{Named list: \code{CI_Median},
#'     \code{CI_Lower}, \code{CI_Upper}, \code{P_Synergy}
#'     (\eqn{P(\text{CI}<1)}), \code{Interpretation} (posterior-median
#'     classification).}
#'   \item{\code{synergy_table}}{Combined summary data frame with one row per
#'     group plus Bliss-expected and Loewe-expected rows.}
#'   \item{\code{posterior_summary}}{Standard \code{brms} posterior summary.}
#'   \item{\code{mcmc_diagnostics}}{Data frame: per-parameter Rhat,
#'     Bulk_ESS, Tail_ESS, Converged (Rhat ≤ 1.01).}
#'   \item{\code{summary}}{Named list of analysis metadata.}
#'   \item{\code{synergy_plot}}{\pkg{ggplot2} bar/point plot comparing
#'     observed and expected fractional effects with 95 % CrI;
#'     \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{posterior_dist_plot}}{\pkg{ggplot2} density of posterior
#'     Bliss excess and Loewe CI draws; \code{NULL} when
#'     \code{plots = FALSE}.}
#' }
#'
#' @importFrom stats quantile median
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar geom_hline
#'   geom_density geom_vline scale_fill_manual facet_wrap
#'   labs theme theme_classic
#' @export
bayesian_synergy <- function(
  df,
  volume_column               = "Volume",
  time_column                 = "Day",
  treatment_column            = "Treatment",
  id_column                   = "ID",
  cage_column                 = NULL,
  drug_a_name,
  drug_b_name,
  combo_name,
  control_name                = "Control",
  endpoint_day                = NULL,
  transform                   = c("log", "sqrt", "none"),
  random_effects_specification = c("intercept_only", "slope"),
  prior_strength              = c(
    "skeptical", "weakly_informative", "informative", "diffuse", "manual"
  ),
  prior_b                     = NULL,
  prior_intercept             = NULL,
  prior_sd                    = NULL,
  prior_sigma                 = NULL,
  n_chains                    = 4L,
  n_warmup                    = 1000L,
  n_iter                      = 500L,
  seed                        = 42L,
  include_cage_effect         = TRUE,
  return_model                = TRUE,
  plots                       = TRUE,
  verbose                     = FALSE,
  backend                     = c("rstan", "cmdstanr"),
  priors                      = NULL,
  mcmc                        = NULL
) {

  # ── Dependency check ───────────────────────────────────────────────────────
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "'brms' is required for bayesian_synergy(). ",
      "Install it with: install.packages('brms')"
    )
  }

  transform  <- match.arg(transform)
  re_spec    <- match.arg(random_effects_specification)
  prior_str  <- match.arg(prior_strength)
  backend    <- resolve_brms_backend(backend)

  # CODE_REVIEW.md Round 2 D.2 — config overrides.
  .p <- .resolve_priors(priors, prior_str, prior_b, prior_intercept,
                        prior_sd, prior_sigma)
  prior_str       <- .p$strength
  prior_b         <- .p$b
  prior_intercept <- .p$intercept
  prior_sd        <- .p$sd
  prior_sigma     <- .p$sigma
  .m <- .resolve_mcmc(mcmc, n_chains, n_warmup, n_iter, seed, backend)
  n_chains <- .m$chains
  n_warmup <- .m$warmup
  n_iter   <- .m$iter
  seed     <- .m$seed
  backend  <- .m$backend

  # ── Input validation ───────────────────────────────────────────────────────
  required_cols <- c(
    volume_column, time_column, treatment_column, id_column
  )
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  all_groups <- unique(df[[treatment_column]])
  required_groups <- c(control_name, drug_a_name, drug_b_name, combo_name)
  missing_groups <- required_groups[!required_groups %in% all_groups]
  if (length(missing_groups) > 0L) {
    stop(
      "The following groups are not found in '", treatment_column, "': ",
      paste(missing_groups, collapse = ", ")
    )
  }

  # ── Subset to the four analysis groups ────────────────────────────────────
  analysis_df <- df[df[[treatment_column]] %in% required_groups, ]
  analysis_df[[treatment_column]] <- factor(
    analysis_df[[treatment_column]],
    levels = c(control_name, drug_a_name, drug_b_name, combo_name)
  )

  # ── Endpoint day ──────────────────────────────────────────────────────────
  all_days <- sort(unique(as.numeric(analysis_df[[time_column]])))
  ep_day   <- if (!is.null(endpoint_day)) {
    as.numeric(endpoint_day)
  } else {
    all_days[length(all_days)]
  }
  if (!ep_day %in% all_days) {
    ep_day <- all_days[which.min(abs(all_days - ep_day))]
    warning("Requested endpoint_day not in data; using closest: ", ep_day)
  }

  if (isTRUE(verbose)) {
    message("bayesian_synergy(): endpoint = day ", ep_day,
            "; transform = ", transform)
  }

  # ── Fit model via shared helper (CODE_REVIEW.md G.6) ──────────────────────
  fit <- bs_fit_synergy_model(
    analysis_df        = analysis_df,
    volume_column      = volume_column,
    treatment_column   = treatment_column,
    time_column        = time_column,
    id_column          = id_column,
    cage_column        = cage_column,
    transform          = transform,
    include_cage_effect = include_cage_effect,
    re_spec            = re_spec,
    prior_str          = prior_str,
    manual_priors      = list(prior_b = prior_b, prior_intercept = prior_intercept,
                              prior_sd = prior_sd, prior_sigma = prior_sigma),
    n_chains = n_chains, n_warmup = n_warmup, n_iter = n_iter,
    seed = seed, verbose = verbose, backend = backend
  )
  model             <- fit$model
  analysis_df       <- fit$analysis_df
  has_cage          <- fit$has_cage
  posterior_summary <- fit$posterior_summary
  mcmc_diagnostics  <- fit$mcmc_diagnostics
  nuts_diagnostics  <- fit$nuts_diagnostics
  loo_diagnostics   <- fit$loo_diagnostics
  bayes_r2          <- fit$bayes_R2
  ppc_coverage      <- fit$ppc_coverage

  # ── Posterior predictive draws at endpoint day ─────────────────────────────
  groups <- c(control_name, drug_a_name, drug_b_name, combo_name)
  nd     <- data.frame(
    g = groups, d = ep_day, stringsAsFactors = FALSE
  )
  names(nd) <- c(treatment_column, time_column)

  epred <- tryCatch(
    brms::posterior_epred(
      model, newdata = nd,
      re_formula = NA, allow_new_levels = TRUE
    ),
    error = function(e) {
      stop("posterior_epred() failed: ", conditionMessage(e))
    }
  )
  colnames(epred) <- groups

  vol_mat <- bayes_backtransform(epred, transform)

  # ── Draw-wise TGI and FE ──────────────────────────────────────────────────
  v_ctrl  <- vol_mat[, control_name]
  fe_mat  <- 1 - sweep(vol_mat, 1L, v_ctrl, FUN = "/")
  colnames(fe_mat) <- groups

  fe_a     <- fe_mat[, drug_a_name]
  fe_b     <- fe_mat[, drug_b_name]
  fe_combo <- fe_mat[, combo_name]

  # ── TGI summary ───────────────────────────────────────────────────────────
  tgi_summary <- do.call(rbind, lapply(groups, function(g) {
    q <- stats::quantile(fe_mat[, g], c(0.025, 0.5, 0.975))
    data.frame(
      Group      = g,
      TGI_Median = round(q["50%"],   3),
      TGI_Lower  = round(q["2.5%"],  3),
      TGI_Upper  = round(q["97.5%"], 3),
      stringsAsFactors = FALSE
    )
  }))
  rownames(tgi_summary) <- NULL

  # ── Draw-wise Bliss Independence ──────────────────────────────────────────
  bliss_expected  <- synergy_bliss_expected(fe_a, fe_b)
  bliss_excess    <- fe_combo - bliss_expected

  bliss_q     <- stats::quantile(bliss_excess,   c(0.025, 0.5, 0.975))
  bliss_fe_q  <- stats::quantile(bliss_expected, c(0.025, 0.5, 0.975))
  obs_fe_med  <- stats::median(fe_combo)

  bliss_summary <- list(
    Excess_Median      = round(bliss_q["50%"],         3),
    Excess_Lower       = round(bliss_q["2.5%"],        3),
    Excess_Upper       = round(bliss_q["97.5%"],       3),
    P_Synergy          = round(mean(bliss_excess > 0), 3),
    Expected_FE_Median = round(bliss_fe_q["50%"],      3),
    Observed_FE_Median = round(obs_fe_med,             3)
  )

  # ── Draw-wise Loewe CI ────────────────────────────────────────────────────
  # Floor is relative to the maximum observed FE to avoid huge CIs when
  # the combo has near-zero effect. A note is added when the floor is applied.
  fe_max          <- max(fe_combo, na.rm = TRUE)
  fe_floor        <- max(fe_max * 1e-4, 1e-4)
  loewe_res       <- synergy_loewe_ci(fe_a, fe_b, fe_combo, fe_floor = fe_floor)
  loewe_num       <- loewe_res$loewe_num
  loewe_ci        <- loewe_res$ci
  floor_applied   <- any(loewe_res$floor_applied)

  loewe_q    <- stats::quantile(loewe_ci, c(0.025, 0.5, 0.975))
  ci_med     <- loewe_q["50%"]

  loewe_interp <- if (ci_med < 0.85) {
    "Synergistic (CI < 0.85)"
  } else if (ci_med <= 1.15) {
    "Additive (0.85 ≤ CI ≤ 1.15)"
  } else {
    "Antagonistic (CI > 1.15)"
  }

  loewe_summary <- list(
    CI_Median      = round(ci_med,               3),
    CI_Lower       = round(loewe_q["2.5%"],      3),
    CI_Upper       = round(loewe_q["97.5%"],     3),
    P_Synergy      = round(mean(loewe_ci < 1),   3),
    Interpretation = loewe_interp,
    Floor_Applied  = floor_applied
  )

  # ── synergy_table ─────────────────────────────────────────────────────────
  vol_ctrl_med  <- stats::median(vol_mat[, control_name])

  synergy_table <- do.call(rbind, lapply(groups, function(g) {
    q   <- stats::quantile(fe_mat[, g], c(0.025, 0.5, 0.975))
    qv  <- stats::quantile(vol_mat[, g], c(0.025, 0.5, 0.975))
    data.frame(
      Group          = g,
      Mean_Volume    = round(qv["50%"],  1),
      TGI_Percent    = round(q["50%"] * 100, 1),
      TGI_Lower      = round(q["2.5%"] * 100, 1),
      TGI_Upper      = round(q["97.5%"] * 100, 1),
      Type           = "Observed",
      stringsAsFactors = FALSE
    )
  }))

  bliss_fe_med <- bliss_fe_q["50%"]
  loewe_fe_med <- stats::median(loewe_num)

  extra_rows <- data.frame(
    Group       = c("Bliss Expected", "Loewe Expected"),
    Mean_Volume = round(vol_ctrl_med * c(
      1 - bliss_fe_med, 1 - loewe_fe_med
    ), 1),
    TGI_Percent = round(c(bliss_fe_med, loewe_fe_med) * 100, 1),
    TGI_Lower   = NA_real_,
    TGI_Upper   = NA_real_,
    Type        = "Expected",
    stringsAsFactors = FALSE
  )
  synergy_table <- rbind(synergy_table, extra_rows)
  rownames(synergy_table) <- NULL

  # ── Plots ──────────────────────────────────────────────────────────────────
  synergy_plot     <- NULL
  post_dist_plot   <- NULL

  if (isTRUE(plots)) {

    # Bar chart: observed TGI per group + Bliss/Loewe expected lines
    obs_rows  <- synergy_table[synergy_table$Type == "Observed", ]
    exp_rows  <- synergy_table[synergy_table$Type == "Expected", ]
    obs_rows$Group <- factor(obs_rows$Group, levels = groups)

    synergy_plot <- ggplot2::ggplot(
      obs_rows,
      ggplot2::aes(
        x    = .data[["Group"]],
        y    = .data[["TGI_Percent"]],
        fill = .data[["Group"]]
      )
    ) +
      ggplot2::geom_col(alpha = 0.75, width = 0.6) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data[["TGI_Lower"]],
          ymax = .data[["TGI_Upper"]]
        ),
        width = 0.2, colour = "grey30"
      ) +
      ggplot2::geom_hline(
        yintercept = exp_rows$TGI_Percent[
          exp_rows$Group == "Bliss Expected"
        ],
        linetype = "dashed", colour = "steelblue", linewidth = 0.8
      ) +
      ggplot2::geom_hline(
        yintercept = exp_rows$TGI_Percent[
          exp_rows$Group == "Loewe Expected"
        ],
        linetype = "dotted", colour = "tomato", linewidth = 0.8
      ) +
      ggplot2::labs(
        title    = "Bayesian Drug Combination Synergy",
        subtitle = paste0(
          "Bars = posterior median TGI ± 95 % CrI; ",
          "dashed = Bliss expected; dotted = Loewe expected"
        ),
        x    = NULL,
        y    = "TGI (%)",
        fill = NULL
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(legend.position = "none")

    # Density overlay: Bliss excess and Loewe CI draws
    n_draws  <- length(bliss_excess)
    dens_df  <- data.frame(
      Value  = c(bliss_excess, loewe_ci),
      Metric = rep(
        c("Bliss Excess (>0 = synergy)",
          "Loewe CI (<1 = synergy)"),
        each = n_draws
      ),
      stringsAsFactors = FALSE
    )
    vlines_df <- data.frame(
      Metric    = c(
        "Bliss Excess (>0 = synergy)",
        "Loewe CI (<1 = synergy)"
      ),
      xintercept = c(0, 1),
      stringsAsFactors = FALSE
    )

    post_dist_plot <- ggplot2::ggplot(
      dens_df,
      ggplot2::aes(x = .data[["Value"]])
    ) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.4) +
      ggplot2::geom_vline(
        data    = vlines_df,
        mapping = ggplot2::aes(xintercept = .data[["xintercept"]]),
        linetype = "dashed", colour = "tomato", linewidth = 0.7
      ) +
      ggplot2::facet_wrap(~ Metric, scales = "free") +
      ggplot2::labs(
        title    = "Posterior Distributions of Synergy Metrics",
        subtitle = paste0(
          "Bliss: P(synergy) = ", bliss_summary$P_Synergy,
          "; Loewe: P(CI<1) = ", loewe_summary$P_Synergy
        ),
        x = "Value",
        y = "Density"
      ) +
      ggplot2::theme_classic(base_size = 14)
  }

  # ── Standard Bayesian diagnostic plots (G.4: match what TG/BW/Survival
  # / DR functions already return). All three are skipped silently if the
  # underlying brms helpers aren't available — keeps the function usable
  # in lighter environments.
  pp_check_plot         <- tryCatch(
    if (requireNamespace("brms", quietly = TRUE)) brms::pp_check(model) else NULL,
    error = function(e) NULL
  )
  prior_posterior_plot  <- tryCatch(
    bayes_prior_posterior_plot(model, treatment_column),
    error = function(e) NULL
  )
  mcmc_trace_plot       <- tryCatch(
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      bayesplot::mcmc_trace(
        brms::as_draws_df(model),
        regex_pars = paste0("^b_", gsub("([.^$*+?()\\[\\]{}|])", "\\\\\\1",
                                         treatment_column))
      )
    } else NULL,
    error = function(e) NULL
  )

  # ── Summary metadata ───────────────────────────────────────────────────────
  analysis_summary <- list(
    analysis_type = paste0(
      "Bayesian Drug Combination Synergy — draw-wise Bliss and Loewe ",
      "metrics (brms LME)"
    ),
    data_description = list(
      control_group = control_name,
      drug_a        = drug_a_name,
      drug_b        = drug_b_name,
      combo         = combo_name,
      endpoint_day  = ep_day,
      n_animals     = nrow(
        analysis_df[analysis_df[[time_column]] == all_days[1L], ]
      )
    ),
    model_specification = list(
      transform       = transform,
      prior_strength  = prior_str,
      n_chains        = n_chains,
      n_warmup        = n_warmup,
      n_iter          = n_iter,
      seed            = seed,
      random_effects  = re_term
    ),
    synergy_interpretation = list(
      bliss = bliss_summary,
      loewe = loewe_summary
    )
  )

  # ── Return ─────────────────────────────────────────────────────────────────
  out <- list(
    model_type_used     = "bayes_synergy",
    model               = if (isTRUE(return_model)) model else NULL,
    transform_used      = transform,
    tgi_summary         = tgi_summary,
    bliss_summary       = bliss_summary,
    loewe_summary       = loewe_summary,
    synergy_table       = synergy_table,
    posterior_summary   = posterior_summary,
    mcmc_diagnostics    = mcmc_diagnostics,
    nuts_diagnostics    = nuts_diagnostics,
    loo_diagnostics     = loo_diagnostics,
    bayes_R2            = bayes_r2,
    ppc_coverage        = ppc_coverage,
    summary               = analysis_summary,
    synergy_plot          = synergy_plot,
    posterior_dist_plot   = post_dist_plot,
    pp_check_plot         = pp_check_plot,
    prior_posterior_plot  = prior_posterior_plot,
    mcmc_trace_plot       = mcmc_trace_plot
  )
  out
}


#' Bayesian Drug Combination Synergy Over Time
#'
#' Fits the same Bayesian linear mixed-effects model as
#' \code{\link{bayesian_synergy}} (a single model on all four treatment groups
#' with a \code{Treatment × Day} interaction term), then evaluates draw-wise
#' Bliss Independence excess and Loewe Combination Index at every observed
#' study day via a single \code{brms::posterior_epred()} call on the full
#' day-by-group grid. Returns per-day posterior summaries that show how synergy
#' evolves across the experiment.
#'
#' @inheritParams bayesian_synergy
#' @param study_days Numeric vector of study days at which to evaluate synergy.
#'   \code{NULL} (default) uses all unique days observed in \code{df}.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{model_type_used}}{Character \code{"bayes_synergy_ot"}.}
#'   \item{\code{model}}{\code{brmsfit} object, or \code{NULL} when
#'     \code{return_model = FALSE}.}
#'   \item{\code{transform_used}}{The transform applied.}
#'   \item{\code{synergy_by_day}}{Data frame with one row per study day:
#'     \code{Day}, \code{Bliss_Median}, \code{Bliss_Lower},
#'     \code{Bliss_Upper}, \code{P_Bliss_Synergy},
#'     \code{Loewe_Median}, \code{Loewe_Lower}, \code{Loewe_Upper},
#'     \code{P_Loewe_Synergy}, \code{Loewe_Floor_Applied}.}
#'   \item{\code{tgi_by_day}}{Data frame with one row per (Group, Day):
#'     \code{Day}, \code{Group}, \code{TGI_Median}, \code{TGI_Lower},
#'     \code{TGI_Upper}.}
#'   \item{\code{peak_bliss_day}}{List: \code{Day} and \code{Bliss_Median}
#'     at the study day showing highest posterior-median Bliss excess.}
#'   \item{\code{peak_loewe_day}}{List: \code{Day} and \code{Loewe_Median}
#'     (lowest posterior-median Loewe CI, i.e. strongest synergy).}
#'   \item{\code{posterior_summary}}{Standard \code{brms} posterior summary.}
#'   \item{\code{mcmc_diagnostics}}{Data frame: per-parameter Rhat,
#'     Bulk_ESS, Tail_ESS, Converged.}
#'   \item{\code{summary}}{Named list of analysis metadata.}
#'   \item{\code{synergy_time_plot}}{\pkg{ggplot2} faceted line plot of
#'     posterior-median Bliss excess and Loewe CI over time with 95 % CrI
#'     ribbons; \code{NULL} when \code{plots = FALSE}.}
#' }
#'
#' @seealso \code{\link{bayesian_synergy}},
#'   \code{\link{analyze_drug_synergy_over_time}}
#'
#' @importFrom stats quantile median
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline
#'   facet_wrap labs theme theme_classic
#' @export
bayesian_synergy_over_time <- function(
  df,
  volume_column                = "Volume",
  time_column                  = "Day",
  treatment_column             = "Treatment",
  id_column                    = "ID",
  cage_column                  = NULL,
  drug_a_name,
  drug_b_name,
  combo_name,
  control_name                 = "Control",
  study_days                   = NULL,
  transform                    = c("log", "sqrt", "none"),
  random_effects_specification = c("intercept_only", "slope"),
  prior_strength               = c(
    "skeptical", "weakly_informative", "informative", "diffuse", "manual"
  ),
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
  verbose                      = FALSE,
  backend                      = c("rstan", "cmdstanr")
) {

  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "'brms' is required for bayesian_synergy_over_time(). ",
      "Install it with: install.packages('brms')"
    )
  }

  transform  <- match.arg(transform)
  re_spec    <- match.arg(random_effects_specification)
  prior_str  <- match.arg(prior_strength)
  backend    <- resolve_brms_backend(backend)

  # ── Input validation ───────────────────────────────────────────────────────
  required_cols <- c(volume_column, time_column, treatment_column, id_column)
  missing_cols  <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  all_groups      <- unique(df[[treatment_column]])
  required_groups <- c(control_name, drug_a_name, drug_b_name, combo_name)
  missing_groups  <- required_groups[!required_groups %in% all_groups]
  if (length(missing_groups) > 0L) {
    stop(
      "The following groups are not found in '", treatment_column, "': ",
      paste(missing_groups, collapse = ", ")
    )
  }

  # ── Subset to the four analysis groups ────────────────────────────────────
  analysis_df <- df[df[[treatment_column]] %in% required_groups, ]
  analysis_df[[treatment_column]] <- factor(
    analysis_df[[treatment_column]],
    levels = c(control_name, drug_a_name, drug_b_name, combo_name)
  )

  all_days   <- sort(unique(as.numeric(analysis_df[[time_column]])))
  if (is.null(study_days)) {
    study_days <- all_days
  } else {
    study_days <- sort(intersect(as.numeric(study_days), all_days))
  }
  if (length(study_days) == 0L) {
    stop("'study_days' contains no days present in the data.")
  }

  if (isTRUE(verbose)) {
    message("bayesian_synergy_over_time(): ", length(study_days),
            " days; transform = ", transform)
  }

  # ── Fit model via shared helper (CODE_REVIEW.md G.6) ──────────────────────
  fit <- bs_fit_synergy_model(
    analysis_df        = analysis_df,
    volume_column      = volume_column,
    treatment_column   = treatment_column,
    time_column        = time_column,
    id_column          = id_column,
    cage_column        = cage_column,
    transform          = transform,
    include_cage_effect = include_cage_effect,
    re_spec            = re_spec,
    prior_str          = prior_str,
    manual_priors      = list(prior_b = prior_b, prior_intercept = prior_intercept,
                              prior_sd = prior_sd, prior_sigma = prior_sigma),
    n_chains = n_chains, n_warmup = n_warmup, n_iter = n_iter,
    seed = seed, verbose = verbose, backend = backend
  )
  model             <- fit$model
  analysis_df       <- fit$analysis_df
  has_cage          <- fit$has_cage
  posterior_summary <- fit$posterior_summary
  mcmc_diagnostics  <- fit$mcmc_diagnostics
  nuts_diagnostics  <- fit$nuts_diagnostics
  loo_diagnostics   <- fit$loo_diagnostics
  bayes_r2          <- fit$bayes_R2
  ppc_coverage      <- fit$ppc_coverage

  # ── Posterior predictive draws for all (group, day) combinations ──────────
  groups  <- c(control_name, drug_a_name, drug_b_name, combo_name)
  n_grps  <- length(groups)
  nd_full <- expand.grid(
    setNames(list(groups, study_days), c(treatment_column, time_column)),
    stringsAsFactors = FALSE
  )

  if (isTRUE(verbose)) {
    message("Computing posterior_epred for ", nrow(nd_full), " conditions …")
  }

  epred_full <- tryCatch(
    brms::posterior_epred(
      model, newdata = nd_full,
      re_formula = NA, allow_new_levels = TRUE
    ),
    error = function(e) {
      stop("posterior_epred() failed: ", conditionMessage(e))
    }
  )

  vol_full <- bayes_backtransform(epred_full, transform)

  # ── Per-day synergy summaries ─────────────────────────────────────────────
  # expand.grid cycles first arg (groups) fastest, so columns in vol_full are:
  #   [group1@day1, group2@day1, ..., group1@day2, ...]
  synergy_rows <- vector("list", length(study_days))
  tgi_rows     <- vector("list", length(study_days))

  for (di in seq_along(study_days)) {
    day     <- study_days[[di]]
    col_idx <- seq_len(n_grps) + (di - 1L) * n_grps
    vm      <- vol_full[, col_idx, drop = FALSE]
    colnames(vm) <- groups

    v_ctrl  <- vm[, control_name]
    fe_mat  <- 1 - sweep(vm, 1L, v_ctrl, FUN = "/")
    colnames(fe_mat) <- groups

    fe_a     <- fe_mat[, drug_a_name]
    fe_b     <- fe_mat[, drug_b_name]
    fe_combo <- fe_mat[, combo_name]

    # Bliss
    bliss_expected <- synergy_bliss_expected(fe_a, fe_b)
    bliss_excess   <- fe_combo - bliss_expected
    bq             <- stats::quantile(bliss_excess, c(0.025, 0.5, 0.975))

    # Loewe
    fe_max        <- max(fe_combo, na.rm = TRUE)
    fe_floor      <- max(fe_max * 1e-4, 1e-4)
    loewe_res     <- synergy_loewe_ci(fe_a, fe_b, fe_combo,
                                       fe_floor = fe_floor)
    loewe_ci      <- loewe_res$ci
    floor_applied <- any(loewe_res$floor_applied)
    lq            <- stats::quantile(loewe_ci, c(0.025, 0.5, 0.975))

    synergy_rows[[di]] <- data.frame(
      Day                = day,
      Bliss_Median       = round(bq["50%"],             3),
      Bliss_Lower        = round(bq["2.5%"],            3),
      Bliss_Upper        = round(bq["97.5%"],           3),
      P_Bliss_Synergy    = round(mean(bliss_excess > 0), 3),
      Loewe_Median       = round(lq["50%"],             3),
      Loewe_Lower        = round(lq["2.5%"],            3),
      Loewe_Upper        = round(lq["97.5%"],           3),
      P_Loewe_Synergy    = round(mean(loewe_ci < 1),    3),
      Loewe_Floor_Applied = floor_applied,
      stringsAsFactors   = FALSE
    )

    tgi_rows[[di]] <- do.call(rbind, lapply(groups, function(g) {
      q <- stats::quantile(fe_mat[, g], c(0.025, 0.5, 0.975))
      data.frame(
        Day            = day,
        Group          = g,
        TGI_Median     = round(q["50%"],   3),
        TGI_Lower      = round(q["2.5%"],  3),
        TGI_Upper      = round(q["97.5%"], 3),
        stringsAsFactors = FALSE
      )
    }))
  }

  synergy_by_day <- do.call(rbind, synergy_rows)
  tgi_by_day     <- do.call(rbind, tgi_rows)
  rownames(synergy_by_day) <- NULL
  rownames(tgi_by_day)     <- NULL

  # ── Peak synergy ──────────────────────────────────────────────────────────
  peak_b_idx <- which.max(synergy_by_day$Bliss_Median)
  peak_l_idx <- which.min(synergy_by_day$Loewe_Median)

  peak_bliss_day <- list(
    Day          = synergy_by_day$Day[peak_b_idx],
    Bliss_Median = synergy_by_day$Bliss_Median[peak_b_idx]
  )
  peak_loewe_day <- list(
    Day          = synergy_by_day$Day[peak_l_idx],
    Loewe_Median = synergy_by_day$Loewe_Median[peak_l_idx]
  )

  # ── Plot ──────────────────────────────────────────────────────────────────
  synergy_time_plot <- NULL

  if (isTRUE(plots)) {
    n_days <- nrow(synergy_by_day)

    plot_df <- data.frame(
      Day    = rep(synergy_by_day$Day, 2L),
      Metric = rep(
        c("Bliss Excess (>0 = synergy)", "Loewe CI (<1 = synergy)"),
        each = n_days
      ),
      Median = c(synergy_by_day$Bliss_Median, synergy_by_day$Loewe_Median),
      Lower  = c(synergy_by_day$Bliss_Lower,  synergy_by_day$Loewe_Lower),
      Upper  = c(synergy_by_day$Bliss_Upper,  synergy_by_day$Loewe_Upper),
      Ref    = c(
        rep(0, n_days),
        rep(1, n_days)
      ),
      stringsAsFactors = FALSE
    )

    synergy_time_plot <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[["Day"]])
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[["Lower"]],
          ymax = .data[["Upper"]]
        ),
        fill = "steelblue", alpha = 0.25
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data[["Median"]]),
        colour = "steelblue", linewidth = 0.9
      ) +
      ggplot2::geom_hline(
        ggplot2::aes(yintercept = .data[["Ref"]]),
        linetype = "dashed", colour = "tomato", linewidth = 0.7
      ) +
      ggplot2::facet_wrap(
        ~ Metric, scales = "free_y", ncol = 1L
      ) +
      ggplot2::labs(
        title    = "Bayesian Drug Combination Synergy Over Time",
        subtitle = paste0(
          drug_a_name, " + ", drug_b_name, " vs. ", combo_name,
          "; ribbon = 95 % CrI"
        ),
        x = "Study Day",
        y = "Posterior Median"
      ) +
      ggplot2::theme_classic(base_size = 14)
  }

  # ── Summary metadata ───────────────────────────────────────────────────────
  analysis_summary <- list(
    analysis_type = paste0(
      "Bayesian Drug Combination Synergy Over Time — draw-wise Bliss and ",
      "Loewe at each study day (brms LME)"
    ),
    data_description = list(
      control_group = control_name,
      drug_a        = drug_a_name,
      drug_b        = drug_b_name,
      combo         = combo_name,
      study_days    = study_days
    ),
    model_specification = list(
      transform      = transform,
      prior_strength = prior_str,
      n_chains       = n_chains,
      n_iter         = n_iter,
      seed           = seed,
      random_effects = re_term
    ),
    peak_synergy = list(
      bliss = peak_bliss_day,
      loewe = peak_loewe_day
    )
  )

  list(
    model_type_used   = "bayes_synergy_ot",
    model             = if (isTRUE(return_model)) model else NULL,
    transform_used    = transform,
    synergy_by_day    = synergy_by_day,
    tgi_by_day        = tgi_by_day,
    peak_bliss_day    = peak_bliss_day,
    peak_loewe_day    = peak_loewe_day,
    posterior_summary = posterior_summary,
    mcmc_diagnostics  = mcmc_diagnostics,
    nuts_diagnostics  = nuts_diagnostics,
    loo_diagnostics   = loo_diagnostics,
    bayes_R2          = bayes_r2,
    ppc_coverage      = ppc_coverage,
    summary           = analysis_summary,
    synergy_time_plot = synergy_time_plot
  )
}
