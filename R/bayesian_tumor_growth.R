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
#'   When supplied and \code{include_cage_effect = TRUE}, a cage-level random
#'   intercept is added to the model (\code{(1|cage/ID)} for intercept-only
#'   random effects, \code{(Day|ID) + (1|cage)} for random slopes).
#'   The column is also used for composite mouse-key construction.
#' @param include_cage_effect Logical. When \code{TRUE} (default) and
#'   \code{cage_column} is supplied, a cage random intercept is included in the
#'   random-effects structure. Set \code{FALSE} to use cage only for key
#'   construction without modelling it.
#' @param dose_column Reserved for future use; currently ignored.
#' @param transform Volume transformation applied before modelling:
#'   \code{"log"} (default, recommended for exponential growth),
#'   \code{"sqrt"}, or \code{"none"}.
#' @param random_effects_specification Random effects structure:
#'   \code{"intercept_only"} (default, \code{(1 | ID)}) or
#'   \code{"slope"} (\code{(Day | ID)}, adds per-animal random slopes).
#' @param reference_group Treatment group used as the reference (Intercept)
#'   level. Auto-detected if \code{NULL}: checks for common control-group
#'   names (\code{Control}, \code{Vehicle}, \code{control}, \code{vehicle},
#'   \code{CTRL}, \code{ctrl}) before falling back to the first level
#'   alphabetically.
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
#' @param n_warmup Warm-up (burn-in) iterations per chain. Default \code{1000}.
#' @param n_iter Post-warmup draws per chain. Default \code{500}.
#' @param seed Integer random seed for reproducibility. Default \code{42}.
#' @param return_model Logical. Return the \code{brmsfit} object? Default
#'   \code{TRUE}. Set \code{FALSE} to reduce memory when only summaries are
#'   needed.
#' @param plots Logical. Pre-compute and return \pkg{ggplot2} / \pkg{bayesplot}
#'   objects? Default \code{TRUE}. Set \code{FALSE} when plots are generated
#'   reactively from the stored model object (e.g. in Shiny).
#' @param verbose Logical. Show Stan compilation and sampling messages? Default
#'   \code{FALSE}.
#' @param backend brms backend: \code{"rstan"} (default) or \code{"cmdstanr"}.
#'   \code{"cmdstanr"} is 3-10x faster to compile and produces better
#'   diagnostic reporting but requires the \pkg{cmdstanr} package plus a
#'   working CmdStan toolchain (install via
#'   \code{cmdstanr::install_cmdstan()}). Recommended for production use
#'   on the VPS; the default keeps \pkg{rstan} so the function works
#'   out-of-the-box.
#' @param priors Optional \code{\link{tg_priors}()} configuration object.
#'   When supplied, its fields override the individual \code{prior_*}
#'   arguments. CODE_REVIEW.md Round 2 D.2 — lets callers bundle the
#'   five prior arguments into a single value.
#' @param mcmc Optional \code{\link{tg_mcmc}()} configuration object. When
#'   supplied, its fields override the individual \code{n_chains},
#'   \code{n_warmup}, \code{n_iter}, \code{seed}, and \code{backend}
#'   arguments.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{model}}{\code{brmsfit} object, or \code{NULL} when
#'     \code{return_model = FALSE}.}
#'   \item{\code{model_type_used}}{Character \code{"bayes_tg"}.}
#'   \item{\code{transform_used}}{Transformation applied.}
#'   \item{\code{summary}}{Named list of analysis metadata (analysis_type,
#'     data_description, model_specification, methods).}
#'   \item{\code{posterior_summary}}{Data frame of fixed-effect posterior
#'     medians, 2.5 %–97.5 % CrI, Rhat, Bulk_ESS, Tail_ESS for all
#'     parameters.}
#'   \item{\code{treatment_effects}}{Data frame with columns
#'     \code{Group}, \code{Adjusted_Mean}, \code{SE}, \code{DF},
#'     \code{Lower_CrI}, \code{Upper_CrI}, \code{Note} — compatible with the
#'     \code{lme4} path. Values are posterior medians and 95 % HPD CrI of
#'     group-level EMMs at the mean study day.}
#'   \item{\code{pairwise_comparisons}}{Data frame of posterior contrasts vs
#'     the reference group (columns: \code{contrast}, \code{estimate},
#'     \code{Lower_CrI}, \code{Upper_CrI}).}
#'   \item{\code{mcmc_diagnostics}}{Data frame with \code{Parameter},
#'     \code{Rhat}, \code{Bulk_ESS}, \code{Tail_ESS}, \code{Converged}
#'     per fixed-effect parameter. Rhat > 1.01 is flagged as not converged.}
#'   \item{\code{pp_check_plot}}{Posterior predictive density overlay ggplot,
#'     or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{posterior_dist_plot}}{Posterior density areas for treatment
#'     parameters (bayesplot), or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{prior_posterior_plot}}{Density overlay of prior (grey) and
#'     posterior (blue) for each treatment-effect coefficient. A vertical
#'     dashed line marks zero (no effect). Useful for assessing how much the
#'     data updated the prior. \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{credible_intervals_plot}}{Forest plot of group EMMs with
#'     95 % CrI, or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{mcmc_trace_plot}}{MCMC trace plot for treatment parameters,
#'     or \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{residuals_plot}}{Posterior mean residuals vs. study day,
#'     faceted by treatment group with a loess smoother. Systematic curvature
#'     indicates that the log-linear time assumption is violated and a GAM or
#'     nonlinear model should be considered. \code{NULL} when
#'     \code{plots = FALSE}.}
#'   \item{\code{growth_rates}}{Per-animal exponential growth rates. When the
#'     LMM was fit with \code{random_effects_specification = "slope"}, these
#'     are derived directly from the brms posterior — each row carries
#'     \code{growth_rate} (posterior median) plus \code{growth_rate_lower}
#'     and \code{growth_rate_upper} (2.5% / 97.5% credible intervals). For
#'     intercept-only and GAM paths, falls back to log-linear OLS per
#'     animal and the credible-interval columns reflect OLS confidence
#'     intervals (Round 2 E.3).}
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
#' @section Separating analysis from visualization (CODE_REVIEW.md D.3):
#' Every \code{bayesian_*} function in this package returns the underlying
#' analysis data as data frames (\code{treatment_effects},
#' \code{posterior_summary}, \code{pairwise_comparisons},
#' \code{mcmc_diagnostics}, \code{growth_rates}, …) and — separately — a
#' set of pre-rendered \pkg{ggplot2} plot objects (\code{pp_check_plot},
#' \code{posterior_dist_plot}, \code{credible_intervals_plot},
#' \code{mcmc_trace_plot}, …). Pass \code{plots = FALSE} to skip plot
#' generation and receive only the data; the plot fields are then
#' \code{NULL} and the returned list is suitable for headless / CI
#' pipelines or for callers that prefer to render plots from the data
#' using the package's \code{plot_*()} helpers (\code{\link{plot_auc}},
#' \code{\link{plot_growth_rate}}, etc.) or their own renderers. The Shiny
#' dashboard uses this pattern: it passes \code{plots = FALSE} and rebuilds
#' every plot reactively from the data frames so user UI changes (colour,
#' group order, …) take effect without re-running the analysis.
#'
#' @references
#' Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models
#' using Stan. \emph{Journal of Statistical Software}, 80(1), 1–28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @importFrom stats as.formula relevel
#' @importFrom ggplot2 ggplot aes geom_density geom_point geom_pointrange
#'   geom_smooth geom_vline geom_hline scale_colour_manual scale_fill_manual
#'   facet_wrap labs theme theme_classic
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
  model_type                   = c("lmm", "gam"),
  random_effects_specification = c("intercept_only", "slope"),
  reference_group              = NULL,
  prior_strength               = c("skeptical", "weakly_informative", "informative", "diffuse", "manual"),
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
  necrotic_column              = NULL,
  necrotic_handling            = c("exclude", "covariate", "none"),
  backend                      = c("rstan", "cmdstanr"),
  priors                       = NULL,
  mcmc                         = NULL
) {

  # ── Dependency checks ──────────────────────────────────────────────────────

  transform                    <- match.arg(transform)
  model_type                   <- match.arg(model_type)
  random_effects_specification <- match.arg(random_effects_specification)
  prior_strength               <- match.arg(prior_strength)
  necrotic_handling            <- match.arg(necrotic_handling)
  backend                      <- resolve_brms_backend(backend)

  # CODE_REVIEW.md Round 2 D.2 — when callers supply `priors = tg_priors()`
  # or `mcmc = tg_mcmc()`, their fields override the legacy individual
  # arguments so the signature stays usable both ways.
  .p <- .resolve_priors(priors, prior_strength, prior_b, prior_intercept,
                        prior_sd, prior_sigma)
  prior_strength  <- .p$strength
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

  # ── Column validation ──────────────────────────────────────────────────────
  required_cols <- c(time_column, volume_column, treatment_column, id_column)
  missing_cols  <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Cage placeholder — mirrors lme4 path so make_mouse_key() always works
  cage_setup   <- setup_cage_column(df, cage_column)
  df           <- cage_setup$df
  cage_column  <- cage_setup$cage_column
  no_cage_mode <- cage_setup$no_cage_mode

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
    vol           <- analysis_df[[volume_column]]
    positive_vals <- vol[is.finite(vol) & vol > 0]
    if (length(positive_vals) == 0L) {
      stop(
        "No positive values in '", volume_column,
        "' after filtering. Cannot apply log transform."
      )
    }
    vol[vol <= 0] <- min(positive_vals, na.rm = TRUE) / 2
    analysis_df[[volume_column]] <- log(vol)
  } else if (transform == "sqrt") {
    analysis_df[[volume_column]] <- sqrt(analysis_df[[volume_column]])
  }

  # ── Formula ────────────────────────────────────────────────────────────────
  use_cage_re <- isTRUE(include_cage_effect) && !no_cage_mode

  re_term <- if (use_cage_re) {
    if (random_effects_specification == "slope") {
      # Random slopes per animal + random intercept per cage (crossed)
      paste0("(", time_column, " | ", id_column, ") + (1 | ", cage_column, ")")
    } else {
      # Animals nested within cages: (1|cage) + (1|cage:ID)
      paste0("(1 | ", cage_column, "/", id_column, ")")
    }
  } else {
    if (random_effects_specification == "slope") {
      paste0("(", time_column, " | ", id_column, ")")
    } else {
      paste0("(1 | ", id_column, ")")
    }
  }

  # Linear interaction (LMM) vs group-specific smoother (GAM via brms s()).
  # The smoother basis dimension is auto-chosen from observed time points.
  if (model_type == "gam") {
    n_days     <- length(unique(analysis_df[[time_column]]))
    k_val      <- max(3L, min(10L, n_days - 1L))
    fixed_part <- paste0(
      volume_column, " ~ ", treatment_column,
      " + s(", time_column, ", by = ", treatment_column,
      ", k = ", k_val, ")"
    )
  } else {
    fixed_part <- paste(volume_column, "~", treatment_column, "*", time_column)
  }
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
    # CODE_REVIEW.md R3.8 / G.5 — data-scaled Intercept and per-coefficient
    # b priors. See R/utils_bayes.R for why the previous fixed
    # normal(0, b_sd * 2.5) Intercept prior sat ~9 SDs from the data on a
    # mm3-scale study, and why one blanket class = "b" prior left the
    # Treatment:Day interaction effectively unconstrained.
    selected_priors <- bayes_scaled_priors(
      brms_formula, analysis_df, volume_column, prior_strength,
      time_column = time_column, include_sd = TRUE
    )
  }

  # ── Fit model ──────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    message("Fitting Bayesian LMM via brms (",
            n_chains, " chains × ", n_iter, " post-warmup draws, ",
            n_warmup, " warmup)...")
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
    backend      = backend,
    silent       = if (isTRUE(verbose)) 0L else 2L,
    refresh      = if (isTRUE(verbose)) 100L else 0L
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
  mcmc_diagnostics <- make_mcmc_diagnostics(posterior_summary,
                                            total_draws = n_chains * n_iter)
  nuts_diagnostics <- make_nuts_diagnostics(model)
  loo_diagnostics  <- bayes_loo(model)
  bayes_r2         <- bayes_r2_summary(model)
  ppc_coverage     <- bayes_ppc_coverage(model)

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
        SE            = if ("SE" %in% names(emm_df))
                          round(emm_df$SE, 3) else NA_real_,
        DF            = NA_real_,
        Lower_CrI     = round(emm_df[[lower_col]],    3),
        Upper_CrI     = round(emm_df[[upper_col]],    3),
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
          contrast    = as.character(pc_df$contrast),
          estimate    = round(pc_df$estimate,    4),
          SE          = if ("SE" %in% names(pc_df)) round(pc_df$SE, 4) else NA_real_,
          Lower_CrI   = round(pc_df[[lower_pc]], 4),
          Upper_CrI   = round(pc_df[[upper_pc]], 4),
          P_direction = emm_p_direction(pc, n_contrasts = nrow(pc_df)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ── Growth rates ──────────────────────────────────────────────────────────
  # When the LMM was fit with random slopes, derive per-animal growth rates
  # from the brms model itself (treatment fixed slope + per-animal random
  # slope draw) rather than falling back to OLS on log-volumes. The brms
  # version gives posterior credible intervals; the OLS version does not.
  # Other paths (intercept_only, GAM) fall back to OLS.
  brms_growth_rates <- NULL
  if (model_type == "lmm" &&
      random_effects_specification == "slope") {
    brms_growth_rates <- tryCatch(
      tg_brms_per_animal_growth_rates(model, treatment_column,
                                       id_column, time_column,
                                       analysis_df),
      error = function(e) NULL
    )
  }
  growth_rates <- if (!is.null(brms_growth_rates)) {
    brms_growth_rates
  } else {
    tgs_compute_growth_rates(
      auc_df, treatment_column, id_column,
      cage_column, time_column, volume_column
    )
  }

  # ── Data summary ───────────────────────────────────────────────────────────
  data_summary <- tgs_compute_summary(
    analysis_df, treatment_column, time_column, volume_column
  )

  # ── Plots ──────────────────────────────────────────────────────────────────
  pp_check_plot           <- NULL
  posterior_dist_plot     <- NULL
  prior_posterior_plot    <- NULL
  credible_intervals_plot <- NULL
  mcmc_trace_plot         <- NULL
  residuals_plot          <- NULL

  if (isTRUE(plots)) {

    # Posterior predictive check
    pp_check_plot <- tryCatch(
      brms::pp_check(model, type = "dens_overlay", ndraws = 50),
      error = function(e) NULL
    )

    # Treatment-parameter posterior densities and trace plots
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      draws_arr <- tryCatch(posterior::as_draws_array(model),
                            error = function(e) NULL)

      if (!is.null(draws_arr)) {
        all_pars <- dimnames(draws_arr)$variable
        tx_pars  <- grep(
          paste0("^b_",
                 gsub("([.^$*+?()\\[\\]{}|])", "\\\\\\1",
                      treatment_column, perl = TRUE)),
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

    # Prior vs posterior overlay
    prior_posterior_plot <- tryCatch(
      bayes_prior_posterior_plot(model, treatment_column),
      error = function(e) NULL
    )

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
            xmin = .data[["Lower_CrI"]],
            xmax = .data[["Upper_CrI"]],
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

    # Residuals vs Day — curvature diagnostic
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
          method   = "loess", se = TRUE, formula = y ~ x,
          colour   = "steelblue", linewidth = 0.8,
          fill     = "steelblue", alpha = 0.15
        ) +
        ggplot2::geom_hline(
          yintercept = 0, linetype = "dashed",
          colour = "grey50", linewidth = 0.5
        ) +
        ggplot2::facet_wrap(~ Treatment) +
        ggplot2::labs(
          title    = "Residuals vs. Study Day",
          subtitle = "Curvature indicates non-linear growth — consider a GAM or nonlinear model",
          x        = time_column,
          y        = paste0("Residual (", transform, " scale)")
        ) +
        ggplot2::theme_classic(base_size = 14)
    }, error = function(e) NULL)
  }

  # ── Analysis summary metadata ──────────────────────────────────────────────
  .prior_desc <- describe_priors(selected_priors)
  analysis_summary <- list(
    analysis_type = if (model_type == "gam") {
      "Bayesian Generalized Additive Mixed Model (brms, group-specific smooths)"
    } else {
      "Bayesian Linear Mixed-Effects Model (brms)"
    },
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
      fixed_effects      = paste(volume_column, "~",
                                 treatment_column, "*", time_column),
      random_effects     = re_term,
      cage_effect_modelled = use_cage_re,
      prior_strength     = prior_strength,
      transform          = transform
    ),
    methods = list(
      engine          = paste0("brms (", n_chains, " chains × ",
                               n_iter, " draws + ", n_warmup,
                               " warmup, seed = ", seed, ")"),
      # CODE_REVIEW.md R3.8 — report the priors actually handed to brms rather
      # than reconstructing them from b_sd, which would now misreport the
      # data-scaled Intercept and the per-coefficient rate priors.
      prior_b         = .prior_desc$prior_b,
      prior_intercept = .prior_desc$prior_intercept,
      prior_sd        = .prior_desc$prior_sd,
      prior_sigma     = .prior_desc$prior_sigma,
      prior_table     = .prior_desc$all,
      prior_scaling   = .prior_desc$scaling,
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
    model_type_used         = if (model_type == "gam") "bayes_tg_gam" else "bayes_tg",
    transform_used          = transform,
    meta = me_result_meta(
      analysis_type   = "Bayesian linear mixed-effects model (brms)",
      model_type_used = "bayes_tg",
      inference       = "bayesian",
      interval_type   = "credible",
      transform_used  = transform,
      estimate_scale  = switch(transform, log = "log volume", sqrt = "sqrt volume", "volume")
    ),
    summary                 = analysis_summary,
    posterior_summary       = posterior_summary,
    treatment_effects       = treatment_effects,
    pairwise_comparisons    = pairwise_comparisons,
    mcmc_diagnostics        = mcmc_diagnostics,
    nuts_diagnostics        = nuts_diagnostics,
    loo_diagnostics         = loo_diagnostics,
    bayes_R2                = bayes_r2,
    ppc_coverage            = ppc_coverage,
    pp_check_plot           = pp_check_plot,
    posterior_dist_plot     = posterior_dist_plot,
    prior_posterior_plot    = prior_posterior_plot,
    credible_intervals_plot = credible_intervals_plot,
    mcmc_trace_plot         = mcmc_trace_plot,
    residuals_plot          = residuals_plot,
    growth_rates            = growth_rates,
    data_summary            = data_summary,
    necrosis_summary        = necrosis_summary
  )
}


#' Per-animal posterior growth rates from a brms LMM
#'
#' Combines treatment-level fixed slope draws (Treatment[i]:Day) with each
#' animal's random slope draw (ranef ID, Day) to produce a posterior
#' distribution of growth rates per mouse. Returns a data frame compatible
#' with the lme4 path's \code{growth_rates} shape plus credible-interval
#' columns.
#'
#' Columns: \code{Treatment}, \code{ID}, \code{Cage} (if available),
#' \code{growth_rate} (posterior median), \code{growth_rate_lower} /
#' \code{growth_rate_upper} (2.5% / 97.5% CrI), \code{R_squared} (NA —
#' the OLS R^2 has no Bayesian counterpart at the per-animal level).
#'
#' Requires \pkg{brms} >= 2.19 and \pkg{posterior}.
#' @noRd
tg_brms_per_animal_growth_rates <- function(model, treatment_column,
                                             id_column, time_column,
                                             analysis_df) {
  # Fixed-effect posterior draws (full matrix): rows = draws, cols = params
  fe <- brms::fixef(model, summary = FALSE)
  fe_names <- colnames(fe)

  # Slope column for the reference treatment level: "<time_column>" alone.
  if (!time_column %in% fe_names) return(NULL)
  ref_slope_draws <- fe[, time_column]

  # Treatment:Day interaction columns: one per non-reference treatment level.
  int_cols <- grep(
    paste0(treatment_column, ".*:", time_column, "|",
           time_column, ":", treatment_column),
    fe_names, value = TRUE
  )

  # Per-animal random slope draws on Day. brms::ranef(summary = FALSE)
  # returns a list of 3D arrays (draws x animals x effects).
  re_all <- tryCatch(brms::ranef(model, summary = FALSE),
                     error = function(e) NULL)
  if (is.null(re_all) || !id_column %in% names(re_all)) return(NULL)
  re_arr <- re_all[[id_column]]
  if (length(dim(re_arr)) != 3L) return(NULL)
  re_dimnames <- dimnames(re_arr)
  if (length(re_dimnames) < 3L) return(NULL)

  # Identify which slice corresponds to the Day random slope
  slope_dim <- which(re_dimnames[[3L]] == time_column)
  if (length(slope_dim) == 0L) return(NULL)
  re_slope <- re_arr[, , slope_dim[[1L]], drop = TRUE]
  # re_slope is now n_draws x n_animals.

  animal_ids <- colnames(re_slope)

  # Build per-animal lookup of Treatment (and Cage if present)
  cage_col_present <- "Cage" %in% colnames(analysis_df)
  id_lookup <- unique(analysis_df[, c(treatment_column, id_column,
                                       intersect("Cage", colnames(analysis_df))),
                                   drop = FALSE])
  id_lookup[[id_column]] <- as.character(id_lookup[[id_column]])

  rows <- lapply(animal_ids, function(aid) {
    meta <- id_lookup[id_lookup[[id_column]] == aid, , drop = FALSE]
    if (nrow(meta) == 0L) return(NULL)
    treat <- as.character(meta[[treatment_column]])[[1L]]
    cage  <- if (cage_col_present) as.character(meta[["Cage"]])[[1L]] else NA_character_

    # Build treatment fixed slope = reference slope + interaction (if any)
    # for this animal's treatment level. brms encodes interaction columns as
    # "<treatment_column><level>:<time>" — drop the reference level.
    int_for_treat <- grep(
      paste0("^", treatment_column,
             gsub("([.^$*+?()\\[\\]{}|])", "\\\\\\1", treat,
                  perl = TRUE),
             ":", time_column, "$"),
      int_cols, value = TRUE
    )
    if (length(int_for_treat) == 1L) {
      treat_slope_draws <- ref_slope_draws + fe[, int_for_treat]
    } else {
      # Treatment is the reference level (or naming differs) — use ref slope.
      treat_slope_draws <- ref_slope_draws
    }
    animal_slope_draws <- treat_slope_draws + re_slope[, aid]

    qq <- stats::quantile(animal_slope_draws,
                          c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(
      Treatment         = treat,
      ID                = aid,
      Cage              = cage,
      growth_rate       = round(qq[[2L]], 4),
      growth_rate_lower = round(qq[[1L]], 4),
      growth_rate_upper = round(qq[[3L]], 4),
      R_squared         = NA_real_,
      stringsAsFactors  = FALSE
    )
  })

  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out) || nrow(out) == 0L) return(NULL)
  if (!cage_col_present) out$Cage <- NULL
  out
}
