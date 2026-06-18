# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Bayesian Dose-Response Analysis (Hill / Emax Model)
#'
#' Fits a Bayesian Hill/Emax inhibition model to dose-response tumour data
#' using \pkg{brms} (Bürkner, 2017). The model is parameterised on the Tumour
#' Growth Inhibition (TGI) scale relative to the vehicle control, providing
#' posterior estimates of EC50, Emax, and the Hill slope with full credible
#' intervals.
#'
#' \strong{Inhibition-only by construction:} the non-linear formula
#' \code{TGI ~ inv_logit(logEmax) / (1 + (exp(logEC50)/Dose)^exp(logHill))}
#' bounds the asymptotic effect in \code{[0, 1]} and constrains the Hill
#' slope to be positive, so the curve is monotonically increasing in dose
#' toward a positive ceiling — i.e. inhibitory only. Data that are actually
#' stimulatory (treatment accelerates growth, negative TGI) cannot be
#' represented; the function emits a warning at fit time and the resulting
#' posterior will not be credible. For stimulatory dose-response, model
#' \eqn{V_{treated}} directly with a stimulatory Hill parameterisation or
#' invert the TGI sign convention.
#'
#' @param df Data frame containing dose-response data.
#' @param dose_column Column name for dose values (numeric, including 0 for
#'   control). Default \code{"Dose"}.
#' @param volume_column Column name for tumour volume at the endpoint day.
#'   Default \code{"Volume"}.
#' @param treatment_column Column name for treatment group. Default
#'   \code{"Treatment"}.
#' @param id_column Column name for individual animal ID. Default \code{"ID"}.
#' @param day_column Column name for study day. Default \code{"Day"}.
#' @param endpoint_day Numeric. Study day to use as the endpoint for TGI
#'   computation. \code{NULL} (default) uses the last observed day.
#' @param reference_group Name of the vehicle/control group (dose = 0).
#'   Auto-detected from common control names if \code{NULL}.
#' @param prior_strength Prior preset. The EC50 prior is \strong{data-aware}:
#'   centred on \code{log(median(non_zero_doses))} so that, regardless of
#'   the units the user supplies doses in (nM, mg/kg, mM…), the prior puts
#'   most of its mass in the range the user actually tested. The width
#'   varies by preset (see below). This is intentionally weakly-informative
#'   on dose magnitude while preserving a flat prior on relative dose
#'   spacing. If the experiment didn't span the true EC50, the posterior
#'   will be pulled toward the median observed dose — supply \code{prior_ec50}
#'   manually if you have external EC50 knowledge.
#'   \describe{
#'     \item{\code{"skeptical"}}{(default) Emax ~ N(1.5, 0.75) on logit scale
#'       (prior median ≈ 0.82); Hill ~ N(0, 0.4) on log scale (prior median
#'       = 1.0); EC50 centred at the geometric mean dose ± 1.5 SD on log
#'       scale.}
#'     \item{\code{"weakly_informative"}}{Emax ~ N(1.5, 1.5), Hill ~ N(0, 0.8),
#'       EC50 centred at geometric mean dose ± 2 SD.}
#'     \item{\code{"informative"}}{Emax ~ N(1.5, 1.0), Hill ~ N(0, 0.6),
#'       EC50 centred at geometric mean dose ± 1.75 SD.}
#'     \item{\code{"diffuse"}}{Very wide priors for exploratory analysis.}
#'     \item{\code{"manual"}}{Use \code{prior_emax}, \code{prior_ec50},
#'       \code{prior_hill}, and \code{prior_sigma} directly.}
#'   }
#' @param prior_emax brms prior string for the logit-transformed Emax
#'   parameter (\code{nlpar = "logEmax"}). Only used when
#'   \code{prior_strength = "manual"}.
#' @param prior_ec50 brms prior string for the log-transformed EC50
#'   (\code{nlpar = "logEC50"}). Only used when
#'   \code{prior_strength = "manual"}.
#' @param prior_hill brms prior string for the log-transformed Hill slope
#'   (\code{nlpar = "logHill"}). Only used when
#'   \code{prior_strength = "manual"}.
#' @param prior_sigma brms prior string for the residual SD (\code{class =
#'   "sigma"}). Only used when \code{prior_strength = "manual"}.
#' @param n_chains Number of MCMC chains. Default \code{4}.
#' @param n_warmup Warm-up (burn-in) iterations per chain. Default \code{1000}.
#' @param n_iter Post-warmup draws per chain. Default \code{500}.
#' @param seed Integer random seed. Default \code{42}.
#' @param return_model Logical. Return the \code{brmsfit} object? Default
#'   \code{TRUE}.
#' @param plots Logical. Pre-compute and return \pkg{ggplot2} plot objects?
#'   Default \code{TRUE}.
#' @param verbose Logical. Show Stan compilation and sampling messages?
#'   Default \code{FALSE}.
#' @param backend brms backend: \code{"rstan"} (default) or \code{"cmdstanr"}.
#'   See \code{\link{bayesian_tumor_growth}} for details.
#'
#' @details
#' The Hill inhibition model fitted is:
#' \deqn{
#'   \mathrm{TGI}_i =
#'   \frac{E_{\max}}{1 + (EC_{50} / \mathrm{Dose}_i)^{\mathrm{Hill}}}
#'   + \varepsilon_i
#' }
#' where \eqn{E_\text{max} = \text{inv\_logit}(\theta_1) \in (0,1)},
#' \eqn{EC_{50} = \exp(\theta_2) > 0}, and
#' \eqn{\text{Hill} = \exp(\theta_3) > 0}.
#' TGI is computed as \eqn{1 - V_i / \bar{V}_\text{control}} using the mean
#' volume of the reference group at the endpoint day. The model is fitted on
#' treated animals only (Dose > 0); the reference group defines TGI = 0 at
#' Dose = 0 by construction.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{model}}{\code{brmsfit} object or \code{NULL}.}
#'   \item{\code{model_type_used}}{Character \code{"bayes_dr"}.}
#'   \item{\code{summary}}{Named list of analysis metadata.}
#'   \item{\code{posterior_summary}}{Data frame of raw (log/logit-scale)
#'     parameter posterior summaries (median, 95 % CrI, Rhat, ESS).}
#'   \item{\code{dr_parameters}}{Data frame with one row per parameter (EC50,
#'     Emax, Hill) on the interpretable back-transformed scale: posterior
#'     median, 2.5 % and 97.5 % CrI.}
#'   \item{\code{mcmc_diagnostics}}{Data frame with Rhat and Converged per
#'     parameter.}
#'   \item{\code{dose_response_summary}}{Data frame: Dose, N, mean and SD of
#'     observed TGI, posterior-predicted median TGI, 95 % CrI — one row per
#'     dose level.}
#'   \item{\code{dose_response_curve_plot}}{Dose-response curve: posterior
#'     median ± 95 % CrI ribbon over a dense dose grid, overlaid with observed
#'     TGI points. \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{prior_posterior_plot}}{Prior vs posterior density overlay for
#'     the three Hill model parameters. \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{pp_check_plot}}{Posterior predictive density check.
#'     \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{mcmc_trace_plot}}{MCMC trace for the three parameters.
#'     \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{tgi_data}}{Data frame used for fitting: one row per treated
#'     animal with columns \code{Dose} and \code{TGI}.}
#'   \item{\code{control_mean_volume}}{Mean control volume at the endpoint day
#'     used to compute TGI.}
#' }
#'
#' @references
#' Bürkner, P.-C. (2017). brms: An R package for Bayesian multilevel models
#' using Stan. \emph{Journal of Statistical Software}, 80(1), 1–28.
#' \doi{10.18637/jss.v080.i01}
#'
#' @importFrom stats as.formula median quantile sd
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point geom_density
#'   geom_vline scale_fill_manual scale_colour_manual facet_wrap labs
#'   theme theme_classic
#' @export
bayesian_dose_response <- function(
  df,
  dose_column      = "Dose",
  volume_column    = "Volume",
  treatment_column = "Treatment",
  id_column        = "ID",
  day_column       = "Day",
  endpoint_day     = NULL,
  reference_group  = NULL,
  prior_strength   = c("skeptical", "weakly_informative",
                       "informative", "diffuse", "manual"),
  prior_emax       = NULL,
  prior_ec50       = NULL,
  prior_hill       = NULL,
  prior_sigma      = NULL,
  n_chains         = 4L,
  n_warmup         = 1000L,
  n_iter           = 500L,
  seed             = 42L,
  return_model     = TRUE,
  plots            = TRUE,
  verbose          = FALSE,
  backend          = c("rstan", "cmdstanr"),
  mcmc             = NULL
) {

  # ── Dependency check ───────────────────────────────────────────────────────
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required for Bayesian analysis.\n",
      "Install it with: install.packages('brms')"
    )
  }

  prior_strength <- match.arg(prior_strength)
  backend        <- resolve_brms_backend(backend)

  # CODE_REVIEW.md Round 2 D.2 — only the MCMC config helper applies here;
  # the Hill-model priors (`prior_emax`, `prior_ec50`, `prior_hill`) are
  # model-specific so they're not bundled by `tg_priors()`.
  .m <- .resolve_mcmc(mcmc, n_chains, n_warmup, n_iter, seed, backend)
  n_chains <- .m$chains
  n_warmup <- .m$warmup
  n_iter   <- .m$iter
  seed     <- .m$seed
  backend  <- .m$backend

  # ── Column validation ──────────────────────────────────────────────────────
  required_cols <- c(dose_column, volume_column, treatment_column,
                     id_column, day_column)
  missing_cols  <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ── Reference group ────────────────────────────────────────────────────────
  treatment_groups <- unique(as.character(df[[treatment_column]]))
  ctrl_candidates  <- c("Control", "Vehicle", "control", "vehicle",
                        "CTRL", "ctrl")
  if (is.null(reference_group)) {
    reference_group <- ctrl_candidates[ctrl_candidates %in% treatment_groups]
    if (length(reference_group) == 0L) {
      reference_group <- sort(treatment_groups)[1L]
    } else {
      reference_group <- reference_group[1L]
    }
  } else if (!reference_group %in% treatment_groups) {
    stop(
      "Reference group '", reference_group,
      "' not found in '", treatment_column, "'."
    )
  }

  # ── Endpoint observations per mouse ────────────────────────────────────────
  # When endpoint_day is unspecified ("All dates") we mirror the frequentist
  # path: take each mouse's LAST observation. This avoids dropping reference
  # animals euthanised before the global max day (typical in dose-escalation
  # studies where controls reach the IACUC volume limit first).
  all_days <- sort(unique(as.numeric(df[[day_column]])))
  ep_label <- if (!is.null(endpoint_day)) {
    paste0("day ", endpoint_day)
  } else {
    "last observation per mouse"
  }

  if (!is.null(endpoint_day)) {
    ep_day <- as.numeric(endpoint_day)
    ep_df  <- df[as.numeric(df[[day_column]]) == ep_day, ]
    if (nrow(ep_df) == 0L) {
      stop(
        "No observations found at endpoint day ", ep_day,
        ". Available days: ", paste(all_days, collapse = ", ")
      )
    }
  } else {
    ep_df <- df %>%
      dplyr::group_by(.data[[id_column]]) %>%
      dplyr::filter(as.numeric(.data[[day_column]]) ==
                      max(as.numeric(.data[[day_column]]), na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }

  # ── Control mean volume ────────────────────────────────────────────────────
  ctrl_vols <- as.numeric(
    ep_df[[volume_column]][ep_df[[treatment_column]] == reference_group]
  )
  if (length(ctrl_vols) == 0L || all(is.na(ctrl_vols))) {
    stop(
      "Reference group '", reference_group,
      "' has no valid volume observations at ", ep_label, "."
    )
  }
  ctrl_mean <- mean(ctrl_vols, na.rm = TRUE)

  # ── TGI data for treated animals ───────────────────────────────────────────
  treated_df <- ep_df[ep_df[[treatment_column]] != reference_group, ]
  if (nrow(treated_df) == 0L) {
    stop("No treated (non-reference) animals found at endpoint day ", ep_day)
  }

  tgi_data <- data.frame(
    ID        = as.character(treated_df[[id_column]]),
    Dose      = as.numeric(treated_df[[dose_column]]),
    TGI       = 1 - as.numeric(treated_df[[volume_column]]) / ctrl_mean,
    Treatment = as.character(treated_df[[treatment_column]]),
    stringsAsFactors = FALSE
  )

  # Remove zero-dose treated rows (should be reference only, but safety check)
  tgi_data <- tgi_data[tgi_data$Dose > 0, ]
  if (nrow(tgi_data) == 0L) {
    stop("All treated animals have Dose = 0. Cannot fit Hill model.")
  }

  non_zero_doses <- tgi_data$Dose[tgi_data$Dose > 0]
  log_median_dose <- median(log(non_zero_doses))

  # ── Direction sanity check ─────────────────────────────────────────────────
  # The Hill formula used below — TGI ~ inv_logit(logEmax) /
  # (1 + (exp(logEC50)/Dose)^exp(logHill)) — has inv_logit(logEmax) in [0, 1]
  # and exp(logHill) > 0, so the curve is *monotonically increasing in dose
  # toward a positive ceiling*. In TGI space (1 - V_treated/V_control), that
  # corresponds to inhibitory dose-response only. If the underlying data show
  # negative-mean TGI across doses (treatment accelerates growth), the model
  # cannot represent it and the fit will be poor with no clear signal as to
  # why. Surface this clearly.
  if (mean(tgi_data$TGI, na.rm = TRUE) < -0.10) {
    warning(
      "bayesian_dose_response(): the supplied data have negative mean TGI ",
      "across doses (treatment appears to accelerate growth). The Hill ",
      "model used here is *inhibition-only by construction* (Emax bounded ",
      "in [0, 1] and Hill > 0) and cannot represent stimulatory dose-",
      "response. The fit will likely be poor with non-credible parameters. ",
      "Consider modelling V_treated directly with a stimulatory Hill ",
      "formula or inverting the sign convention.",
      call. = FALSE
    )
  }

  # ── Prior specification ────────────────────────────────────────────────────
  if (prior_strength == "manual") {
    if (any(is.null(c(prior_emax, prior_ec50, prior_hill, prior_sigma)))) {
      stop(
        "When prior_strength = 'manual', ",
        "all four prior_* arguments must be supplied."
      )
    }
    selected_priors <- c(
      brms::prior_string(prior_emax,  nlpar = "logEmax"),
      brms::prior_string(prior_ec50,  nlpar = "logEC50"),
      brms::prior_string(prior_hill,  nlpar = "logHill"),
      brms::prior_string(prior_sigma, class = "sigma")
    )
  } else {
    emax_sd  <- switch(prior_strength,
      skeptical          = 0.75,
      weakly_informative = 1.5,
      informative        = 1.0,
      diffuse            = 3.0
    )
    ec50_sd  <- switch(prior_strength,
      skeptical          = 1.5,
      weakly_informative = 2.0,
      informative        = 1.75,
      diffuse            = 4.0
    )
    hill_sd  <- switch(prior_strength,
      skeptical          = 0.4,
      weakly_informative = 0.8,
      informative        = 0.6,
      diffuse            = 2.0
    )
    sig_rate <- switch(prior_strength,
      skeptical          = 4,
      weakly_informative = 2,
      informative        = 3,
      diffuse            = 0.5
    )
    selected_priors <- c(
      brms::prior_string(
        paste0("normal(1.5, ", emax_sd, ")"),
        nlpar = "logEmax"
      ),
      brms::prior_string(
        paste0("normal(", round(log_median_dose, 3), ", ", ec50_sd, ")"),
        nlpar = "logEC50"
      ),
      brms::prior_string(
        paste0("normal(0, ", hill_sd, ")"),
        nlpar = "logHill"
      ),
      brms::prior_string(
        paste0("exponential(", sig_rate, ")"),
        class = "sigma"
      )
    )
  }

  # ── brms nonlinear formula ─────────────────────────────────────────────────
  brms_formula <- brms::bf(
    TGI ~ inv_logit(logEmax) / (1 + (exp(logEC50) / Dose)^exp(logHill)),
    logEmax ~ 1,
    logEC50 ~ 1,
    logHill ~ 1,
    nl = TRUE
  )

  # ── Fit ────────────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) {
    message(
      "Fitting Bayesian Hill model via brms (",
      n_chains, " chains × ", n_iter, " post-warmup draws, ",
      n_warmup, " warmup)..."
    )
  }

  model <- brms::brm(
    formula      = brms_formula,
    data         = tgi_data,
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

  # ── Posterior summary (raw log/logit scale) ────────────────────────────────
  brms_smry  <- summary(model)
  fixed_df   <- as.data.frame(brms_smry$fixed)
  fixed_df   <- cbind(
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
  ppc_coverage     <- bayes_ppc_coverage(model)

  # ── Back-transformed parameter posteriors ──────────────────────────────────
  draws <- tryCatch(brms::as_draws_df(model), error = function(e) NULL)

  dr_parameters <- NULL
  if (!is.null(draws)) {
    log_emax_draws <- draws[["b_logEmax_Intercept"]]
    log_ec50_draws <- draws[["b_logEC50_Intercept"]]
    log_hill_draws <- draws[["b_logHill_Intercept"]]

    inv_logit <- function(x) 1 / (1 + exp(-x))

    make_param_row <- function(name, vals, trans_fn) {
      bv  <- trans_fn(vals)
      q   <- stats::quantile(bv, c(0.025, 0.5, 0.975), na.rm = TRUE)
      data.frame(
        Parameter  = name,
        Median     = round(q["50%"],   4),
        Lower_CrI  = round(q["2.5%"],  4),
        Upper_CrI  = round(q["97.5%"], 4),
        stringsAsFactors = FALSE
      )
    }

    dr_parameters <- rbind(
      make_param_row("EC50", log_ec50_draws, exp),
      make_param_row("Emax", log_emax_draws, inv_logit),
      make_param_row("Hill", log_hill_draws, exp)
    )
    rownames(dr_parameters) <- NULL
  }

  # ── Dose-response summary per dose level ───────────────────────────────────
  dose_levels <- sort(unique(tgi_data$Dose))
  dose_grid_fit <- data.frame(
    Dose = dose_levels, stringsAsFactors = FALSE
  )

  epred_fit <- tryCatch(
    brms::posterior_epred(
      model, newdata = dose_grid_fit, allow_new_levels = TRUE
    ),
    error = function(e) NULL
  )

  dose_response_summary <- do.call(rbind, lapply(
    seq_along(dose_levels),
    function(i) {
      d      <- dose_levels[i]
      obs    <- tgi_data$TGI[tgi_data$Dose == d]
      pred_q <- if (!is.null(epred_fit)) {
        stats::quantile(epred_fit[, i], c(0.025, 0.5, 0.975))
      } else {
        c("2.5%" = NA, "50%" = NA, "97.5%" = NA)
      }
      data.frame(
        Dose              = d,
        N                 = length(obs),
        Observed_TGI_Mean = round(mean(obs, na.rm = TRUE), 3),
        Observed_TGI_SD   = round(stats::sd(obs, na.rm = TRUE), 3),
        Predicted_Median  = round(pred_q["50%"],   3),
        Predicted_Lower   = round(pred_q["2.5%"],  3),
        Predicted_Upper   = round(pred_q["97.5%"], 3),
        stringsAsFactors  = FALSE
      )
    }
  ))
  rownames(dose_response_summary) <- NULL

  # ── Plots ──────────────────────────────────────────────────────────────────
  dose_response_curve_plot <- NULL
  prior_posterior_plot_dr  <- NULL
  pp_check_plot            <- NULL
  mcmc_trace_plot          <- NULL

  if (isTRUE(plots)) {

    pp_check_plot <- tryCatch(
      brms::pp_check(model, type = "dens_overlay", ndraws = 50),
      error = function(e) NULL
    )

    if (requireNamespace("bayesplot", quietly = TRUE)) {
      draws_arr <- tryCatch(posterior::as_draws_array(model),
                            error = function(e) NULL)
      if (!is.null(draws_arr)) {
        all_pars  <- dimnames(draws_arr)$variable
        hill_pars <- grep("^b_log", all_pars, value = TRUE)
        if (length(hill_pars) > 0) {
          mcmc_trace_plot <- tryCatch(
            bayesplot::mcmc_trace(draws_arr, pars = hill_pars),
            error = function(e) NULL
          )
        }
      }
    }

    # Dose-response curve: posterior credible band
    dose_response_curve_plot <- tryCatch({
      max_dose   <- max(tgi_data$Dose)
      dose_seq   <- seq(max_dose / 200, max_dose * 1.1, length.out = 200)
      curve_nd   <- data.frame(Dose = dose_seq)

      ep_curve   <- brms::posterior_epred(
        model, newdata = curve_nd, allow_new_levels = TRUE
      )
      curve_nd$Median <- apply(ep_curve, 2L, stats::median)
      curve_nd$Lower  <- apply(
        ep_curve, 2L, function(x) stats::quantile(x, 0.025)
      )
      curve_nd$Upper  <- apply(
        ep_curve, 2L, function(x) stats::quantile(x, 0.975)
      )

      ec50_med <- if (!is.null(dr_parameters)) {
        dr_parameters$Median[dr_parameters$Parameter == "EC50"]
      } else {
        NA_real_
      }

      ggplot2::ggplot() +
        ggplot2::geom_ribbon(
          data = curve_nd,
          ggplot2::aes(
            x    = .data[["Dose"]],
            ymin = .data[["Lower"]],
            ymax = .data[["Upper"]]
          ),
          fill = "steelblue", alpha = 0.2
        ) +
        ggplot2::geom_line(
          data = curve_nd,
          ggplot2::aes(x = .data[["Dose"]], y = .data[["Median"]]),
          colour = "steelblue", linewidth = 1
        ) +
        ggplot2::geom_point(
          data = tgi_data,
          ggplot2::aes(x = .data[["Dose"]], y = .data[["TGI"]]),
          alpha = 0.6, size = 2.5, colour = "grey30"
        ) +
        ggplot2::geom_vline(
          xintercept = ec50_med, linetype = "dashed",
          colour = "tomato", linewidth = 0.7
        ) +
        ggplot2::geom_hline(
          yintercept = c(0, 0.5), linetype = c("solid", "dashed"),
          colour = c("grey70", "tomato"), linewidth = 0.5
        ) +
        ggplot2::labs(
          title    = "Bayesian Dose-Response Curve (Hill / Emax Model)",
          subtitle = paste0(
            "Posterior median ± 95 % CrI;",
            " dashed red = EC50 (", round(ec50_med, 2), ")"
          ),
          x = paste0("Dose (", dose_column, ")"),
          y = "TGI (fraction of control)"
        ) +
        ggplot2::theme_classic(base_size = 14)
    }, error = function(e) NULL)

    # Prior vs posterior for Hill parameters
    prior_posterior_plot_dr <- tryCatch({
      if (is.null(draws)) {
        return(NULL)
      }
      params <- list(
        list(raw = "b_logEmax_Intercept",
             prior = "prior_b_logEmax_Intercept",
             label = "Emax (logit scale)", trans = "none"),
        list(raw = "b_logEC50_Intercept",
             prior = "prior_b_logEC50_Intercept",
             label = "EC50 (log scale)", trans = "none"),
        list(raw = "b_logHill_Intercept",
             prior = "prior_b_logHill_Intercept",
             label = "Hill (log scale)", trans = "none")
      )
      rows <- lapply(params, function(p) {
        post_col  <- if (p$raw   %in% names(draws)) draws[[p$raw]]   else NULL
        prior_col <- if (p$prior %in% names(draws)) draws[[p$prior]] else NULL
        if (is.null(post_col)) {
          return(NULL)
        }
        post_df <- data.frame(
          Parameter = p$label, Value = post_col,
          Source = "Posterior", stringsAsFactors = FALSE
        )
        if (!is.null(prior_col)) {
          rbind(
            post_df,
            data.frame(
              Parameter = p$label, Value = prior_col,
              Source = "Prior", stringsAsFactors = FALSE
            )
          )
        } else {
          post_df
        }
      })
      plot_df <- do.call(rbind, rows[!sapply(rows, is.null)])
      if (nrow(plot_df) == 0L) {
        return(NULL)
      }
      plot_df$Source <- factor(
        plot_df$Source, levels = c("Prior", "Posterior")
      )

      ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x      = .data[["Value"]],
          fill   = .data[["Source"]],
          colour = .data[["Source"]]
        )
      ) +
        ggplot2::geom_density(alpha = 0.35, linewidth = 0.6) +
        ggplot2::geom_vline(
          xintercept = 0, linetype = "dashed",
          colour = "grey30", linewidth = 0.5
        ) +
        ggplot2::facet_wrap(~ Parameter, scales = "free") +
        ggplot2::scale_fill_manual(
          values = c(Prior = "grey60", Posterior = "steelblue")
        ) +
        ggplot2::scale_colour_manual(
          values = c(Prior = "grey40", Posterior = "steelblue4")
        ) +
        ggplot2::labs(
          title    = "Prior vs. Posterior Distributions",
          subtitle = "Hill model parameters on log / logit scale",
          x        = "Parameter value",
          y        = "Density",
          fill     = NULL,
          colour   = NULL
        ) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::theme(legend.position = "top")
    }, error = function(e) NULL)
  }

  # ── Analysis summary metadata ──────────────────────────────────────────────
  analysis_summary <- list(
    analysis_type = paste0(
      "Bayesian Hill/Emax Dose-Response Model (brms); ",
      "TGI = Emax / (1 + (EC50/Dose)^Hill)"
    ),
    data_description = list(
      endpoint_day      = ep_day,
      reference_group   = reference_group,
      control_mean_vol  = round(ctrl_mean, 2),
      n_dose_levels     = length(dose_levels),
      dose_levels       = dose_levels,
      n_treated_animals = nrow(tgi_data)
    ),
    model_specification = list(
      formula        = paste0(
        "TGI ~ inv_logit(logEmax)",
        " / (1 + (exp(logEC50)/Dose)^exp(logHill))"
      ),
      prior_strength = prior_strength
    ),
    methods = list(
      engine = paste0(
        "brms (", n_chains, " chains × ",
        n_iter, " draws + ", n_warmup, " warmup, seed = ", seed, ")"
      ),
      prior_logEmax = if (prior_strength == "manual") {
        prior_emax
      } else {
        paste0("normal(1.5, ", emax_sd, ")")
      },
      prior_logEC50 = if (prior_strength == "manual") {
        prior_ec50
      } else {
        paste0(
          "normal(", round(log_median_dose, 3), ", ", ec50_sd, ")"
        )
      },
      prior_logHill = if (prior_strength == "manual") {
        prior_hill
      } else {
        paste0("normal(0, ", hill_sd, ")")
      }
    )
  )

  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    model                    = if (isTRUE(return_model)) model else NULL,
    model_type_used          = "bayes_dr",
    summary                  = analysis_summary,
    posterior_summary        = posterior_summary,
    dr_parameters            = dr_parameters,
    mcmc_diagnostics         = mcmc_diagnostics,
    nuts_diagnostics         = nuts_diagnostics,
    loo_diagnostics          = loo_diagnostics,
    bayes_R2                 = bayes_r2,
    ppc_coverage             = ppc_coverage,
    dose_response_summary    = dose_response_summary,
    dose_response_curve_plot = dose_response_curve_plot,
    prior_posterior_plot     = prior_posterior_plot_dr,
    pp_check_plot            = pp_check_plot,
    mcmc_trace_plot          = mcmc_trace_plot,
    tgi_data                 = tgi_data,
    control_mean_volume      = ctrl_mean
  )
}
