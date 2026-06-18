# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Shared helper utilities used internally by all six Bayesian analysis
# functions: bayesian_tumor_growth(), bayesian_body_weight(),
# bayesian_survival(), bayesian_dose_response(), bayesian_synergy(), and
# bayesian_therapeutic_window().

# ── Back-transform ─────────────────────────────────────────────────────────────

#' Back-transform posterior predictions to the original measurement scale
#'
#' @param x Numeric vector or matrix of predictions on the modelling scale.
#' @param transform Character: \code{"log"}, \code{"sqrt"}, or \code{"none"}.
#' @noRd
bayes_backtransform <- function(x, transform) {
  switch(transform, log = exp(x), sqrt = x ^ 2, x)
}


# ── MCMC diagnostics ───────────────────────────────────────────────────────────

#' Build a standardised MCMC diagnostics data frame from a posterior summary
#'
#' Expects a data frame with columns \code{Parameter}, \code{Rhat},
#' \code{Bulk_ESS}, and \code{Tail_ESS} as produced by the
#' \code{summary(brmsfit)$fixed} workflow used throughout this package.
#' Adds an \code{ESS_per_draw} column (Bulk_ESS / total_draws) when the total
#' number of draws is supplied; ratios < 0.1 indicate an inefficient
#' parameterisation per Stan recommendations.
#' @noRd
make_mcmc_diagnostics <- function(posterior_summary_df, total_draws = NULL) {
  d <- data.frame(
    Parameter = posterior_summary_df$Parameter,
    Rhat      = round(posterior_summary_df$Rhat,     4),
    Bulk_ESS  = round(posterior_summary_df$Bulk_ESS, 0),
    Tail_ESS  = round(posterior_summary_df$Tail_ESS, 0),
    Converged = posterior_summary_df$Rhat <= 1.01,
    # CODE_REVIEW.md DIAGNOSTICS gap (9) — ESS adequacy flag mirrors the
    # Converged flag for ESS. 400 is the de-facto minimum from the Stan
    # team's guidance; below it the posterior is noisier than the chain
    # length suggests.
    ESS_Adequate = (posterior_summary_df$Bulk_ESS >= 400) &
                   (posterior_summary_df$Tail_ESS >= 400),
    stringsAsFactors = FALSE
  )
  if (!is.null(total_draws) && is.numeric(total_draws) && total_draws > 0) {
    d$ESS_per_draw <- round(d$Bulk_ESS / total_draws, 3)
  }
  d
}


#' Extract NUTS sampler-level diagnostics from a brmsfit
#'
#' Surfaces the three first-line Stan/NUTS sampler diagnostics that Rhat /
#' ESS alone cannot detect:
#' \itemize{
#'   \item \code{n_divergent} — number of divergent transitions across all
#'     post-warmup iterations. Any non-zero count signals untrustworthy
#'     posterior near the divergent points. \code{> 0} \emph{should} be
#'     investigated even when Rhat = 1.00.
#'   \item \code{n_max_treedepth} — number of iterations that hit
#'     \code{max_treedepth}. Indicates inefficient exploration (not
#'     incorrect, but slow); reparameterise or raise treedepth.
#'   \item \code{ebfmi_min} — minimum per-chain E-BFMI (energy Bayesian
#'     fraction of missing information). Values < 0.3 indicate the sampler
#'     is struggling to explore the typical set.
#' }
#'
#' \code{clean = TRUE} iff all three pass (no divergences, no max-treedepth
#' hits, E-BFMI \eqn{\ge} 0.3).
#'
#' Returns a one-row data frame so it can be displayed in a dashboard table
#' alongside the per-parameter diagnostics.
#' @noRd
make_nuts_diagnostics <- function(model) {
  default <- data.frame(
    n_divergent     = NA_integer_,
    n_max_treedepth = NA_integer_,
    ebfmi_min       = NA_real_,
    clean           = NA,
    stringsAsFactors = FALSE
  )
  if (!requireNamespace("brms", quietly = TRUE) || is.null(model)) return(default)

  np <- tryCatch(brms::nuts_params(model), error = function(e) NULL)
  if (is.null(np) || !is.data.frame(np) || !"Parameter" %in% colnames(np)) {
    return(default)
  }

  div <- np$Value[np$Parameter == "divergent__"]
  td  <- np$Value[np$Parameter == "treedepth__"]
  max_td_val <- if (length(td) > 0) max(td, na.rm = TRUE) else NA_real_
  n_div   <- if (length(div) > 0) as.integer(sum(div > 0, na.rm = TRUE)) else NA_integer_
  n_maxtd <- if (length(td) > 0 && is.finite(max_td_val)) {
    as.integer(sum(td >= max_td_val, na.rm = TRUE))
  } else NA_integer_

  ebfmi <- tryCatch({
    if (utils::packageVersion("brms") >= "2.20.0" &&
        exists("ebfmi", where = asNamespace("brms"))) {
      vals <- brms::ebfmi(model)
      if (length(vals) > 0) min(vals, na.rm = TRUE) else NA_real_
    } else {
      # Fall back to rstan if available
      if (requireNamespace("rstan", quietly = TRUE) && !is.null(model$fit)) {
        vals <- rstan::get_bfmi(model$fit)
        if (length(vals) > 0) min(vals, na.rm = TRUE) else NA_real_
      } else NA_real_
    }
  }, error = function(e) NA_real_)

  clean_flag <- isTRUE(n_div == 0L) &&
                isTRUE(is.finite(ebfmi) && ebfmi >= 0.3)

  data.frame(
    n_divergent     = n_div,
    n_max_treedepth = n_maxtd,
    ebfmi_min       = if (is.finite(ebfmi)) round(ebfmi, 3) else NA_real_,
    clean           = clean_flag,
    stringsAsFactors = FALSE
  )
}


# ── Prior presets ──────────────────────────────────────────────────────────────

#' Return prior hyperparameters for a named prior-strength preset
#'
#' Returns a list with \code{b_sd} (normal SD for fixed-effect coefficients)
#' and \code{exp_rate} (rate for Exponential priors on SD and sigma parameters).
#' Used by bayesian_tumor_growth, bayesian_body_weight, bayesian_survival, and
#' bayesian_synergy.
#' @noRd
bayes_prior_params <- function(prior_strength) {
  switch(prior_strength,
    skeptical          = list(b_sd = 0.25, exp_rate = 2),
    weakly_informative = list(b_sd = 1.0,  exp_rate = 1),
    informative        = list(b_sd = 0.5,  exp_rate = 2),
    diffuse            = list(b_sd = 2.5,  exp_rate = 0.5),
    stop("Unknown prior_strength: '", prior_strength,
         "'. Expected one of skeptical, weakly_informative, ",
         "informative, diffuse, or manual (manual is handled by the ",
         "caller).")
  )
}


# ── Cage column setup ──────────────────────────────────────────────────────────

#' Validate and resolve a brms backend choice
#'
#' Returns the validated backend choice. \code{"rstan"} (default) is always
#' available because brms imports rstan. \code{"cmdstanr"} requires the
#' cmdstanr package plus a working CmdStan installation; if either is
#' missing we \code{stop()} with a pointer to install instructions rather
#' than silently fall back, so users who chose cmdstanr know why their
#' choice didn't apply.
#'
#' \code{cmdstanr} is the actively-maintained Stan interface (3–10× faster
#' compilation, better diagnostic reporting) and is the recommended backend
#' for production / VPS use.
#' @noRd
resolve_brms_backend <- function(backend = c("rstan", "cmdstanr")) {
  backend <- match.arg(backend)
  if (backend == "cmdstanr") {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      stop(
        "backend = 'cmdstanr' requires the cmdstanr package. ",
        "Install it with: ",
        "install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos'))) ",
        "and then a CmdStan toolchain via cmdstanr::install_cmdstan()."
      )
    }
    cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(),
                             error = function(e) NULL)
    if (is.null(cmdstan_path) || !nzchar(cmdstan_path)) {
      stop(
        "backend = 'cmdstanr' is installed but no CmdStan toolchain was ",
        "found. Run cmdstanr::install_cmdstan() once to install it."
      )
    }
  }
  backend
}


#' Empirical coverage of posterior predictive intervals
#'
#' For a fitted brmsfit, simulates from the posterior predictive
#' distribution at the training data points and computes the empirical
#' coverage of nominal 50\%, 80\%, and 95\% intervals against the observed
#' response. Good coverage means \code{empirical / nominal ~ 1}; large
#' deviations indicate model mis-specification.
#'
#' Returns a one-row data frame:
#' \code{cov_50}, \code{cov_80}, \code{cov_95} (empirical coverage as a
#' proportion of training observations falling within the interval), and
#' an \code{n_obs} count. Returns \code{NULL} on any failure so the
#' analysis still completes.
#'
#' Limitation: this is \emph{in-sample} coverage. It catches obvious
#' mis-specification (e.g. residuals far heavier-tailed than the assumed
#' family) but does not substitute for out-of-sample validation
#' (\code{\link{bayes_loo}} for that).
#' @noRd
bayes_ppc_coverage <- function(model, response_var = NULL) {
  if (!requireNamespace("brms", quietly = TRUE) || is.null(model)) return(NULL)
  yrep <- tryCatch(brms::posterior_predict(model),
                   error = function(e) NULL,
                   warning = function(w) NULL)
  if (is.null(yrep) || !is.matrix(yrep)) return(NULL)

  y_obs <- tryCatch({
    if (is.null(response_var)) {
      brms::standata(model)$Y
    } else {
      model$data[[response_var]]
    }
  }, error = function(e) NULL)
  if (is.null(y_obs) || length(y_obs) != ncol(yrep)) return(NULL)

  # Per-observation predictive quantile bounds
  q_bounds <- apply(yrep, 2L, function(col) {
    stats::quantile(col, c(0.025, 0.10, 0.25, 0.75, 0.90, 0.975),
                    names = FALSE, na.rm = TRUE)
  })
  # q_bounds is 6 x n_obs (rows: q025, q10, q25, q75, q90, q975)
  cov50 <- mean(y_obs >= q_bounds[3L, ] & y_obs <= q_bounds[4L, ], na.rm = TRUE)
  cov80 <- mean(y_obs >= q_bounds[2L, ] & y_obs <= q_bounds[5L, ], na.rm = TRUE)
  cov95 <- mean(y_obs >= q_bounds[1L, ] & y_obs <= q_bounds[6L, ], na.rm = TRUE)

  data.frame(
    cov_50  = round(cov50, 4),
    cov_80  = round(cov80, 4),
    cov_95  = round(cov95, 4),
    n_obs   = as.integer(length(y_obs)),
    stringsAsFactors = FALSE
  )
}


#' Posterior probability of direction for emmeans contrasts on a brmsfit
#'
#' Returns a numeric vector of length \code{ncol(contrast_draws)} with the
#' probability that each contrast's effect is in its dominant direction.
#' Reports \code{max(P(effect > 0), P(effect < 0))} — the directional
#' posterior probability that replaces the awkward "does the 95\% CrI exclude
#' zero?" interpretation with a quantitative posterior statement.
#'
#' Implementation pulls posterior draws via \code{emmeans::as.mcmc.emmGrid}
#' (which works for brms emmGrid objects). Falls back to \code{NA_real_}
#' when draws can't be extracted.
#' @noRd
emm_p_direction <- function(emm_or_contrast, n_contrasts = NULL) {
  default <- if (is.null(n_contrasts)) NA_real_ else rep(NA_real_, n_contrasts)
  if (!requireNamespace("emmeans", quietly = TRUE) ||
      is.null(emm_or_contrast)) return(default)
  draws_mat <- tryCatch({
    mc <- emmeans::as.mcmc.emmGrid(emm_or_contrast)
    as.matrix(mc)
  }, error = function(e) NULL)
  if (is.null(draws_mat) || !is.matrix(draws_mat)) return(default)
  vapply(seq_len(ncol(draws_mat)), function(j) {
    x <- draws_mat[, j]
    p_pos <- mean(x > 0, na.rm = TRUE)
    round(max(p_pos, 1 - p_pos), 4)
  }, numeric(1L))
}


#' Posterior Bayesian R^2 summary for a brmsfit
#'
#' Returns the Estimate, Est.Error, and 95\% CrI of \code{brms::bayes_R2}
#' as a one-row data frame. The Bayesian analogue of OLS R^2: a posterior
#' distribution of variance explained, not a point estimate. Returns
#' \code{NULL} when the model can't be evaluated.
#' @noRd
bayes_r2_summary <- function(model) {
  if (!requireNamespace("brms", quietly = TRUE) || is.null(model)) return(NULL)
  r2 <- tryCatch(brms::bayes_R2(model, summary = TRUE),
                 error = function(e) NULL,
                 warning = function(w) NULL)
  if (is.null(r2) || !is.matrix(r2)) return(NULL)
  data.frame(
    Estimate     = round(unname(r2["R2", "Estimate"]),  4),
    Est_Error    = round(unname(r2["R2", "Est.Error"]), 4),
    Lower_95_CrI = round(unname(r2["R2", "Q2.5"]),      4),
    Upper_95_CrI = round(unname(r2["R2", "Q97.5"]),     4),
    stringsAsFactors = FALSE
  )
}


#' PSIS-LOO cross-validation + Pareto-k diagnostics for a brmsfit
#'
#' Returns a one-row data frame with the standard summary plus a per-mouse
#' \code{pareto_k} vector (\code{NA} when LOO failed). The Bayesian
#' counterparts of AIC (\code{elpd_loo}) and Cook's distance
#' (\code{pareto_k > 0.7} flags influential observations). Returns
#' \code{NULL} on any error so the analysis still returns a result list.
#'
#' \itemize{
#'   \item \code{elpd_loo} — expected log pointwise predictive density.
#'         Use for model comparison via \code{loo::loo_compare()}.
#'   \item \code{p_loo} — effective number of parameters; if larger than the
#'         actual parameter count, the model may be mis-specified.
#'   \item \code{looic} — \code{-2 * elpd_loo}; on the AIC/BIC scale.
#'   \item \code{n_high_k} — count of observations with Pareto-k > 0.7.
#'         These are the mice / data points the LOO approximation can't
#'         reliably estimate; investigate them individually.
#'   \item \code{pareto_k} — list-column with the full per-observation
#'         Pareto-k vector.
#' }
#'
#' @noRd
bayes_loo <- function(model) {
  if (!requireNamespace("brms", quietly = TRUE) || is.null(model)) return(NULL)
  loo_obj <- tryCatch(brms::loo(model, save_psis = TRUE),
                      error = function(e) NULL,
                      warning = function(w) NULL)
  if (is.null(loo_obj) || !is.list(loo_obj)) return(NULL)

  ests <- tryCatch(loo_obj$estimates, error = function(e) NULL)
  pk   <- tryCatch(loo_obj$diagnostics$pareto_k, error = function(e) NULL)
  if (is.null(ests) || is.null(pk)) return(NULL)

  data.frame(
    elpd_loo = round(unname(ests["elpd_loo", "Estimate"]), 3),
    se_elpd  = round(unname(ests["elpd_loo", "SE"]),       3),
    p_loo    = round(unname(ests["p_loo",    "Estimate"]), 3),
    looic    = round(unname(ests["looic",    "Estimate"]), 3),
    n_high_k = as.integer(sum(pk > 0.7, na.rm = TRUE)),
    pareto_k = I(list(pk)),
    stringsAsFactors = FALSE
  )
}


#' Resolve cage column, inserting a placeholder when none is supplied
#'
#' Returns a named list with elements \code{df} (possibly modified data frame),
#' \code{cage_column} (resolved column name), and \code{no_cage_mode} (logical).
#' When \code{cage_column} is \code{NULL} or absent from \code{df}, a column
#' named \code{".cage_placeholder"} with all values \code{"1"} is added so that
#' make_mouse_key() and RE formulas always have a valid cage column.
#' @noRd
setup_cage_column <- function(df, cage_column) {
  no_cage_mode <- is.null(cage_column) || !cage_column %in% colnames(df)
  if (no_cage_mode) {
    cage_column       <- ".cage_placeholder"
    df[[cage_column]] <- "1"
  }
  list(df = df, cage_column = cage_column, no_cage_mode = no_cage_mode)
}


# ── Prior vs posterior plot ────────────────────────────────────────────────────

#' Prior vs posterior density overlay for treatment-effect coefficients
#'
#' Used by bayesian_tumor_growth(), bayesian_body_weight(), and
#' bayesian_survival(). Requires \code{sample_prior = "yes"} at fit time.
#' @noRd
bayes_prior_posterior_plot <- function(model, treatment_column) {
  post <- tryCatch(brms::as_draws_df(model), error = function(e) NULL)
  if (is.null(post)) return(NULL)

  safe_tx <- gsub("([.^$*+?()\\[\\]{}|])", "\\\\\\1", treatment_column,
                  perl = TRUE)
  tx_cols <- grep(paste0("^b_", safe_tx), names(post), value = TRUE)
  if (length(tx_cols) == 0) return(NULL)

  prior_col <- if ("prior_b" %in% names(post)) post$prior_b else NULL

  clean <- function(x) sub(paste0("^b_", treatment_column), "", x)

  post_long <- do.call(rbind, lapply(tx_cols, function(col) {
    data.frame(Parameter = clean(col), Value = post[[col]],
               Source = "Posterior", stringsAsFactors = FALSE)
  }))

  plot_df <- if (!is.null(prior_col)) {
    prior_long <- do.call(rbind, lapply(tx_cols, function(col) {
      data.frame(Parameter = clean(col), Value = prior_col,
                 Source = "Prior", stringsAsFactors = FALSE)
    }))
    rbind(post_long, prior_long)
  } else {
    post_long
  }

  plot_df$Source <- factor(plot_df$Source, levels = c("Prior", "Posterior"))

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
    ggplot2::facet_wrap(~ Parameter, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = c(Prior = "grey60", Posterior = "steelblue")
    ) +
    ggplot2::scale_colour_manual(
      values = c(Prior = "grey40", Posterior = "steelblue4")
    ) +
    ggplot2::labs(
      title    = "Prior vs. Posterior Distributions",
      subtitle = paste0(
        "Treatment-effect coefficients; dashed line = no effect (0)"
      ),
      x      = "Coefficient value",
      y      = "Density",
      fill   = NULL,
      colour = NULL
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(legend.position = "top")
}


# ── Posterior summary builder ──────────────────────────────────────────────────

#' Extract and standardise fixed-effects posterior summary from a brmsfit
#'
#' Returns a data frame with columns: Parameter, Estimate, Est.Error,
#' Lower_95_CrI, Upper_95_CrI, Rhat, Bulk_ESS, Tail_ESS.
#' @noRd
build_posterior_summary <- function(model) {
  brms_smry <- summary(model)
  fixed_df  <- as.data.frame(brms_smry$fixed)
  fixed_df  <- cbind(
    Parameter = rownames(fixed_df), fixed_df, stringsAsFactors = FALSE
  )
  rownames(fixed_df) <- NULL
  names(fixed_df)[names(fixed_df) == "l-95% CI"] <- "Lower_95_CrI"
  names(fixed_df)[names(fixed_df) == "u-95% CI"] <- "Upper_95_CrI"
  fixed_df
}
