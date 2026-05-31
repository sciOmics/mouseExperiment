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
    diffuse            = list(b_sd = 2.5,  exp_rate = 0.5)
  )
}


# ── Cage column setup ──────────────────────────────────────────────────────────

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

  safe_tx <- gsub("([.^$*+?()\\[\\]{}|])", "\\\\\\1", treatment_column)
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
