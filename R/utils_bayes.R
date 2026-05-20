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
#' @noRd
make_mcmc_diagnostics <- function(posterior_summary_df) {
  data.frame(
    Parameter = posterior_summary_df$Parameter,
    Rhat      = round(posterior_summary_df$Rhat,     4),
    Bulk_ESS  = round(posterior_summary_df$Bulk_ESS, 0),
    Tail_ESS  = round(posterior_summary_df$Tail_ESS, 0),
    Converged = posterior_summary_df$Rhat <= 1.01,
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
