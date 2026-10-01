#' Bayesian-fit display helpers (CODE_REVIEW.md C.1/C.3).
#'
#' Split out of `helpers.R` (was 1,417 LOC) so the Bayesian-diagnostics
#' surface — used by every Bayesian module — lives in its own file rather
#' than buried among 38 unrelated utility functions.

# ---- Elapsed-time formatter -------------------------------------------------

#' Human-readable elapsed time string.
#'
#' Used in the Summary tab of each Bayesian module to show how long the
#' fit took. attach_bayes_poller stamps `elapsed_seconds` onto every
#' completed-Bayesian result based on the job registry's started_at /
#' completed_at.
#'
#' @param seconds Numeric seconds. NA / NULL returns NULL.
#' @return Character string like "47 s", "3 min 12 s", or NULL.
#' @keywords internal
format_elapsed <- function(seconds) {
  if (is.null(seconds) || is.na(seconds) || !is.finite(seconds) ||
      seconds < 0) return(NULL)
  if (seconds < 60) return(sprintf("%.0f s", seconds))
  mins <- floor(seconds / 60)
  secs <- round(seconds - mins * 60)
  if (mins < 60) return(sprintf("%d min %d s", mins, secs))
  hrs  <- floor(mins / 60)
  rem  <- mins - hrs * 60
  sprintf("%d h %d min", hrs, rem)
}

# ---- Bayesian diagnostics panel (CODE_REVIEW.md Round 2 D-series) -----------

#' Standard Bayesian diagnostics panel (model-fit + sampler banner)
#'
#' Renders a single block summarising every Bayesian-fit diagnostic added
#' to mouseExperiment v0.4.5 / v0.4.6:
#' \itemize{
#'   \item NUTS sampler banner — divergent transitions, max-treedepth hits,
#'     and minimum E-BFMI. Surfaces as a warning alert when any pathology
#'     is detected, success alert otherwise.
#'   \item Bayes R^2 — posterior median + 95\% credible interval.
#'   \item PPC coverage — empirical 50 / 80 / 95 \% interval coverage.
#'   \item PSIS-LOO — \code{elpd_loo} ± SE and the count of observations
#'     with Pareto-k above the backend's threshold (influential / unreliable LOO).
#' }
#' Fields that are missing from the result list are silently skipped, so
#' the helper works for both v0.4.5 results (no PPC coverage) and v0.4.6+.
#'
#' @param result A Bayesian-analysis result list (e.g. from
#'   \code{bayesian_tumor_growth()}).
#' @return A \code{shiny::tagList()} ready to render in a \code{uiOutput}.
#' @keywords internal
bayes_diagnostics_panel <- function(result) {
  if (is.null(result)) return(NULL)
  bayes_R2 <- result$bayes_R2
  ppc      <- result$ppc_coverage
  nuts     <- result$nuts_diagnostics
  loo      <- result$loo_diagnostics

  parts <- list()

  # NUTS sampler banner
  if (!is.null(nuts) && nrow(nuts) > 0L) {
    n_div  <- as.integer(nuts$n_divergent[1])
    n_td   <- as.integer(nuts$n_max_treedepth[1])
    ebfmi  <- as.numeric(nuts$ebfmi_min[1])
    issues <- character(0)
    if (isTRUE(n_div > 0L)) {
      issues <- c(issues, sprintf("%d divergent transition(s)", n_div))
    }
    if (isTRUE(n_td > 0L)) {
      issues <- c(issues, sprintf("%d iteration(s) hit max_treedepth", n_td))
    }
    if (isTRUE(is.finite(ebfmi) && ebfmi < 0.3)) {
      issues <- c(issues, sprintf("E-BFMI low (%.2f)", ebfmi))
    }
    if (length(issues) > 0L) {
      parts$nuts <- shiny::div(
        class = "alert alert-warning",
        shiny::icon("triangle-exclamation"),
        shiny::strong(" Sampler diagnostic warnings: "),
        paste(issues, collapse = "; "), ".",
        shiny::br(),
        shiny::tags$small(
          "Posterior estimates may be unreliable for affected parameters. ",
          "Investigate (e.g. reparameterise, raise adapt_delta, ",
          "or increase max_treedepth) before reporting."
        )
      )
    } else {
      parts$nuts <- shiny::div(
        class = "alert alert-success",
        shiny::icon("check-circle"),
        shiny::strong(" NUTS sampler diagnostics clean: "),
        sprintf(
          "%d divergence(s), %d max-treedepth hit(s), E-BFMI = %.2f.",
          n_div, n_td, ebfmi
        )
      )
    }
  }

  # Stats panel
  rows <- list()
  if (!is.null(bayes_R2) && nrow(bayes_R2) > 0L) {
    rows$bayes_R2 <- sprintf(
      "%.3f (95%% CrI: %.3f – %.3f)",
      bayes_R2$Estimate[1], bayes_R2$Lower_95_CrI[1], bayes_R2$Upper_95_CrI[1]
    )
  }
  if (!is.null(ppc) && nrow(ppc) > 0L) {
    rows$ppc <- sprintf(
      "50%%: %d%% | 80%%: %d%% | 95%%: %d%% (n_obs = %d)",
      round(100 * ppc$cov_50[1]), round(100 * ppc$cov_80[1]),
      round(100 * ppc$cov_95[1]), as.integer(ppc$n_obs[1])
    )
  }
  if (!is.null(loo) && nrow(loo) > 0L) {
    rows$loo <- sprintf(
      "elpd_loo = %.2f (SE %.2f); %d observation(s) with Pareto-k > %s",
      loo$elpd_loo[1], loo$se_elpd[1], as.integer(loo$n_high_k[1]),
      format(pareto_k_threshold(list(loo_diagnostics = loo)))
    )
  }

  if (length(rows) > 0L) {
    dl_items <- list()
    if (!is.null(rows$bayes_R2)) {
      dl_items <- c(dl_items,
        list(shiny::tags$dt("Bayes R²"), shiny::tags$dd(rows$bayes_R2)))
    }
    if (!is.null(rows$ppc)) {
      dl_items <- c(dl_items,
        list(shiny::tags$dt("PPC coverage"), shiny::tags$dd(rows$ppc)))
    }
    if (!is.null(rows$loo)) {
      dl_items <- c(dl_items,
        list(shiny::tags$dt("PSIS-LOO"), shiny::tags$dd(rows$loo)))
    }
    parts$stats <- shiny::div(
      class = "panel panel-default",
      shiny::div(class = "panel-heading",
        shiny::h5(class = "panel-title", "Model fit & diagnostics")),
      shiny::div(class = "panel-body",
        do.call(shiny::tags$dl,
                c(class = "dl-horizontal", dl_items))
      )
    )
  }

  if (length(parts) == 0L) return(NULL)
  do.call(shiny::tagList, unname(parts))
}

#' List of influential observations from a Bayesian fit's PSIS-LOO output
#'
#' Filters \code{result$loo_diagnostics$pareto_k} (a list-column carrying the
#' per-observation Pareto-k vector) for values > the threshold and returns a
#' small data frame with observation index + Pareto-k value, intended for
#' display in a DT table on the dashboard's MCMC Diagnostics tab.
#'
#' @param result A Bayesian-analysis result list.
#' @param threshold Pareto-k cutoff. Default: the backend's
#'   \code{k_threshold}, which depends on the number of draws (0.667 at
#'   1,000; backend R20.35), or 0.7 for an older result.
#' @return A data frame with columns \code{Obs_Index} and \code{Pareto_k},
#'   sorted descending. Returns an empty data frame when no influential
#'   observations exist or the LOO field is missing.
#' @keywords internal
bayes_influential_obs <- function(result, threshold = pareto_k_threshold(result)) {
  loo <- result$loo_diagnostics
  if (is.null(loo) || nrow(loo) == 0L) {
    return(data.frame(Obs_Index = integer(0), Pareto_k = numeric(0)))
  }
  pk <- loo$pareto_k[[1]]
  if (is.null(pk)) {
    return(data.frame(Obs_Index = integer(0), Pareto_k = numeric(0)))
  }
  high <- which(pk > threshold)
  if (length(high) == 0L) {
    return(data.frame(Obs_Index = integer(0), Pareto_k = numeric(0)))
  }
  df <- data.frame(
    Obs_Index = as.integer(high),
    Pareto_k  = round(as.numeric(pk[high]), 4),
    stringsAsFactors = FALSE
  )
  df[order(-df$Pareto_k), , drop = FALSE]
}

#' Build a posterior-predictive coverage bar chart.
#'
#' Plots empirical vs nominal coverage at 50 / 80 / 95 %. Bars near the
#' dashed nominal line = well-calibrated; below = under-coverage (too
#' narrow); above = over-coverage (too wide).
#'
#' Closes DIAGNOSTICS.md gap (10) — turning the dl row "50%: 48% | 80%:
#' 78% | 95%: 94%" into a quick-scan visual.
#'
#' @param result A Bayesian-analysis result list (with `ppc_coverage`).
#' @return ggplot, or NULL if `ppc_coverage` isn't populated.
#' @keywords internal
build_ppc_coverage_plot <- function(result) {
  if (is.null(result)) return(NULL)
  ppc <- result$ppc_coverage
  if (is.null(ppc) || nrow(ppc) == 0L) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  cov_df <- data.frame(
    Nominal   = c(50, 80, 95),
    Empirical = c(as.numeric(ppc$cov_50[1]),
                  as.numeric(ppc$cov_80[1]),
                  as.numeric(ppc$cov_95[1])) * 100
  )
  ggplot2::ggplot(cov_df, ggplot2::aes(x = factor(.data[["Nominal"]]),
                                       y = .data[["Empirical"]])) +
    ggplot2::geom_col(fill = "#2c7bb6", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%",
                                                   .data[["Empirical"]])),
                       vjust = -0.4, size = 4) +
    ggplot2::geom_segment(
      data = data.frame(x = c(0.5, 1.5, 2.5), xend = c(1.5, 2.5, 3.5),
                        y = c(50, 80, 95)),
      ggplot2::aes(x = .data[["x"]], xend = .data[["xend"]],
                   y = .data[["y"]],  yend = .data[["y"]]),
      colour = "red", linetype = "dashed", inherit.aes = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 105)) +
    ggplot2::labs(title = "Posterior predictive coverage",
                  subtitle = "Red dashed = nominal target; bars = empirical",
                  x = "Nominal interval (%)",
                  y = "Empirical coverage (%)") +
    ggplot2::theme_classic()
}

#' Build a Pareto-k histogram with the backend's cut-off line.
#'
#' Closes DIAGNOSTICS.md gap (11). The dashboard already shows the table
#' of high-k observations; the histogram gives the full distribution so
#' users can see whether 5 high-k observations have nearest neighbours at
#' k = 0.65 or k = 0.3.
#'
#' @param result A Bayesian-analysis result list.
#' @return ggplot, or NULL if Pareto-k isn't populated.
#' @keywords internal
build_pareto_k_plot <- function(result) {
  if (is.null(result)) return(NULL)
  loo <- result$loo_diagnostics
  if (is.null(loo) || nrow(loo) == 0L) return(NULL)
  pk <- loo$pareto_k[[1]]
  if (is.null(pk) || length(pk) == 0L) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  k_thr <- pareto_k_threshold(result)
  df <- data.frame(k = as.numeric(pk))
  ggplot2::ggplot(df, ggplot2::aes(x = .data[["k"]])) +
    ggplot2::geom_histogram(bins = 30L, fill = "#2c7bb6",
                            colour = "white") +
    ggplot2::geom_vline(xintercept = k_thr, colour = "red",
                        linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0.5, colour = "orange",
                        linetype = "dotted", linewidth = 0.6) +
    ggplot2::annotate("text", x = k_thr, y = Inf, label = sprintf("  k = %s", format(k_thr)),
                      colour = "red", hjust = 0, vjust = 1.4, size = 3.5) +
    ggplot2::annotate("text", x = 0.5, y = Inf, label = "  k = 0.5",
                      colour = "orange", hjust = 0, vjust = 1.4, size = 3.5) +
    ggplot2::labs(title = "PSIS-LOO Pareto-k distribution",
                  subtitle = sprintf("n = %d; cut-off at %s (red; depends on the number of draws), 0.5 (orange)",
                                     length(pk), format(k_thr)),
                  x = "Pareto k", y = "Count of observations") +
    ggplot2::theme_classic()
}

#' The Pareto-k threshold of a result
#'
#' Backend 0.28.0 (R20.35) reports `k_threshold`, which depends on the number
#' of draws; an older result falls back to 0.7.
#' @keywords internal
pareto_k_threshold <- function(result) {
  k <- result$loo_diagnostics$k_threshold
  if (length(k) >= 1L && is.finite(k[1])) as.numeric(k[1]) else 0.7
}

#' Build an MCMC rank-overlay plot.
#'
#' Modern convergence diagnostic that replaces Rhat for finer-grained
#' detection. Each chain's draws are rank-normalised and plotted as a
#' histogram per chain; chains that have converged should look like
#' uniform draws (flat histograms) and resemble each other.
#'
#' Closes DIAGNOSTICS.md gap (12).
#'
#' @param brmsfit A fitted brms model (or any model that
#'   `bayesplot::mcmc_rank_overlay` accepts).
#' @param pars Optional character vector of parameter names to show.
#'   Default uses brms's `pars = NA` which shows all non-`lp__`
#'   parameters.
#' @return ggplot, or NULL on failure.
#' @keywords internal
build_mcmc_rank_plot <- function(draws, pars = NULL) {
  if (is.null(draws)) return(NULL)
  if (!requireNamespace("bayesplot", quietly = TRUE)) return(NULL)
  if (!requireNamespace("posterior", quietly = TRUE)) return(NULL)
  tryCatch({
    # R20.D18: the result's headline draws (backend >= 0.28.0), or a fit.
    draws <- posterior::as_draws_array(draws)
    # Drop lp__ which adds noise; keep treatment-relevant params.
    if (is.null(pars)) {
      keep <- dimnames(draws)$variable
      keep <- keep[!grepl("^lp__$", keep)]
      # Limit to first 12 to keep the plot legible
      if (length(keep) > 12L) keep <- keep[seq_len(12L)]
      draws <- draws[, , keep, drop = FALSE]
    }
    bayesplot::mcmc_rank_overlay(draws)
  }, error = function(e) NULL)
}
