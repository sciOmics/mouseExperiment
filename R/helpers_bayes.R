#' Bayesian-fit display helpers (CODE_REVIEW.md C.1/C.3).
#'
#' Split out of `helpers.R` (was 1,417 LOC) so the Bayesian-diagnostics
#' surface — used by every Bayesian module — lives in its own file rather
#' than buried among 38 unrelated utility functions.

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
#'     with Pareto-k > 0.7 (influential / unreliable LOO).
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
      "elpd_loo = %.2f (SE %.2f); %d observation(s) with Pareto-k > 0.7",
      loo$elpd_loo[1], loo$se_elpd[1], as.integer(loo$n_high_k[1])
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
#' @param threshold Pareto-k cutoff (default 0.7).
#' @return A data frame with columns \code{Obs_Index} and \code{Pareto_k},
#'   sorted descending. Returns an empty data frame when no influential
#'   observations exist or the LOO field is missing.
#' @keywords internal
bayes_influential_obs <- function(result, threshold = 0.7) {
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
