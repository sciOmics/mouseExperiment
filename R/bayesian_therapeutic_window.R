# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Bayesian Therapeutic Window Metric
#'
#' Propagates full posterior uncertainty through the Therapeutic Window
#' Metric (TWM = TGI / |weight-loss %|) by performing draw-wise arithmetic
#' on posterior predictive samples from a \code{\link{bayesian_tumor_growth}}
#' result and a \code{\link{bayesian_body_weight}} result. The resulting
#' posterior distribution of TWM per treatment group has a direct probability
#' interpretation: e.g. "there is 95 % posterior probability that TWM lies
#' between X and Y."
#'
#' @param tg_result A list returned by \code{\link{bayesian_tumor_growth}} with
#'   \code{return_model = TRUE}. The \code{brmsfit} object stored at
#'   \code{tg_result$model} is used to generate posterior volume predictions.
#' @param bw_result A list returned by \code{\link{bayesian_body_weight}} with
#'   \code{return_model = TRUE}. The \code{brmsfit} object stored at
#'   \code{bw_result$model} is used to generate posterior weight predictions.
#' @param treatment_column Name of the treatment-group column that appears in
#'   both fitted models. Default \code{"Treatment"}.
#' @param time_column Name of the study-day column in both models. Default
#'   \code{"Day"}.
#' @param reference_group Name of the vehicle/control group. Auto-detected
#'   from the factor levels of the TG model if \code{NULL}.
#' @param noise_floor Minimum absolute weight-loss percentage used as the TWM
#'   denominator. When \eqn{|WL\%| \leq} \code{noise_floor}, the denominator
#'   is clamped to \code{noise_floor} so that TWM = TGI / (noise_floor / 100).
#'   Default \code{1.0} (i.e. 1 %).
#' @param tg_endpoint_day Study day for TGI computation. \code{NULL} (default)
#'   uses the last observed day in the TG model data.
#' @param plots Logical. Return \pkg{ggplot2} plot objects? Default \code{TRUE}.
#' @param verbose Logical. Print progress messages? Default \code{FALSE}.
#'
#' @details
#' The joint posterior TWM for group \eqn{g} at draw \eqn{d} is:
#' \deqn{
#'   \mathrm{TWM}_g^{(d)} =
#'   \frac{\mathrm{TGI}_g^{(d)}}{
#'     \max(|WL_g^{(d)}\%|, \text{noise\_floor}) / 100}
#' }
#' where
#' \eqn{
#'   \mathrm{TGI}_g^{(d)} =
#'   1 - \hat{V}_g^{(d)} / \hat{V}_{\text{ctrl}}^{(d)}
#' }
#' (population-level posterior predictive volume at the endpoint day,
#' back-transformed from the model scale), and
#' \eqn{WL_g^{(d)}} is the posterior predictive percentage weight change from
#' the first to the last study day in the BW model.
#'
#' Draws from the two independently fitted models are treated as independent
#' samples from their respective marginal posteriors and paired draw-by-draw
#' (using the minimum number of available draws if the draw counts differ).
#'
#' @return A named list:
#' \describe{
#'   \item{\code{model_type_used}}{Character \code{"bayes_twm"}.}
#'   \item{\code{twm_table}}{Data frame with one row per treated group:
#'     \code{Group}, \code{TWM_Median}, \code{TWM_Lower} (2.5 % CrI),
#'     \code{TWM_Upper} (97.5 % CrI), \code{TGI_Median}, \code{WL_Pct_Median},
#'     \code{N_Draws}.}
#'   \item{\code{tgi_summary}}{Data frame of per-group TGI (including the
#'     reference): \code{Group}, \code{TGI_Median}, \code{TGI_Lower},
#'     \code{TGI_Upper}.}
#'   \item{\code{wl_summary}}{Data frame of per-group weight-loss %:
#'     \code{Group}, \code{WL_Median}, \code{WL_Lower}, \code{WL_Upper}.}
#'   \item{\code{summary}}{Named list of metadata.}
#'   \item{\code{twm_plot}}{Forest plot of TWM with 95 % CrI; vertical dashed
#'     line at TWM = 1 (break-even). \code{NULL} when \code{plots = FALSE}.}
#'   \item{\code{tgi_wl_plot}}{Scatter of posterior-median TGI vs weight-loss %
#'     per group with TWM = 1 isoline. \code{NULL} when
#'     \code{plots = FALSE}.}
#' }
#'
#' @importFrom stats median quantile
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_vline geom_point
#'   geom_text geom_abline labs theme theme_classic
#' @export
bayesian_therapeutic_window <- function(
  tg_result,
  bw_result,
  treatment_column = "Treatment",
  time_column      = "Day",
  reference_group  = NULL,
  noise_floor      = 1.0,
  tg_endpoint_day  = NULL,
  plots            = TRUE,
  verbose          = FALSE
) {

  # ── Validate inputs ────────────────────────────────────────────────────────
  if (!inherits(tg_result$model, "brmsfit")) {
    stop(
      "'tg_result$model' must be a brmsfit object. ",
      "Re-run bayesian_tumor_growth() with return_model = TRUE."
    )
  }
  if (!inherits(bw_result$model, "brmsfit")) {
    stop(
      "'bw_result$model' must be a brmsfit object. ",
      "Re-run bayesian_body_weight() with return_model = TRUE."
    )
  }

  # ── Extract model data ─────────────────────────────────────────────────────
  tg_data <- tg_result$model$data
  bw_data <- bw_result$model$data

  if (!treatment_column %in% colnames(tg_data)) {
    stop(
      "'", treatment_column, "' not found in TG model data. ",
      "Check the treatment_column argument."
    )
  }
  if (!treatment_column %in% colnames(bw_data)) {
    stop(
      "'", treatment_column, "' not found in BW model data. ",
      "Check the treatment_column argument."
    )
  }

  tg_groups <- levels(tg_data[[treatment_column]])
  bw_groups <- levels(bw_data[[treatment_column]])

  # ── Reference group ────────────────────────────────────────────────────────
  if (is.null(reference_group)) {
    reference_group <- tg_groups[1L]
  } else if (!reference_group %in% tg_groups) {
    stop(
      "Reference group '", reference_group,
      "' not found in TG model treatment levels: ",
      paste(tg_groups, collapse = ", ")
    )
  }

  # ── Endpoint and study-day ranges ──────────────────────────────────────────
  tg_days <- sort(unique(as.numeric(tg_data[[time_column]])))
  bw_days <- sort(unique(as.numeric(bw_data[[time_column]])))

  ep_day     <- if (!is.null(tg_endpoint_day)) {
    as.numeric(tg_endpoint_day)
  } else {
    tg_days[length(tg_days)]
  }
  bw_first   <- bw_days[1L]
  bw_last    <- bw_days[length(bw_days)]

  if (isTRUE(verbose)) {
    message(
      "Computing posterior TWM: TG endpoint = day ", ep_day,
      "; BW change = day ", bw_first, " to day ", bw_last
    )
  }

  # ── Helper: make newdata for population-level predictions ──────────────────
  make_nd <- function(groups, day) {
    nd <- data.frame(
      g = groups, d = day, stringsAsFactors = FALSE
    )
    names(nd) <- c(treatment_column, time_column)
    nd
  }

  # ── Posterior predictive draws (population-level) ──────────────────────────
  # TG model: volume at endpoint day per group
  epred_tg <- tryCatch(
    brms::posterior_epred(
      tg_result$model,
      newdata       = make_nd(tg_groups, ep_day),
      re_formula    = NA,
      allow_new_levels = TRUE
    ),
    error = function(e) {
      stop("posterior_epred() failed on TG model: ", conditionMessage(e))
    }
  )

  # BW model: weight at first and last BW day per group
  epred_bw_first <- tryCatch(
    brms::posterior_epred(
      bw_result$model,
      newdata       = make_nd(bw_groups, bw_first),
      re_formula    = NA,
      allow_new_levels = TRUE
    ),
    error = function(e) {
      stop("posterior_epred() failed on BW model (first day): ",
           conditionMessage(e))
    }
  )
  epred_bw_last <- tryCatch(
    brms::posterior_epred(
      bw_result$model,
      newdata       = make_nd(bw_groups, bw_last),
      re_formula    = NA,
      allow_new_levels = TRUE
    ),
    error = function(e) {
      stop("posterior_epred() failed on BW model (last day): ",
           conditionMessage(e))
    }
  )

  # ── Back-transform to original scale ──────────────────────────────────────
  bt <- function(x, transform) {
    if (transform == "log")  return(exp(x))
    if (transform == "sqrt") return(x ^ 2)
    x
  }

  vols_mat   <- bt(epred_tg,        tg_result$transform_used)
  wts_first  <- bt(epred_bw_first,  bw_result$transform_used)
  wts_last   <- bt(epred_bw_last,   bw_result$transform_used)

  # ── TGI matrix [n_draws × n_tg_groups] ────────────────────────────────────
  ref_col_tg <- which(tg_groups == reference_group)
  ctrl_vol   <- vols_mat[, ref_col_tg]
  tgi_mat    <- 1 - sweep(vols_mat, 1L, ctrl_vol, FUN = "/")
  colnames(tgi_mat) <- tg_groups

  # ── Weight-loss % matrix [n_draws × n_bw_groups] ──────────────────────────
  wl_mat <- (wts_last - wts_first) / wts_first * 100
  colnames(wl_mat) <- bw_groups

  # ── Align draws and compute TWM for treated groups ────────────────────────
  treated_tg <- tg_groups[tg_groups != reference_group]
  treated_bw <- bw_groups[bw_groups != reference_group]
  common     <- intersect(treated_tg, treated_bw)

  if (length(common) == 0L) {
    stop(
      "No common treated groups between TG and BW models after excluding '",
      reference_group, "'. Check that both models were fitted on the same ",
      "treatment groups."
    )
  }

  n_draws <- min(nrow(tgi_mat), nrow(wl_mat))

  twm_rows <- lapply(common, function(g) {
    tgi_draws <- tgi_mat[seq_len(n_draws), g]
    wl_draws  <- wl_mat[seq_len(n_draws), g]
    denom     <- pmax(abs(wl_draws), noise_floor) / 100
    twm_draws <- tgi_draws / denom

    tgi_q <- stats::quantile(tgi_draws, c(0.025, 0.5, 0.975))
    twm_q <- stats::quantile(twm_draws, c(0.025, 0.5, 0.975))

    data.frame(
      Group         = g,
      TWM_Median    = round(twm_q["50%"],   3),
      TWM_Lower     = round(twm_q["2.5%"],  3),
      TWM_Upper     = round(twm_q["97.5%"], 3),
      TGI_Median    = round(tgi_q["50%"],   3),
      WL_Pct_Median = round(
        stats::median(wl_draws), 1
      ),
      N_Draws       = n_draws,
      stringsAsFactors = FALSE
    )
  })
  twm_table <- do.call(rbind, twm_rows)
  rownames(twm_table) <- NULL

  # ── TGI summary (all groups including reference) ───────────────────────────
  tgi_summary <- do.call(rbind, lapply(tg_groups, function(g) {
    q <- stats::quantile(tgi_mat[, g], c(0.025, 0.5, 0.975))
    data.frame(
      Group     = g,
      TGI_Median = round(q["50%"],   3),
      TGI_Lower  = round(q["2.5%"],  3),
      TGI_Upper  = round(q["97.5%"], 3),
      stringsAsFactors = FALSE
    )
  }))
  rownames(tgi_summary) <- NULL

  # ── Weight-loss summary (all BW groups) ───────────────────────────────────
  wl_summary <- do.call(rbind, lapply(bw_groups, function(g) {
    q <- stats::quantile(wl_mat[, g], c(0.025, 0.5, 0.975))
    data.frame(
      Group      = g,
      WL_Median  = round(q["50%"],   1),
      WL_Lower   = round(q["2.5%"],  1),
      WL_Upper   = round(q["97.5%"], 1),
      stringsAsFactors = FALSE
    )
  }))
  rownames(wl_summary) <- NULL

  # ── Plots ──────────────────────────────────────────────────────────────────
  twm_plot    <- NULL
  tgi_wl_plot <- NULL

  if (isTRUE(plots) && nrow(twm_table) > 0L) {

    te       <- twm_table
    te$Group <- factor(te$Group, levels = rev(te$Group))

    twm_plot <- ggplot2::ggplot(
      te,
      ggplot2::aes(
        x    = .data[["TWM_Median"]],
        y    = .data[["Group"]],
        xmin = .data[["TWM_Lower"]],
        xmax = .data[["TWM_Upper"]]
      )
    ) +
      ggplot2::geom_pointrange(
        colour = "steelblue", size = 0.8, linewidth = 0.8
      ) +
      ggplot2::geom_vline(
        xintercept = 1, linetype = "dashed",
        colour = "grey50", linewidth = 0.5
      ) +
      ggplot2::labs(
        title    = "Bayesian Therapeutic Window Metric",
        subtitle = paste0(
          "Posterior median ± 95 % CrI; TWM > 1 indicates favourable ",
          "efficacy-safety balance"
        ),
        x = "TWM = TGI / |Weight Loss %|",
        y = NULL
      ) +
      ggplot2::theme_classic(base_size = 14)

    # Scatter: TGI vs weight loss with TWM=1 isoline
    scatter_df <- merge(
      data.frame(
        Group      = twm_table$Group,
        TGI_Median = twm_table$TGI_Median,
        stringsAsFactors = FALSE
      ),
      data.frame(
        Group     = twm_table$Group,
        WL_Median = twm_table$WL_Pct_Median,
        stringsAsFactors = FALSE
      ),
      by = "Group"
    )

    tgi_wl_plot <- tryCatch({
      ggplot2::ggplot(
        scatter_df,
        ggplot2::aes(
          x      = .data[["WL_Median"]],
          y      = .data[["TGI_Median"]],
          label  = .data[["Group"]]
        )
      ) +
        ggplot2::geom_point(
          colour = "steelblue", size = 3
        ) +
        ggplot2::geom_text(
          vjust = -0.8, size = 3.5, colour = "grey30"
        ) +
        ggplot2::geom_vline(
          xintercept = -noise_floor, linetype = "dashed",
          colour = "grey60", linewidth = 0.5
        ) +
        ggplot2::geom_abline(
          slope     = -1 / 100,
          intercept = noise_floor / 100,
          linetype  = "dotted",
          colour    = "tomato",
          linewidth = 0.7
        ) +
        ggplot2::labs(
          title    = "Efficacy vs. Safety (Posterior Medians)",
          subtitle = "Dotted red line = TWM 1 isoline; points above = TWM > 1",
          x        = "Mean weight change (%)",
          y        = "TGI (fraction)"
        ) +
        ggplot2::theme_classic(base_size = 14)
    }, error = function(e) NULL)
  }

  # ── Summary metadata ───────────────────────────────────────────────────────
  analysis_summary <- list(
    analysis_type = paste0(
      "Bayesian Therapeutic Window Metric — draw-wise posterior arithmetic ",
      "(TWM = TGI / |WL%| per MCMC draw)"
    ),
    tg_model = list(
      transform   = tg_result$transform_used,
      endpoint_day = ep_day,
      groups      = tg_groups
    ),
    bw_model = list(
      transform   = bw_result$transform_used,
      first_day   = bw_first,
      last_day    = bw_last,
      groups      = bw_groups
    ),
    parameters = list(
      reference_group  = reference_group,
      noise_floor      = noise_floor,
      n_draws_used     = n_draws,
      common_groups    = common
    )
  )

  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    model_type_used = "bayes_twm",
    twm_table       = twm_table,
    tgi_summary     = tgi_summary,
    wl_summary      = wl_summary,
    summary         = analysis_summary,
    twm_plot        = twm_plot,
    tgi_wl_plot     = tgi_wl_plot
  )
}
