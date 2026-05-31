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
#' @section Independence assumption — caveat:
#' Because the TG and BW models are fitted independently, any biological
#' correlation between within-animal efficacy and toxicity (an animal that
#' does poorly on both, or unusually well on both) is \strong{not} captured
#' in the joint posterior. The paired-draw construction implicitly assumes
#' the two posteriors are independent. The TWM CrI will be:
#' \itemize{
#'   \item too \emph{wide} when TG and BW responses are positively correlated
#'     within animal (the joint uncertainty would partially cancel);
#'   \item too \emph{narrow} when they are negatively correlated.
#' }
#' A joint multivariate brms model (\code{brms::brm(mvbind(Volume, Weight) ~
#' ...)}) would resolve this. Until then, treat the TWM CrI as a reasonable
#' first-order propagation, not a calibrated joint credible interval.
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
  vols_mat   <- bayes_backtransform(epred_tg,       tg_result$transform_used)
  wts_first  <- bayes_backtransform(epred_bw_first, bw_result$transform_used)
  wts_last   <- bayes_backtransform(epred_bw_last,  bw_result$transform_used)

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
  # Normalise for whitespace/case differences before intersection (B3.4)
  common <- intersect(
    trimws(treated_tg), trimws(treated_bw)
  )
  # Resolve back to original names from the TG model
  common <- treated_tg[trimws(treated_tg) %in% common]

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


#' Fit and Compute Bayesian Therapeutic Window from Raw Data
#'
#' Convenience wrapper that fits a \code{\link{bayesian_tumor_growth}} model
#' and a \code{\link{bayesian_body_weight}} model from raw data frames in a
#' single call, then passes both fitted results to
#' \code{\link{bayesian_therapeutic_window}} to compute the posterior TWM
#' distribution. All arguments that control the fitting of both sub-models are
#' exposed directly; see the individual function documentation for details.
#'
#' @param tg_df Data frame of tumor-growth observations. Must contain
#'   \code{time_column}, \code{treatment_column}, \code{id_column}, and
#'   \code{volume_column}.
#' @param bw_df Data frame of body-weight observations. Must contain
#'   \code{time_column}, \code{treatment_column}, \code{id_column}, and
#'   \code{weight_column}.
#' @param time_column Name of the study-day column present in both data frames.
#'   Default \code{"Day"}.
#' @param treatment_column Name of the treatment-group column present in both
#'   data frames. Default \code{"Treatment"}.
#' @param id_column Name of the animal-ID column present in both data frames.
#'   Default \code{"ID"}.
#' @param cage_column Name of the cage column present in both data frames.
#'   \code{NULL} (default) omits cage random effects.
#' @param reference_group Name of the vehicle/control group. \code{NULL}
#'   (default) auto-detects from factor levels.
#' @param volume_column Name of the tumor-volume column in \code{tg_df}.
#'   Default \code{"Volume"}.
#' @param tg_transform Variance-stabilising transform applied to tumor volume
#'   before fitting. One of \code{"log"} (default), \code{"sqrt"},
#'   \code{"none"}.
#' @param necrotic_column Optional column in \code{tg_df} flagging necrotic
#'   observations. Passed to \code{bayesian_tumor_growth()}.
#' @param weight_column Name of the body-weight column in \code{bw_df}.
#'   Default \code{"Weight"}.
#' @param bw_transform Variance-stabilising transform applied to body weight
#'   before fitting. One of \code{"none"} (default), \code{"log"},
#'   \code{"sqrt"}.
#' @param prior_strength Prior strength used for both sub-models. One of
#'   \code{"skeptical"} (default), \code{"weakly_informative"},
#'   \code{"informative"}, \code{"diffuse"}, \code{"manual"}.
#' @param n_chains Number of MCMC chains. Default \code{4L}.
#' @param n_warmup Warm-up (burn-in) iterations per chain. Default \code{1000L}.
#' @param n_iter Post-warmup draws per chain. Default \code{500L}.
#' @param seed Random seed for reproducibility. Default \code{42L}.
#' @param include_cage_effect Logical. Include cage-level random intercept in
#'   both sub-models when \code{cage_column} is supplied? Default \code{TRUE}.
#' @param noise_floor Minimum absolute weight-loss percentage used as the TWM
#'   denominator (passed to \code{bayesian_therapeutic_window()}). Default
#'   \code{1.0}.
#' @param tg_endpoint_day Study day for TGI computation. \code{NULL} (default)
#'   uses the last observed day in the TG model data.
#' @param plots Logical. Return \pkg{ggplot2} plot objects? Default \code{TRUE}.
#' @param verbose Logical. Print progress messages? Default \code{FALSE}.
#'
#' @return The same named list as \code{\link{bayesian_therapeutic_window}},
#'   with an additional element \code{tg_result} and \code{bw_result} holding
#'   the full outputs of the two sub-model fits.
#'
#' @seealso \code{\link{bayesian_therapeutic_window}},
#'   \code{\link{bayesian_tumor_growth}}, \code{\link{bayesian_body_weight}}
#'
#' @export
bayesian_twm_from_data <- function(
  tg_df,
  bw_df,
  time_column         = "Day",
  treatment_column    = "Treatment",
  id_column           = "ID",
  cage_column         = NULL,
  reference_group     = NULL,
  volume_column       = "Volume",
  tg_transform        = c("log", "sqrt", "none"),
  necrotic_column     = NULL,
  weight_column       = "Weight",
  bw_transform        = c("none", "log", "sqrt"),
  prior_strength      = c("skeptical", "weakly_informative",
                          "informative", "diffuse", "manual"),
  n_chains            = 4L,
  n_warmup            = 1000L,
  n_iter              = 500L,
  seed                = 42L,
  include_cage_effect = TRUE,
  noise_floor         = 1.0,
  tg_endpoint_day     = NULL,
  plots               = TRUE,
  verbose             = FALSE
) {
  tg_transform   <- match.arg(tg_transform)
  bw_transform   <- match.arg(bw_transform)
  prior_strength <- match.arg(prior_strength)

  if (isTRUE(verbose)) message("Fitting Bayesian tumor-growth model...")
  tg_result <- bayesian_tumor_growth(
    df               = tg_df,
    time_column      = time_column,
    volume_column    = volume_column,
    treatment_column = treatment_column,
    id_column        = id_column,
    cage_column      = cage_column,
    transform        = tg_transform,
    reference_group  = reference_group,
    prior_strength   = prior_strength,
    n_chains         = n_chains,
    n_warmup         = n_warmup,
    n_iter           = n_iter,
    seed             = seed,
    include_cage_effect = include_cage_effect,
    return_model     = TRUE,
    plots            = plots,
    verbose          = verbose,
    necrotic_column  = necrotic_column
  )

  if (isTRUE(verbose)) message("Fitting Bayesian body-weight model...")
  bw_result <- bayesian_body_weight(
    df               = bw_df,
    weight_column    = weight_column,
    time_column      = time_column,
    treatment_column = treatment_column,
    id_column        = id_column,
    cage_column      = cage_column,
    transform        = bw_transform,
    reference_group  = reference_group,
    prior_strength   = prior_strength,
    n_chains         = n_chains,
    n_warmup         = n_warmup,
    n_iter           = n_iter,
    seed             = seed,
    include_cage_effect = include_cage_effect,
    return_model     = TRUE,
    plots            = plots,
    verbose          = verbose
  )

  if (isTRUE(verbose)) message("Computing therapeutic window metric...")
  twm_result <- bayesian_therapeutic_window(
    tg_result        = tg_result,
    bw_result        = bw_result,
    treatment_column = treatment_column,
    time_column      = time_column,
    reference_group  = reference_group,
    noise_floor      = noise_floor,
    tg_endpoint_day  = tg_endpoint_day,
    plots            = plots,
    verbose          = verbose
  )

  c(twm_result, list(tg_result = tg_result, bw_result = bw_result))
}
