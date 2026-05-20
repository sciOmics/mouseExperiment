# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

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
#' @param n_iter Total iterations per chain (including warmup). Default
#'   \code{2000}.
#' @param seed Integer random seed. Default \code{42}.
#' @param include_cage_effect Logical. Include cage random intercept when
#'   \code{cage_column} is supplied? Default \code{TRUE}.
#' @param return_model Logical. Return the \code{brmsfit} object? Default
#'   \code{TRUE}.
#' @param plots Logical. Return \pkg{ggplot2} plot objects? Default
#'   \code{TRUE}.
#' @param verbose Logical. Print progress messages? Default \code{FALSE}.
#'
#' @details
#' The model formula is \code{Volume_transformed ~ Treatment * Day + re},
#' where \code{re} is the random-effects term for animal-level intercepts (or
#' slopes). A single model is fitted on all four groups simultaneously so
#' group-specific uncertainty is fully propagated.
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
  n_iter                      = 2000L,
  seed                        = 42L,
  include_cage_effect         = TRUE,
  return_model                = TRUE,
  plots                       = TRUE,
  verbose                     = FALSE
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

  # ── Apply transform ───────────────────────────────────────────────────────
  analysis_df$Response <- switch(
    transform,
    log  = log(analysis_df[[volume_column]]),
    sqrt = sqrt(analysis_df[[volume_column]]),
    analysis_df[[volume_column]]
  )

  # ── Cage placeholder (keeps model data consistent with bayesian_tumor_growth)
  has_cage <- !is.null(cage_column) && cage_column %in% colnames(analysis_df)
  if (has_cage) {
    analysis_df$cage_var <- analysis_df[[cage_column]]
  } else {
    analysis_df$cage_var <- analysis_df[[id_column]]
  }
  analysis_df$id_var <- paste0(
    analysis_df$cage_var, "__", analysis_df[[id_column]]
  )

  # ── Random-effects formula ─────────────────────────────────────────────────
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

  # ── Priors ────────────────────────────────────────────────────────────────
  prior_params <- switch(prior_str,
    skeptical          = list(b = 0.25, sd_rate = 2, sigma_rate = 2),
    weakly_informative = list(b = 1.0,  sd_rate = 1, sigma_rate = 1),
    informative        = list(b = 0.5,  sd_rate = 2, sigma_rate = 2),
    diffuse            = list(b = 2.5,  sd_rate = 0.5, sigma_rate = 0.5),
    manual             = NULL
  )

  if (prior_str == "manual") {
    if (is.null(prior_b) || is.null(prior_intercept) ||
          is.null(prior_sd) || is.null(prior_sigma)) {
      stop(
        "prior_strength = 'manual' requires prior_b, prior_intercept, ",
        "prior_sd, and prior_sigma."
      )
    }
    priors <- c(
      brms::set_prior(prior_b,         class = "b"),
      brms::set_prior(prior_intercept, class = "Intercept"),
      brms::set_prior(prior_sd,        class = "sd"),
      brms::set_prior(prior_sigma,     class = "sigma")
    )
  } else {
    pp <- prior_params
    priors <- c(
      brms::set_prior(
        paste0("normal(0, ", pp$b, ")"), class = "b"
      ),
      brms::set_prior(
        paste0("normal(0, ", pp$b * 2, ")"), class = "Intercept"
      ),
      brms::set_prior(
        paste0("exponential(", pp$sd_rate, ")"), class = "sd"
      ),
      brms::set_prior(
        paste0("exponential(", pp$sigma_rate, ")"), class = "sigma"
      )
    )
  }

  # ── Fit model ─────────────────────────────────────────────────────────────
  if (isTRUE(verbose)) message("Fitting brms model …")

  model <- suppressWarnings(
    brms::brm(
      formula      = brms_formula,
      data         = analysis_df,
      prior        = priors,
      chains       = n_chains,
      iter         = n_iter,
      seed         = seed,
      sample_prior = "yes",
      silent       = if (isTRUE(verbose)) 0L else 2L,
      refresh      = if (isTRUE(verbose)) 100L else 0L
    )
  )

  # ── MCMC diagnostics ──────────────────────────────────────────────────────
  post_sum <- brms::posterior_summary(model)
  post_sum_df <- as.data.frame(post_sum)
  post_sum_df$Parameter <- rownames(post_sum_df)
  rownames(post_sum_df) <- NULL

  rhat_vals <- tryCatch(
    brms::rhat(model), error = function(e) NULL
  )
  if (!is.null(rhat_vals)) {
    diag_df <- data.frame(
      Parameter = names(rhat_vals),
      Rhat      = as.numeric(rhat_vals),
      stringsAsFactors = FALSE
    )
    ess_vals <- tryCatch(brms::neff_ratio(model), error = function(e) NULL)
    if (!is.null(ess_vals)) {
      diag_df$ESS_Bulk <- as.numeric(ess_vals)[
        match(diag_df$Parameter, names(ess_vals))
      ]
    }
    diag_df$Converged <- !is.na(diag_df$Rhat) & diag_df$Rhat <= 1.01
  } else {
    diag_df <- data.frame(
      Parameter = character(0), Rhat = numeric(0),
      Converged = logical(0), stringsAsFactors = FALSE
    )
  }

  # ── Back-transform helper ──────────────────────────────────────────────────
  bt <- function(x) {
    if (transform == "log")  return(exp(x))
    if (transform == "sqrt") return(x ^ 2)
    x
  }

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

  vol_mat <- bt(epred)

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
  bliss_expected  <- fe_a + fe_b - fe_a * fe_b
  bliss_excess    <- fe_combo - bliss_expected

  bliss_q     <- stats::quantile(bliss_excess,   c(0.025, 0.5, 0.975))
  bliss_fe_q  <- stats::quantile(bliss_expected, c(0.025, 0.5, 0.975))
  obs_fe_med  <- stats::median(fe_combo)

  bliss_summary <- list(
    Excess_Median      = round(bliss_q["50%"],         4),
    Excess_Lower       = round(bliss_q["2.5%"],        4),
    Excess_Upper       = round(bliss_q["97.5%"],       4),
    P_Synergy          = round(mean(bliss_excess > 0), 3),
    Expected_FE_Median = round(bliss_fe_q["50%"],      3),
    Observed_FE_Median = round(obs_fe_med,             3)
  )

  # ── Draw-wise Loewe CI ────────────────────────────────────────────────────
  loewe_num  <- pmin(fe_a + fe_b, 1)
  loewe_ci   <- loewe_num / pmax(fe_combo, 1e-6)

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
    Interpretation = loewe_interp
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
    model_type_used  = "bayes_synergy",
    model            = if (isTRUE(return_model)) model else NULL,
    transform_used   = transform,
    tgi_summary      = tgi_summary,
    bliss_summary    = bliss_summary,
    loewe_summary    = loewe_summary,
    synergy_table    = synergy_table,
    posterior_summary = post_sum_df,
    mcmc_diagnostics = diag_df,
    summary          = analysis_summary,
    synergy_plot     = synergy_plot,
    posterior_dist_plot = post_dist_plot
  )
  out
}
