# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Bayesian A Priori Power Analysis via Simulation
#'
#' Estimates the Bayesian power — the probability that a study of size N
#' yields a posterior probability of a meaningful treatment effect exceeding a
#' specified target — using simulation. For each candidate sample size a brms
#' longitudinal model is fitted to \code{n_simulations} synthetic datasets; the
#' fraction of fits where \eqn{P(\beta_{\text{trt:Day}} < -\delta \mid y) >}
#' \code{target_prob} is the estimated Bayesian power at that N.
#'
#' \strong{Computational note:} Fitting brms models is slow. The compiled Stan
#' program is reused across simulations via \code{brms::update()}, but even so
#' \code{n_simulations} should be kept small (\eqn{\leq 50}) unless run on a
#' compute cluster. For quick exploratory use the defaults
#' (\code{n_chains = 2, n_iter = 1000, n_simulations = 20}) produce rough
#' estimates in a few minutes per sample-size candidate.
#'
#' @param n_per_group Integer vector of animals-per-group values to evaluate.
#'   Default \code{c(5L, 8L, 10L, 12L, 15L)}.
#' @param control_growth_rate Exponential growth rate of the control group
#'   on the log-volume scale (per day). Default \code{0.15}.
#' @param treatment_effect True reduction in log-scale growth rate for the
#'   treatment group relative to control (positive = inhibition). Default
#'   \code{0.10}.
#' @param animal_intercept_sd Between-animal random-intercept SD on the
#'   log-volume scale. Default \code{0.20}.
#' @param residual_sd Within-animal residual SD on the log-volume scale.
#'   Default \code{0.15}.
#' @param baseline_log_volume Mean log-volume at day 0. Default \code{log(100)}.
#' @param timepoints Numeric vector of study days. Default \code{0:14}.
#' @param n_simulations Number of synthetic datasets to simulate per N.
#'   Default \code{20L}. Keep small; see Computational note above.
#' @param target_prob Posterior probability threshold: a simulation is counted
#'   as a "success" when \eqn{P(\beta < -\delta \mid y) >} \code{target_prob}.
#'   Default \code{0.95}.
#' @param effect_threshold Minimum log-scale slope difference that constitutes
#'   a meaningful effect (\eqn{\delta} above). \code{0} (default) means any
#'   treatment-induced reduction in growth is scored as a success.
#' @param prior_strength Prior preset passed to the brms sub-models. One of
#'   \code{"skeptical"} (default), \code{"weakly_informative"},
#'   \code{"informative"}, \code{"diffuse"}.
#' @param n_chains Number of MCMC chains per fit. Default \code{2L} for speed.
#' @param n_iter Iterations per chain (including warm-up). Default \code{1000L}
#'   for speed.
#' @param seed Integer random seed for reproducibility. Default \code{42L}.
#' @param progress_fn Optional callback \code{function(n, sim_i, n_sims)}
#'   called after each simulation fit. Useful for Shiny progress bars.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{power_table}}{Data frame: \code{N_Per_Group},
#'     \code{Bayes_Power} (proportion of successes), \code{SE_Power},
#'     \code{N_Success}, \code{N_Simulations}.}
#'   \item{\code{power_curve_data}}{Same as \code{power_table} — for direct
#'     \pkg{ggplot2} use.}
#'   \item{\code{params}}{Named list of simulation parameters used.}
#'   \item{\code{power_curve_plot}}{\pkg{ggplot2} power curve with 95 %
#'     Wilson CI error bars. \code{NULL} if \pkg{ggplot2} is unavailable.}
#' }
#'
#' @importFrom stats rnorm qnorm
#' @export
bayesian_power_analysis <- function(
  n_per_group         = c(5L, 8L, 10L, 12L, 15L),
  control_growth_rate = 0.15,
  treatment_effect    = 0.10,
  animal_intercept_sd = 0.20,
  residual_sd         = 0.15,
  baseline_log_volume = log(100),
  timepoints          = 0:14,
  n_simulations       = 20L,
  target_prob         = 0.95,
  effect_threshold    = 0,
  prior_strength      = c(
    "skeptical", "weakly_informative", "informative", "diffuse"
  ),
  n_chains            = 2L,
  n_iter              = 1000L,
  seed                = 42L,
  progress_fn         = NULL
) {

  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "'brms' is required for bayesian_power_analysis(). ",
      "Install it with: install.packages('brms')"
    )
  }

  prior_str     <- match.arg(prior_strength)
  n_per_group   <- as.integer(n_per_group)
  n_simulations <- as.integer(n_simulations)
  n_chains      <- as.integer(n_chains)
  n_iter        <- as.integer(n_iter)
  timepoints    <- as.numeric(timepoints)

  if (any(n_per_group < 2L)) stop("All 'n_per_group' values must be >= 2.")
  if (n_simulations < 5L)    stop("'n_simulations' must be >= 5.")
  if (length(timepoints) < 2L) {
    stop("'timepoints' must have at least 2 values.")
  }
  if (treatment_effect < 0) {
    stop("'treatment_effect' must be >= 0 (positive = inhibition).")
  }

  pp       <- bayes_prior_params(prior_str)
  b_sd     <- pp$b_sd
  exp_rate <- pp$exp_rate
  priors   <- c(
    brms::prior_string(paste0("normal(0, ", b_sd, ")"),
                       class = "b"),
    brms::prior_string(paste0("normal(0, ", b_sd * 2.5, ")"),
                       class = "Intercept"),
    brms::prior_string(paste0("exponential(", exp_rate, ")"),
                       class = "sd"),
    brms::prior_string(paste0("exponential(", exp_rate, ")"),
                       class = "sigma")
  )

  # ── Simulation helper ──────────────────────────────────────────────────────
  sim_data <- function(n, rng_seed) {
    set.seed(rng_seed)
    n_t       <- length(timepoints)
    n_animals <- n * 2L
    ids       <- rep(seq_len(n_animals), each = n_t)
    grp       <- rep(
      c(rep("Control", n), rep("Treatment", n)),
      each = n_t
    )
    days      <- rep(timepoints, times = n_animals)
    rates     <- ifelse(
      grp == "Control",
      control_growth_rate,
      control_growth_rate - treatment_effect
    )
    u_i       <- rep(
      stats::rnorm(n_animals, 0, animal_intercept_sd),
      each = n_t
    )
    eps       <- stats::rnorm(n_animals * n_t, 0, residual_sd)
    log_vol   <- baseline_log_volume + u_i + rates * days + eps

    data.frame(
      ID        = factor(ids),
      Treatment = factor(grp, levels = c("Control", "Treatment")),
      Day       = days,
      Log_Vol   = log_vol,
      stringsAsFactors = FALSE
    )
  }

  # ── Per-N power estimation ─────────────────────────────────────────────────
  power_rows <- vector("list", length(n_per_group))

  for (ni in seq_along(n_per_group)) {
    n          <- n_per_group[[ni]]
    n_success  <- 0L
    base_model <- NULL

    for (si in seq_len(n_simulations)) {
      sim_seed <- seed * 1000L + (ni - 1L) * 100L + si
      df_sim   <- sim_data(n, sim_seed)

      fit <- tryCatch({
        if (is.null(base_model)) {
          brms::brm(
            formula = Log_Vol ~ Treatment * Day + (1 | ID),
            data    = df_sim,
            prior   = priors,
            chains  = n_chains,
            iter    = n_iter,
            seed    = seed,
            silent  = 2L,
            refresh = 0L
          )
        } else {
          brms::update(
            base_model,
            newdata  = df_sim,
            recompile = FALSE,
            silent   = 2L,
            refresh  = 0L
          )
        }
      }, error = function(e) NULL)

      if (is.null(fit)) next
      if (si == 1L) base_model <- fit

      fe_draws <- tryCatch(
        brms::fixef(fit, summary = FALSE),
        error = function(e) NULL
      )
      if (is.null(fe_draws)) next

      int_col <- grep("Treatment.*Day|Day.*Treatment",
                      colnames(fe_draws), value = TRUE)
      if (length(int_col) == 0L) next
      int_col <- int_col[[1L]]

      p_effect <- mean(fe_draws[, int_col] < -effect_threshold)
      if (p_effect > target_prob) n_success <- n_success + 1L

      if (!is.null(progress_fn)) progress_fn(n, si, n_simulations)
    }

    bayes_power <- n_success / n_simulations
    se_power    <- sqrt(bayes_power * (1 - bayes_power) / n_simulations)

    power_rows[[ni]] <- data.frame(
      N_Per_Group   = n,
      Bayes_Power   = round(bayes_power, 3),
      SE_Power      = round(se_power,    3),
      N_Success     = n_success,
      N_Simulations = n_simulations,
      stringsAsFactors = FALSE
    )
  }

  power_table <- do.call(rbind, power_rows)
  rownames(power_table) <- NULL

  # ── Power curve plot ───────────────────────────────────────────────────────
  power_curve_plot <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    z   <- stats::qnorm(0.975)
    pt  <- power_table
    pt$Lower <- pmax(0, pt$Bayes_Power - z * pt$SE_Power)
    pt$Upper <- pmin(1, pt$Bayes_Power + z * pt$SE_Power)

    power_curve_plot <- ggplot2::ggplot(
      pt,
      ggplot2::aes(
        x = .data[["N_Per_Group"]],
        y = .data[["Bayes_Power"]]
      )
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[["Lower"]],
          ymax = .data[["Upper"]]
        ),
        fill = "steelblue", alpha = 0.2
      ) +
      ggplot2::geom_line(colour = "steelblue", linewidth = 0.9) +
      ggplot2::geom_point(size = 2.5, colour = "steelblue") +
      ggplot2::geom_hline(
        yintercept = target_prob,
        linetype = "dashed", colour = "tomato", linewidth = 0.7
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        title    = "Bayesian Power Analysis",
        subtitle = paste0(
          "P(β < -", effect_threshold, ") > ", target_prob,
          " as success criterion; prior = ", prior_str
        ),
        x = "Animals per Group",
        y = paste0("Bayesian Power (", n_simulations, " sims / N)")
      ) +
      ggplot2::theme_classic(base_size = 14)
  }

  list(
    power_table      = power_table,
    power_curve_data = power_table,
    params           = list(
      n_per_group         = n_per_group,
      control_growth_rate = control_growth_rate,
      treatment_effect    = treatment_effect,
      animal_intercept_sd = animal_intercept_sd,
      residual_sd         = residual_sd,
      baseline_log_volume = baseline_log_volume,
      timepoints          = timepoints,
      n_simulations       = n_simulations,
      target_prob         = target_prob,
      effect_threshold    = effect_threshold,
      prior_strength      = prior_str,
      n_chains            = n_chains,
      n_iter              = n_iter,
      seed                = seed
    ),
    power_curve_plot = power_curve_plot
  )
}
