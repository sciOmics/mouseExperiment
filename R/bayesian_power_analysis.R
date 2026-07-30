# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Bayesian A Priori Power Analysis via Simulation
#'
#' Estimates the Bayesian power — the probability that a study of size N
#' yields a posterior probability of a meaningful treatment effect exceeding a
#' specified target — using simulation. For each candidate sample size a brms
#' longitudinal model is fitted to \code{n_simulations} synthetic datasets; the
#' fraction of fits where the posterior probability that at least one treatment
#' group's Treatment:Day interaction is more negative than \eqn{-\delta}
#' exceeds \code{target_prob} is the estimated Bayesian power at that N.
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
#' @param n_groups Total number of treatment groups including control. Default
#'   \code{2L} (single treatment vs control). For multi-group designs (e.g. a
#'   dose-escalation with three treated arms vs control, \code{n_groups = 4L})
#'   the simulator replicates \code{treatment_effect} across all non-control
#'   groups unless a vector of length \code{n_groups - 1} is supplied. The
#'   power success criterion under multi-group is "at least one Treatment:Day
#'   interaction has P(β < −δ | y) > target_prob" — the natural multi-arm
#'   analogue of the single-treatment case.
#' @param control_growth_rate Exponential growth rate of the control group
#'   on the log-volume scale (per day). Default \code{0.15}.
#' @param treatment_effect True reduction in log-scale growth rate for each
#'   treated group relative to control (positive = inhibition). Scalar (same
#'   effect for every non-control group) or numeric vector of length
#'   \code{n_groups - 1L}. Default \code{0.10}.
#' @param random_effects_specification One of \code{"intercept_only"} (default,
#'   \code{(1 | ID)}) or \code{"slope"} (\code{(Day | ID)}, per-animal random
#'   slopes). Match the spec used by your downstream
#'   \code{\link{bayesian_tumor_growth}} fit so the power estimate reflects
#'   the model you'll actually run.
#' @param animal_intercept_sd Between-animal random-intercept SD on the
#'   log-volume scale. Default \code{0.20}.
#' @param animal_slope_sd Between-animal random-slope SD on the log-volume
#'   scale (only used when \code{random_effects_specification = "slope"}).
#'   Default \code{0.05}.
#' @param residual_sd Within-animal residual SD on the log-volume scale.
#'   Default \code{0.15}.
#' @param baseline_log_volume Mean log-volume at day 0. Default \code{log(100)}.
#' @param timepoints Numeric vector of study days. Default \code{0:14}.
#' @param n_simulations Number of synthetic datasets to simulate per N.
#'   Default \code{20L}. Keep small; see Computational note above.
#' @param target_prob Posterior probability threshold: a simulation is counted
#'   as a "success" when \eqn{P(\beta < -\delta \mid y) >} \code{target_prob}
#'   for at least one Treatment:Day coefficient. Default \code{0.95}.
#' @param effect_threshold Minimum log-scale slope difference that constitutes
#'   a meaningful effect (\eqn{\delta} above). \code{0} (default) means any
#'   treatment-induced reduction in growth is scored as a success.
#' @param prior_strength Prior preset passed to the brms sub-models. One of
#'   \code{"skeptical"} (default), \code{"weakly_informative"},
#'   \code{"informative"}, \code{"diffuse"}.
#' @param null_calibration Logical. When \code{TRUE}, runs an additional
#'   simulation with \code{treatment_effect = 0} (everything else held
#'   constant) and reports the empirical false-positive rate at the same
#'   target_prob threshold. Under a calibrated decision rule the false
#'   positive rate should be at most \code{1 - target_prob}. Substantially
#'   doubles wall-clock time. Default \code{FALSE}.
#' @param n_chains Number of MCMC chains per fit. Default \code{2L} for speed.
#' @param n_warmup Warm-up (burn-in) iterations per chain. Default \code{500L}.
#' @param n_iter Post-warmup draws per chain. Default \code{500L} for speed.
#' @param seed Integer random seed for reproducibility. Default \code{42L}.
#' @param progress_fn Optional callback \code{function(n, sim_i, n_sims)}
#'   called after each simulation fit. Useful for Shiny progress bars.
#' @param backend brms backend: \code{"rstan"} (default) or \code{"cmdstanr"}.
#'   With cmdstanr the compile-once-and-reuse pattern via
#'   \code{brms::update()} runs noticeably faster across many simulations.
#'   See \code{\link{bayesian_tumor_growth}} for installation details.
#'
#' @return A named list:
#' \describe{
#'   \item{\code{power_table}}{Data frame: \code{N_Per_Group},
#'     \code{Bayes_Power} (proportion of successes), \code{SE_Power},
#'     \code{N_Success}, \code{N_Simulations}. When
#'     \code{null_calibration = TRUE}, also includes \code{Type_I_Rate}
#'     and \code{SE_Type_I} for the same Ns.}
#'   \item{\code{power_curve_data}}{Same as \code{power_table} — for direct
#'     \pkg{ggplot2} use.}
#'   \item{\code{params}}{Named list of simulation parameters used.}
#'   \item{\code{power_curve_plot}}{\pkg{ggplot2} power curve with 95 %
#'     Wilson CI error bars. When \code{null_calibration = TRUE}, also
#'     overlays the empirical Type-I curve. \code{NULL} if \pkg{ggplot2}
#'     is unavailable.}
#' }
#'
#' @importFrom stats rnorm qnorm
#' @export
bayesian_power_analysis <- function(
  n_per_group         = c(5L, 8L, 10L, 12L, 15L),
  n_groups            = 2L,
  control_growth_rate = 0.15,
  treatment_effect    = 0.10,
  random_effects_specification = c("intercept_only", "slope"),
  animal_intercept_sd = 0.20,
  animal_slope_sd     = 0.05,
  residual_sd         = 0.15,
  baseline_log_volume = log(100),
  timepoints          = 0:14,
  n_simulations       = 20L,
  target_prob         = 0.95,
  effect_threshold    = 0,
  prior_strength      = c(
    "skeptical", "weakly_informative", "informative", "diffuse"
  ),
  null_calibration    = FALSE,
  n_chains            = 2L,
  n_warmup            = 500L,
  n_iter              = 500L,
  seed                = 42L,
  progress_fn         = NULL,
  backend             = c("rstan", "cmdstanr")
) {


  prior_str     <- match.arg(prior_strength)
  re_spec       <- match.arg(random_effects_specification)
  backend       <- resolve_brms_backend(backend)
  n_per_group   <- as.integer(n_per_group)
  n_groups      <- as.integer(n_groups)
  n_simulations <- as.integer(n_simulations)
  n_chains      <- as.integer(n_chains)
  n_warmup      <- as.integer(n_warmup)
  n_iter        <- as.integer(n_iter)
  timepoints    <- as.numeric(timepoints)

  if (any(n_per_group < 2L)) stop("All 'n_per_group' values must be >= 2.")
  if (n_simulations < 5L)    stop("'n_simulations' must be >= 5.")
  if (length(timepoints) < 2L) {
    stop("'timepoints' must have at least 2 values.")
  }
  if (n_groups < 2L) {
    stop("'n_groups' must be >= 2 (one control + at least one treatment).")
  }

  # Treatment effect: scalar → replicate; vector → must match n_groups - 1.
  n_trt <- n_groups - 1L
  if (length(treatment_effect) == 1L) {
    treatment_effects <- rep(treatment_effect, n_trt)
  } else if (length(treatment_effect) == n_trt) {
    treatment_effects <- as.numeric(treatment_effect)
  } else {
    stop(sprintf(
      "'treatment_effect' must be a scalar or length-%d vector (n_groups - 1).",
      n_trt
    ))
  }
  if (any(treatment_effects < 0)) {
    stop("'treatment_effect' values must be >= 0 (positive = inhibition).")
  }

  group_names <- c("Control",
                   if (n_trt == 1L) "Treatment"
                   else sprintf("Treatment_%d", seq_len(n_trt)))

  re_formula_term <- if (re_spec == "slope") "(Day | ID)" else "(1 | ID)"
  bf_formula <- stats::as.formula(
    paste("Log_Vol ~ Treatment * Day +", re_formula_term)
  )
  # Priors are built below, from a pilot simulated dataset — see the note after
  # sim_data(). They cannot be constructed here because the data-scaled priors
  # (CODE_REVIEW.md R3.8) need a response to scale against, and no data exist
  # yet at this point.
  priors <- NULL

  # ── Simulation helper ──────────────────────────────────────────────────────
  sim_data <- function(n, rng_seed, trt_effects) {
    set.seed(rng_seed)
    n_t       <- length(timepoints)
    n_animals <- n * n_groups
    ids       <- rep(seq_len(n_animals), each = n_t)
    grp       <- rep(
      unlist(lapply(group_names, function(g) rep(g, n))),
      each = n_t
    )
    days      <- rep(timepoints, times = n_animals)

    # Per-group true slope: control = control_growth_rate, treated groups =
    # control_growth_rate - effect[g-1]
    group_slope <- c(control_growth_rate,
                     control_growth_rate - trt_effects)
    names(group_slope) <- group_names
    rates       <- group_slope[grp]

    u_i       <- rep(
      stats::rnorm(n_animals, 0, animal_intercept_sd),
      each = n_t
    )
    v_i_per_animal <- if (re_spec == "slope") {
      stats::rnorm(n_animals, 0, animal_slope_sd)
    } else rep(0, n_animals)
    v_i       <- rep(v_i_per_animal, each = n_t)
    eps       <- stats::rnorm(n_animals * n_t, 0, residual_sd)
    log_vol   <- baseline_log_volume + u_i + (rates + v_i) * days + eps

    data.frame(
      ID        = factor(ids),
      Treatment = factor(grp, levels = group_names),
      Day       = days,
      Log_Vol   = log_vol,
      stringsAsFactors = FALSE
    )
  }

  # ── Priors, from a pilot simulated dataset ─────────────────────────────────
  #
  # CODE_REVIEW.md R3.8 / J.6 — the analysis pipeline now uses data-scaled,
  # per-coefficient priors, so the power simulation must use the same ones or its
  # estimates do not describe the pipeline the user will actually run. Build them
  # from one pilot dataset at the largest N: the simulation parameters are fixed
  # across iterations, so the scale constants are too, which keeps the prior
  # structure and values constant and lets brms compile the Stan model once.
  priors <- bayes_scaled_priors(
    bf_formula,
    sim_data(max(n_per_group), seed, treatment_effects),
    response = "Log_Vol", prior_strength = prior_str,
    time_column = "Day", include_sd = TRUE
  )

  # ── Per-N power estimation ─────────────────────────────────────────────────
  # Returns a list with n_success and (used count of) fitted sims.
  estimate_per_N <- function(n, trt_effects, label) {
    n_success  <- 0L
    base_model <- NULL
    n_used     <- 0L

    for (si in seq_len(n_simulations)) {
      sim_seed <- seed * 1000L + (which(n_per_group == n) - 1L) * 100L + si
      # Different seed offset for null calibration so it's an independent draw
      if (identical(label, "null")) sim_seed <- sim_seed + 7919L
      df_sim   <- sim_data(n, sim_seed, trt_effects)

      fit <- tryCatch({
        if (is.null(base_model)) {
          brms::brm(
            formula = bf_formula,
            data    = df_sim,
            prior   = priors,
            chains  = n_chains,
            cores   = n_chains,
            iter    = n_warmup + n_iter,
            warmup  = n_warmup,
            seed    = seed,
            backend = backend,
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
      if (n_used == 0L) base_model <- fit
      n_used <- n_used + 1L

      fe_draws <- tryCatch(
        brms::fixef(fit, summary = FALSE),
        error = function(e) NULL
      )
      if (is.null(fe_draws)) next

      # Find every Treatment:Day interaction coefficient (one per treated group)
      int_cols <- grep("Treatment.*:Day|Day:Treatment",
                       colnames(fe_draws), value = TRUE)
      if (length(int_cols) == 0L) next

      # Per-coefficient posterior probability of "effect more negative than -δ"
      p_effects <- vapply(int_cols, function(cc) {
        mean(fe_draws[, cc] < -effect_threshold)
      }, numeric(1L))

      # Multi-group success criterion: at least one treated group passes.
      if (max(p_effects) > target_prob) n_success <- n_success + 1L

      if (!is.null(progress_fn)) progress_fn(n, si, n_simulations)
    }

    list(n_success = n_success, n_used = n_used)
  }

  zero_effects <- rep(0, n_trt)

  power_rows <- vector("list", length(n_per_group))
  for (ni in seq_along(n_per_group)) {
    n <- n_per_group[[ni]]
    pwr <- estimate_per_N(n, treatment_effects, label = "alt")
    bayes_power <- pwr$n_success / n_simulations
    se_power    <- sqrt(bayes_power * (1 - bayes_power) / n_simulations)

    row <- data.frame(
      N_Per_Group   = n,
      Bayes_Power   = round(bayes_power, 3),
      SE_Power      = round(se_power,    3),
      N_Success     = pwr$n_success,
      N_Simulations = n_simulations,
      stringsAsFactors = FALSE
    )

    if (isTRUE(null_calibration)) {
      null_pwr  <- estimate_per_N(n, zero_effects, label = "null")
      type_i    <- null_pwr$n_success / n_simulations
      se_type_i <- sqrt(type_i * (1 - type_i) / n_simulations)
      row$Type_I_Rate <- round(type_i,    3)
      row$SE_Type_I   <- round(se_type_i, 3)
    }

    power_rows[[ni]] <- row
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
        title    = paste0(
          "Bayesian Power Analysis",
          if (n_groups > 2L) sprintf(" (%d groups)", n_groups) else ""
        ),
        subtitle = paste0(
          "max P(β < -", effect_threshold, ") > ", target_prob,
          " across treated groups; prior = ", prior_str,
          if (re_spec == "slope") "; RE: (Day | ID)" else "; RE: (1 | ID)"
        ),
        x = "Animals per Group",
        y = paste0("Bayesian Power (", n_simulations, " sims / N)")
      ) +
      ggplot2::theme_classic(base_size = 14)

    if (isTRUE(null_calibration)) {
      # Overlay empirical Type-I curve in tomato
      power_curve_plot <- power_curve_plot +
        ggplot2::geom_line(
          ggplot2::aes(y = .data[["Type_I_Rate"]]),
          colour = "tomato", linewidth = 0.8, linetype = "dotdash"
        ) +
        ggplot2::geom_point(
          ggplot2::aes(y = .data[["Type_I_Rate"]]),
          colour = "tomato", size = 2.2
        ) +
        ggplot2::labs(
          caption = sprintf(
            "Dot-dash line = empirical Type-I rate under H0 (target = %.2f).",
            1 - target_prob
          )
        )
    }
  }

  list(
    power_table      = power_table,
    power_curve_data = power_table,
    params           = list(
      n_per_group         = n_per_group,
      n_groups            = n_groups,
      control_growth_rate = control_growth_rate,
      treatment_effect    = treatment_effects,
      random_effects_specification = re_spec,
      animal_intercept_sd = animal_intercept_sd,
      animal_slope_sd     = animal_slope_sd,
      residual_sd         = residual_sd,
      baseline_log_volume = baseline_log_volume,
      timepoints          = timepoints,
      n_simulations       = n_simulations,
      target_prob         = target_prob,
      effect_threshold    = effect_threshold,
      prior_strength      = prior_str,
      null_calibration    = null_calibration,
      n_chains            = n_chains,
      n_warmup            = n_warmup,
      n_iter              = n_iter,
      seed                = seed
    ),
    power_curve_plot = power_curve_plot
  )
}
