#' A Priori Power Analysis via LMM Simulation (Option 3)
#'
#' Estimates prospective power by simulating longitudinal tumour-growth
#' datasets from user-supplied growth parameters, fitting the same linear
#' mixed-effects model used by \code{tumor_growth_statistics()}, and
#' recording the proportion of simulations where the Treatment × Day
#' interaction is significant.
#'
#' @section Model:
#' Simulated data: \eqn{y_{ij} = (\mu_i + b_{0,j}) + (\beta_i + b_{1,j}) \cdot t_{ij} + \varepsilon_{ij}}
#' where \eqn{\mu_i} = group baseline, \eqn{\beta_i} = group growth rate,
#' \eqn{b_0 \sim N(0, \sigma_{b0}^2)}, \eqn{b_1 \sim N(0, \sigma_{b1}^2)},
#' \eqn{\varepsilon \sim N(0, \sigma_e^2)}.
#' The fitted model is \code{Volume ~ Treatment * Day + (Day | ID)};
#' the p-value for the \code{Treatment:Day} interaction is extracted.
#'
#' @param n_per_group Integer vector. Sample sizes per group to evaluate.
#'   Power is computed for each value independently. Default \code{c(5, 8, 10, 12, 15)}.
#' @param n_groups Integer ≥ 2. Number of treatment groups (default 2). When
#'   > 2, all non-control groups share the same \code{treatment_effect}.
#' @param baseline_volume Numeric > 0. Mean tumour volume at Day 0 (default 100).
#' @param baseline_sd Numeric > 0. SD of baseline volume (default 20).
#' @param control_growth_rate Numeric. Daily growth rate for the control group
#'   in the same units as \code{baseline_volume} (default 15).
#' @param treatment_effect Numeric. Reduction in daily growth rate for treated
#'   groups relative to control (default 10; i.e. treated rate = control_rate −
#'   treatment_effect). Positive values mean tumour suppression.
#' @param random_slope_sd Numeric ≥ 0. SD of the per-mouse random slope
#'   (default 3).
#' @param random_intercept_sd Numeric ≥ 0. SD of the per-mouse random
#'   intercept (default 10).
#' @param residual_sd Numeric > 0. Residual (measurement) SD (default 10).
#' @param timepoints Integer vector. Study days at which measurements are
#'   taken (default \code{0:14}).
#' @param alpha Numeric vector of significance levels (default \code{c(0.05)}).
#' @param n_simulations Integer. Number of LMM fits per (N, alpha) combination
#'   (default 500). Values ≥ 1000 give more stable estimates but are slow.
#' @param seed Integer or \code{NULL}. RNG seed for reproducibility (default
#'   \code{NULL}).
#' @param progress_fn Function or \code{NULL}. Called with a single numeric
#'   argument (proportion complete 0–1) after each completed N scenario; use
#'   for progress reporting in Shiny.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{power_table}}{Data frame: N_Per_Group, N_Groups, Alpha,
#'     Power, SE_Power, N_Significant, N_Simulations.}
#'   \item{\code{power_curve_data}}{Same as \code{power_table} — suitable
#'     for direct ggplot2 use.}
#'   \item{\code{params}}{Named list of the simulation parameters used.}
#' }
#'
#' @export
apriori_power_simulation <- function(n_per_group          = c(5L, 8L, 10L, 12L, 15L),
                                     n_groups             = 2L,
                                     baseline_volume      = 100,
                                     baseline_sd          = 20,
                                     control_growth_rate  = 15,
                                     treatment_effect     = 10,
                                     random_slope_sd      = 3,
                                     random_intercept_sd  = 10,
                                     residual_sd          = 10,
                                     timepoints           = 0:14,
                                     alpha                = c(0.05),
                                     n_simulations        = 500L,
                                     seed                 = NULL,
                                     progress_fn          = NULL) {

  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required for simulation-based a priori power analysis.")

  n_groups     <- as.integer(n_groups)
  n_simulations <- as.integer(n_simulations)
  n_per_group  <- as.integer(n_per_group)
  if (n_groups < 2L)   stop("n_groups must be >= 2.")
  if (n_simulations < 10L) stop("n_simulations must be >= 10.")
  if (length(timepoints) < 2) stop("timepoints must have at least 2 values.")
  if (!is.null(seed)) set.seed(seed)

  group_names   <- c("Control", paste0("Treatment", seq_len(n_groups - 1L)))
  growth_rates  <- c(control_growth_rate,
                     rep(control_growth_rate - treatment_effect, n_groups - 1L))

  n_timepoints <- length(timepoints)
  total_steps  <- length(n_per_group) * length(alpha)
  step_i       <- 0L

  rows <- vector("list", total_steps)
  idx  <- 0L

  for (n in n_per_group) {
    for (a in alpha) {
      n_sig <- 0L
      n_fit <- 0L

      for (sim in seq_len(n_simulations)) {
        # ---- Simulate data ---------------------------------------------------
        sim_list <- vector("list", n_groups)
        for (g in seq_len(n_groups)) {
          ids <- paste0(group_names[g], "_", seq_len(n))
          b0  <- stats::rnorm(n, 0, random_intercept_sd)
          b1  <- stats::rnorm(n, 0, random_slope_sd)
          for (j in seq_len(n)) {
            vol <- (baseline_volume + b0[j]) +
              (growth_rates[g] + b1[j]) * timepoints +
              stats::rnorm(n_timepoints, 0, residual_sd)
            # Clamp to avoid negative volumes (log-scale artifacts)
            vol <- pmax(vol, 0.1)
            sim_list[[(g - 1L) * n + j]] <- data.frame(
              ID        = ids[j],
              Treatment = group_names[g],
              Day       = timepoints,
              Volume    = vol,
              stringsAsFactors = FALSE
            )
          }
        }
        df_sim <- do.call(rbind, sim_list)
        df_sim$Treatment <- factor(df_sim$Treatment, levels = group_names)

        # ---- Fit LMM ---------------------------------------------------------
        fit <- tryCatch(
          lme4::lmer(Volume ~ Treatment * Day + (Day | ID),
                     data    = df_sim,
                     REML    = FALSE,
                     control = lme4::lmerControl(optimizer = "bobyqa",
                                                 calc.derivs = FALSE)),
          error   = function(e) NULL,
          warning = function(w) {
            tryCatch(
              lme4::lmer(Volume ~ Treatment * Day + (1 | ID),
                         data    = df_sim,
                         REML    = FALSE,
                         control = lme4::lmerControl(optimizer = "bobyqa",
                                                     calc.derivs = FALSE)),
              error = function(e2) NULL
            )
          }
        )
        if (is.null(fit)) next
        n_fit <- n_fit + 1L

        # ---- Extract p-value for Treatment:Day ------------------------------
        p_val <- tryCatch({
          coefs  <- as.data.frame(stats::coef(summary(fit)))
          # Rows whose name contains ":"  — the interaction terms
          int_rows <- coefs[grep(":", rownames(coefs), fixed = TRUE), , drop = FALSE]
          if (nrow(int_rows) == 0) {
            # Fall back: Satterthwaite via lmerTest if available
            if (requireNamespace("lmerTest", quietly = TRUE)) {
              fit2 <- lmerTest::as_lmerModLmerTest(fit)
              coefs2 <- as.data.frame(stats::coef(summary(fit2)))
              int_rows2 <- coefs2[grep(":", rownames(coefs2), fixed = TRUE), , drop = FALSE]
              if (nrow(int_rows2) > 0 && "Pr(>|t|)" %in% colnames(int_rows2))
                min(int_rows2[["Pr(>|t|)"]], na.rm = TRUE)
              else NA_real_
            } else NA_real_
          } else {
            # Use likelihood-ratio test between full and reduced model
            fit_red <- tryCatch(
              lme4::lmer(Volume ~ Treatment + Day + (Day | ID),
                         data    = df_sim,
                         REML    = FALSE,
                         control = lme4::lmerControl(optimizer = "bobyqa",
                                                     calc.derivs = FALSE)),
              error = function(e) NULL
            )
            if (!is.null(fit_red)) {
              lrt <- stats::anova(fit_red, fit)
              lrt$`Pr(>Chisq)`[2]
            } else NA_real_
          }
        }, error = function(e) NA_real_)

        if (!is.na(p_val) && p_val < a) n_sig <- n_sig + 1L
      }

      power_est <- if (n_fit > 0) n_sig / n_fit else NA_real_
      se_est    <- if (n_fit > 1 && !is.na(power_est))
        sqrt(power_est * (1 - power_est) / n_fit) else NA_real_

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        N_Per_Group   = n,
        N_Groups      = n_groups,
        Alpha         = a,
        Power         = power_est,
        SE_Power      = se_est,
        N_Significant = n_sig,
        N_Simulations = n_fit,
        stringsAsFactors = FALSE
      )

      step_i <- step_i + 1L
      if (is.function(progress_fn))
        progress_fn(step_i / total_steps)
    }
  }

  power_table <- do.call(rbind, rows[seq_len(idx)])

  list(
    power_table      = power_table,
    power_curve_data = power_table,   # same structure; alias for plot code
    params = list(
      n_per_group          = n_per_group,
      n_groups             = n_groups,
      baseline_volume      = baseline_volume,
      baseline_sd          = baseline_sd,
      control_growth_rate  = control_growth_rate,
      treatment_effect     = treatment_effect,
      random_slope_sd      = random_slope_sd,
      random_intercept_sd  = random_intercept_sd,
      residual_sd          = residual_sd,
      timepoints           = timepoints,
      alpha                = alpha,
      n_simulations        = n_simulations
    )
  )
}
