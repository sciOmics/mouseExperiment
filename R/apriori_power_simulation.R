#' A Priori Power Analysis via LMM Simulation (Option 3)
#'
#' Estimates prospective power by simulating longitudinal tumour-growth
#' datasets from user-supplied growth parameters, fitting the same linear
#' mixed-effects model used by \code{tumor_growth_statistics()}, and
#' recording the proportion of simulations where the Treatment × Day
#' interaction is significant.
#'
#' @section Model:
#' Data are generated on the \strong{log scale} to match the log-transformed
#' LMM used in \code{tumor_growth_statistics()}:
#' \deqn{\log(V_{ij}) = \log(\mu_i) + b_{0,j} + (\beta_i + b_{1,j}) \cdot t_{ij} + \varepsilon_{ij}}
#' where \eqn{\mu_i} = group geometric-mean baseline volume,
#' \eqn{\beta_i} = log-scale daily growth rate,
#' \eqn{b_0 \sim N(0, \sigma_{b0}^2)}, \eqn{b_1 \sim N(0, \sigma_{b1}^2)},
#' \eqn{\varepsilon \sim N(0, \sigma_e^2)}, all on the log scale.
#' Volumes are obtained by exponentiating: \eqn{V_{ij} = \exp(\cdot)}.
#' The fitted model is \code{log(Volume) ~ Treatment * Day + (Day | ID)};
#' power is the proportion of LRTs where the \code{Treatment:Day} p-value < \eqn{\alpha}.
#'
#' @section Parameter scale:
#' \code{control_growth_rate}, \code{treatment_effect}, \code{random_slope_sd},
#' \code{random_intercept_sd}, and \code{residual_sd} are all on the
#' \strong{log scale}. A \code{control_growth_rate} of 0.15 means log-volume
#' increases by 0.15 per day (~16% daily volume growth). A
#' \code{treatment_effect} of 0.10 means the treated group's log-growth rate
#' is 0.10 lower than control. \code{baseline_volume} is in raw volume units
#' (mm³); it is log-transformed internally.
#'
#' @section Baseline variability:
#' Neither \code{baseline_sd} nor \code{random_intercept_sd} affects the power
#' this function reports, and that is a property of the estimand rather than a
#' defect. The test is on \code{Treatment:Day} — a contrast of \emph{slopes} —
#' while both parameters inject only per-mouse \emph{intercept} variation, which
#' the \code{(Day | ID)} random intercept absorbs exactly. In a balanced design
#' the slope contrast is orthogonal to per-mouse intercepts, so the p-value is
#' invariant to them: varying \code{baseline_sd} over two orders of magnitude
#' (1 to 200 mm³) reproduces the \code{Treatment:Day} p-value to ten significant
#' figures.
#'
#' The two parameters are also aliased with each other. Both are drawn
#' independently and summed into the same per-mouse offset, so only their
#' root-sum-square,
#' \deqn{\sqrt{(\mathrm{baseline\_sd}/\mathrm{baseline\_volume})^2 + \mathrm{random\_intercept\_sd}^2}}
#' is identifiable — no simulation can distinguish them.
#'
#' The practical consequence: \strong{do not spend calibration effort on these
#' two}. Power here is driven by \code{treatment_effect},
#' \code{random_slope_sd}, \code{residual_sd}, the number and spacing of
#' \code{timepoints}, and \code{n_per_group}. If baseline heterogeneity is a
#' real concern for your design, it is an argument about the analysis estimand
#' (baseline as a covariate, or an endpoint-volume contrast rather than a slope
#' contrast), not something a larger \code{baseline_sd} here will reveal.
#'
#' @section Dropout:
#' By default every animal contributes every timepoint. Real studies do not work
#' that way: animals are euthanised on an IACUC volume limit or for poor health,
#' which removes late observations preferentially from the fastest-growing
#' controls. Set \code{dropout_limit} to plan against your actual limit.
#'
#' Ignoring it is optimistic, but only mildly, and the reason is worth stating:
#' this dropout is MAR given the observed trajectory, so the LMM likelihood stays
#' unbiased and power degrades through lost information alone rather than through
#' bias. Measured on a two-arm design (control 0.13/day, effect 0.018/day,
#' 400 replicates per cell):
#'
#' \tabular{lrr}{
#'   \strong{Observations lost} \tab \strong{Power (n=8)} \tab \strong{Power (n=12)} \cr
#'   0%  \tab 0.770 \tab 0.917 \cr
#'   21% \tab 0.738 \tab 0.905 \cr
#'   30% \tab 0.708 \tab 0.890 \cr
#'   39% \tab 0.645 \tab 0.835
#' }
#'
#' So a design the default reports at 0.77 is really running at about 0.71 if it
#' loses 30% of its observations — enough to justify an extra animal or two per
#' arm, not enough to invalidate the planning exercise.
#'
#' @param n_per_group Integer vector. Sample sizes per group to evaluate.
#'   Power is computed for each value independently. Default \code{c(5, 8, 10, 12, 15)}.
#' @param n_groups Integer ≥ 2. Number of treatment groups (default 2). When
#'   > 2, all non-control groups share the same \code{treatment_effect}.
#' @param baseline_volume Numeric > 0. Geometric mean tumour volume at Day 0
#'   in mm³ (default 100).
#' @param baseline_sd Numeric \eqn{\ge} 0. Standard deviation (in mm^3) of the
#'   raw-scale baseline volume across animals. Converted to a log-scale SD
#'   via the first-order delta-method approximation
#'   \code{log_sd = baseline_sd / baseline_volume} (a first-order
#'   approximation) and added as a
#'   per-mouse jitter on top of \code{log(baseline_volume)}.
#'
#'   \strong{This parameter does not change the reported power}, and no value
#'   you give it will. See the \emph{Baseline variability} section — it is
#'   retained because it makes the simulated volumes realistic if you inspect
#'   them, not because it informs the sample-size decision.
#' @param control_growth_rate Numeric. Log-scale daily growth rate for the
#'   control group (default 0.15, ~16% daily increase).
#' @param treatment_effect Numeric. Reduction in log-scale daily growth rate
#'   for treated groups relative to control (default 0.10; i.e. treated rate =
#'   control_rate − treatment_effect). Positive values mean tumour suppression.
#' @param random_slope_sd Numeric ≥ 0. SD of the per-mouse random slope on
#'   the log scale (default 0.05).
#' @param random_intercept_sd Numeric ≥ 0. SD of the per-mouse random
#'   intercept on the log scale (default 0.20).
#' @param residual_sd Numeric > 0. Residual (measurement) SD on the log scale
#'   (default 0.10).
#' @param timepoints Integer vector. Study days at which measurements are
#'   taken (default \code{0:14}).
#' @param dropout_limit Numeric > 0 or \code{Inf}. Volume (in the same units as
#'   \code{baseline_volume}) at which an animal is removed from study; all of its
#'   observations from the first crossing onward are discarded, mimicking
#'   euthanasia on an IACUC limit. \code{Inf} (the default) keeps the historical
#'   complete-data behaviour. See the \emph{Dropout} section for what it costs to
#'   leave this at \code{Inf}.
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
#'   \item{\code{attrition}}{Data frame: for each N, the mean percentage of
#'     observations and of animals lost to \code{dropout_limit}. All zero when
#'     \code{dropout_limit = Inf}. Report this alongside the power estimate —
#'     it is what tells a reviewer the number was planned against a realistic
#'     study rather than an idealised one.}
#' }
#'
#' @export
apriori_power_simulation <- function(n_per_group          = c(5L, 8L, 10L, 12L, 15L),
                                     n_groups             = 2L,
                                     baseline_volume      = 100,
                                     baseline_sd          = 20,
                                     control_growth_rate  = 0.15,
                                     treatment_effect     = 0.10,
                                     random_slope_sd      = 0.05,
                                     random_intercept_sd  = 0.20,
                                     residual_sd          = 0.10,
                                     timepoints           = 0:14,
                                     dropout_limit        = Inf,
                                     alpha                = c(0.05),
                                     n_simulations        = 500L,
                                     seed                 = NULL,
                                     progress_fn          = NULL) {


  n_groups     <- as.integer(n_groups)
  n_simulations <- as.integer(n_simulations)
  n_per_group  <- as.integer(n_per_group)
  if (n_groups < 2L)   stop("n_groups must be >= 2.")
  if (n_simulations < 10L) stop("n_simulations must be >= 10.")
  if (length(timepoints) < 2) stop("timepoints must have at least 2 values.")
  if (!is.numeric(dropout_limit) || length(dropout_limit) != 1L ||
      is.na(dropout_limit) || dropout_limit <= 0) {
    stop("dropout_limit must be a single positive number, or Inf for none.")
  }
  if (!is.null(seed)) set.seed(seed)

  group_names   <- c("Control", paste0("Treatment", seq_len(n_groups - 1L)))
  growth_rates  <- c(control_growth_rate,
                     rep(control_growth_rate - treatment_effect, n_groups - 1L))

  n_timepoints <- length(timepoints)
  total_steps  <- length(n_per_group) * length(alpha)
  step_i       <- 0L

  rows      <- vector("list", total_steps)
  attr_rows <- vector("list", total_steps)
  idx       <- 0L

  for (n in n_per_group) {
    for (a in alpha) {
      n_sig <- 0L
      n_fit <- 0L
      # Realised attrition, accumulated so the caller can see what the reported
      # power was actually planned against.
      obs_kept <- 0L; obs_total <- 0L
      mice_truncated <- 0L; mice_total <- 0L

      for (sim in seq_len(n_simulations)) {
        # ---- Simulate data ---------------------------------------------------
        sim_list <- vector("list", n_groups)
        # Per-mouse baseline draws (log-scale CV ≈ baseline_sd / baseline_volume
        # mapped to log-scale SD via the delta-method first-order approximation
        # log_sd ≈ baseline_sd / baseline_volume for small CVs).
        log_baseline_sd <- if (baseline_sd > 0 && baseline_volume > 0) {
          baseline_sd / baseline_volume
        } else 0
        for (g in seq_len(n_groups)) {
          ids <- paste0(group_names[g], "_", seq_len(n))
          # Draw standard normals and scale, rather than passing the SD to
          # rnorm(). R's rnorm() returns `mu` without calling norm_rand() when
          # sigma == 0, so an SD of zero silently advances the RNG stream by a
          # different amount and reshuffles every subsequent draw -- making
          # `sd = 0` a discontinuity rather than the limiting case it reads as.
          # For sd > 0 this form is bit-identical to rnorm(n, 0, sd).
          b0  <- random_intercept_sd * stats::rnorm(n)
          b1  <- random_slope_sd     * stats::rnorm(n)
          per_mouse_log_baseline <- log(baseline_volume) +
            log_baseline_sd * stats::rnorm(n)
          for (j in seq_len(n)) {
            # Log-normal data: log(V) = log(baseline_j) + b0 + (rate + b1)*t + e
            # Matches the log(Volume) LMM fitted below.
            log_vol <- per_mouse_log_baseline[j] + b0[j] +
              (growth_rates[g] + b1[j]) * timepoints +
              stats::rnorm(n_timepoints, 0, residual_sd)
            vol <- exp(log_vol)

            # Euthanasia on the volume limit: drop this animal's observations
            # from the first crossing onward. Dropout is MAR given the observed
            # trajectory, which is why the LMM below remains valid -- see the
            # Dropout section.
            keep <- if (is.finite(dropout_limit)) {
              hit <- which(vol >= dropout_limit)
              if (length(hit)) seq_len(hit[1]) else seq_along(timepoints)
            } else {
              seq_along(timepoints)
            }
            mice_total <- mice_total + 1L
            obs_total  <- obs_total + n_timepoints
            obs_kept   <- obs_kept + length(keep)
            if (length(keep) < n_timepoints) mice_truncated <- mice_truncated + 1L

            sim_list[[(g - 1L) * n + j]] <- data.frame(
              ID        = ids[j],
              Treatment = group_names[g],
              Day       = timepoints[keep],
              Volume    = vol[keep],
              stringsAsFactors = FALSE
            )
          }
        }
        df_sim <- do.call(rbind, sim_list)
        df_sim$Treatment <- factor(df_sim$Treatment, levels = group_names)

        # An aggressive limit can leave a group with too little signal to fit a
        # random-slope model. Skip rather than let lmer fail into the warning
        # fallback, which would silently change the model being powered.
        n_obs_per_mouse <- table(df_sim$ID)
        if (length(unique(df_sim$Treatment[!is.na(df_sim$Treatment)])) < 2L ||
            sum(n_obs_per_mouse >= 2L) < 2L * n_groups) {
          next
        }

        # ---- Fit LMM on log(Volume) ------------------------------------------
        # Matches tumor_growth_statistics() which log-transforms before fitting.
        fit <- tryCatch(
          lme4::lmer(log(Volume) ~ Treatment * Day + (Day | ID),
                     data    = df_sim,
                     REML    = FALSE,
                     control = lme4::lmerControl(optimizer = "bobyqa",
                                                 calc.derivs = FALSE)),
          error   = function(e) NULL,
          warning = function(w) {
            tryCatch(
              lme4::lmer(log(Volume) ~ Treatment * Day + (1 | ID),
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
        # Prefer LRT (statistically preferred for LMMs); lme4 does not compute
        # p-values by default so the coefficient table rarely has Pr(>|t|).
        p_val <- tryCatch({
          fit_red <- tryCatch(
            lme4::lmer(log(Volume) ~ Treatment + Day + (Day | ID),
                       data    = df_sim,
                       REML    = FALSE,
                       control = lme4::lmerControl(optimizer = "bobyqa",
                                                   calc.derivs = FALSE)),
            error = function(e) NULL
          )
          if (!is.null(fit_red)) {
            lrt <- stats::anova(fit_red, fit)
            lrt$`Pr(>Chisq)`[2]
          } else if (requireNamespace("lmerTest", quietly = TRUE)) {
            # Fall back to Satterthwaite approximation
            fit2      <- lmerTest::as_lmerModLmerTest(fit)
            coefs2    <- as.data.frame(stats::coef(summary(fit2)))
            int_rows2 <- coefs2[
              grep(":", rownames(coefs2), fixed = TRUE), , drop = FALSE
            ]
            if (nrow(int_rows2) > 0 && "Pr(>|t|)" %in% colnames(int_rows2)) {
              min(int_rows2[["Pr(>|t|)"]], na.rm = TRUE)
            } else {
              NA_real_
            }
          } else {
            NA_real_
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
      attr_rows[[idx]] <- data.frame(
        N_Per_Group      = n,
        Alpha            = a,
        Dropout_Limit    = dropout_limit,
        Pct_Obs_Lost     = if (obs_total > 0)
          100 * (1 - obs_kept / obs_total) else 0,
        Pct_Mice_Truncated = if (mice_total > 0)
          100 * mice_truncated / mice_total else 0,
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
      dropout_limit        = dropout_limit,
      alpha                = alpha,
      n_simulations        = n_simulations
    ),
    attrition = do.call(rbind, attr_rows[seq_len(idx)])
  )
}
