#' A Priori Power Analysis (Analytical)
#'
#' Computes prospective power / required sample size from user-supplied effect
#' size parameters, without needing experimental data. Supports two-sample
#' t-test (two groups) and one-way ANOVA omnibus F-test (k ≥ 3 groups).
#'
#' @param effect_size Numeric scalar. Cohen's d (for two groups) or Cohen's f
#'   (for ANOVA; f = sd_means / pooled_within_sd). If \code{delta} and
#'   \code{pooled_sd} are supplied instead, \code{effect_size} is computed as
#'   \code{abs(delta) / pooled_sd} and the two-group t-test path is used.
#' @param delta Numeric scalar. Raw mean difference between groups. Used only
#'   when \code{pooled_sd} is also supplied. Ignored if \code{effect_size} is
#'   provided directly.
#' @param pooled_sd Numeric scalar > 0. Pooled within-group standard deviation.
#'   Used to derive Cohen's d from \code{delta}. Ignored if \code{effect_size}
#'   is provided.
#' @param n_groups Integer ≥ 2. Number of treatment groups (default 2).
#' @param alpha Numeric vector of significance levels (default \code{c(0.05,
#'   0.01)}).
#' @param target_power Numeric vector of desired power levels (default
#'   \code{c(0.80, 0.90)}).
#' @param n_per_group Integer vector or \code{NULL}. If provided, power is
#'   computed at these sample sizes instead of solving for required N.
#' @param mode Character. \code{"find_n"} (default) — for each (alpha,
#'   target_power) combination return the N per group required to achieve it.
#'   \code{"find_power"} — requires \code{n_per_group} to be set; returns
#'   achieved power at those N values.
#' @param sd_sensitivity_range Numeric scalar 0–1 (default 0.4). If > 0, a
#'   sensitivity table is appended showing how required N changes when the
#'   assumed SD is varied by ±(sensitivity_range × 100)% and
#'   ±(sensitivity_range/2 × 100)%.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{scenario_table}}{Data frame: Effect_Size, N_Groups, Alpha,
#'     Target_Power / N_Per_Group, Required_N / Achieved_Power (depending on
#'     mode).}
#'   \item{\code{power_curve_data}}{Data frame suitable for plotting power
#'     vs N per group for each alpha level.}
#'   \item{\code{sensitivity_table}}{Data frame showing required N across SD
#'     perturbations (only when \code{sd_sensitivity_range > 0} and
#'     \code{pooled_sd} is known).}
#'   \item{\code{effect_size}}{The Cohen's d / f used.}
#'   \item{\code{pooled_sd}}{The pooled SD used (NA when not applicable).}
#'   \item{\code{n_groups}}{Number of groups.}
#'   \item{\code{mode}}{Analysis mode.}
#' }
#'
#' @export
apriori_power_analysis <- function(effect_size   = NULL,
                                   delta         = NULL,
                                   pooled_sd     = NULL,
                                   n_groups      = 2L,
                                   alpha         = c(0.05, 0.01),
                                   target_power  = c(0.80, 0.90),
                                   n_per_group   = NULL,
                                   mode          = c("find_n", "find_power"),
                                   sd_sensitivity_range = 0.4) {

  mode <- match.arg(mode)
  n_groups <- as.integer(n_groups)
  if (n_groups < 2L) stop("n_groups must be >= 2")

  # ---- Resolve effect size ---------------------------------------------------
  if (is.null(effect_size)) {
    if (is.null(delta) || is.null(pooled_sd))
      stop("Provide either 'effect_size' (Cohen's d/f) or both 'delta' and 'pooled_sd'.")
    if (!is.numeric(pooled_sd) || pooled_sd <= 0)
      stop("'pooled_sd' must be a positive number.")
    effect_size <- abs(as.numeric(delta)) / as.numeric(pooled_sd)
  }
  effect_size <- as.numeric(effect_size)
  if (length(effect_size) != 1 || !is.finite(effect_size) || effect_size <= 0)
    stop("'effect_size' must be a single positive finite number.")

  stored_pooled_sd <- if (!is.null(pooled_sd)) as.numeric(pooled_sd) else NA_real_

  # ---- Build scenario table --------------------------------------------------
  power_fn <- if (n_groups == 2L) {
    # Two-sample t-test
    function(n, a) {
      tryCatch(
        stats::power.t.test(n = n, delta = effect_size, sd = 1,
                            sig.level = a, type = "two.sample")$power,
        error = function(e) NA_real_
      )
    }
  } else {
    # One-way ANOVA: convert d → f = d / sqrt(2)  (equal group means, two
    # extreme groups), then use pwr::pwr.anova.test if available, otherwise
    # approximate via non-central F.
    f_val <- effect_size / sqrt(2)
    function(n, a) {
      if (requireNamespace("pwr", quietly = TRUE)) {
        tryCatch(
          pwr::pwr.anova.test(k = n_groups, n = n, f = f_val,
                              sig.level = a)$power,
          error = function(e) NA_real_
        )
      } else {
        # Non-central F approximation
        tryCatch({
          df1 <- n_groups - 1L
          df2 <- n_groups * (n - 1L)
          ncp  <- n * n_groups * f_val^2
          crit <- stats::qf(1 - a, df1, df2)
          1 - stats::pf(crit, df1, df2, ncp = ncp)
        }, error = function(e) NA_real_)
      }
    }
  }

  n_fn <- function(a, pwr) {
    if (n_groups == 2L) {
      tryCatch(
        ceiling(stats::power.t.test(power = pwr, delta = effect_size, sd = 1,
                                    sig.level = a, type = "two.sample")$n),
        error = function(e) NA_integer_
      )
    } else {
      f_val <- effect_size / sqrt(2)
      if (requireNamespace("pwr", quietly = TRUE)) {
        tryCatch(
          ceiling(pwr::pwr.anova.test(k = n_groups, f = f_val,
                                      sig.level = a, power = pwr)$n),
          error = function(e) NA_integer_
        )
      } else {
        # Binary search fallback
        lo <- 2L; hi <- 10000L
        while (lo < hi) {
          mid <- (lo + hi) %/% 2L
          if (isTRUE(power_fn(mid, a) >= pwr)) hi <- mid else lo <- mid + 1L
        }
        lo
      }
    }
  }

  if (mode == "find_n") {
    rows <- lapply(alpha, function(a) {
      lapply(target_power, function(tp) {
        data.frame(Effect_Size = effect_size, N_Groups = n_groups,
                   Alpha = a, Target_Power = tp,
                   Required_N = n_fn(a, tp),
                   stringsAsFactors = FALSE)
      })
    })
    scenario_table <- do.call(rbind, unlist(rows, recursive = FALSE))
  } else {
    if (is.null(n_per_group) || length(n_per_group) == 0)
      stop("'n_per_group' must be supplied when mode = 'find_power'.")
    rows <- lapply(as.integer(n_per_group), function(n) {
      lapply(alpha, function(a) {
        data.frame(Effect_Size = effect_size, N_Groups = n_groups,
                   Alpha = a, N_Per_Group = n,
                   Achieved_Power = power_fn(n, a),
                   stringsAsFactors = FALSE)
      })
    })
    scenario_table <- do.call(rbind, unlist(rows, recursive = FALSE))
  }

  # ---- Power curve data ------------------------------------------------------
  n_seq <- seq(2L, 80L, by = 1L)
  curve_rows <- lapply(alpha, function(a) {
    data.frame(
      N_Per_Group = n_seq,
      Power       = vapply(n_seq, function(n) power_fn(n, a), numeric(1)),
      Alpha       = a,
      stringsAsFactors = FALSE
    )
  })
  power_curve_data <- do.call(rbind, curve_rows)

  # ---- Sensitivity table (SD perturbations) ---------------------------------
  sensitivity_table <- NULL
  if (sd_sensitivity_range > 0 && !is.na(stored_pooled_sd) && !is.null(delta)) {
    raw_delta  <- abs(as.numeric(delta))
    sd_factors <- c(1 - sd_sensitivity_range,
                    1 - sd_sensitivity_range / 2,
                    1,
                    1 + sd_sensitivity_range / 2,
                    1 + sd_sensitivity_range)
    sd_labels  <- paste0(round((sd_factors - 1) * 100), "%")
    ref_a  <- min(alpha)
    ref_tp <- max(target_power)
    sens_rows <- lapply(seq_along(sd_factors), function(i) {
      sd_i  <- stored_pooled_sd * sd_factors[i]
      d_i   <- raw_delta / sd_i
      n_i   <- if (n_groups == 2L) {
        tryCatch(ceiling(stats::power.t.test(power = ref_tp, delta = d_i, sd = 1,
                                             sig.level = ref_a, type = "two.sample")$n),
                 error = function(e) NA_integer_)
      } else {
        n_fn_local <- function(a2, pwr2) {
          f2 <- d_i / sqrt(2)
          if (requireNamespace("pwr", quietly = TRUE))
            tryCatch(ceiling(pwr::pwr.anova.test(k = n_groups, f = f2,
                                                 sig.level = a2, power = pwr2)$n),
                     error = function(e) NA_integer_)
          else NA_integer_
        }
        n_fn_local(ref_a, ref_tp)
      }
      data.frame(SD_Change = sd_labels[i], Assumed_SD = round(sd_i, 3),
                 Cohens_d = round(d_i, 3), Required_N = n_i,
                 Alpha = ref_a, Target_Power = ref_tp,
                 stringsAsFactors = FALSE)
    })
    sensitivity_table <- do.call(rbind, sens_rows)
  }

  list(
    scenario_table    = scenario_table,
    power_curve_data  = power_curve_data,
    sensitivity_table = sensitivity_table,
    effect_size       = effect_size,
    pooled_sd         = stored_pooled_sd,
    n_groups          = n_groups,
    mode              = mode
  )
}
