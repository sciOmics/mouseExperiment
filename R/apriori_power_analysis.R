#' A Priori Power Analysis (Analytical)
#'
#' Computes prospective power / required sample size from user-supplied effect
#' size parameters, without needing experimental data. Power is for a
#' two-sample t-test between a treated arm and the control arm. With two groups
#' that is the one comparison the study makes. With \eqn{k \ge 3} groups the
#' study is analysed as the \eqn{k - 1} treated-vs-control comparisons, so the
#' power reported is for each of those comparisons at the per-comparison alpha
#' (see \code{p_adjust_method}).
#'
#' @param effect_size Numeric scalar: Cohen's d, the standardised difference to
#'   detect between a treated arm and the control arm. If \code{delta} and
#'   \code{pooled_sd} are supplied instead, \code{effect_size} is computed as
#'   \code{abs(delta) / pooled_sd}.
#'
#'   \strong{Three or more groups} (CODE_REVIEW.md R20.7). Each treated arm
#'   that differs from control by \code{d} is detected with the reported power;
#'   the chance of detecting several such arms at once is lower. Versions before
#'   0.24.0 powered an omnibus ANOVA with \eqn{f = d / \sqrt{2}} (the conversion
#'   for its own two-extreme-groups configuration is \eqn{d / \sqrt{2k}}), which
#'   returned sample sizes 2.6--4 times too small, and applied the Bonferroni
#'   per-comparison alpha to that omnibus test.
#'
#'   \strong{Effect-size scale:} Cohen's d here is the standardised mean
#'   difference on the \emph{modelling scale} (\code{log(Volume)} when the
#'   downstream LMM uses \code{transform = "log"}, the default in this
#'   package). A user thinking "d = 0.5 corresponds to a 0.5 mm^3 difference"
#'   will badly mis-specify the analysis — log-scale d = 0.5 corresponds to
#'   approximately a \code{exp(0.5 * sigma_log)} fold-difference between
#'   group means, where \code{sigma_log} is the within-group SD of
#'   log-volume (often ~0.3-0.5 in preclinical TG data, so d = 0.5
#'   ≈ 1.2-1.3× fold-difference). When in doubt, compute d from pilot data
#'   on the log scale and pass \code{delta} + \code{pooled_sd} explicitly.
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
#' @param n_comparisons Integer. Number of treatment-vs-control comparisons the
#'   analysis will make. Used with \code{p_adjust_method = "bonferroni"} to
#'   derive the effective per-comparison alpha. Defaults to \code{n_groups - 1}
#'   (the vs-control family these designs normally use).
#' @param p_adjust_method \code{"bonferroni"} (default) or \code{"none"}.
#'   \code{"bonferroni"} computes power at \code{alpha / n_comparisons}, which
#'   matches an analysis that adjusts over the \eqn{k - 1} vs-control
#'   comparisons, as \code{tumor_growth_statistics()} does by default. With two
#'   groups there is one comparison, so it changes nothing. Powering at an
#'   unadjusted alpha and then analysing with an adjustment delivers materially
#'   less power than the nominal target. Dunnett's exact correction is slightly
#'   less conservative than Bonferroni, so for a Dunnett analysis the Bonferroni
#'   N is a small overestimate. The default was \code{"none"} before 0.24.0.
#' @param dropout_rate Numeric in [0, 1). Expected proportion of enrolled
#'   animals that will not be analysable (euthanasia for tumour burden,
#'   technical failure). \code{Required_N} is the number that must be
#'   analysable; \code{Enroll_N = ceiling(Required_N / (1 - dropout_rate))} is
#'   the number to enrol. Default 0 for backward compatibility, but in
#'   preclinical oncology it is rarely 0.
#' @export
apriori_power_analysis <- function(effect_size   = NULL,
                                   delta         = NULL,
                                   pooled_sd     = NULL,
                                   n_groups      = 2L,
                                   alpha         = c(0.05, 0.01),
                                   target_power  = c(0.80, 0.90),
                                   n_per_group   = NULL,
                                   mode          = c("find_n", "find_power"),
                                   n_comparisons = NULL,
                                   p_adjust_method = c("bonferroni", "none"),
                                   dropout_rate  = 0,
                                   sd_sensitivity_range = 0.4) {

  mode <- match.arg(mode)
  p_adjust_method <- match.arg(p_adjust_method)

  # CODE_REVIEW.md R3.15 — two omissions made Required_N an underestimate of
  # what a study needs.
  #
  # (1) Multiplicity. These designs are analysed as k-1 vs-control contrasts
  #     with an adjustment (the package's own default is Bonferroni), but the
  #     power calculation used `alpha` as given. A 4-arm study powered at
  #     alpha = 0.05 is actually analysed at an effective 0.0167 per
  #     comparison, so the delivered power falls well short of the target.
  # (2) Attrition. Required_N is the number of *analysable* animals. Animals
  #     are lost to euthanasia and technical failure, differentially by arm
  #     (controls first), so the number to enrol is larger than the number the
  #     calculation assumed would complete.
  if (!is.numeric(dropout_rate) || length(dropout_rate) != 1L ||
      dropout_rate < 0 || dropout_rate >= 1) {
    stop("'dropout_rate' must be a single number in [0, 1).", call. = FALSE)
  }
  if (!is.null(n_comparisons)) {
    n_comparisons <- as.integer(n_comparisons)
    if (is.na(n_comparisons) || n_comparisons < 1L) {
      stop("'n_comparisons' must be a positive integer.", call. = FALSE)
    }
  }
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

  # Effective per-comparison alpha. Defaults to the k-1 vs-control family when
  # the user does not say otherwise, which is what these studies actually run.
  n_comp_used <- if (!is.null(n_comparisons)) {
    n_comparisons
  } else if (p_adjust_method == "none") {
    1L
  } else {
    max(1L, n_groups - 1L)
  }
  alpha_requested <- alpha
  alpha_effective <- if (p_adjust_method == "bonferroni") {
    alpha / n_comp_used
  } else {
    alpha
  }

  # ---- Power for one treated-vs-control comparison --------------------------
  # CODE_REVIEW.md R20.7 -- k enters only through the per-comparison alpha: a
  # k-arm study is analysed as k - 1 treated-vs-control comparisons, each a
  # two-sample contrast. The former k >= 3 path powered an omnibus F-test with
  # f = d / sqrt(2) (correct for its own configuration is d / sqrt(2k)) and
  # returned N 2.6-4x too small; it also applied the per-comparison Bonferroni
  # alpha to that omnibus test.
  power_fn <- function(n, a) {
    tryCatch(
      stats::power.t.test(n = n, delta = effect_size, sd = 1,
                          sig.level = a, type = "two.sample")$power,
      error = function(e) NA_real_
    )
  }

  n_fn <- function(a, pwr, d = effect_size) {
    tryCatch(
      ceiling(stats::power.t.test(power = pwr, delta = d, sd = 1,
                                  sig.level = a, type = "two.sample")$n),
      error = function(e) NA_integer_
    )
  }

  if (mode == "find_n") {
    rows <- lapply(seq_along(alpha), function(i) {
      a_req <- alpha_requested[i]
      a_eff <- alpha_effective[i]
      lapply(target_power, function(tp) {
        n_analysable <- n_fn(a_eff, tp)
        data.frame(Effect_Size = effect_size, N_Groups = n_groups,
                   Alpha = a_req,
                   Alpha_Per_Comparison = a_eff,
                   N_Comparisons = n_comp_used,
                   Target_Power = tp,
                   # The number that must be analysable to hit the target...
                   Required_N = n_analysable,
                   # ...and the number to actually enrol to end up with it.
                   Enroll_N = ceiling(n_analysable / (1 - dropout_rate)),
                   Dropout_Rate = dropout_rate,
                   stringsAsFactors = FALSE)
      })
    })
    scenario_table <- do.call(rbind, unlist(rows, recursive = FALSE))
  } else {
    if (is.null(n_per_group) || length(n_per_group) == 0)
      stop("'n_per_group' must be supplied when mode = 'find_power'.")
    rows <- lapply(as.integer(n_per_group), function(n) {
      lapply(seq_along(alpha), function(i) {
        a <- alpha_effective[i]; a_req <- alpha_requested[i]
        # n is what the user has; attrition reduces it to n_analysable.
        n_analysable <- max(2L, floor(n * (1 - dropout_rate)))
        data.frame(Effect_Size = effect_size, N_Groups = n_groups,
                   Alpha = a_req,
                   Alpha_Per_Comparison = a,
                   N_Comparisons = n_comp_used,
                   N_Per_Group = n,
                   N_Analysable = n_analysable,
                   Dropout_Rate = dropout_rate,
                   Achieved_Power = power_fn(n_analysable, a),
                   stringsAsFactors = FALSE)
      })
    })
    scenario_table <- do.call(rbind, unlist(rows, recursive = FALSE))
  }

  # ---- Power curve data ------------------------------------------------------
  n_seq <- seq(2L, 80L, by = 1L)
  curve_rows <- lapply(seq_along(alpha), function(i) {
    a_eff <- alpha_effective[i]
    data.frame(
      Effect_Size = effect_size,
      N_Per_Group = n_seq,
      Power       = vapply(n_seq, function(n) power_fn(n, a_eff), numeric(1)),
      Alpha       = alpha_requested[i],
      Alpha_Per_Comparison = a_eff,
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
    # The same per-comparison alpha as the scenario table (R20.7 / R20.71).
    ref_i  <- which.min(alpha)
    ref_a  <- alpha_requested[ref_i]
    ref_a_eff <- alpha_effective[ref_i]
    ref_tp <- max(target_power)
    sens_rows <- lapply(seq_along(sd_factors), function(i) {
      sd_i  <- stored_pooled_sd * sd_factors[i]
      d_i   <- raw_delta / sd_i
      n_i   <- n_fn(ref_a_eff, ref_tp, d = d_i)
      data.frame(SD_Change = sd_labels[i], Assumed_SD = round(sd_i, 3),
                 Cohens_d = round(d_i, 3), Required_N = n_i,
                 Alpha = ref_a, Alpha_Per_Comparison = ref_a_eff,
                 Target_Power = ref_tp,
                 stringsAsFactors = FALSE)
    })
    sensitivity_table <- do.call(rbind, sens_rows)
  }

  method_note <- if (n_groups >= 3L) {
    paste0(
      "Power is for one treated-vs-control comparison (two-sample t-test) at ",
      "alpha ", paste(signif(alpha_effective, 3), collapse = " / "),
      " per comparison",
      if (p_adjust_method == "bonferroni") {
        paste0(" (Bonferroni over ", n_comp_used, " comparisons)")
      } else " (no multiplicity adjustment)",
      ". Each treated arm that differs from control by d = ",
      signif(effect_size, 3), " is detected with this power; the chance of ",
      "detecting several such arms at once is lower."
    )
  } else NULL

  list(
    scenario_table    = scenario_table,
    power_curve_data  = power_curve_data,
    sensitivity_table = sensitivity_table,
    effect_size       = effect_size,
    pooled_sd         = stored_pooled_sd,
    n_groups          = n_groups,
    mode              = mode,
    method_note       = method_note
  )
}
