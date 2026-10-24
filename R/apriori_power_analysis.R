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
#'   analysis will make. Used with \code{p_adjust_method = "bonferroni"} or
#'   \code{"dunnett"} to derive the effective per-comparison alpha. Defaults to \code{n_groups - 1}
#'   (the vs-control family these designs normally use).
#' @param p_adjust_method \code{"bonferroni"} (default), \code{"dunnett"} or
#'   \code{"none"}.
#'   \code{"bonferroni"} computes power at \code{alpha / n_comparisons}, which
#'   matches an analysis that adjusts over the \eqn{k - 1} vs-control
#'   comparisons, as \code{tumor_growth_statistics()} does by default. With two
#'   groups there is one comparison, so it changes nothing. Powering at an
#'   unadjusted alpha and then analysing with an adjustment delivers materially
#'   less power than the nominal target. \code{"dunnett"} (v0.44.0) powers each
#'   comparison at Dunnett's per-comparison alpha: the two-sided critical value
#'   of the \eqn{k - 1} comparisons with a shared control (equal group sizes,
#'   correlation 0.5) on the pooled error degrees of freedom \eqn{k (n - 1)},
#'   so it depends on \eqn{n}. It is slightly less conservative than
#'   Bonferroni, and matches the Tumour Growth default. The default was
#'   \code{"none"} before 0.24.0.
#' @param dropout_rate Numeric in [0, 1). Expected proportion of enrolled
#'   animals that will not be analysable (euthanasia for tumour burden,
#'   technical failure). \code{Required_N} is the number that must be
#'   analysable; \code{Enroll_N = ceiling(Required_N / (1 - dropout_rate))} is
#'   the number to enrol. Default 0 for backward compatibility, but in
#'   preclinical oncology it is rarely 0.
#' @examples
#' # Animals per group for a standardised effect of 1.2 with three groups
#' apriori_power_analysis(effect_size = 1.2, n_groups = 3)
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
                                   p_adjust_method = c("bonferroni", "dunnett", "none"),
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
  # R20.30: n_comparisons was silently ignored without an adjustment, yet
  # echoed in the table as if it had been applied.
  if (!is.null(n_comparisons) && p_adjust_method == "none" && n_comparisons > 1L) {
    warning("n_comparisons = ", n_comparisons, " has no effect with p_adjust_method = ",
            "\"none\": each comparison is tested at alpha = ",
            paste(alpha, collapse = ", "), ".", call. = FALSE)
  }
  n_comp_used <- if (!is.null(n_comparisons)) {
    n_comparisons
  } else if (p_adjust_method == "none") {
    1L
  } else {
    max(1L, n_groups - 1L)
  }
  alpha_requested <- alpha
  # F6: Dunnett's per-comparison alpha depends on n through the error degrees
  # of freedom, so alpha_at() gives it at each n; alpha_effective holds its
  # large-sample value for the notes.
  alpha_effective <- switch(p_adjust_method,
    bonferroni = alpha / n_comp_used,
    dunnett    = vapply(alpha, me_dunnett_alpha, numeric(1), m = n_comp_used),
    none       = alpha)
  alpha_at <- function(i, n) {
    if (p_adjust_method != "dunnett") return(alpha_effective[i])
    me_dunnett_alpha(alpha_requested[i], n_comp_used, df = n_groups * (max(n, 2L) - 1L))
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

  # The analysable n for power `pwr` at the i-th alpha. Under Dunnett the
  # answer lies between the unadjusted and the Bonferroni ones, so those
  # bracket the search.
  n_for <- function(i, pwr, d = effect_size) {
    if (p_adjust_method != "dunnett") return(n_fn(alpha_effective[i], pwr, d))
    lo <- n_fn(alpha_requested[i], pwr, d)
    hi <- n_fn(alpha_requested[i] / n_comp_used, pwr, d)
    if (is.na(lo) || is.na(hi)) return(NA_integer_)
    for (n in seq(max(2L, lo), max(2L, lo, hi))) {
      p <- tryCatch(stats::power.t.test(n = n, delta = d, sd = 1, sig.level = alpha_at(i, n),
                                        type = "two.sample")$power, error = function(e) NA_real_)
      if (isTRUE(p >= pwr)) return(as.integer(n))
    }
    as.integer(max(lo, hi))
  }

  # R20.71: one definition of n. Enrolment from an analysable n rounds before
  # ceiling(): 21 / (1 - 0.3) is 30.000000000000004 in floating point, which
  # ceiling() made 31. Analysable from enrolled: floor(n * (1 - dropout)).
  enroll_n     <- function(n_a) ceiling(round(n_a / (1 - dropout_rate), 8))
  analysable_n <- function(n_e) pmax(2L, floor(round(n_e * (1 - dropout_rate), 8)))

  if (mode == "find_n") {
    rows <- lapply(seq_along(alpha), function(i) {
      a_req <- alpha_requested[i]
      a_eff <- alpha_effective[i]
      lapply(target_power, function(tp) {
        n_analysable <- n_for(i, tp)
        data.frame(Effect_Size = effect_size, N_Groups = n_groups,
                   Alpha = a_req,
                   Alpha_Per_Comparison = if (is.na(n_analysable)) a_eff else alpha_at(i, n_analysable),
                   N_Comparisons = n_comp_used,
                   Target_Power = tp,
                   # The number that must be analysable to hit the target...
                   Required_N = n_analysable,
                   # ...and the number to actually enrol to end up with it.
                   Enroll_N = enroll_n(n_analysable),
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
        a_req <- alpha_requested[i]
        # n is what the user has; attrition reduces it to n_analysable.
        n_analysable <- analysable_n(n)
        a <- alpha_at(i, n_analysable)
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
  # R20.71: N_Per_Group is the number enrolled, as in the find_power table;
  # the curve plotted power at that many analysable animals (0.562 where the
  # table gave 0.406 at 30 % dropout).
  n_seq <- seq(2L, 80L, by = 1L)
  curve_rows <- lapply(seq_along(alpha), function(i) {
    a_eff <- vapply(analysable_n(n_seq), function(n) alpha_at(i, n), numeric(1))
    data.frame(
      Effect_Size = effect_size,
      N_Per_Group = n_seq,
      N_Analysable = analysable_n(n_seq),
      Power       = mapply(power_fn, analysable_n(n_seq), a_eff),
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
      n_i   <- n_for(ref_i, ref_tp, d = d_i)
      # R20.71: with the enrolment that dropout implies, like the scenario table.
      data.frame(SD_Change = sd_labels[i], Assumed_SD = round(sd_i, 3),
                 Cohens_d = round(d_i, 3), Required_N = n_i,
                 Enroll_N = enroll_n(n_i), Dropout_Rate = dropout_rate,
                 Alpha = ref_a,
                 Alpha_Per_Comparison = if (is.na(n_i)) ref_a_eff else alpha_at(ref_i, n_i),
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
      switch(p_adjust_method,
        bonferroni = paste0(" (Bonferroni over ", n_comp_used, " comparisons)"),
        dunnett    = paste0(" (Dunnett over ", n_comp_used, " comparisons with the control; ",
                            "at the study's N the value is a little smaller)"),
        none       = " (no multiplicity adjustment)"),
      ". Each treated arm that differs from control by d = ",
      signif(effect_size, 3), " is detected with this power; the chance of ",
      "detecting several such arms at once is lower."
    )
  } else NULL

  me_result(list(
    scenario_table    = scenario_table,
    power_curve_data  = power_curve_data,
    sensitivity_table = sensitivity_table,
    effect_size       = effect_size,
    pooled_sd         = stored_pooled_sd,
    n_groups          = n_groups,
    mode              = mode,
    method_note       = method_note
  ), "me_power")
}
