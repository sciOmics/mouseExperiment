# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Shared synergy formulas used by analyze_drug_synergy() (scalar form on
# point estimates) and bayesian_synergy() / bayesian_synergy_over_time()
# (vectorised over posterior draws). Both functions previously had local
# copies of the same algebra; this file consolidates them. The formulas
# work elementwise on numeric vectors of any length, so the same
# implementation serves both call sites.

#' Bliss-independence expected fraction-of-effect
#'
#' Computes the Bliss-independence expected combined effect from two
#' individual fraction-of-effect values:
#' \deqn{ \mathrm{FE}_{\mathrm{Bliss}} = \mathrm{FE}_A + \mathrm{FE}_B
#'        - \mathrm{FE}_A \cdot \mathrm{FE}_B }
#' Bliss assumes the two drugs act independently. Excess over Bliss
#' (\code{fe_combo - synergy_bliss_expected(...)}) is positive for
#' synergy and negative for antagonism. Vectorised — works on scalars
#' or numeric vectors of any length.
#' @noRd
synergy_bliss_expected <- function(fe_a, fe_b) {
  fe_a + fe_b - fe_a * fe_b
}


#' Mouse-level bootstrap CIs for the synergy metrics
#'
#' CODE_REVIEW.md R3.6 / R3.7 / G.6 Option A — `analyze_drug_synergy()` reduced
#' each arm to a single group mean and then made a categorical scientific call
#' ("Strong Synergy" / "Additive" / "Antagonistic") from hard thresholds, with no
#' standard error, interval, or test anywhere in the path. At the group sizes
#' these studies use, the sampling SD of a single TGI is easily 10 percentage
#' points, so the threshold crossings were well inside noise.
#'
#' Resampling **mice within group** — including the control arm — and recomputing
#' the whole statistic per resample fixes two findings at once: it yields a
#' percentile interval for each metric (R3.7), and because the control mean is
#' recomputed on every resample it propagates the uncertainty of the TGI
#' denominator that was previously treated as a known constant (R3.6).
#'
#' Note the resampling unit is the mouse, not the observation: the arms are
#' independent groups of animals, and the mouse is the experimental unit.
#'
#' @param vols Named list of numeric vectors — per-mouse volumes at the
#'   evaluation timepoint, one element per arm
#'   (`control`, `drug_a`, `drug_b`, `combo`).
#' @param n_boot Number of resamples.
#' @param seed Optional RNG seed; the caller's RNG state is restored on exit.
#' @param ci_floor_frac Relative floor applied to the combination fractional
#' @return A data frame with one row per metric and columns `Metric`,
#'   `Estimate`, `CI_Lower`, `CI_Upper`, `n_boot`, or `NULL` when any arm has
#'   fewer than two mice.
#' @noRd
#' @keywords internal
synergy_bootstrap <- function(vols, n_boot = 2000L, seed = NULL,
                              ci_floor_frac = 1e-4) {
  need <- c("control", "drug_a", "drug_b", "combo")
  if (!all(need %in% names(vols))) return(NULL)
  vols <- lapply(vols[need], function(v) v[is.finite(v)])
  if (any(vapply(vols, length, integer(1)) < 2L)) return(NULL)
  if (n_boot < 2L) return(NULL)

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }

  # One resample: draw mice with replacement within each arm, then recompute
  # every derived quantity from scratch.
  one <- function() {
    m <- vapply(vols, function(v) mean(sample(v, length(v), replace = TRUE)),
                numeric(1))
    ctrl <- m[["control"]]
    if (!is.finite(ctrl) || ctrl <= 0) return(rep(NA_real_, 4L))

    fe_a     <- 1 - m[["drug_a"]] / ctrl
    fe_b     <- 1 - m[["drug_b"]] / ctrl
    fe_combo <- 1 - m[["combo"]]  / ctrl

    # R17.2: the Loewe_CI draw was dropped with the rest of the response-additivity
    # path. Bliss excess is the quantity the verdict now rests on, so it is the
    # one that needs an interval.
    bliss_exp <- synergy_bliss_expected(fe_a, fe_b)
    c(fe_a * 100, fe_b * 100, fe_combo * 100,
      fe_combo - bliss_exp)
  }

  draws <- vapply(seq_len(n_boot), function(i) one(), numeric(4L))
  metrics <- c("TGI_A_pct", "TGI_B_pct", "TGI_Combo_pct",
               "Bliss_Excess_FE")

  do.call(rbind, lapply(seq_along(metrics), function(i) {
    d <- draws[i, ]
    d <- d[is.finite(d)]
    if (length(d) < 2L) {
      return(data.frame(Metric = metrics[i], Estimate = NA_real_,
                        CI_Lower = NA_real_, CI_Upper = NA_real_,
                        n_boot = 0L, stringsAsFactors = FALSE))
    }
    q <- stats::quantile(d, c(0.025, 0.975), names = FALSE)
    data.frame(Metric = metrics[i], Estimate = stats::median(d),
               CI_Lower = q[1], CI_Upper = q[2], n_boot = length(d),
               stringsAsFactors = FALSE)
  }))
}
