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


#' Loewe combination index (single-dose approximation)
#'
#' Computes the single-dose Loewe combination index:
#' \deqn{ \mathrm{CI}_{\mathrm{Loewe}} =
#'        \frac{\min(\mathrm{FE}_A + \mathrm{FE}_B,\ 1)}
#'             {\max(\mathrm{FE}_{\mathrm{combo}},\ \mathrm{floor})} }
#' \code{CI < 1} indicates synergy; \code{CI > 1} indicates antagonism.
#' The denominator is floored at \code{fe_floor} (default \code{1e-4}) to
#' avoid pathological CI values when combo efficacy is near zero.
#' Returns a list of \code{ci} (the index, same length as the inputs) and
#' \code{floor_applied} (logical of equal length: TRUE where the floor
#' was active for that element). Vectorised.
#'
#' This is the *single-dose* form of the Loewe additivity criterion, which
#' assumes a linear dose-response relationship between fraction-of-effect
#' and dose. Use with caution when actual dose-response curvature is
#' substantial; consider \code{drc::isobole()} or a full Loewe surface for
#' more rigorous combination analysis.
#' @noRd
synergy_loewe_ci <- function(fe_a, fe_b, fe_combo, fe_floor = 1e-4) {
  num            <- pmin(fe_a + fe_b, 1)
  denom_floored  <- pmax(fe_combo, fe_floor)
  ci             <- num / denom_floored
  list(ci = ci, floor_applied = fe_combo < fe_floor, loewe_num = num)
}
