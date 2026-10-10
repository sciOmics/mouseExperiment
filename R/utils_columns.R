# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

#' Quote a column name for use in a pasted formula
#'
#' CODE_REVIEW.md R20.43 -- formulas built with paste() broke on names such as
#' "Study Day" ("unexpected symbol"). Backticks make any name safe for lm(),
#' aov(), drm() and survival formulas. brms is stricter (syntactic names, no
#' double underscores), so the Bayesian entry points copy their columns into
#' fixed internal names instead.
#'
#' @param x Character vector of column names.
#' @return The names wrapped in backticks.
#' @noRd
#' @keywords internal
me_bt <- function(x) paste0("`", x, "`")
