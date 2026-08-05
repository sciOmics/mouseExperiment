# Copyright (c) 2026 mouseExperiment Contributors
# Licensed under the MIT License - see LICENSE file

# Canonical result provenance — CODE_REVIEW.md B7.1 / B7.2.
#
# B7.1 asked for the treatment-effects schema to be harmonised between the
# frequentist and Bayesian functions. Renaming the columns to match would be the
# wrong fix: `Lower_CL` is a confidence limit and `Lower_CrI` is a credible
# interval, and Round 1 §B1.4 deliberately separated them *because* they are
# different things. Forcing one name back onto both would trade a discoverability
# problem for a correctness problem.
#
# The actual defect is that a consumer cannot tell which it is holding without
# guessing. The dashboard currently does exactly that:
#
#   intersect(c("lower.CL", "Lower_CrI", "Lower", "lower"), colnames(df))[1]
#
# which is the §K.2 pattern — if a name changes the guess silently returns NA and
# the panel goes blank rather than failing.
#
# So: keep the semantically-correct column names, and add a small uniform block
# that says what the object *is*. Consumers branch on `meta$interval_type` instead
# of sniffing column names.
#
# B7.2 is the same shape. `bayesian_survival()` omits `transform_used` because an
# AFT model applies no volume transform — reasonable, but it leaves a consumer
# unable to distinguish "no transform" from "field absent". Every function now
# reports it, with `"none"` meaning none was applied.
#
# This is purely additive: no existing field is renamed or removed.

#' Build a canonical result-provenance block
#'
#' Attached as `$meta` by every analysis entry point so downstream consumers can
#' interrogate a result object rather than infer its nature from column names.
#'
#' @param analysis_type Human-readable description, e.g. "Linear mixed-effects
#'   model".
#' @param model_type_used Machine-readable model identifier, e.g. `"lme4"`,
#'   `"gam"`, `"auc"`, `"bayes_tg"`.
#' @param inference `"frequentist"` or `"bayesian"`. Determines how the interval
#'   columns should be read and reported.
#' @param interval_type What the interval columns mean: `"confidence"` (columns
#'   suffixed `_CL`), `"credible"` (`_CrI`), or `"none"`.
#' @param transform_used Transform applied to the response on the modelling scale:
#'   `"log"`, `"sqrt"`, or `"none"`. `"none"` means none was applied — never
#'   absent, so a consumer can always tell.
#' @param estimate_scale What the point estimates are on, e.g.
#'   `"log volume"`, `"volume (mm3)"`, `"log time ratio"`. This is what a caller
#'   needs in order to back-transform correctly.
#' @param comparison_family,p_adjust_method Multiplicity provenance, or `NA` when
#'   the function reports no pairwise comparisons.
#' @param extra Optional named list of additional fields.
#' @return A list of class `me_meta`.
#' @noRd
#' @keywords internal
me_result_meta <- function(analysis_type,
                           model_type_used,
                           inference        = c("frequentist", "bayesian"),
                           interval_type    = c("confidence", "credible", "none"),
                           transform_used   = "none",
                           estimate_scale   = NA_character_,
                           comparison_family = NA_character_,
                           p_adjust_method   = NA_character_,
                           extra             = NULL) {
  inference     <- match.arg(inference)
  interval_type <- match.arg(interval_type)

  out <- list(
    analysis_type     = analysis_type,
    model_type_used   = model_type_used,
    inference         = inference,
    interval_type     = interval_type,
    # Column names a consumer should read for the interval, derived from
    # interval_type so nobody has to hard-code the mapping.
    interval_columns  = switch(interval_type,
      confidence = c(lower = "Lower_CL",  upper = "Upper_CL"),
      credible   = c(lower = "Lower_CrI", upper = "Upper_CrI"),
      none       = c(lower = NA_character_, upper = NA_character_)
    ),
    transform_used    = transform_used,
    estimate_scale    = estimate_scale,
    comparison_family = comparison_family,
    p_adjust_method   = p_adjust_method,
    package_version   = as.character(utils::packageVersion("mouseExperiment"))
  )
  if (!is.null(extra)) out <- c(out, extra)
  structure(out, class = c("me_meta", "list"))
}

#' @export
print.me_meta <- function(x, ...) {
  cat("<me_meta>\n")
  cat("  analysis      :", x$analysis_type, "\n")
  cat("  model_type    :", x$model_type_used, "\n")
  cat("  inference     :", x$inference, "\n")
  cat("  intervals     :", x$interval_type,
      if (!is.na(x$interval_columns[["lower"]])) {
        paste0("(", x$interval_columns[["lower"]], " / ",
               x$interval_columns[["upper"]], ")")
      } else "", "\n")
  cat("  transform     :", x$transform_used, "\n")
  if (!is.na(x$estimate_scale)) cat("  estimate scale:", x$estimate_scale, "\n")
  if (!is.na(x$comparison_family)) {
    cat("  comparisons   :", x$comparison_family,
        paste0("(", x$p_adjust_method, ")"), "\n")
  }
  cat("  built by      : mouseExperiment", x$package_version, "\n")
  invisible(x)
}

#' Read a result's interval columns without guessing at their names
#'
#' Replaces the `intersect(c("lower.CL", "Lower_CrI", ...))` sniffing that the
#' dashboard used, which silently yields `NA` when a name changes and blanks the
#' panel rather than failing (§K.2).
#'
#' @param result Any analysis result carrying `$meta`.
#' @param df The data frame within it to read, e.g. `result$treatment_effects`.
#' @return A list with `lower` and `upper` numeric vectors, or `NULL` when the
#'   result reports no intervals.
#' @noRd
#' @keywords internal
me_interval_cols <- function(result, df) {
  meta <- result$meta
  if (is.null(meta) || identical(meta$interval_type, "none")) return(NULL)
  lo <- meta$interval_columns[["lower"]]
  hi <- meta$interval_columns[["upper"]]
  if (is.na(lo) || !all(c(lo, hi) %in% names(df))) {
    warning("`meta` declares interval columns '", lo, "' / '", hi,
            "' but they are not present in the supplied data frame. ",
            "This is a bug in the analysis function, not the caller.",
            call. = FALSE)
    return(NULL)
  }
  list(lower = df[[lo]], upper = df[[hi]])
}
