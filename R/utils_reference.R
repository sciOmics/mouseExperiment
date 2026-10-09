# Choosing the reference (control) arm -- CODE_REVIEW.md R20.46
#
# Tumour growth, survival, body weight and both Bayesian entry points used to
# default to the alphabetically first arm. That choice depends on the locale,
# and it is often a treated arm: the shipped combination data gave "aPD1"
# under en_US. The toxicity functions matched six exact names, so
# "Vehicle (10% DMSO)" was not recognised, and a reference_group that matched
# no arm was silently replaced by another.
#
# Decided 2026-10-09: one shared detector. If no control is recognised, stop
# and ask for reference_group. A named reference must exist. The arm used is
# always reported.

# Control names, matched as whole words anywhere in an arm's name,
# case-insensitively ("Vehicle (10% DMSO)", "IgG2a isotype", "PBS control").
me_control_pattern <- paste0(
  "\\b(control|ctrl|vehicle|dmso|pbs|saline|untreated|placebo|isotype|",
  "igg[0-9a-z]*)\\b")

#' Resolve the reference arm for comparisons
#'
#' @param arms Treatment labels (any vector; NAs ignored).
#' @param reference_group The caller's reference_group, or `NULL` to detect one.
#' @param fn Name of the calling function, for messages.
#' @param quiet Suppress the "using ... as the reference" message.
#' @return The reference arm, a single string.
#' @keywords internal
#' @noRd
me_resolve_reference <- function(arms, reference_group = NULL, fn = "analysis",
                                 quiet = FALSE) {
  arms <- unique(as.character(arms[!is.na(arms)]))
  shown <- paste0("\"", arms, "\"", collapse = ", ")

  if (!is.null(reference_group)) {
    reference_group <- as.character(reference_group)
    if (length(reference_group) != 1L || is.na(reference_group) ||
        !reference_group %in% arms) {
      stop(fn, "(): reference_group \"", paste(reference_group, collapse = ", "),
           "\" is not one of the treatment arms: ", shown, ".", call. = FALSE)
    }
    return(reference_group)
  }

  hits <- arms[grepl(me_control_pattern, arms, ignore.case = TRUE, perl = TRUE)]
  if (length(hits) > 1L) {
    # "Vehicle" over "Drug A + Vehicle": a combination is not the control.
    single <- hits[!grepl("\\+|&|\\band\\b|\\bwith\\b", hits, ignore.case = TRUE)]
    if (length(single) >= 1L) hits <- single
  }
  if (length(hits) == 0L) {
    stop(fn, "(): no control arm recognised among ", shown,
         ". Set reference_group to the arm comparisons should be made against.",
         call. = FALSE)
  }
  if (length(hits) > 1L) {
    stop(fn, "(): several arms look like a control (",
         paste0("\"", hits, "\"", collapse = ", "),
         "). Set reference_group to the one comparisons should be made against.",
         call. = FALSE)
  }
  if (!isTRUE(quiet)) {
    message(fn, "(): using \"", hits, "\" as the reference group (recognised as ",
            "the control). Set reference_group to choose another.")
  }
  hits
}
