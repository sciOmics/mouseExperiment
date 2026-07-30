# Coverage baseline runner (CODE_REVIEW.md Round 2 K.10).
#
# Usage from package root:
#   Rscript coverage.R           # print per-file coverage summary to stdout
#   Rscript coverage.R --html    # also write HTML report (covr_report.html)
#
# Output is a `covr::package_coverage()` result with line-, branch-, and
# function-level coverage rolled up per file plus the overall total.
# Functions with Bayesian fits are excluded by default because `brms::brm()`
# requires Stan compilation that won't fit in a quick coverage run — they're
# present in the test suite but the long path isn't exercised on each run.
# Drop the `function_exclusions` argument to measure them too (slow).
#
# CI tip: piping the printed result through `grep "Coverage:"` gives a
# single-line summary suitable for a build status check.
#
# The 2026-06-02 caveat is resolved: the stale tests that blocked this script
# (`test-post_power_analysis.R`, and the `"log_cell_kill"` assertions in
# `test-toxicity_functions.R`) were removed or repaired in v0.5.0, so
# `covr::package_coverage()` now produces a summary.
#
# Read the Bayesian exclusion carefully. Those seven functions are the package's
# heaviest and most defect-prone surface -- CODE_REVIEW.md R3-L found two Critical
# bugs in them that had survived five releases -- so a headline number computed
# with them excluded flatters the package. Treat the figure this script prints as
# "coverage of the non-Bayesian surface", and run without `function_exclusions`
# (slow: every Stan model compiles) before making any claim about overall
# coverage.

stopifnot(requireNamespace("covr", quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
want_html <- "--html" %in% args

cov <- covr::package_coverage(
  type                = "tests",
  function_exclusions = c(
    "^bayesian_tumor_growth$",
    "^bayesian_body_weight$",
    "^bayesian_survival$",
    "^bayesian_synergy(_over_time)?$",
    "^bayesian_therapeutic_window$",
    "^bayesian_dose_response$",
    "^bayesian_power_analysis$"
  ),
  quiet = FALSE
)

print(cov)

if (want_html) {
  out <- normalizePath("covr_report.html", mustWork = FALSE)
  covr::report(cov, file = out, browse = FALSE)
  message("Wrote HTML report: ", out)
}
