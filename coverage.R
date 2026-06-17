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
# Known caveat (2026-06-02): the test files `test-post_power_analysis.R` and
# parts of `test-toxicity_functions.R` reference removed / renamed functions
# (`post_power_analysis()` was deleted in v0.3.4; `efficacy_toxicity_bivariate`
# no longer accepts the `"log_cell_kill"` metric). Until those stale tests
# are fixed under CODE_REVIEW.md K.11, `covr::package_coverage()` exits with
# a test-failure error before producing a summary. Either:
#   (a) Fix the stale tests first (recommended), or
#   (b) Use `covr::file_coverage()` on a per-file basis to bypass the test
#       harness entirely.

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
