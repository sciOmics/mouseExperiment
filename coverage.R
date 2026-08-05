# Coverage baseline runner (CODE_REVIEW.md Round 2 K.10).
#
# Usage from package root:
#   Rscript coverage.R           # non-Bayesian surface (fast)
#   Rscript coverage.R --bayes   # include the Bayesian paths (slow: Stan compiles)
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
suppressMessages({
  library(covr); library(testthat); library(mouseExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
want_html     <- "--html"  %in% args
# Bayesian coverage requires compiling every Stan model -- slow, but the only
# honest way to quote an overall figure. Off by default.
include_bayes <- "--bayes" %in% args

# `covr::package_coverage(function_exclusions = ...)` fails inside covr's own
# exclude() machinery on this package ("operations are possible only for numeric,
# logical or complex types", from Ops.data.frame on the exclusion list). Use
# file_coverage(), which pairs source files to test files directly and bypasses
# that path entirely.
# ALL source files must be listed: file_coverage() only makes the listed files
# available, and the tests call package internals across file boundaries
# (make_mouse_key, synergy_bliss_expected, me_result_meta, ...). Excluding a
# source file therefore breaks the tests rather than merely omitting it from the
# report. Only the TEST files are filtered.
src_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
tst_files <- list.files("tests/testthat", pattern = "^test-.*[.]R$",
                        full.names = TRUE)

if (!include_bayes) {
  tst_files <- tst_files[!grepl("/test-bayesian", tst_files)]
}

# file_coverage() sys.source()s each test file directly; unlike test_dir() it does
# NOT load helper-*.R first, so every fixture (make_tg_simple, make_bw_simple, ...)
# would be missing. Prepend them.
helper_files <- list.files("tests/testthat", pattern = "^helper-.*[.]R$",
                           full.names = TRUE)
tst_files <- c(helper_files, tst_files)

cov <- covr::file_coverage(source_files = src_files, test_files = tst_files)

cat(sprintf("\n%s LINE COVERAGE: %.1f%%\n",
            if (include_bayes) "OVERALL" else "NON-BAYESIAN TESTS ONLY",
            covr::percent_coverage(cov)))
if (!include_bayes) {
  cat("NOTE: the bayesian_*.R files are measured but their tests were not run,\n",
      "      so their low figures are an artefact of this lane, not a finding.\n",
      "      Use --bayes for a number you can quote.\n", sep = "")
}
print(cov)

if (want_html) {
  out <- normalizePath("covr_report.html", mustWork = FALSE)
  covr::report(cov, file = out, browse = FALSE)
  message("Wrote HTML report: ", out)
}
