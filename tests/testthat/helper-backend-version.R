# CODE_REVIEW.md R4.D2 — fail loudly on backend version skew.
#
# `devtools::load_all()` resolves `Imports: mouseExperiment` against the INSTALLED
# library, not the sibling source tree. During the Round 3 review the installed
# backend was 0.4.10 while the source under review was 0.9.0 — six releases and the
# whole of Round 3 apart. Every dashboard test that exercises backend behaviour was
# testing code nobody had changed.
#
# Re-running against a freshly-installed backend produced the SAME 210 passing, so
# nothing was hidden that time. That identical result is the real finding: a suite
# whose outcome does not move across six releases of its primary dependency is not
# exercising it, which corroborates §I.1 (`expect_true(TRUE)` placeholders) and
# §I.2 (three largest server modules untested) from a different direction.
#
# This helper makes the skew impossible to miss. It deliberately does NOT skip —
# skipping is the failure mode Round 3 spent its time removing (§R3-L). If the
# installed backend is older than the sibling source, that is a setup error the
# developer must fix, and a test run against the wrong backend is worse than no
# test run at all because it reports green.

.me_sibling_backend_version <- function() {
  # tests/testthat -> tests -> pkg root -> repos root -> mouseExperiment
  candidates <- c(
    file.path("..", "..", "..", "mouseExperiment", "DESCRIPTION"),
    file.path("..", "..", "..", "..", "mouseExperiment", "DESCRIPTION")
  )
  for (p in candidates) {
    if (file.exists(p)) {
      return(list(version = as.character(read.dcf(p, "Version")[1, 1]), path = p))
    }
  }
  NULL
}
