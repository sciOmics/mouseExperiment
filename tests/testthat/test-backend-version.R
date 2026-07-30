# CODE_REVIEW.md R4.D2 — the installed-backend / sibling-source skew check.
# The helper that locates the sibling DESCRIPTION lives in
# helper-backend-version.R; the assertion lives here so it runs as a real test
# rather than at source time.

test_that("R4.D2: the installed backend is not older than the sibling source", {
  installed <- tryCatch(utils::packageVersion("mouseExperiment"),
                        error = function(e) NULL)
  expect_false(is.null(installed),
               info = "mouseExperiment is not installed at all.")

  sibling <- .me_sibling_backend_version()
  if (is.null(sibling)) {
    # Not a side-by-side checkout (e.g. CI installing from GitHub). Nothing to
    # compare against, so only assert the DESCRIPTION floor is satisfied — which
    # load_all() has already enforced by this point.
    succeed()
    return(invisible(NULL))
  }

  src <- package_version(sibling$version)
  if (installed < src) {
    fail(paste0(
      "Installed mouseExperiment is ", installed, " but the sibling source at ",
      normalizePath(sibling$path, mustWork = FALSE), " is ", src, ".\n",
      "These tests are therefore exercising a stale backend, which is how ",
      "CODE_REVIEW.md R4.D2 arose: the suite reported 210 passing against a ",
      "backend six releases behind the code under review.\n\n",
      "Install the sibling source first:\n",
      "  devtools::install('../mouseExperiment', upgrade = FALSE)\n\n",
      "This is a hard failure rather than a skip on purpose — a green run against ",
      "the wrong backend is worse than no run."
    ))
  } else {
    succeed()
  }
})
