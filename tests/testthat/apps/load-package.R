# Sourced by the shinytest2 app directories next to this file. AppDriver runs
# each app in a fresh background R process, so the app has to load the package
# itself (CODE_REVIEW.md R20.D11). From a source checkout that is the source
# tree (devtools::test()); under R CMD check, where the tests are copied out of
# the source tree, it is the installed package being checked.
local({
  root <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
  from_source <- file.exists(file.path(root, "DESCRIPTION")) &&
    any(grepl("^Package: mouseExperiment$",
              readLines(file.path(root, "DESCRIPTION"), warn = FALSE))) &&
    requireNamespace("pkgload", quietly = TRUE)
  if (from_source) {
    pkgload::load_all(root, quiet = TRUE, helpers = FALSE, export_all = TRUE)
  } else {
    library(mouseExperiment)
  }
})
