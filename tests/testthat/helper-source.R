# Tests that read the package's own R sources can only run from the source
# tree. R CMD check runs the tests against the installed package, where R/ does
# not exist: there a file list came back empty and such tests passed without
# checking anything, or failed to open a file. They now skip instead, and CI
# runs them from the source tree (CODE_REVIEW.md R20.53, R20.55).
pkg_source_files <- function(...) {
  dir <- testthat::test_path("..", "..", "R")
  testthat::skip_if_not(dir.exists(dir), "package R/ sources not available (R CMD check)")
  if (length(list(...))) file.path(dir, ...)
  else list.files(dir, pattern = "[.]R$", full.names = TRUE)
}
