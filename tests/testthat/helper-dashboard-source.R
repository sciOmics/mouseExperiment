# =============================================================================
# Tests that read the dashboard's own sources (R/, DESCRIPTION, .Rbuildignore)
# can only run from the source tree. The VPS staging gate runs the installed
# tests (--install-tests), and R CMD check runs a copy; neither has R/ beside
# tests/. There such reads skip instead of erroring (CODE_REVIEW.md R20-D.T).
# Within a test, read sources after the behavioural checks, so those still run.
# =============================================================================

dash_source_dir <- function() {
  dir <- testthat::test_path("..", "..")
  testthat::skip_if_not(dir.exists(file.path(dir, "R")),
                        "dashboard sources not available (installed tests)")
  dir
}

# The text of a source file, e.g. dash_source("R", "app.R").
dash_source <- function(...) {
  paste(readLines(file.path(dash_source_dir(), ...), warn = FALSE), collapse = "\n")
}

# Every file in R/, concatenated.
dash_all_sources <- function() {
  files <- list.files(file.path(dash_source_dir(), "R"), pattern = "[.]R$")
  paste(vapply(files, function(f) dash_source("R", f), ""), collapse = "\n")
}
