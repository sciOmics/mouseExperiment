# =============================================================================
# README usage examples run (CODE_REVIEW.md R20.51)
#
# Round 20 found most README examples broken: wrong argument names, fields that
# do not exist, arm names that are not in the data, an undefined data frame.
# Every R block in the README's "Usage" section is extracted and run here, so a
# README that drifts from the API fails the suite. Blocks calling a Bayesian
# entry point are skipped: they compile Stan, and the Bayesian tests cover them.
#
# R CMD check runs the tests from a copy that has no README, so this file
# skips there; CI runs it from the source tree (.github/workflows).
# =============================================================================

readme_path <- function() {
  p <- testthat::test_path("..", "..", "README.md")
  skip_if_not(file.exists(p), "README.md not available (R CMD check copy)")
  p
}

readme_usage_blocks <- function() {
  lines <- readLines(readme_path(), warn = FALSE)
  start <- grep("^## Usage", lines)
  end   <- grep("^## ", lines); end <- end[end > start][1] - 1L
  sec   <- lines[start:end]
  open  <- grep("^```r\\s*$", sec); close <- grep("^```\\s*$", sec)
  lapply(open, function(o) sec[(o + 1L):(close[close > o][1] - 1L)])
}

test_that("R20.51: the README's Usage section has R examples to run", {
  blocks <- readme_usage_blocks()
  expect_gte(length(blocks), 8L)
})

test_that("R20.51: every non-Bayesian README example runs", {
  blocks <- readme_usage_blocks()
  ran <- 0L
  for (b in blocks) {
    code <- paste(b, collapse = "\n")
    if (grepl("bayesian_", code)) next
    env <- new.env(parent = globalenv())
    res <- tryCatch({
      suppressMessages(suppressWarnings(eval(parse(text = code), envir = env)))
      "ok"
    }, error = function(e) conditionMessage(e))
    expect_identical(res, "ok", info = paste0("README block failed:\n", code))
    ran <- ran + 1L
  }
  expect_gte(ran, 7L)
  grDevices::graphics.off()
})
