# =============================================================================
# Round 16 — the max-treedepth diagnostic cried wolf
#
# Unusual for this review: most findings are silent failures, this one is a false
# alarm. The cost is the same. A diagnostic that always reports a large number is
# one nobody can act on, and it masks the genuine saturation it exists to detect.
# =============================================================================

test_that("R16.1: the configured treedepth limit is read, with Stan's default", {
  # Stan's default max_treedepth is 10. The fallback matters because the count is
  # meaningless if the limit is guessed from the draws themselves.
  expect_equal(mouseExperiment:::.me_max_treedepth(NULL), 10)
  expect_equal(mouseExperiment:::.me_max_treedepth(list()), 10)
})

test_that("R16.1: hits are counted against the limit, not the observed maximum", {
  skip_if_not_installed("brms")

  set.seed(53)
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) {
    lam <- if (tx == "Control") 0.10 else 0.05
    for (m in 1:25) {
      k <- k + 1L
      t <- stats::rexp(1, lam)
      rows[[k]] <- data.frame(ID = paste0(tx, m), Treatment = tx,
                              Time = min(t, 60),
                              Event = as.integer(t <= 60),
                              stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)

  r <- suppressWarnings(suppressMessages(bayesian_survival(
    df, time_column = "Time", event_column = "Event",
    treatment_column = "Treatment", id_column = "ID",
    reference_group = "Control", family = "weibull",
    n_iter = me_test_niter(), n_warmup = me_test_niter(), n_chains = 2L,
    seed = 5L, plots = FALSE, verbose = FALSE)))

  np <- brms::nuts_params(r$model)
  td <- np$Value[np$Parameter == "treedepth__"]
  limit <- mouseExperiment:::.me_max_treedepth(r$model)

  # The reported count must equal the true count against the configured limit.
  expect_identical(r$nuts_diagnostics$n_max_treedepth,
                   as.integer(sum(td >= limit, na.rm = TRUE)))

  # And the specific regression: a well-behaved low-dimensional model whose
  # treedepths sit far below the limit must report 0, not "most draws". The old
  # code counted draws at max(td), which was 1425 of 1500 on this model.
  if (max(td, na.rm = TRUE) < limit) {
    expect_identical(r$nuts_diagnostics$n_max_treedepth, 0L)
    expect_gt(sum(td == max(td, na.rm = TRUE)), 0L)   # the old count was non-zero
  }
})

test_that("R16.1: clean reflects validity, not efficiency", {
  skip_if_not_installed("brms")
  # Documented behaviour: divergences and E-BFMI decide `clean`; treedepth is
  # reported but excluded, because hitting the limit makes sampling slow rather
  # than wrong. The roxygen used to claim all three were required, which the code
  # has never done.
  doc <- paste(readLines(testthat::test_path("..", "..", "R", "utils_bayes.R"),
                         warn = FALSE), collapse = "\n")
  expect_false(grepl("iff all three pass", doc, fixed = TRUE))
  expect_match(doc, "deliberately excluded", fixed = TRUE)
})
