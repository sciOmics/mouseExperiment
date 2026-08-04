# =============================================================================
# Cross-repo call resolution — CODE_REVIEW.md R17.3
#
# Three `mouseExperiment::` calls pointed at functions that do not exist:
#
#   plot_survival                 removed from the backend in ace5b94
#   prepare_dose_data             internal helper, never exported
#   generate_summary_statistics   internal helper, never exported
#
# All three sat inside tryCatch, so every survival and dose-response run threw,
# logged a caught error and silently took a fallback path. Users saw output, so
# nothing looked wrong. The KM plot case is the sharpest: a comment said "use the
# exported plot_survival() function (matches demo notebook approach)", so a
# specific intended plot had never once been drawn and everyone was looking at
# the fallback.
#
# This is the same cross-repo staleness as R4.D2, where the suite passed against
# a backend six releases old. The check is four lines and could have caught all
# three the day they broke, so it is now permanent.
# =============================================================================

test_that("R17.3: every mouseExperiment:: call resolves to an exported function", {
  skip_if_not_installed("mouseExperiment")

  r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "R/ not locatable")

  # Strip comments before scanning. Comments naming a removed function are fine
  # and in fact wanted -- they explain what replaced it. What must not survive is
  # a live call, which errors at runtime. (Getting this wrong is how the earlier
  # plot_treatments guard produced false positives on its own documentation.)
  code <- unlist(lapply(list.files(r_dir, pattern = "[.]R$", full.names = TRUE),
                        function(f) sub("#.*$", "", readLines(f, warn = FALSE))))

  hits <- unlist(regmatches(
    code, gregexpr("mouseExperiment::[A-Za-z0-9_.]+", code)))
  called <- unique(sub("mouseExperiment::", "", hits))
  expect_gt(length(called), 0L)   # the scan itself must not silently find nothing

  exported <- getNamespaceExports("mouseExperiment")
  unresolved <- setdiff(called, exported)

  if (length(unresolved)) {
    cat("\nCalls to non-existent backend functions:\n  ",
        paste(unresolved, collapse = "\n  "), "\n")
  }
  expect_length(unresolved, 0L)
})

test_that("R17.3: the scan would actually catch a broken call", {
  # A guard that cannot fail is not a guard. Confirm the matcher fires on a
  # fabricated call and that comment-stripping does not defeat it.
  fake <- c('  x <- mouseExperiment::this_does_not_exist(1)',
            '  # mouseExperiment::also_fake() -- only a comment')
  code <- sub("#.*$", "", fake)
  called <- unique(sub("mouseExperiment::", "",
                       unlist(regmatches(code,
                         gregexpr("mouseExperiment::[A-Za-z0-9_.]+", code)))))
  expect_true("this_does_not_exist" %in% called)
  expect_false("also_fake" %in% called)   # comments must not trip it
})

test_that("R17.3: the three previously-dead calls are gone from live code", {
  r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "R/ not locatable")
  code <- unlist(lapply(list.files(r_dir, pattern = "[.]R$", full.names = TRUE),
                        function(f) sub("#.*$", "", readLines(f, warn = FALSE))))
  for (fn in c("plot_survival", "prepare_dose_data",
               "generate_summary_statistics")) {
    expect_false(any(grepl(paste0("mouseExperiment::", fn, "\\b"), code)),
                 info = paste("live call to removed/internal backend fn:", fn))
  }
})

test_that("R17.3: dose-response plots have the data they need", {
  # The consequence of the prepare_dose_data call failing: result$analysis_data
  # was never set, three renderers hit a silent req(ad), and every dose-response
  # plot came out blank. The backend now returns the frame it analysed.
  skip_if_not_installed("mouseExperiment")
  set.seed(3)
  doses <- c(0, 5, 10, 25, 50)
  rows <- list(); k <- 0L
  for (dz in doses) for (m in 1:6) {
    k <- k + 1L
    mu <- 100 + 900 / (1 + (max(dz, 1e-9) / 25)^1.5)
    rows[[k]] <- data.frame(ID = paste0("d", dz, "_", m),
                            Treatment = paste0("D", dz), Dose = dz, Day = 21,
                            Volume = max(1, stats::rnorm(1, mu, 60)),
                            stringsAsFactors = FALSE)
  }
  res <- suppressWarnings(suppressMessages(
    mouseExperiment::dose_response_statistics(
      do.call(rbind, rows), dose_column = "Dose", treatment_column = "Treatment",
      volume_column = "Volume", day_column = "Day", id_column = "ID",
      control_group_name = "D0", verbose = FALSE)))

  expect_false(is.null(res$analysis_data))
  expect_s3_class(as.data.frame(res$analysis_data), "data.frame")
  expect_gt(nrow(res$analysis_data), 0L)
  expect_true(all(c("Dose", "Volume") %in% names(res$analysis_data)))
})
