# =============================================================================
# Demo data reachability — CODE_REVIEW.md D.12
#
# Two demo datasets shipped by the backend (`master_synthetic_data`,
# `synthetic_data`) had no CSV here and no entry in the loader, so neither could
# ever be selected. Two more files (`combo_treatment_schedule.csv`,
# `dose_levels_treatment_schedule.csv`) sit in this package's own
# `inst/sample_data/` with nothing in the app able to open them.
#
# Nothing detected any of it, because a demo dataset that is never offered
# produces no error — it simply is not there. That is the §K.2 shape again: the
# absence is invisible.
#
# These tests close the loop in both directions: every file that ships must be
# reachable or explicitly classified, and every loader entry must resolve to a
# file that actually parses.
# =============================================================================

sample_dir <- function() {
  p <- system.file("sample_data", package = "mouseExperimentDashboard")
  if (!nzchar(p) || !dir.exists(p)) p <- testthat::test_path("..", "..", "inst",
                                                             "sample_data")
  skip_if_not(dir.exists(p), "sample_data directory not locatable")
  p
}

# The loader's table, read out of the source so the test cannot drift from it.
loader_keys <- function() {
  src <- paste(readLines(testthat::test_path("..", "..", "R", "mod_data_upload.R"),
                         warn = FALSE), collapse = "\n")
  blk <- regmatches(src, regexpr("demo_files <- list\\([\\s\\S]*?\\n    \\)",
                                 src, perl = TRUE))
  skip_if(length(blk) != 1L, "demo_files block not locatable")
  files <- regmatches(blk, gregexpr('data\\s*=\\s*"[^"]+"', blk))[[1]]
  sub('.*"([^"]+)".*', "\\1", files)
}

# Schedules are dosing tables (Treatment / Day / Dose), not measurement data.
# Loading one as the main dataset would break every analysis, so they are
# deliberately not in the demo selector. Listed here so the completeness test
# below distinguishes "intentionally not offered" from "forgotten".
SCHEDULE_FILES <- c("combo_treatment_schedule.csv",
                    "dose_levels_treatment_schedule.csv")

test_that("D.12: every shipped sample file is either loadable or a known schedule", {
  shipped <- list.files(sample_dir(), pattern = "\\.csv$")
  offered <- loader_keys()
  unreachable <- setdiff(shipped, c(offered, SCHEDULE_FILES))
  expect_length(unreachable, 0L)
  if (length(unreachable)) {
    cat("\nShipped but unreachable:\n  ", paste(unreachable, collapse = "\n  "), "\n")
  }
})

test_that("D.12: every loader entry resolves to a file that parses", {
  # The other direction: an entry naming a file that is not shipped fails at
  # runtime with "Demo file not found", which users see and developers do not.
  dir <- sample_dir()
  for (f in loader_keys()) {
    p <- file.path(dir, f)
    expect_true(file.exists(p), info = paste("loader names a missing file:", f))
    if (file.exists(p)) {
      res <- validate_uploaded_file(p, "csv")
      expect_true(res$valid, info = paste(f, "-", res$message %||% ""))
    }
  }
})

test_that("D.12: the master demo exercises every optional column", {
  # Its whole reason for being added: one dataset that reaches every module, so
  # a new user can see the app work without supplying their own data. If a
  # column is ever dropped from it, the demo quietly stops covering a module.
  p <- file.path(sample_dir(), "master_synthetic_data.csv")
  skip_if_not(file.exists(p), "master demo not present")
  d <- validate_uploaded_file(p, "csv")$data
  dc <- detect_columns(colnames(d))

  for (role in c("id", "day", "volume", "treatment", "cage", "weight",
                 "dose", "necrotic", "survival_event", "length", "width")) {
    expect_false(is.null(dc[[role]]), info = paste("master demo lost:", role))
  }
  # Height is what makes the 3-axis ellipsoid formula demonstrable; it has no
  # detect_columns role, so assert the column directly.
  expect_true("Height" %in% colnames(d))
  expect_gt(length(unique(d$Treatment)), 4L)
})

test_that("D.12: the simple demo is genuinely minimal but sufficient", {
  p <- file.path(sample_dir(), "synthetic_data.csv")
  skip_if_not(file.exists(p), "simple demo not present")
  d <- validate_uploaded_file(p, "csv")$data
  dc <- detect_columns(colnames(d))
  for (role in c("id", "day", "volume", "treatment")) {
    expect_false(is.null(dc[[role]]), info = paste("simple demo lost:", role))
  }
  expect_gte(length(unique(d$Treatment)), 2L)
})

test_that("D.12: the UI selector offers exactly what the loader knows", {
  # A choice in the dropdown with no matching loader entry silently does nothing
  # when clicked; a loader entry with no choice is unreachable. Both directions.
  src <- paste(readLines(testthat::test_path("..", "..", "R", "ui_modules.R"),
                         warn = FALSE), collapse = "\n")
  blk <- regmatches(src, regexpr('"-- Select a demo dataset --" = ""[\\s\\S]*?\\n                \\)',
                                 src, perl = TRUE))
  skip_if(length(blk) != 1L, "demo_dataset choices not locatable")
  ui_keys <- sub('.*=\\s*"([^"]+)"$', "\\1",
                 regmatches(blk, gregexpr('=\\s*"[a-z_]+"', blk))[[1]])
  ui_keys <- setdiff(ui_keys, "")

  src2 <- paste(readLines(testthat::test_path("..", "..", "R", "mod_data_upload.R"),
                          warn = FALSE), collapse = "\n")
  blk2 <- regmatches(src2, regexpr("demo_files <- list\\([\\s\\S]*?\\n    \\)",
                                   src2, perl = TRUE))
  loader_names <- gsub("[^a-z_]", "",
                       regmatches(blk2, gregexpr("\\n      [a-z_]+ = list", blk2))[[1]])
  loader_names <- sub("list$", "", loader_names)

  expect_setequal(ui_keys, loader_names)
})
