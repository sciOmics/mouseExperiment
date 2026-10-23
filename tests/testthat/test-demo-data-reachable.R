# =============================================================================
# Demo data reachability — CODE_REVIEW.md D.12
#
# Two demo datasets shipped by the backend (`master_synthetic_data`,
# `synthetic_data`) had no CSV here and no entry in the loader, so neither could
# ever be selected. Two further files — the treatment schedules — sat in this
# package's own `inst/sample_data/` with nothing in the app able to open them;
# they were removed under D.15 along with plot_treatments().
#
# Nothing detected any of it, because a demo dataset that is never offered
# produces no error — it simply is not there. That is the §K.2 shape again: the
# absence is invisible.
#
# These tests close the loop in both directions: every file that ships must be
# reachable, and every loader entry must resolve to a file that actually parses.
#
# 2026-10-09: the "Simple" and "Necrotic" demos were removed, the others were
# given short names, and "Master" became "All Features" (listed last), with a
# removal-reason column so every module, including the reason-dependent
# options, has data to show.
# =============================================================================

sample_dir <- function() {
  p <- system.file("sample_data", package = "mouseExperiment")
  if (!nzchar(p) || !dir.exists(p)) p <- testthat::test_path("..", "..", "inst",
                                                             "sample_data")
  skip_if_not(dir.exists(p), "sample_data directory not locatable")
  p
}

# The loader's table, read out of the source so the test cannot drift from it.
loader_keys <- function() {
  src <- dash_source("R", "mod_data_upload.R")   # skips outside the source tree
  blk <- regmatches(src, regexpr("demo_files <- list\\([\\s\\S]*?\\n    \\)",
                                 src, perl = TRUE))
  skip_if(length(blk) != 1L, "demo_files block not locatable")
  files <- regmatches(blk, gregexpr('data\\s*=\\s*"[^"]+"', blk))[[1]]
  sub('.*"([^"]+)".*', "\\1", files)
}

# D.15/D.16: the two schedule CSVs that used to sit here unreachable were removed
# along with plot_treatments(). The dosing-annotation feature that briefly replaced
# it was withdrawn too (D.16), but the files stay gone either way — nothing ships
# that cannot be loaded, so the exemption list is empty and should stay so.
SCHEDULE_FILES <- character(0)

test_that("D.12: every shipped sample file is loadable from the demo selector", {
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

test_that("D.12: the All Features demo exercises every optional column", {
  # Its whole reason for being: one dataset that reaches every module, so a new
  # user can see the app work without supplying their own data. If a column is
  # ever dropped from it, the demo quietly stops covering a module.
  p <- file.path(sample_dir(), "all_features_synthetic_data.csv")
  skip_if_not(file.exists(p), "All Features demo not present")
  d <- validate_uploaded_file(p, "csv")$data
  dc <- detect_columns(colnames(d))

  for (role in c("id", "day", "volume", "treatment", "cage", "weight", "dose",
                 "necrotic", "survival_event", "length", "width", "removal_reason")) {
    expect_false(is.null(dc[[role]]), info = paste("All Features demo lost:", role))
  }
  # Height is what makes the 3-axis ellipsoid formula demonstrable; it has no
  # detect_columns role, so assert the column directly.
  expect_true("Height" %in% colnames(d))
  expect_gt(length(unique(d$Treatment)), 4L)

  # Each reason-dependent option needs at least one animal to act on: deaths
  # (Survival's "count deaths as events"), weight-loss removals (Toxicity), and
  # the limit and planned-end removals they are told apart from.
  reasons <- removal_reason_values(d, dc$removal_reason)
  expect_setequal(unique(guess_removal_category(reasons)),
                  c("limit", "death", "weight", "planned"))
})

test_that("Demo selector: four short names, All Features last", {
  # 2026-10-09: "Simple" and "Necrotic" were removed, the parenthetical
  # descriptions dropped, and "Master" renamed "All Features" and moved last.
  for (f in c("synthetic_data.csv", "necrotic_synthetic_data.csv", "master_synthetic_data.csv"))
    expect_false(file.exists(file.path(sample_dir(), f)), info = paste("still shipped:", f))
  src <- dash_source("R", "ui_data.R")
  blk <- regmatches(src, regexpr('"-- Select a demo dataset --" = ""[\\s\\S]*?\\n                \\)',
                                 src, perl = TRUE))
  skip_if(length(blk) != 1L, "demo_dataset choices not locatable")
  labels <- sub('^"([^"]+)".*', "\\1",
                regmatches(blk, gregexpr('"[^"]+"\\s*=\\s*"[a-z_]+"', blk))[[1]])
  expect_identical(labels, c("Combo Treatment", "Dose Levels", "Weight", "All Features"))
})

test_that("D.12: the UI selector offers exactly what the loader knows", {
  # A choice in the dropdown with no matching loader entry silently does nothing
  # when clicked; a loader entry with no choice is unreachable. Both directions.
  src <- dash_source("R", "ui_data.R")
  blk <- regmatches(src, regexpr('"-- Select a demo dataset --" = ""[\\s\\S]*?\\n                \\)',
                                 src, perl = TRUE))
  skip_if(length(blk) != 1L, "demo_dataset choices not locatable")
  ui_keys <- sub('.*=\\s*"([^"]+)"$', "\\1",
                 regmatches(blk, gregexpr('=\\s*"[a-z_]+"', blk))[[1]])
  ui_keys <- setdiff(ui_keys, "")

  src2 <- dash_source("R", "mod_data_upload.R")
  blk2 <- regmatches(src2, regexpr("demo_files <- list\\([\\s\\S]*?\\n    \\)",
                                   src2, perl = TRUE))
  loader_names <- gsub("[^a-z_]", "",
                       regmatches(blk2, gregexpr("\\n      [a-z_]+ = list", blk2))[[1]])
  loader_names <- sub("list$", "", loader_names)

  expect_setequal(ui_keys, loader_names)
})
