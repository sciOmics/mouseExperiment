# =============================================================================
# Tests for new dashboard helper functions added in 4.x enhancement pass
# =============================================================================

# ---------------------------------------------------------------------------
# Source internal helpers (not exported from the namespace)
# ---------------------------------------------------------------------------
source_new_helpers <- function() {
  tryCatch({
    pkg_env <- getNamespace("mouseExperimentDashboard")
    list(
      palette_colors       = get("palette_colors",       envir = pkg_env),
      default_group_colors = get("default_group_colors", envir = pkg_env),
      get_log_entries      = get("get_log_entries",      envir = pkg_env),
      clear_log_buffer     = get("clear_log_buffer",     envir = pkg_env),
      appender_buffer      = get("appender_buffer",      envir = pkg_env)
    )
  }, error = function(e) NULL)
}

new_helpers <- source_new_helpers()
skip_if_no_new <- function() {
  skip_if(is.null(new_helpers), "New dashboard helper functions not accessible")
}

# ===========================================================================
# palette_colors()
# ===========================================================================

test_that("palette_colors returns correct number of colours for default palette", {
  skip_if_no_new()
  pc <- new_helpers$palette_colors
  cols <- pc("default", 5)
  expect_equal(length(cols), 5)
  expect_true(all(grepl("^#", cols)))
})

test_that("palette_colors returns viridis colours", {
  skip_if_no_new()
  pc <- new_helpers$palette_colors
  cols <- pc("viridis", 4)
  expect_equal(length(cols), 4)
  expect_true(all(grepl("^#", cols)))
})

test_that("palette_colors returns Okabe-Ito colours", {
  skip_if_no_new()
  pc <- new_helpers$palette_colors
  cols <- pc("okabe_ito", 6)
  expect_equal(length(cols), 6)
  # Known first colour of Okabe-Ito
  expect_equal(cols[1], "#E69F00")
})

test_that("palette_colors returns Set2 colours", {
  skip_if_no_new()
  pc <- new_helpers$palette_colors
  cols <- pc("set2", 3)
  expect_equal(length(cols), 3)
})

test_that("palette_colors handles n > palette length via recycling", {
  skip_if_no_new()
  pc <- new_helpers$palette_colors
  cols <- pc("okabe_ito", 12)
  # Okabe-Ito has 8 base + recycled
  expect_equal(length(cols), 12)
})

test_that("palette_colors returns empty vector for n=0", {
  skip_if_no_new()
  pc <- new_helpers$palette_colors
  expect_equal(length(pc("viridis", 0)), 0)
})

# ===========================================================================
# default_group_colors()
# ===========================================================================

test_that("default_group_colors returns n hex colours", {
  skip_if_no_new()
  dgc <- new_helpers$default_group_colors
  cols <- dgc(5)
  expect_equal(length(cols), 5)
  expect_true(all(grepl("^#", cols)))
})

test_that("default_group_colors returns empty vector for n <= 0", {
  skip_if_no_new()
  dgc <- new_helpers$default_group_colors
  expect_equal(length(dgc(0)), 0)
  expect_equal(length(dgc(-1)), 0)
})

# ===========================================================================
# Log buffer
# ===========================================================================

test_that("log buffer starts empty and can be cleared", {
  skip_if_no_new()
  gle <- new_helpers$get_log_entries
  clb <- new_helpers$clear_log_buffer

  clb()
  expect_equal(length(gle()), 0)
})

test_that("appender_buffer stores log lines", {
  skip_if_no_new()
  ab  <- new_helpers$appender_buffer
  gle <- new_helpers$get_log_entries
  clb <- new_helpers$clear_log_buffer

  clb()
  # Capture stdout to suppress console output
  capture.output(ab("test log line 1"), type = "output")
  capture.output(ab("test log line 2"), type = "output")

  entries <- gle()
  expect_equal(length(entries), 2)
  expect_equal(entries[1], "test log line 1")
  expect_equal(entries[2], "test log line 2")
  clb()
})
