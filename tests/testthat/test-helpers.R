# =============================================================================
# Dashboard helper tests
# =============================================================================
#
# Source the helpers directly (they're not exported).
# In devtools::test() / R CMD check these run from the package root.
# =============================================================================

# ---------------------------------------------------------------------------
# Source internal helpers (not exported from the namespace)
# ---------------------------------------------------------------------------
source_helpers <- function() {
  tryCatch({
    pkg_env <- getNamespace("mouseExperiment")
    list(
      validate_uploaded_file  = get("validate_uploaded_file",  envir = pkg_env),
      detect_columns          = get("detect_columns",          envir = pkg_env),
      format_df_for_display   = get("format_df_for_display",   envir = pkg_env),
      format_pval             = get("format_pval",             envir = pkg_env)
    )
  }, error = function(e) NULL)
}

helpers <- source_helpers()
skip_if_no_helpers <- function() {
  skip_if(is.null(helpers), "Dashboard helper functions not accessible")
}

# Helper: write a small CSV to a temp file and return its path
write_temp_csv <- function(df) {
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(df, tmp, row.names = FALSE)
  tmp
}

# ===========================================================================
# validate_uploaded_file()
# ===========================================================================

test_that("validate_uploaded_file: valid CSV returns valid=TRUE and a data frame", {
  skip_if_no_helpers()
  vuf <- helpers$validate_uploaded_file

  df  <- data.frame(ID = 1:3, Value = c(1.1, 2.2, 3.3))
  tmp <- write_temp_csv(df)
  on.exit(unlink(tmp))

  res <- vuf(tmp, "csv")
  expect_true(res$valid)
  expect_true(is.data.frame(res$data))
  expect_equal(nrow(res$data), 3L)
})

test_that("validate_uploaded_file: non-existent file returns valid=FALSE", {
  skip_if_no_helpers()
  vuf <- helpers$validate_uploaded_file

  res <- vuf("/tmp/this_file_does_not_exist_abc123.csv", "csv")
  expect_false(res$valid)
  expect_true(nzchar(res$message))
})

test_that("validate_uploaded_file: empty CSV (header only) returns valid=FALSE", {
  skip_if_no_helpers()
  vuf <- helpers$validate_uploaded_file

  tmp <- tempfile(fileext = ".csv")
  writeLines("ID,Value", tmp)
  on.exit(unlink(tmp))

  res <- vuf(tmp, "csv")
  expect_false(res$valid)
})

# ===========================================================================
# detect_columns()
# ===========================================================================

test_that("detect_columns: detects standard column names", {
  skip_if_no_helpers()
  dc <- helpers$detect_columns

  cols <- c("ID", "Cage", "Treatment", "Day", "Volume", "Survival_Censor")
  res  <- dc(cols)

  expect_equal(res$id,       "ID",              ignore_attr = TRUE)
  expect_equal(res$cage,     "Cage",            ignore_attr = TRUE)
  expect_equal(res$treatment,"Treatment",        ignore_attr = TRUE)
})

test_that("detect_columns: detects Date column (not Day) when Day absent", {
  skip_if_no_helpers()
  dc <- helpers$detect_columns

  cols <- c("ID", "Cage", "Treatment", "Date", "Length", "Width")
  res  <- dc(cols)

  expect_true(!is.null(res$date) || !is.null(res$day) || nzchar(res$date %||% ""),
              info = "Expected date/day detection when Date column present")
})

test_that("detect_columns: detects Dose column", {
  skip_if_no_helpers()
  dc <- helpers$detect_columns

  cols <- c("ID", "Cage", "Treatment", "Dose", "Day", "Volume")
  res  <- dc(cols)

  dose_val <- res$dose %||% res[["Dose"]] %||% ""
  expect_true(nzchar(dose_val) || "dose" %in% tolower(names(res)),
              info = "Expected Dose column to be detected")
})

test_that("detect_columns: returns a named list", {
  skip_if_no_helpers()
  dc <- helpers$detect_columns

  cols <- c("ID", "Treatment", "Day", "Volume")
  res  <- dc(cols)

  expect_true(is.list(res))
  expect_true(length(names(res)) > 0)
})

# ===========================================================================
# format_pval()
# ===========================================================================

test_that("format_pval: formats numeric p-value as string", {
  skip_if_no_helpers()
  fp <- helpers$format_pval

  out <- fp(0.0123)
  expect_true(is.character(out))
  expect_true(nzchar(out))
})

test_that("format_pval: small p-value uses scientific or < notation", {
  skip_if_no_helpers()
  fp <- helpers$format_pval

  out <- fp(1e-9)
  expect_true(is.character(out))
  # Should not just be "0"
  expect_false(out == "0")
})

test_that("format_pval: NA input returns NA_character_ or empty string, not error", {
  skip_if_no_helpers()
  fp <- helpers$format_pval

  expect_no_error({
    out <- fp(NA_real_)
    expect_true(is.na(out) || is.character(out))
  })
})

test_that("format_pval: p ≈ 1 returns a non-empty string", {
  skip_if_no_helpers()
  fp <- helpers$format_pval

  out <- fp(1.0)
  expect_true(is.character(out))
  expect_true(nzchar(out))
})

# ===========================================================================
# format_df_for_display()
# ===========================================================================

test_that("format_df_for_display: returns a data frame", {
  skip_if_no_helpers()
  fdf <- helpers$format_df_for_display

  df  <- data.frame(A = c(1.23456, 2.34567), B = c("x", "y"),
                    stringsAsFactors = FALSE)
  out <- fdf(df)
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), nrow(df))
})

test_that("format_df_for_display: numeric columns are rounded to specified digits", {
  skip_if_no_helpers()
  fdf <- helpers$format_df_for_display

  df  <- data.frame(Value = c(1.23456789))
  out <- fdf(df, digits = 3)

  # After formatting, the column should be character or a rounded numeric
  val <- as.numeric(out$Value)
  expect_equal(val, 1.235, tolerance = 0.0001)
})

test_that("format_df_for_display: handles table/array class columns without error", {
  skip_if_no_helpers()
  fdf <- helpers$format_df_for_display

  # Simulate what survival_statistics returns: table/array-class columns
  df <- data.frame(Group = c("A", "B"), stringsAsFactors = FALSE)
  df$Events <- structure(c(3L, 5L), class = "table", dim = 2L,
                         dimnames = list(c("A", "B")))
  df$Total  <- structure(c(6L, 6L), class = "table", dim = 2L,
                         dimnames = list(c("A", "B")))

  expect_no_error({
    out <- fdf(df)
    expect_true(is.data.frame(out))
  })
})

test_that("format_df_for_display: NA values in numeric column do not cause error", {
  skip_if_no_helpers()
  fdf <- helpers$format_df_for_display

  df  <- data.frame(Value = c(1.5, NA, 3.14))
  expect_no_error(fdf(df))
})

# Pipe-friendly NULL coalescing for tests above
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
