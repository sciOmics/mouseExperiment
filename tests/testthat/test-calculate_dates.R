# =============================================================================
# Tests for calculate_dates()
# Verifies date-to-study-day arithmetic for all supported date formats.
# =============================================================================

# Construct a small data frame with known dates relative to a start date.
# start_date = "2025-03-24" (ISO), measurement dates: same day, +7, +14, +21
make_date_df <- function(date_str, fmt) {
  data.frame(Date = date_str, stringsAsFactors = FALSE)
}

test_that("ISO YYYY-MM-DD: day 0 when date equals start_date", {
  df  <- data.frame(Date = "2025-03-24", stringsAsFactors = FALSE)
  out <- suppressMessages(calculate_dates(df, start_date = "2025-03-24"))

  expect_true("Day" %in% colnames(out))
  expect_equal(out$Day, 0L)
})

test_that("ISO YYYY-MM-DD: correct day count for multiple rows", {
  df <- data.frame(
    Date = c("2025-03-24", "2025-03-31", "2025-04-07", "2025-04-14"),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(calculate_dates(df, start_date = "2025-03-24"))

  expect_equal(out$Day, c(0, 7, 14, 21))
})

test_that("MM/DD/YYYY format: correct day count", {
  df <- data.frame(
    Date = c("03/24/2025", "03/31/2025", "04/07/2025"),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(
    calculate_dates(df, start_date = "03/24/2025", date_format = "%m/%d/%Y")
  )

  expect_equal(out$Day, c(0, 7, 14))
})

test_that("DD-Mon format with explicit year: correct day count", {
  df <- data.frame(
    Date = c("24-Mar", "31-Mar", "07-Apr"),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(
    calculate_dates(df, start_date = "24-Mar",
                    date_format = "%d-%b", year = 2025)
  )

  expect_equal(out$Day, c(0, 7, 14))
})

test_that("Day column is non-negative for dates >= start_date", {
  df <- data.frame(
    Date = c("2025-03-24", "2025-04-14"),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(calculate_dates(df, start_date = "2025-03-24"))

  expect_true(all(out$Day >= 0))
})

test_that("original columns are preserved after date conversion", {
  df <- data.frame(
    ID   = c("M1", "M2"),
    Date = c("2025-03-24", "2025-03-31"),
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(calculate_dates(df, start_date = "2025-03-24"))

  expect_true(all(c("ID", "Date", "Day") %in% colnames(out)))
})

test_that("missing date column raises an error", {
  df <- data.frame(Day = 0)
  expect_error(calculate_dates(df, start_date = "2025-03-24"), "Date")
})

test_that("demo CSV dates round-trip correctly (combo dataset)", {
  # The demo combo CSV has dates in the form "03/24/2025"
  demo_path <- system.file(
    "sample_data", "combo_treatment_synthetic_data.csv",
    package = "mouseExperiment"
  )
  skip_if_not(nzchar(demo_path) && file.exists(demo_path),
              "Demo data not installed")

  raw <- utils::read.csv(demo_path, stringsAsFactors = FALSE)
  skip_if_not("Date" %in% colnames(raw), "Date column not present")

  out <- suppressMessages(
    calculate_dates(raw, start_date = "03/24/2025", date_format = "%m/%d/%Y")
  )

  expect_true("Day" %in% colnames(out))
  expect_true(all(out$Day >= 0))
  expect_equal(min(out$Day), 0)
})
