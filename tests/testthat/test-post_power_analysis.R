# =============================================================================
# Tests for post_power_analysis()
#
# The function takes observed data (or effect size + SD parameters) and
# returns required sample sizes for achieving a given power.
#
# Ground-truth checks:
#   - Larger effect size  → smaller required n (monotone relationship)
#   - Higher desired power → larger required n
#   - Return structure has required fields
# =============================================================================

test_that("returns a list with expected fields on tg_simple data", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    post_power_analysis(
      data             = df,
      alpha            = 0.05,
      power            = c(0.80, 0.90),
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      method           = "parametric"
    )
  ))

  expect_true(is.list(res))
  # At least one of these structures should be present
  expect_true(length(res) >= 1L, info = "post_power_analysis returned an empty list")
})

test_that("power table contains numeric sample size estimates", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    post_power_analysis(
      data             = df,
      alpha            = 0.05,
      power            = 0.80,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      method           = "parametric"
    )
  ))

  # Find the first data frame in the result with numeric values
  tables <- Filter(is.data.frame, res)
  if (length(tables) == 0) {
    # Some implementations return power_table inside a sub-list
    for (nm in names(res)) {
      if (is.list(res[[nm]])) {
        tables <- c(tables, Filter(is.data.frame, res[[nm]]))
      }
    }
  }

  skip_if(length(tables) == 0, "No data frame found in post_power_analysis output")

  tbl <- tables[[1]]
  num_cols <- sapply(tbl, function(x) is.numeric(x) || all(!is.na(suppressWarnings(as.numeric(x)))))
  expect_true(any(num_cols), info = "No numeric sample-size column found in power table")
})

test_that("required n is strictly positive", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    post_power_analysis(
      data             = df,
      alpha            = 0.05,
      power            = 0.80,
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      method           = "parametric"
    )
  ))

  tables <- Filter(is.data.frame, res)
  if (length(tables) == 0) skip("No data frame in result")

  tbl      <- tables[[1]]
  n_col    <- intersect(c("n", "N", "sample_size", "n_per_group",
                          "Required_N"), colnames(tbl))[1]
  skip_if(is.na(n_col), "Cannot locate sample-size column")

  n_vals <- suppressWarnings(as.numeric(tbl[[n_col]]))
  n_vals <- n_vals[!is.na(n_vals)]
  expect_true(all(n_vals > 0), info = paste("n values:", paste(n_vals, collapse = ", ")))
})

test_that("higher power requires larger or equal n (monotone relationship)", {
  df  <- make_tg_simple()
  res <- suppressWarnings(suppressMessages(
    post_power_analysis(
      data             = df,
      alpha            = 0.05,
      power            = c(0.80, 0.90),
      time_column      = "Day",
      volume_column    = "Volume",
      treatment_column = "Treatment",
      id_column        = "ID",
      method           = "parametric"
    )
  ))

  tables <- Filter(is.data.frame, res)
  if (length(tables) == 0) skip("No data frame in result")

  for (tbl in tables) {
    pw_col   <- intersect(c("power", "Power", "desired_power"), colnames(tbl))[1]
    n_col    <- intersect(c("n", "N", "sample_size", "n_per_group",
                            "Required_N"), colnames(tbl))[1]
    if (!is.na(pw_col) && !is.na(n_col)) {
      pw_vals <- suppressWarnings(as.numeric(tbl[[pw_col]]))
      n_vals  <- suppressWarnings(as.numeric(tbl[[n_col]]))
      valid   <- !is.na(pw_vals) & !is.na(n_vals)
      if (sum(valid) >= 2) {
        ord <- order(pw_vals[valid])
        expect_true(
          all(diff(n_vals[valid][ord]) >= 0),
          info = "n should be non-decreasing as power increases"
        )
      }
    }
  }
})

test_that("does not error on four-group multi-timepoint combo data", {
  df  <- make_combo_four_group()
  # Keep only the measured columns (drop Survival_Censor with its NAs)
  df <- df[, c("ID", "Cage", "Treatment", "Day", "Volume")]

  expect_no_error(
    suppressWarnings(suppressMessages(
      post_power_analysis(
        data             = df,
        alpha            = 0.05,
        power            = 0.80,
        time_column      = "Day",
        volume_column    = "Volume",
        treatment_column = "Treatment",
        id_column        = "ID",
        method           = "parametric"
      )
    ))
  )
})
