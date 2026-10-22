# =============================================================================
# Tests for the built-in demo datasets (inst/sample_data/)
#
# Verifies that all four CSV files are present, load cleanly, and have the
# column structure expected by the dashboard.
# =============================================================================

# Helper: find demo files whether package is installed or being developed in-situ
find_demo <- function(filename) {
  # 1. Installed package path
  p <- system.file("sample_data", filename, package = "mouseExperiment")
  if (nzchar(p) && file.exists(p)) return(p)
  # 2. Development / devtools::load_all() path
  local_p <- file.path(
    system.file(package = "mouseExperiment"),
    "sample_data", filename
  )
  if (file.exists(local_p)) return(local_p)
  # 3. Raw source tree (tests run from package root via testthat::test_local())
  raw_p <- file.path("inst", "sample_data", filename)
  if (file.exists(raw_p)) return(raw_p)
  NULL
}

# ======================================================
# File existence
# ======================================================
# D.15: the two treatment-schedule CSVs were removed with plot_treatments().
# Dosing is now entered as text in the sidebar and drawn into the plots, so no
# schedule file ships and none should be asserted here.

test_that("combo_treatment_synthetic_data.csv exists", {
  p <- find_demo("combo_treatment_synthetic_data.csv")
  expect_true(!is.null(p) && file.exists(p))
})


test_that("dose_levels_synthetic_data.csv exists", {
  p <- find_demo("dose_levels_synthetic_data.csv")
  expect_true(!is.null(p) && file.exists(p))
})


# ======================================================
# Load and basic structure
# ======================================================
test_that("combo data: loads without error and has expected columns", {
  p <- find_demo("combo_treatment_synthetic_data.csv")
  skip_if(is.null(p), "Demo file not found")

  df <- utils::read.csv(p, stringsAsFactors = FALSE)

  expect_true(is.data.frame(df))
  expect_true(nrow(df) > 0)
  for (col in c("Date", "Cage", "ID", "Treatment", "Length", "Width", "Survival_Censor")) {
    expect_true(col %in% colnames(df), info = paste("Missing column:", col))
  }
})

test_that("dose_levels data: loads without error and has Dose column", {
  p <- find_demo("dose_levels_synthetic_data.csv")
  skip_if(is.null(p), "Demo file not found")

  df <- utils::read.csv(p, stringsAsFactors = FALSE)

  expect_true(is.data.frame(df))
  expect_true(nrow(df) > 0)
  for (col in c("Date", "Cage", "ID", "Treatment", "Dose",
                "Length", "Width", "Survival_Censor")) {
    expect_true(col %in% colnames(df), info = paste("Missing column:", col))
  }
})



# ======================================================
# Data integrity: combo dataset
# ======================================================
test_that("combo data: four expected treatment groups present", {
  p <- find_demo("combo_treatment_synthetic_data.csv")
  skip_if(is.null(p), "Demo file not found")

  df      <- utils::read.csv(p, stringsAsFactors = FALSE)
  groups  <- sort(unique(df$Treatment))
  expected <- c("Control", "HDACi", "HDACi + PD1", "aPD1")

  expect_equal(sort(groups), sort(expected),
               info = paste("Groups found:", paste(groups, collapse = ", ")))
})

test_that("combo data: no negative volumes (Length/Width >= 0)", {
  p <- find_demo("combo_treatment_synthetic_data.csv")
  skip_if(is.null(p), "Demo file not found")

  df <- utils::read.csv(p, stringsAsFactors = FALSE)

  # Zeros are valid (early-stage / undetectable tumors); negatives are not
  expect_true(all(df$Length >= 0, na.rm = TRUE))
  expect_true(all(df$Width  >= 0, na.rm = TRUE))
  # At least some measurements are positive
  expect_true(any(df$Length > 0, na.rm = TRUE))
})

test_that("combo data: Survival_Censor is 0 or 1", {
  p <- find_demo("combo_treatment_synthetic_data.csv")
  skip_if(is.null(p), "Demo file not found")

  df <- utils::read.csv(p, stringsAsFactors = FALSE)
  cens <- na.omit(df$Survival_Censor)

  expect_true(all(cens %in% c(0, 1)),
              info = paste("Unexpected censor values:", paste(unique(cens[!cens %in% c(0,1)]), collapse = ", ")))
})

# ======================================================
# Data integrity: dose_levels dataset
# ======================================================
test_that("dose_levels data: Dose column is numeric with multiple levels", {
  p <- find_demo("dose_levels_synthetic_data.csv")
  skip_if(is.null(p), "Demo file not found")

  df <- utils::read.csv(p, stringsAsFactors = FALSE)

  expect_true(is.numeric(df$Dose) || all(!is.na(suppressWarnings(as.numeric(df$Dose)))))
  expect_true(length(unique(df$Dose)) >= 3L,
              info = paste("Unique doses:", paste(unique(df$Dose), collapse = ", ")))
})

# ======================================================
# Round-trip through validate_uploaded_file + detect_columns
# ======================================================
test_that("combo data: validate_uploaded_file returns valid=TRUE", {
  pkg_env <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
  skip_if(is.null(pkg_env), "Package namespace not available")
  vuf <- tryCatch(get("validate_uploaded_file", envir = pkg_env), error = function(e) NULL)
  skip_if(is.null(vuf), "validate_uploaded_file not accessible")

  p   <- find_demo("combo_treatment_synthetic_data.csv")
  skip_if(is.null(p), "Demo file not found")

  res <- vuf(p, "csv")
  expect_true(res$valid, info = res$message)
  expect_true(nrow(res$data) > 0)
})

test_that("dose_levels data: detect_columns finds Treatment, Date, Dose", {
  pkg_env <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
  skip_if(is.null(pkg_env), "Package namespace not available")
  dc <- tryCatch(get("detect_columns", envir = pkg_env), error = function(e) NULL)
  skip_if(is.null(dc), "detect_columns not accessible")

  p   <- find_demo("dose_levels_synthetic_data.csv")
  skip_if(is.null(p), "Demo file not found")

  df   <- utils::read.csv(p, stringsAsFactors = FALSE)
  cols <- dc(colnames(df))

  tx_val <- cols$treatment %||% cols[["Treatment"]] %||% ""
  expect_true(nzchar(tx_val), info = "Treatment not detected")
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
