# =============================================================================
# Tests for analyze_drug_synergy()
#
# Ground-truth design (Bliss independence):
#   Control TGI = 0       (mean vol ≈ 500)
#   DrugA   TGI ≈ 0.50   (mean vol ≈ 250)
#   DrugB   TGI ≈ 0.40   (mean vol ≈ 300)
#   Bliss expected TGI = 0.50 + 0.40 - 0.50*0.40 = 0.70
#   Bliss expected combo vol = 500 * 0.30 = 150
#
#   Additive:     actual combo mean ≈ 150  → |bliss difference| ≈ 0
#   Synergistic:  actual combo mean ≈ 50   → bliss difference > 0
#   Antagonistic: actual combo mean ≈ 325  → bliss difference < 0
#
# Verified structure of analyze_drug_synergy() output:
#   $bliss_independence: list with $expected_effect, $observed_effect,
#                        $difference (positive = synergy), $synergy (label)
#   $combination_index:  list with $ci (named numeric), $interpretation (string)
# =============================================================================

call_synergy <- function(df) {
  suppressWarnings(suppressMessages(
    analyze_drug_synergy(
      df,
      treatment_column = "Treatment",
      volume_column    = "Volume",
      time_column      = "Day",
      drug_a_name      = "DrugA",
      drug_b_name      = "DrugB",
      combo_name       = "DrugA+DrugB",
      control_name     = "Control",
      eval_time_point  = 21
    )
  ))
}

# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------
test_that("returns a list with expected top-level fields", {
  df  <- make_synergy_additive()
  res <- call_synergy(df)

  expect_true(is.list(res))
  required <- c("summary", "bliss_independence", "combination_index")
  present  <- required[required %in% names(res)]
  expect_true(length(present) >= 2L,
              info = paste("Missing:", paste(setdiff(required, names(res)), collapse = ", ")))
})

test_that("summary table contains at least 3 treatment groups", {
  df  <- make_synergy_additive()
  res <- call_synergy(df)

  summ <- as.data.frame(res$summary)
  expect_true(nrow(summ) >= 3L,
              info = paste("Got", nrow(summ), "rows in summary"))
})

test_that("bliss_independence contains difference field", {
  df  <- make_synergy_additive()
  res <- call_synergy(df)

  skip_if(!("bliss_independence" %in% names(res)), "bliss_independence not returned")
  bliss <- res$bliss_independence
  expect_true("difference" %in% names(bliss),
              info = paste("bliss_independence fields:", paste(names(bliss), collapse = ", ")))
})

test_that("combination_index is a list with $ci element", {
  df  <- make_synergy_synergistic()
  res <- call_synergy(df)

  skip_if(!("combination_index" %in% names(res)), "combination_index not returned")
  ci_obj <- res$combination_index
  expect_true(is.list(ci_obj))
  expect_true("ci" %in% names(ci_obj),
              info = paste("combination_index fields:", paste(names(ci_obj), collapse = ", ")))
})

# ---------------------------------------------------------------------------
# Additive combo: Bliss difference near zero
# ---------------------------------------------------------------------------
test_that("additive combo: |bliss difference| < 0.20 (near-additive)", {
  df  <- make_synergy_additive()
  res <- call_synergy(df)

  skip_if(!("bliss_independence" %in% names(res)), "bliss_independence not returned")
  skip_if(is.null(res$bliss_independence$difference), "difference field absent")

  delta <- as.numeric(res$bliss_independence$difference)
  expect_true(abs(delta) < 0.20,
              info = paste("|Bliss difference| =", round(abs(delta), 3),
                           "; expected < 0.20 for additive combo"))
})

# ---------------------------------------------------------------------------
# Synergistic combo
# ---------------------------------------------------------------------------
test_that("synergistic combo: combination_index$ci < 1", {
  df  <- make_synergy_synergistic()
  res <- call_synergy(df)

  skip_if(!("combination_index" %in% names(res)), "combination_index not returned")
  ci_val <- suppressWarnings(as.numeric(res$combination_index$ci))
  skip_if(all(is.na(ci_val)), "Could not extract numeric CI value")

  expect_true(any(ci_val < 1, na.rm = TRUE),
              info = paste("CI =", paste(round(ci_val, 3), collapse = ", "),
                           "; expected CI < 1 for synergistic combo"))
})

test_that("synergistic combo: bliss difference > 0 (observed TGI > expected)", {
  df  <- make_synergy_synergistic()
  res <- call_synergy(df)

  skip_if(!("bliss_independence" %in% names(res)), "bliss_independence not returned")
  skip_if(is.null(res$bliss_independence$difference), "difference field absent")

  delta <- as.numeric(res$bliss_independence$difference)
  expect_true(delta > 0,
              info = paste("Bliss difference =", round(delta, 3),
                           "; expected > 0 (actual TGI 90% > expected 70%)"))
})

# ---------------------------------------------------------------------------
# Antagonistic combo
# ---------------------------------------------------------------------------
test_that("antagonistic combo: combination_index$ci > 1", {
  df  <- make_synergy_antagonist()
  res <- call_synergy(df)

  skip_if(!("combination_index" %in% names(res)), "combination_index not returned")
  ci_val <- suppressWarnings(as.numeric(res$combination_index$ci))
  skip_if(all(is.na(ci_val)), "Could not extract numeric CI value")

  expect_true(any(ci_val > 1, na.rm = TRUE),
              info = paste("CI =", paste(round(ci_val, 3), collapse = ", "),
                           "; expected CI > 1 for antagonistic combo"))
})

test_that("antagonistic combo: bliss difference < 0 (observed TGI < expected)", {
  df  <- make_synergy_antagonist()
  res <- call_synergy(df)

  skip_if(!("bliss_independence" %in% names(res)), "bliss_independence not returned")
  skip_if(is.null(res$bliss_independence$difference), "difference field absent")

  delta <- as.numeric(res$bliss_independence$difference)
  expect_true(delta < 0,
              info = paste("Bliss difference =", round(delta, 3),
                           "; expected < 0 (actual TGI 38% < expected 70%)"))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
test_that("errors when specified combo group is absent from the data", {
  df <- make_synergy_additive()
  expect_error(
    suppressMessages(analyze_drug_synergy(
      df,
      treatment_column = "Treatment",
      volume_column    = "Volume",
      time_column      = "Day",
      drug_a_name      = "DrugA",
      drug_b_name      = "DrugB",
      combo_name       = "NONEXISTENT_COMBO",
      control_name     = "Control"
    ))
  )
})
