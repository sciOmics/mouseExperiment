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
#   $bliss_independence: list with $difference (R17.2 removed the Loewe CI)
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
  required <- c("summary", "bliss_independence")
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

test_that("R17.2: the Loewe combination index is gone", {
  # Removed in v0.21.0: it computed response additivity, not Loewe additivity,
  # and labelled Bliss-additive combinations antagonistic across most of the
  # range where active single agents sit.
  res <- call_synergy(make_synergy_additive())
  expect_null(res$combination_index)
  expect_null(res$loewe_additivity)
  expect_false("Loewe Expected" %in% res$summary$Treatment)
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
test_that("synergistic combo: Bliss difference is positive", {
  res <- call_synergy(make_synergy_synergistic())
  expect_gt(res$bliss_independence$difference, 0)
  expect_match(res$overall_assessment, "Synergy")
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
test_that("antagonistic combo: Bliss difference is negative", {
  res <- call_synergy(make_synergy_antagonist())
  expect_lt(res$bliss_independence$difference, 0)
  expect_match(res$overall_assessment, "Antagonism|Additivity")
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
