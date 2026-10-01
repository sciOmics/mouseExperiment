# =============================================================================
# Regression tests for CODE_REVIEW.md Round 1 and Round 2 fixes — closes K.11.
#
# K.11 flagged that the ~40 Round 1/2 fixes shipped without regression tests.
# Round 3 then demonstrated the cost twice over:
#
#   R3.35  bayesian_synergy() was non-functional from v0.4.6 to v0.9.0 -- five
#          releases -- because the G.6 refactor moved a variable and left both
#          consumers pointing at the old local.
#   R3.37  bayesian_therapeutic_window()'s scatter plot was always NULL, because
#          the v0.4.5 fix for J.2 introduced a reference to a data frame that
#          does not exist.
#
# Both were *fixes* that stopped working. Neither had a test. (Both functions
# were removed in v0.23.0, CODE_REVIEW.md R20-K.)
#
# These tests target the Round 1/2 fixes whose failure mode is SILENT — a wrong
# number or an empty result rather than an error. A fix that fails loudly does not
# need a regression test as urgently; one that degrades quietly does.
# =============================================================================

# ---- Round 1 §1.5 / §3.1: composite mouse keys ------------------------------

test_that("R1-1.5: composite keys survive separators that appear in real labels", {
  # The original bug used "_" as the separator, so a treatment named "Drug_A"
  # split into the wrong fields. The package standardised on "|||".
  k <- make_mouse_key("Drug_A", "01", "Cage_1")
  parts <- split_mouse_key(k)
  expect_equal(length(parts), 3L)
  expect_identical(parts[1], "Drug_A")
  expect_identical(parts[2], "01")
  expect_identical(parts[3], "Cage_1")

  # Round-trip through every character class that occurs in real study labels.
  for (lab in c("Drug A", "Drug-A", "Drug.A", "Drug_A", "aPD1 + HDACi", "A/B")) {
    expect_identical(split_mouse_key(make_mouse_key(lab, "1", "C1"))[1], lab,
                     info = paste("failed round-trip for:", lab))
  }
})

test_that("R1-1.8: same ID in different arms stays two distinct animals", {
  # The recurring bug class: aggregating by ID alone collapses mice that share a
  # numeric ear-tag across treatment groups. Fixed in body_weight_auc (v0.3.1),
  # weight_corrected_tgi (1.8), therapeutic_window_metric (J.14),
  # efficacy_toxicity_bivariate (J.18) and analyze_body_weight (R3.11); all but
  # therapeutic_window_metric and analyze_body_weight were removed in v0.23.0.
  k1 <- make_mouse_key("Control", "1", "C1")
  k2 <- make_mouse_key("DrugA",   "1", "C1")
  expect_false(identical(k1, k2))

  # And the same ID in different cages within one arm.
  k3 <- make_mouse_key("DrugA", "1", "C2")
  expect_false(identical(k2, k3))
  expect_equal(length(unique(c(k1, k2, k3))), 3L)
})

# ---- Round 1 §3.6: one trapezoidal AUC implementation -----------------------

test_that("R1-3.6: calculate_auc is the single AUC implementation and is correct", {
  # Three files had copy-pasted trap_auc() variants with different edge-case
  # handling. All now call calculate_auc(); check it against a hand computation.
  t <- c(0, 2, 4);  v <- c(0, 10, 20)      # two trapezoids: 10 + 30
  expect_equal(calculate_auc(t, v), 10 + 30)

  # Unsorted input must give the same answer -- the local copies did not all sort.
  expect_equal(calculate_auc(c(4, 0, 2), c(20, 0, 10)), 40)
  # NA pairs dropped, not propagated.
  expect_equal(calculate_auc(c(0, 2, 4, 6), c(0, 10, 20, NA)), 40)
  # Fewer than two usable points is NA, not an error or 0.
  expect_true(is.na(calculate_auc(0, 5)))
  expect_error(calculate_auc(c(0, 1), 5), regexp = "same length")
})

# ---- Round 1 §2.5: LOCF AUC is an area, not a volume ------------------------


# ---- Round 1 §1.7: baselines must come from the earliest day ----------------

test_that("R1-1.7: baseline is the earliest day regardless of row order", {
  # aggregate(x[1]) took whatever row happened to be first. Shuffling the rows
  # must not change the answer.
  # Three animals: the evaluable-day rule (v0.26.0) needs at least 3 per arm.
  mk <- function(order_idx) {
    d <- data.frame(
      ID = rep(c("m1", "m2", "m3"), each = 3), Treatment = "A", Cage = "C1",
      Day = rep(c(0, 7, 14), 3),
      Weight = c(20, 19, 18, 22, 21, 20, 21, 20, 19),
      Volume = rep(c(100, 200, 300), 3),
      stringsAsFactors = FALSE)
    d[order_idx, , drop = FALSE]
  }
  ordered_res <- suppressWarnings(suppressMessages(therapeutic_window_metric(
    mk(1:9), reference_group = "A", adjust_tumor_weight = FALSE, n_boot = 0)))
  shuffled_res <- suppressWarnings(suppressMessages(therapeutic_window_metric(
    mk(c(3, 6, 9, 1, 4, 7, 2, 5, 8)), reference_group = "A",
    adjust_tumor_weight = FALSE, n_boot = 0)))

  expect_equal(ordered_res$weight_loss_data$Baseline_Weight,
               shuffled_res$weight_loss_data$Baseline_Weight)
  # Baselines must be the day-0 weights, not some later row.
  expect_setequal(round(ordered_res$weight_loss_data$Baseline_Weight, 6), c(20, 22, 21))
})

# ---- Round 1 §2.4 / §I.1: one Bliss implementation --------------------------

test_that("R1-2.4: Bliss expectation is shared and behaves at its ceiling", {
  # The formula was duplicated between the frequentist and Bayesian paths.
  expect_equal(synergy_bliss_expected(0.5, 0.5), 0.75)
  expect_equal(synergy_bliss_expected(0, 0.4), 0.4)
  # Vectorised, because the Bayesian path applies it across posterior draws.
  expect_equal(synergy_bliss_expected(c(0.5, 0.2), c(0.5, 0.5)), c(0.75, 0.6))

  # The documented ceiling effect: two strong single agents leave almost no room
  # to demonstrate synergy. This is the limitation §2.4 required be documented,
  # pinned so the formula cannot drift away from it.
  expect_gt(synergy_bliss_expected(0.9, 0.9), 0.98)
})

# ---- Round 2 §J.2 / R3.37: plots must actually be produced ------------------

test_that("R2-J.2/R3.37: plot slots hold plots, not silent NULLs", {
  # R3.37: a tryCatch swallowed an undefined-variable error and the plot was
  # always NULL. The general lesson -- assert that a documented plot field is a
  # ggplot, not merely that the call returned.
  set.seed(5)
  d <- data.frame(
    ID = rep(1:8, each = 4), Cage = rep(1:4, each = 8),
    Treatment = rep(c("Control", "DrugA"), each = 16),
    Day = rep(c(0, 7, 14, 21), 8), stringsAsFactors = FALSE)
  d$Volume <- exp(log(200) + 0.08 * d$Day + stats::rnorm(32, 0, 0.2))

  res <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    d, plots = FALSE, verbose = FALSE, include_diagnostics = TRUE)))

  for (slot in c("diag_qq_plot", "diag_resid_fitted_plot",
                 "diag_scale_location_plot")) {
    expect_s3_class(res[[slot]], "ggplot")
  }
})

# ---- Round 1 §3.3 / J.20: no parent-frame magic -----------------------------

test_that("R1-3.3/J.20: in_place does not mutate the caller's environment", {
  # Both calculate_volume() and calculate_dates() used assign() into
  # parent.frame(), which silently did the wrong thing when called from inside
  # another function. Both were deprecated to plain returns.
  df_local <- data.frame(Length = c(10, 12), Width = c(5, 6),
                         stringsAsFactors = FALSE)
  before <- df_local

  out <- suppressWarnings(calculate_volume(df_local, in_place = TRUE))

  # The caller's object must be untouched; the result comes back as a value.
  expect_identical(df_local, before)
  expect_true("Volume" %in% names(out))
  expect_false("Volume" %in% names(df_local))
})

# ---- Round 2 §A.1: the PH assumption is actually tested ---------------------

test_that("R2-A.1: survival_statistics returns a proportional-hazards check", {
  set.seed(8)
  sdf <- data.frame(
    ID = 1:24, Cage = rep(1:6, each = 4),
    Treatment = rep(c("Control", "DrugA"), each = 12),
    Day = c(stats::rnorm(12, 16, 3), stats::rnorm(12, 26, 4)),
    Survival_Censor = c(rep(1L, 11), 0L, rep(1L, 11), 0L),
    stringsAsFactors = FALSE)

  res <- suppressWarnings(suppressMessages(
    survival_statistics(sdf, reference_group = "Control", verbose = FALSE)))

  # A.1 added cox.zph(); E.2 added concordance. Both must survive on the Cox path.
  expect_identical(res$method_used, "cox")
  expect_false(is.null(res$ph_test))
  expect_true("GLOBAL" %in% rownames(res$ph_test$table))
  expect_false(is.null(res$c_index))
})
