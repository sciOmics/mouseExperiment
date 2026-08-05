# =============================================================================
# Round 15 — survival robustness and per-animal baselines
#
# Both findings are reachable from the dashboard with ordinary data. Neither
# produced an error message that pointed anywhere near its cause: one failed with
# "argument is of length zero", the other silently returned a smaller number.
# =============================================================================

surv_df <- function(with_cage = TRUE, seed = 31, n = 30) {
  set.seed(seed)
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) {
    lam <- if (tx == "Control") 0.10 else 0.05
    for (m in seq_len(n)) {
      k <- k + 1L
      t_ev <- stats::rexp(1, lam)
      row <- data.frame(ID = paste0(tx, m), Treatment = tx,
                        Time = min(t_ev, 60),
                        Event = as.integer(t_ev <= 60),
                        stringsAsFactors = FALSE)
      if (with_cage) row$Cage <- paste0("C", ceiling(m / 4))
      rows[[k]] <- row
    }
  }
  do.call(rbind, rows)
}
surv_run <- function(df, cage) {
  suppressWarnings(suppressMessages(survival_statistics(
    df, time_column = "Time", censor_column = "Event",
    treatment_column = "Treatment", id_column = "ID", cage_column = cage,
    reference_group = "Control", verbose = FALSE)))
}

# ---- R15.1 survival could not run without a cage column ---------------------

test_that("R15.1: survival runs when no cage column is mapped", {
  # The dashboard passes cage_column = NULL whenever the user has not mapped a
  # cage. That hit `if (NULL %in% colnames(df))`, which is `if (logical(0))` and
  # errors with "argument is of length zero" -- so the Survival tab failed for
  # every upload without a cage column.
  expect_s3_class(surv_run(surv_df(with_cage = FALSE), NULL)$results, "data.frame")
})

test_that("R15.1: a cage_column naming an absent column degrades, not errors", {
  # The default is "Cage", which plenty of real uploads lack. A bare data-frame
  # subset then gave "undefined columns selected".
  r <- surv_run(surv_df(with_cage = FALSE), "Cage")
  expect_s3_class(r$results, "data.frame")
  expect_identical(r$cage_structure$structure, "none")
})

test_that("R15.1: results are identical however the absent cage is expressed", {
  a <- surv_run(surv_df(with_cage = FALSE), NULL)$results
  b <- surv_run(surv_df(with_cage = FALSE), "Cage")$results
  expect_equal(a$HR, b$HR)
})

test_that("R15.1: a real cage column is still used, not ignored by the guard", {
  # The fix must not silently discard cage information that IS present.
  r <- surv_run(surv_df(with_cage = TRUE), "Cage")
  expect_false(identical(r$cage_structure$structure, "none"))
})

test_that("R15.1: the hazard ratio recovers a known truth", {
  # Exponential survival at 0.10/day against 0.05/day => true HR 0.5.
  r <- surv_run(surv_df(with_cage = TRUE, n = 60), "Cage")
  row <- r$results[r$results$Group == "DrugA", ]
  expect_gt(row$HR, 0.3); expect_lt(row$HR, 0.8)
  expect_lt(row$CI_Lower, 0.5); expect_gt(row$CI_Upper, 0.5)
})

# ---- R15.2 baselines anchored to the global first day -----------------------

# Two DrugA animals enrol late (no day-0 weight) and are far more toxic. If they
# are dropped, the reported mean loss collapses to the value from the other four.
stagger_df <- function() {
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) for (m in 1:6) {
    k <- k + 1L
    late <- (tx == "DrugA" && m <= 2)
    days <- if (late) seq(3, 21, 3) else seq(0, 21, 3)
    drop <- if (tx != "DrugA") 0.02 else if (late) 0.30 else 0.10
    rows[[k]] <- data.frame(
      ID = paste0(tx, m), Treatment = tx, Day = days,
      Weight = 25 * (1 - drop * days / 21),
      Volume = exp(log(150) + (if (tx == "Control") 0.12 else 0.07) * days),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

test_that("R15.2: late-enrolling animals are not dropped from weight loss", {
  # The defect returned exactly 10.0 -- the mean of the four animals that had a
  # day-0 weight -- while silently excluding the two most-toxic ones. The bias
  # has a direction: excluded animals leave the toxicity denominator, so the
  # therapeutic window reads safer than it is.
  r <- suppressWarnings(suppressMessages(therapeutic_window_metric(
    stagger_df(), volume_column = "Volume", weight_column = "Weight",
    time_column = "Day", treatment_column = "Treatment", id_column = "ID",
    reference_group = "Control", adjust_tumor_weight = FALSE, n_boot = 0)))
  got <- r$twm_table$Mean_Pct_Weight_Loss[r$twm_table$Treatment == "DrugA"]
  expect_gt(got, 13)          # the old, wrong answer was exactly 10.0
  expect_lt(got, 18)
})

test_that("R15.2: baseline is each animal's own first observation", {
  # Not the group's first day. The late animals' true observable loss is measured
  # from day 3 (26.9%), not from a day-0 weight they never had (30%) -- so the
  # correct group mean is ~15.6, and a test expecting 16.7 would be asserting an
  # unobservable quantity.
  b <- me_per_mouse_baseline(stagger_df(), c("ID", "Treatment"), "Weight")
  late <- b$Baseline_Weight[b$ID == "DrugA1"]
  expect_equal(late, 25 * (1 - 0.30 * 3 / 21), tolerance = 1e-8)
  # And an animal present at day 0 still baselines there.
  expect_equal(b$Baseline_Weight[b$ID == "DrugA3"], 25, tolerance = 1e-8)
  expect_identical(nrow(b), 12L)
})

test_that("R15.2: complete data is unaffected by the change", {
  # Regression guard: when every animal has the same first day, per-animal and
  # global baselines coincide, so nothing should move.
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) for (m in 1:6) {
    k <- k + 1L
    days <- seq(0, 21, 3)
    drop <- if (tx == "DrugA") 0.10 else 0.02
    rows[[k]] <- data.frame(ID = paste0(tx, m), Treatment = tx, Day = days,
      Weight = 25 * (1 - drop * days / 21),
      Volume = exp(log(150) + (if (tx == "Control") 0.12 else 0.07) * days),
      stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, rows)
  b <- me_per_mouse_baseline(df, c("ID", "Treatment"), "Weight")
  expect_true(all(abs(b$Baseline_Weight - 25) < 1e-8))
})

test_that("R15.2: the helper averages duplicate rows on the first day", {
  # Matching the previous behaviour for repeated measurements on one day.
  df <- data.frame(ID = c("m1", "m1", "m1"), Treatment = "A",
                   Day = c(0, 0, 7), Weight = c(24, 26, 20),
                   stringsAsFactors = FALSE)
  b <- me_per_mouse_baseline(df, c("ID", "Treatment"), "Weight")
  expect_equal(b$Baseline_Weight, 25, tolerance = 1e-8)
})
