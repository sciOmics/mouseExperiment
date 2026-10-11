# =============================================================================
# Round 20 addendum Majors closed in v0.31.0 (CODE_REVIEW.md R20-U): the
# findings of the dashboard audit that live in this package. Each test fails on
# v0.30.0 for the reason its item describes.
# =============================================================================

m <- master_synthetic_data
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
tg_fit <- function(d = m, ...) {
  quiet(tumor_growth_statistics(d, reference_group = "Vehicle", plots = FALSE, ...))
}

# ---- R20.84: Q-Q plots on standardised values ---------------------------------

test_that("R20.84: the residual Q-Q is standardised, so y = x is the right reference", {
  tg <- tg_fit()
  qd <- ggplot2::ggplot_build(tg$diag_qq_plot)$data[[1]]
  expect_equal(stats::sd(qd$y), 1, tolerance = 1e-8)
  # v0.30.0 plotted raw residuals (SD about 0.16 here) against y = x: almost
  # no point lay near the line. Standardised, most do.
  near <- mean(abs(qd$y - qd$x) < 0.25)
  expect_gt(near, 0.5)
})

test_that("R20.84: the random-effects Q-Q is standardised, and skipped with no spread", {
  tg <- tg_fit()
  rd <- ggplot2::ggplot_build(tg$diag_re_qq_plot)$data[[1]]
  sds <- tapply(rd$y, rd$PANEL, stats::sd)
  expect_true(all(abs(sds - 1) < 1e-8))
  expect_null(me_qq_frame(rep(1e-17, 10)))       # a singular fit's BLUPs
  expect_null(me_qq_frame(c(1, NA)))
})

# ---- R20.85: intervals adjusted like the p-values ------------------------------

test_that("R20.85: lme4 contrasts carry intervals adjusted like their p-values", {
  for (adj in c("dunnett", "bonferroni", "none")) {
    pc <- tg_fit(p_adjust_method = adj, include_diagnostics = FALSE)$pairwise_comparisons
    expect_true(all(c("CI_Adjusted_Lower", "CI_Adjusted_Upper", "CI_Adjustment") %in% names(pc)),
                info = adj)
    # An interval excludes 0 exactly when the adjusted p is below 0.05.
    excl <- pc$CI_Adjusted_Lower > 0 | pc$CI_Adjusted_Upper < 0
    expect_identical(excl, pc$P_Value_Adjusted < 0.05, info = adj)
  }
  bonf <- tg_fit(p_adjust_method = "bonferroni", include_diagnostics = FALSE)$pairwise_comparisons
  none <- tg_fit(p_adjust_method = "none", include_diagnostics = FALSE)$pairwise_comparisons
  expect_true(all(bonf$CI_Adjusted_Upper - bonf$CI_Adjusted_Lower >
                    none$CI_Adjusted_Upper - none$CI_Adjusted_Lower))
  expect_match(bonf$CI_Adjustment[1], "bonferroni")
})

test_that("R20.85: AUC contrasts carry adjusted intervals that agree with p_adjusted", {
  pc <- tg_fit(model_type = "auc", p_adjust_method = "bonferroni",
               include_diagnostics = FALSE)$pairwise_comparisons
  expect_true(all(c("CI_Adjusted_Lower", "Ratio_CI_Adjusted_Lower") %in% names(pc)))
  # On the scale of the test the agreement is exact.
  excl <- pc$Ratio_CI_Adjusted_Lower > 1 | pc$Ratio_CI_Adjusted_Upper < 1
  expect_identical(excl, pc$P_Value_Adjusted < 0.05)
  # Wider than the unadjusted 95% interval.
  expect_true(all(pc$CI_Adjusted_Upper - pc$CI_Adjusted_Lower > pc$ci_upper - pc$ci_lower))
  expect_identical(me_auc_ci_adjustment("holm", 5)$level, 0.99)
  expect_match(me_auc_ci_adjustment("fdr", 5)$label, "no simultaneous interval")
})

# ---- R20.57: the data summary on the volume scale -------------------------------

test_that("R20.57: data_summary is on the volume scale, with a geometric mean under log", {
  tg <- tg_fit(include_diagnostics = FALSE)
  ds <- tg$data_summary
  raw <- m$Volume[m$Treatment == "Vehicle" & m$Day == max(m$Day[m$Treatment == "Vehicle"])]
  row <- ds[ds$Treatment == "Vehicle" & ds$Day == max(ds$Day[ds$Treatment == "Vehicle"]), ]
  expect_equal(row$mean_volume, mean(raw, na.rm = TRUE))         # v0.30.0: log scale
  expect_equal(row$geometric_mean_volume, exp(mean(log(pmax(raw[!is.na(raw)], 2)))))
  sq <- tg_fit(transform = "sqrt", include_diagnostics = FALSE)$data_summary
  expect_false("geometric_mean_volume" %in% names(sq))
  expect_equal(sq$mean_volume, ds$mean_volume)
})

test_that("R20.57: n counts measured volumes, not rows", {
  d <- m
  first_vehicle <- which(d$Treatment == "Vehicle" & d$Day == 0)[1:2]
  d$Volume[first_vehicle] <- NA
  ds <- tg_fit(d, include_diagnostics = FALSE)$data_summary
  n_rows <- sum(m$Treatment == "Vehicle" & m$Day == 0)
  expect_identical(as.integer(ds$n[ds$Treatment == "Vehicle" & ds$Day == 0]), n_rows - 2L)
})
