# =============================================================================
# Round 20 Minor and efficiency items closed in v0.32.0 (CODE_REVIEW.md R20-V).
# Each test fails on v0.31.0 for the reason its item describes.
# =============================================================================

m <- master_synthetic_data
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
tg_fit <- function(d = m, ...) {
  quiet(tumor_growth_statistics(d, reference_group = "Vehicle", plots = FALSE,
                                include_diagnostics = FALSE, ...))
}

# ---- R20.59: Dunnett p-values are reproducible and leave the RNG alone ----------

test_that("R20.59: Dunnett (mvt) summaries repeat exactly, and the caller's RNG is untouched", {
  a <- tg_fit(p_adjust_method = "dunnett")$pairwise_comparisons
  b <- tg_fit(p_adjust_method = "dunnett")$pairwise_comparisons
  expect_identical(a$P_Value_Adjusted, b$P_Value_Adjusted)
  expect_true(all(a$P_Value_Adjusted >= a$P_Value_Raw))
  set.seed(11); before <- .Random.seed
  emm <- emmeans::emmeans(lm(Sepal.Length ~ Species, iris), ~ Species)
  pc  <- emmeans::contrast(emm, "trt.vs.ctrl", adjust = "mvt")
  s1 <- me_contrast_summary(pc)
  expect_identical(.Random.seed, before)
  expect_identical(me_contrast_summary(pc)$p.value, s1$p.value)
})

# ---- R20.65: weighings before the baseline are not at risk -------------------------

test_that("R20.65: a low weight before the baseline day is not a weight-loss event", {
  d <- data.frame(
    ID = rep(1:6, each = 4), Treatment = rep(c("Vehicle", "Drug"), each = 12),
    Day = rep(c(-7, 0, 7, 14), 6),
    Weight = c(15, 20, 20.5, 21,  20, 20, 20, 20,  20, 20.2, 20.1, 20,
               20, 20, 19.9, 20,  20, 20.1, 20, 20.2,  20, 20, 20, 15.0))
  r <- quiet(weight_loss_threshold(d, threshold = 0.2, baseline_day = 0,
                                   reference_group = "Vehicle"))
  ev <- r$event_data
  # v0.31.0: animal 1 had an "event" at day -7, before its day-0 baseline.
  expect_equal(ev$Event[ev$ID == 1 & ev$Treatment == "Vehicle"], 0L)
  expect_equal(sum(ev$Event), 1L)                  # animal 6, at day 14
  expect_true(all(ev$Time[ev$Event == 1L] >= 0))
  # The cumulative incidence lists the reference first, as the KM fit does.
  expect_match(names(r$cuminc$strata)[1], "Vehicle")
})

# ---- R20.68: the cage ICC does not depend on the time units -------------------------

test_that("R20.68: the cage ICC is the same with time in days or weeks", {
  wk <- m; wk$Day <- wk$Day / 7
  a <- tg_fit()$cage_analysis$icc
  b <- tg_fit(wk)$cage_analysis$icc
  expect_false(is.null(a))
  expect_equal(a$ICC, b$ICC, tolerance = 1e-3)   # 0.513 vs 0.395 before
})

# ---- R20.71: one definition of n in the analytic power analysis ----------------------

test_that("R20.71: the power curve and the find_power table agree; enrolment is rounded", {
  fp <- quiet(apriori_power_analysis(mode = "find_power", effect_size = 1,
                                     n_per_group = 10, n_groups = 2, alpha = 0.05,
                                     dropout_rate = 0.3))
  curve <- fp$power_curve_data
  expect_equal(curve$Power[curve$N_Per_Group == 10], fp$scenario_table$Achieved_Power)
  fn <- quiet(apriori_power_analysis(mode = "find_n", delta = 300, pooled_sd = 265,
                                     n_groups = 2, alpha = 0.05, target_power = 0.8,
                                     dropout_rate = 0.3))
  expect_true("Enroll_N" %in% names(fn$sensitivity_table))
  expect_equal(fn$sensitivity_table$Enroll_N,
               ceiling(round(fn$sensitivity_table$Required_N / 0.7, 8)))
})

# ---- R20.66 / R20.75 / R20.69 / R20.58: counts, messages and metadata -----------------

test_that("R20.66: the dose-response summary counts measured volumes", {
  d <- data.frame(Dose = rep(c(0, 10), each = 5),
                  Volume = c(100, 120, NA, 90, 110,  60, NA, NA, 70, 65))
  st <- generate_summary_statistics(d, dose_column = "Dose", volume_column = "Volume",
                                    verbose = FALSE)
  expect_equal(as.integer(st$n), c(4L, 3L))          # 5 and 5 before
  expect_equal(st$sem_volume[1], stats::sd(c(100, 120, 90, 110)) / 2)
})

test_that("R20.75: the survivors warning says how many animals are still on study", {
  d <- m
  gone <- d$Treatment == "Vehicle" & d$ID == d$ID[d$Treatment == "Vehicle"][1] &
    d$Cage == d$Cage[d$Treatment == "Vehicle"][1] & d$Day > 8
  d <- d[!gone, ]                                  # one control animal leaves at day 8
  expect_warning(endpoint_volumes(d, cage_column = "Cage", endpoint_method = "survivors",
                                  endpoint_day = 16),
                 "Vehicle: 7 of 8 still on study")
})

test_that("R20.69: repeated animal-day rows draw a warning", {
  d <- rbind(m, m[m$ID == m$ID[1] & m$Treatment == m$Treatment[1] & m$Day == 0, ])
  expect_warning(suppressMessages(tumor_growth_statistics(
    d, reference_group = "Vehicle", plots = FALSE, include_diagnostics = FALSE)),
    "repeat an animal and day")
})

test_that("R20.58: the summary describes the model that was fitted", {
  tg <- tg_fit()
  ms <- tg$summary$model_specification
  expect_match(ms$random_effects, "\\(1 \\| Cage\\)")         # cage was left out
  expect_identical(tg$summary$data_description$subjects,
                   as.integer(lme4::ngrps(tg$model)[["ID"]]))
  d <- m; d$Cage[d$Treatment == "Drug_B"] <- NA
  expect_warning(suppressMessages(tumor_growth_statistics(
    d, reference_group = "Vehicle", plots = FALSE, include_diagnostics = FALSE)),
    "left out of the model")
})

# ---- R20.61: plots ---------------------------------------------------------------------

test_that("R20.61: the mean line ends when an arm is mostly gone; units and pre-treatment days", {
  d <- m
  d$Day[d$Day == 0] <- -3                          # a pre-treatment day
  p <- quiet(plot_tumor_growth(d, volume_units = "cm3"))
  b <- ggplot2::ggplot_build(p)
  expect_true(min(b$layout$panel_params[[1]]$x.range) <= -3)   # not dropped
  expect_match(paste(deparse(p$labels$y), collapse = ""), "cm\\^3")
  mean_layer <- which(vapply(p$layers, function(l) "Mean" %in% names(l$data), logical(1)))
  md <- p$layers[[mean_layer]]$data
  enrolled <- tapply(d$ID[d$Treatment == "Vehicle"], d$Cage[d$Treatment == "Vehicle"],
                     function(x) length(unique(x)))
  n_veh <- sum(enrolled)
  last_day <- max(md$Day[md$Group == "Vehicle"])
  n_left <- length(unique(paste(d$ID, d$Cage)[d$Treatment == "Vehicle" & d$Day == last_day &
                                                  is.finite(d$Volume)]))
  expect_gte(n_left, ceiling(n_veh / 2))
  expect_match(p$labels$caption, "Mean lines end")
})

test_that("R20.61: plot_auc's interval uses t; the caterpillar plot uses Satterthwaite df", {
  ci_fun <- function(x) {
    n <- sum(!is.na(x)); stats::qt(0.975, n - 1) * stats::sd(x) / sqrt(n)
  }
  tg <- quiet(tumor_growth_statistics(m, reference_group = "Vehicle", model_type = "auc",
                                      plots = FALSE, include_diagnostics = FALSE))
  ind <- tg$auc_analysis$individual
  p <- plot_auc(ind, error_bar_type = "CI")
  bars <- ggplot2::ggplot_build(p)$data
  eb <- bars[[which(vapply(bars, function(b) all(c("ymin", "ymax") %in% names(b)), logical(1)))[1]]]
  v <- ind$AUC[ind$Group == sort(unique(ind$Group))[1] & is.finite(ind$AUC)]
  expect_equal(eb$ymax[1] - eb$y[1], ci_fun(v), tolerance = 1e-6)

  fit <- tg_fit(return_model = TRUE)$model
  cp <- plot_caterpillar(fit)
  cd <- cp$data
  ct <- stats::coef(summary(lmerTest::as_lmerModLmerTest(fit)))
  k <- rownames(ct)[2]
  half <- ct[k, "Std. Error"] * stats::qt(0.975, ct[k, "df"])
  expect_equal(cd$ci_upper[cd$term == k] - cd$estimate[cd$term == k], half, tolerance = 1e-6)
})

# ---- R20.67: synergy over time ---------------------------------------------------------

test_that("R20.67: plot_bliss labels the analysed arms; day p-values are adjusted too", {
  tr <- quiet(analyze_drug_synergy_over_time(m, drug_a_name = "Drug_A Mid",
    drug_b_name = "Drug_B", combo_name = "Drug_A Mid + Drug_B", reference_group = "Vehicle",
    min_day = 2, n_boot = 50, seed = 1, verbose = FALSE))
  p <- plot_bliss(tr)
  lv <- levels(factor(p$data$Group))
  expect_true(all(c("Drug_A Mid", "Drug_B", "Drug_A Mid + Drug_B") %in% p$data$Group))
  ss <- tr$synergy_summary
  expect_equal(ss$P_Value_vs_Drug_A_Holm, stats::p.adjust(ss$P_Value_vs_Drug_A, "holm"))
  expect_false("ci_floor_frac" %in% names(formals(synergy_bootstrap)))
})

# ---- R20.82: exact randomisation tests refit each partition once -------------------------

test_that("R20.82: exhaustive enumeration refits partitions, not relabellings", {
  a <- enumerate_label_perms(rep(c("A", "B", "C"), each = 2))
  expect_length(a, 90L)
  expect_length(me_partition_representatives(a), 15L)
  d <- m[m$Treatment %in% c("Vehicle", "Drug_A High", "Drug_B"), ]
  r <- quiet(trajectory_permutation_test(d, spec = perm_spec("cage", n_perm = 999, seed = 1)))
  expect_true(r$exhaustive)
  expect_identical(r$n_perm, 15L)
  expect_equal(r$p_value, 1 / 15)                  # the same p as all 90 relabellings
})
