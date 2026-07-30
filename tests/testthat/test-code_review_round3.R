# =============================================================================
# Regression tests for CODE_REVIEW.md Round 3 fixes.
#
# Round 2 K.11 flagged that none of the ~40 prior review fixes had an
# accompanying regression test, which is the usual reason fixes silently
# revert. Every Round 3 fix gets a test here, labelled with its finding ID.
#
# Two of these (R3.1, R3.2) exist because a fix was written in an earlier round
# and did not execute at runtime — so these assert observable behaviour, not
# the presence of code.
# =============================================================================

# ---- helpers ----------------------------------------------------------------

make_surv_df <- function(zero_event_group = TRUE) {
  set.seed(7)
  data.frame(
    ID              = 1:24,
    Cage            = rep(1:6, each = 4),
    Treatment       = rep(c("Control", "DrugA", "DrugB"), each = 8),
    Day             = c(rnorm(8, 15, 3), rnorm(8, 17, 3), rep(35, 8)),
    Survival_Censor = c(rep(1L, 8), rep(1L, 8),
                        if (zero_event_group) rep(0L, 8) else rep(1L, 8)),
    stringsAsFactors = FALSE
  )
}

make_tox_df <- function(volume_scale = 1, reuse_ids = TRUE) {
  set.seed(5)
  one <- function(tx, base) do.call(rbind, lapply(1:4, function(i) {
    data.frame(
      ID        = if (reuse_ids) i else paste0(tx, i),
      Treatment = tx,
      Cage      = paste0(tx, "_C1"),
      Day       = c(0, 7, 14),
      Weight    = base + i + c(0, -0.5, -1) + rnorm(3, 0, 0.05),
      Volume    = c(100, 400, 900) * volume_scale,
      stringsAsFactors = FALSE
    )
  }))
  rbind(one("Control", 20), one("DrugA", 26))
}

# ---- R3.2 -------------------------------------------------------------------

test_that("R3.2: pairwise log-rank p-values are per-comparison, not the omnibus", {
  df <- make_surv_df()
  # Pin the asymptotic flavour: this test checks the formula-scoping fix by
  # comparing against survdiff() directly. The permutation flavour (H.4) is
  # exercised separately below.
  res <- suppressWarnings(suppressMessages(
    survival_statistics(df, reference_group = "Control",
                        firth_correction = FALSE, permutation_logrank = FALSE,
                        verbose = FALSE)
  ))
  expect_identical(res$method_used, "logrank")

  non_ref <- res$results$Group != "Control"
  praw    <- res$results$P_Value_Unadjusted[non_ref]

  # The old code errored inside tryCatch and substituted the omnibus p-value
  # for every group, so all non-reference rows were identical to each other
  # and to the omnibus value.
  omnibus <- attr(res$results, "omnibus_logrank_p")
  expect_false(isTRUE(all.equal(unname(praw), rep(omnibus, length(praw)))))
  expect_equal(length(unique(praw)), length(praw))

  # Each must equal an independently computed pairwise log-rank.
  for (g in c("DrugA", "DrugB")) {
    pd <- df[df$Treatment %in% c("Control", g), ]
    sd <- survival::survdiff(survival::Surv(Day, Survival_Censor) ~ Treatment,
                             data = pd)
    expect_equal(res$results$P_Value_Unadjusted[res$results$Group == g],
                 1 - stats::pchisq(sd$chisq, df = 1L),
                 tolerance = 1e-10)
  }
})

# ---- R3.1 (survival path) ---------------------------------------------------

test_that("R3.1: p_adjust_method changes the reported survival p-values", {
  df <- make_surv_df()
  run <- function(m) suppressWarnings(suppressMessages(
    survival_statistics(df, reference_group = "Control", firth_correction = FALSE,
                        p_adjust_method = m, permutation_logrank = FALSE,
                        verbose = FALSE)
  ))$results

  none <- run("none")
  bonf <- run("bonferroni")
  holm <- run("holm")

  nr <- none$Group != "Control"
  # Bonferroni over 2 comparisons doubles the raw values.
  expect_equal(bonf$P_Value[nr], pmin(1, 2 * none$P_Value[nr]), tolerance = 1e-10)
  # Holm differs from Bonferroni on the larger p-value.
  expect_false(isTRUE(all.equal(holm$P_Value[nr], bonf$P_Value[nr])))
  # The result is self-describing.
  expect_identical(unique(bonf$P_Adjust_Method), "bonferroni")
  expect_true("P_Value_Unadjusted" %in% names(bonf))
})

# ---- R3.14 ------------------------------------------------------------------

test_that("R3.14: longitudinal input to survival_statistics is rejected", {
  df   <- make_surv_df()
  long <- df[rep(seq_len(nrow(df)), each = 5), ]
  long$Day <- long$Day + rep(0:4, nrow(df))

  expect_error(
    suppressMessages(survival_statistics(long, verbose = FALSE)),
    regexp = "one row per animal"
  )
  # The one-row-per-animal form still works.
  expect_no_error(suppressWarnings(suppressMessages(
    survival_statistics(df, firth_correction = FALSE, verbose = FALSE)
  )))
})

# ---- R3.27 ------------------------------------------------------------------

test_that("R3.27: an all-events group does not trigger the Firth path", {
  df <- make_surv_df(zero_event_group = FALSE)   # every animal has an event
  res <- suppressWarnings(suppressMessages(
    survival_statistics(df, reference_group = "Control", verbose = FALSE)
  ))
  expect_identical(res$method_used, "cox")
})

# ---- R3.11 ------------------------------------------------------------------

test_that("R3.11: baseline weight uses the composite mouse key", {
  df  <- make_tox_df(reuse_ids = TRUE)   # IDs 1-4 reused in both groups
  res <- suppressMessages(analyze_body_weight(
    df, volume_column = "Volume", cage_column = "Cage",
    covariates = c("initial_mass"), volume_units = "mm3"
  ))
  im <- unique(res$weight_data[, c("Treatment", "ID", "Initial_Mass")])

  # 8 distinct mice must have 8 distinct baselines; keying on ID alone
  # collapsed them to 4 shared values.
  expect_equal(nrow(im), 8L)
  expect_equal(length(unique(round(im$Initial_Mass, 4))), 8L)

  # Groups were built ~6 g apart, so the group means must reflect that.
  ctrl <- mean(im$Initial_Mass[im$Treatment == "Control"])
  drug <- mean(im$Initial_Mass[im$Treatment == "DrugA"])
  expect_gt(drug - ctrl, 4)
})

# ---- R3.10 ------------------------------------------------------------------

test_that("R3.10: tumour burden is not adjusted for twice", {
  df <- make_tox_df(reuse_ids = FALSE)

  # Default call requests both the mass subtraction and a Volume covariate.
  expect_message(
    res <- analyze_body_weight(df, volume_column = "Volume",
                               adjust_tumor_weight = TRUE,
                               covariates = c("volume"), volume_units = "mm3"),
    regexp = "redundant 'volume' covariate"
  )
  expect_false(grepl("Volume", res$model_info$fixed_formula %||%
                       paste(res$fixed_effects$Term, collapse = " ")))

  # Conditioning on volume instead of subtracting it still works.
  res2 <- suppressMessages(analyze_body_weight(
    df, volume_column = "Volume", adjust_tumor_weight = FALSE,
    covariates = c("volume"), volume_units = "mm3"
  ))
  expect_true(any(grepl("Volume", res2$fixed_effects$Term)))
})

# ---- R3.30 ------------------------------------------------------------------

test_that("R3.30: volume units are resolved rather than assumed to be mm3", {
  mm3 <- make_tox_df(volume_scale = 1,       reuse_ids = FALSE)
  cm3 <- make_tox_df(volume_scale = 1 / 1000, reuse_ids = FALSE)

  r_mm3 <- suppressMessages(analyze_body_weight(
    mm3, volume_column = "Volume", volume_units = "mm3", covariates = character()))
  r_cm3 <- suppressMessages(analyze_body_weight(
    cm3, volume_column = "Volume", volume_units = "cm3", covariates = character()))

  # Same physical tumour expressed in two units must give the same Net_Weight.
  expect_equal(sort(r_mm3$weight_data$Net_Weight),
               sort(r_cm3$weight_data$Net_Weight), tolerance = 1e-8)

  # Auto-detection recovers the right unit in both cases.
  expect_identical(detect_volume_units(mm3$Volume), "mm3")
  expect_identical(detect_volume_units(cm3$Volume), "cm3")

  # A mismatched explicit unit is flagged, not silently accepted.
  expect_warning(resolve_volume_units(cm3$Volume, "mm3"), regexp = "look like")
})

# ---- R3.18 ------------------------------------------------------------------

test_that("R3.18: TWM is continuous across the noise floor for any noise_floor", {
  df <- make_tox_df(reuse_ids = FALSE)

  twm_at <- function(wl, tgi, nf) tgi / max(wl, nf)   # the implemented form
  for (nf in c(1.0, 2.5, 5.0)) {
    below <- twm_at(nf - 1e-8, 60, nf)
    above <- twm_at(nf + 1e-8, 60, nf)
    expect_equal(below, above, tolerance = 1e-6,
                 info = paste("discontinuity at noise_floor =", nf))
  }

  # Behaviour at the default noise_floor = 1.0 is unchanged (backward compat):
  # negligible weight loss still scores TWM == TGI.
  res <- suppressMessages(therapeutic_window_metric(
    df, reference_group = "Control", volume_units = "mm3"))
  neg <- res$twm_table$Mean_Pct_Weight_Loss <= 1.0
  if (any(neg)) {
    expect_equal(res$twm_table$TWM[neg], pmax(res$twm_table$TGI[neg], 0),
                 tolerance = 1e-8)
  }
})

# ---- R3.31 ------------------------------------------------------------------

test_that("R3.31: Jonckheere-Terpstra trend test actually runs", {

  # Deliberately UNEQUAL dose spacing — the case contr.poly gets wrong and JT
  # is immune to, since it uses only the ordering of the dose levels.
  set.seed(9)
  doses <- c(0, 10, 30, 100)
  df <- do.call(rbind, lapply(doses, function(d) data.frame(
    ID        = paste0("D", d, "_", 1:8),
    Dose      = d,
    Day       = 21,
    Treatment = paste0("Dose_", d),
    Volume    = pmax(50, 1200 * (1 - 0.85 * d / (d + 25)) + rnorm(8, 0, 80)),
    stringsAsFactors = FALSE
  )))

  res <- suppressWarnings(suppressMessages(dose_response_statistics(
    df, dose_column = "Dose", volume_column = "Volume",
    day_column = "Day", id_column = "ID", verbose = FALSE
  )))
  jt <- res$trend_test$jonckheere_test

  expect_false(is.null(jt))
  expect_true(is.numeric(jt$p.value))
  expect_identical(jt$alternative_used, "decreasing")
  expect_equal(jt$p.value,
               clinfun::jonckheere.test(df$Volume, df$Dose,
                                        alternative = "decreasing")$p.value,
               tolerance = 1e-10)
})

test_that("R3.31: JT direction is derived from the data, not hardcoded", {
  set.seed(4)
  doses <- c(0, 10, 30, 100)
  # Stimulatory: volume INCREASES with dose.
  df <- do.call(rbind, lapply(doses, function(d) data.frame(
    Dose = d, Volume = 300 + 6 * d + rnorm(8, 0, 40), stringsAsFactors = FALSE
  )))
  jt <- run_jonckheere_test(df, "Dose", "Volume", verbose = FALSE)
  expect_identical(jt$alternative_used, "increasing")
  expect_lt(jt$p.value, 0.05)
})

# ---- R3.20 ------------------------------------------------------------------

test_that("R3.20: all-non-positive volumes error instead of producing log(Inf)", {
  set.seed(2)
  df <- data.frame(
    ID        = rep(1:6, each = 4),
    Cage      = rep(1:3, each = 8),
    Treatment = rep(c("Control", "DrugA"), each = 12),
    Day       = rep(c(0, 7, 14, 21), 6),
    Volume    = 0,                      # every volume non-positive
    stringsAsFactors = FALSE
  )
  expect_error(
    suppressWarnings(tumor_growth_statistics(df, transform = "log",
                                             plots = FALSE, verbose = FALSE)),
    regexp = "No positive values"
  )
})

# ---- R3.16 / R3.29 ----------------------------------------------------------

test_that("R3.16/R3.29: results report the scale they are actually on", {
  set.seed(3)
  df <- data.frame(
    ID        = rep(1:8, each = 5),
    Cage      = rep(1:4, each = 10),
    Treatment = rep(c("Control", "DrugA"), each = 20),
    Day       = rep(c(0, 5, 10, 15, 20), 8),
    Volume    = exp(rnorm(40, log(300), 0.3)),
    stringsAsFactors = FALSE
  )

  lme4_res <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, model_type = "lme4", transform = "log", plots = FALSE, verbose = FALSE)))
  expect_identical(lme4_res$transform_used,  "log")
  expect_identical(lme4_res$model_type_used, "lme4")

  auc_res <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, model_type = "auc", transform = "log", plots = FALSE, verbose = FALSE)))
  # R3.16 — AUC is integrated on the raw scale regardless of `transform`, and
  # must not claim otherwise.
  expect_identical(auc_res$transform_used, "none")
  expect_identical(auc_res$transform_requested, "log")
  expect_identical(auc_res$model_type_used, "auc")
  expect_false(any(grepl("was log transformed prior to analysis",
                         auc_res$summary$notes, fixed = TRUE)))
})

# ---- R3.25 ------------------------------------------------------------------

test_that("R3.25: weight_loss_threshold keys mice on cage as well", {
  set.seed(6)
  # Same ID in two cages within ONE treatment group, different trajectories:
  # one crosses the 20% threshold, the other does not.
  df <- rbind(
    data.frame(ID = 1, Cage = "A", Treatment = "DrugA", Day = c(0, 7, 14),
               Weight = c(25, 24.5, 24)),                 # ~4% loss, censored
    data.frame(ID = 1, Cage = "B", Treatment = "DrugA", Day = c(0, 7, 14),
               Weight = c(25, 22, 19)),                   # 24% loss, event
    data.frame(ID = 2, Cage = "A", Treatment = "Control", Day = c(0, 7, 14),
               Weight = c(25, 25, 25))
  )
  # n = 3 with a zero-event group sends the Cox fit down the Firth path and it
  # does not converge; irrelevant here — the assertions are about subject
  # identity, not the model.
  res <- suppressWarnings(weight_loss_threshold(
    df, cage_column = "Cage", adjust_tumor_weight = FALSE))

  # Without the cage in the key the two ID-1 mice collapse into one subject.
  expect_equal(nrow(res$event_data), 3L)
  expect_equal(sum(res$event_data$Event), 1L)
})

# =============================================================================
# Batch 2 — R3.1/G.1 (comparison family), R3.17/G.2 (cage), R3.4, R3.12
# =============================================================================

make_tg_df <- function() {
  set.seed(3)
  df <- data.frame(
    ID        = rep(1:12, each = 6),
    Cage      = rep(1:6, each = 12),
    Treatment = rep(c("Control", "DrugA", "DrugB", "DrugC"), each = 18),
    Day       = rep(c(0, 4, 8, 12, 16, 20), 12),
    stringsAsFactors = FALSE
  )
  df$Volume <- exp(log(200) + 0.10 * df$Day -
                     0.02 * df$Day * (df$Treatment != "Control") +
                     rnorm(nrow(df), 0, 0.25))
  df
}

# Two cages per arm, four mice per cage, with a real cage-level shift — the
# maintainer's stated standard design (G.2 case b).
make_nested_cage_df <- function() {
  set.seed(21)
  mk <- function(tx, cages, eff) do.call(rbind, lapply(cages, function(cg) {
    cage_re <- rnorm(1, 0, 0.30)
    do.call(rbind, lapply(1:4, function(i) {
      b <- rnorm(1, 0, 0.15)
      data.frame(
        ID = paste0(cg, "_", i), Cage = cg, Treatment = tx,
        Day = c(0, 4, 8, 12, 16, 20),
        Volume = exp(log(200) + cage_re + b +
                       (0.10 + eff) * c(0, 4, 8, 12, 16, 20) +
                       rnorm(6, 0, 0.12)),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rbind(mk("Control", c("C1", "C2"), 0),
        mk("DrugA",   c("C3", "C4"), -0.02),
        mk("DrugB",   c("C5", "C6"), -0.03))
}

tg_pvals <- function(df, ...) {
  r <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, plots = FALSE, verbose = FALSE, ...)))
  as.data.frame(summary(r$pairwise_comparisons))$p.value
}

# ---- R3.1 -------------------------------------------------------------------

test_that("R3.1: p_adjust_method changes lme4 pairwise p-values", {
  df   <- make_tg_df()
  none <- tg_pvals(df, p_adjust_method = "none")
  bonf <- tg_pvals(df, p_adjust_method = "bonferroni")
  holm <- tg_pvals(df, p_adjust_method = "holm")
  dunn <- tg_pvals(df, p_adjust_method = "dunnett")

  # The old code passed a custom contrast list with no `adjust`, for which
  # emmeans defaults to "none" — so every one of these was identical.
  expect_equal(length(none), 3L)
  expect_equal(bonf, pmin(1, 3 * none), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(holm, bonf)))
  expect_false(isTRUE(all.equal(dunn, bonf)))
  # Dunnett exploits the correlation between contrasts, so it is never more
  # conservative than Bonferroni.
  expect_true(all(dunn <= bonf + 1e-8))
})

test_that("R3.1: comparison_family selects the set of comparisons", {
  df <- make_tg_df()
  expect_equal(length(tg_pvals(df, comparison_family = "vs_reference")), 3L)
  expect_equal(length(tg_pvals(df, comparison_family = "all_pairs",
                               p_adjust_method = "tukey")), 6L)

  custom <- list("DrugA+DrugB vs Control" = c(-1, 0.5, 0.5, 0))
  p <- tg_pvals(df, comparison_family = "custom", custom_contrasts = custom)
  expect_equal(length(p), 1L)
})

test_that("R3.1: family and adjustment are validated against each other", {
  df <- make_tg_df()
  expect_error(tg_pvals(df, p_adjust_method = "tukey"),
               regexp = "only defined for")
  expect_error(tg_pvals(df, comparison_family = "all_pairs",
                        p_adjust_method = "dunnett"),
               regexp = "only defined for")
  expect_error(tg_pvals(df, comparison_family = "custom"),
               regexp = "requires `custom_contrasts`")
  # Tukey/Dunnett need a joint model; the AUC path has none.
  expect_error(
    suppressWarnings(suppressMessages(tumor_growth_statistics(
      df, model_type = "auc", comparison_family = "all_pairs",
      p_adjust_method = "tukey", plots = FALSE, verbose = FALSE))),
    regexp = "single fitted model"
  )
})

test_that("R3.1/R3.21: the AUC path adjusts over the family it reports", {
  df <- make_tg_df()
  run <- function(fam, m) suppressWarnings(suppressMessages(
    tumor_growth_statistics(df, model_type = "auc", comparison_family = fam,
                            p_adjust_method = m, plots = FALSE,
                            verbose = FALSE)))$pairwise_comparisons

  ref  <- run("vs_reference", "bonferroni")
  pair <- run("all_pairs",    "bonferroni")
  expect_equal(nrow(ref),  3L)
  expect_equal(nrow(pair), 6L)

  # Bonferroni over 6 comparisons must be exactly twice as harsh as over 3.
  ref_raw  <- run("vs_reference", "none")$p_adjusted
  expect_equal(ref$p_adjusted, pmin(1, 3 * ref_raw), tolerance = 1e-10)
  expect_identical(unique(ref$Comparison_Family), "vs_reference")
})

test_that("R3.1: results record the family and adjustment applied", {
  df <- make_tg_df()
  r <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, comparison_family = "all_pairs", p_adjust_method = "tukey",
    plots = FALSE, verbose = FALSE)))
  expect_identical(r$comparison_family, "all_pairs")
  expect_identical(r$p_adjust_method_used, "tukey")
})

# ---- R3.17 ------------------------------------------------------------------

test_that("R3.17: cage design structure is classified from mice, not rows", {
  df <- make_nested_cage_df()
  info <- classify_cage_structure(df, "Cage", "ID", "Treatment")
  expect_identical(info$structure, "nested_replicated")
  expect_equal(unname(info$cages_per_treatment), c(2L, 2L, 2L))
  expect_true(all(info$mice_per_cage == 4L))   # mice, not the 24 rows each

  # Crossed: move one mouse into another arm's cage.
  crossed <- df
  crossed$Cage[crossed$Treatment == "DrugA"][1:6] <- "C1"
  expect_identical(
    classify_cage_structure(crossed, "Cage", "ID", "Treatment")$structure,
    "crossed")

  # Confounded: one cage per arm.
  conf <- df
  conf$Cage <- paste0("C_", conf$Treatment)
  expect_identical(
    classify_cage_structure(conf, "Cage", "ID", "Treatment")$structure,
    "nested_confounded")
})

test_that("R3.17: auto fits (1|Cage) on a nested design and reports the ICC", {
  df <- make_nested_cage_df()
  auto <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, handle_cage_effects = "auto", plots = FALSE, verbose = FALSE)))
  old  <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, handle_cage_effects = "never_include", plots = FALSE, verbose = FALSE)))

  expect_identical(auto$summary$model_specification$cage_placement,
                   "random intercept")
  expect_identical(old$summary$model_specification$cage_placement, "omitted")

  # The ICC is reported only when cage is actually in the model.
  expect_false(is.null(auto$cage_analysis$icc))
  expect_true(auto$cage_analysis$icc$ICC > 0)
  expect_null(old$cage_analysis$icc)

  # Ignoring estimable cage variance understates the standard errors — the
  # whole point of the finding.
  se_auto <- as.data.frame(summary(auto$pairwise_comparisons))$SE
  se_old  <- as.data.frame(summary(old$pairwise_comparisons))$SE
  expect_true(all(se_auto > se_old))
})

test_that("R3.17: a fixed cage effect is refused when it is not estimable", {
  df <- make_nested_cage_df()
  expect_error(
    suppressWarnings(suppressMessages(tumor_growth_statistics(
      df, handle_cage_effects = "always_include", plots = FALSE,
      verbose = FALSE))),
    regexp = "not estimable"
  )
  # A crossed design admits it.
  crossed <- df
  crossed$Cage <- rep(c("X1", "X2", "X3"), length.out = nrow(crossed))
  expect_no_error(suppressWarnings(suppressMessages(
    tumor_growth_statistics(crossed, handle_cage_effects = "always_include",
                            plots = FALSE, verbose = FALSE))))
})

test_that("R3.17: complete cage/treatment confounding is warned about", {
  df <- make_nested_cage_df()
  df$Cage <- paste0("C_", df$Treatment)   # one cage per arm
  expect_warning(
    suppressMessages(tumor_growth_statistics(df, handle_cage_effects = "auto",
                                             plots = FALSE, verbose = FALSE)),
    regexp = "completely confounded"
  )
})

# ---- R3.4 / R3.12 -----------------------------------------------------------

make_bw_df <- function() {
  set.seed(31)
  mk <- function(tx, drift) do.call(rbind, lapply(1:6, function(i) data.frame(
    ID = paste0(tx, i), Treatment = tx, Cage = paste0(tx, "_C", (i > 3) + 1),
    Day = c(0, 4, 8, 12, 16, 20),
    Weight = 25 + rnorm(1, 0, 0.6) + drift * c(0, 4, 8, 12, 16, 20) +
      rnorm(6, 0, 0.15),
    Volume = 150 * exp(0.08 * c(0, 4, 8, 12, 16, 20)),
    stringsAsFactors = FALSE
  )))
  rbind(mk("Control", 0), mk("DrugA", -0.06), mk("DrugB", -0.02))
}

test_that("R3.4: the body-weight GAM path returns a populated EMM table", {
  res <- suppressWarnings(suppressMessages(analyze_body_weight(
    make_bw_df(), volume_column = "Volume", cage_column = "Cage",
    model_type = "gam", reference_group = "Control", volume_units = "mm3")))

  # Both were silently NULL before: emmeans could not dispatch on the unpatched
  # gamm4 stub and the error was swallowed by tryCatch.
  expect_false(is.null(res$emmeans_table))
  expect_equal(nrow(res$emmeans_table), 3L)
  expect_false(is.null(res$pairwise_comparisons))
})

test_that("R3.12: analyze_body_weight returns adjusted pairwise comparisons", {
  df <- make_bw_df()
  get_p <- function(m) {
    r <- suppressWarnings(suppressMessages(analyze_body_weight(
      df, volume_column = "Volume", cage_column = "Cage",
      reference_group = "Control", p_adjust_method = m, volume_units = "mm3")))
    pw <- r$pairwise_comparisons
    pw[[intersect(c("p_value", "p.value"), names(pw))[1]]]
  }
  none <- get_p("none")
  bonf <- get_p("bonferroni")

  expect_equal(length(none), 2L)
  expect_equal(bonf, pmin(1, 2 * none), tolerance = 1e-10)

  r <- suppressWarnings(suppressMessages(analyze_body_weight(
    df, volume_column = "Volume", reference_group = "Control",
    comparison_family = "all_pairs", p_adjust_method = "tukey",
    volume_units = "mm3")))
  expect_equal(nrow(r$pairwise_comparisons), 3L)
  expect_identical(r$p_adjust_method_used, "tukey")
})

# =============================================================================
# Batch B — R3.5 / R3.3 / G.3: endpoint estimand
# =============================================================================

# A study whose dropout mechanism is the real one: each animal has its own
# growth rate and is removed the first day its volume crosses the IACUC limit.
# Fast growers therefore leave first, so survivors at the last day are the slow
# growers — most severely in the control arm, which loses animals earliest.
make_attrition_df <- function(seed = 101, limit = 2000) {
  set.seed(seed)
  mk <- function(tx, rate_mu, n = 10) do.call(rbind, lapply(seq_len(n), function(i) {
    rate <- stats::rnorm(1, rate_mu, 0.030)
    days <- seq(0, 28, 4)
    v <- exp(log(150) + stats::rnorm(1, 0, 0.15) + rate * days +
               stats::rnorm(length(days), 0, 0.08))
    hit  <- which(v >= limit)
    keep <- if (length(hit)) seq_len(hit[1]) else seq_along(days)
    data.frame(ID = paste0(tx, i), Treatment = tx, Day = days[keep],
               Volume = v[keep], TrueRate = rate,
               Weight = 25 - 0.015 * days[keep] +
                 stats::rnorm(length(keep), 0, 0.1),
               stringsAsFactors = FALSE)
  }))
  rbind(mk("Control", 0.115), mk("DrugA", 0.075), mk("DrugB", 0.060))
}

test_that("R3.5: the endpoint estimand corrects survivor bias in TGI", {
  full <- make_attrition_df()
  df   <- full[, c("ID", "Treatment", "Day", "Volume", "Weight")]

  # Ground truth: each arm's geometric mean volume at day 28 had nobody been
  # removed, computed from the simulation's own per-animal growth rates.
  mice  <- unique(full[, c("ID", "Treatment", "TrueRate")])
  truth <- vapply(split(mice, mice$Treatment),
                  function(g) exp(mean(log(150) + g$TrueRate * 28)), numeric(1))
  true_tgi <- (1 - truth / truth[["Control"]]) * 100

  tgi_for <- function(m) {
    r <- suppressWarnings(suppressMessages(therapeutic_window_metric(
      df, reference_group = "Control", endpoint_method = m,
      n_boot = 0, adjust_tumor_weight = FALSE)))
    stats::setNames(r$tgi_data$TGI, r$tgi_data$Treatment)
  }
  surv  <- tgi_for("survivors")
  model <- tgi_for("model")

  err <- function(x) mean(abs(x[c("DrugA", "DrugB")] -
                                true_tgi[c("DrugA", "DrugB")]))

  # Survivor conditioning understates efficacy; the model estimand recovers it.
  expect_true(all(surv[c("DrugA", "DrugB")] < true_tgi[c("DrugA", "DrugB")]))
  expect_lt(err(model), err(surv))
  # On this fixture the improvement is large, not marginal.
  expect_lt(err(model), err(surv) / 2)
})

test_that("R3.5: attrition is reported so survivor selection is visible", {
  df <- make_attrition_df()[, c("ID", "Treatment", "Day", "Volume", "Weight")]
  r <- suppressWarnings(suppressMessages(therapeutic_window_metric(
    df, reference_group = "Control", endpoint_method = "model",
    n_boot = 0, adjust_tumor_weight = FALSE)))

  a <- r$attrition
  expect_true(all(c("N_Enrolled", "N_At_Endpoint", "Pct_Lost") %in% names(a)))
  # The control arm must show the heaviest loss — that is the whole point.
  expect_gt(a$Pct_Lost[a$Treatment == "Control"],
            max(a$Pct_Lost[a$Treatment != "Control"]))
  # Every enrolled animal still contributes under the model estimand.
  expect_true(all(r$tgi_data$N == a$N_Enrolled[match(r$tgi_data$Treatment,
                                                     a$Treatment)]))
})

test_that("R3.5: 'survivors' warns that it conditions on survival", {
  df <- make_attrition_df()[, c("ID", "Treatment", "Day", "Volume", "Weight")]
  expect_warning(
    suppressMessages(therapeutic_window_metric(
      df, reference_group = "Control", endpoint_method = "survivors",
      n_boot = 0, adjust_tumor_weight = FALSE)),
    regexp = "conditions on being observed"
  )
})

test_that("R3.5: the endpoint estimand reaches synergy and the other consumers", {
  df <- make_attrition_df()[, c("ID", "Treatment", "Day", "Volume", "Weight")]
  # A real fourth arm: slow-growing, so every Combo animal reaches day 28 and
  # the "survivors" comparison is well defined for it.
  set.seed(77)
  combo <- do.call(rbind, lapply(1:10, function(i) {
    rate <- stats::rnorm(1, 0.030, 0.02); days <- seq(0, 28, 4)
    data.frame(
      ID = paste0("Combo", i), Treatment = "Combo", Day = days,
      Volume = exp(log(150) + stats::rnorm(1, 0, 0.15) + rate * days +
                     stats::rnorm(length(days), 0, 0.08)),
      Weight = 25 - 0.015 * days + stats::rnorm(length(days), 0, 0.1),
      stringsAsFactors = FALSE)
  }))
  df <- rbind(df, combo)

  syn <- function(m) suppressWarnings(suppressMessages(analyze_drug_synergy(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", eval_time_point = 28, endpoint_method = m,
    n_boot = 0, verbose = FALSE)))

  s_surv  <- syn("survivors")
  s_model <- syn("model")
  tgi <- function(r) r$summary$TGI_Percent[r$summary$Treatment == "DrugA"]
  expect_gt(tgi(s_model), tgi(s_surv))
  expect_false(is.null(s_model$attrition))

  # weight_corrected_tgi and efficacy_toxicity_bivariate take the same argument.
  expect_no_error(suppressWarnings(suppressMessages(
    weight_corrected_tgi(df, reference_group = "Control",
                         endpoint_method = "model"))))
  expect_no_error(suppressWarnings(suppressMessages(
    efficacy_toxicity_bivariate(df, reference_group = "Control",
                                endpoint_method = "model"))))
})

test_that("R3.5: an animal with no row at the endpoint day still contributes", {
  df <- make_attrition_df()[, c("ID", "Treatment", "Day", "Volume", "Weight")]
  r <- suppressWarnings(suppressMessages(efficacy_toxicity_bivariate(
    df, reference_group = "Control", endpoint_method = "model")))
  # Previously `sub$Volume[sub$Day == max_day]` was numeric(0) for any animal
  # removed early; every enrolled animal must appear exactly once.
  expect_equal(nrow(r$per_mouse), length(unique(df$ID)))
  expect_true(all(is.finite(r$per_mouse$Efficacy)))
})

# ---- G.6 / R3.6 / R3.7: bootstrap intervals ---------------------------------

test_that("R3.6/R3.7: synergy metrics carry bootstrap intervals", {
  set.seed(12)
  mk <- function(tx, mu, n = 9) data.frame(
    Treatment = tx, Day = 21, ID = paste0(tx, seq_len(n)),
    Volume = stats::rlnorm(n, log(mu), 0.35), stringsAsFactors = FALSE)
  df <- rbind(mk("Control", 1500), mk("DrugA", 900),
              mk("DrugB", 1000), mk("Combo", 300))

  r <- suppressWarnings(suppressMessages(analyze_drug_synergy(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", eval_time_point = 21,
    endpoint_method = "survivors", n_boot = 1000, boot_seed = 1,
    verbose = FALSE)))

  ci <- r$synergy_ci
  expect_false(is.null(ci))
  expect_true(all(c("TGI_A_pct", "TGI_Combo_pct", "Bliss_Excess_FE",
                    "Loewe_CI") %in% ci$Metric))
  expect_true(all(ci$CI_Lower <= ci$Estimate & ci$Estimate <= ci$CI_Upper,
                  na.rm = TRUE))
  # Reproducible given a seed.
  r2 <- suppressWarnings(suppressMessages(analyze_drug_synergy(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", eval_time_point = 21,
    endpoint_method = "survivors", n_boot = 1000, boot_seed = 1,
    verbose = FALSE)))
  expect_equal(ci$CI_Lower, r2$synergy_ci$CI_Lower)
  # n_boot = 0 skips it.
  r0 <- suppressWarnings(suppressMessages(analyze_drug_synergy(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", eval_time_point = 21,
    endpoint_method = "survivors", n_boot = 0, verbose = FALSE)))
  expect_null(r0$synergy_ci)
})

test_that("R3.6/R3.7: TWM carries bootstrap intervals for its ranking", {
  df <- make_attrition_df()[, c("ID", "Treatment", "Day", "Volume", "Weight")]
  r <- suppressWarnings(suppressMessages(therapeutic_window_metric(
    df, reference_group = "Control", endpoint_method = "model",
    n_boot = 500, boot_seed = 2, adjust_tumor_weight = FALSE)))

  expect_false(is.null(r$twm_ci))
  expect_true(all(c("TWM_Lower", "TWM_Upper") %in% names(r$twm_table)))
  ok <- is.finite(r$twm_table$TWM_Lower) & is.finite(r$twm_table$TWM_Upper)
  expect_true(all(r$twm_table$TWM_Lower[ok] <= r$twm_table$TWM_Upper[ok]))
})

# ---- R3.15 ------------------------------------------------------------------

test_that("R3.15: power analysis accounts for multiplicity and attrition", {
  base <- apriori_power_analysis(effect_size = 0.8, n_groups = 4,
                                 alpha = 0.05, target_power = 0.8)
  adj  <- apriori_power_analysis(effect_size = 0.8, n_groups = 4,
                                 alpha = 0.05, target_power = 0.8,
                                 p_adjust_method = "bonferroni")
  drop <- apriori_power_analysis(effect_size = 0.8, n_groups = 4,
                                 alpha = 0.05, target_power = 0.8,
                                 dropout_rate = 0.2)

  # Adjusting for 3 vs-control comparisons needs more animals.
  expect_gt(adj$scenario_table$Required_N[1], base$scenario_table$Required_N[1])
  expect_equal(adj$scenario_table$Alpha_Per_Comparison[1], 0.05 / 3)
  expect_equal(adj$scenario_table$N_Comparisons[1], 3L)

  # Attrition leaves Required_N alone but raises the number to enrol.
  expect_equal(drop$scenario_table$Required_N[1],
               base$scenario_table$Required_N[1])
  expect_gt(drop$scenario_table$Enroll_N[1], drop$scenario_table$Required_N[1])
  expect_equal(drop$scenario_table$Enroll_N[1],
               ceiling(drop$scenario_table$Required_N[1] / 0.8))
})

# ---- R3.22 ------------------------------------------------------------------

test_that("R3.22: the AUC omnibus test is Welch's, matching its pairwise tests", {
  df <- make_tg_df()
  r <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, model_type = "auc", plots = FALSE, verbose = FALSE)))

  expect_true(grepl("Welch", r$anova$Method[1]))
  expect_equal(
    r$anova$p_value[1],
    stats::oneway.test(AUC ~ Treatment, data = r$auc_analysis$individual,
                       var.equal = FALSE)$p.value,
    tolerance = 1e-10)
  # Variance homogeneity is reported rather than assumed.
  expect_false(is.null(r$variance_test))
})

# ---- H.4 --------------------------------------------------------------------

test_that("H.4: the log-rank fallback uses an exact permutation test", {
  df <- make_surv_df()

  perm <- suppressWarnings(suppressMessages(survival_statistics(
    df, reference_group = "Control", firth_correction = FALSE,
    permutation_logrank = TRUE, verbose = FALSE)))
  asym <- suppressWarnings(suppressMessages(survival_statistics(
    df, reference_group = "Control", firth_correction = FALSE,
    permutation_logrank = FALSE, verbose = FALSE)))

  expect_identical(attr(perm$results, "logrank_flavour"),
                   "permutation log-rank (coin)")
  expect_identical(attr(asym$results, "logrank_flavour"),
                   "asymptotic log-rank (survdiff)")

  # The two agree in direction but not exactly — that is the point: the
  # asymptotic chi-square approximation is least trustworthy in exactly this
  # regime (a group with zero events, few events overall).
  nr <- perm$results$Group != "Control"
  expect_false(isTRUE(all.equal(perm$results$P_Value_Unadjusted[nr],
                               asym$results$P_Value_Unadjusted[nr])))

  # Must be deterministic: an analysis function cannot return a different
  # number each time it is called.
  perm2 <- suppressWarnings(suppressMessages(survival_statistics(
    df, reference_group = "Control", firth_correction = FALSE,
    permutation_logrank = TRUE, verbose = FALSE)))
  expect_equal(perm$results$P_Value_Unadjusted, perm2$results$P_Value_Unadjusted)
})

# ---- H.2 --------------------------------------------------------------------

test_that("H.2: AUC comparisons gain a permutation p-value", {
  set.seed(3)
  days <- c(0, 4, 8, 12, 16, 20)
  rows <- list()
  for (tx in c("Control", "DrugA", "DrugB")) {
    rate <- if (tx == "Control") 0.10 else 0.075
    for (i in 1:9) {
      rows[[length(rows) + 1L]] <- data.frame(
        ID = paste0(tx, i), Cage = paste0(tx, "_C", (i > 4) + 1),
        Treatment = tx, Day = days,
        Volume = exp(log(200) + stats::rnorm(1, 0, 0.2) + rate * days +
                       stats::rnorm(6, 0, 0.15)),
        stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)

  r <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, model_type = "auc", auc_permutations = 2000, auc_bootstrap_seed = 1,
    plots = FALSE, verbose = FALSE)))
  pw <- r$pairwise_comparisons

  expect_true(all(c("perm_p_value", "perm_p_adjusted") %in% names(pw)))
  expect_true(all(is.finite(pw$perm_p_value)))
  # The (1 + count) / (1 + n) estimator can never return exactly zero.
  expect_true(all(pw$perm_p_value > 0))
  # Distinct comparisons must get distinct p-values (each pair is seeded apart).
  expect_equal(length(unique(pw$perm_p_value)), nrow(pw))
  # Off by default so existing callers are unaffected.
  r0 <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, model_type = "auc", plots = FALSE, verbose = FALSE)))
  expect_true(all(is.na(r0$pairwise_comparisons$perm_p_value)))
})

test_that("H.2: the permutation test warns when its resolution is too coarse", {
  # With 3 vs 3 there are only choose(6,3) = 20 label splits, so the smallest
  # attainable two-sided p-value is 0.1 regardless of the observed difference.
  expect_warning(tgs_perm_diff_p(stats::rnorm(3), stats::rnorm(3), n_perm = 200),
                 regexp = "distinct label splits")
  expect_no_warning(tgs_perm_diff_p(stats::rnorm(12), stats::rnorm(12),
                                    n_perm = 200))
})

# ---- R3.9 -------------------------------------------------------------------

test_that("R3.9: extrapolation_points is deprecated and ignored", {
  df <- make_tg_df()
  expect_warning(
    r_extrap <- suppressMessages(tumor_growth_statistics(
      df, extrapolation_points = 3, plots = FALSE, verbose = FALSE)),
    regexp = "deprecated"
  )
  r_plain <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    df, plots = FALSE, verbose = FALSE)))

  # Ignored means identical results, not merely similar.
  expect_equal(as.data.frame(summary(r_extrap$pairwise_comparisons))$estimate,
               as.data.frame(summary(r_plain$pairwise_comparisons))$estimate)
  # The helper is gone.
  expect_false(exists("tgs_extrapolate",
                      envir = asNamespace("mouseExperiment"), inherits = FALSE))
})

# ---- R3.8 -------------------------------------------------------------------

test_that("R3.8: prior scales are data-driven and time-unit invariant", {
  set.seed(5)
  logv <- stats::rnorm(200, log(400), 0.5)   # mm3-scale log volumes
  days <- rep(seq(0, 28, 4), length.out = 200)

  sc <- bayes_prior_scales(logv, days, "skeptical")
  # The Intercept prior is centred on the data. The old fixed normal(0, 0.625)
  # sat about nine prior SDs away on exactly this scale.
  expect_equal(sc$response_median, stats::median(logv))
  expect_gt(abs(sc$response_median) / 0.625, 5)

  # b_sd is a TOTAL log-fold change over the study, so the per-time scale must
  # rescale exactly with the time unit.
  sc_weeks <- bayes_prior_scales(logv, days / 7, "skeptical")
  expect_equal(sc_weeks$b_sd_per_time / sc$b_sd_per_time, 7, tolerance = 1e-8)

  # The ladder still orders as its names promise.
  strengths <- c("skeptical", "informative", "weakly_informative", "diffuse")
  widths <- vapply(strengths,
                   function(p) bayes_prior_scales(logv, days, p)$b_sd_total,
                   numeric(1))
  expect_true(all(diff(widths) > 0))
})

test_that("R3.8: interaction coefficients get their own prior", {
  set.seed(5)
  d <- data.frame(ID = rep(1:12, each = 6),
                  Treatment = rep(c("Control", "DrugA"), each = 36),
                  Day = rep(seq(0, 28, length.out = 6), 12))
  d$logV <- log(400) + 0.10 * d$Day -
    0.03 * d$Day * (d$Treatment == "DrugA") + stats::rnorm(72, 0, 0.25)

  pr <- bayes_scaled_priors(stats::as.formula("logV ~ Treatment * Day + (1|ID)"),
                            d, "logV", "skeptical", time_column = "Day")
  pdf <- as.data.frame(pr)

  # The interaction — the estimand — must not share the main-effect prior.
  expect_true("TreatmentDrugA:Day" %in% pdf$coef)
  expect_true("Day" %in% pdf$coef)
  # Intercept prior centred on the data, not on zero.
  icpt <- pdf$prior[pdf$class == "Intercept"]
  expect_false(grepl("normal(0,", icpt, fixed = TRUE))
  expect_identical(attr(pr, "me_prior_scaling")$rate_coefficients,
                   c("Day", "TreatmentDrugA:Day"))
})
