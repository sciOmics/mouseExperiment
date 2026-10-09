# =============================================================================
# Round 20 statistical Majors closed in v0.30.0 (CODE_REVIEW.md R20-T).
# Each test fails on v0.29.0 for the reason its item describes.
# =============================================================================

m <- master_synthetic_data
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# ---- R20.46: one reference rule ---------------------------------------------

test_that("R20.46: the control is recognised by name, and inside longer names", {
  expect_identical(me_resolve_reference(m$Treatment, quiet = TRUE), "Vehicle")
  expect_identical(me_resolve_reference(c("Active", "Vehicle (10% DMSO)"), quiet = TRUE),
                   "Vehicle (10% DMSO)")
  expect_identical(me_resolve_reference(c("aPD1", "IgG2a isotype"), quiet = TRUE),
                   "IgG2a isotype")
  # A combination containing the vehicle is not the control.
  expect_identical(me_resolve_reference(c("Drug A + Vehicle", "Vehicle", "Drug A"),
                                        quiet = TRUE), "Vehicle")
  expect_message(me_resolve_reference(m$Treatment, fn = "f"),
                 "using \"Vehicle\" as the reference group")
})

test_that("R20.46: no recognisable control, or an unknown reference, stops", {
  expect_error(me_resolve_reference(c("A", "B")), "no control arm recognised")
  expect_error(me_resolve_reference(c("Vehicle", "A"), "vehicle"),
               "is not one of the treatment arms")
})

test_that("R20.46: tumour growth defaults to the control, not the first arm alphabetically", {
  # v0.29.0 sorted the arms and took "Drug_A High".
  tg <- quiet(tumor_growth_statistics(m, plots = FALSE, include_diagnostics = FALSE))
  expect_true(all(grepl("Vehicle", tg$pairwise_comparisons$contrast)))
})

test_that("R20.46: body weight stops on an unknown reference instead of using another arm", {
  # v0.29.0 ignored "vehicle" (wrong case) and compared against "Drug_A High".
  expect_error(quiet(analyze_body_weight(m, volume_column = "Volume", volume_units = "mm3",
                                         reference_group = "vehicle")),
               "is not one of the treatment arms")
})

test_that("R20.46: the toxicity functions return the reference they used", {
  tw <- quiet(therapeutic_window_metric(m, volume_units = "mm3", n_boot = 50, boot_seed = 1))
  expect_identical(tw$reference_group, "Vehicle")
  wl <- quiet(weight_loss_threshold(m, threshold = 0.1))
  expect_identical(wl$reference_group, "Vehicle")
})

# ---- R20.48: synergy trend legend -----------------------------------------------

test_that("R20.48: each synergy-trend line is labelled with its own name", {
  tr <- quiet(analyze_drug_synergy_over_time(m, drug_a_name = "Drug_A Mid",
    drug_b_name = "Drug_B", combo_name = "Drug_A Mid + Drug_B", control_name = "Vehicle",
    min_time_point = 2, n_boot = 50, boot_seed = 1, verbose = FALSE))
  sc <- ggplot2::ggplot_build(plot_synergy_trend(tr))$plot$scales$get_scales("colour")
  br <- sc$get_breaks()
  # v0.29.0 labelled the red Drug A line "Drug_B" and the grey Bliss line "Drug_A Mid".
  expect_identical(as.character(sc$get_labels(br)), as.character(br))
  expect_identical(unname(sc$map("Drug_A Mid")), "red")
})

# ---- R20.24 / R20.63 / R20.64: survival ---------------------------------------

with_warnings <- function(expr) {
  ws <- character()
  val <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = ws)
}
per_animal <- local({
  mm <- m[order(m$ID, m$Day), ]
  mm[!duplicated(mm$ID, fromLast = TRUE), ]
})

test_that("R20.24: Firth intervals are coxphf's profile limits and agree with its p-values", {
  s <- quiet(survival_statistics(per_animal, reference_group = "Vehicle", verbose = FALSE))
  expect_identical(s$method_used, "coxphf")
  r <- s$results[s$results$Group != "Vehicle", ]
  # v0.29.0 used confint() (Wald) beside coxphf's profile-likelihood p-value.
  expect_equal(unname(r$CI_Lower), unname(s$model$ci.lower), tolerance = 1e-10)
  expect_equal(unname(r$CI_Upper), unname(s$model$ci.upper), tolerance = 1e-10)
  expect_true(all((r$CI_Upper < 1 | r$CI_Lower > 1) == (r$P_Value_Unadjusted < 0.05)))
})

test_that("R20.24: a Firth failure falls through to the log-rank tests instead of crashing", {
  # An event at time 0 with a zero-event arm makes coxphf fail; v0.29.0 then
  # returned an empty table and crashed computing medians.
  d <- data.frame(ID = 1:12, Treatment = rep(c("Vehicle", "A", "B"), each = 4), Cage = 1,
                  Day = c(0, 5, 8, 10, 20, 20, 20, 20, 6, 9, 20, 20),
                  Survival_Censor = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0))
  out <- with_warnings(survival_statistics(d, reference_group = "Vehicle", verbose = FALSE))
  expect_true(any(grepl("using pairwise log-rank tests instead", out$warnings)))
  expect_identical(out$value$method_used, "logrank")
  r <- out$value$results
  expect_true(all(is.finite(r$P_Value_Unadjusted[r$Group != "Vehicle"])))
})

test_that("R20.24 / R20.63: the event column must be 0/1; a missing event is dropped with a warning", {
  bad12 <- per_animal; bad12$Survival_Censor <- bad12$Survival_Censor + 1
  expect_error(survival_statistics(bad12, reference_group = "Vehicle", verbose = FALSE),
               "must be 0 \\(censored\\) or 1")
  txt <- per_animal
  txt$Survival_Censor <- ifelse(txt$Survival_Censor == 1, "Dead", "Alive")
  expect_error(survival_statistics(txt, reference_group = "Vehicle", verbose = FALSE),
               "found \"Dead\"")
  na1 <- per_animal; na1$Survival_Censor[3] <- NA
  out <- with_warnings(survival_statistics(na1, reference_group = "Vehicle", verbose = FALSE))
  expect_true(any(grepl("were left out", out$warnings)))
  expect_equal(sum(out$value$results$Total), nrow(per_animal) - 1L)
})

test_that("R20.64: verbose = FALSE prints nothing, and data with no events stop", {
  expect_no_message(suppressWarnings(
    survival_statistics(per_animal, reference_group = "Vehicle", verbose = FALSE)))
  none <- per_animal; none$Survival_Censor <- 0
  expect_error(survival_statistics(none, reference_group = "Vehicle", verbose = FALSE),
               "no animal has an event")
})

# ---- R20.25: cage handling and large data -------------------------------------

has_cage_term <- function(tg) grepl("Cage", paste(deparse(stats::formula(tg$model)), collapse = ""))
tg_fit <- function(d, ...) quiet(tumor_growth_statistics(d, reference_group = "Vehicle",
                                plots = FALSE, include_diagnostics = FALSE, ...))

test_that("R20.25: an unused factor level does not drop the cage random effect", {
  d <- m
  d$Treatment <- factor(d$Treatment, levels = c(unique(m$Treatment), "Unused arm"))
  tg <- tg_fit(d)
  # v0.29.0 counted the unused level as an arm with no cages: "confounded".
  expect_identical(tg$cage_analysis$structure$structure, "nested_replicated")
  expect_true(has_cage_term(tg))
})

test_that("R20.25: one single-cage arm keeps (1|Cage) for the study, with a warning", {
  d <- m; d$Cage[d$Treatment == "Drug_B"] <- "C09"
  out <- with_warnings(tumor_growth_statistics(d, reference_group = "Vehicle",
                                               plots = FALSE, include_diagnostics = FALSE))
  expect_identical(out$value$cage_analysis$structure$structure, "nested_partial")
  expect_true(has_cage_term(out$value))
  expect_true(any(grepl("Drug_B is housed in a single cage", out$warnings)))
})

test_that("R20.25: a request that cannot be honoured is refused out loud", {
  d <- m; d$Cage <- paste0("K", as.integer(factor(d$Treatment)))   # every arm one cage
  out <- with_warnings(tumor_growth_statistics(d, reference_group = "Vehicle",
                         handle_cage_effects = "as_random_effect", plots = FALSE,
                         include_diagnostics = FALSE))
  expect_true(any(grepl("'as_random_effect' cannot be honoured", out$warnings)))
  expect_false(has_cage_term(out$value))
})

test_that("R20.25: the default path works above 3,000 rows (asymptotic emmeans columns)", {
  big <- do.call(rbind, lapply(1:6, function(k) {
    d <- m; d$ID <- d$ID + 1000 * k; d$Cage <- paste0(d$Cage, "_", k); d
  }))
  expect_gt(nrow(big), 3000)
  tg <- tg_fit(big)                       # v0.29.0: "differing number of rows: 6, 0"
  expect_true(all(is.finite(tg$treatment_effects$Lower_CL)))
  expect_true(all(is.finite(tg$treatment_effects$Upper_CL)))
})

# ---- R20.30: arguments that changed nothing -----------------------------------

test_that("R20.30: necrotic exclusion reaches AUC, growth rates and the permutation test", {
  fit <- function(...) quiet(tumor_growth_statistics(m, reference_group = "Vehicle",
                                                     plots = FALSE, ...))
  ex <- list(necrotic_column = "Necrotic", necrotic_handling = "exclude")
  a0 <- fit(model_type = "auc"); a1 <- do.call(fit, c(list(model_type = "auc"), ex))
  expect_false(isTRUE(all.equal(a0$treatment_effects, a1$treatment_effects)))
  g0 <- fit(); g1 <- do.call(fit, ex)
  expect_false(isTRUE(all.equal(g0$growth_rates, g1$growth_rates)))
  ps <- perm_spec(unit = "mouse", n_perm = 49, seed = 1)
  p0 <- fit(permutation_test = ps); p1 <- do.call(fit, c(list(permutation_test = ps), ex))
  expect_false(isTRUE(all.equal(p0$permutation_test$observed_statistic,
                                p1$permutation_test$observed_statistic)))
})

test_that("R20.30: dose_column analyses one arm per dose, and leaves dosed labels alone", {
  dd <- quiet(calculate_volume(calculate_dates(dose_levels_synthetic_data,
                                               start_date = "24-Mar", year = 2025)))
  t1 <- quiet(tumor_growth_statistics(dd, reference_group = "Control", plots = FALSE,
                                      dose_column = "Dose"))
  # v0.29.0 pooled the three HDACi doses into one arm.
  expect_setequal(unique(t1$growth_rates$Treatment),
                  c("Control", "HDACi 10", "HDACi 50", "HDACi 100"))
  # Labels that already include the dose (as the dashboard passes them) are kept.
  expect_identical(me_arms_with_dose(c("HDACi 10", "HDACi 50"), data.frame(D = c(10, 50)), "D"),
                   c("HDACi 10", "HDACi 50"))
  pa <- dd[order(dd$Cage, dd$ID, dd$Day), ]
  pa <- pa[!duplicated(paste(pa$Treatment, pa$Cage, pa$ID), fromLast = TRUE), ]
  s1 <- quiet(survival_statistics(pa, reference_group = "Control", dose_column = "Dose",
                                  verbose = FALSE))
  expect_setequal(s1$results$Group, c("Control", "HDACi 10", "HDACi 50", "HDACi 100"))
})

test_that("R20.30: calculate_volume replaces an existing Volume and accepts tibbles", {
  cv <- quiet(calculate_volume(m[1:3, ], formula = "sphere"))
  expect_equal(sum(names(cv) == "Volume"), 1L)       # v0.29.0 appended a second one
  expect_equal(cv$Volume, pi / 6 * pmin(m$Length[1:3], m$Width[1:3])^3)
  tb <- dplyr::as_tibble(combo_treatment_synthetic_data[1:5, ])
  expect_true(is.numeric(calculate_volume(tb)$Volume))   # v0.29.0 errored
})

test_that("R20.30: n_comparisons without an adjustment warns instead of being ignored", {
  expect_warning(apriori_power_analysis(effect_size = 1, n_groups = 4, n_comparisons = 5,
                                        p_adjust_method = "none"),
                 "has no effect with p_adjust_method")
})

# ---- R20.32: the zero floor ------------------------------------------------------

test_that("R20.32: zeros take a fixed floor of 2 mm3 (0.002 cm3), reported in the result", {
  expect_identical(me_zero_floor("mm3"), 2)
  expect_identical(me_zero_floor("cm3"), 0.002)
  expect_identical(me_zero_floor("mm3", 0.5), 0.5)
  v <- suppressWarnings(me_floor_zeros(c(0, 150, 300, 0), 2))
  expect_equal(as.numeric(v), c(2, 150, 300, 2))
  expect_identical(attr(v, "n_floored"), 2L)
  dd <- quiet(calculate_volume(calculate_dates(dose_levels_synthetic_data,
                                               start_date = "24-Mar", year = 2025)))
  t1 <- quiet(tumor_growth_statistics(dd, reference_group = "Control", plots = FALSE))
  expect_identical(t1$zero_floor$value, 2)
  expect_gt(t1$zero_floor$n_floored, 0L)
})

test_that("R20.32: one small reading elsewhere no longer moves the floor", {
  d <- m
  d$Volume[d$Treatment == "Drug_A High" & d$Day >= 20] <- 0     # complete regressions
  floor_seen <- function(x) {
    t <- quiet(tumor_growth_statistics(x, reference_group = "Vehicle", plots = FALSE,
                                       include_diagnostics = FALSE))
    exp(min(stats::model.response(stats::model.frame(t$model))))
  }
  expect_equal(floor_seen(d), 2)
  d2 <- d; d2$Volume[which(d2$Treatment == "Vehicle")[1]] <- 2.5
  # v0.29.0 used half the smallest positive volume: 1.25 after this one edit.
  expect_equal(floor_seen(d2), 2)
})

test_that("R20.32: volumes that look like cm3 declared as mm3 draw a warning", {
  cm <- m; cm$Volume <- cm$Volume / 1000
  expect_warning(quiet_msg <- suppressMessages(tumor_growth_statistics(cm,
                   reference_group = "Vehicle", plots = FALSE, include_diagnostics = FALSE)),
                 "Are the volumes in mm3")
})

# ---- R20.31: growth rate first; Treatment at the mean day ---------------------

test_that("R20.31: the growth-rate test leads, and each row says what it tests", {
  tg <- tg_fit(m)
  expect_identical(rownames(tg$anova)[1], "Treatment:Day")
  expect_true(grepl("^Treatment \\(at day", rownames(tg$anova)[2]))
  expect_identical(tg$anova$Tests[1], "Do growth rates differ between arms?")
  expect_match(tg$anova_method, "Treatment compares the arms at the mean study day")
})

test_that("R20.31: the Treatment row no longer depends on where time is counted from", {
  # Random intercepts only, so the fit itself is identical under a time shift.
  f <- function(d) tg_fit(d, random_effects_specification = "intercept_only")$anova
  a <- f(m); m14 <- m; m14$Day <- m14$Day + 14
  b <- f(m14)
  # v0.29.0 tested the arms at day 0: p = 0.86 here, and 7e-26 from day 14.
  expect_equal(a[["Pr(>F)"]], b[["Pr(>F)"]], tolerance = 1e-6)
  expect_equal(a[["F value"]], b[["F value"]], tolerance = 1e-6)
})
