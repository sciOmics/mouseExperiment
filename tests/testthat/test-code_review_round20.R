# Regression tests for CODE_REVIEW.md Round 20 fixes, labelled by finding ID.
# Each asserts observable behaviour, not the presence of code: two earlier
# "fixes" in this package existed in the source and never executed.

# ---- R20.7 ------------------------------------------------------------------

test_that("R20.7: k >= 3 powers each treated-vs-control comparison", {
  n_for <- function(k, d) {
    apriori_power_analysis(effect_size = d, n_groups = k, alpha = 0.05,
                           target_power = 0.8)$scenario_table$Required_N
  }
  # power.t.test at alpha / (k - 1), the reference column in R20.7. The old
  # omnibus path returned 27, 8, 10 and 7.
  expect_equal(n_for(3, 0.5), 78)
  expect_equal(n_for(3, 1.0), 21)
  expect_equal(n_for(4, 0.8), 35)
  expect_equal(n_for(4, 1.0), 23)

  # At the old recommendation (n = 8 for k = 3, d = 1) each comparison has about
  # a third of the intended power; a 20,000-study simulation of Bonferroni
  # t-tests against control gave 0.338.
  fp <- apriori_power_analysis(effect_size = 1, n_groups = 3, alpha = 0.05,
                               n_per_group = 8, mode = "find_power")
  expect_equal(fp$scenario_table$Achieved_Power, 0.338, tolerance = 0.005)
  expect_equal(fp$scenario_table$Alpha_Per_Comparison, 0.025)
  expect_match(fp$method_note, "treated-vs-control comparison")
})

test_that("R20.7: the default correction changes nothing with two groups", {
  a <- apriori_power_analysis(effect_size = 0.8, n_groups = 2, alpha = 0.05,
                              target_power = 0.8)
  b <- apriori_power_analysis(effect_size = 0.8, n_groups = 2, alpha = 0.05,
                              target_power = 0.8, p_adjust_method = "none")
  expect_identical(a$scenario_table$Required_N, b$scenario_table$Required_N)
  expect_null(a$method_note)
})

test_that("R20.7: the SD sensitivity table uses the per-comparison alpha", {
  r <- apriori_power_analysis(delta = 0.5, pooled_sd = 0.5, n_groups = 4,
                              alpha = 0.05, target_power = 0.8)
  mid <- r$sensitivity_table[r$sensitivity_table$SD_Change == "0%", ]
  # d = 1 at the assumed SD, so the middle row must match the scenario table.
  expect_equal(mid$Required_N, r$scenario_table$Required_N[1])
  expect_equal(mid$Alpha_Per_Comparison, 0.05 / 3)
})

# ---- R20.8 ------------------------------------------------------------------

# Four arms with log-linear growth; `rates` are per-day growth rates.
r20_synergy_df <- function(rates, n = 8, days = seq(0, 28, by = 3.5), seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(names(rates), function(arm) {
    do.call(rbind, lapply(seq_len(n), function(i) {
      b0 <- stats::rnorm(1, log(100), 0.1)
      r  <- rates[[arm]] + stats::rnorm(1, 0, 0.005)
      data.frame(ID = paste(arm, i), Treatment = arm, Day = days,
                 Volume = exp(b0 + r * days + stats::rnorm(length(days), 0, 0.08)))
    }))
  }))
}

test_that("R20.8: an agent that accelerates growth yields no Bliss quantities", {
  # The R20.8 scenario: DrugA accelerates growth, DrugB is inert, the
  # combination equals control.
  df <- r20_synergy_df(c(Control = 0.10, DrugA = 0.13, DrugB = 0.10, Combo = 0.10))
  res <- suppressWarnings(analyze_drug_synergy(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE, n_boot = 200, boot_seed = 1))

  expect_false(res$evaluable)
  expect_true(is.na(res$bliss_independence$synergy))      # was TRUE
  expect_true(is.na(res$bliss_independence$difference))
  expect_true(is.na(res$bliss_independence$expected_effect))
  excess <- res$synergy_ci[res$synergy_ci$Metric == "Bliss_Excess_FE", ]
  expect_true(is.na(excess$CI_Lower) && is.na(excess$CI_Upper))  # was [0.80, 3.79]
  expect_match(res$overall_assessment, "Not evaluable")

  ot <- suppressWarnings(suppressMessages(analyze_drug_synergy_over_time(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE)))
  expect_false(any(ot$synergy_summary$Evaluable))
  expect_equal(nrow(ot$peak_bliss_synergy), 0L)            # was day 28, 173.8 %

  p <- suppressWarnings(plot_synergy_trend(ot))
  built <- suppressWarnings(ggplot2::ggplot_build(p))
  ribbon_rows <- vapply(seq_along(p$layers), function(i) {
    if (inherits(p$layers[[i]]$geom, "GeomRibbon")) nrow(built$data[[i]]) else 0L
  }, integer(1))
  expect_equal(sum(ribbon_rows), 0L)                       # was drawn on all 9 days
})

test_that("R20.8: inhibitory agents are unaffected", {
  df <- r20_synergy_df(c(Control = 0.10, DrugA = 0.07, DrugB = 0.075, Combo = 0.03))
  res <- suppressWarnings(analyze_drug_synergy(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE, n_boot = 200, boot_seed = 1))
  expect_true(res$evaluable)
  expect_true(is.logical(res$bliss_independence$synergy) &&
                !is.na(res$bliss_independence$synergy))
  excess <- res$synergy_ci[res$synergy_ci$Metric == "Bliss_Excess_FE", ]
  expect_true(is.finite(excess$CI_Lower) && is.finite(excess$CI_Upper))

  ot <- suppressWarnings(suppressMessages(analyze_drug_synergy_over_time(
    df, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE)))
  expect_equal(nrow(ot$peak_bliss_synergy), 1L)
  expect_true(ot$synergy_summary$Evaluable[
    ot$synergy_summary$Time_Point == ot$peak_bliss_synergy$Time_Point])
})

# ---- T1 / R20.3 / R20.15 / R20.43: animal identity and column names ----------

# Four arms x 2 cages x 4 animals, 7 days. Three ID schemes for the same animals:
#   UID     unique across the study;
#   RID_arm 1..8 within each arm (reused across arms);
#   RID_cg  1..4 within each cage (reused across cages and arms).
r20_t1_df <- function(seed = 11) {
  set.seed(seed)
  rates <- c(Control = 0.12, DrugA = 0.09, DrugB = 0.095, Combo = 0.05)
  doses <- c(Control = 0, DrugA = 10, DrugB = 20, Combo = 30)
  days  <- c(0, 3, 7, 10, 14, 17, 21)
  out <- list(); uid <- 0L
  for (arm in names(rates)) for (cg in 1:2) for (k in 1:4) {
    uid <- uid + 1L
    b0  <- stats::rnorm(1, log(100), 0.15)
    r   <- rates[[arm]] + stats::rnorm(1, 0, 0.02)
    w0  <- stats::rnorm(1, 22, 1)
    out[[uid]] <- data.frame(
      UID = uid, RID_arm = (cg - 1L) * 4L + k, RID_cg = k,
      Treatment = arm, Cage = paste0(arm, "-C", cg), Day = days,
      Volume = exp(b0 + r * days + stats::rnorm(length(days), 0, 0.1)),
      Weight = w0 - 0.05 * days * (arm == "Combo") +
        stats::rnorm(length(days), 0, 0.2),
      Dose = doses[[arm]], stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# One row per animal: time to 1,000 mm3, censored at day 21.
r20_t1_surv <- function(d) {
  keys <- unique(d[, c("UID", "RID_arm", "RID_cg", "Treatment", "Cage")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    a <- d[d$UID == keys$UID[i], ]
    hit <- a$Day[a$Volume >= 1000]
    cbind(keys[i, ], Time = if (length(hit)) min(hit) else 21,
          Event = as.integer(length(hit) > 0))
  }))
}

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("T1: make_mouse_key() refuses a missing component", {
  expect_error(make_mouse_key(c("A", "B"), NULL), "NULL or empty")
  expect_error(make_mouse_key(c("A", "B"), character(0)), "NULL or empty")
  expect_error(make_mouse_key(c("A", "B"), c("1", "2", "3")), "different lengths")
  expect_equal(make_mouse_key(c("A", "B"), "1"), c("A|||1", "B|||1"))
  # A frame with no rows is not a missing column: no animals, no keys. The
  # weight-loss and body-weight functions reached this on an empty subset and
  # reported missing columns.
  expect_identical(make_mouse_key(character(0), character(0), character(0)),
                   character(0))
  expect_error(make_mouse_key(character(0), NULL), "NULL or empty")
})

test_that("R20.3: tumour-growth and body-weight models give the same answer with reused IDs", {
  d <- r20_t1_df()
  tg <- function(id) quiet(tumor_growth_statistics(
    d, id_column = id, cage_column = NULL, plots = FALSE, verbose = FALSE,
    reference_group = "Control"))
  a <- tg("UID"); b <- tg("RID_arm")
  expect_equal(lme4::fixef(b$model), lme4::fixef(a$model), tolerance = 1e-8)
  # Same partition of the data, so the fits agree up to optimiser noise (the
  # grouping levels are ordered differently): ~1e-9 on this fixture.
  expect_equal(as.matrix(stats::vcov(b$model)), as.matrix(stats::vcov(a$model)),
               tolerance = 1e-6)
  expect_equal(lme4::ngrps(b$model)[["RID_arm"]], 32L)   # was 8

  bw <- function(id) quiet(analyze_body_weight(
    d, weight_column = "Weight", id_column = id, volume_column = "Volume",
    adjust_tumor_weight = FALSE, reference_group = "Control"))
  x <- bw("UID"); y <- bw("RID_arm")
  expect_equal(y$fixed_effects, x$fixed_effects, tolerance = 1e-8)
  expect_equal(y$model_info$n_subjects, 32L)               # was 8
})

test_that("R20.3: synergy counts animals by treatment + ID + cage", {
  d <- r20_t1_df()
  syn <- function(id, cage) quiet(analyze_drug_synergy(
    d, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", id_column = id, cage_column = cage,
    n_boot = 0, verbose = FALSE))
  a <- syn("UID", NULL); b <- syn("RID_cg", "Cage")
  expect_equal(unname(b$group_n), rep(8L, 4L))            # was 4 per arm
  expect_equal(b$bliss_independence, a$bliss_independence, tolerance = 1e-8)
})

test_that("R20.15: over-time synergy takes an ID column that is not called ID", {
  d <- r20_t1_df()
  names(d)[names(d) == "UID"] <- "Animal"
  ot <- quiet(analyze_drug_synergy_over_time(
    d, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", id_column = "Animal", n_boot = 0, verbose = FALSE))
  expect_equal(nrow(ot$synergy_summary), 7L)               # was an error
  d$Animal <- NULL
  expect_error(quiet(analyze_drug_synergy_over_time(
    d, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE)), "Missing required columns")
})

test_that("R20.3 / R20.15: survival works with reused IDs and a differently named ID column", {
  s <- r20_t1_surv(r20_t1_df())
  sv <- function(id) quiet(survival_statistics(
    s, time_column = "Time", censor_column = "Event", treatment_column = "Treatment",
    id_column = id, cage_column = NULL, reference_group = "Control", verbose = FALSE))
  a <- sv("UID"); b <- sv("RID_arm")
  expect_equal(b$results$HR, a$results$HR, tolerance = 1e-8)
  expect_error(sv("NoSuchColumn"), "ID column 'NoSuchColumn' not found")
})

test_that("R20.3: dose-response groups animals by the key", {
  d <- r20_t1_df()
  dr <- function(id, cage = NULL) quiet(dose_response_statistics(
    d, dose_column = "Dose", treatment_column = "Treatment", volume_column = "Volume",
    day_column = "Day", id_column = id, cage_column = cage, verbose = FALSE))
  # IDs restart in every cage, so only treatment + ID + cage identifies an animal.
  a <- dr("UID"); b <- dr("RID_cg", "Cage")
  expect_equal(nrow(b$analysis_data), 32L)                  # 16 before: cage-mates merged
  expect_equal(stats::coef(b$linear_model), stats::coef(a$linear_model), tolerance = 1e-8)
})

test_that("R20.43: entry points accept column names with spaces", {
  d <- r20_t1_df()
  s <- r20_t1_surv(d)
  ren <- c(Day = "Study Day", Volume = "Tumor Volume", UID = "Animal ID",
           Treatment = "Treatment Group", Cage = "Cage No", Weight = "Body Weight",
           Dose = "Dose mg", Time = "Days On Study", Event = "Death Event")
  rn <- function(x) { names(x)[names(x) %in% names(ren)] <- ren[names(x)[names(x) %in% names(ren)]]; x }
  dd <- rn(d); ss <- rn(s)

  tg0 <- quiet(tumor_growth_statistics(d, id_column = "UID", cage_column = "Cage",
    plots = FALSE, verbose = FALSE, reference_group = "Control"))
  tg1 <- quiet(tumor_growth_statistics(dd, time_column = "Study Day",
    volume_column = "Tumor Volume", treatment_column = "Treatment Group",
    id_column = "Animal ID", cage_column = "Cage No", plots = FALSE,
    verbose = FALSE, reference_group = "Control"))
  expect_equal(unname(lme4::fixef(tg1$model)), unname(lme4::fixef(tg0$model)),
               tolerance = 1e-8)

  sv0 <- quiet(survival_statistics(s, time_column = "Time", censor_column = "Event",
    treatment_column = "Treatment", id_column = "UID", cage_column = NULL,
    reference_group = "Control", verbose = FALSE))
  sv1 <- quiet(survival_statistics(ss, time_column = "Days On Study",
    censor_column = "Death Event", treatment_column = "Treatment Group",
    id_column = "Animal ID", cage_column = NULL, reference_group = "Control",
    verbose = FALSE))                                       # "unexpected symbol" before
  expect_equal(sv1$results$HR, sv0$results$HR, tolerance = 1e-8)

  dr0 <- quiet(dose_response_statistics(d, id_column = "UID", verbose = FALSE))
  dr1 <- quiet(dose_response_statistics(dd, dose_column = "Dose mg",
    treatment_column = "Treatment Group", volume_column = "Tumor Volume",
    day_column = "Study Day", id_column = "Animal ID", verbose = FALSE))
  expect_equal(unname(stats::coef(dr1$linear_model)), unname(stats::coef(dr0$linear_model)),
               tolerance = 1e-8)                             # "unexpected symbol" before

  ot1 <- quiet(analyze_drug_synergy_over_time(dd, treatment_column = "Treatment Group",
    volume_column = "Tumor Volume", time_column = "Study Day", drug_a_name = "DrugA",
    drug_b_name = "DrugB", combo_name = "Combo", control_name = "Control",
    id_column = "Animal ID", n_boot = 0, verbose = FALSE))
  expect_equal(nrow(ot1$synergy_summary), 7L)
})

test_that("R20.3 / R20.43 / R20.62: Bayesian fits group by animal and accept any column names", {
  skip_on_cran()
  d <- r20_t1_df()
  d <- d[d$Day %in% c(0, 7, 14, 21), ]
  names(d)[names(d) == "RID_arm"] <- "Ear Tag"
  names(d)[names(d) == "Day"]     <- "Study Day"
  # R20.37: labels as the dashboard builds them with a dose column -- with
  # spaces and "+", which brms renames.
  labs <- c(Control = "Control", DrugA = "Drug A 10", DrugB = "Drug B 20",
            Combo = "Drug A + B")
  d$Treatment <- unname(labs[d$Treatment])
  bt <- quiet(bayesian_tumor_growth(
    d, time_column = "Study Day", volume_column = "Volume", id_column = "Ear Tag",
    reference_group = "Control", random_effects_specification = "slope",
    plots = FALSE, verbose = FALSE,
    mcmc = tg_mcmc(chains = 1, warmup = 150, iter = 150, seed = 1)))
  expect_equal(brms::ngrps(bt$model)$Animal, 32L)          # was 8 ear tags
  expect_setequal(unique(bt$growth_rates$ID), as.character(1:8))  # original IDs shown
  expect_equal(nrow(bt$growth_rates), 32L)
  # R20.37: each arm's animals get their arm's slope, not the control's. The
  # fixture's true rates are 0.12 (Control) and 0.05 (Drug A + B).
  gr <- tapply(bt$growth_rates$growth_rate, bt$growth_rates$Treatment, mean)
  expect_lt(gr[["Drug A + B"]], gr[["Control"]] - 0.04)    # equal before
  expect_equal(gr[["Drug A + B"]], 0.05, tolerance = 0.3)

  s <- r20_t1_surv(r20_t1_df())
  s$`.__DerivedEvent` <- factor(as.character(s$Event), levels = c("0", "1"))
  names(s)[names(s) == "Time"] <- ".__DerivedTime"
  bs <- quiet(bayesian_survival(
    s, time_column = ".__DerivedTime", event_column = ".__DerivedEvent",
    treatment_column = "Treatment", id_column = "RID_arm", reference_group = "Control",
    family = "weibull", plots = FALSE, verbose = FALSE,
    mcmc = tg_mcmc(chains = 1, warmup = 150, iter = 150, seed = 1)))
  expect_equal(bs$summary$data_description$subjects, 32L)  # was 8
  expect_equal(bs$summary$data_description$total_events, sum(s$Event))
})

# ---- R20.83: declared volume units -------------------------------------------

test_that("R20.83: a small-tumour mm3 study is not read as cm3", {
  # The dashboard's weight demo: 0.1-1,092.6 mm3, median 6.21 -- read as cm3
  # before, which made every net weight 1000x wrong.
  set.seed(3)
  v <- c(stats::runif(300, 0.1, 10), stats::runif(60, 200, 1092.6))
  expect_lt(stats::median(v), 20)
  expect_identical(detect_volume_units(v), "mm3")
  expect_identical(detect_volume_units(v / 1000), "cm3")
})

test_that("R20.83: mass adjustment needs declared units and refuses implausible masses", {
  d <- r20_t1_df()
  bw <- function(...) analyze_body_weight(d, weight_column = "Weight",
    id_column = "UID", volume_column = "Volume", reference_group = "Control", ...)
  expect_error(bw(), "volume_units is required")
  expect_error(quiet(bw(volume_units = "cm3")), "exceeds 50% of body weight")
  ok <- quiet(bw(volume_units = "mm3"))
  expect_true(all(is.finite(ok$fixed_effects$Estimate)))
  # Declaring mm3 for cm3 data is plausible by mass but flagged by the data.
  d_cm3 <- d; d_cm3$Volume <- d_cm3$Volume / 1000
  # Collect every warning, in case the fit adds unrelated ones.
  msgs <- character(0)
  withCallingHandlers(
    analyze_body_weight(d_cm3, weight_column = "Weight", id_column = "UID",
      volume_column = "Volume", reference_group = "Control", volume_units = "mm3"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("look like 'cm3'", msgs, fixed = TRUE)))
  expect_error(weight_loss_threshold(d, weight_column = "Weight", id_column = "UID",
    volume_column = "Volume"), "volume_units is required")
  expect_error(therapeutic_window_metric(d, weight_column = "Weight", id_column = "UID",
    volume_column = "Volume", reference_group = "Control"), "volume_units is required")
})

# ---- R20.6: survival by declared unit of randomisation ------------------------

# Two arms x `n_cages` cages x `per_cage` mice, one row per animal. A shared cage
# frailty makes cage-mates correlated.
r20_cage_surv <- function(n_cages = 2, per_cage = 5, hr = 1, seed = 1) {
  set.seed(seed)
  out <- list()
  for (arm in c("Control", "Drug")) for (cg in seq_len(n_cages)) {
    frail <- stats::rnorm(1, 0, 0.8)
    rate  <- 0.08 * exp(frail) * (if (arm == "Drug") hr else 1)
    t     <- stats::rexp(per_cage, rate)
    out[[length(out) + 1]] <- data.frame(
      ID = seq_len(per_cage), Treatment = arm, Cage = paste0(arm, "-", cg),
      Time = pmin(t, 30), Event = as.integer(t <= 30))
  }
  do.call(rbind, out)
}

test_that("R20.6: randomisation by animal uses ordinary standard errors, no cluster()", {
  s <- r20_cage_surv()
  r <- quiet(survival_statistics(s, time_column = "Time", censor_column = "Event",
    treatment_column = "Treatment", cage_column = "Cage", id_column = "ID",
    reference_group = "Control", verbose = FALSE))
  expect_identical(r$randomisation_unit, "mouse")
  expect_false(r$cage_cluster_used)
  expect_false(any(grepl("cluster", deparse(stats::formula(r$model)))))
  expect_true(is.finite(r$results$CI_Lower[r$results$Group == "Drug"]))
  # Housing by arm is reported, not hidden.
  expect_match(r$cage_caveat, "cage-mates as independent")
})

test_that("R20.6: randomisation by cage permutes whole cages and reports the floor", {
  s <- r20_cage_surv(n_cages = 2, hr = 0.2)
  expect_warning(
    r <- survival_statistics(s, time_column = "Time", censor_column = "Event",
      treatment_column = "Treatment", cage_column = "Cage", id_column = "ID",
      randomisation_unit = "cage", reference_group = "Control", verbose = FALSE),
    "cannot reach p <")
  cp <- r$cage_permutation
  expect_equal(cp$Assignments, 6)                    # choose(4, 2)
  expect_equal(cp$Min_Attainable_P, 1 / 3)
  drug <- r$results[r$results$Group == "Drug", ]
  expect_gte(drug$P_Value_Unadjusted, 1 / 3)         # however strong the effect
  expect_true(is.na(drug$CI_Lower) && is.na(drug$CI_Upper))
  expect_identical(drug$P_Method, "cage-level permutation log-rank")

  # Three cages per arm: floor 0.1.
  s3 <- r20_cage_surv(n_cages = 3, hr = 0.2)
  r3 <- quiet(survival_statistics(s3, time_column = "Time", censor_column = "Event",
    treatment_column = "Treatment", cage_column = "Cage", id_column = "ID",
    randomisation_unit = "cage", reference_group = "Control", verbose = FALSE))
  expect_equal(r3$cage_permutation$Min_Attainable_P, 0.1)
})

test_that("R20.6: cage randomisation is refused when cages hold several treatments", {
  s <- r20_cage_surv()
  s$Cage <- rep(c("C1", "C2"), length.out = nrow(s))   # every cage mixes arms
  expect_error(quiet(survival_statistics(s, time_column = "Time", censor_column = "Event",
    treatment_column = "Treatment", cage_column = "Cage", id_column = "ID",
    randomisation_unit = "cage", reference_group = "Control", verbose = FALSE)),
    "hold more than one treatment")
})

# ---- Warnings the dashboard now shows (dashboard R20.D34) ----------------------

test_that("R20.76 / R20-N: an ordinary tumour-growth fit raises no spurious warnings", {
  # Two warnings reached dashboard users on every lme4 run once its warnings
  # were surfaced: "disregarded additional arguments" (influence() was given
  # an obs argument it does not have) and "Chi-squared approximation may be
  # incorrect" (an unused cage-treatment chi-square on measurement rows).
  d <- r20_t1_df()
  msgs <- character(0)
  res <- withCallingHandlers(
    tumor_growth_statistics(d, time_column = "Day", volume_column = "Volume",
                            id_column = "UID", treatment_column = "Treatment",
                            cage_column = "Cage", reference_group = "Control",
                            model_type = "lme4", include_diagnostics = TRUE,
                            plots = FALSE, verbose = FALSE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("disregarded additional arguments", msgs, fixed = TRUE)))
  expect_false(any(grepl("Chi-squared approximation", msgs, fixed = TRUE)))
  expect_null(res$cage_analysis$collinearity_test)
  expect_false(is.null(res$diag_cooks_distance))    # influence still computed
})

# ---- Step 4 (v0.26.0): evaluable days, the endpoint model, AUC, random slopes --

# Gompertz growth, V(t) = V0 exp((a / b)(1 - exp(-b t))), with removal once a
# tumour passes `limit`: the shape the old log-linear endpoint model could not
# follow (R20.1).
r20_gompertz <- function(seed = 1, n = 10, limit = 2000,
                         a = c(Control = 0.41, A = 0.30, B = 0.32, Combo = 0.20),
                         b = 0.1, days = c(0, 3, 7, 10, 14, 17, 21, 24, 28)) {
  set.seed(seed)
  out <- list()
  for (arm in names(a)) for (i in seq_len(n)) {
    lv <- log(100) + stats::rnorm(1, 0, 0.2) +
      (a[[arm]] * exp(stats::rnorm(1, 0, 0.08)) / b) * (1 - exp(-b * days)) +
      stats::rnorm(length(days), 0, 0.1)
    stop_at <- which(exp(lv) > limit)[1]
    keep <- if (is.na(stop_at)) seq_along(days) else seq_len(stop_at)
    out[[length(out) + 1L]] <- data.frame(ID = i, Treatment = arm, Day = days[keep],
                                          Volume = exp(lv[keep]))
  }
  do.call(rbind, out)
}
r20_gompertz_truth <- function(arm, t,
                               a = c(Control = 0.41, A = 0.30, B = 0.32, Combo = 0.20),
                               b = 0.1) {
  100 * exp((a[[arm]] / b) * (1 - exp(-b * t)))
}

test_that("R20-K: evaluable_days() applies the 50 % and 3-animal rule", {
  # 8 animals; the last measurements fall on days 7 (x2), 14 (x3) and 21 (x3).
  d <- do.call(rbind, lapply(1:8, function(i) {
    last <- c(7, 7, 14, 14, 14, 21, 21, 21)[i]
    data.frame(ID = i, Treatment = "A", Day = seq(0, last, 7))
  }))
  ev <- evaluable_days(d)
  tab <- ev$table
  expect_equal(tab$N_On_Study, c(8, 8, 6, 3))      # days 0, 7, 14, 21
  expect_equal(tab$Evaluable, c(TRUE, TRUE, TRUE, FALSE))  # 3/8 < 50 %
  expect_equal(ev$last_day, 14)
  expect_equal(ev$excluded$Day, 21)
  expect_match(ev$excluded$Arms, "A 3/8")
  expect_match(ev$rule, "at least 50%")
  # An arm of two animals is never evaluable (fewer than 3).
  two <- d[d$ID %in% 1:2, ]
  expect_true(is.na(evaluable_days(two)$last_day))
})

test_that("R20.1: synergy defaults to the last evaluable day and refuses later days", {
  d <- r20_gompertz(1)
  ev <- evaluable_days(d)
  r <- quiet(analyze_drug_synergy(d, drug_a_name = "A", drug_b_name = "B",
    combo_name = "Combo", control_name = "Control", n_boot = 0, verbose = FALSE))
  expect_equal(r$eval_time_point, ev$last_day)
  expect_lt(r$eval_time_point, 28)                  # the last study day before
  expect_equal(r$evaluability$last_day, ev$last_day)
  expect_error(quiet(analyze_drug_synergy(d, drug_a_name = "A", drug_b_name = "B",
    combo_name = "Combo", control_name = "Control", eval_time_point = 28,
    n_boot = 0, verbose = FALSE)), "is not evaluable")
})

test_that("R20.1: with Gompertz growth the endpoint model recovers TGI and the Bliss excess", {
  # The straight-line model put the control at ~10x its true volume on the
  # last day and read TGIs of 81 / 75 / 94 against a truth of 28 / 22 / 60.
  est <- vapply(1:6, function(s) {
    d <- r20_gompertz(s)
    r <- quiet(analyze_drug_synergy(d, drug_a_name = "A", drug_b_name = "B",
      combo_name = "Combo", control_name = "Control", n_boot = 0, verbose = FALSE))
    t <- r$eval_time_point
    tru <- vapply(c("Control", "A", "B", "Combo"), r20_gompertz_truth, numeric(1), t = t)
    fe  <- 1 - tru[2:4] / tru[["Control"]]
    ctrl_est <- r$summary$Mean_Volume[1] / (1 - r$summary$TGI_Percent[1] / 100)
    c(ctrl_ratio = ctrl_est / tru[["Control"]],
      err_a = r$summary$TGI_Percent[1] - 100 * fe[[1]],
      err_combo = r$summary$TGI_Percent[3] - 100 * fe[[3]],
      err_bliss = r$bliss_independence$difference -
        (fe[[3]] - synergy_bliss_expected(fe[[1]], fe[[2]])))
  }, numeric(4))
  expect_true(all(abs(est["ctrl_ratio", ] - 1) < 0.25))
  expect_lt(abs(mean(est["err_a", ])), 6)
  expect_lt(abs(mean(est["err_combo", ])), 6)
  expect_lt(abs(mean(est["err_bliss", ])), 0.05)
})

test_that("R20.2: synergy intervals, tests and point estimates describe one estimand", {
  d <- r20_gompertz(2)
  r <- quiet(analyze_drug_synergy(d, drug_a_name = "A", drug_b_name = "B",
    combo_name = "Combo", control_name = "Control", n_boot = 500, boot_seed = 1,
    verbose = FALSE))
  ci <- r$synergy_ci
  # The point estimate is reported (not the median of the draws), and it lies
  # inside its own interval.
  expect_equal(ci$Estimate[ci$Metric == "TGI_A_pct"], r$summary$TGI_Percent[1])
  expect_equal(ci$Estimate[ci$Metric == "Bliss_Excess_FE"],
               r$bliss_independence$difference)
  expect_true(all(ci$CI_Lower <= ci$Estimate & ci$Estimate <= ci$CI_Upper))
  expect_match(r$interval_method, "endpoint model")
  expect_match(r$statistical_tests$Method[1], "endpoint model")
  # The per-animal estimand pairs with a per-animal bootstrap, and also reports
  # its point estimate.
  s <- quiet(analyze_drug_synergy(d, drug_a_name = "A", drug_b_name = "B",
    combo_name = "Combo", control_name = "Control", endpoint_method = "last_obs",
    n_boot = 300, boot_seed = 1, verbose = FALSE))
  expect_equal(s$synergy_ci$Estimate[1], s$summary$TGI_Percent[1])
  expect_match(s$interval_method, "bootstrap of animals")
})

test_that("R20.1 / R20.29: the endpoint model uses random slopes and leaves out pre-palpable zeros", {
  d <- r20_gompertz(3)
  # Two animals not yet palpable at day 0, one regressing to 0 at its last day.
  d$Volume[d$Treatment == "A" & d$ID %in% 1:2 & d$Day == 0] <- 0
  last_combo <- max(d$Day[d$Treatment == "Combo" & d$ID == 1])
  d$Volume[d$Treatment == "Combo" & d$ID == 1 & d$Day == last_combo] <- 0
  r <- quiet(analyze_drug_synergy(d, drug_a_name = "A", drug_b_name = "B",
    combo_name = "Combo", control_name = "Control", n_boot = 0, verbose = FALSE))
  info <- r$endpoint_model
  expect_true(info$random_effects %in% c("correlated", "uncorrelated"))
  expect_identical(info$time_basis, "natural spline, 3 df")
  expect_equal(info$n_prepalpable_excluded, 2L)
  expect_equal(info$n_zero_after_positive, 1L)
})

test_that("R20.1: over-time synergy analyses evaluable days only, and lists the others", {
  d <- r20_gompertz(4)
  ev <- evaluable_days(d)
  r <- quiet(analyze_drug_synergy_over_time(d, drug_a_name = "A", drug_b_name = "B",
    combo_name = "Combo", control_name = "Control", n_boot = 200, boot_seed = 1,
    verbose = FALSE))
  expect_equal(r$synergy_summary$Time_Point, ev$days)
  expect_equal(r$evaluability$excluded$Day, ev$excluded$Day)
  s <- r$synergy_summary
  expect_true(all(s$TGI_Combo_Lower <= s$TGI_Combo & s$TGI_Combo <= s$TGI_Combo_Upper))
})

test_that("R20.17: dose-response analyses the last evaluable day and reports TGI per dose", {
  set.seed(5)
  days <- c(0, 4, 8, 12, 16, 20)
  d <- do.call(rbind, lapply(c(0, 10, 30), function(dose) do.call(rbind, lapply(1:6, function(i) {
    lv <- log(100) + stats::rnorm(1, 0, 0.15) + (0.22 - 0.005 * dose) * days +
      stats::rnorm(length(days), 0, 0.08)
    stop_at <- which(exp(lv) > 2000)[1]
    keep <- if (is.na(stop_at)) seq_along(days) else seq_len(stop_at)
    data.frame(ID = i, Treatment = if (dose == 0) "Control" else "Drug", Dose = dose,
               Day = days[keep], Volume = exp(lv[keep]))
  }))))
  ev <- evaluable_days(transform(d, Treatment = as.character(Dose)))
  r <- quiet(dose_response_statistics(d, verbose = FALSE))
  expect_equal(r$endpoint_day, ev$last_day)
  expect_lt(r$endpoint_day, 20)       # controls are removed before day 20
  expect_true(all(r$analysis_data$Day == r$endpoint_day))
  tt <- r$tgi_table
  expect_equal(tt$Dose, c(0, 10, 30))
  expect_true(tt$Control[tt$Dose == 0])
  expect_true(all(diff(tt$TGI) > 0))
  expect_true(all(tt$TGI_Lower <= tt$TGI & tt$TGI <= tt$TGI_Upper))
  # A day on which the control has thinned out is refused, not silently
  # analysed without it.
  expect_error(quiet(dose_response_statistics(d, time_point = 20, verbose = FALSE)),
               "is not evaluable")
})

test_that("R20.4: the AUC compares fitted curves over a common window", {
  # The review's case: control 0.15/day, drug 0.08/day, removal at 2,000.
  # Per-animal trapezoids integrate each animal over its own follow-up, so the
  # early-removed controls get small AUCs.
  set.seed(8)
  days <- seq(0, 30, 3)
  d <- do.call(rbind, lapply(names(c(Control = 0.15, Drug = 0.08)), function(arm) {
    rate <- c(Control = 0.15, Drug = 0.08)[[arm]]
    do.call(rbind, lapply(1:10, function(i) {
      lv <- log(150) + stats::rnorm(1, 0, 0.2) + (rate + stats::rnorm(1, 0, 0.01)) * days +
        stats::rnorm(length(days), 0, 0.08)
      stop_at <- which(exp(lv) > 2000)[1]
      keep <- if (is.na(stop_at)) seq_along(days) else seq_len(stop_at)
      data.frame(ID = i, Treatment = arm, Day = days[keep], Volume = exp(lv[keep]))
    }))
  }))
  r <- quiet(tumor_growth_statistics(d, cage_column = NULL, model_type = "auc",
    reference_group = "Control", p_adjust_method = "holm", plots = FALSE))
  # The descriptive trapezoids reverse the effect here...
  ind <- r$auc_analysis$individual
  expect_gt(mean(ind$AUC[ind$Treatment == "Drug"]),
            mean(ind$AUC[ind$Treatment == "Control"]))
  # ...the model-based AUC does not.
  pw <- r$posthoc$pairwise
  expect_lt(pw$ratio, 1)
  expect_lt(pw$ratio_upper, 1)
  expect_lt(pw$p_adjusted, 0.05)
  expect_equal(unname(r$auc_window["end"]), evaluable_days(d)$last_day)
  te <- r$treatment_effects
  expect_true(all(te$Lower_CL <= te$AUC & te$AUC <= te$Upper_CL))
})

test_that("R20.5: tumour growth uses random slopes and F-tests by default", {
  d <- r20_t1_df()
  r <- quiet(tumor_growth_statistics(d, id_column = "UID", cage_column = "Cage",
    reference_group = "Control", plots = FALSE, include_diagnostics = FALSE))
  expect_identical(r$random_effects$requested, "slope")
  expect_true(r$random_effects$used %in% c("slope", "slope_uncorrelated"))
  expect_true(all(c("F value", "DenDF") %in% names(r$anova)))
  expect_match(r$anova_method, "Satterthwaite")
  # Intercept-only stays available on request.
  r0 <- quiet(tumor_growth_statistics(d, id_column = "UID", cage_column = "Cage",
    reference_group = "Control", random_effects_specification = "intercept_only",
    plots = FALSE, include_diagnostics = FALSE))
  expect_identical(r0$random_effects$used, "intercept_only")
})
