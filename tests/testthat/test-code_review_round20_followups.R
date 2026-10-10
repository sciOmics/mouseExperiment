# =============================================================================
# Round 20 follow-ups closed in v0.33.0 (CODE_REVIEW.md R20-W): the endpoint
# model's zero floor, and correlated-slope fits started from centred time.
# Each test fails on v0.32.0.
# =============================================================================

m <- master_synthetic_data
quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# ---- Centred time: the fit does not depend on the time origin ------------------

test_that("R20-W: the default (correlated-slope) fit is the same counted from any origin", {
  fit <- function(d) quiet(tumor_growth_statistics(d, reference_group = "Vehicle",
                                                   plots = FALSE, include_diagnostics = FALSE))
  a <- fit(m)
  m14 <- m; m14$Day <- m14$Day + 14                 # from implantation, not randomisation
  b <- fit(m14)
  expect_identical(a$random_effects$used, "slope")
  # v0.32.0 stopped short on the shifted data: 365.93 against 366.45, and the
  # Treatment row's F was 112.7 against 92.1.
  expect_equal(as.numeric(logLik(b$model)), as.numeric(logLik(a$model)), tolerance = 1e-6)
  expect_equal(b$anova[["F value"]], a$anova[["F value"]], tolerance = 1e-3)
  slopes <- function(r) { f <- lme4::fixef(r$model); f[grepl("Day:|:Day", names(f))] }
  expect_equal(unname(slopes(b)), unname(slopes(a)), tolerance = 1e-5)
})

test_that("R20-W: the starting theta maps a centred fit back to the original origin", {
  d <- data.frame(ID = rep(1:12, each = 6), Day = rep(seq(10, 35, by = 5), 12))
  set.seed(4)
  b0 <- rnorm(12, 0, 0.3); b1 <- rnorm(12, 0, 0.04)
  d$y <- 4 + b0[d$ID] + (0.1 + b1[d$ID]) * d$Day + rnorm(nrow(d), 0, 0.1)
  d$ID <- factor(d$ID)
  f <- y ~ Day + (Day | ID)
  ctrl <- lme4::lmerControl()
  th <- me_centred_start(f, d, "Day", "ID", shift = mean(d$Day), control = ctrl)
  plain <- lme4::lmer(f, data = d, control = ctrl)
  started <- lme4::lmer(f, data = d, control = ctrl, start = list(theta = th))
  expect_equal(as.numeric(logLik(started)), as.numeric(logLik(plain)), tolerance = 1e-6)
  expect_equal(th, unname(lme4::getME(started, "theta")), tolerance = 1e-3)
})

# ---- The endpoint model's fixed zero floor ---------------------------------------

endpoint_df <- function(d) {
  data.frame(MouseKey = paste(d$Treatment, d$ID, d$Cage), Treatment = d$Treatment,
             Day = d$Day, Volume = d$Volume, stringsAsFactors = FALSE)
}
regressed <- function() {
  d <- m
  d$Volume[d$Treatment == "Drug_A High" & d$Day >= 20] <- 0     # complete regressions
  d
}

test_that("R20-W: a regression to zero takes the fixed floor, not the smallest volume", {
  d <- regressed()
  small <- which(d$Treatment == "Vehicle" & d$Day == 0)[1]
  d$Volume[small] <- 0.5            # one small reading in another arm
  em <- quiet(me_endpoint_model(endpoint_df(d)))
  y  <- stats::model.response(stats::model.frame(em$fit))
  n_zero <- sum(d$Treatment == "Drug_A High" & d$Day >= 20 & is.finite(d$Volume))
  # v0.32.0 set them to 0.5, the smallest positive volume in the study.
  expect_equal(sum(abs(y - log(2)) < 1e-12), n_zero)
  expect_identical(em$detection_limit, 2)
  expect_true(any(abs(y - log(0.5)) < 1e-12))    # the measured 0.5 is kept as measured
  em_cm3 <- quiet(me_endpoint_model(endpoint_df(transform(d, Volume = Volume / 1000)),
                                    zero_floor = 0.002, volume_units = "cm3"))
  expect_identical(em_cm3$detection_limit, 0.002)
})

test_that("R20-W: synergy, dose-response and the therapeutic window take the floor's units", {
  d <- regressed()
  syn <- function(...) quiet(analyze_drug_synergy(d, drug_a_name = "Drug_A High",
    drug_b_name = "Drug_B", combo_name = "Drug_A Mid + Drug_B", control_name = "Vehicle",
    cage_column = "Cage", n_boot = 0, verbose = FALSE, ...))
  a <- syn(volume_units = "mm3")
  b <- syn(zero_floor = 20)
  tgi <- function(r) r$summary$TGI_Percent[r$summary$Treatment == "Drug_A High"]
  expect_false(isTRUE(all.equal(tgi(a), tgi(b))))  # the floor reaches the estimate
  expect_true("volume_units" %in% names(formals(analyze_drug_synergy_over_time)))
  expect_true(all(c("volume_units", "zero_floor") %in% names(formals(dose_response_statistics))))
  expect_true("zero_floor" %in% names(formals(therapeutic_window_metric)))
})

# ---- Dashboard R20.D60: the Bliss difference in percentage points ---------------

test_that("R20.D60: the synergy metrics report the Bliss difference in percentage points", {
  r <- quiet(analyze_drug_synergy(m, drug_a_name = "Drug_A Mid", drug_b_name = "Drug_B",
    combo_name = "Drug_A Mid + Drug_B", control_name = "Vehicle", cage_column = "Cage",
    n_boot = 0, verbose = FALSE))
  sm <- r$synergy_metrics
  expect_identical(sm$Metric[1], "Bliss Difference (percentage points)")
  # v0.32.0 put the raw fraction here ("0.0179841522745299").
  expect_equal(as.numeric(sm$Value[1]), round(100 * r$bliss_independence$difference, 1))
})
