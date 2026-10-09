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
