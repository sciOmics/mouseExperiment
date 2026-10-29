# Round 22 proposal F8 (CODE_REVIEW.md R22-I, 0.47.0): result classes with
# print() and summary() methods.

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
q <- function(e) suppressWarnings(suppressMessages(e))
m <- master_synthetic_data

results <- function() {
  ev <- survival_events(m, volume_column = "Volume", volume_limit = 1500)
  list(
    me_tumor_growth = q(tumor_growth_statistics(m, reference_group = "Vehicle", plots = FALSE)),
    me_survival = q(survival_statistics(ev, time_column = "Time", event_column = "Event",
                                        reference_group = "Vehicle", verbose = FALSE)),
    me_body_weight = q(analyze_body_weight(m, reference_group = "Vehicle", volume_units = "mm3")),
    me_weight_loss = q(weight_loss_threshold(m, threshold = 0.1, reference_group = "Vehicle")),
    me_therapeutic_window = q(therapeutic_window_metric(m, reference_group = "Vehicle",
                                                        volume_units = "mm3", n_boot = 20, seed = 1)),
    me_synergy = q(analyze_drug_synergy(m, drug_a_name = "Drug_A Mid", drug_b_name = "Drug_B",
                                        combo_name = "Drug_A Mid + Drug_B", reference_group = "Vehicle",
                                        n_boot = 20, seed = 1, verbose = FALSE)),
    me_dose_response = q(dose_response_statistics(m, treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
                                                  reference_group = "Vehicle", verbose = FALSE)),
    me_power = apriori_power_analysis(effect_size = 1),
    me_survival_power = apriori_survival_power(0.5, 0.9),
    me_permutation = q(trajectory_permutation_test(m, spec = perm_spec(unit = "mouse", n_perm = 19, seed = 1))),
    me_growth_delay = growth_delay(m, volume = 1000, reference_group = "Vehicle", n_boot = 20, seed = 1),
    me_data_check = check_study_data(m))
}

test_that("F8: every analysis returns a classed list that prints what it found", {
  for (cls in names(r <- results())) {
    x <- r[[cls]]
    expect_s3_class(x, c(cls, "me_result"), exact = TRUE)
    expect_type(x, "list")                                  # still a list
    out <- capture.output(print(x))
    title <- g(".me_result_spec")[[cls]]$title(x)
    expect_identical(out[1], paste0("<mouseExperiment: ", title, ">"), info = cls)
    expect_true(any(grepl("names(x)", out, fixed = TRUE)), info = cls)
    expect_lt(length(out), 200)                             # not the wall of fields
    s <- summary(x)
    expect_s3_class(s, "summary.me_result")
    expect_true(all(vapply(s, is.data.frame, NA)), info = cls)
    expect_silent(capture.output(print(s)))
  }
})

test_that("F8: the tables print their first rows and their names; subsets are plain lists", {
  tg <- q(tumor_growth_statistics(m, reference_group = "Vehicle", plots = FALSE))
  out <- capture.output(print(tg, n = 2))
  expect_true(any(grepl("^Treatment:Day ", out)))         # an ANOVA's terms are its row names
  expect_true(any(grepl("... and 4 more rows", out, fixed = TRUE)))
  expect_identical(class(tg[c("anova", "summary")]), "list")
  expect_identical(class(unclass(tg)), "list")
  # Bayesian results, shaped as the fits return them.
  b <- g("me_result")(list(family_used = "weibull",
                           treatment_effects = data.frame(Group = c("Vehicle", "A"), Time_Ratio = c(1, 2))),
                      "me_bayes_survival")
  expect_match(capture.output(print(b))[1], "Bayesian AFT, weibull", fixed = TRUE)
  expect_error(g("me_result")(list(), "me_unknown"))
})
