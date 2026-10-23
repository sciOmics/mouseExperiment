# Round 22's Major fixes (CODE_REVIEW.md Round 22). Each test fails on 0.38.0.

# ---- R22.1: a seeded result does not depend on the caller's RNG kind -------
# The dashboard's worker runs jobs with L'Ecuyer-CMRG (future's seed = TRUE);
# an R session uses Mersenne-Twister. set.seed(seed) kept the active kind, so
# the same seed drew different numbers in the two.

under_kind <- function(kind, expr) {
  old <- RNGkind(kind)
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  set.seed(1)
  force(expr)
}

test_that("R22.1: the LMM power simulation gives the same power under either RNG kind", {
  sim <- function() suppressWarnings(apriori_power_simulation(
    n_per_group = 6, n_groups = 2, n_simulations = 20, timepoints = seq(0, 21, 3), seed = 20260730L))
  a <- under_kind("Mersenne-Twister", sim())
  b <- under_kind("L'Ecuyer-CMRG", sim())
  expect_equal(a$power_table, b$power_table)
})

test_that("R22.1: the randomisation test draws the same permutations under either RNG kind", {
  perm <- function() trajectory_permutation_test(master_synthetic_data,
    spec = perm_spec(unit = "mouse", n_perm = 49, seed = 20260730L))
  a <- under_kind("Mersenne-Twister", perm())
  b <- under_kind("L'Ecuyer-CMRG", perm())
  expect_identical(a$perm_statistics, b$perm_statistics)
  expect_identical(a$p_value, b$p_value)
})

test_that("R22.1: the synergy bootstrap and the AUC draws do not depend on the RNG kind", {
  syn <- function() analyze_drug_synergy(master_synthetic_data, drug_a_name = "Drug_A Mid",
    drug_b_name = "Drug_B", combo_name = "Drug_A Mid + Drug_B", control_name = "Vehicle",
    n_boot = 100, boot_seed = 20260730L, verbose = FALSE)
  expect_equal(under_kind("Mersenne-Twister", syn())$synergy_ci,
               under_kind("L'Ecuyer-CMRG", syn())$synergy_ci)
  auc <- function() suppressWarnings(tumor_growth_statistics(master_synthetic_data,
    model_type = "auc", reference_group = "Vehicle", p_adjust_method = "holm", plots = FALSE))
  expect_equal(under_kind("Mersenne-Twister", auc())$treatment_effects,
               under_kind("L'Ecuyer-CMRG", auc())$treatment_effects)
})

test_that("R22.1: the caller's RNG kind and stream are left as they were", {
  under_kind("L'Ecuyer-CMRG", {
    before <- .Random.seed
    invisible(suppressWarnings(apriori_power_simulation(
      n_per_group = 6, n_groups = 2, n_simulations = 10, timepoints = seq(0, 9, 3), seed = 1L)))
    expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
    expect_identical(.Random.seed, before)
  })
})

# ---- R22.2: a Bayesian survival run's effects are downloaded and drawn --------

bayes_surv_result <- function() list(
  model_type_used = "bayes", error = FALSE, reference_group = "Vehicle",
  treatment_effects = data.frame(Group = c("Vehicle", "A", "B"), Time_Ratio = c(1, 1.8, 1.2),
                                 Lower_CrI = c(1, 1.2, 0.8), Upper_CrI = c(1, 2.6, 1.7),
                                 HR = c(1, 0.4, 0.8), stringsAsFactors = FALSE),
  survival_data = data.frame(Time = c(5, 8, 9), Event = c(1, 0, 1), Treatment = c("Vehicle", "A", "B")))

test_that("R22.2: Bayesian survival offers its time ratios in Download Results and draws them", {
  res <- bayes_surv_result()
  ef <- surv_effects(res)
  expect_identical(ef$measure, "Time Ratio")
  expect_identical(ef$table, res$treatment_effects)
  items <- surv_dl_items(res, km_plot = function() NULL,
                         forest_plot = function() surv_forest_plot_styled(res, list()))
  ids <- vapply(items, `[[`, "", "id")
  expect_true(all(c("results", "plot_forest") %in% ids))
  expect_identical(items[[match("results", ids)]]$value(), res$treatment_effects)
  p <- surv_forest_plot_styled(res, list())
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "Time Ratios vs Vehicle", fixed = TRUE)
  expect_match(p$labels$x, "Time Ratio", fixed = TRUE)
  # A Cox run is unchanged: hazard ratios.
  cox <- list(model_type_used = "cox", error = FALSE, reference_group = "Vehicle",
              results = data.frame(Group = c("A", "Vehicle"), HR = c(0.5, 1), CI_Lower = c(0.2, NA),
                                   CI_Upper = c(0.9, NA), P_Value = c(0.02, NA)))
  expect_match(surv_forest_plot_styled(cox, list())$labels$title, "Hazard Ratios", fixed = TRUE)
})

test_that("R22.2: the CSV button beside the table writes a Bayesian run's effects", {
  rv <- shiny::reactiveValues(surv_results = bayes_surv_result())
  shiny::testServer(function(id) shiny::moduleServer(id, function(input, output, session)
    add_surv_download_handlers(input, output, session, rv)), {
      out <- utils::read.csv(output$download_surv_cox)
      expect_identical(out$Time_Ratio, c(1, 1.8, 1.2))
    })
})

test_that("R22.2: a Bayesian survival script prints the treatment effects, not NULL", {
  res <- bayes_surv_result()
  res$call_record <- record_call("bayesian_survival",
    list(df = res$survival_data, time_column = "Time", event_column = "Event"))
  s <- surv_script_bundle(res)$lines
  expect_true(any(s == "print(result$treatment_effects)"))
  expect_false(any(s == "print(result$results)"))
})

# ---- R22.4: weight_loss_threshold()'s Firth fit converges --------------------

wl_study <- function(seed = 1) {
  set.seed(seed)
  arms <- c("Vehicle", "A", "B", "C", "D", "E"); loss <- c(0.012, 0.008, 0.0015, 0, 0.009, 0)
  do.call(rbind, lapply(seq_along(arms), function(a) do.call(rbind, lapply(1:8, function(i) {
    w0 <- stats::rnorm(1, 22, 1); r <- max(0, stats::rnorm(1, loss[a], 0.004)) * (loss[a] > 0)
    data.frame(ID = i, Treatment = arms[a], Day = seq(0, 28, 4),
               Weight = w0 * (1 - r * seq(0, 28, 4)) + stats::rnorm(8, 0, 0.15))
  }))))
}

test_that("R22.4: the weight-loss threshold's Firth p-values converge for event-free arms", {
  skip_if_not_installed("coxphf")
  d <- wl_study()
  w <- character()
  wl <- withCallingHandlers(weight_loss_threshold(d, weight_column = "Weight",
    treatment_column = "Treatment", id_column = "ID", cage_column = NULL,
    reference_group = "Vehicle", threshold = 0.2),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_identical(wl$cox_method, "coxphf")
  expect_false(any(grepl("profile likelihood", w, fixed = TRUE)))
  # 0.38.0 (coxphf's defaults): arm A's p-value was NA.
  expect_false(anyNA(wl$cox_summary$prob))
  m <- coxphf::coxphf(survival::Surv(Time, Event) ~ Treatment, data = wl$event_data,
                      maxit = 1000, maxstep = 0.05)
  expect_equal(unname(wl$cox_summary$prob), unname(m$prob), tolerance = 1e-8)
})
