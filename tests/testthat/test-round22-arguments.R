# R22.7: the exported functions name the same thing the same way (0.40.0).
# There were no users yet, so the old names were removed, not aliased.

exported_formals <- function() {
  ns <- asNamespace("mouseExperiment")
  ex <- getNamespaceExports("mouseExperiment")
  ex <- ex[vapply(ex, function(f) is.function(get(f, envir = ns)), NA)]
  stats::setNames(lapply(ex, function(f) names(formals(get(f, envir = ns)))), ex)
}

test_that("R22.7: no exported function takes one of the retired argument names", {
  retired <- c("ID_column", "day_column", "survival_column", "censor_column", "control_name",
               "control_group_name", "boot_seed", "auc_bootstrap_seed", "auc_bootstrap_n",
               "auc_permutations", "eval_time_point", "time_point", "min_time_point",
               "max_time_point", "start_date", "custom_title", "custom_colors")
  fm <- exported_formals()
  hits <- Filter(length, lapply(fm, intersect, retired))
  expect_identical(hits, stats::setNames(list(), character()))
})

test_that("R22.7: the shared roles use one name each", {
  fm <- exported_formals()
  takes <- function(nm) sort(names(Filter(function(a) nm %in% a, fm)))
  expect_true(all(c("analyze_drug_synergy", "analyze_drug_synergy_over_time",
                    "dose_response_statistics", "survival_statistics",
                    "tumor_growth_statistics", "therapeutic_window_metric") %in% takes("reference_group")))
  expect_true(all(c("plot_tumor_growth", "survival_statistics", "bayesian_survival") %in%
                  takes("event_column")))
  expect_true(all(c("analyze_drug_synergy", "dose_response_statistics", "therapeutic_window_metric") %in%
                  takes("endpoint_day")))
  expect_true(all(c("plot_drug_synergy", "plot_synergy_trend", "plot_auc", "plot_bliss") %in% takes("title")))
  expect_true(all(c("analyze_drug_synergy", "therapeutic_window_metric", "tumor_growth_statistics") %in%
                  takes("seed")))
  expect_true("reference_date" %in% fm$calculate_dates)
})

test_that("R22.7: the dashboard's calls use the new names, so the R script does too", {
  # The calls are do.call() argument lists, which R's parser cannot tie to
  # the function, so they are checked as text.
  src <- paste(vapply(c("mod_dose_response.R", "mod_synergy.R", "mod_survival.R", "mod_toxicity.R"),
                      function(f) dash_source("R", f), ""), collapse = "\n")
  for (old in c("control_group_name =", "control_name =", "boot_seed =", "eval_time_point =",
                "censor_column =", "day_column ="))
    expect_false(grepl(old, src, fixed = TRUE), info = old)
})
