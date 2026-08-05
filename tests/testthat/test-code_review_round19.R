# =============================================================================
# Round 19 — the five exported functions with no test of any kind
#
# Found by scanning exports against the test corpus: 37 of 42 exports were
# referenced somewhere, five were not. That list is not academic --
# `bayesian_power_analysis()` contained a live `brms::update` call that could
# never resolve (R18.1), and nothing would have caught it or the next one.
#
# The three MCMC entry points are exercised at minimal depth: the point is that
# they run end to end and return the documented shape, not that they sample well.
# `tg_mcmc()` / `tg_priors()` get real behavioural tests, because their whole job
# is to override other arguments and a silent failure to do so is the K.4 bug
# class -- the user's settings would simply be ignored.
# =============================================================================

# ---- tg_mcmc() / tg_priors(): configuration objects -------------------------

test_that("R19.1: tg_mcmc() returns a tagged list with the documented defaults", {
  m <- tg_mcmc()
  expect_s3_class(m, "tg_mcmc")
  expect_true(all(c("chains", "warmup", "iter", "seed", "backend") %in% names(m)))
  expect_identical(m$chains, 4L)
  expect_identical(m$backend, "rstan")
})

test_that("R19.1: tg_priors() returns a tagged list with the documented defaults", {
  p <- tg_priors()
  expect_s3_class(p, "tg_priors")
  expect_true(all(c("strength", "b", "intercept", "sd", "sigma") %in% names(p)))
  expect_identical(p$strength, "skeptical")
})

test_that("R19.1: tg_mcmc() actually overrides the individual arguments", {
  # The entire reason the object exists. If the override silently failed, a user
  # who set chains/iter/seed through it would get the defaults instead and have
  # no way to tell -- the sampler would still run and still return results.
  none <- mouseExperiment:::.resolve_mcmc(NULL, 4L, 1000L, 500L, 42L, "rstan")
  expect_identical(none$chains, 4L)
  expect_identical(none$iter, 500L)

  over <- mouseExperiment:::.resolve_mcmc(
    tg_mcmc(chains = 2L, warmup = 250L, iter = 100L, seed = 7L, backend = "rstan"),
    4L, 1000L, 500L, 42L, "rstan")
  expect_identical(over$chains, 2L)
  expect_identical(over$warmup, 250L)
  expect_identical(over$iter,   100L)
  expect_identical(over$seed,   7L)
})

test_that("R19.1: tg_priors() actually overrides prior_strength", {
  none <- mouseExperiment:::.resolve_priors(NULL, "skeptical", NULL, NULL, NULL, NULL)
  expect_identical(none$strength, "skeptical")

  over <- mouseExperiment:::.resolve_priors(
    tg_priors(strength = "diffuse", sd = 1.5),
    "skeptical", NULL, NULL, NULL, NULL)
  expect_identical(over$strength, "diffuse")
  expect_equal(over$sd, 1.5)
})

test_that("R19.1: the resolvers reject a wrongly-typed config", {
  # A bare list looks close enough to be passed by mistake; taking it silently
  # would apply some settings and drop others.
  expect_error(mouseExperiment:::.resolve_mcmc(list(chains = 2L), 4L, 1000L, 500L,
                                               42L, "rstan"),
               "tg_mcmc")
  expect_error(mouseExperiment:::.resolve_priors(list(strength = "diffuse"),
                                                 "skeptical", NULL, NULL, NULL, NULL),
               "tg_priors")
})

# ---- bayesian_power_analysis() ----------------------------------------------

test_that("R19.2: bayesian_power_analysis returns the documented shape", {
  skip_if_not_installed("brms")
  # This is the function that carried the unresolvable brms::update call. Even a
  # shape test would have failed loudly on that, which is the point.
  r <- suppressWarnings(suppressMessages(bayesian_power_analysis(
    n_per_group = 4L, n_groups = 2L, n_simulations = 5L,  # 5 is the floor
    timepoints = c(0, 7, 14), n_chains = 1L, n_warmup = 100L, n_iter = 100L,
    seed = 3L)))

  expect_type(r, "list")
  expect_true(all(c("power_table", "power_curve_data", "params") %in% names(r)))
  expect_s3_class(r$power_table, "data.frame")
  expect_gt(nrow(r$power_table), 0L)
  pw <- r$power_table[[grep("[Pp]ower", names(r$power_table))[1]]]
  expect_true(all(pw >= 0 & pw <= 1, na.rm = TRUE))
  expect_identical(r$params$n_groups, 2L)
})

# ---- bayesian_synergy_over_time() -------------------------------------------

syn_ot_df <- function(seed = 21) {
  set.seed(seed)
  rates <- c(Control = 0.13, DrugA = 0.10, DrugB = 0.105, Combo = 0.06)
  rows <- list(); k <- 0L
  for (tx in names(rates)) for (m in 1:5) {
    k <- k + 1L; d <- seq(0, 14, 7)
    rows[[k]] <- data.frame(
      ID = paste0(tx, m), Treatment = tx, Day = d,
      Volume = exp(log(150) + stats::rnorm(1, 0, 0.08) + rates[[tx]] * d +
                     stats::rnorm(length(d), 0, 0.05)),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

test_that("R19.3: bayesian_synergy_over_time returns per-day Bliss with intervals", {
  skip_if_not_installed("brms")
  r <- suppressWarnings(suppressMessages(bayesian_synergy_over_time(
    syn_ot_df(), volume_column = "Volume", time_column = "Day",
    treatment_column = "Treatment", id_column = "ID",
    drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control",
    n_chains = 1L, n_warmup = 150L, n_iter = 150L, seed = 4L,
    plots = FALSE, verbose = FALSE)))

  expect_true(all(c("synergy_by_day", "tgi_by_day", "peak_bliss_day") %in% names(r)))
  sbd <- r$synergy_by_day
  expect_s3_class(sbd, "data.frame")
  expect_true(all(c("Day", "Bliss_Median", "Bliss_Lower", "Bliss_Upper",
                    "P_Bliss_Synergy") %in% names(sbd)))
  # R17.2 removed the Loewe columns; they must not come back.
  expect_false(any(grepl("Loewe", names(sbd))))
  # Intervals must bracket the median, and the probability must be a probability.
  ok <- is.finite(sbd$Bliss_Median)
  expect_true(all(sbd$Bliss_Lower[ok] <= sbd$Bliss_Median[ok]))
  expect_true(all(sbd$Bliss_Upper[ok] >= sbd$Bliss_Median[ok]))
  expect_true(all(sbd$P_Bliss_Synergy[ok] >= 0 & sbd$P_Bliss_Synergy[ok] <= 1))
  expect_true(all(c("Day", "Bliss_Median") %in% names(as.list(r$peak_bliss_day))))
})

# ---- bayesian_twm_from_data() -----------------------------------------------

test_that("R19.4: bayesian_twm_from_data composes a TWM from two posteriors", {
  skip_if_not_installed("brms")
  set.seed(31)
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) for (m in 1:5) {
    k <- k + 1L; d <- seq(0, 14, 7)
    drop <- if (tx == "DrugA") 0.10 else 0.02
    rows[[k]] <- data.frame(
      ID = paste0(tx, m), Treatment = tx, Day = d,
      Volume = exp(log(150) + (if (tx == "Control") 0.12 else 0.07) * d +
                     stats::rnorm(length(d), 0, 0.05)),
      Weight = 25 * (1 - drop * d / 14) + stats::rnorm(length(d), 0, 0.05),
      stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, rows)

  r <- suppressWarnings(suppressMessages(bayesian_twm_from_data(
    tg_df = df, bw_df = df, time_column = "Day", treatment_column = "Treatment",
    id_column = "ID", reference_group = "Control", volume_column = "Volume",
    weight_column = "Weight",
    n_chains = 1L, n_warmup = 150L, n_iter = 150L, seed = 6L,
    plots = FALSE, verbose = FALSE)))

  expect_identical(r$model_type_used, "bayes_twm")
  expect_true(all(c("twm_table", "tgi_summary", "wl_summary") %in% names(r)))
  expect_s3_class(r$twm_table, "data.frame")
  expect_gt(nrow(r$twm_table), 0L)
  # The composite scale must be declared -- TWM mixes a log-volume posterior with
  # a grams posterior, so a consumer cannot infer it (B7.1).
  expect_identical(r$meta$inference, "bayesian")
  expect_identical(r$meta$interval_type, "credible")
  expect_match(r$meta$estimate_scale, "TWM")
})

# ---- the scan that found them ----------------------------------------------

test_that("R19.5: every exported function is referenced by at least one test", {
  # The check that produced this round. 37 of 42 exports were covered; the five
  # that were not are tested above. Keeping the scan means a new export cannot be
  # added without a test, which is how bayesian_power_analysis went four rounds
  # of auditing with an unresolvable call in it.
  tdir <- testthat::test_path(".")
  corpus <- paste(unlist(lapply(
    list.files(tdir, pattern = "[.]R$", full.names = TRUE),
    function(f) readLines(f, warn = FALSE))), collapse = "\n")

  ex <- getNamespaceExports("mouseExperiment")
  ex <- ex[!grepl("^(print|summary|plot|format)\\.", ex)]   # S3 methods
  untested <- ex[!vapply(ex, function(f) grepl(f, corpus, fixed = TRUE), logical(1))]

  if (length(untested)) {
    cat("\nExported but untested:\n  ", paste(untested, collapse = "\n  "), "\n")
  }
  expect_length(untested, 0L)
})
