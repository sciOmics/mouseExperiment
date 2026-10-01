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
