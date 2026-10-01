# =============================================================================
# Round 20, step 6 (v0.28.0): Bayesian robustness
#
# The compiled-model cache (R20.78), plots that do not hold the fitted model
# (R20.80), LOO kept despite warnings (R20.35), the shape prior (R20.38),
# one row per animal for Bayesian survival (R20.39), the prior layer (R20.40),
# raw-scale prior widths (R20.41) and smaller items (R20.74). Each test fails
# on v0.27.0.
#
# The fits reuse one compiled model per structure, so after the first fit of
# a structure each one samples without compiling.
# =============================================================================

quiet6 <- function(expr) suppressWarnings(suppressMessages(expr))

tg6 <- function(d, ref = "Control", plots = FALSE, ...) {
  quiet6(bayesian_tumor_growth(d, reference_group = ref, n_chains = 2L,
                               n_iter = 300L, n_warmup = 300L, seed = 3L,
                               plots = plots, ...))
}

# TRUE when anything reachable from `x` (lists, attributes, closures and
# non-namespace environments) is a brmsfit.
holds_brmsfit <- function(x) {
  seen  <- new.env()
  found <- FALSE
  walk <- function(o, depth = 0L) {
    if (found || depth > 40L || is.null(o)) return(invisible())
    if (inherits(o, "brmsfit")) { found <<- TRUE; return(invisible()) }
    if (is.environment(o)) {
      if (isNamespace(o) || identical(o, globalenv()) || identical(o, baseenv()) ||
          identical(o, emptyenv())) return(invisible())
      key <- format(o)
      if (!is.null(seen[[key]])) return(invisible())
      assign(key, TRUE, envir = seen)
      for (nm in ls(o, all.names = TRUE)) {
        walk(tryCatch(get(nm, envir = o, inherits = FALSE), error = function(e) NULL),
             depth + 1L)
      }
      walk(parent.env(o), depth + 1L)
    } else if (is.function(o)) {
      walk(environment(o), depth + 1L)
    } else if (is.list(o)) {
      # By index: `for (el in o)` binds el to an empty symbol when a list holds
      # one (ggplot2 4's S7 objects do), and forcing it errors.
      for (i in seq_along(o)) walk(o[[i]], depth + 1L)
    }
    at <- attributes(o)
    for (i in seq_along(at)) walk(at[[i]], depth + 1L)
    invisible()
  }
  walk(x)
  found
}

test_that("R20.78: a fit of the same structure reuses the compiled model and draws the same", {
  d1 <- make_tg_simple()
  d2 <- d1
  d2$Treatment <- ifelse(d2$Treatment == "Control", "Vehicle", "Drug X+Y")  # brms renames these
  d2$Volume    <- d2$Volume * 1.2
  a <- tg6(d1)
  b <- tg6(d2, ref = "Vehicle")
  expect_true(b$model_reused)                     # compiled again before (about 20 s)
  expect_match(b$summary$methods$compiled_model, "reused")
  old <- options(mouseExperiment.cache_compiled_models = FALSE)
  on.exit(options(old), add = TRUE)
  fresh <- tg6(d2, ref = "Vehicle")
  expect_false(fresh$model_reused)
  # Same compiled code, data and seed: the same draws.
  expect_equal(b$posterior_summary$Estimate, fresh$posterior_summary$Estimate)
  expect_equal(b$posterior_summary$Parameter, fresh$posterior_summary$Parameter)
})

test_that("R20.80 / R20.40 / R20.74: plots do not hold the model; the prior layer is drawn; residuals survive a missing volume", {
  d <- make_tg_simple()
  d$Volume[3] <- NA
  r <- tg6(d, plots = TRUE)
  plots <- r[grepl("_plot$", names(r))]
  plots <- plots[!vapply(plots, is.null, logical(1L))]
  expect_true(all(c("pp_check_plot", "posterior_dist_plot", "prior_posterior_plot",
                    "credible_intervals_plot", "mcmc_trace_plot", "residuals_plot") %in%
                    names(plots)))                # residuals_plot was NULL with an NA volume
  for (nm in names(plots)) {
    expect_false(holds_brmsfit(plots[[nm]]), info = nm)   # every plot held it
  }
  # The prior layer (R20.40): brms names per-coefficient prior draws prior_b_<coef>.
  expect_true(all(c("Prior", "Posterior") %in% r$prior_posterior_plot$data$Source))
  # Small headline draws replace the model for trace and rank plots.
  expect_s3_class(r$posterior_draws, "draws_array")
  expect_true(all(c("b_Intercept", "sigma") %in% posterior::variables(r$posterior_draws)))
  expect_false(any(grepl("^r_", posterior::variables(r$posterior_draws))))
})

test_that("R20.35: LOO survives an influential observation and counts it", {
  d <- make_tg_simple()
  d$Volume[5] <- d$Volume[5] * 50                 # a planted outlier
  loo <- tg6(d)$loo_diagnostics
  expect_false(is.null(loo))                       # was NULL on the warning
  expect_gte(loo$n_high_k, 1L)                     # could only ever be 0
  expect_equal(loo$k_threshold, round(1 - 1 / log10(600), 3))   # 600 draws
  expect_true(nzchar(loo$warnings))
})

# Weibull survival with shape 6 and HR 0.088: 8 animals per arm, followed to
# day 30, so a few treated animals are censored.
weib6 <- function(seed = 5, n = 8, shape = 6, scale_c = 20, hr = 0.088, end = 30) {
  set.seed(seed)
  scale_t <- scale_c * hr^(-1 / shape)
  d <- data.frame(ID = seq_len(2 * n), Treatment = rep(c("Control", "Drug"), each = n),
                  Cage = 1)
  t <- c(scale_c * (-log(stats::runif(n)))^(1 / shape),
         scale_t * (-log(stats::runif(n)))^(1 / shape))
  d$Day   <- pmin(t, end)
  d$Event <- as.integer(t <= end)
  d
}

test_that("R20.38 / R20.74: the shape prior does not pull the hazard ratio toward 1", {
  d <- weib6()
  r <- quiet6(bayesian_survival(d, time_column = "Day", event_column = "Event",
                                reference_group = "Control", include_cage_effect = FALSE,
                                n_chains = 2L, n_iter = 1000L, n_warmup = 1000L,
                                seed = 11L))
  shape <- posterior::as_draws_df(r$model)$shape
  q <- stats::quantile(shape, c(0.025, 0.975), names = FALSE)
  expect_true(q[1] < 6 && 6 < q[2])               # was [1.66, 5.45]
  expect_lt(r$treatment_effects$HR[r$treatment_effects$Group == "Drug"], 0.15)  # was 0.242
  # R20.74: the metadata reports the priors used, with the shape prior.
  expect_identical(r$summary$methods$prior_aux, "lognormal(1, 1)")
  expect_false(grepl("^normal\\(0,", r$summary$methods$prior_intercept))
  # R20.74: coverage counts uncensored times only.
  expect_equal(r$ppc_coverage$n_obs, sum(d$Event))
})

test_that("R20.39: bayesian_survival() refuses longitudinal data", {
  d <- weib6()
  expect_error(bayesian_survival(rbind(d, d), time_column = "Day", event_column = "Event"),
               "requires one row per animal")
})

test_that("R20.41: on a raw volume scale the treatment prior width scales with the response", {
  y  <- c(100, 400, 900, 1600, 2500, 3600)
  tt <- c(0, 7, 14, 21, 28, 35)
  expect_equal(bayes_prior_scales(y, tt, "skeptical", log_scale = FALSE)$b_sd_total,
               0.25 * stats::mad(y))
  expect_equal(bayes_prior_scales(log(y), tt, "skeptical")$b_sd_total, 0.25)
})

test_that("R20.74: a manual prior names what is missing, and a backend from tg_mcmc() is checked", {
  d <- make_tg_simple()
  expect_error(bayesian_tumor_growth(d, prior_strength = "manual", prior_b = "normal(0, 1)"),
               "missing: prior_intercept, prior_sd, prior_sigma")
  skip_if(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr is installed")
  expect_error(bayesian_tumor_growth(d, mcmc = tg_mcmc(backend = "cmdstanr")),
               "requires the cmdstanr package")
})

# Last in the file: it empties the cache the fits above filled.
test_that("R20.78: clear_compiled_model_cache() empties the cache", {
  expect_gte(clear_compiled_model_cache(), 1L)
  expect_identical(clear_compiled_model_cache(), 0L)
})
