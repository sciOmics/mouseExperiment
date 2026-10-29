# Posterior distributions with a highest-density interval and the prior
# (CODE_REVIEW.md Round 23, R23.12, R23-G).

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

fake_draws <- function() {
  set.seed(7)
  rbind(data.frame(Parameter = "Drug", Value = stats::rnorm(4000, -0.4, 0.1), Source = "Posterior"),
        data.frame(Parameter = "Drug", Value = stats::rnorm(4000, 0, 2), Source = "Prior"),
        data.frame(Parameter = "Drug:Day", Value = stats::rexp(4000, 20), Source = "Posterior"),
        data.frame(Parameter = "Drug:Day", Value = stats::rnorm(4000, 0, 1), Source = "Prior"))
}

test_that("R23.12: the HDI is the narrowest interval holding the probability", {
  h <- g("me_hdi")(stats::qnorm(stats::ppoints(20000)), 0.95)
  expect_equal(h, c(-1.96, 1.96), tolerance = 0.01)
  # A skewed posterior: the HDI starts at the mode, narrower than the
  # equal-tailed interval.
  x <- stats::qexp(stats::ppoints(20000))
  h <- g("me_hdi")(x, 0.9)
  expect_lt(h[1], 0.01)
  expect_lt(diff(h), diff(stats::quantile(x, c(0.05, 0.95), names = FALSE)))
  expect_equal(mean(x >= h[1] & x <= h[2]), 0.9, tolerance = 0.001)
  tab <- posterior_hdi(fake_draws(), 0.8)
  expect_identical(names(tab), c("Parameter", "Lower", "Upper", "Level"))
  expect_identical(tab$Parameter, c("Drug", "Drug:Day"))
  expect_error(posterior_hdi(fake_draws(), 95), "between 0 and 1")
  expect_error(posterior_hdi(fake_draws(), parameters = "Other"), "Unknown parameter")
  expect_error(posterior_hdi(list(model_type_used = "bayes")), "0.55.0")
})

test_that("R23.12: plot_posterior() shades the HDI, names its limits, and adds the prior on request", {
  d <- fake_draws()
  p <- plot_posterior(d, hdi = 0.9)
  expect_s3_class(p, "ggplot")
  lim <- posterior_hdi(d, 0.9)
  shaded <- p$layers[[1]]$data
  for (k in seq_len(nrow(lim))) {
    s <- shaded[shaded$Parameter == lim$Parameter[k], ]
    expect_true(all(s$x >= lim$Lower[k] & s$x <= lim$Upper[k]))
  }
  expect_true(all(grepl("90% HDI \\[", levels(shaded$Label))))
  expect_length(p$layers, 3L)                       # area, posterior, zero line
  pp <- plot_posterior(d, show_prior = TRUE)
  expect_length(pp$layers, 4L)
  prior <- pp$layers[[4]]$data
  post <- pp$layers[[2]]$data
  # The prior is drawn over the posterior's range only.
  for (par in unique(post$Parameter))
    expect_equal(range(prior$x[prior$Parameter == par]), range(post$x[post$Parameter == par]))
  expect_match(pp$labels$subtitle, "grey dashed: the prior", fixed = TRUE)
  expect_identical(levels(plot_posterior(d, hdi = 0.9, parameters = "Drug")$layers[[2]]$data$Label),
                   grep("^Drug\n", levels(shaded$Label), value = TRUE))
})

test_that("R23.12: a Bayesian fit keeps its prior and posterior draws", {
  skip_if_not_installed("brms")
  skip_on_cran()
  fit <- suppressWarnings(suppressMessages(bayesian_tumor_growth(
    master_synthetic_data, reference_group = "Vehicle", cage_column = "Cage",
    n_chains = 2, n_warmup = 300, n_iter = 300, seed = 1, plots = FALSE, verbose = FALSE)))
  d <- fit$prior_posterior_data
  expect_true(is.data.frame(d))
  expect_setequal(unique(d$Source), c("Posterior", "Prior"))
  # The arms as labelled, not brms's names (Drug_AMidPDrug_B:Day).
  expect_true("Drug_A Mid + Drug_B" %in% d$Parameter)
  expect_true(any(grepl(" Day$", unique(d$Parameter))))
  expect_false(any(grepl(":", d$Parameter, fixed = TRUE)))
  expect_s3_class(plot_posterior(fit, show_prior = TRUE), "ggplot")
})

test_that("R23.12: Tumor Growth and Survival draw it in Plots; Survival has no Prior vs Posterior tab", {
  tg <- paste(as.character(g("tumorGrowthUI")("m")), collapse = "")
  sv <- paste(as.character(g("survivalUI")("m")), collapse = "")
  expect_match(tg, '<option value="posterior">Posterior Distributions (Bayesian)</option>', fixed = TRUE)
  expect_match(sv, '<option value="posterior">Posterior Distributions (Bayesian)</option>', fixed = TRUE)
  expect_match(tg, 'id="m-tg_opt_post_hdi"', fixed = TRUE)
  expect_false(grepl('data-value="Prior vs Posterior"', sv, fixed = TRUE))
  # The R script redraws it on the script's result.
  rec <- g("record_call")("plot_posterior", list(x = NULL, show_prior = TRUE, hdi = 0.9),
                          data_arg = NULL, exprs = list(x = quote(result)))
  steps <- g("tg_plot_steps")(NULL, rec)
  expect_length(steps, 1L)
  code <- paste(g("script_call")(steps[[1]]$rec, steps[[1]]$name), collapse = "\n")
  expect_match(code, "p_posterior <- plot_posterior(", fixed = TRUE)
  expect_match(code, "x = result", fixed = TRUE)
  expect_match(code, "hdi = 0.9", fixed = TRUE)
})

test_that("R23.12: the Tumor Growth server draws and offers it after a Bayesian run", {
  rv <- shiny::reactiveValues(
    uploaded_data = master_synthetic_data, data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
  shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    rv$tg_results <- list(error = FALSE, model_type_used = "bayes", prior_posterior_data = fake_draws())
    session$setInputs(tg_plot_selector = "posterior", tg_opt_post_prior = TRUE, tg_opt_post_hdi = 80)
    p <- tg_plots$posterior()
    expect_s3_class(p, "ggplot")
    expect_true(all(grepl("80% HDI", levels(p$layers[[1]]$data$Label), fixed = TRUE)))
    rec <- tg_plots$posterior_record()
    expect_identical(rec$args$hdi, 0.8)
    expect_true(rec$args$show_prior)
  })
})
