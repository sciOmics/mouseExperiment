# The maintainer's manual review: explanations (CODE_REVIEW.md Round 23,
# R23-C). The Bayesian diagnostics say what they mean (R23.6), the influential
# observations name their animal (R23.6), and the mixed models' Info tabs
# explain fixed and random effects (R23.21).

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
ui_html <- function(f) paste(as.character(g(f)("m")), collapse = "\n")

test_that("R23.6: each diagnostic is explained under the panel, in Tumor Growth and Survival", {
  help <- paste(as.character(g("bayes_diagnostics_help")()), collapse = "")
  for (term in c("Divergent transitions", "Max-treedepth hits", "E-BFMI", "Bayes R", "PPC coverage",
                 "PSIS-LOO", "Pareto-k"))
    expect_match(help, term, fixed = TRUE)
  expect_match(help, "efficiency, not validity", fixed = TRUE)
  for (f in c("tumorGrowthUI", "survivalUI")) {
    h <- ui_html(f)
    diag <- substr(h, regexpr('data-value="Diagnostics"', h), nchar(h))
    expect_match(diag, "What these diagnostics mean", fixed = TRUE, info = f)
    expect_false(grepl("Investigate these mice individually", h, fixed = TRUE), info = f)
  }
})

test_that("R23.6: an influential observation names its animal and measurement", {
  obs <- data.frame(Treatment = c("A", "A", "B"), ID = c("1", "2", "1"), Cage = c("C1", "C1", "C2"),
                    Day = c(0, 7, 14), Volume = c(100, 250, 900))
  res <- list(loo_diagnostics = data.frame(k_threshold = 0.7, n_high_k = 2L,
                                           pareto_k = I(list(c(0.2, 0.9, 1.3))),
                                           observations = I(list(obs))))
  inf <- g("bayes_influential_obs")(res)
  expect_identical(names(inf), c("Treatment", "ID", "Cage", "Day", "Volume", "Obs_Index", "Pareto_k"))
  expect_identical(inf$ID, c("1", "2"))              # by Pareto-k, highest first
  expect_identical(inf$Treatment, c("B", "A"))
  expect_identical(inf$Obs_Index, c(3L, 2L))
  # A result from before 0.51.0 has no observations: the index only, as before.
  res$loo_diagnostics$observations <- NULL
  expect_identical(names(g("bayes_influential_obs")(res)), c("Obs_Index", "Pareto_k"))
  # Observations that do not match the Pareto-k vector are not used.
  res$loo_diagnostics$observations <- I(list(obs[1:2, ]))
  expect_identical(names(g("bayes_influential_obs")(res)), c("Obs_Index", "Pareto_k"))
})

test_that("R23.6: a Bayesian tumour-growth fit records each observation's animal, day and volume", {
  skip_if_not_installed("brms")
  skip_on_cran()
  fit <- suppressWarnings(suppressMessages(bayesian_tumor_growth(
    master_synthetic_data, reference_group = "Vehicle", cage_column = "Cage",
    n_chains = 2, n_warmup = 300, n_iter = 300, seed = 1, plots = FALSE, verbose = FALSE)))
  loo <- fit$loo_diagnostics
  obs <- loo$observations[[1]]
  expect_identical(names(obs), c("Treatment", "ID", "Cage", "Day", "Volume"))
  expect_identical(nrow(obs), length(loo$pareto_k[[1]]))
  # The volumes are the data's, back on the volume scale.
  key <- paste(obs$Treatment, obs$ID, obs$Cage, obs$Day)
  d <- master_synthetic_data[is.finite(master_synthetic_data$Volume) & master_synthetic_data$Volume > 0, ]
  dk <- paste(d$Treatment, d$ID, d$Cage, d$Day)
  expect_true(all(key %in% dk))
  m <- match(key, dk)
  expect_equal(obs$Volume, d$Volume[m], tolerance = 1e-8)
})

test_that("R23.21: the mixed models' Info tabs explain fixed and random effects", {
  help <- paste(as.character(g("effects_help")("animal")), collapse = "")
  for (term in c("Fixed effects", "Random effects", "Why both", "random slope"))
    expect_match(help, term, fixed = TRUE)
  for (f in c("tumorGrowthUI", "toxicityUI"))
    expect_match(ui_html(f), "Fixed and random effects", fixed = TRUE, info = f)
  # The Info tab's formulas show the default random slopes.
  expect_false(grepl("Day + (1 | ID)", ui_html("tumorGrowthUI"), fixed = TRUE))
})
