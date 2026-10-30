# Compiled Stan models kept on disk, so they outlive the R process
# (CODE_REVIEW.md Round 23, R23.29, R23-K).

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
fake_fit <- function() structure(list(fit = NULL), class = "brmsfit")

test_that("R23.29: the disk cache is off unless a directory is given", {
  withr::local_options(mouseExperiment.model_cache_dir = NULL)
  withr::local_envvar(MOUSEEXPERIMENT_MODEL_CACHE = "")
  expect_null(g("me_model_cache_dir")())
  d <- withr::local_tempdir()
  withr::local_envvar(MOUSEEXPERIMENT_MODEL_CACHE = d)
  expect_identical(g("me_model_cache_dir")(), d)
  # The option wins over the environment.
  d2 <- file.path(withr::local_tempdir(), "made")
  withr::local_options(mouseExperiment.model_cache_dir = d2)
  expect_identical(g("me_model_cache_dir")(), d2)
  expect_true(dir.exists(d2))
})

test_that("R23.29: the key covers the package versions, R, the platform and the compiler flags", {
  env <- g("me_model_build_env")()
  expect_true(all(c("brms", "rstan", "StanHeaders", "Rcpp") %in% names(env$packages)))
  expect_identical(env$platform, R.version$platform)
  f <- withr::local_tempfile(lines = "CXX17FLAGS = -O1")
  withr::local_envvar(R_MAKEVARS_USER = f)
  a <- g("me_model_build_env")()
  writeLines("CXX17FLAGS = -O2", f)
  b <- g("me_model_build_env")()
  expect_false(identical(rlang::hash(a), rlang::hash(b)))
})

test_that("R23.29: models are written whole, read back, dropped when unreadable, and pruned", {
  d <- withr::local_tempdir()
  withr::local_options(mouseExperiment.model_cache_dir = d, mouseExperiment.model_cache_disk_size = 2L)
  cache <- g(".me_model_cache")
  withr::defer(clear_compiled_model_cache())
  expect_true(g("me_model_cache_write")("k1", fake_fit()))
  expect_identical(list.files(d), "k1.rds")              # no temporary file left
  rm(list = ls(cache, all.names = TRUE), envir = cache)
  e <- g("me_model_cache_read")("k1")
  expect_s3_class(e$template, "brmsfit")
  expect_true(e$from_disk)
  expect_true(exists("k1", envir = cache))
  expect_null(g("me_model_cache_read")("absent"))
  # An unreadable file is deleted, not trusted.
  writeLines("not an RDS", file.path(d, "bad.rds"))
  expect_null(g("me_model_cache_read")("bad"))
  expect_false(file.exists(file.path(d, "bad.rds")))
  # The least recently used are pruned to the size set.
  Sys.setFileTime(file.path(d, "k1.rds"), Sys.time() - 3600)
  g("me_model_cache_write")("k2", fake_fit())
  g("me_model_cache_write")("k3", fake_fit())
  expect_setequal(list.files(d), c("k2.rds", "k3.rds"))
  # clear_compiled_model_cache(disk = TRUE) empties both.
  assign("k9", list(template = fake_fit()), envir = cache)
  # k1 (read back above) and k9 in the session, k2 and k3 on disk.
  expect_identical(clear_compiled_model_cache(disk = TRUE), 4L)
  expect_length(list.files(d), 0L)
})

test_that("R23.29: a Bayesian fit reuses a model from the disk cache after the session's cache is gone", {
  skip_if_not_installed("brms")
  skip_on_cran()
  d <- withr::local_tempdir()
  withr::local_options(mouseExperiment.model_cache_dir = d)
  clear_compiled_model_cache()
  withr::defer(clear_compiled_model_cache())
  fit <- function() suppressWarnings(suppressMessages(bayesian_tumor_growth(
    master_synthetic_data, reference_group = "Vehicle", cage_column = "Cage",
    n_chains = 2, n_warmup = 300, n_iter = 300, seed = 1, plots = FALSE, verbose = FALSE)))
  a <- fit()
  expect_identical(a$summary$methods$compiled_model, "compiled for this fit")
  expect_length(list.files(d, pattern = "[.]rds$"), 1L)
  clear_compiled_model_cache()                            # as after a restart
  b <- fit()
  expect_identical(b$summary$methods$compiled_model, "reused from an earlier fit")
  expect_equal(b$posterior_summary, a$posterior_summary, tolerance = 0)
})

test_that("R23.29: the server image turns the cache on, and CI proves it across containers", {
  dock <- dash_source("docker", "Dockerfile")
  expect_match(dock, "ENV MOUSEEXPERIMENT_MODEL_CACHE=/var/cache/mouseexperiment/models", fixed = TRUE)
  ci <- dash_source(".github", "workflows", "image.yaml")
  expect_match(ci, "-e REQUIRE_REUSE=true", fixed = TRUE)
  expect_match(dash_source("docker", "shiny-gate-fit.R"), "the compiled model was not reused", fixed = TRUE)
})
