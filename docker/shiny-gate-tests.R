# Image gate: the package's own test suite, run against the candidate image
# (sciOmicsVPS CODE_REVIEW.md V3-Q). Since mouseExperiment 0.34.0 that suite
# holds the dashboard's tests too (V3-S).
#
# CI (.github/workflows/image.yaml) pipes this into `R --no-save` inside
# throwaway containers of the candidate, after docker/shiny-gate.R, as the
# image's own user and without a network. A failed or erroring test fails the
# job, which a pull request into production must pass. Until 0.49.2 the
# server ran it, in the staging environment's build (R23-A); production
# builds skip it, as before: their commits passed it in CI.
#
# The Dockerfile installs the package with --install-tests, so the suite is
# that of the exact commit in the image. It runs from a writable copy (tests
# write plots and snapshots into their own directory), with NOT_CRAN=true as
# devtools::test() sets it.
#
# Bayesian fits are skipped. brms::brm() is traced to signal a testthat skip,
# so a test that reaches a fit stops there and counts as skipped. The package
# starts every fit with brm() (later fits reuse its compiled model), and
# compiling them would add about six minutes; gate 1 already compiles and
# samples with both Stan backends.
#
# Environment (optional):
#   GATE_TESTS_SHARD     "k/n": run the k-th of n shares of the test files
#                        (default "1/1", all of them), balanced by file size.
#                        CI runs the shards side by side, one container
#                        each: a suite uses one core.
#   GATE_TESTS_BAYES=true runs the Bayesian fits too.

options(warn = 1)
t0 <- Sys.time()
step <- function(...) cat(sprintf("[tests %6.1fs] ", as.numeric(Sys.time() - t0, units = "secs")), ..., "\n", sep = "")

Sys.setenv(NOT_CRAN = "true")
setwd(tempdir())

if (!identical(Sys.getenv("GATE_TESTS_BAYES"), "true")) {
  suppressMessages(trace("brm", where = asNamespace("brms"), print = FALSE,
    tracer = quote(testthat::skip("Bayesian fit; not run in the deploy gate"))))
  step("Bayesian fits will be skipped (brms::brm traced)")
}

shard <- suppressWarnings(as.integer(strsplit(Sys.getenv("GATE_TESTS_SHARD", "1/1"), "/", fixed = TRUE)[[1]]))
if (length(shard) != 2 || anyNA(shard) || shard[1] < 1 || shard[1] > shard[2])
  stop("GATE_TESTS_SHARD must be k/n with 1 <= k <= n, got '", Sys.getenv("GATE_TESTS_SHARD"), "'")

run_suite <- function(pkg) {
  src <- system.file("tests", "testthat", package = pkg)
  if (!nzchar(src)) stop(pkg, " was installed without its tests (--install-tests)")
  dir <- file.path(tempdir(), pkg)
  dir.create(dir)
  stopifnot(file.copy(src, dir, recursive = TRUE))
  # The shards are balanced by file size, a rough stand-in for run time: the
  # largest file first, each to the shard with the least so far, ties by name.
  # Every container computes the same assignment. On 0.34.0 the two shards ran
  # in 123 s and 94 s on a Mac (by name alternately: 158 s and 59 s).
  files <- sort(list.files(file.path(dir, "testthat"), pattern = "^test[-_].*[.][rR]$"))
  size <- file.size(file.path(dir, "testthat", files))
  load <- numeric(shard[2])
  owner <- integer(length(files))
  for (i in order(-size, files)) {
    k <- which.min(load)
    load[k] <- load[k] + size[i]
    owner[i] <- k
  }
  mine <- files[owner == shard[1]]
  # This shard's files, as a test_dir() filter. The other files stay on disk:
  # a test may read them (R19.5 checks that every export is named in some test
  # file). test_dir() matches the filter against each file's name without
  # "test-" and ".R"; characters other than letters, digits and _ are escaped.
  stems <- gsub("([^A-Za-z0-9_])", "[\\1]", sub("[.][rR]$", "", sub("^test[-_]", "", mine)))
  step(pkg, " ", as.character(utils::packageVersion(pkg)), ": running test shard ",
       shard[1], "/", shard[2], " (", length(mine), " of ", length(files), " files)")
  if (!length(mine)) return(0L)
  res <- testthat::test_dir(file.path(dir, "testthat"), package = pkg,
                            filter = paste0("^(", paste(stems, collapse = "|"), ")$"),
                            load_package = "installed", reporter = "check",
                            stop_on_failure = FALSE, stop_on_warning = FALSE)
  df <- as.data.frame(res)
  bad <- df[df$failed > 0 | df$error, c("file", "test"), drop = FALSE]
  step(pkg, ": ", nrow(df), " tests, ", nrow(bad), " failed, ",
       sum(df$skipped), " skipped")
  if (nrow(bad)) print(bad, row.names = FALSE)
  nrow(bad)   # an integer
}

failed <- run_suite("mouseExperiment")
if (failed > 0) stop("failing tests: ", failed)

cat("TESTS_OK\n")
