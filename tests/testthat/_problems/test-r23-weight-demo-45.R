# Extracted from test-r23-weight-demo.R:45

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "mouseExperiment", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
skip_if_not_installed("shiny", minimum_version = "1.5.0")
.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

# test -------------------------------------------------------------------------
h <- paste(as.character(g("toxicityUI")("m")), collapse = "")
expect_match(h, "Whether an arm differs from the reference is tested in the next table", fixed = TRUE)
