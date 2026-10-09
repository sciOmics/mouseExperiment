# =============================================================================
# Power module parameter plumbing — CODE_REVIEW.md D.4 / K.15
#
# D.4: mod_power.R passed `%||%` fallbacks left over from the percent-scale
# parameterisation the backend abandoned when its data-generating model moved to
# the log scale (§2.2). `control_growth_rate %||% 15` where the parameter is a
# log rate means exp(15 * t), which overflows to Inf by day 14 and kills the fit
# with an opaque error. The fallbacks only fire when an input reads NULL, so this
# was latent -- exactly the kind of thing a test should hold down rather than a
# reader be expected to notice.
#
# These tests read the fallbacks out of the source rather than driving the UI,
# because the failure mode is a literal in the call, not a behaviour. That makes
# them cheap and precise: they fail on the specific mistake and nothing else.
# =============================================================================

power_src <- function() {
  p <- testthat::test_path("..", "..", "R", "mod_power.R")
  if (!file.exists(p)) p <- "R/mod_power.R"
  skip_if_not(file.exists(p), "mod_power.R not locatable")
  paste(readLines(p, warn = FALSE), collapse = "\n")
}

test_that("D.4: LMM fallbacks are on the log scale the backend expects", {
  src <- power_src()

  # Each pair is (input id, the value the backend's own default uses).
  expected <- list(
    pw_lmm_ctrl_rate     = "0.15",
    pw_lmm_tx_effect     = "0.10",
    pw_lmm_rand_slope_sd = "0.05",
    pw_lmm_rand_int_sd   = "0.20",
    pw_lmm_resid_sd      = "0.10"
  )

  for (id in names(expected)) {
    m <- regmatches(src, regexpr(
      paste0("input\\$", id, "\\s*%\\|\\|%\\s*[0-9.]+"), src))
    expect_length(m, 1L)
    got <- sub(".*%\\|\\|%\\s*", "", m)
    expect_identical(got, expected[[id]],
                     info = paste0(id, ": fallback must match the backend default; ",
                                   "a percent-scale value here overflows exp()"))
  }
})

test_that("D.4: the fallbacks agree with the UI defaults next to them", {
  # The UI input and its fallback should not be able to drift apart silently.
  src <- power_src()
  for (id in c("pw_lmm_ctrl_rate", "pw_lmm_tx_effect", "pw_lmm_rand_slope_sd",
               "pw_lmm_rand_int_sd", "pw_lmm_resid_sd")) {
    # Labels contain parentheses ("(/day, log-scale)"), so the span between the
    # id and its value must be matched permissively but bounded.
    ui <- regmatches(src, regexpr(
      paste0('ns\\("', id, '"\\)[\\s\\S]{0,200}?value = [0-9.]+'),
      src, perl = TRUE))
    fb <- regmatches(src, regexpr(
      paste0("input\\$", id, "\\s*%\\|\\|%\\s*[0-9.]+"), src))
    expect_length(ui, 1L)
    expect_length(fb, 1L)
    ui_val <- as.numeric(sub(".*value = ", "", ui))
    fb_val <- as.numeric(sub(".*%\\|\\|%\\s*", "", fb))
    expect_equal(ui_val, fb_val, info = paste(id, "UI default vs fallback"))
  }
})

test_that("K.15: the dropout control exists and degrades to Inf when blank", {
  src <- power_src()
  expect_match(src, 'ns\\("pw_lmm_dropout_limit"\\)',
               info = "euthanasia limit input must exist")
  expect_match(src, "dropout_limit\\s*=",
               info = "and must be passed to the backend")
  # A blank numericInput yields NA; the backend rejects non-positive values, so
  # the module has to map blank/NA/<=0 to Inf rather than pass it through.
  expect_match(src, "if \\(length\\(dl\\) != 1L \\|\\| is\\.na\\(dl\\) \\|\\| dl <= 0\\) Inf else dl",
               fixed = FALSE,
               info = "blank input must become Inf, not NA")
})

test_that("K.14: the UI no longer claims random intercept SD drives power", {
  # It provably does not: power is identical across a 30-fold range, because the
  # Treatment:Day slope contrast is orthogonal to per-mouse intercepts. The old
  # text sent users to calibrate the one dial that cannot move the answer.
  src <- power_src()
  expect_false(grepl("Higher values require larger N", src, fixed = TRUE),
               info = "the false claim about random intercept SD must stay gone")
  expect_match(src, "Does not affect the power estimate")
  expect_match(src, "dominant driver of power",
               info = "random_slope_sd is what actually matters")
})

test_that("K.15: attrition is carried into the results object", {
  src <- power_src()
  # R20.D6: the simulation's result arrives in .finalize_power_success().
  expect_match(src, "lmm_attrition\\s*=\\s*result(_lmm)?\\$attrition")
  expect_match(src, "Dropout modelled:",
               info = "Summary tab must report realised attrition")
  expect_match(src, "none modelled",
               info = "and must state the assumption when there is no dropout")
})
