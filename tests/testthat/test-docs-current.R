# =============================================================================
# Info / Help currency — CODE_REVIEW.md D.13
#
# The Info and Help tabs had drifted badly from the app: the Tumor Growth Info
# announced "two statistical modeling approaches" when there were five, described
# the LME4 ANOVA as an F-test when it reports chi-square, and documented
# extrapolation as feeding the statistical model after that was deliberately
# removed. The Help tab described a Bayesian Power mode that does not exist and
# claimed the AUC path used Kruskal-Wallis when it uses a Welch omnibus.
#
# Nothing detected any of this because prose has no failure mode. It is not
# executed, so it cannot be wrong in a way a test suite notices — it just quietly
# stops being true, which is worse than an error because users act on it.
#
# These tests do what can reasonably be automated: assert that documented option
# counts match the actual `selectInput` choices, that terms describing removed or
# never-present functionality do not reappear, and that newly added features are
# mentioned somewhere. They cannot verify that the prose is *good*, only that its
# checkable claims are still true.
# =============================================================================

read_src <- function(f) {
  p <- testthat::test_path("..", "..", "R", f)
  skip_if_not(file.exists(p), paste(f, "not locatable"))
  paste(readLines(p, warn = FALSE), collapse = "\n")
}

# Pull the values out of a selectInput's `choices = c(...)` block.
choice_values <- function(src, input_id) {
  blk <- regmatches(src, regexpr(
    paste0('ns\\("', input_id, '"\\)[\\s\\S]*?choices\\s*=\\s*c\\(([\\s\\S]*?)\\)\\s*,'),
    src, perl = TRUE))
  if (!length(blk)) return(character(0))
  vals <- regmatches(blk, gregexpr('=\\s*"[a-z_0-9]+"', blk))[[1]]
  unique(sub('.*"([^"]+)".*', "\\1", vals))
}

# ---- Tumor Growth ----------------------------------------------------------

test_that("D.13: the Info tab's model count matches the model selector", {
  src <- read_src("ui_modules.R")
  n_models <- length(choice_values(src, "tg_model_type"))
  expect_gt(n_models, 0L)

  # The Info tab states the count in words. Keep them in step.
  words <- c("2" = "two", "3" = "three", "4" = "four", "5" = "five",
             "6" = "six", "7" = "seven")
  expected <- words[[as.character(n_models)]]
  expect_match(src, paste0(expected, " statistical modeling approaches"),
               info = paste("model selector offers", n_models, "options"))
})

test_that("D.13: every model type in the selector is named in the Info tab", {
  src <- read_src("ui_modules.R")
  # A mode users can pick but that the documentation never mentions is the gap
  # this whole round was about.
  labels <- c(lme4 = "LME4", auc = "Area Under Curve", bayes = "Bayesian LMM")
  for (key in choice_values(src, "tg_model_type")) {
    lbl <- labels[[key]]
    if (is.null(lbl)) next  # a new mode: add it to this map and to the Info tab
    expect_match(src, lbl, fixed = TRUE,
                 info = paste("model", key, "is selectable but undocumented"))
  }
})

test_that("D.13: the ANOVA tab is not described as an F-test", {
  # The LME4 path reports Chisq / Df / Pr(>Chisq) — a Wald chi-square. Calling it
  # an F-test misstates both the statistic and why it is used.
  src <- read_src("ui_modules.R")
  expect_false(grepl("Type III ANOVA F-test", src, fixed = TRUE))
  expect_match(src, "Wald chi-square", fixed = TRUE)
})

test_that("D.13: extrapolation is documented as affecting plots only", {
  # It was removed from the statistical path deliberately: appending an imputed
  # observation and treating it as measured made standard errors ignore the
  # imputation. The docs described the old behaviour for several releases.
  src <- read_src("ui_modules.R")
  expect_false(grepl("LOCF (last observation carried forward) can extrapolate AUC",
                     src, fixed = TRUE))
  expect_match(src, "growth-trajectory plots only", fixed = TRUE)
})

test_that("D.13: the randomisation test and correction options are documented", {
  src <- read_src("ui_modules.R")
  # Both were added to the sidebar without ever reaching the Info tab.
  expect_match(src, "Randomisation Test", fixed = TRUE)
  expect_match(src, "unit that was randomised", fixed = TRUE)
  for (corr in c("Tukey HSD", "Bonferroni", "Holm", "Benjamini-Hochberg FDR")) {
    expect_match(src, corr, fixed = TRUE, info = corr)
  }
})

test_that("D.13: structure-aware cage handling is documented, not the legacy rule", {
  src <- read_src("ui_modules.R")
  for (term in c("Crossed", "Nested, replicated", "Nested, confounded")) {
    expect_match(src, term, fixed = TRUE, info = term)
  }
})

# ---- Help ------------------------------------------------------------------

test_that("D.13: Help does not describe a Power mode the app does not offer", {
  # The Help tab listed "Bayesian simulation" among the Power modes. The selector
  # offers exactly two. Documenting a third sends users looking for a control
  # that is not there.
  ui  <- read_src("ui_modules.R")
  pw  <- read_src("mod_power.R")
  modes <- choice_values(pw, "pw_mode")
  expect_setequal(modes, c("apriori_analytic", "apriori_lmm"))
  expect_false(grepl("Bayesian simulation", ui, fixed = TRUE),
               label = "Help must not claim a Bayesian power mode")
})

test_that("D.13 / R20.4: Help states the AUC comparison method the backend actually uses", {
  # It said Kruskal-Wallis / pairwise Wilcoxon, then Welch's ANOVA on per-animal
  # AUCs. Since backend v0.26.0 the AUC is the area under each arm's fitted
  # curve, compared by AUC ratios.
  ui <- read_src("ui_modules.R")
  expect_false(grepl("Kruskal-Wallis", ui, fixed = TRUE))
  expect_false(grepl("Welch one-way ANOVA", ui, fixed = TRUE))
  expect_match(ui, "ratio of their AUCs", fixed = TRUE)
  expect_match(ui, "fitted growth curve", fixed = TRUE)
})

test_that("R20.D38: the Synergy Info tab states Bliss, not response additivity", {
  # It documented "Expected TGI = min(FE_A + FE_B, 1.0) x 100", the null R17.2
  # removed, while the backend computes FE_A + FE_B - FE_A * FE_B.
  syn <- read_src("mod_synergy.R")
  expect_false(grepl("min(FE_A + FE_B", syn, fixed = TRUE))
  expect_match(syn, "FE_A + FE_B \\u2212 FE_A \\u00d7 FE_B", fixed = TRUE)
  expect_false(grepl("CI threshold", syn, fixed = TRUE))
})

test_that("R20-K: the Synergy, Toxicity and Dose Response Info tabs state the evaluable-day rule", {
  for (f in c("mod_synergy.R", "mod_toxicity.R", "mod_dose_response.R")) {
    expect_match(read_src(f), "evaluable_rule_text()", fixed = TRUE, info = f)
  }
})

test_that("D.13: Help no longer claims Tukey is the default correction", {
  # The default is an exact Dunnett; the old text said "emmeans with Tukey
  # adjustment", which is the wrong family as well as the wrong correction.
  ui <- read_src("ui_modules.R")
  expect_false(grepl("emmeans with Tukey adjustment", ui, fixed = TRUE))
  expect_match(ui, "exact Dunnett correction", fixed = TRUE)
})

test_that("D.13: every auto-detected column role appears in the Help table", {
  # detect_columns() is the authority. A role it detects but Help omits leaves
  # users unaware the app can use a column they already have — which is how the
  # necrosis, sex and height columns went unmentioned.
  ui <- read_src("ui_modules.R")
  roles <- c("Necrotic", "Sex", "Date", "Survival_Time", "Height", "Removal_Reason")
  for (r in roles) {
    expect_match(ui, r, fixed = TRUE, info = paste("Help omits role:", r))
  }
})

test_that("D.13: Help documents the volume formula choice, not just one formula", {
  # Six formulas are selectable in Data Upload; Help hard-coded the ellipsoid.
  ui <- read_src("ui_modules.R")
  expect_match(ui, "modified ellipsoid", fixed = TRUE)
  expect_match(ui, "six formulas", fixed = TRUE)
})

test_that("D.13: Help lists the demo datasets the Data tab offers", {
  ui  <- read_src("ui_modules.R")
  upl <- read_src("mod_data_upload.R")
  blk <- regmatches(upl, regexpr("demo_files <- list\\([\\s\\S]*?\\n    \\)",
                                 upl, perl = TRUE))
  skip_if(length(blk) != 1L, "demo_files block not locatable")
  labels <- sub('.*"([^"]+)".*', "\\1",
                regmatches(blk, gregexpr('label\\s*=\\s*"[^"]+"', blk))[[1]])
  # Each demo's leading word should appear in the Help demo list.
  for (lbl in labels) {
    key <- sub(" .*", "", lbl)
    expect_match(ui, key, fixed = TRUE, info = paste("Help omits demo:", lbl))
  }
})

# ---- Toxicity / Synergy ----------------------------------------------------

test_that("D.13: the endpoint estimand is documented wherever it is offered", {
  # It changes the headline number, so a user who does not know it exists cannot
  # know which estimand their result reports.
  for (f in c("mod_toxicity.R", "mod_synergy.R")) {
    src <- read_src(f)
    if (!grepl('ns\\("[a-z]+_endpoint_method"\\)', src)) next
    expect_match(src, "Endpoint estimand", fixed = TRUE, info = f)
    expect_match(src, "Observed", fixed = TRUE, info = f)
  }
})

test_that("D.13: interval draws are documented wherever they are offered", {
  # Since backend v0.26.0 the draws are of the fitted model under the model
  # estimand and of the animals otherwise, so the control is "Interval draws".
  for (f in c("mod_toxicity.R", "mod_synergy.R")) {
    src <- read_src(f)
    if (!grepl('ns\\("[a-z]+_n_boot"\\)', src)) next
    expect_match(src, "Interval draws", fixed = TRUE, info = f)
  }
})

test_that("D.13: competing risks are documented in the toxicity KM section", {
  src <- read_src("mod_toxicity.R")
  expect_match(src, "Competing risks", fixed = TRUE)
  expect_match(src, "Aalen-Johansen", fixed = TRUE)
  # R20.23: which removals compete depends on the removal-reason column.
  expect_match(src, "removal-reason column", fixed = TRUE)
})

test_that("R20-K: the therapeutic window is documented as two axes, not the TWM ratio", {
  # The Info and Help tabs described TWM = TGI / weight loss and its ranking
  # for as long as the ratio existed; v26.09.005 withdrew it.
  tox <- read_src("mod_toxicity.R")
  ui  <- read_src("ui_modules.R")
  for (src in list(tox, ui)) {
    expect_false(grepl("TWM = TGI", src, fixed = TRUE))
    expect_false(grepl("Therapeutic Window Metric (TWM)", src, fixed = TRUE))
  }
  for (flag in c("Tolerated", "Not tolerated", "Unclear")) {
    expect_match(tox, flag, fixed = TRUE, info = flag)
  }
  expect_match(tox, "worst weight loss", fixed = TRUE)
})

test_that("R20.19: the Synergy docs describe the interval-based verdict, not the bands", {
  syn <- read_src("mod_synergy.R")
  expect_false(grepl("Strong Synergy above +0.1", syn, fixed = TRUE))
  expect_match(syn, "Additive (no departure from Bliss detected)", fixed = TRUE)
})

test_that("R20.D47 / R20.D49 / R20.D13: Dose Response and Survival document what changed", {
  dr <- read_src("mod_dose_response.R")
  expect_match(dr, "One agent at a time", fixed = TRUE)
  expect_match(dr, "It is two-sided", fixed = TRUE)
  expect_match(dr, "computed on log dose", fixed = TRUE)
  surv <- read_src("mod_survival.R")
  expect_match(surv, "removal-reason column", fixed = TRUE)
  expect_match(surv, "Count deaths from other causes as events", fixed = TRUE)
})

test_that("R14.3: the Tumor Growth tab does not claim to report TGI", {
  # tumor_growth_statistics() returns no TGI field and the tab renders none, but
  # the Info panel stated "All modes compute Tumor Growth Inhibition (TGI) ...".
  # Promising a number the tab does not produce sends users looking for it.
  ui <- read_src("ui_modules.R")
  expect_false(grepl("All modes compute Tumor Growth Inhibition", ui, fixed = TRUE))
  expect_false(grepl("Both modes compute Tumor Growth Inhibition", ui, fixed = TRUE))
  expect_match(ui, "does ", fixed = TRUE)
})

test_that("R14.7: the Bayesian summary is not described as medians or blanket HPD", {
  # brms summary() reports the posterior MEAN with an EQUAL-TAILED interval;
  # emmeans-derived tables report a genuine highest-density interval. Both appear
  # under Lower_95_CrI / Upper_95_CrI, so the tab has to say which is which.
  ui <- read_src("ui_modules.R")
  expect_false(grepl("posterior medians with", ui, fixed = TRUE))
  expect_match(ui, "posterior means", fixed = TRUE)
  expect_match(ui, "equal-tailed", fixed = TRUE)
  expect_match(ui, "highest-density", fixed = TRUE)
})
