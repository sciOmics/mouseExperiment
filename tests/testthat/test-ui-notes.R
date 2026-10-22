# The maintainer's notes from using the dashboard (2026-10-22; 26.10.007,
# CODE_REVIEW.md R20-D.W). Each test fails on 26.10.006.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
html <- function(x) paste(as.character(x), collapse = "\n")

test_that("the About tab gives a contact address", {
  expect_match(html(g("citationUI")()), 'href="mailto:david@sciomics.bio"', fixed = TRUE)
})

test_that("the Survival model menu says Cox Proportional Hazards", {
  ui <- html(g("survivalUI")("surv"))
  expect_match(ui, ">Cox Proportional Hazards</option>", fixed = TRUE)
  expect_false(grepl("Frequentist (Cox PH)", ui, fixed = TRUE))
})

test_that("Toxicity's estimation method explains REML and ML", {
  ui <- html(g("toxicityUI")("tox"))
  expect_match(ui, "Estimation method[^<]*<span class=\"help-tip[^>]*title=\"REML \\(restricted")
})

test_that("each group's colour and shape fit the sidebar: two lines, a symbol menu", {
  ui <- html(g("group_order_and_color_ui")("tg", c("Vehicle - Dose: 0", "Drug_A Mid + Drug_B - Dose: 15")))
  expect_match(ui, "group-style-controls", fixed = TRUE)
  # The shape menu filled 125 px with the shape names; it is the symbol now.
  expect_false(grepl("width: 125px", ui, fixed = TRUE))
  expect_false(grepl("Triangle", ui, fixed = TRUE))
  expect_match(ui, ">\u25b2</option>", fixed = TRUE)
})

test_that("every 'How to read this plot' note sits above its plot, in the plot_help style", {
  for (ui in list(html(g("tumorGrowthUI")("tg")), html(g("toxicityUI")("tox")))) {
    expect_false(grepl("<summary>\\s*<small>How to read", ui))
    # The residual Q-Q note comes before the plot, not a screen below it.
    qq_note <- regexpr("How to read this plot", ui, fixed = TRUE)
    qq_plot <- regexpr("diag_qq\"", ui, fixed = TRUE)
    expect_true(qq_note > 0 && qq_plot > 0 && qq_note < qq_plot)
  }
  # The MCMC trace plots, which had none, have one.
  for (ui in list(html(g("tumorGrowthUI")("tg")), html(g("survivalUI")("surv")))) {
    expect_match(ui, "MCMC trace plots</strong>[\\s\\S]*fuzzy", perl = TRUE)
  }
})

test_that("the Tumour Growth results have no Compare tab", {
  # It showed the previous and current runs' tables side by side, only when the
  # model type changed, and so put LME4 rate differences beside AUC ratios.
  ui <- html(g("tumorGrowthUI")("tg"))
  expect_false(grepl('data-value="Compare"', ui, fixed = TRUE))
  expect_false(grepl("tg_compare_panel", ui, fixed = TRUE))
})

test_that("the AUC help says how animals that die early are handled", {
  tg <- html(g("tumorGrowthUI")("tg"))
  # The sidebar tooltip, and the Info tab's fuller account with its limits.
  expect_match(tg, "animals that die early still count through all their measurements", fixed = TRUE)
  expect_match(tg, "Animals that die or are removed early", fixed = TRUE)
  expect_match(tg, "read the Survival analysis alongside the AUC", fixed = TRUE)
  expect_match(html(g("helpUI")()), "early deaths of large tumours do not pull the AUC down", fixed = TRUE)
})

test_that("Power Analysis comes after Drug Synergy in the navbar", {
  app <- suppressMessages(g("dashboard_app")())
  ui  <- get("ui", envir = environment(app$httpHandler))
  if (is.function(ui)) ui <- ui(list())
  h   <- paste(as.character(htmltools::renderTags(ui)$html), collapse = "\n")
  nav <- regmatches(h, regexpr('<ul class="nav navbar-nav[^"]*"[^>]*>[\\s\\S]*?</ul>', h, perl = TRUE))
  tabs <- gsub('^data-value="|"$', "", regmatches(nav, gregexpr('data-value="[^"]+"', nav))[[1]])
  # 26.10.006 put Power Analysis second, between Data and Tumor Growth.
  expect_identical(tabs, c("Data", "Tumor Growth", "Survival", "Toxicity", "Dose Response",
                           "Drug Synergy", "Power Analysis", "Help", "About"))
})

test_that("26.10.009: the Help tab matches the app", {
  h <- html(g("helpUI")())
  # Features it did not mention.
  for (s in c("Download Results", "Save Preset", "resume or cancel", "How to read this plot",
              "up to 30 MB")) expect_match(h, s, fixed = TRUE, info = s)
  # The endpoint estimands by the names the sidebars use.
  for (s in c("Each animal's last observation", "Survivors")) expect_match(h, s, fixed = TRUE)
  expect_false(grepl("‘Observed’", h, fixed = TRUE))
  # A survival-time column is detected by name but never used.
  expect_false(grepl("Survival_Time, Time_To_Event", h, fixed = TRUE))
  # Bliss is independent action, not additivity.
  expect_false(grepl("expected additive effect", h, fixed = TRUE))
  expect_match(h, "FE_A + FE_B", fixed = TRUE)
  # The Synergy Info tab: one reference model, described correctly, and no
  # need to run Tumor Growth first.
  syn <- html(g("synergyUI")("syn"))
  expect_false(grepl("sum of individual TGIs", syn, fixed = TRUE))
  expect_false(grepl("Two reference models", syn, fixed = TRUE))
  expect_false(grepl("must be run first", syn, fixed = TRUE))
  # Power's analytic mode is a t-test per comparison (no ANOVA since 26.09.002).
  expect_false(grepl("one-way ANOVA from an effect", html(g("powerAnalysisUI")("pw")), fixed = TRUE))
})

test_that("26.10.009: Help has a quick start, the evaluable-day rule and a glossary", {
  h <- html(g("helpUI")())
  expect_match(h, "Quick start", fixed = TRUE)
  expect_match(h, "Load Demo Data", fixed = TRUE)
  # The rule is the same text the Study Design panel shows.
  expect_match(h, g("evaluable_rule_text")(), fixed = TRUE)
  for (term in c("Reference group", "Geometric mean", "TGI (tumour growth inhibition)",
                 "Estimand", "Random effects (random slopes)", "Adjusted and unadjusted p-values",
                 "Confidence and credible intervals", "Hazard ratio and time ratio",
                 "Bliss difference")) {
    expect_match(h, paste0("<dt>", term, "</dt>"), fixed = TRUE, info = term)
  }
  # Release history is not user guidance.
  expect_false(grepl("Versions before 26.09.002", html(g("powerAnalysisUI")("pw")), fixed = TRUE))
})
