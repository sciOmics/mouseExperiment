# =============================================================================
# Dosing-schedule annotation, dashboard side — CODE_REVIEW.md D.15
#
# `plot_treatments()` and its two schedule datasets were removed. They required a
# separate uploaded file and a second panel whose x-axis had to be kept aligned
# with the growth plot by hand — and the files shipped here for months with
# nothing in the app able to open them (D.12).
#
# The replacement draws dosing days into the plots themselves. These tests cover
# the dashboard half: the controls exist on all three time-axis modules, the
# toggle actually gates the annotation, blank input is harmless, and parse
# warnings reach the user rather than the log.
# =============================================================================

suppressMessages(library(ggplot2))

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")

src_of <- function(f) {
  p <- testthat::test_path("..", "..", "R", f)
  skip_if_not(file.exists(p), paste(f, "not locatable"))
  paste(readLines(p, warn = FALSE), collapse = "\n")
}

demo_plot <- function() {
  set.seed(4)
  df <- do.call(rbind, lapply(c("Control", "DrugA"), function(tx) {
    d <- rep(seq(0, 21, 3), 3)
    data.frame(Day = d, Treatment = tx,
               Volume = exp(log(120) + 0.1 * d + stats::rnorm(length(d), 0, 0.1)),
               stringsAsFactors = FALSE)
  }))
  list(df = df,
       p = ggplot(df, aes(Day, Volume, colour = Treatment)) +
         geom_point() + theme_classic())
}
n_mark_layers <- function(p) {
  g <- ggplot_build(p)
  sum(vapply(g$data, function(L)
    if (!is.null(L$shape) && any(L$shape == 25)) 1L else 0L, integer(1)))
}

# ---- the removal -----------------------------------------------------------

test_that("D.15: nothing still calls the removed plot_treatments()", {
  # Comments naming it are fine and in fact wanted — they explain what replaced
  # it. What must not survive is a call, which would now be an error at runtime.
  for (f in list.files(testthat::test_path("..", "..", "R"), full.names = TRUE)) {
    code <- sub("#.*$", "", readLines(f, warn = FALSE))
    expect_false(any(grepl("plot_treatments\\s*\\(", code)),
                 info = paste("live call remains in", basename(f)))
  }
  expect_false("plot_treatments" %in% getNamespaceExports("mouseExperiment"))
})

test_that("D.15: the schedule CSVs are gone from sample_data", {
  d <- system.file("sample_data", package = "mouseExperimentDashboard")
  if (!nzchar(d) || !dir.exists(d)) {
    d <- testthat::test_path("..", "..", "inst", "sample_data")
  }
  skip_if_not(dir.exists(d), "sample_data not locatable")
  expect_length(list.files(d, pattern = "schedule"), 0L)
})

# ---- controls present on every time-axis module ----------------------------

test_that("D.15: all three time-axis modules expose the dosing controls", {
  # The user asked for growth, body weight and survival. A module that renders a
  # day axis but has no control is the gap this test exists to catch.
  expect_match(src_of("ui_modules.R"),  'dose_schedule_ui\\(ns, "tg"\\)')
  expect_match(src_of("mod_toxicity.R"), 'dose_schedule_ui\\(ns, "tox"\\)')
  expect_match(src_of("mod_survival.R"), 'dose_schedule_ui\\(ns, "surv"\\)')
})

test_that("D.15: each module applies the annotation to its plot", {
  expect_match(src_of("mod_tumor_growth_plots.R"),
               'apply_dose_schedule\\(p, input, "tg"')
  # Toxicity has two trajectory branches (spaghetti and mean +/- error); both
  # must be annotated or the marks vanish when the user switches plot type.
  tox <- src_of("mod_toxicity_plots.R")
  expect_equal(length(gregexpr('apply_dose_schedule\\(p, input, "tox"',
                               tox)[[1]]), 2L)
  expect_match(src_of("mod_survival.R"), 'apply_dose_schedule\\(km_obj\\$plot')
})

test_that("D.15: the controls are gated behind a toggle", {
  ui <- as.character(get("dose_schedule_ui", envir = .pkg)(shiny::NS("m"), "tg"))
  expect_match(ui, "m-tg_show_dosing", fixed = TRUE)
  expect_match(ui, "m-tg_dosing_days", fixed = TRUE)
  expect_match(ui, "m-tg_dosing_style", fixed = TRUE)
  # The field only appears once the box is ticked.
  expect_match(ui, "input.tg_show_dosing == true", fixed = TRUE)
})

# ---- gating behaviour ------------------------------------------------------

test_that("D.15: the toggle actually gates the annotation", {
  apply_fn <- get("apply_dose_schedule", envir = .pkg)
  d <- demo_plot()
  grps <- unique(d$df$Treatment)

  off <- list(tg_show_dosing = FALSE, tg_dosing_days = "1,5,9,13")
  expect_identical(n_mark_layers(apply_fn(d$p, off, "tg", grps)), 0L)

  on <- list(tg_show_dosing = TRUE, tg_dosing_days = "1,5,9,13")
  expect_gt(n_mark_layers(apply_fn(d$p, on, "tg", grps)), 0L)
})

test_that("D.15: toggled on with a blank field leaves the plot untouched", {
  # The obvious way for a user to reach a broken state: tick the box, type
  # nothing. It must render normally rather than error or blank the panel.
  apply_fn <- get("apply_dose_schedule", envir = .pkg)
  d <- demo_plot()
  for (v in list("", "   ", NULL)) {
    p <- apply_fn(d$p, list(tg_show_dosing = TRUE, tg_dosing_days = v), "tg",
                  unique(d$df$Treatment))
    expect_identical(length(p$layers), length(d$p$layers))
  }
})

test_that("D.15: the style selector overrides the automatic choice", {
  apply_fn <- get("apply_dose_schedule", envir = .pkg)
  d <- demo_plot(); grps <- unique(d$df$Treatment)
  dense <- paste(0:20, collapse = ",")

  auto <- apply_fn(d$p, list(tg_show_dosing = TRUE, tg_dosing_days = dense,
                             tg_dosing_style = "auto"), "tg", grps)
  expect_identical(n_mark_layers(auto), 0L)          # auto -> band

  forced <- apply_fn(d$p, list(tg_show_dosing = TRUE, tg_dosing_days = dense,
                               tg_dosing_style = "rug"), "tg", grps)
  expect_gt(n_mark_layers(forced), 0L)
})

test_that("D.15: annotation failure never costs the user their plot", {
  # A malformed schedule or an internal error must degrade to the unannotated
  # plot, not an error panel where the analysis used to be.
  apply_fn <- get("apply_dose_schedule", envir = .pkg)
  d <- demo_plot()
  p <- apply_fn(d$p, list(tg_show_dosing = TRUE,
                          tg_dosing_days = "not a number at all"), "tg",
                unique(d$df$Treatment))
  expect_s3_class(p, "ggplot")
  expect_identical(length(p$layers), length(d$p$layers))
})

# ---- feedback --------------------------------------------------------------

test_that("D.15: parse warnings reach the user instead of the log", {
  # A mistyped arm name silently dropping its row is the K.2 failure this review
  # has found five times. resolve_dose_schedule captures the warning; the module
  # renders it beside the field.
  resolve <- get("resolve_dose_schedule", envir = .pkg)
  res <- resolve(list(tg_show_dosing = TRUE,
                      tg_dosing_days = "Nope: 1,5\nControl: 1,5"), "tg",
                 c("Control", "DrugA"))
  expect_gt(length(res$messages), 0L)
  expect_match(paste(res$messages, collapse = " "), "not found in the data")
  # The valid line still produces a schedule.
  expect_false(is.null(res$schedule))

  msg_ui <- get("dose_schedule_msg_ui", envir = .pkg)
  html <- as.character(msg_ui(list(tg_show_dosing = TRUE,
                                   tg_dosing_days = "Nope: 1,5"), "tg",
                              c("Control")))
  expect_match(html, "alert", fixed = TRUE)
})

test_that("D.15: a clean schedule produces no warning banner", {
  msg_ui <- get("dose_schedule_msg_ui", envir = .pkg)
  expect_null(msg_ui(list(tg_show_dosing = TRUE, tg_dosing_days = "1,5,9,13"),
                     "tg", c("Control", "DrugA")))
  # And nothing at all when the toggle is off.
  expect_null(msg_ui(list(tg_show_dosing = FALSE, tg_dosing_days = "Nope: 1"),
                     "tg", c("Control")))
})

test_that("D.15: all three modules render the message slot", {
  expect_match(src_of("mod_tumor_growth.R"), "tg_dosing_msg", fixed = TRUE)
  expect_match(src_of("mod_toxicity.R"),     "tox_dosing_msg", fixed = TRUE)
  expect_match(src_of("mod_survival.R"),     "surv_dosing_msg", fixed = TRUE)
})

# ---- documentation ---------------------------------------------------------

test_that("D.15: the Info tab documents the feature and its caveat", {
  ui <- src_of("ui_modules.R")
  expect_match(ui, "Dosing schedule", fixed = TRUE)
  expect_match(ui, "Dosing days cannot be read from your data", fixed = TRUE)
  # The over-reading caveat is the part that matters statistically: separation
  # near a dose mark is not evidence the dose caused it.
  expect_match(ui, "Reading dose marks", fixed = TRUE)
  expect_match(ui, "descriptive context, not evidence", fixed = TRUE)
})
