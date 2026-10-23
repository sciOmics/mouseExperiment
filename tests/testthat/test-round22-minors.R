# Round 22's Minor items (CODE_REVIEW.md R22-D, 0.42.0). Each test fails on 0.41.0.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
html_of <- function(x) paste(as.character(x), collapse = "\n")

# ---- R22.8: no hazard ratios the data cannot support ------------------------------

# Six arms of 8. Weight falls only where `loss` is positive; Vehicle animals
# leave at `vehicle_end` (the volume limit, say).
wl_data <- function(loss, vehicle_end = 28, seed = 1) {
  set.seed(seed)
  arms <- names(loss)
  do.call(rbind, lapply(arms, function(a) do.call(rbind, lapply(1:8, function(i) {
    days <- seq(0, if (a == "Vehicle") vehicle_end else 28, 4)
    w0 <- 22 + stats::rnorm(1, 0, 0.5)
    data.frame(ID = i, Treatment = a, Day = days,
               Weight = w0 * (1 - loss[[a]][i] * days / 28))
  }))))
}
no_loss <- rep(0, 8)

test_that("R22.8: one late weight-loss event, with no Vehicle animal left, gives no Cox table", {
  loss <- list(Vehicle = no_loss, A = c(0.3, rep(0, 7)), B = no_loss, C = no_loss)
  wl <- suppressWarnings(weight_loss_threshold(wl_data(loss, vehicle_end = 12),
    weight_column = "Weight", reference_group = "Vehicle", threshold = 0.2))
  expect_identical(sum(wl$event_data$Event), 1L)
  expect_null(wl$cox_model)
  expect_null(wl$cox_summary)
  expect_match(wl$cox_skipped, "1 weight-loss event", fixed = TRUE)
  expect_false(is.null(wl$km_fit))
})

test_that("R22.8: events with no reference animal at risk give no Cox table either", {
  loss <- list(Vehicle = no_loss, A = c(0.3, 0.3, 0.3, 0.3, rep(0, 4)), B = no_loss)
  wl <- suppressWarnings(weight_loss_threshold(wl_data(loss, vehicle_end = 8),
    weight_column = "Weight", reference_group = "Vehicle", threshold = 0.2))
  expect_identical(sum(wl$event_data$Event), 4L)
  expect_null(wl$cox_model)
  expect_match(wl$cox_skipped, "no Vehicle animal was still on study", fixed = TRUE)
})

test_that("R22.8: with events while the reference is at risk, the Cox model is fitted", {
  loss <- list(Vehicle = c(0.3, rep(0, 7)), A = c(0.3, 0.3, 0.3, rep(0, 5)), B = no_loss)
  wl <- suppressWarnings(weight_loss_threshold(wl_data(loss),
    weight_column = "Weight", reference_group = "Vehicle", threshold = 0.2))
  expect_null(wl$cox_skipped)
  expect_false(is.null(wl$cox_model))
})

# ---- R22.10, R22.15: one animal key; doubling-time rows ---------------------------

test_that("R22.10: animal_key() is make_mouse_key(), with its check for a missing column", {
  expect_identical(g("animal_key")("Vehicle", c(1, 2), "C1"), c("Vehicle|||1|||C1", "Vehicle|||2|||C1"))
  expect_identical(g("animal_key")(c("A", "B"), c(1, 2)), g("make_mouse_key")(c("A", "B"), c(1, 2)))
  # A missing ID column used to give every animal of an arm one key.
  expect_error(g("animal_key")(c("A", "B"), NULL), "a key component is NULL or empty")
})

test_that("R22.10: calculate_dates() points to study_days() for an unknown day/month order", {
  expect_false("in_place" %in% names(formals(calculate_dates)))
  expect_match(dash_source("man", "calculate_dates.Rd"), "link[=study_days]", fixed = TRUE)
})

test_that("R22.15: tumor_doubling_time() has plain row numbers, not 'times'", {
  dt <- suppressWarnings(tumor_doubling_time(master_synthetic_data))
  expect_identical(rownames(dt), as.character(seq_len(nrow(dt))))
  expect_null(names(dt$growth_rate))
})

# ---- R22.11, R22.18: the analytics helpers ------------------------------------------

test_that("R22.11/R22.18: every analytics helper the server calls exists, and no other", {
  js <- paste(readLines(system.file("www", "js", "umami-events.js", package = "mouseExperiment")),
              collapse = "\n")
  defined <- regmatches(js, gregexpr("window\\.track[A-Za-z]+(?= = function)", js, perl = TRUE))[[1]]
  defined <- sub("window.", "", defined, fixed = TRUE)
  # Without the Umami script it neither polls for it nor logs each event.
  expect_false(grepl("setInterval", js, fixed = TRUE))
  expect_false(grepl("console.log", js, fixed = TRUE))
  src <- dash_all_sources()
  called <- unique(sub("window.", "", regmatches(src, gregexpr("window\\.track[A-Za-z]+", src))[[1]],
                       fixed = TRUE))
  expect_setequal(setdiff(defined, "trackEvent"), setdiff(called, "trackEvent"))
})

# ---- R22.13: arm labels in the Bayesian posterior summary ----------------------------

test_that("R22.13: the posterior summary names each arm as labelled, not as brms renamed it", {
  arms <- c("Vehicle", "Drug A + Drug B 10", "Drug-C (5 mg/kg)", "D's")
  par <- c("Intercept", "TreatmentDrugAPDrugB10", "Day", "TreatmentDrugAPDrugB10:Day",
           "TreatmentDrugMC5mgDkg", "TreatmentDs")
  expect_identical(g("bayes_term_labels")(par, arms),
                   c("Intercept", "Treatment: Drug A + Drug B 10", "Day",
                     "Treatment: Drug A + Drug B 10 × Day", "Treatment: Drug-C (5 mg/kg)",
                     "Treatment: D's"))
  # The same renaming as brms's own.
  skip_if_not(exists("rename", envir = asNamespace("brms"), inherits = FALSE))
  expect_identical(g("brms_level_names")(arms), get("rename", envir = asNamespace("brms"))(arms))
})

# ---- R22.14: no version history in the Help ------------------------------------------

test_that("R22.14: the Help explains methods without citing package versions", {
  txt <- html_of(g("helpUI")())
  expect_false(grepl("\\bv0\\.[0-9]+", txt))
  expect_match(txt, "Why there is no Combination Index", fixed = TRUE)
})

# ---- R22.16: one dose of a treatment ---------------------------------------------

test_that("R22.16: apply_exclusions() removes one dose arm, and the Subset Filter offers each", {
  d <- data.frame(ID = 1:6, Treatment = c("Vehicle", "Vehicle", "Drug A", "Drug A", "Drug A", "Drug B"),
                  Dose = c(0, 0, 10, 20, 20, 5), Day = 0, Volume = 100)
  out <- apply_exclusions(d, arms = "Drug A 20", dose_column = "Dose")
  expect_identical(out$ID, c(1L, 2L, 3L, 6L))
  s <- list(col_treatment = "Treatment", col_dose = "Dose", col_id = "ID")
  ch <- g("treatment_filter_choices")(d, s)
  # Each treatment, then the doses of the one given at several.
  expect_identical(unname(ch), c("Drug A", "Drug B", "Vehicle",
                                 animal_key("Drug A", "10"), animal_key("Drug A", "20")))
  expect_identical(names(ch)[4:5], c("Drug A, dose 10", "Drug A, dose 20"))
  expect_identical(unname(g("treatment_filter_choices")(d, utils::modifyList(s, list(col_dose = "")))),
                   c("Drug A", "Drug B", "Vehicle"))
  rv <- shiny::reactiveValues(settings = s, exclude_treatments = c("Drug B", animal_key("Drug A", "20")))
  args <- shiny::isolate(g("subset_filter_args")(d, rv))
  expect_identical(args$treatments, "Drug B")
  expect_identical(args$arms, "Drug A 20")
  expect_identical(do.call(apply_exclusions, c(list(df = d), args))$ID, c(1L, 2L, 3L))
})

# ---- R22.18: .Rbuildignore -------------------------------------------------------

test_that("R22.18: every .Rbuildignore entry names a file that exists or is made locally", {
  src <- dash_source_dir()
  ign <- readLines(file.path(src, ".Rbuildignore"))
  literal <- gsub("\\\\.", ".", sub("^\\^(.*)\\$$", "\\1", ign))
  local <- sub("/$", "", trimws(readLines(file.path(src, ".gitignore"))))
  stale <- literal[!file.exists(file.path(src, literal)) & !literal %in% local]
  expect_identical(stale, character(0))
})
