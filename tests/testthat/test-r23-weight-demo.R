# The Weight demo shows toxicity (CODE_REVIEW.md R23.24), and the Toxicity
# rates table says where the test is (R23.22).

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

test_that("R23.24: the combination is severely toxic and the single agents mild", {
  d <- suppressMessages(calculate_dates(weight_synthetic_data, reference_date = "March 1", year = 2025))
  names(d) <- make.names(names(d))
  d <- d[order(d$Treatment, d$Ear.Tag, d$Day), ]
  d$pct <- 100 * (d$Mouse.Weight / ave(d$Mouse.Weight, d$Treatment, d$Ear.Tag,
                                       FUN = function(w) w[1]) - 1)
  worst <- aggregate(pct ~ Treatment + Ear.Tag, d, min)
  by_arm <- split(worst$pct, worst$Treatment)
  expect_identical(sum(by_arm[["Drug A + B"]] <= -20), 3L)
  expect_true(all(by_arm[["Drug A"]] > -5))
  expect_true(all(by_arm[["Drug B"]] > -10 & by_arm[["Drug B"]] < -5))
  expect_true(all(by_arm[["DMSO"]] > -3))
  # Two animals leave the study once past 20 %: no 15 March rows.
  expect_identical(nrow(weight_synthetic_data), 98L)
  last <- aggregate(Day ~ Treatment + Ear.Tag, d, max)
  expect_identical(sort(last$Ear.Tag[last$Day < 14]), c("N", "R"))
  expect_true(all(last$Treatment[last$Day < 14] == "Drug A + B"))
  # The package's data and the dashboard's demo file are the same.
  csv <- utils::read.csv(system.file("sample_data", "weight_synthetic_data.csv",
                                     package = "mouseExperiment"), check.names = FALSE)
  expect_equal(csv, as.data.frame(weight_synthetic_data), ignore_attr = TRUE)
})

test_that("R23.24: the threshold analysis finds the three, in the combination", {
  d <- suppressMessages(calculate_dates(weight_synthetic_data, reference_date = "March 1", year = 2025))
  wl <- suppressWarnings(weight_loss_threshold(d, threshold = 0.2, time_column = "Day",
                                               id_column = "Ear Tag", treatment_column = "Treatment",
                                               weight_column = "Mouse Weight", reference_group = "DMSO"))
  ev <- wl$event_data
  expect_identical(sum(ev$Event), 3L)
  expect_identical(unique(as.character(ev$Treatment[ev$Event == 1])), "Drug A + B")
})

test_that("R23.22: the rates table points to the comparison that tests them", {
  h <- paste(as.character(g("toxicityUI")("m")), collapse = "")
  expect_match(h, "differs from the reference is tested in the next table", fixed = TRUE)
})
