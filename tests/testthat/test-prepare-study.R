# The dashboard's data preparation as package functions (CODE_REVIEW.md
# R20-AB): read_study_csv(), apply_exclusions(), study_days(), label_arms()
# and survival_events(), which the dashboard now calls itself.

test_that("R20-AB: read_study_csv() reads semicolons with decimal commas", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("ID;Day;Volume;Treatment", "1;0;54,5;Vehicle", "1;3;98,25;Vehicle"), f)
  d <- read_study_csv(f)
  expect_identical(names(d), c("ID", "Day", "Volume", "Treatment"))
  expect_equal(d$Volume, c(54.5, 98.25))
  writeLines(c("ID,Day,Treatment", "1,0,Vehicle", "2,0,Drug, A"), f)
  expect_error(read_study_csv(f), "Line 3 has 4 fields")
})

test_that("R20-AB: apply_exclusions() removes animals by arm, ID and cage, and arms and cages", {
  d <- data.frame(Treatment = rep(c("Ctl", "DrugA"), each = 4),
                  Dose = rep(c(0, 10), each = 4),
                  ID = rep(c("1", "2"), 4), Cage = rep(c("C1", "C1", "C2", "C2"), 2),
                  Day = 1:8, stringsAsFactors = FALSE)
  # Ear tag 1 in cage C1 of the control arm only: the same tag elsewhere stays.
  out <- apply_exclusions(d, animals = data.frame(arm = "Ctl", id = "1", cage = "C1"),
                          cage_column = "Cage")
  expect_identical(nrow(out), 7L)
  expect_false(any(out$Treatment == "Ctl" & out$ID == "1" & out$Cage == "C1"))
  # With a dose column, the arm is "Treatment Dose".
  out <- apply_exclusions(d, animals = data.frame(arm = "DrugA 10", id = "2", cage = "C2"),
                          cage_column = "Cage", dose_column = "Dose")
  expect_identical(nrow(out), 7L)
  expect_identical(nrow(apply_exclusions(d, treatments = "DrugA")), 4L)
  expect_identical(nrow(apply_exclusions(d, cages = "C2", cage_column = "Cage")), 4L)
  expect_error(apply_exclusions(d, animals = data.frame(arm = "Ctl", id = "1"),
                                cage_column = "Cage"), "and cage")
  expect_identical(apply_exclusions(d), d)
})

test_that("R20-AB: the subset filter is apply_exclusions() with the Data Upload settings", {
  d <- utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                   package = "mouseExperiment"), stringsAsFactors = FALSE)
  key <- animal_key(paste(d$Treatment, d$Dose), d$ID, d$Cage)[c(1, 200)]
  rv <- shiny::reactiveValues(
    exclude_ids = key, exclude_treatments = "Drug_B", exclude_cages = "C03",
    settings = list(col_treatment = "Treatment", col_id = "ID", col_cage = "Cage",
                    col_dose = "Dose"))
  args <- shiny::isolate(subset_filter_args(d, rv))
  expect_equal(args$animals$arm, paste(d$Treatment[c(1, 200)], d$Dose[c(1, 200)]))
  out <- shiny::isolate(apply_subset_filter(d, rv))
  expect_identical(out, do.call(apply_exclusions, c(list(df = d), args)))
  gone <- animal_key(paste(out$Treatment, out$Dose), out$ID, out$Cage)
  expect_false(any(key %in% gone))
  expect_false(any(out$Treatment == "Drug_B") || any(out$Cage == "C03"))
})

test_that("R20-AB: study_days() returns the data, and warns about what it could not read", {
  d <- data.frame(ID = 1:3, Date = c("24/03/2025", "27/03/2025", "31/03/2025"))
  out <- study_days(d, "Date", "2025-03-24")
  expect_equal(out$Day, c(0, 3, 7))
  expect_identical(out$ID, d$ID)
  d2 <- data.frame(ID = 1:3, Date = c("2025-03-24", "not a date", "2025-03-31"))
  expect_warning(out2 <- study_days(d2, "Date", "2025-03-24"), "1 ")
  expect_true(is.na(out2$Day[2]))
  # Ambiguous d/d/yyyy dates stop unless the order is given.
  d3 <- data.frame(Date = c("01/03/2025", "04/03/2025"))
  expect_error(study_days(d3, "Date", "2025-03-01"), "date order")
  expect_equal(study_days(d3, "Date", "2025-03-01", "dmy")$Day, c(0, 3))
})

test_that("R20-AB: label_arms() labels arms by treatment and dose", {
  d <- data.frame(Treatment = c("A", "A", "B"), Dose = c(1, 10, 0))
  expect_identical(label_arms(d, "Treatment", "Dose")$Arm, c("A 1", "A 10", "B 0"))
  expect_identical(label_arms(d, "Treatment", NULL, name = "G")$G, c("A", "A", "B"))
  expect_error(label_arms(d, "Nope"), "not in the data")
})

test_that("R20-AB: survival_events() gives each animal its time and event", {
  d <- data.frame(Treatment = rep(c("Ctl", "Drug"), each = 6),
                  ID = rep(rep(c("1", "2"), each = 3), 2),
                  Day = rep(c(0, 5, 10), 4),
                  Volume = c(100, 900, 1600,  100, 300, 500,  100, 200, 300,  100, 1200, 2100),
                  Ev = c(0, 0, 1,  0, 0, 0,  0, 0, 0,  0, 1, 1),
                  Reason = c("", "", "Volume limit", "", "", "Found dead", "", "", "", "", "", ""),
                  stringsAsFactors = FALSE)
  # From the event column: the first flagged day.
  e <- survival_events(d, event_column = "Ev")
  e <- e[order(e$Treatment, e$ID), ]
  expect_equal(e$Time, c(10, 10, 10, 5))
  expect_equal(e$Event, c(1, 0, 0, 1))
  expect_true(is.na(attr(e, "death_events_added")))
  # From a volume limit: the first day at or over it; others censored at their last day.
  v <- survival_events(d, volume_column = "Volume", volume_limit = 1500)
  v <- v[order(v$Treatment, v$ID), ]
  expect_equal(v$Time, c(10, 10, 10, 10))
  expect_equal(v$Event, c(1, 0, 0, 1))
  # A death from another cause becomes an event at the animal's last day.
  w <- survival_events(d, event_column = "Ev", removal_reason_column = "Reason",
                       death_reasons = "Found dead")
  w <- w[order(w$Treatment, w$ID), ]
  expect_equal(w$Event, c(1, 1, 0, 1))
  expect_identical(attr(w, "death_events_added"), 1L)
  expect_error(survival_events(d), "event column")
  expect_error(survival_events(d, volume_limit = 1500), "volume column")
})
