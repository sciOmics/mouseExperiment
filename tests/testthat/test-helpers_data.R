# =============================================================================
# helpers_data.R — CODE_REVIEW.md D.1 / D.9
#
# These functions decide what the data means before any model sees it. A defect
# here is invisible downstream, because every number the user reads is already
# conditioned on it — and 29 of the 39 functions in the former helpers.R had no
# test referencing them at all.
#
# The two that matter most are `pick_reference_group()`, which determines what
# every contrast is measured against, and `synthesize_censor_column()`, which
# fabricates the event indicator the survival analysis consumes. Both are tested
# here for the first time.
# =============================================================================

# ---- pick_reference_group --------------------------------------------------

test_that("D.9: the control arm is found regardless of case or wording", {
  # The reference group is a statistical input, not a display choice: it is the
  # denominator of every TGI and the baseline of every contrast.
  expect_identical(pick_reference_group(c("DrugA", "Vehicle", "DrugB")), "Vehicle")
  expect_identical(pick_reference_group(c("DrugA", "control")), "control")
  expect_identical(pick_reference_group(c("DrugA", "IgG")), "IgG")
  expect_identical(pick_reference_group(c("DrugA", "DMSO")), "DMSO")
  expect_identical(pick_reference_group(c("DrugA", "VEHICLE-1")), "VEHICLE-1")
})

test_that("D.9: it falls back to the first group when nothing matches", {
  expect_identical(pick_reference_group(c("ArmX", "ArmY")), "ArmX")
  expect_identical(pick_reference_group(character(0)), "")
})

test_that("D.9: group order wins over pattern order", {
  # Worth pinning because it is arbitrary but load-bearing: with both a vehicle
  # and a DMSO arm present, which becomes the reference is decided by the order
  # the groups arrive in, not by the order of the internal pattern list.
  expect_identical(pick_reference_group(c("DMSO", "Vehicle")), "DMSO")
  expect_identical(pick_reference_group(c("Vehicle", "DMSO")), "Vehicle")
})

test_that("D.9: matching is on substrings, which can misfire", {
  # `grepl("control", ...)` matches anywhere in the name, so an arm whose name
  # merely contains a pattern is silently promoted to reference. Documented
  # rather than silently fixed: renaming is the user's call, and a "fix" that
  # required exact equality would break the "VEHICLE-1" case above, which is a
  # real and common naming style.
  expect_identical(pick_reference_group(c("Uncontrolled", "Vehicle")),
                   "Uncontrolled")
  # The practical guard is that a real study almost always has a genuine control
  # arm, and it is listed first. The fallback is the same either way.
})

# ---- synthesize_censor_column ----------------------------------------------

# Two animals reach study end, one is removed early.
censor_df <- function() {
  data.frame(
    ID  = c(rep("m1", 4), rep("m2", 4), rep("m3", 2)),
    Day = c(0, 7, 14, 21, 0, 7, 14, 21, 0, 7),
    stringsAsFactors = FALSE
  )
}

test_that("D.9: an animal whose record stops early is flagged, once", {
  out <- synthesize_censor_column(censor_df(), "ID", "Day")
  d <- out$data
  expect_identical(out$n_flagged, 1L)
  expect_identical(sum(d[[out$col]]), 1L)
  # Flagged on the animal's LAST row, not its first, and not on every row.
  m3 <- d[d$ID == "m3", ]
  expect_identical(m3[[out$col]], c(0L, 1L))
  # Animals present at study end are untouched.
  expect_identical(sum(d[d$ID %in% c("m1", "m2"), out$col]), 0L)
})

test_that("D.9: nothing is flagged when every animal reaches study end", {
  df <- censor_df()
  df <- df[df$ID != "m3", ]
  out <- synthesize_censor_column(df, "ID", "Day")
  expect_identical(out$n_flagged, 0L)
  expect_identical(sum(out$data[[out$col]]), 0L)
})

test_that("D.9: study end is taken from the data, not assumed", {
  # If the longest-followed animal stops at day 14, day 14 is study end and
  # nobody reaching it is flagged — the function must not assume a fixed grid.
  df <- data.frame(ID = c("m1", "m1", "m2", "m2"), Day = c(0, 14, 0, 14),
                   stringsAsFactors = FALSE)
  out <- synthesize_censor_column(df, "ID", "Day")
  expect_identical(out$n_flagged, 0L)
})

test_that("D.9: a benign gap in follow-up is indistinguishable from an event", {
  # This is the substantive limitation, asserted so it is a known property rather
  # than a surprise. The rule is "last observation before study end => event". An
  # animal that survived but whose final measurement was skipped — a missed
  # session, a data-entry gap — is coded identically to one that was euthanised.
  #
  # It cannot be fixed here: the information needed to tell the two apart is not
  # in the data frame. What can be done is refuse to let it pass unnoticed, which
  # is why the function reports `n_flagged` for the caller to display. This test
  # pins that contract.
  df <- data.frame(ID = c("alive", "alive", "died", "died"),
                   Day = c(0, 14, 0, 14), stringsAsFactors = FALSE)
  df <- rbind(df, data.frame(ID = "alive", Day = 21))  # only 'alive' reaches 21
  out <- synthesize_censor_column(df, "ID", "Day")
  # 'died' is flagged. So would a survivor with a missed day-21 measurement.
  expect_identical(out$n_flagged, 1L)
  expect_true(out$n_flagged > 0,
              info = "n_flagged must be surfaced so the user can sanity-check it")
})

test_that("D.9: the synthesized column does not collide with existing names", {
  df <- censor_df()
  out <- synthesize_censor_column(df, "ID", "Day")
  expect_true(out$col %in% names(out$data))
  # Original columns survive untouched.
  expect_true(all(c("ID", "Day") %in% names(out$data)))
  expect_identical(nrow(out$data), nrow(df))
})

# ---- column detection ------------------------------------------------------

test_that("D.9: detect_column prefers an exact match over a substring", {
  # Order matters: a dataset with both "Day" and "Day_Of_Week" must resolve to
  # "Day". The two-pass structure (exact, then contains) is the thing being
  # asserted, since a single contains-pass would return whichever came first.
  cn <- c("Day_Of_Week", "Day")
  expect_identical(detect_column(cn, tolower(cn), "day"), "Day")
})

test_that("D.9: detect_column falls back to a substring match", {
  cn <- c("Tumor Length (mm)", "Mouse")
  expect_identical(detect_column(cn, tolower(cn), "length"), "Tumor Length (mm)")
})

test_that("D.9: detect_column returns NULL when nothing matches", {
  cn <- c("Alpha", "Beta")
  expect_null(detect_column(cn, tolower(cn), "volume"))
})

# ---- required-column checks ------------------------------------------------

test_that("D.9: has_col and check_required_cols agree on NULL and empty", {
  df <- data.frame(ID = 1, Day = 2)
  expect_true(has_col(df, "ID"))
  expect_false(has_col(df, "Missing"))
  expect_false(has_col(df, NULL))
  expect_false(has_col(df, ""))

  # check_required_cols returns the ones that failed, so an empty return is the
  # success signal. NULL and "" must count as missing, not as satisfied — this is
  # the branch that stops an unmapped column reaching the backend.
  expect_length(check_required_cols(df, c("ID", "Day")), 0L)
  expect_identical(unname(check_required_cols(df, c("ID", "Nope"))), "Nope")
  expect_length(check_required_cols(df, c("ID", "")), 1L)
})

# ---- apply_subset_filter ---------------------------------------------------

# Exclusions are a statistical act: dropping animals changes every estimate that
# follows. These assert the filter removes exactly what was asked for and nothing
# else, and that a mis-mapped column fails safe (drops nothing) rather than
# silently emptying the dataset.

subset_df <- function() data.frame(
  ID = c("m1", "m2", "m3", "m4"), Treatment = c("Ctl", "Ctl", "DrugA", "DrugA"),
  Cage = c("C1", "C1", "C2", "C2"), Volume = c(100, 110, 90, 95),
  stringsAsFactors = FALSE)

subset_rv <- function(...) {
  args <- list(...)
  shiny::reactiveValues(
    exclude_ids        = args$ids %||% character(0),
    exclude_treatments = args$tx  %||% character(0),
    exclude_cages      = args$cg  %||% character(0),
    settings = list(col_id = args$id_col %||% "ID",
                    col_treatment = args$tx_col %||% "Treatment",
                    col_cage = args$cg_col %||% "Cage"))
}

test_that("D.9: excluding animals, arms and cages each remove exactly their rows", {
  expect_identical(nrow(apply_subset_filter(subset_df(), subset_rv())), 4L)

  by_id <- apply_subset_filter(subset_df(), subset_rv(ids = "m2"))
  expect_identical(by_id$ID, c("m1", "m3", "m4"))

  by_tx <- apply_subset_filter(subset_df(), subset_rv(tx = "DrugA"))
  expect_identical(unique(by_tx$Treatment), "Ctl")

  by_cage <- apply_subset_filter(subset_df(), subset_rv(cg = "C2"))
  expect_identical(unique(by_cage$Cage), "C1")
})

test_that("D.9: exclusions compose rather than override one another", {
  out <- apply_subset_filter(subset_df(), subset_rv(ids = "m1", cg = "C2"))
  expect_identical(out$ID, "m2")
})

test_that("D.9: an unmapped or wrong column drops nothing", {
  # Fail-safe direction matters. If a mis-mapped ID column caused the %in% test
  # to match everything, the filter would silently return zero rows and the
  # analysis would fail somewhere far away from the cause.
  out <- apply_subset_filter(subset_df(), subset_rv(ids = "m1", id_col = "NoSuchCol"))
  expect_identical(nrow(out), 4L)
  out2 <- apply_subset_filter(subset_df(), subset_rv(ids = "m1", id_col = ""))
  expect_identical(nrow(out2), 4L)
})

test_that("D.9: IDs are compared as character, so numeric IDs still match", {
  df <- data.frame(ID = c(1, 2, 3), Treatment = "Ctl", Cage = "C1",
                   stringsAsFactors = FALSE)
  out <- apply_subset_filter(df, subset_rv(ids = "2"))
  expect_identical(as.character(out$ID), c("1", "3"))
})

test_that("D.9: empty and NULL input survive the filter", {
  expect_null(apply_subset_filter(NULL, subset_rv()))
  empty <- subset_df()[0, ]
  expect_identical(nrow(apply_subset_filter(empty, subset_rv())), 0L)
})

# ---- detect_date_format ----------------------------------------------------

# A misread date format silently produces wrong elapsed-day values, which shifts
# every growth rate in the study. Each recognised layout is pinned.

test_that("D.9: each supported date layout is recognised", {
  mk <- function(v) data.frame(Date = v, stringsAsFactors = FALSE)
  expect_identical(detect_date_format(mk("03/24/2025"), "Date", "2025-03-24")$date_format,
                   "%m/%d/%Y")
  expect_identical(detect_date_format(mk("2025-03-24"), "Date", "2025-03-24")$date_format,
                   "%Y-%m-%d")
  expect_identical(detect_date_format(mk("24-Mar"), "Date", "2025-03-24")$date_format,
                   "%d-%b")
  expect_identical(detect_date_format(mk("March 24"), "Date", "2025-03-24")$date_format,
                   "%B %d %Y")
})

test_that("D.9: an unrecognised layout returns NULL rather than guessing", {
  # Guessing would be worse than refusing: a wrong format yields plausible but
  # incorrect day numbers, which no downstream check would catch.
  out <- detect_date_format(data.frame(Date = "24.03.2025", stringsAsFactors = FALSE),
                            "Date", "2025-03-24")
  expect_null(out$date_format)
})

test_that("D.9: the year-less layout is completed from the reference date", {
  # "March 24" carries no year, so the function appends the reference year to the
  # data. Getting this wrong shifts every date by a whole year.
  out <- detect_date_format(data.frame(Date = c("March 24", "March 31"),
                                       stringsAsFactors = FALSE),
                            "Date", "2025-03-24")
  expect_match(out$data$Date[1], "2025")
  expect_identical(out$start_date, "March 24 2025")
})

test_that("D.9: the start date is returned in the same format as the data", {
  # calculate_dates() is then called with both, so a mismatch fails to parse.
  out <- detect_date_format(data.frame(Date = "03/24/2025", stringsAsFactors = FALSE),
                            "Date", "2025-03-24")
  expect_identical(out$start_date, "03/24/2025")
  expect_identical(out$date_format, "%m/%d/%Y")
})
