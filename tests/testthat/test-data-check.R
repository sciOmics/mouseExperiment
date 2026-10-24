# Round 22 proposal F5 (CODE_REVIEW.md R22-G, 0.45.0): a check of the data on
# upload, with the flagged rows.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}

# Three arms of four animals, measured on days 0-20 every 4 days.
clean_study <- function() {
  days <- seq(0, 20, 4)
  do.call(rbind, lapply(c("Vehicle", "A", "B"), function(arm) do.call(rbind, lapply(1:4, function(i) {
    r <- switch(arm, Vehicle = 0.12, A = 0.08, B = 0.05) + i / 200
    data.frame(ID = paste0(substr(arm, 1, 1), i), Treatment = arm, Cage = paste0(arm, "-C", (i > 2) + 1),
               Day = days, Volume = round(100 * exp(r * days), 1),
               Weight = round(22 - 0.02 * days, 2), Event = 0)
  }))))
}
checks_of <- function(qc) sort(unique(qc$flagged$Check))

test_that("F5: a clean study passes, and each planted error is found on its row", {
  expect_identical(nrow(check_study_data(clean_study(), cage_column = "Cage",
                                         weight_column = "Weight", event_column = "Event")$summary), 0L)
  d <- clean_study()
  d$Line <- seq_len(nrow(d)) + 1L
  d$Volume[3] <- -d$Volume[3]                               # sign error
  d <- rbind(d, transform(d[8, ], Line = 999L))             # entered twice
  d$Volume[20] <- d$Volume[20] * 6                          # a decimal slip
  d$Weight[30] <- d$Weight[30] * 0.7                        # 30 % weight drop
  d$ID[40] <- NA                                            # unplaceable
  d$Event[45] <- 1                                          # measured after its event
  d <- d[-50, ]                                             # a missed day (B3 day 4)
  qc <- check_study_data(d, cage_column = "Cage", weight_column = "Weight", event_column = "Event",
                         row_column = "Line")
  expect_true(all(c("Negative volume", "Duplicate measurement", "Large jump", "Weight drop",
                    "Missing ID, arm or day", "Measured after its event", "Missed measurement") %in%
                  checks_of(qc)))
  f <- qc$flagged
  expect_identical(sort(f$Row[f$Check == "Duplicate measurement"]), c(9L, 999L))
  expect_identical(f$Row[f$Check == "Negative volume"], 4L)
  expect_true(21L %in% f$Row[f$Check == "Large jump"])
  expect_identical(f$Row[f$Check == "Weight drop"], 31L)
  expect_identical(f$Severity[f$Check == "Negative volume"], "Problem")
  missed <- f[f$Check == "Missed measurement", ]
  # B1's deleted day 4, and A3's day 12, whose row lost its ID.
  expect_identical(sort(paste(missed$ID, missed$Day)), c("A3 12", "B1 4"))
  expect_true(all(is.na(missed$Row)))
  expect_identical(qc$summary$Severity[1], "Problem")   # problems first
})

test_that("F5: design notes, calipers as measured, and rows without a volume", {
  d <- clean_study()
  d$ID <- sub("^[A-Z]", "M", d$ID)                         # tags reused in every arm
  d$Cage <- sub("^[A-Za-z]+-", "", d$Cage)                 # cages numbered within each arm
  d$Length <- 10; d$Width <- 8; d$Width[7] <- 12           # a reversed pair
  d$Volume[12] <- NA
  qc <- check_study_data(d, cage_column = "Cage", length_column = "Length", width_column = "Width")
  expect_true(all(c("IDs reused across arms", "Cage labels reused across arms",
                    "Width above length", "No volume") %in% checks_of(qc)))
  expect_identical(qc$summary$Severity[qc$summary$Check == "IDs reused across arms"], "Note")
  # A row with no volume is not also a missed measurement.
  expect_false("Missed measurement" %in% checks_of(qc))
})

test_that("F5: the Data tab shows the check and downloads the flagged rows with their lines", {
  d <- clean_study()
  d$Volume[3] <- -d$Volume[3]
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    rv$uploaded_data <- d
    rv$data_valid <- TRUE
    session$setInputs(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage",
                      col_day = "Day", time_method = "direct_day", volume_method = "direct",
                      col_volume = "Volume", col_weight = "Weight", col_survival_event = "Event")
    ui <- html_of(output$data_check_ui)
    expect_match(ui, "Data check", fixed = TRUE)
    expect_match(ui, "Negative volume", fixed = TRUE)
    expect_match(ui, "Possible errors", fixed = TRUE)
    csv <- utils::read.csv(output$download_data_check)
    expect_identical(csv$Line[csv$Check == "Negative volume"], 4L)   # the header is line 1
  })
})
