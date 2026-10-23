# Regression tests for step 3 of the Round 20/21 implementation (CODE_REVIEW.md
# R20-D.K): the Study Design panel, one animal key, declared volume units and the
# day/month date fix. Each test fails on v26.09.002. Where a test drives a
# module, it also sets the inputs that version read (tg_reference,
# surv_reference_group, ...), so that it fails on the defect rather than on a
# missing input.

skip_if_not_installed("shiny", minimum_version = "1.5.0")
skip_if_not_installed("mouseExperiment", minimum_version = "0.25.2")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}

demo_csv <- function(f) {
  utils::read.csv(system.file("sample_data", f, package = "mouseExperiment"),
                  stringsAsFactors = FALSE)
}

master_settings <- function(...) {
  utils::modifyList(list(
    col_id = "ID", col_day = "Day", time_method = "day",
    col_volume = "Volume", volume_method = "precalculated",
    col_treatment = "Treatment", col_cage = "Cage", col_weight = "Weight",
    col_dose = "Dose", col_survival_event = "Survival_Censor",
    control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse"),
    list(...))
}

master_rv <- function(data = demo_csv("all_features_synthetic_data.csv"), ...) {
  shiny::reactiveValues(uploaded_data = data, data_valid = TRUE, session_id = NULL,
                        settings = master_settings(...))
}

# Two arms whose ear tags restart at 1, no cage column. Vehicle animals reach
# the endpoint on days 8-20; DrugA animals later or not at all. DrugA animal 1
# is flagged on day 20 but measured to day 28.
reused_tag_surv <- function() {
  days <- seq(0, 28, by = 4)
  rows <- list()
  for (arm in c("Vehicle", "DrugA")) for (i in 1:10) {
    ev_day <- if (arm == "Vehicle") c(8, 12, 12, 16, 16, 20, 20, 20, 24, 24)[i]
              else c(20, 28, NA, NA, 28, NA, NA, 24, NA, NA)[i]
    last <- if (is.na(ev_day)) 28 else ev_day
    keep <- if (arm == "DrugA" && i == 1) days else days[days <= last]
    rows[[length(rows) + 1L]] <- data.frame(
      ID = i, Treatment = arm, Day = keep,
      Volume = 100 * exp(0.1 * keep),
      Event = as.integer(!is.na(ev_day) & keep == ev_day))
  }
  do.call(rbind, rows)
}

# ---- R20.D1: dates --------------------------------------------------------------

test_that("R20.D1: DD/MM/YYYY dates are read day-first", {
  # The review's case: a study starting 1 March 2025. The old reader gave
  # -57, 33, 155, 247, NA, NA.
  df <- data.frame(Date = c("01/03/2025", "04/03/2025", "08/03/2025",
                            "11/03/2025", "15/03/2025", "18/03/2025"))
  sd <- g("study_days_info")(df, "Date", "2025-03-01")
  expect_equal(sd$data$Day, c(0, 3, 7, 10, 14, 17))
  expect_equal(sd$layout, "DD/MM/YYYY")
  expect_match(sd$reason, "15/03/2025")
})

test_that("R20.D1: an ambiguous column asks for the order, and a declared order is used", {
  df <- data.frame(Date = c("01/03/2025", "04/03/2025", "08/03/2025"))
  expect_error(g("study_days_info")(df, "Date", "2025-03-01"), "Choose the date order")
  expect_equal(g("study_days_info")(df, "Date", "2025-03-01", "dmy")$data$Day, c(0, 3, 7))
  expect_equal(g("study_days_info")(df, "Date", "2025-01-03", "mdy")$data$Day, c(0, 90, 212))
  # A declared order the data contradict is an error, not a silent NA.
  df2 <- data.frame(Date = c("01/03/2025", "15/03/2025"))
  expect_error(g("study_days_info")(df2, "Date", "2025-03-01", "mdy"), "month 15")
})

test_that("R20.D1: year-less dates across New Year do not go negative", {
  df <- data.frame(Date = c("20-Dec", "27-Dec", "03-Jan", "10-Jan"))
  sd <- g("study_days_info")(df, "Date", "2025-12-20")
  expect_equal(sd$data$Day, c(0, 7, 14, 21))
  expect_match(paste(sd$notes, collapse = " "), "crosses New Year")
})

test_that("R20.D1: unreadable dates and days before Day 0 are reported", {
  df <- data.frame(Date = c("2025-03-20", "2025-03-24", "2025-03-31", "n/a"))
  sd <- g("study_days_info")(df, "Date", "2025-03-24")
  expect_equal(sd$data$Day, c(-4, 0, 7, NA))
  notes <- paste(sd$notes, collapse = " ")
  expect_match(notes, "1 row\\(s\\) have no usable date")
  expect_match(notes, "1 measurement\\(s\\) fall before the reference date")
})

test_that("R20.D1: every tab's preprocessing honours the declared date order", {
  df <- data.frame(ID = rep(1:3, each = 3), Treatment = rep(c("A", "B", "C"), each = 3),
                   Date = rep(c("01/03/2025", "04/03/2025", "08/03/2025"), 3),
                   Volume = 100)
  rv <- shiny::reactiveValues(
    uploaded_data = df, data_valid = TRUE,
    settings = list(col_id = "ID", col_treatment = "Treatment", time_method = "from_date",
                    col_date_for_day = "Date", reference_date = "2025-03-01",
                    date_order = "dmy", volume_method = "direct", col_volume = "Volume"))
  pp <- shiny::isolate(g("ensure_preprocessed")(rv))
  expect_equal(unique(pp$data$Day), c(0, 3, 7))     # was 0, 31, 59 read month-first
})

test_that("R20.D1: Data Upload says how the dates were read", {
  df <- data.frame(ID = rep(1:2, each = 2), Treatment = "A",
                   Date = rep(c("01/03/2025", "15/03/2025"), 2), Volume = 100)
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    rv$uploaded_data <- df
    rv$data_valid <- TRUE
    session$setInputs(time_method = "from_date", col_date_for_day = "Date",
                      reference_date = "2025-03-01", date_order = "auto",
                      col_id = "ID", col_treatment = "Treatment",
                      volume_method = "direct", col_volume = "Volume")
    invisible(output$data_summary_ui)          # evaluates the preview
    expect_match(html_of(output$reference_date_status), "DD/MM/YYYY")
    expect_equal(rv$settings$date_order, "auto")
  })
})

# ---- Study Design panel (R20-D.H) -----------------------------------------------

test_that("Study Design: a demo declares mm3 and the control; settings carry the design", {
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(demo_dataset = "all_features", load_demo_data = 1)
    ui <- html_of(output$study_design_ui)
    expect_match(ui, 'value="mm3"[^>]*checked')
    expect_match(ui, '<option value="Vehicle" selected>')
    session$setInputs(col_treatment = "Treatment", col_id = "ID", col_cage = "Cage",
                      volume_units = "mm3", control_group = "Vehicle",
                      randomisation_unit = "cage")
    expect_equal(rv$settings$volume_units, "mm3")
    expect_equal(rv$settings$control_group, "Vehicle")
    expect_equal(rv$settings$randomisation_unit, "cage")
    # Without a cage column the unit is individual animals.
    session$setInputs(col_cage = "")
    expect_equal(rv$settings$randomisation_unit, "mouse")
  })
})

test_that("Study Design: an uploaded volume column starts with units undeclared", {
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    f <- tempfile(fileext = ".csv")
    utils::write.csv(demo_csv("all_features_synthetic_data.csv"), f, row.names = FALSE)
    session$setInputs(data_file = list(name = "study.csv", datapath = f))
    expect_true(rv$data_valid)
    ui <- html_of(output$study_design_ui)
    expect_false(grepl('value="mm3"[^>]*checked', ui))
    expect_false(grepl('value="cm3"[^>]*checked', ui))
    expect_match(html_of(output$volume_units_status), "Not declared")
  })
})

test_that("Study Design: the preset carries the design and the subset filter", {
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(demo_dataset = "all_features", load_demo_data = 1)
    session$setInputs(col_treatment = "Treatment", col_id = "ID", col_cage = "Cage",
                      volume_units = "cm3", control_group = "Drug_B",
                      randomisation_unit = "cage", date_order = "dmy")
    rv$exclude_cages <- "C03"
    saved <- jsonlite::read_json(output$save_preset, simplifyVector = TRUE)
    expect_equal(saved$settings$volume_units, "cm3")
    expect_equal(saved$settings$control_group, "Drug_B")
    expect_equal(saved$settings$randomisation_unit, "cage")
    expect_equal(saved$settings$date_order, "dmy")
    expect_equal(saved$filters$exclude_cages, "C03")

    # Loading it restores the filter.
    rv$exclude_cages <- character(0)
    f <- tempfile(fileext = ".json")
    jsonlite::write_json(saved, f, auto_unbox = TRUE)
    session$setInputs(load_preset = list(name = "p.json", datapath = f))
    expect_equal(rv$exclude_cages, "C03")
  })
})

test_that("Study Design: a changed reference date or date order clears the cache", {
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(demo_dataset = "combo", load_demo_data = 1)
    session$setInputs(col_treatment = "Treatment", col_id = "ID", col_cage = "Cage",
                      time_method = "from_date", col_date_for_day = "Date",
                      reference_date = "2025-03-24", volume_method = "calculate",
                      col_length = "Length", col_width = "Width")
    rv$processed_data <- data.frame(stale = TRUE)
    session$setInputs(reference_date = "2025-03-20")
    expect_false(identical(rv$processed_data, data.frame(stale = TRUE)))
    rv$processed_data <- data.frame(stale = TRUE)
    session$setInputs(date_order = "mdy")
    expect_false(identical(rv$processed_data, data.frame(stale = TRUE)))
  })
})

# ---- Tumour growth ----------------------------------------------------------------

test_that("Study Design: Tumour Growth compares against the study's control group", {
  suppressWarnings(shiny::testServer(g("tumor_growth_server"),
                                     args = list(rv = master_rv(control_group = "Drug_B")), {
    session$setInputs(tg_model_type = "lme4", tg_reference = "Vehicle 0",
                      tg_p_adjust = "dunnett", run_tg_analysis = 1)
    te <- rv$tg_results$treatment_effects
    expect_equal(as.character(te$Group[te$Note == "Reference group"]), "Drug_B 10")
    expect_equal(rv$tg_reference_used, "Drug_B 10")
    expect_match(html_of(output$tg_reference_note), "Drug_B 10")
  }))
})

test_that("R20.D57: the drill-down keeps cage-mates that share an ear tag apart", {
  df <- data.frame(ID = "L", Treatment = "A", Cage = rep(c("C1", "C2"), each = 3),
                   Day = rep(c(0, 3, 7), 2), Volume = 100)
  rows <- g("drilldown_rows")(df, data.frame(Treatment = "A", ID = "L", Cage = "C2"),
                              "ID", "Treatment", "Cage")
  expect_equal(nrow(rows), 3L)                     # 6 before: both cages
  expect_equal(unique(rows$Cage), "C2")
})

# ---- Survival ---------------------------------------------------------------------

test_that("R20.D2: survival keeps ear tags reused across arms apart and uses the first event day", {
  rv <- shiny::reactiveValues(
    uploaded_data = reused_tag_surv(), data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "", col_dose = "",
                    col_survival_event = "Event", control_group = "Vehicle",
                    randomisation_unit = "mouse"))
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = rv), {
    session$setInputs(surv_reference_group = "Vehicle", surv_event_source = "column",
                      run_surv_analysis = 1)
    expect_false(isTRUE(rv$surv_results$error))
    sd <- rv$surv_results$survival_data
    expect_equal(nrow(sd), 20L)                     # 10 before: IDs merged across arms
    expect_equal(sum(sd$Event), 14L)
    expect_equal(as.vector(table(sd$Treatment)[c("Vehicle", "DrugA")]), c(10L, 10L))
    # DrugA animal 1 was flagged on day 20 and measured to day 28: its event
    # time is 20, not 28.
    expect_equal(sort(sd$Time[sd$Treatment == "DrugA" & sd$Event == 1]),
                 c(20, 24, 28, 28))
  }))
})

test_that("R20.D2: the per-animal reduction handles row-level and animal-level flags", {
  df <- data.frame(ID = c(1, 1, 1, 1, 2, 2, 2), Treatment = c(rep("A", 4), rep("B", 3)),
                   Day = c(0, 10, 20, 28, 0, 10, 20),
                   Ev = c(0, 0, 1, 0, 1, 1, 1))
  out <- g("survival_animal_rows")(df, "Day", "Treatment", "ID", NULL, "Ev")
  expect_equal(nrow(out), 2L)
  expect_equal(out$Day[out$Treatment == "A"], 20)    # first flagged day, not 28
  expect_equal(out$Day[out$Treatment == "B"], 20)    # every row flagged: last day
  expect_equal(out$Ev, c(1, 1))
  ev <- g("volume_threshold_events")(
    data.frame(ID = c(1, 1, 1, 1, 1, 1), Treatment = rep(c("A", "B"), each = 3),
               Day = c(0, 5, 10, 0, 5, 10), V = c(10, 60, 80, 10, 20, 30)),
    "Day", "Treatment", "ID", NULL, "V", 50)
  expect_equal(ev$Time, c(5, 10))
  expect_equal(ev$Event, c(1L, 0L))
})

test_that("R20.D14: survival runs with an ID column not named 'ID'", {
  d <- demo_csv("all_features_synthetic_data.csv")
  names(d)[names(d) == "ID"] <- "Animal"
  suppressWarnings(shiny::testServer(g("survival_server"),
                                     args = list(rv = master_rv(d, col_id = "Animal")), {
    session$setInputs(surv_reference_group = "Vehicle 0", run_surv_analysis = 1)
    expect_false(isTRUE(rv$surv_results$error))     # "undefined columns selected" before
    expect_equal(nrow(rv$surv_results$survival_data), 48L)
  }))
})

test_that("R20.6: whole-cage randomisation gives the cage permutation and says so", {
  suppressWarnings(shiny::testServer(g("survival_server"),
                                     args = list(rv = master_rv(randomisation_unit = "cage")), {
    session$setInputs(surv_reference_group = "Vehicle 0", run_surv_analysis = 1)
    res <- rv$surv_results
    expect_false(isTRUE(res$error))
    expect_equal(res$randomisation_unit, "cage")
    expect_true(is.data.frame(res$cage_permutation) && nrow(res$cage_permutation) == 5L)
    panel <- html_of(output$surv_model_quality_panel)
    expect_match(panel, "Whole cages were randomised")
    expect_match(panel, "Smallest attainable p")
    # The floor warning reaches the user.
    expect_match(html_of(output$surv_status_msg), "cannot reach p")
  }))
})

test_that("R20.D13: no unread natural-death checkbox; the assumption is stated", {
  ui <- html_of(g("survivalUI")("s"))
  expect_false(grepl("natural_death_as_event", ui, fixed = TRUE))
  shiny::testServer(g("survival_server"), args = list(rv = master_rv()), {
    session$setInputs(surv_event_source = "volume")
    expect_match(html_of(output$surv_event_assumption), "censored")
    session$setInputs(surv_event_source = "column")
    expect_match(html_of(output$surv_event_assumption), "coded 1")
  })
})

# ---- Toxicity (R20.D34, R20.D31) --------------------------------------------------

weight_rv <- function(units = "mm3") {
  shiny::reactiveValues(
    uploaded_data = demo_csv("weight_synthetic_data.csv"), data_valid = TRUE,
    settings = list(col_id = "Ear.Tag", col_treatment = "Treatment",
                    time_method = "from_date", col_date_for_day = "Date",
                    reference_date = "2025-03-01", date_order = "auto",
                    volume_method = "direct", col_volume = "Tumor.Volume",
                    col_weight = "Mouse.Weight", col_cage = "", col_dose = "",
                    control_group = "DMSO", volume_units = units))
}

run_tox_bw <- function(session) {
  session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE,
                    tox_tumor_density = 1, tox_covariates = "volume",
                    tox_estimation = "REML", tox_threshold = 20,
                    tox_baseline_day = "first", tox_reference_group = "DMSO",
                    run_tox_analysis = 1)
}

test_that("R20.D34: the Weight demo's tumour-mass correction uses the declared mm3", {
  rv <- weight_rv("mm3")
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    run_tox_bw(session)
    res <- rv$tox_results
    expect_false(isTRUE(res$mixed_model$error))
    expect_false(isTRUE(res$weight_loss$error))
    # 20 of 20 animals crossed 20 % loss and net weights reached -1,070 g before.
    expect_equal(sum(res$weight_loss$event_data$Event), 0L)
    td <- res$trajectory_data
    expect_true(all(td$Weight > 10))
    # The plot shows the model's net weight, not a re-derivation.
    wd <- res$mixed_model$weight_data
    expect_equal(sort(td$Weight), sort(wd$Net_Weight[!is.na(wd$Net_Weight)]))
    expect_equal(attr(td, "weight_label"), "Net body weight (g, tumour mass removed)")
  }))
})

test_that("R20.D34: undeclared units block the tumour-mass correction", {
  rv <- weight_rv("")
  shiny::isolate(g("ensure_preprocessed")(rv))
  shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    run_tox_bw(session)
    expect_null(rv$tox_results)
    expect_match(html_of(output$tox_status_msg), "Declare the tumour volume units")
  })
})

test_that("R20.D31: the therapeutic window gets the cage column and the declared units", {
  seen <- NULL
  testthat::local_mocked_bindings(
    therapeutic_window_metric = function(...) { seen <<- list(...); stop("stop here") },
    .package = "mouseExperiment")
  rv <- master_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "et", tox_adjust_tumor = TRUE,
                      tox_tumor_density = 1, tox_endpoint_method = "model",
                      tox_n_boot = 0, tox_reference_group = "Vehicle",
                      run_tox_analysis = 1)
  })
  expect_equal(seen$cage_column, "Cage")
  expect_equal(seen$volume_units, "mm3")
  # R20-AB: with the dose column mapped, arms are "Treatment Dose"; the
  # study's control "Vehicle" is the arm "Vehicle 0".
  expect_equal(seen$reference_group, "Vehicle 0")
  expect_equal(seen$treatment_column, "TreatmentGroup")
})

# ---- Dose-response and synergy -----------------------------------------------------

test_that("Study Design: dose-response uses the study's control and the cage key", {
  seen <- NULL
  testthat::local_mocked_bindings(
    dose_response_statistics = function(...) { seen <<- list(...); stop("stop here") },
    .package = "mouseExperiment")
  rv <- master_rv(control_group = "Vehicle")
  shiny::isolate(g("ensure_preprocessed")(rv))
  shiny::isolate(rv$dose_col <- "Dose")
  shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    # Since v26.09.005 the arms of one agent's series are chosen (R20.D47).
    session$setInputs(dr_control_group = "Drug_B", dr_time_point = "all",
                      dr_treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
                      run_dr_analysis = 1)
  })
  expect_equal(seen$control_group_name, "Vehicle")
  expect_equal(seen$cage_column, "Cage")
})

synergy_rv <- function(id_name = "ID") {
  set.seed(1)
  arms <- c(Control = 0.10, DrugA = 0.07, DrugB = 0.07, Combo = 0.03)
  df <- do.call(rbind, lapply(names(arms), function(arm) do.call(rbind, lapply(1:6, function(i) {
    days <- seq(0, 21, by = 3.5)
    data.frame(ID = paste(arm, i), Treatment = arm, Day = days,
               Volume = exp(log(100) + arms[[arm]] * days + stats::rnorm(length(days), 0, 0.08)))
  }))))
  names(df)[names(df) == "ID"] <- id_name
  shiny::reactiveValues(
    uploaded_data = df, data_valid = TRUE,
    settings = list(col_id = id_name, col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "", col_dose = "",
                    control_group = "Control"))
}

test_that("R20.D14: over-time synergy runs with an ID column not named 'ID'", {
  suppressWarnings(shiny::testServer(g("synergy_server"),
                                     args = list(rv = synergy_rv("Animal")), {
    session$setInputs(syn_drug_a = "DrugA", syn_drug_b = "DrugB", syn_combo = "Combo",
                      syn_control = "Control", syn_time_point = "all",
                      syn_endpoint_method = "model", syn_n_boot = 0,
                      run_syn_analysis = 1)
    # "Could not calculate synergy for any time points" before.
    expect_false(isTRUE(rv$syn_results$error))
    expect_gt(nrow(rv$syn_results$synergy_summary), 3L)
    expect_equal(rv$syn_plot_inputs$control, "Control")
  }))
})

test_that("R20.D36: after a Tumour Growth run with a dose column, Synergy's own choices run", {
  rv <- master_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_reference = "Vehicle 0",
                      tg_p_adjust = "dunnett", run_tg_analysis = 1)
  }))
  # R20-AB: Tumor Growth labels arms on its own copy; the shared data keep the
  # treatment column, so other tabs no longer depend on whether it ran.
  expect_equal(shiny::isolate(rv$tx_col), "Treatment")
  expect_false("TreatmentGroup" %in% names(shiny::isolate(rv$processed_data)))
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = rv), {
    ui <- html_of(output$syn_group_selectors_ui)
    opts <- unique(regmatches(ui, gregexpr('(?<=<option value=")[^"]+', ui, perl = TRUE))[[1]])
    pick <- function(p) grep(p, opts, value = TRUE)[1]
    session$setInputs(syn_drug_a = pick("^Drug_A Mid [0-9]"), syn_drug_b = pick("^Drug_B"),
                      syn_combo = pick("\\+"), syn_control = pick("^Vehicle"),
                      syn_time_point = "14", syn_endpoint_method = "model",
                      syn_n_boot = 0, run_syn_analysis = 1)
    # "The following specified groups do not exist" before.
    expect_false(isTRUE(rv$syn_results$error))
    expect_equal(rv$syn_plot_inputs$control, "Vehicle 0")
  }))
})

# ---- Warnings reach the user (R20.D34) ---------------------------------------------

test_that("R20.D34: safe_analysis() keeps the backend's warnings and the status shows them", {
  r <- g("safe_analysis")({
    warning("volumes look like cm3")
    list(ok = TRUE)
  }, module = "test")
  expect_equal(r$analysis_warnings, "volumes look like cm3")
  s <- html_of(g("completion_status")(list(r)))
  expect_match(s, "Analysis completed with 1 warning")
  expect_match(s, "volumes look like cm3")
  expect_match(html_of(g("completion_status")(list(list(ok = TRUE)))), "completed successfully")
})
