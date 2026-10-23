# Regression tests for step 5 of the Round 20/21 implementation (CODE_REVIEW.md
# R20-D.M): the two-axis therapeutic window, the interval-based synergy verdict,
# one agent's dose series with its trend test and EC50 interval, and the
# optional removal-reason column. Each test fails on v26.09.004 run against
# backend 0.27.0.

skip_if_not_installed("shiny", minimum_version = "1.5.0")
skip_if_not_installed("mouseExperiment", minimum_version = "0.27.0")

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

master_rv <- function(...) {
  shiny::reactiveValues(
    uploaded_data = demo_csv("all_features_synthetic_data.csv"), data_valid = TRUE,
    session_id = NULL,
    settings = utils::modifyList(list(
      col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
      volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
      col_weight = "Weight", col_dose = "Dose", col_survival_event = "Survival_Censor",
      control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse"),
      list(...)))
}

combo_rv <- function() {
  shiny::reactiveValues(
    uploaded_data = demo_csv("combo_treatment_synthetic_data.csv"), data_valid = TRUE,
    session_id = NULL,
    settings = list(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", time_method = "from_date",
                    col_date_for_day = "Date", reference_date = "2025-03-24",
                    date_order = "auto", volume_method = "calculate",
                    col_length = "Length", col_width = "Width",
                    volume_formula = "modified_ellipsoid",
                    control_group = "Control", volume_units = "mm3",
                    randomisation_unit = "mouse"))
}

# Two arms of 8 animals, days 0-20, with a removal reason on each animal's last
# row. Vehicle tumours reach 2,000 mm3 near day 18 ("Tumour volume limit").
# Drug animals 1-2 die on day 8 with small tumours ("Found dead"); animal 3 is
# removed on day 12 for "Body condition", having lost about 14 %; the rest
# reach the end ("End of study").
reason_df <- function() {
  set.seed(21)
  days <- seq(0, 20, 2)
  rows <- list()
  for (arm in c("Vehicle", "Drug")) for (i in 1:8) {
    rate <- (if (arm == "Vehicle") 0.15 else 0.05) + stats::rnorm(1, 0, 0.005)
    v <- 150 * exp(rate * days)
    w <- 22 + stats::rnorm(length(days), 0, 0.1)
    end <- 20; why <- "End of study"
    if (arm == "Vehicle" && any(v >= 2000)) {
      end <- days[which(v >= 2000)[1]]; why <- "Tumour volume limit"
    }
    if (arm == "Drug" && i %in% 1:2) { end <- 8;  why <- "Found dead" }
    if (arm == "Drug" && i == 3)     { end <- 12; why <- "Body condition"
                                       w <- w * (1 - 0.012 * days) }
    keep <- days <= end
    rr <- rep("", sum(keep)); rr[sum(keep)] <- why
    rows[[length(rows) + 1L]] <- data.frame(
      ID = i, Treatment = arm, Cage = paste0(arm, "-C1"), Day = days[keep],
      Volume = v[keep], Weight = w[keep], Removal_Reason = rr,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

reason_map <- list("Body condition" = "weight", "End of study" = "planned",
                   "Found dead" = "death", "Tumour volume limit" = "limit")

reason_rv <- function(mapped = TRUE) {
  shiny::reactiveValues(
    uploaded_data = reason_df(), data_valid = TRUE, session_id = NULL,
    settings = list(
      col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
      volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
      col_weight = "Weight", control_group = "Vehicle", volume_units = "mm3",
      randomisation_unit = "mouse",
      col_removal_reason = if (mapped) "Removal_Reason" else "",
      removal_reason_map = if (mapped) reason_map else list()))
}

# ---- R20-K: the therapeutic window, two axes ------------------------------------

test_that("R20-K: Toxicity sets TGI beside worst weight loss, with intervals and a flag", {
  rv <- master_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "et", tox_adjust_tumor = TRUE,
                      tox_tumor_density = 1, tox_endpoint_method = "model",
                      tox_n_boot = 200, tox_tol_threshold = 12, run_tox_analysis = 1)
    tw <- rv$tox_results$twm
    expect_false(isTRUE(tw$error))
    expect_equal(tw$tolerability_threshold, 12)                  # the sidebar's threshold
    expect_true(all(c("Worst_Loss", "Tolerability") %in% names(tw$window_table)))
    txt <- paste(output$tox_summary, collapse = " ")
    expect_match(txt, "Tolerability threshold: 12% weight loss")
    expect_match(txt, "Vehicle:  TGI \\(reference\\)")
    expect_false(grepl("TWM", txt, fixed = TRUE))
    expect_match(html_of(output$tox_twm_tolerability_note), "at or above 12%")
    expect_true(nzchar(output$tox_et_twm_plot$src))
  }))
})

test_that("R20-K: the displayed window table reads as two axes, not a ranking", {
  w <- data.frame(Treatment = c("Vehicle", "A"), N_Animals = c(8L, 8L),
                  TGI = c(0, 61.23), TGI_Lower = c(0, 50.01), TGI_Upper = c(0, 70.4),
                  Worst_Loss = c(3.2, 21.7), Worst_Loss_Lower = c(2.1, 19.2),
                  Worst_Loss_Upper = c(4.4, NA), N_Over_Threshold = c(0L, 5L),
                  Tolerability = c("Tolerated", NA), stringsAsFactors = FALSE)
  d <- g("window_display_table")(w, "Vehicle")
  expect_identical(d$`TGI 95% interval`, c("reference", "[50.0, 70.4]"))
  expect_true(is.na(d$`TGI (%)`[1]))
  expect_identical(d$`Weight loss 95% interval`[2], "—")
  expect_identical(d$`Over threshold`, c("0 of 8", "5 of 8"))
  expect_identical(d$Tolerability, c("Tolerated", "No interval"))
})

# ---- R20.19: the synergy verdict ---------------------------------------------------

test_that("R20.19: the Synergy verdict rests on the interval, and says so without one", {
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = combo_rv()), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "28", syn_endpoint_method = "model",
                      syn_n_boot = 200, run_syn_analysis = 1)
    res <- rv$syn_results
    expect_true(res$overall_assessment %in%
                  c("Synergy", "Antagonism", "Additive (no departure from Bliss detected)"))
    note <- html_of(output$syn_consensus_note)
    expect_match(note, "95% interval for the Bliss excess")        # the rule is shown
    session$setInputs(syn_n_boot = 0, run_syn_analysis = 2)
    expect_match(rv$syn_results$overall_assessment, "no interval")
    expect_match(html_of(output$syn_consensus_note), "alert-info")  # not a green verdict
  }))
})

# ---- R20.D47 / R20.D49: dose-response ----------------------------------------------

test_that("R20.D47: a dose series is suggested from the arm names, or left to the user", {
  m <- demo_csv("all_features_synthetic_data.csv")
  arms <- setdiff(sort(unique(m$Treatment)), "Vehicle")
  expect_setequal(g("guess_dose_series")(m, "Treatment", "Dose", arms),
                  c("Drug_A Low", "Drug_A Mid", "Drug_A High"))
  expect_identical(g("guess_dose_series")(m, "Treatment", "Dose", "Drug X"), "Drug X")
  # Two agents with one arm each: nothing to suggest.
  expect_identical(g("guess_dose_series")(m, "Treatment", "Dose", c("Drug_B", "Other")),
                   character(0))
})

test_that("R20.D47 / R20.D49: Dose Response fits one agent and shows the JT test and the EC50 interval", {
  rv <- master_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    ui <- html_of(output$dr_series_ui)
    expect_match(ui, "Dose series")
    for (a in c("Drug_A Low", "Drug_A Mid", "Drug_A High"))
      expect_match(ui, paste0('<option value="', a, '" selected>'), fixed = TRUE)
    expect_false(grepl('<option value="Drug_B" selected>', ui, fixed = TRUE))
    # Nothing chosen: no analysis.
    session$setInputs(dr_time_point = "default", run_dr_analysis = 1)
    expect_null(rv$dr_results)
    expect_match(html_of(output$dr_status_msg), "one agent")
    session$setInputs(dr_treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
                      run_dr_analysis = 2)
    res <- rv$dr_results
    expect_false(isTRUE(res$error))
    expect_setequal(res$series$treatments, c("Drug_A Low", "Drug_A Mid", "Drug_A High"))
    expect_setequal(unique(res$analysis_data$Treatment),
                    c("Vehicle", "Drug_A Low", "Drug_A Mid", "Drug_A High"))
    txt <- paste(output$dr_summary, collapse = " ")
    expect_match(txt, "Jonckheere-Terpstra, two-sided")
    expect_match(txt, "on log dose")
    expect_match(txt, "outside the tested doses")                   # EC50 2.4 < 5
    expect_match(paste(output$dr_stat_tests, collapse = " "), "JONCKHEERE-TERPSTRA")
    # Two agents on one axis are refused, with the backend's reason.
    session$setInputs(dr_treatments = c("Drug_A Mid", "Drug_A Mid + Drug_B"),
                      run_dr_analysis = 3)
    expect_true(isTRUE(rv$dr_results$error))
    expect_match(rv$dr_results$message, "Dose 15 holds more than one arm")
    expect_match(rv$dr_results$message, "under 'Dose series' in the sidebar", fixed = TRUE)
  }))
})

test_that("R20.D47: with one treatment name at every dose, the dose-0 rows are the control", {
  # One agent named at every dose, its dose-0 animals the control.
  set.seed(4)
  d <- do.call(rbind, lapply(c(0, 10, 30, 100), function(dz) do.call(rbind, lapply(1:6, function(i) {
    days <- seq(0, 21, 3)
    data.frame(ID = paste0(dz, "-", i), Treatment = "Drug X", Dose = dz, Day = days,
               Volume = 100 * exp((0.15 - 0.001 * dz) * days + stats::rnorm(8, 0, 0.05)))
  }))))
  rv <- shiny::reactiveValues(
    uploaded_data = d, data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_dose = "Dose",
                    control_group = "Drug X", volume_units = "mm3",
                    randomisation_unit = "mouse"))
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    expect_match(html_of(output$dr_reference_note), "the rows with dose 0")
    expect_match(html_of(output$dr_series_ui), '<option value="Drug X" selected>', fixed = TRUE)
    session$setInputs(dr_treatments = "Drug X", dr_time_point = "default",
                      run_dr_analysis = 1)
    res <- rv$dr_results
    expect_false(isTRUE(res$error))               # "the control must have dose 0" otherwise
    expect_equal(res$tgi_table$Dose[res$tgi_table$Control], 0)
    expect_equal(res$tgi_table$Dose, c(0, 10, 30, 100))
  }))
})

# ---- R20-D.H: the removal-reason column ---------------------------------------------

test_that("R20-D.H: a removal-reason column is mapped, and each value given a meaning", {
  expect_identical(g("guess_removal_category")(
    c("Found dead", "Moribund", "Body condition score 2", "Weight loss >20%",
      "Tumour volume limit", "Tumor ulcerated", "End of study", "Data cut", "Escaped")),
    c("death", "death", "weight", "weight", "limit", "limit", "planned", "planned", "other"))
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    f <- tempfile(fileext = ".csv")
    utils::write.csv(reason_df(), f, row.names = FALSE)
    session$setInputs(data_file = list(name = "study.csv", datapath = f))
    expect_true(rv$data_valid)
    expect_identical(rv$detected_cols$removal_reason, "Removal_Reason")
    session$setInputs(col_treatment = "Treatment", col_id = "ID",
                      col_removal_reason = "Removal_Reason")
    ui <- html_of(output$removal_reasons_ui)
    for (v in names(reason_map)) expect_match(ui, v, fixed = TRUE)
    expect_equal(rv$settings$removal_reason_map, reason_map)       # guessed from the words
    session$setInputs(reason_cat_4 = "other")                      # "Tumour volume limit"
    expect_equal(rv$settings$removal_reason_map[["Tumour volume limit"]], "other")
    saved <- jsonlite::read_json(output$save_preset, simplifyVector = TRUE)
    expect_equal(saved$settings$col_removal_reason, "Removal_Reason")
    expect_equal(saved$settings$removal_reason_map[["Tumour volume limit"]], "other")
  })
})

test_that("R20.D13: Survival counts deaths as events only when the reasons say who died", {
  # Without a reason column there is no checkbox, and deaths are censored.
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = reason_rv(FALSE)), {
    expect_identical(html_of(output$surv_death_option_ui), "")
    session$setInputs(surv_event_source = "volume", surv_vol_threshold = 2000,
                      run_surv_analysis = 1)
    sd <- rv$surv_results$survival_data
    expect_equal(sum(sd$Event[sd$Treatment == "Drug"]), 0)
    expect_true(is.na(rv$surv_results$death_events_added))
  }))
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = reason_rv(TRUE)), {
    expect_match(html_of(output$surv_death_option_ui), "Count deaths from other causes")
    session$setInputs(surv_event_source = "volume", surv_vol_threshold = 2000,
                      surv_death_as_event = TRUE, run_surv_analysis = 1)
    expect_match(html_of(output$surv_event_assumption), "Deaths from other causes are events")
    sd <- rv$surv_results$survival_data
    expect_equal(sum(sd$Event[sd$Treatment == "Drug"]), 2)        # the two found dead
    expect_equal(sort(sd$Time[sd$Treatment == "Drug" & sd$Event == 1]), c(8, 8))
    expect_equal(rv$surv_results$death_events_added, 2L)
    expect_match(paste(output$surv_summary, collapse = " "),
                 "Deaths from other causes counted as events: 2")
    # Unticked: censored again.
    session$setInputs(surv_death_as_event = FALSE, run_surv_analysis = 2)
    sd <- rv$surv_results$survival_data
    expect_equal(sum(sd$Event[sd$Treatment == "Drug"]), 0)
  }))
})

test_that("R20.23: the weight-loss threshold uses the removal reasons, and says how", {
  run_bw <- function(rv) {
    shiny::isolate(g("ensure_preprocessed")(rv))
    res <- NULL
    suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
      session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE,
                        tox_tumor_density = 1, tox_covariates = "volume",
                        tox_estimation = "REML", tox_threshold = 20,
                        tox_baseline_day = "first", run_tox_analysis = 1)
      res <<- list(wl = rv$tox_results$weight_loss,
                   note = html_of(output$tox_km_competing_note))
    }))
    res
  }
  with_r <- run_bw(reason_rv(TRUE))
  ed <- with_r$wl$event_data
  body <- ed[ed$Treatment == "Drug" & ed$ID == "3", ]
  expect_equal(body$Event, 1L)                           # removed for body condition
  expect_equal(body$Time, 12)
  expect_equal(with_r$wl$n_competing_risk, 10L)          # 8 tumour limit + 2 deaths
  expect_match(with_r$note, "Removal reasons")
  without <- run_bw(reason_rv(FALSE))
  expect_equal(without$wl$n_competing_risk, 0L)
  expect_match(without$note, "No removal-reason column")
})
