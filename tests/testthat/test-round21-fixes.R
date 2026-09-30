# Regression tests for CODE_REVIEW.md Round 21 fixes (IDs R20.D*), labelled by
# finding ID. They drive the code path through testServer() or call the helper;
# none of them only reads the source.

skip_if_not_installed("shiny", minimum_version = "1.5.0")
skip_if_not_installed("mouseExperiment", minimum_version = "0.24.0")

.pkg <- tryCatch(getNamespace("mouseExperimentDashboard"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperimentDashboard namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

render_text <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}

master_rv <- function(...) {
  path <- system.file("sample_data", "master_synthetic_data.csv",
                      package = "mouseExperimentDashboard")
  shiny::reactiveValues(
    uploaded_data = utils::read.csv(path, stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day",
                    col_volume = "Volume", volume_method = "precalculated",
                    col_treatment = "Treatment", col_cage = "Cage",
                    col_weight = "Weight", col_dose = "Dose",
                    col_survival_event = "Survival_Censor"),
    ...)
}

# Four arms with log-linear growth; `rates` are per-day growth rates.
synergy_df <- function(rates, n = 8, days = seq(0, 28, by = 3.5), seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(names(rates), function(arm) {
    do.call(rbind, lapply(seq_len(n), function(i) {
      b0 <- stats::rnorm(1, log(100), 0.1)
      r  <- rates[[arm]] + stats::rnorm(1, 0, 0.005)
      data.frame(ID = paste(arm, i), Treatment = arm, Day = days,
                 Volume = exp(b0 + r * days + stats::rnorm(length(days), 0, 0.08)))
    }))
  }))
}

# ---- R20.D3 -----------------------------------------------------------------

test_that("R20.D3: a session without the cookie gets a key no browser can present", {
  k <- g("make_local_sid")()
  expect_match(k, "^local-[a-f0-9]{64}$")
  # get_sid() takes exactly 64 hex characters from the Cookie header.
  expect_null(g("get_sid")(list(request = list(HTTP_COOKIE = paste0("mE_session=", k)))))
  expect_false(exists(".legacy_inline_tg_bayes",   envir = .pkg, inherits = FALSE))
  expect_false(exists(".legacy_inline_surv_bayes", envir = .pkg, inherits = FALSE))
})

test_that("R20.D3: a Bayesian launch without a cookie runs through the registry", {
  g(".reset_job_registry")()
  on.exit(g(".reset_job_registry")(), add = TRUE)
  testthat::local_mocked_bindings(
    .tg_bayes_promise = function(bayes_params) promises::promise_resolve(list(ok = TRUE)),
    .package = "mouseExperimentDashboard")
  srv <- function(id, rv) shiny::moduleServer(id, function(input, output, session) {
    g("launch_tg_bayes_job")(g("SID_NONE"), list(data = NULL), input, output, session, rv)
  })
  rv <- shiny::reactiveValues(session_id = g("SID_NONE"))
  shiny::testServer(srv, args = list(rv = rv), {
    sid <- shiny::isolate(rv$session_id)
    expect_match(sid, "^local-")
    jobs <- g("list_jobs")(sid)
    expect_length(jobs, 1L)
    expect_equal(jobs[[1L]]$kind, "bayes_tg")
  })
})

# ---- R20.D4 -----------------------------------------------------------------

test_that("R20.D4: an unusable reference date is reported, not a missing symbol", {
  f  <- g("detect_date_format")
  df <- data.frame(Date = c("03/24/2025", "03/28/2025"))
  expect_error(f(df, "Date", "24/03/2025"), "not a recognised date")
  expect_error(f(df, "Date", ""), "No reference date is set")
  ok <- f(df, "Date", "2025-03-24")
  expect_equal(ok$date_format, "%m/%d/%Y")
  expect_equal(ok$start_date, "03/24/2025")
})

test_that("R20.D4: Tumor Growth reports a bad reference date and stops running", {
  df <- data.frame(
    ID = rep(1:6, each = 3), Treatment = rep(c("Control", "Drug"), each = 9),
    Cage = "C1", Date = rep(c("03/24/2025", "03/28/2025", "04/01/2025"), 6),
    Volume = c(100, 150, 220, 95, 140, 210, 105, 160, 230,
               100, 120, 140, 98, 115, 138, 102, 125, 150))
  rv <- shiny::reactiveValues(
    uploaded_data = df, data_valid = TRUE,
    settings = list(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage",
                    time_method = "from_date", col_date_for_day = "Date",
                    reference_date = "24/03/2025", volume_method = "direct",
                    col_volume = "Volume"))
  shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_reference = "Control",
                      run_tg_analysis = 1)
    expect_false(isTRUE(rv$tg_running))
    expect_match(render_text(output$tg_status_msg), "not a recognised date")
    expect_null(rv$tg_results)
  })
})

# ---- R20.D19 ----------------------------------------------------------------

test_that("R20.D19: the survival forest-plot download produces a PNG", {
  expect_true(exists("build_surv_forest_plot", envir = .pkg, inherits = FALSE))
  # The master demo has arms with few events, so coxph warns about convergence.
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = master_rv()), {
    shiny::isolate(g("ensure_preprocessed")(rv))
    session$setInputs(run_surv_analysis = 1)
    expect_false(isTRUE(rv$surv_results$error))
    path <- output$download_surv_forest_plot          # HTTP 500 before
    expect_true(file.exists(path) && file.size(path) > 0)
    expect_identical(readBin(path, "raw", 4L), as.raw(c(0x89, 0x50, 0x4e, 0x47)))
  }))
})

# ---- R20.D21 ----------------------------------------------------------------

test_that("R20.D21: info alerts use a known icon", {
  # The unknown-icon notice is a message, not a warning.
  expect_no_message(x <- g("status_alert")("Running", "info"))
  expect_match(as.character(x), "fa-spin", fixed = TRUE)
})

# ---- R20.D39 ----------------------------------------------------------------

test_that("R20.D39: AUC is not offered Dunnett or Tukey, and a stale Dunnett still runs", {
  ch <- g("tg_p_adjust_choices")
  expect_false(any(c("dunnett", "tukey") %in% ch("auc")))
  expect_true(all(c("dunnett", "tukey") %in% ch("lme4")))
  # The AUC path's equal-variance test warns on small cells.
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = master_rv()), {
    # With a dose column mapped, arms are labelled "Treatment Dose".
    session$setInputs(tg_model_type = "auc", tg_reference = "Vehicle 0",
                      tg_p_adjust = "dunnett", run_tg_analysis = 1)
    expect_false(isTRUE(rv$tg_results$error))          # errored before
    expect_equal(rv$tg_results$p_adjust_method_used, "holm")
  }))
})

# ---- R20.D17 ----------------------------------------------------------------

test_that("R20.D17: the analytic Power tab plans k - 1 corrected comparisons and attrition", {
  shiny::testServer(g("power_server"), args = list(rv = shiny::reactiveValues()), {
    session$setInputs(pw_mode = "apriori_analytic", pw_ap_es_mode = "cohens_d",
                      pw_ap_cohens_d = "1", pw_ap_n_groups = "4",
                      pw_ap_submode = "find_n", pw_alpha = "0.05", pw_power = "0.80",
                      pw_ap_p_adjust = "bonferroni", pw_ap_n_comparisons = 3,
                      pw_ap_dropout = 20, run_pw_analysis = 1)
    st <- rv$ap_results$scenario_table
    expect_equal(st$Alpha_Per_Comparison, 0.05 / 3)
    expect_equal(st$Required_N, 23)                     # the tab said 7 (R20.7 + R20.D17)
    expect_equal(st$Enroll_N, ceiling(23 / 0.8))
    txt <- render_text(output$pw_ap_summary_ui)
    expect_match(txt, "Bonferroni over 3")
    expect_match(txt, "20 %")
  })
})

# ---- R20.8 (dashboard side) -------------------------------------------------

test_that("R20.8: the synergy summary says 'not evaluable' rather than printing NA", {
  harm <- synergy_df(c(Control = 0.10, DrugA = 0.13, DrugB = 0.10, Combo = 0.10))
  res <- suppressWarnings(mouseExperiment::analyze_drug_synergy(
    harm, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE, n_boot = 200, boot_seed = 1))
  rv <- shiny::reactiveValues(syn_results = res, syn_over_time = FALSE,
                              uploaded_data = harm, data_valid = TRUE)
  # The bar plot drops the NA "Bliss Expected" bar, with a ggplot warning.
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = rv), {
    txt <- render_text(output$syn_summary)
    expect_match(txt, "Not evaluable")
    expect_false(grepl("Synergistic:", txt, fixed = TRUE))  # printed "TRUE" before
    expect_equal(render_text(output$syn_ci_verdict_note), "") # green note before
  }))

  ot <- suppressWarnings(suppressMessages(mouseExperiment::analyze_drug_synergy_over_time(
    harm, drug_a_name = "DrugA", drug_b_name = "DrugB", combo_name = "Combo",
    control_name = "Control", verbose = FALSE)))
  rv2 <- shiny::reactiveValues(syn_results = ot, syn_over_time = TRUE,
                               uploaded_data = harm, data_valid = TRUE)
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = rv2), {
    rv$syn_over_time <- TRUE    # the module initialises it to FALSE
    session$flushReact()
    txt <- render_text(output$syn_summary)
    expect_match(txt, "Evaluable timepoints: 0 of 9")
    expect_match(txt, "No Bliss peak")
  }))
})
