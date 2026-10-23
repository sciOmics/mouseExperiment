# Scripts that start from the uploaded file (CODE_REVIEW.md R20-AC). With
# "Include my original file" ticked, the ZIP holds the file and the script
# rebuilds the analysed data from it with the package's preparation functions
# (R20-AB), checks the two agree, and runs the analysis on the rebuilt data.
# Each test runs a tab on an upload read as the dashboard reads one, downloads
# the script, runs it, and requires the dashboard's tables.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

upload_rv <- function(file, ...) {
  path <- system.file("sample_data", file, package = "mouseExperiment")
  shiny::reactiveValues(
    uploaded_data = read_study_csv(path), source_file = g("source_file_record")(path),
    data_valid = TRUE, session_id = NULL, settings = list(...))
}
features_rv <- function(...) upload_rv("all_features_synthetic_data.csv",
  col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
  volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
  col_weight = "Weight", col_dose = "Dose", col_survival_event = "Survival_Censor",
  control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse", ...)
combo_rv <- function() upload_rv("combo_treatment_synthetic_data.csv",
  col_id = "ID", col_treatment = "Treatment", col_cage = "Cage", col_dose = "",
  time_method = "from_date", col_date_for_day = "Date", reference_date = "2025-03-24",
  date_order = "auto", volume_method = "calculate", col_length = "Length", col_width = "Width",
  volume_formula = "modified_ellipsoid", control_group = "Control", volume_units = "mm3",
  randomisation_unit = "mouse")
# Exclude one animal (by its arm label, ID and cage) and one cage.
exclude_some <- function(rv) {
  d <- shiny::isolate(rv$uploaded_data)
  shiny::isolate({
    rv$exclude_ids   <- animal_key(paste(d$Treatment[20], d$Dose[20]), d$ID[20], d$Cage[20])
    rv$exclude_cages <- d$Cage[nrow(d)]
  })
  invisible(rv)
}

# Tick only the R script and "Include my original file", download, unzip, run.
# `edit` changes the unzipped original file first.
run_from_upload <- function(output, session, prefix, edit = NULL) {
  groups <- c("tables", "plots", "diagnostics", "data", "report")
  ins <- stats::setNames(rep(list(character(0)), length(groups)),
                         paste0(prefix, "_dl_group_", groups))
  ins[[paste0(prefix, "_dl_group_code")]] <- "script"
  ins[[paste0(prefix, "_dl_include_upload")]] <- TRUE
  do.call(session$setInputs, ins)
  z <- output[[paste0(prefix, "_dl_zip")]]
  dir <- tempfile("upload_")
  utils::unzip(z, exdir = dir)
  if (is.function(edit)) edit(file.path(dir, "data", "original_upload.csv"))
  env <- new.env(parent = globalenv())
  old <- setwd(dir)
  on.exit(setwd(old))
  err <- tryCatch({
    suppressMessages(suppressWarnings(utils::capture.output(
      sys.source("analysis.R", envir = env, keep.source = FALSE))))
    NULL
  }, error = function(e) conditionMessage(e))
  list(env = env, error = err, files = utils::unzip(z, list = TRUE)$Name,
       script = readLines(file.path(dir, "analysis.R")))
}
same_table <- function(a, b, label) {
  expect_false(is.null(a), info = label)
  expect_equal(as.data.frame(a), as.data.frame(b), tolerance = 1e-8, ignore_attr = TRUE,
               info = label)
}
has_line <- function(script, text) any(grepl(text, script, fixed = TRUE))

test_that("R20-AC: Tumor Growth's script rebuilds the data from the upload, exclusions and arm labels included", {
  rv <- exclude_some(features_rv())
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference",
                      tg_p_adjust = "holm")
    session$setInputs(run_tg_analysis = 1)
    res <- rv$tg_results
    expect_false(isTRUE(res$error))
    run <- run_from_upload(output, session, "tg")
    expect_null(run$error)
    expect_true("data/original_upload.csv" %in% run$files)
    expect_true(has_line(run$script, 'prepared <- read_study_csv("data/original_upload.csv")'))
    expect_true(has_line(run$script, "prepared <- apply_exclusions("))
    expect_true(has_line(run$script, "animals = data.frame(arm = "))
    expect_true(has_line(run$script, "prepared <- label_arms("))
    expect_true(has_line(run$script, "df <- prepared"))
    for (f in c("anova", "pairwise_comparisons", "treatment_effects"))
      same_table(run$env$result[[f]], res[[f]], f)

    # A file that does not give the analysed data stops the script.
    bad <- run_from_upload(output, session, "tg", edit = function(path) {
      x <- readLines(path); x[5] <- sub(",[0-9.]+,", ",999,", x[5]); writeLines(x, path)
    })
    expect_match(bad$error, "differ from data/analysed_data.rds", fixed = TRUE)
  }))
})

test_that("R20-AC: without the tick, the ZIP holds no original file and the script no rebuild", {
  rv <- features_rv()
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference",
                      tg_p_adjust = "holm")
    session$setInputs(run_tg_analysis = 1)
    session$setInputs(tg_dl_group_tables = character(0), tg_dl_group_plots = character(0),
                      tg_dl_group_diagnostics = character(0), tg_dl_group_data = character(0),
                      tg_dl_group_report = character(0), tg_dl_group_code = "script",
                      tg_dl_include_upload = FALSE)
    files <- utils::unzip(output$tg_dl_zip, list = TRUE)$Name
    expect_false(any(grepl("original_upload", files, fixed = TRUE)))
  }))
})

test_that("R20-AC: Survival's script rebuilds the per-animal events from the upload", {
  rv <- exclude_some(features_rv())
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = rv), {
    session$setInputs(surv_event_source = "volume", surv_vol_threshold = 1500)
    session$setInputs(run_surv_analysis = 1)
    res <- rv$surv_results
    expect_false(isTRUE(res$error))
    run <- run_from_upload(output, session, "surv")
    expect_null(run$error)
    expect_true(has_line(run$script, "prepared <- survival_events("))
    expect_true(has_line(run$script, "volume_limit = 1500"))
    same_table(run$env$result$results, res$results, "results")
  }))
})

test_that("R20-AC: Toxicity's script rebuilds the data, with the arms labelled by dose", {
  rv <- features_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("toxicity_server"), args = list(rv = rv), {
    session$setInputs(tox_model_type = "bw", tox_adjust_tumor = TRUE, tox_tumor_density = 1,
                      tox_covariates = "volume", tox_estimation = "REML", tox_threshold = 20,
                      tox_baseline_day = "first", run_tox_analysis = 1)
    res <- rv$tox_results
    skip_if(is.null(res$mixed_model), "the body-weight model did not finish here")
    run <- run_from_upload(output, session, "tox")
    expect_null(run$error)
    expect_true(has_line(run$script, 'name = "TreatmentGroup"'))
    same_table(run$env$bw$fixed_effects, res$mixed_model$fixed_effects, "fixed effects")
    same_table(run$env$wl$event_data, res$weight_loss$event_data, "weight-loss events")
  }))
})

test_that("R20-AC: Drug Synergy's script rebuilds days from dates and volumes from calipers", {
  rv <- combo_rv()
  suppressWarnings(shiny::testServer(g("synergy_server"), args = list(rv = rv), {
    session$setInputs(syn_drug_a = "aPD1", syn_drug_b = "HDACi", syn_combo = "HDACi + PD1",
                      syn_time_point = "28", syn_endpoint_method = "model", syn_n_boot = 50,
                      run_syn_analysis = 1)
    res <- rv$syn_results
    expect_false(isTRUE(res$error))
    run <- run_from_upload(output, session, "syn")
    expect_null(run$error)
    expect_true(has_line(run$script, "prepared <- study_days("))
    expect_true(has_line(run$script, "prepared <- calculate_volume("))
    expect_true(has_line(run$script, 'prepared[["Day"]] <- suppressWarnings(as.numeric(prepared[["Day"]]))'))
    same_table(run$env$result$synergy_metrics, res$synergy_metrics, "metrics")
    same_table(run$env$result$synergy_ci, res$synergy_ci, "intervals")
  }))
})

test_that("R20-AC: Dose Response's script rebuilds the data, the dose made numeric", {
  rv <- features_rv()
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("dose_response_server"), args = list(rv = rv), {
    session$setInputs(dr_treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"),
                      dr_time_point = "default", run_dr_analysis = 1)
    res <- rv$dr_results
    expect_false(isTRUE(res$error))
    run <- run_from_upload(output, session, "dr")
    expect_null(run$error)
    expect_true(has_line(run$script, 'prepared[["Dose"]] <- suppressWarnings(as.numeric(prepared[["Dose"]]))'))
    same_table(run$env$result$tgi_table, res$tgi_table, "TGI")
  }))
})

test_that("R20-AC: Power has no data, so it offers no original file", {
  ui <- paste(as.character(g("powerAnalysisUI")("pw")), collapse = " ")
  expect_false(grepl("pw_dl_include_upload", ui, fixed = TRUE))
  ui_tg <- paste(as.character(g("tumorGrowthUI")("tg")), collapse = " ")
  expect_true(grepl("tg-tg_dl_include_upload", ui_tg, fixed = TRUE))
})
