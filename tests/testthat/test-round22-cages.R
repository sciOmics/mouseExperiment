# R22.3: cages numbered within each arm (CODE_REVIEW.md R22-C, 0.41.0). With
# C1 and C2 in every arm, each label was read as one cage holding every arm:
# cage became a fixed effect estimated across arms, and the p-values were too
# small. label_cages() and Study Design's "Numbered within each arm" make the
# labels unique; the analyses warn when the arms reuse one set of labels.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)
html_of <- function(x) {
  if (is.null(x)) return("")
  if (is.list(x) && !is.null(x$html)) return(paste(as.character(x$html), collapse = " "))
  paste(as.character(x), collapse = " ")
}

# master_synthetic_data has two cages per arm, labelled across the study
# (C01-C12). Here each arm's are C1 and C2.
within_arm <- function(d = master_synthetic_data, arm = d$Treatment) {
  d$Cage <- stats::ave(as.character(d$Cage), arm, FUN = function(x) paste0("C", match(x, unique(x))))
  d
}
warnings_of <- function(expr) {
  w <- character()
  withCallingHandlers(suppressMessages(expr), warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  w
}
reuse_warning <- "cage labels come from one set"

test_that("R22.3: label_cages() puts each cage's arm in its label", {
  d <- data.frame(ID = 1:6, Treatment = rep(c("Vehicle", "Drug A"), each = 3),
                  Dose = c(0, 0, 0, 10, 10, 20), Cage = c("C1", "C1", NA, "C1", "C1", "C1"))
  expect_identical(label_cages(d)$Cage,
                   c("Vehicle C1", "Vehicle C1", NA, "Drug A C1", "Drug A C1", "Drug A C1"))
  # With a dose, the arm is "Treatment Dose", as label_arms() writes it.
  expect_identical(label_cages(d, dose_column = "Dose")$Cage[4:6],
                   c("Drug A 10 C1", "Drug A 10 C1", "Drug A 20 C1"))
  expect_error(label_cages(d, cage_column = "Box"), "'Box' is not in the data")
})

test_that("R22.3: the reuse pattern is every arm's labels coming from one set", {
  f <- g("reused_cage_labels")
  arm <- rep(c("V", "A", "B"), each = 4)
  expect_identical(f(rep(c("C1", "C2"), 6), arm), c("C1", "C2"))
  # Vehicle in three cages, the others in two: still numbered within arm.
  expect_identical(f(c("C1", "C2", "C3", "C3", rep(c("C1", "C2"), 4)), arm), c("C1", "C2", "C3"))
  expect_identical(f(c("C10", "C2", "C1", "C1", "C2", "C10"), rep(c("V", "A"), each = 3)),
                   c("C1", "C2", "C10"))
  # Unique labels, and cages that each hold a different pair of arms, are not.
  expect_identical(f(paste0("C", rep(1:6, each = 2)), arm), character(0))
  expect_identical(f(c("C1", "C1", "C2", "C2", "C1", "C1", "C3", "C3", "C2", "C2", "C3", "C3"), arm),
                   character(0))
  expect_identical(f(c("C1", "C2"), c("V", "V")), character(0))   # one arm
})

test_that("R22.3: Tumor Growth warns on reused labels, and label_cages() makes cage a random effect", {
  w <- within_arm()
  fit <- function(d) suppressWarnings(suppressMessages(tumor_growth_statistics(
    d, reference_group = "Vehicle", cage_column = "Cage", plots = FALSE)))
  expect_true(any(grepl(reuse_warning, warnings_of(tumor_growth_statistics(
    w, reference_group = "Vehicle", cage_column = "Cage", plots = FALSE)), fixed = TRUE)))
  reused <- fit(w)
  fixed  <- fit(label_cages(w))
  expect_identical(reused$cage_analysis$structure$structure, "crossed")
  expect_identical(fixed$cage_analysis$structure$structure, "nested_replicated")
  # With the labels unique the model is the one the original C01-C12 give
  # (to the optimiser's precision: the cage levels are in another order).
  orig <- fit(master_synthetic_data)
  for (col in c("estimate", "SE"))
    expect_equal(fixed$pairwise_comparisons[[col]], orig$pairwise_comparisons[[col]],
                 tolerance = 1e-4, info = col)
  # Read as shared cages, every p-value was smaller.
  expect_true(all(reused$pairwise_comparisons$p.value < fixed$pairwise_comparisons$p.value))
  expect_false(any(grepl(reuse_warning, warnings_of(tumor_growth_statistics(
    label_cages(w), reference_group = "Vehicle", cage_column = "Cage", plots = FALSE)), fixed = TRUE)))
})

test_that("R22.3: Survival and the body-weight model warn on reused labels too", {
  w <- within_arm()
  ev <- survival_events(w, cage_column = "Cage", volume_column = "Volume", volume_limit = 1500)
  expect_true(any(grepl(reuse_warning, warnings_of(survival_statistics(
    ev, time_column = "Time", event_column = "Event", cage_column = "Cage",
    reference_group = "Vehicle", verbose = FALSE)), fixed = TRUE)))
  expect_true(any(grepl(reuse_warning, warnings_of(analyze_body_weight(
    w, cage_column = "Cage", reference_group = "Vehicle", volume_units = "mm3")), fixed = TRUE)))
  expect_false(any(grepl(reuse_warning, warnings_of(survival_statistics(
    survival_events(label_cages(w), cage_column = "Cage", volume_column = "Volume", volume_limit = 1500),
    time_column = "Time", event_column = "Event", cage_column = "Cage",
    reference_group = "Vehicle", verbose = FALSE)), fixed = TRUE)))
})

test_that("R22.3: the cage-level randomisation test names reused labels when it stops", {
  expect_error(trajectory_permutation_test(within_arm(), cage_column = "Cage",
                                           spec = perm_spec(unit = "cage", n_perm = 19, seed = 1)),
               "label_cages()", fixed = TRUE)
})

test_that("R22.3: apply_exclusions() removes one arm's cage, given arm and cage", {
  w <- within_arm()
  out <- apply_exclusions(w, cages = data.frame(arm = "Vehicle", cage = "C1"), cage_column = "Cage")
  expect_identical(nrow(w) - nrow(out), sum(w$Treatment == "Vehicle" & w$Cage == "C1"))
  expect_true(all(out$Cage[out$Treatment == "Drug_B"] %in% c("C1", "C2")))
  expect_true("C1" %in% out$Cage[out$Treatment == "Drug_B"])
  # A label alone is still the whole cage, in every arm.
  expect_false("C1" %in% apply_exclusions(w, cages = "C1", cage_column = "Cage")$Cage)
  expect_error(apply_exclusions(w, cages = data.frame(arm = "Vehicle"), cage_column = "Cage"),
               "columns arm and cage")
})

# ---- the dashboard ---------------------------------------------------------------

within_settings <- function(...) utils::modifyList(list(
  col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
  volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
  col_dose = "", control_group = "Vehicle", volume_units = "mm3",
  randomisation_unit = "mouse", cage_labels = "arm"), list(...))

test_that("R22.3: the preparation labels the cages after the exclusions, and records the step", {
  w <- within_arm()
  rv <- shiny::reactiveValues(uploaded_data = w, data_valid = TRUE, settings = within_settings(),
                              exclude_ids = animal_key("Vehicle", w$ID[1], w$Cage[1]),
                              exclude_cages = animal_key("Drug_B", "C2"))
  pu <- shiny::isolate(g("prepare_upload")(rv))
  fns <- vapply(pu$steps, `[[`, "", "fun")
  expect_identical(fns[1:2], c("apply_exclusions", "label_cages"))
  expect_true(all(startsWith(pu$data$Cage, paste0(pu$data$Treatment, " C"))))
  expect_false(w$ID[1] %in% pu$data$ID[pu$data$Treatment == "Vehicle"])
  expect_false("Drug_B C2" %in% pu$data$Cage)
  expect_true("Drug_A Low C2" %in% pu$data$Cage)
  # Labels unique across the study: no relabelling step.
  shiny::isolate(rv$settings <- within_settings(cage_labels = "study"))
  pu <- shiny::isolate(g("prepare_upload")(rv))
  expect_false("label_cages" %in% vapply(pu$steps, `[[`, "", "fun"))
  expect_true(all(pu$data$Cage %in% c("C1", "C2")))
})

test_that("R22.3: the Subset Filter offers each arm's cage, and keeps exclusions across the readings", {
  w <- within_arm()
  f <- g("cage_filter_choices")
  arm <- f(w, within_settings(), stored = c("C2", animal_key("Vehicle", "C1")))
  expect_identical(names(arm$choices)[1:2], c("Drug_A High C1", "Drug_A High C2"))
  # A whole cage C2 shows as each arm's C2; Vehicle's C1 as itself.
  expect_setequal(arm$selected, c(animal_key(unique(w$Treatment), "C2"), animal_key("Vehicle", "C1")))
  study <- f(w, within_settings(cage_labels = "study"), stored = c("C2", animal_key("Vehicle", "C1")))
  expect_identical(unname(study$choices), c("C1", "C2", animal_key("Vehicle", "C1")))
  expect_identical(names(study$choices)[3], "Vehicle C1")
  expect_setequal(study$selected, c("C2", animal_key("Vehicle", "C1")))
  # Either way, what is excluded is the same rows.
  ce <- g("cage_exclusions")(w, within_settings(), c("C2", animal_key("Vehicle", "C1")))
  out <- apply_exclusions(w, cages = ce, cage_column = "Cage")
  expect_identical(nrow(out), sum(w$Cage == "C1" & w$Treatment != "Vehicle"))
  expect_identical(g("cage_exclusions")(w, within_settings(), "C2"), "C2")
})

test_that("R22.3: Study Design asks how cages are labelled, guessing from the labels", {
  ui <- html_of(g("cage_labels_ui")(c("C1", "C2")))
  expect_match(ui, 'value="arm"[^>]*checked')
  expect_match(html_of(g("cage_labels_ui")(character(0))), 'value="study"[^>]*checked')
  expect_match(html_of(g("cage_labels_ui")(c("C1", "C2"), selected = "study")), 'value="study"[^>]*checked')
  expect_match(html_of(g("cage_labels_note")(c("C1", "C2"), "study")), "text-warning")
  w <- within_arm()
  shiny::testServer(g("data_upload_server"), args = list(rv = shiny::reactiveValues()), {
    rv$uploaded_data <- w
    rv$data_valid <- TRUE
    session$setInputs(col_id = "ID", col_treatment = "Treatment", col_cage = "Cage",
                      col_day = "Day", time_method = "direct_day", volume_method = "direct",
                      col_volume = "Volume")
    expect_match(html_of(output$cage_labels_ui), 'value="arm"[^>]*checked')
    expect_match(html_of(output$cage_labels_ui), "Cage labels are", fixed = TRUE)
    session$setInputs(cage_labels = "arm")
    expect_identical(rv$settings$cage_labels, "arm")
    expect_match(html_of(output$cage_labels_status), "Every arm uses the labels C1, C2", fixed = TRUE)
    session$setInputs(cage_labels = "study")
    expect_match(html_of(output$cage_labels_status), "text-warning", fixed = TRUE)
    # Without a cage column the question is not asked, and nothing is relabelled.
    session$setInputs(col_cage = "")
    expect_identical(rv$settings$cage_labels, "study")
  })
})

test_that("R22.3: the R script from the original file relabels the cages as the dashboard did", {
  path <- file.path(tempfile("cages_"), "study.csv")
  dir.create(dirname(path))
  utils::write.csv(within_arm()[, c("ID", "Day", "Volume", "Treatment", "Cage")], path, row.names = FALSE)
  rv <- shiny::reactiveValues(uploaded_data = read_study_csv(path),
                              source_file = g("source_file_record")(path), data_valid = TRUE,
                              session_id = NULL, settings = within_settings(),
                              exclude_cages = animal_key("Drug_B", "C2"))
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference", tg_p_adjust = "holm")
    session$setInputs(run_tg_analysis = 1)
    res <- rv$tg_results
    expect_false(isTRUE(res$error))
    expect_identical(res$cage_analysis$structure$structure, "nested_partial")
    groups <- c("tables", "plots", "diagnostics", "data", "report")
    ins <- stats::setNames(rep(list(character(0)), length(groups)), paste0("tg_dl_group_", groups))
    do.call(session$setInputs, c(ins, list(tg_dl_group_code = "script", tg_dl_include_upload = TRUE)))
    dir <- tempfile("script_")
    utils::unzip(output$tg_dl_zip, exdir = dir)
    script <- readLines(file.path(dir, "analysis.R"))
    expect_true(any(grepl("prepared <- label_cages(", script, fixed = TRUE)))
    expect_true(any(grepl('cages = data.frame(arm = "Drug_B", cage = "C2")', script, fixed = TRUE)))
    env <- new.env(parent = globalenv())
    old <- setwd(dir); on.exit(setwd(old))
    suppressMessages(suppressWarnings(utils::capture.output(
      sys.source("analysis.R", envir = env, keep.source = FALSE))))
    expect_equal(as.data.frame(env$result$pairwise_comparisons), as.data.frame(res$pairwise_comparisons),
                 tolerance = 1e-8, ignore_attr = TRUE)
  }))
})
