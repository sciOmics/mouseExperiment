# Every run of the session is kept, can be opened again, downloaded with
# others, and an identical Bayesian fit is not refitted (CODE_REVIEW.md
# Round 23, R23.15, R23-I).

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

fake_run <- function(transform, data = master_synthetic_data[1:20, ]) {
  rec <- g("record_call")("tumor_growth_statistics", list(df = data, transform = transform))
  list(error = FALSE, model_type_used = "lme4", transform_used = transform, reference_group = "Vehicle",
       anova = data.frame(term = "Day", p = 0.01), call_record = rec)
}
# Stands in for Tumor Growth's Download Results registration.
fake_registry <- function(rv) {
  list(tumor_growth = list(
    title = function() "Tumor Growth",
    items = function() list(g("dl_item")("anova", "tables", "ANOVA", "anova", "table",
                                         function() data.frame(transform = rv$tg_results$transform_used),
                                         "The model's ANOVA")),
    header = function() "mouseExperiment dashboard: Tumor Growth results",
    size = function() list()))
}

test_that("R23.15: a run's fingerprint is its call, not the time it was made", {
  fp <- g("history_fingerprint")
  a <- fake_run("log")$call_record
  b <- a; b$when <- b$when + 3600
  expect_identical(fp(a), fp(b))
  expect_false(identical(fp(a), fp(fake_run("sqrt")$call_record)))
  expect_false(identical(fp(a), fp(fake_run("log", master_synthetic_data[1:21, ])$call_record)))
  expect_null(fp(NULL))
  # A tab that makes two calls (Toxicity) is identified by both.
  expect_type(fp(list(model = a, threshold = b)), "character")
})

test_that("R23.15: runs are saved, identical ones replace the last, and the table lists them", {
  rv <- shiny::reactiveValues(tg_results = NULL, tg_plot_data = NULL)
  shiny::testServer(g("run_history_server"), args = list(rv = rv), {
    h <- session$userData$run_history
    expect_length(h$entries, 0L)
    rv$tg_plot_data <- data.frame(Day = 1:3)
    rv$tg_results <- fake_run("log"); session$flushReact()
    expect_length(h$entries, 1L)
    expect_identical(h$entries[[1]]$model, "LME4 mixed model, log")
    expect_identical(h$entries[[1]]$rows, 20L)
    expect_true(file.exists(h$entries[[1]]$file))
    rv$tg_results <- fake_run("sqrt"); session$flushReact()
    expect_length(h$entries, 2L)
    # The same analysis again replaces the latest entry rather than adding one.
    rv$tg_results <- fake_run("sqrt"); session$flushReact()
    expect_length(h$entries, 2L)
    # A failed run is not kept.
    rv$tg_results <- list(error = TRUE, message = "no"); session$flushReact()
    expect_length(h$entries, 2L)
    expect_match(paste(as.character(output$history_table), collapse = ""), "<th>Model", fixed = TRUE)
    # Opening run 1 shows it, with its plot data, and is not saved again.
    rv$tg_plot_data <- data.frame(Day = 9)
    g("history_restore")(session, rv, h$entries[[1]])
    session$flushReact()
    expect_identical(rv$tg_results$transform_used, "log")
    expect_identical(rv$tg_plot_data, data.frame(Day = 1:3))
    expect_identical(rv$tg_results$me_history_id, 1L)
    expect_length(h$entries, 2L)
    # Found by fingerprint, for the Bayesian reruns.
    hit <- g("history_find")(session, "tumor_growth", g("history_fingerprint")(fake_run("sqrt")$call_record))
    expect_identical(hit$id, 2L)
    expect_null(g("history_find")(session, "survival", hit$fingerprint))
    expect_null(g("history_find")(session, "tumor_growth", "nope"))
  })
})

test_that("R23.15: several runs download together, a folder each, and the tab gets its run back", {
  rv <- shiny::reactiveValues(tg_results = NULL)
  shiny::testServer(g("run_history_server"), args = list(rv = rv), {
    session$userData$dl_registry <- fake_registry(rv)
    rv$tg_results <- fake_run("log"); session$flushReact()
    rv$tg_results <- fake_run("sqrt"); session$flushReact()
    rv$tg_results <- fake_run("none"); session$flushReact()
    h <- session$userData$run_history
    expect_length(h$entries, 3L)
    # Newest first: rows 2 and 3 are runs 2 and 1.
    session$setInputs(history_table_rows_selected = c(2, 3))
    before <- shiny::isolate(h$changed())
    when3 <- h$entries[[3]]$when
    dir <- tempfile("runs_")
    utils::unzip(output$history_download, exdir = dir)
    files <- list.files(dir, recursive = TRUE)
    expect_true("README.txt" %in% files)
    folders <- unique(sub("/.*", "", files[grepl("/", files)]))
    expect_identical(sort(folders), c("run_01_tumor_growth_lme4_mixed_model_log",
                                      "run_02_tumor_growth_lme4_mixed_model_sqrt"))
    a1 <- utils::read.csv(file.path(dir, "run_01_tumor_growth_lme4_mixed_model_log", "tables", "anova.csv"))
    a2 <- utils::read.csv(file.path(dir, "run_02_tumor_growth_lme4_mixed_model_sqrt", "tables", "anova.csv"))
    expect_identical(c(a1$transform, a2$transform), c("log", "sqrt"))
    # The tab shows what it showed before the download, and the list is as it
    # was: putting run 3 back is not a new run (the table, and its selection,
    # are not redrawn).
    session$flushReact()
    expect_identical(rv$tg_results$transform_used, "none")
    expect_length(h$entries, 3L)
    expect_identical(h$entries[[3]]$when, when3)
    expect_identical(shiny::isolate(h$changed()), before)
    # Removing a run deletes its file.
    f <- h$entries[[1]]$file
    session$setInputs(history_table_rows_selected = 3, history_remove = 1)
    expect_length(h$entries, 2L)
    expect_false(file.exists(f))
  })
})

test_that("R23.15: the Results tab is in the app, and the Bayesian launches look in the history first", {
  ui <- paste(as.character(g("runHistoryUI")("history")), collapse = "")
  expect_match(ui, "history-history_download", fixed = TRUE)
  for (f in c("tg_resume.R", "surv_resume.R")) {
    src <- dash_source("R", f)
    expect_match(src, "hit <- history_find(session", fixed = TRUE, info = f)
    # The check comes before the job is launched.
    expect_lt(regexpr("history_find(", src, fixed = TRUE), regexpr("launch_job(", src, fixed = TRUE))
  }
  expect_match(dash_source("R", "app.R"), 'tabPanel("Results",', fixed = TRUE)
})
