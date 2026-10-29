# Every run of the session, kept (CODE_REVIEW.md Round 23, R23.15, R23-I).
#
# Each analysis tab used to keep only its latest run: a new run replaced it,
# and Download Results offered only that one. Now every finished run is saved
# to a folder of the session's own, with the shared fields its tab needs to
# show it again. The Results tab lists them: one can be opened in its tab
# again, and several downloaded together, each with its tab's own Download
# Results (tables, plots with the tab's current settings, the R script). A
# Bayesian fit with the same data and settings as a saved one is shown again
# instead of refitted (history_find(), used by tg_resume.R and
# surv_resume.R).

#' The analysis tabs a run is kept for
#'
#' Per tab: its title and navbar tab, the shared field holding its result,
#' and the other shared fields its run handler sets and its tabs read (the
#' run's plot data and groups). Run-state fields (running, job ids) are not
#' kept.
#' @keywords internal
history_modules <- function() {
  list(
    tumor_growth  = list(title = "Tumor Growth", tab = "Tumor Growth", result = "tg_results",
                         fields = c("tg_col_names", "tg_plot_data", "tg_plot_groups", "tg_reference_used")),
    survival      = list(title = "Survival", tab = "Survival", result = "surv_results",
                         fields = "surv_plot_groups"),
    toxicity      = list(title = "Toxicity", tab = "Toxicity", result = "tox_results",
                         fields = c("tox_plot_groups", "tox_et_plot_groups")),
    synergy       = list(title = "Drug Synergy", tab = "Drug Synergy", result = "syn_results",
                         fields = c("syn_over_time", "syn_plot_inputs")),
    dose_response = list(title = "Dose Response", tab = "Dose Response", result = "dr_results",
                         fields = character(0)),
    power         = list(title = "Power Analysis", tab = "Power Analysis", result = "ap_results",
                         fields = character(0)))
}

#' What identifies a run: its call, without the time it was made
#'
#' Two runs with the same fingerprint called the same function with the same
#' data and arguments. NULL without a call record.
#' @param rec A [record_call()] (a result's `call_record`), or a named list of
#'   them (`call_records`).
#' @keywords internal
history_fingerprint <- function(rec) {
  if (is.null(rec)) return(NULL)
  strip <- function(r) if (is.list(r) && !is.null(r$fun)) r[c("fun", "args", "data_arg", "data", "exprs")] else r
  rec <- if (!is.null(rec$fun)) strip(rec) else lapply(rec, strip)
  rlang::hash(rec)
}

# x, or y where x is NA (a lookup that missed).
`%|NA|%` <- function(x, y) if (length(x) && !is.na(x[1])) unname(x[1]) else y

#' A one-line description of a run's model, for the Results table
#' @keywords internal
history_model_label <- function(module, res) {
  tg <- c(lme4 = "LME4 mixed model", auc = "AUC", bayes = "Bayesian LMM")
  cox <- c(coxph = "Cox", coxphf = "Cox (Firth)", logrank = "Log-rank")
  mt <- res$model_type_used %||% "lme4"
  lab <- switch(module,
    tumor_growth = paste0(tg[mt] %|NA|% mt,
                          if (!is.null(res$transform_used) && !identical(mt, "auc"))
                            paste0(", ", res$transform_used)),
    survival = if (identical(res$model_type_used, "aft")) sprintf("AFT (%s)", res$distribution)
               else if (identical(res$model_type_used, "bayes")) "Bayesian AFT"
               else cox[res$method_used %||% "coxph"] %|NA|% (res$method_used %||% "Cox"),
    toxicity = if (!is.null(res$twm)) "Therapeutic window" else "Body weight",
    synergy = if (!is.null(res$synergy_summary)) "Bliss, every evaluable day" else "Bliss, one day",
    dose_response = "Dose response",
    power = res$mode %||% "Power",
    "")
  as.character(lab %||% "")
}

#' The session's run history: its folder and entries
#'
#' Kept in `session$userData`, which the tabs' module sessions share; the
#' entries are also a reactiveVal, for the Results tab.
#' @keywords internal
history_store <- function(session) {
  h <- session$userData$run_history
  if (is.null(h)) {
    h <- new.env(parent = emptyenv())
    h$dir <- tempfile("mouseExperiment-runs-")
    dir.create(h$dir, recursive = TRUE)
    h$entries <- list()
    h$next_id <- 1L
    h$changed <- shiny::reactiveVal(0L)
    session$userData$run_history <- h
    session$onSessionEnded(function() unlink(h$dir, recursive = TRUE))
  }
  h
}

#' Save a finished run
#'
#' A run with the same fingerprint as the tab's latest saved one replaces it:
#' the same analysis repeated, or completed in steps (Tumor Growth's
#' randomisation test arrives after the model).
#' @param rv The shared reactiveValues, read for the tab's fields.
#' @return The entry, invisibly.
#' @keywords internal
history_save <- function(session, rv, module) {
  m <- history_modules()[[module]]
  res <- shiny::isolate(rv[[m$result]])
  if (is.null(res) || isTRUE(res$error) || !is.null(res$me_history_id)) return(invisible(NULL))
  h <- history_store(session)
  fp <- history_fingerprint(res$call_record %||% res$call_records)
  mine <- Filter(function(e) identical(e$module, module), h$entries)
  last <- if (length(mine)) mine[[length(mine)]]
  replace <- !is.null(last) && !is.null(fp) && identical(last$fingerprint, fp)
  id <- if (replace) last$id else h$next_id
  snap <- list(result = res,
               fields = shiny::isolate(stats::setNames(lapply(m$fields, function(f) rv[[f]]), m$fields)))
  file <- file.path(h$dir, sprintf("run_%03d.rds", id))
  saveRDS(snap, file)
  rec <- res$call_record %||% res$call_records[[1]]
  # The same run seen again (put back after a download of several runs, or
  # completed by its randomisation test): its file is updated, its entry
  # and the table are not.
  if (replace && identical(last$run_when, rec$when)) return(invisible(last))
  data <- rec$data
  ref <- res$reference_group %||% res$reference_used %||%
    if ("tg_reference_used" %in% m$fields) shiny::isolate(rv$tg_reference_used)
  entry <- list(id = id, module = module, title = m$title, model = history_model_label(module, res),
                when = Sys.time(), run_when = rec$when,
                rows = if (is.data.frame(data)) nrow(data) else NA_integer_,
                reference = as.character(ref %||% "")[1], fingerprint = fp, file = file)
  if (replace) {
    h$entries[[which(vapply(h$entries, `[[`, 0L, "id") == id)]] <- entry
  } else {
    h$entries[[length(h$entries) + 1L]] <- entry
    h$next_id <- h$next_id + 1L
  }
  h$changed(h$changed() + 1L)
  invisible(entry)
}

#' A saved run of `module` with this fingerprint, or NULL
#' @keywords internal
history_find <- function(session, module, fingerprint) {
  if (is.null(fingerprint)) return(NULL)
  h <- session$userData$run_history
  if (is.null(h)) return(NULL)
  hit <- Filter(function(e) identical(e$module, module) && identical(e$fingerprint, fingerprint),
                h$entries)
  if (length(hit) && file.exists(hit[[length(hit)]]$file)) hit[[length(hit)]] else NULL
}

#' Show a saved run in its tab again
#'
#' Its shared fields are written back, then its result, marked so that it is
#' not saved a second time.
#' @keywords internal
history_restore <- function(session, rv, entry) {
  m <- history_modules()[[entry$module]]
  snap <- readRDS(entry$file)
  res <- snap$result
  res$me_history_id <- entry$id
  shiny::isolate({
    for (f in names(snap$fields)) rv[[f]] <- snap$fields[[f]]
    rv[[m$result]] <- res
  })
  invisible(res)
}

#' Write saved runs into one ZIP, a folder each
#'
#' Each run is shown in its tab for the moment it takes to collect what that
#' tab's Download Results offers (tables, plots drawn with the tab's current
#' settings, diagnostics, data and the R script), then the tab gets back what
#' it showed.
#' @keywords internal
history_export_zip <- function(zipfile, session, rv, entries, fmt = "png") {
  outdir <- tempfile("runs_")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  mods <- history_modules()
  reg <- session$userData$dl_registry %||% list()
  index <- character(0)
  for (e in entries) {
    m <- mods[[e$module]]
    keep <- shiny::isolate(lapply(stats::setNames(nm = c(m$result, m$fields)), function(f) rv[[f]]))
    folder <- sprintf("run_%02d_%s", e$id, gsub("[^A-Za-z0-9]+", "_", tolower(paste(e$title, e$model))))
    ok <- tryCatch({
      history_restore(session, rv, e)
      r <- reg[[e$module]]
      it <- shiny::isolate(r$items())
      header <- shiny::isolate(tryCatch(r$header(), error = function(err) character(0)))
      size <- shiny::isolate(tryCatch(r$size(), error = function(err) list()))
      write_results_dir(file.path(outdir, folder), it, fmt = fmt, size = size,
                        header = c(header, sprintf("Run %d of the session, at %s.", e$id,
                                                   format(e$when, "%Y-%m-%d %H:%M"))))
      TRUE
    }, error = function(err) {
      log_caught_error("run_history", "export", err)
      FALSE
    }, finally = shiny::isolate(for (f in names(keep)) rv[[f]] <- keep[[f]]))
    index <- c(index, sprintf("  %-50s %s, %s, run at %s%s", paste0(folder, "/"), e$title, e$model,
                              format(e$when, "%H:%M"), if (ok) "" else " (could not be written)"))
  }
  writeLines(c("mouseExperiment dashboard: results of several runs", "",
               "Each folder holds one run's Download Results: its tables, plots (drawn with the",
               "tab's settings at the time of this download), diagnostics, data and R script.", "",
               "Runs", index), file.path(outdir, "README.txt"))
  old <- setwd(outdir)
  on.exit(setwd(old), add = TRUE, after = FALSE)
  utils::zip(zipfile, list.files(".", recursive = TRUE), flags = "-q9X")
  invisible(zipfile)
}

#' The Results tab
#' @keywords internal
runHistoryUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    h2("Results"),
    p("Every analysis run in this session, newest first. Open one to show it again in its tab, ",
      "with its tables, plots and downloads; select several to download them together, each in ",
      "a folder of its own. Runs are kept until the session ends."),
    uiOutput(ns("history_note")),
    DT::DTOutput(ns("history_table")),
    br(),
    actionButton(ns("history_open"), "Open in its tab", icon = icon("arrow-up-right-from-square")),
    downloadButton(ns("history_download"), "Download selected (ZIP)"),
    actionButton(ns("history_remove"), "Remove selected", icon = icon("trash"))
  )
}

#' Server side of the Results tab, and the saving of every run
#' @param id Module id.
#' @param rv The shared reactiveValues.
#' @keywords internal
run_history_server <- function(id, rv) {
  moduleServer(id, function(input, output, session) {
    h <- history_store(session)
    mods <- history_modules()
    # Save each tab's run as it finishes.
    lapply(names(mods), function(mod) {
      observeEvent(rv[[mods[[mod]]$result]], history_save(session, rv, mod), ignoreNULL = TRUE)
    })

    entries <- reactive({
      h$changed()
      rev(h$entries)
    })

    output$history_note <- renderUI({
      if (length(entries())) return(NULL)
      div(class = "alert alert-info", icon("circle-info"),
          " No analysis has been run yet in this session; each run will be listed here.")
    })

    output$history_table <- DT::renderDT({
      e <- entries()
      if (!length(e)) return(NULL)
      df <- data.frame(
        Run = vapply(e, `[[`, 0L, "id"),
        Analysis = vapply(e, `[[`, "", "title"),
        Model = vapply(e, `[[`, "", "model"),
        `Run at` = vapply(e, function(x) format(x$when, "%H:%M:%S"), ""),
        Reference = vapply(e, function(x) x$reference %||% "", ""),
        `Data rows` = vapply(e, function(x) as.integer(x$rows), 0L),
        check.names = FALSE, stringsAsFactors = FALSE)
      dt_table(df, rownames = FALSE, selection = list(mode = "multiple", target = "row"),
               options = list(dom = "t", paging = FALSE, ordering = FALSE))
    })

    selected <- function() {
      e <- entries()
      i <- input$history_table_rows_selected
      if (!length(i)) return(list())
      e[i[i <= length(e)]]
    }

    observeEvent(input$history_open, {
      sel <- selected()
      if (length(sel) != 1L) {
        showNotification("Select one run to open.", type = "warning", duration = 5)
        return()
      }
      history_restore(session, rv, sel[[1]])
      updateNavbarPage(session$rootScope(), "main_nav", selected = mods[[sel[[1]]$module]]$tab)
    })

    observeEvent(input$history_remove, {
      sel <- selected()
      if (!length(sel)) return()
      ids <- vapply(sel, `[[`, 0L, "id")
      for (s in sel) unlink(s$file)
      h$entries <- Filter(function(x) !x$id %in% ids, h$entries)
      h$changed(h$changed() + 1L)
    })

    output$history_download <- downloadHandler(
      filename = function() paste0("mouseExperiment_runs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
      content = function(file) {
        sel <- selected()
        if (!length(sel)) sel <- entries()
        history_export_zip(file, session, rv, rev(sel))
      },
      contentType = "application/zip")
  })
}
