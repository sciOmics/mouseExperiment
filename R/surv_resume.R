# Survival — F.7 resume helpers (CODE_REVIEW.md F.7).
# Same shape as `tg_resume.R`; see that file for the design notes.

#' @keywords internal
launch_surv_bayes_job <- function(sid, bayes_params,
                                  input, output, session, rv) {
  # R20.D3: every session now has a registry key -- app.R gives a session
  # without the cookie a local one -- so the inline fallback, which could
  # never run, is gone. This guard covers callers that bypass app.R.
  if (is.null(sid) || identical(sid, SID_NONE)) {
    sid <- make_local_sid()
    rv$session_id <- sid
  }
  # R20.D7: one Bayesian fit per browser at a time.
  if (!is.null(busy <- active_job(sid))) {
    rv$surv_running <- FALSE
    shiny::showNotification(sprintf(paste0(
      "A Bayesian %s fit from this browser is already %s. Wait for it, or ",
      "cancel it from the banner at the top of the page."),
      if (identical(busy$kind, "bayes_surv")) "survival" else "tumour-growth",
      if (identical(busy$status, "queued")) "waiting" else "running"),
      type = "warning", duration = 8)
    return(invisible(NULL))
  }
  # R23.15: the same data and settings as a fit saved this session: show it
  # again instead of refitting (the seed is fixed).
  task <- .surv_bayes_task(bayes_params)
  call_record <- record_call("bayesian_survival", task$args)
  hit <- history_find(session, "survival", history_fingerprint(call_record))
  if (!is.null(hit)) {
    history_restore(session, rv, hit)
    rv$surv_running <- FALSE
    shiny::showNotification(sprintf(paste0(
      "This fit has the same data and settings as run %d (at %s), so its results are shown ",
      "again without refitting. Results lists every run."), hit$id, format(hit$when, "%H:%M")),
      type = "message", duration = 10)
    return(invisible(NULL))
  }
  show_analysis_modal(
    "Bayesian survival fit running\u2026",
    # R23.20: how compiling works is in Help, not in the progress window.
    "Fitting the Bayesian AFT model with brms.",
    cancel_id = session$ns("surv_cancel_bayes"))
  # R20.D12: no study strings (the reference arm's name) go to analytics.
  shinyjs::runjs("window.trackAnalysisStart('survival', 'bayes', 'study_control');")
  # R20-AA: the call exactly as the worker makes it, for the R script.
  job_id  <- launch_job(sid, "bayes_surv", task = task,
                        context = list(bayes_params = bayes_params, call_record = call_record))
  rv$surv_job_id <- job_id
  attach_bayes_poller(input, output, session, rv,
                      job_id_field          = "surv_job_id",
                      finalize_success      = .finalize_surv_success,
                      finalize_failure      = .finalize_surv_failure,
                      finalize_cancellation = .finalize_surv_cancellation,
                      finalize_lost         = .finalize_surv_lost)
  invisible(NULL)
}

#' @keywords internal
register_surv_resume <- function(input, output, session, rv) {
  register_bayes_resume(
    kind                  = "bayes_surv",
    tab_name              = "Survival",
    modal_title           = "Bayesian survival fit running\u2026",
    modal_detail          = "Results will appear automatically when complete.",
    input = input, output = output, session = session, rv = rv,
    job_id_field          = "surv_job_id",
    finalize_success      = .finalize_surv_success,
    finalize_failure      = .finalize_surv_failure,
    finalize_cancellation = .finalize_surv_cancellation,
    finalize_lost         = .finalize_surv_lost,
    cancel_input          = "surv_cancel_bayes")
}

# The fit, as a task for the job registry's background worker.
.surv_bayes_task <- function(p) {
  list(
    fun  = mouseExperiment::bayesian_survival,
    args = list(
      df               = p$data,
      time_column      = p$time_col,
      event_column     = p$event_col,
      treatment_column = p$tx_col,
      id_column        = p$id_col,
      cage_column      = p$cage_col,
      family           = p$family,
      reference_group  = p$ref_group,
      prior_strength   = p$prior_str,
      n_chains         = p$n_chains,
      n_warmup         = p$n_warmup,
      n_iter           = p$n_iter,
      backend          = p$backend,
      # R20.D18: the model is not sent back (see tg_resume.R).
      return_model     = FALSE,
      plots            = TRUE),
    packages = "mouseExperiment")
}

.finalize_surv_success <- function(result, context,
                                   input, output, session, rv) {
  bayes_params <- context$bayes_params
  result$model_type_used    <- "bayes"
  if (!is.null(bayes_params)) {
    result$event_src_used     <- bayes_params$event_src
    result$death_events_added <- bayes_params$death_events_added %||% NA_integer_
    result$vol_threshold_used <- bayes_params$vol_threshold
    result$time_method        <- bayes_params$time_method
    result$time_column_used   <- bayes_params$time_col_orig
    result$reference_date     <- bayes_params$ref_date
  }
  if (!isTRUE(result$error) && !is.null(result$survival_data)) {
    tryCatch({
      sd <- result$survival_data
      result$km_fit <- survival::survfit(
        survival::Surv(Time, Event) ~ Treatment, data = sd)
    }, error = function(e) NULL)
  }
  result$model <- NULL   # not requested (R20.D18); a no-op kept for safety
  result <- with_provenance(result, context$call_record, bayes_params$prep_steps,
                            bayes_params$source_file)
  rv$surv_results <- result
  shiny::updateTabsetPanel(session, "surv_results_tabs", selected = "Summary")
  rv$surv_running <- FALSE
  rv$surv_plot_groups <- if (!is.null(result$survival_data) &&
                             has_col(result$survival_data, "Treatment")) {
    unique(as.character(result$survival_data$Treatment))
  } else NULL
  shiny::removeModal()
  if (isTRUE(result$error)) {
    err_msg <- result$message
    shinyjs::runjs(sprintf(
      "window.trackAnalysisComplete('survival', false, %s);",
      js_string(err_msg)))
    output$surv_status_msg <- shiny::renderUI(
      status_alert(paste("Bayesian analysis failed:", translate_error(err_msg)), "danger"))
  } else {
    shinyjs::runjs("window.trackAnalysisComplete('survival', true);")
    output$surv_status_msg <- shiny::renderUI(
      status_alert("Bayesian survival analysis completed!", "success"))
  }
}

.finalize_surv_failure <- function(err_msg, input, output, session, rv) {
  rv$surv_running <- FALSE
  shiny::removeModal()
  shinyjs::runjs(sprintf(
    "window.trackAnalysisComplete('survival', false, %s);",
    js_string(err_msg)))
  output$surv_status_msg <- shiny::renderUI(
    status_alert(paste("Bayesian analysis failed:",
                       translate_error(err_msg)), "danger"))
}

.finalize_surv_cancellation <- function(input, output, session, rv) {
  rv$surv_running <- FALSE
  shiny::removeModal()
  output$surv_status_msg <- shiny::renderUI(
    status_alert("Bayesian survival analysis was cancelled.", "warning"))
}

# R20.D5: the job left the registry before its result was shown.
.finalize_surv_lost <- function(input, output, session, rv) {
  rv$surv_running <- FALSE
  shiny::removeModal()
  output$surv_status_msg <- shiny::renderUI(status_alert(paste(
    "This Bayesian analysis was dismissed or expired before its result was",
    "shown. Run it again."), "warning"))
}
