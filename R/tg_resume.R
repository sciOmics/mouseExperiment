#' Tumor Growth — F.7 resume helpers (CODE_REVIEW.md F.7).
#'
#' Routes the Bayesian TG fit through the package-level job registry
#' so the fit survives tab disconnect / NPM proxy timeout / browser
#' sleep. Both code paths share the same finalizers:
#'   * Fresh "Run analysis" click → launch_tg_bayes_job() registers the
#'     promise, sets rv$tg_job_id, attaches the polling observer.
#'   * Resume from banner → register_tg_resume() handler fires, sets
#'     rv$tg_job_id to the existing job, re-attaches the same poller.
#' Each path ends in .finalize_tg_success / _failure / _cancellation.

#' Launch the Bayesian TG fit through the job registry.
#' @keywords internal
launch_tg_bayes_job <- function(sid, bayes_params,
                                input, output, session, rv) {
  # R20.D3: every session now has a registry key -- app.R gives a session
  # without the cookie a local one -- so the inline fallback, which could
  # never run, is gone. This guard covers callers that bypass app.R.
  if (is.null(sid) || identical(sid, SID_NONE)) {
    sid <- make_local_sid()
    rv$session_id <- sid
  }
  # R20.D7: one Bayesian fit per browser at a time. Repeated clicks used to
  # queue several full fits on the one worker.
  if (!is.null(busy <- active_job(sid))) {
    rv$tg_running <- FALSE
    shiny::showNotification(sprintf(paste0(
      "A Bayesian %s fit from this browser is already %s. Wait for it, or ",
      "cancel it from the banner at the top of the page."),
      if (identical(busy$kind, "bayes_surv")) "survival" else "tumour-growth",
      if (identical(busy$status, "queued")) "waiting" else "running"),
      type = "warning", duration = 8)
    return(invisible(NULL))
  }
  show_analysis_modal(
    "Bayesian tumor growth fit running…",
    paste("Fitting the Bayesian LMM with brms. The first fit of a model",
          "compiles it (about 20 s); later fits of the same structure reuse it."),
    cancel_id = session$ns("tg_cancel_bayes"))

  job_id <- launch_job(sid, "bayes_tg", task = .tg_bayes_task(bayes_params),
                       context = list(bayes_params = bayes_params))
  rv$tg_job_id <- job_id
  attach_bayes_poller(input, output, session, rv,
                      job_id_field          = "tg_job_id",
                      finalize_success      = .finalize_tg_success,
                      finalize_failure      = .finalize_tg_failure,
                      finalize_cancellation = .finalize_tg_cancellation,
                      finalize_lost         = .finalize_tg_lost)
  invisible(NULL)
}

#' Register the bayes_tg resume handler. Called from `tumor_growth_server()`.
#' @keywords internal
register_tg_resume <- function(input, output, session, rv) {
  register_bayes_resume(
    kind                  = "bayes_tg",
    tab_name              = "Tumor Growth",
    modal_title           = "Bayesian tumor growth fit running…",
    modal_detail          = "Results will appear automatically when complete.",
    input = input, output = output, session = session, rv = rv,
    job_id_field          = "tg_job_id",
    finalize_success      = .finalize_tg_success,
    finalize_failure      = .finalize_tg_failure,
    finalize_cancellation = .finalize_tg_cancellation,
    finalize_lost         = .finalize_tg_lost,
    cancel_input          = "tg_cancel_bayes")
}

# ---- Per-kind specifics (package-private) ---------------------------------

# The fit, as a task for the job registry's background worker.
.tg_bayes_task <- function(bayes_params) {
  list(
    fun  = mouseExperiment::bayesian_tumor_growth,
    args = list(
      df                           = bayes_params$data,
      time_column                  = bayes_params$col_names$time_column,
      volume_column                = bayes_params$col_names$volume_column,
      id_column                    = bayes_params$col_names$id_column,
      treatment_column             = bayes_params$col_names$treatment_column,
      cage_column                  = bayes_params$col_names$cage_column,
      dose_column                  = bayes_params$dose_col,
      random_effects_specification = bayes_params$random_effects,
      reference_group              = bayes_params$reference_group,
      prior_strength               = bayes_params$prior_strength,
      prior_b                      = bayes_params$manual_priors$prior_b,
      prior_intercept              = bayes_params$manual_priors$prior_intercept,
      prior_sd                     = bayes_params$manual_priors$prior_sd,
      prior_sigma                  = bayes_params$manual_priors$prior_sigma,
      n_chains                     = bayes_params$n_chains,
      n_warmup                     = bayes_params$n_warmup,
      n_iter                       = bayes_params$n_iter,
      backend                      = bayes_params$backend,
      necrotic_column              = bayes_params$necrotic_col,
      necrotic_handling            = bayes_params$necrotic_handling,
      # R20.D18: the model is not sent back. The diagnostics tab draws from
      # the result's plots (built without the model since backend 0.28.0)
      # and its posterior_draws; keeping the model cost ~31 MB per session.
      return_model                 = FALSE,
      plots                        = TRUE,
      verbose                      = FALSE),
    packages = "mouseExperiment")
}

.finalize_tg_success <- function(result, context,
                                 input, output, session, rv) {
  bayes_params <- context$bayes_params
  result$model_type_used <- "bayes"
  result$model <- NULL   # not requested (R20.D18); a no-op kept for safety
  if (!is.null(rv$tg_results) &&
      !identical(rv$tg_results$model_type_used %||% "",
                 result$model_type_used)) {
    rv$tg_results_prev <- shiny::isolate(rv$tg_results)
  }
  rv$tg_results <- result
  shiny::updateTabsetPanel(session, "tg_results_tabs", selected = "Summary")
  rv$tg_running <- FALSE
  if (!is.null(bayes_params)) {
    rv$tg_plot_data <- bayes_params$data
    rv$tg_col_names <- bayes_params$col_names
  }
  shiny::removeModal()
  output$tg_status_msg <- shiny::renderUI(
    status_alert("Bayesian analysis completed successfully!", "success"))
  shinyjs::runjs("window.trackAnalysisComplete('tumor_growth', true);")
}

.finalize_tg_failure <- function(err_msg, input, output, session, rv) {
  rv$tg_running <- FALSE
  shiny::removeModal()
  output$tg_status_msg <- shiny::renderUI(
    status_alert(paste("Bayesian analysis failed:", err_msg), "danger"))
  shinyjs::runjs(sprintf(
    "window.trackAnalysisComplete('tumor_growth', false, %s);",
    js_string(err_msg)))
}

.finalize_tg_cancellation <- function(input, output, session, rv) {
  rv$tg_running <- FALSE
  shiny::removeModal()
  output$tg_status_msg <- shiny::renderUI(
    status_alert("Bayesian analysis was cancelled.", "warning"))
}

# R20.D5: the job left the registry before its result was shown.
.finalize_tg_lost <- function(input, output, session, rv) {
  rv$tg_running <- FALSE
  shiny::removeModal()
  output$tg_status_msg <- shiny::renderUI(status_alert(paste(
    "This Bayesian analysis was dismissed or expired before its result was",
    "shown. Run it again."), "warning"))
}
