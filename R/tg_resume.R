# Tumor Growth — F.7 resume helpers (CODE_REVIEW.md F.7).
#
# Routes the Bayesian TG fit through the package-level job registry
# so the fit survives tab disconnect / NPM proxy timeout / browser
# sleep. Both code paths share the same finalizers:
#   * Fresh "Run analysis" click → launch_tg_bayes_job() registers the
#     promise, sets rv$tg_job_id, attaches the polling observer.
#   * Resume from banner → register_tg_resume() handler fires, sets
#     rv$tg_job_id to the existing job, re-attaches the same poller.
# Each path ends in .finalize_tg_success / _failure / _cancellation.

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
    "Bayesian tumor growth fit running\u2026",
    paste("Fitting the Bayesian LMM with brms. The first fit of a model",
          "compiles it (about 20 s); later fits of the same structure reuse it."),
    cancel_id = session$ns("tg_cancel_bayes"))

  task <- .tg_bayes_task(bayes_params)
  # R20-Z: the call exactly as the worker makes it, for the R script.
  job_id <- launch_job(sid, "bayes_tg", task = task,
                       context = list(bayes_params = bayes_params,
                                      call_record = record_call("bayesian_tumor_growth", task$args)))
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
    modal_title           = "Bayesian tumor growth fit running\u2026",
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
      volume_units                 = bayes_params$volume_units %||% "mm3",
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
  # R20.D46: what this run was asked for, for the report.
  result$run_settings <- list(
    model_type      = "bayes",
    random_effects  = bayes_params$random_effects %||% "",
    reference_group = bayes_params$reference_group %||% "",
    necrotic_handling = bayes_params$necrotic_handling %||% "",
    n_rows          = if (!is.null(bayes_params$data)) nrow(bayes_params$data) else NA)
  result$model <- NULL   # not requested (R20.D18); a no-op kept for safety
  result$call_record <- context$call_record   # R20-Z: for the R script
  # R20.D59: the posterior contrasts vs the reference feed the forest plot,
  # which Bayesian runs never had.
  if (is.data.frame(result$pairwise_comparisons) && nrow(result$pairwise_comparisons)) {
    result$reference_comparisons_dunnett <- result$pairwise_comparisons
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

# ---- The randomisation test on the worker (CODE_REVIEW.md R20.D6) ----------

# Runs tumor_growth_statistics() with the test and returns only the test.
# Base environment, so the worker needs only the backend.
.tg_perm_only <- function(...) mouseExperiment::tumor_growth_statistics(...)$permutation_test
environment(.tg_perm_only) <- baseenv()

# The result joins the run it belongs to; a later run discards it.
.finalize_perm_success <- function(result, context, input, output, session, rv) {
  cur <- shiny::isolate(rv$tg_results)
  if (is.null(cur) || !identical(cur$run_token, context$run_token)) return(invisible(NULL))
  if (is.null(result$p_value)) {
    rv$tg_perm_status <- paste(c("The randomisation test could not be run on these data.",
                                 result$analysis_warnings), collapse = " ")
    return(invisible(NULL))
  }
  cur$permutation_test <- result
  rv$tg_results <- cur
  rv$tg_perm_status <- NULL
}

.finalize_perm_failure <- function(err_msg, input, output, session, rv) {
  res <- background_failure(err_msg, "tg_perm")
  rv$tg_perm_status <- paste("The randomisation test failed:", res$message)
}

.finalize_perm_cancellation <- function(input, output, session, rv) {
  rv$tg_perm_status <- "cancelled"
}

.finalize_perm_lost <- function(input, output, session, rv) {
  rv$tg_perm_status <- paste("The randomisation test was dismissed or expired",
                             "before its result was shown. Run the analysis again.")
}

