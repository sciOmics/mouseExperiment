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
  show_analysis_modal(
    "Bayesian tumor growth fit running…",
    "Fitting Bayesian LMM with brms. This may take 3–12 minutes.")

  promise <- .tg_bayes_promise(bayes_params)
  job_id <- launch_job(sid, "bayes_tg", promise,
                       context = list(bayes_params = bayes_params))
  rv$tg_job_id <- job_id
  attach_bayes_poller(input, output, session, rv,
                      job_id_field          = "tg_job_id",
                      finalize_success      = .finalize_tg_success,
                      finalize_failure      = .finalize_tg_failure,
                      finalize_cancellation = .finalize_tg_cancellation)
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
    finalize_cancellation = .finalize_tg_cancellation)
}

# ---- Per-kind specifics (package-private) ---------------------------------

.tg_bayes_promise <- function(bayes_params) {
  promises::future_promise({
    mouseExperiment::bayesian_tumor_growth(
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
      return_model                 = TRUE,
      # plots = TRUE so the MCMC Diagnostics tab can render
      # posterior_dist_plot + residuals_plot (CODE_REVIEW.md DIAGNOSTICS
      # gap 2). Small extra fit-time cost; the dashboard's other plots
      # (growth curves etc.) are still built reactively from data.
      plots                        = TRUE,
      verbose                      = FALSE
    )
  })
}

.finalize_tg_success <- function(result, context,
                                 input, output, session, rv) {
  bayes_params <- context$bayes_params
  result$model_type_used <- "bayes"
  if (!is.null(rv$bayes_models)) rv$bayes_models$tg <- result$model
  else rv$bayes_models <- list(tg = result$model)
  result$model <- NULL
  gc()
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
