#' Dose-Response — F.7 resume helpers (CODE_REVIEW.md F.7).
#' Same shape as `tg_resume.R`; see that file for the design notes.

#' @keywords internal
launch_dr_bayes_job <- function(sid, dr_bayes_params,
                                input, output, session, rv) {
  if (is.null(sid) || identical(sid, SID_NONE)) {
    .legacy_inline_dr_bayes(dr_bayes_params, input, output, session, rv)
    return(invisible(NULL))
  }
  show_analysis_modal(
    "Bayesian dose-response fit running…",
    "Fitting Bayesian Hill/Emax model with brms. This may take 3–12 minutes.")
  shinyjs::runjs(
    "window.trackAnalysisStart('dose_response', 'bayes', 'hill_emax');")
  promise <- .dr_bayes_promise(dr_bayes_params)
  job_id  <- launch_job(sid, "bayes_dr", promise,
                        context = list(bayes_params = dr_bayes_params))
  rv$dr_job_id <- job_id
  attach_bayes_poller(input, output, session, rv,
                      job_id_field          = "dr_job_id",
                      finalize_success      = .finalize_dr_success,
                      finalize_failure      = .finalize_dr_failure,
                      finalize_cancellation = .finalize_dr_cancellation)
  invisible(NULL)
}

#' @keywords internal
register_dr_resume <- function(input, output, session, rv) {
  register_bayes_resume(
    kind                  = "bayes_dr",
    tab_name              = "Dose Response",
    modal_title           = "Bayesian dose-response fit running…",
    modal_detail          = "Results will appear automatically when complete.",
    input = input, output = output, session = session, rv = rv,
    job_id_field          = "dr_job_id",
    finalize_success      = .finalize_dr_success,
    finalize_failure      = .finalize_dr_failure,
    finalize_cancellation = .finalize_dr_cancellation)
}

.dr_bayes_promise <- function(p) {
  promises::future_promise({
    mouseExperiment::bayesian_dose_response(
      df               = p$data,
      dose_column      = p$dose_col,
      volume_column    = p$vol_col,
      treatment_column = p$tx_col,
      day_column       = p$time_col,
      id_column        = p$id_col,
      endpoint_day     = p$time_point,
      reference_group  = p$ref_group,
      prior_strength   = p$prior_str,
      n_chains         = p$n_chains,
      n_warmup         = p$n_warmup,
      n_iter           = p$n_iter,
      backend          = p$backend,
      return_model     = TRUE,
      plots            = TRUE)
  })
}

.finalize_dr_success <- function(result, context,
                                 input, output, session, rv) {
  result$model_type_used <- "bayes_dr"
  if (!is.null(rv$bayes_models)) rv$bayes_models$dr <- result$model
  else rv$bayes_models <- list(dr = result$model)
  result$model <- NULL
  gc()
  rv$dr_results <- result
  shiny::updateTabsetPanel(session, "dr_results_tabs", selected = "Summary")
  shiny::removeModal()
  if (isTRUE(result$error)) {
    err_msg <- result$message
    shinyjs::runjs(sprintf(
      "window.trackAnalysisComplete('dose_response', false, %s);",
      js_string(err_msg)))
    output$dr_status_msg <- shiny::renderUI(
      status_alert(paste("Bayesian DR failed:", err_msg), "danger"))
  } else {
    shinyjs::runjs("window.trackAnalysisComplete('dose_response', true);")
    output$dr_status_msg <- shiny::renderUI(
      status_alert("Bayesian dose-response analysis completed!", "success"))
  }
}

.finalize_dr_failure <- function(err_msg, input, output, session, rv) {
  shiny::removeModal()
  shinyjs::runjs(sprintf(
    "window.trackAnalysisComplete('dose_response', false, %s);",
    js_string(err_msg)))
  output$dr_status_msg <- shiny::renderUI(
    status_alert(paste("Bayesian DR failed:",
                       translate_error(err_msg)), "danger"))
}

.finalize_dr_cancellation <- function(input, output, session, rv) {
  shiny::removeModal()
  output$dr_status_msg <- shiny::renderUI(
    status_alert("Bayesian DR analysis was cancelled.", "warning"))
}

.legacy_inline_dr_bayes <- function(p, input, output, session, rv) {
  show_analysis_modal(
    "Bayesian dose-response fit running…",
    "Fitting Bayesian Hill/Emax model with brms. This may take 3–12 minutes.")
  .dr_bayes_promise(p) %...>% (function(result) {
    .finalize_dr_success(result, list(bayes_params = p),
                         input, output, session, rv)
  }) %...!% (function(err) {
    .finalize_dr_failure(conditionMessage(err),
                         input, output, session, rv)
  })
  invisible(NULL)
}
