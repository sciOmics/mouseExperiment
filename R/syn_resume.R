#' Synergy — F.7 resume helpers (CODE_REVIEW.md F.7).
#' Same shape as `tg_resume.R`; see that file for the design notes.

#' @keywords internal
launch_syn_bayes_job <- function(sid, syn_bayes_params,
                                 input, output, session, rv) {
  if (is.null(sid) || identical(sid, SID_NONE)) {
    .legacy_inline_syn_bayes(syn_bayes_params, input, output, session, rv)
    return(invisible(NULL))
  }
  show_analysis_modal(
    "Bayesian synergy fit running…",
    "Fitting synergy model with brms. This may take 3–12 minutes.")
  shinyjs::runjs(
    "window.trackAnalysisStart('drug_synergy', 'bayes', 'combo_vs_monotherapy');")
  promise <- .syn_bayes_promise(syn_bayes_params)
  job_id  <- launch_job(sid, "bayes_syn", promise,
                        context = list(bayes_params = syn_bayes_params))
  rv$syn_job_id <- job_id
  attach_bayes_poller(input, output, session, rv,
                      job_id_field          = "syn_job_id",
                      finalize_success      = .finalize_syn_success,
                      finalize_failure      = .finalize_syn_failure,
                      finalize_cancellation = .finalize_syn_cancellation)
  invisible(NULL)
}

#' @keywords internal
register_syn_resume <- function(input, output, session, rv) {
  register_bayes_resume(
    kind                  = "bayes_syn",
    tab_name              = "Drug Synergy",
    modal_title           = "Bayesian synergy fit running…",
    modal_detail          = "Results will appear automatically when complete.",
    input = input, output = output, session = session, rv = rv,
    job_id_field          = "syn_job_id",
    finalize_success      = .finalize_syn_success,
    finalize_failure      = .finalize_syn_failure,
    finalize_cancellation = .finalize_syn_cancellation)
}

.syn_bayes_promise <- function(p) {
  promises::future_promise({
    mouseExperiment::bayesian_synergy(
      df                           = p$data,
      volume_column                = p$vol_col,
      time_column                  = p$time_col,
      treatment_column             = p$tx_col,
      id_column                    = p$id_col,
      cage_column                  = p$cage_col,
      drug_a_name                  = p$drug_a,
      drug_b_name                  = p$drug_b,
      combo_name                   = p$combo,
      control_name                 = p$control,
      endpoint_day                 = p$endpoint_day,
      random_effects_specification = p$re_spec,
      prior_strength               = p$prior_str,
      n_chains                     = p$n_chains,
      n_warmup                     = p$n_warmup,
      n_iter                       = p$n_iter,
      backend                      = p$backend,
      return_model                 = TRUE,
      plots                        = TRUE)
  })
}

.finalize_syn_success <- function(result, context,
                                  input, output, session, rv) {
  p <- context$bayes_params
  if (!is.null(rv$bayes_models)) rv$bayes_models$syn <- result$model
  else rv$bayes_models <- list(syn = result$model)
  result$model <- NULL
  gc()
  rv$syn_results <- result
  shiny::updateTabsetPanel(session, "syn_results_tabs", selected = "Summary")
  rv$syn_over_time <- FALSE
  if (!is.null(p)) {
    rv$syn_plot_inputs <- list(
      drug_a = p$drug_a, drug_b = p$drug_b,
      combo  = p$combo,  control = p$control)
  }
  shiny::removeModal()
  if (isTRUE(result$error)) {
    err_msg <- result$message
    shinyjs::runjs(sprintf(
      "window.trackAnalysisComplete('drug_synergy', false, %s);",
      js_string(err_msg)))
    output$syn_status_msg <- shiny::renderUI(
      status_alert(paste("Bayesian synergy failed:", err_msg), "danger"))
  } else {
    shinyjs::runjs("window.trackAnalysisComplete('drug_synergy', true);")
    output$syn_status_msg <- shiny::renderUI(
      status_alert("Bayesian synergy analysis completed!", "success"))
  }
}

.finalize_syn_failure <- function(err_msg, input, output, session, rv) {
  shiny::removeModal()
  shinyjs::runjs(sprintf(
    "window.trackAnalysisComplete('drug_synergy', false, %s);",
    js_string(err_msg)))
  output$syn_status_msg <- shiny::renderUI(
    status_alert(paste("Bayesian synergy failed:",
                       translate_error(err_msg)), "danger"))
}

.finalize_syn_cancellation <- function(input, output, session, rv) {
  shiny::removeModal()
  output$syn_status_msg <- shiny::renderUI(
    status_alert("Bayesian synergy analysis was cancelled.", "warning"))
}

.legacy_inline_syn_bayes <- function(p, input, output, session, rv) {
  show_analysis_modal(
    "Bayesian synergy fit running…",
    "Fitting synergy model with brms. This may take 3–12 minutes.")
  .syn_bayes_promise(p) %...>% (function(result) {
    .finalize_syn_success(result, list(bayes_params = p),
                          input, output, session, rv)
  }) %...!% (function(err) {
    .finalize_syn_failure(conditionMessage(err),
                          input, output, session, rv)
  })
  invisible(NULL)
}
