#' Therapeutic Window — F.7 resume helpers (CODE_REVIEW.md F.7).
#'
#' TWM combines two pre-existing brms fits (Bayesian TG + Bayesian BW)
#' rather than running brms itself. The combine step takes a few
#' minutes — still worth routing through the registry so a tab
#' disconnect mid-combine doesn't lose the result.

#' @keywords internal
launch_twm_bayes_job <- function(sid, twm_params,
                                 input, output, session, rv) {
  if (is.null(sid) || identical(sid, SID_NONE)) {
    .legacy_inline_twm_bayes(twm_params, input, output, session, rv)
    return(invisible(NULL))
  }
  show_analysis_modal(
    "Bayesian TWM running…",
    "Combining TG and BW posterior samples. This may take a few minutes.")
  promise <- .twm_bayes_promise(twm_params)
  job_id  <- launch_job(sid, "bayes_twm", promise,
                        context = list(bayes_params = twm_params))
  rv$twm_job_id <- job_id
  attach_bayes_poller(input, output, session, rv,
                      job_id_field          = "twm_job_id",
                      finalize_success      = .finalize_twm_success,
                      finalize_failure      = .finalize_twm_failure,
                      finalize_cancellation = .finalize_twm_cancellation)
  invisible(NULL)
}

#' @keywords internal
register_twm_resume <- function(input, output, session, rv) {
  register_bayes_resume(
    kind                  = "bayes_twm",
    tab_name              = "Toxicity",
    modal_title           = "Bayesian TWM running…",
    modal_detail          = "Results will appear automatically when complete.",
    input = input, output = output, session = session, rv = rv,
    job_id_field          = "twm_job_id",
    finalize_success      = .finalize_twm_success,
    finalize_failure      = .finalize_twm_failure,
    finalize_cancellation = .finalize_twm_cancellation)
}

.twm_bayes_promise <- function(p) {
  promises::future_promise({
    mouseExperiment::bayesian_therapeutic_window(
      tg_result        = p$tg_result,
      bw_result        = p$bw_result,
      treatment_column = p$tx_col,
      time_column      = p$time_col,
      reference_group  = p$ref_group,
      plots            = TRUE)
  })
}

.finalize_twm_success <- function(result, context,
                                  input, output, session, rv) {
  result$model_type <- "bayes_twm"
  rv$tox_results <- list(model_type = "bayes_twm", twm = result)
  shiny::updateTabsetPanel(session, "tox_results_tabs", selected = "Summary")
  shiny::removeModal()
  if (isTRUE(result$error)) {
    output$tox_status_msg <- shiny::renderUI(
      status_alert(paste("Bayesian TWM failed:", result$message), "danger"))
  } else {
    output$tox_status_msg <- shiny::renderUI(
      status_alert("Bayesian Therapeutic Window analysis completed!",
                   "success"))
  }
}

.finalize_twm_failure <- function(err_msg, input, output, session, rv) {
  shiny::removeModal()
  output$tox_status_msg <- shiny::renderUI(
    status_alert(paste("Bayesian TWM failed:",
                       translate_error(err_msg)), "danger"))
}

.finalize_twm_cancellation <- function(input, output, session, rv) {
  shiny::removeModal()
  output$tox_status_msg <- shiny::renderUI(
    status_alert("Bayesian TWM analysis was cancelled.", "warning"))
}

.legacy_inline_twm_bayes <- function(p, input, output, session, rv) {
  show_analysis_modal(
    "Bayesian TWM running…",
    "Combining TG and BW posterior samples. This may take a few minutes.")
  .twm_bayes_promise(p) %...>% (function(result) {
    .finalize_twm_success(result, list(bayes_params = p),
                          input, output, session, rv)
  }) %...!% (function(err) {
    .finalize_twm_failure(conditionMessage(err),
                          input, output, session, rv)
  })
  invisible(NULL)
}
