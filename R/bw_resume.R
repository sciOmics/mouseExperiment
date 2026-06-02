#' Body Weight (toxicity) — F.7 resume helpers (CODE_REVIEW.md F.7).
#' Same shape as `tg_resume.R`; see that file for the design notes.

#' @keywords internal
launch_bw_bayes_job <- function(sid, bw_bayes_params,
                                input, output, session, rv) {
  if (is.null(sid) || identical(sid, SID_NONE)) {
    .legacy_inline_bw_bayes(bw_bayes_params, input, output, session, rv)
    return(invisible(NULL))
  }
  show_analysis_modal(
    "Bayesian body weight fit running…",
    "Fitting body weight model with brms. This may take 3–12 minutes.")
  promise <- .bw_bayes_promise(bw_bayes_params)
  job_id  <- launch_job(sid, "bayes_bw", promise,
                        context = list(bayes_params = bw_bayes_params))
  rv$tox_bw_job_id <- job_id
  attach_bayes_poller(input, output, session, rv,
                      job_id_field          = "tox_bw_job_id",
                      finalize_success      = .finalize_bw_success,
                      finalize_failure      = .finalize_bw_failure,
                      finalize_cancellation = .finalize_bw_cancellation)
  invisible(NULL)
}

#' @keywords internal
register_bw_resume <- function(input, output, session, rv) {
  register_bayes_resume(
    kind                  = "bayes_bw",
    tab_name              = "Toxicity",
    modal_title           = "Bayesian body weight fit running…",
    modal_detail          = "Results will appear automatically when complete.",
    input = input, output = output, session = session, rv = rv,
    job_id_field          = "tox_bw_job_id",
    finalize_success      = .finalize_bw_success,
    finalize_failure      = .finalize_bw_failure,
    finalize_cancellation = .finalize_bw_cancellation)
}

.bw_bayes_promise <- function(p) {
  promises::future_promise({
    mouseExperiment::bayesian_body_weight(
      df                           = p$data,
      weight_column                = p$wt_col,
      time_column                  = p$time_col,
      treatment_column             = p$tx_col,
      id_column                    = p$id_col,
      cage_column                  = p$cage_col,
      volume_column                = p$vol_col,
      adjust_tumor_weight          = p$adjust,
      tumor_density                = p$density,
      random_effects_specification = p$re_spec,
      reference_group              = p$ref_group,
      prior_strength               = p$prior_str,
      n_chains                     = p$n_chains,
      n_warmup                     = p$n_warmup,
      n_iter                       = p$n_iter,
      backend                      = p$backend,
      return_model                 = TRUE,
      plots                        = TRUE)
  })
}

.finalize_bw_success <- function(result, context,
                                 input, output, session, rv) {
  result$model_type <- "bw_bayes"
  if (!is.null(rv$bayes_models)) rv$bayes_models$bw <- result$model
  else rv$bayes_models <- list(bw = result$model)
  result$model <- NULL
  gc()
  rv$tox_results <- list(
    model_type  = "bw_bayes",
    mixed_model = result,
    bw_auc      = NULL,
    weight_loss = NULL)
  shiny::updateTabsetPanel(session, "tox_results_tabs", selected = "Summary")
  rv$tox_plot_groups <- if (!is.null(result$treatment_effects) &&
                            has_col(as.data.frame(result$treatment_effects),
                                    "Group")) {
    unique(as.character(as.data.frame(result$treatment_effects)$Group))
  } else NULL
  shiny::removeModal()
  if (isTRUE(result$error)) {
    output$tox_status_msg <- shiny::renderUI(
      status_alert(paste("Bayesian BW failed:", result$message), "danger"))
  } else {
    output$tox_status_msg <- shiny::renderUI(
      status_alert("Bayesian body weight analysis completed!", "success"))
  }
}

.finalize_bw_failure <- function(err_msg, input, output, session, rv) {
  shiny::removeModal()
  output$tox_status_msg <- shiny::renderUI(
    status_alert(paste("Bayesian BW failed:",
                       translate_error(err_msg)), "danger"))
}

.finalize_bw_cancellation <- function(input, output, session, rv) {
  shiny::removeModal()
  output$tox_status_msg <- shiny::renderUI(
    status_alert("Bayesian BW analysis was cancelled.", "warning"))
}

.legacy_inline_bw_bayes <- function(p, input, output, session, rv) {
  show_analysis_modal(
    "Bayesian body weight fit running…",
    "Fitting body weight model with brms. This may take 3–12 minutes.")
  .bw_bayes_promise(p) %...>% (function(result) {
    .finalize_bw_success(result, NULL, input, output, session, rv)
  }) %...!% (function(err) {
    .finalize_bw_failure(conditionMessage(err),
                         input, output, session, rv)
  })
  invisible(NULL)
}
