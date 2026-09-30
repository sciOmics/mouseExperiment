#' Survival — F.7 resume helpers (CODE_REVIEW.md F.7).
#' Same shape as `tg_resume.R`; see that file for the design notes.

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
  show_analysis_modal(
    "Bayesian survival fit running…",
    "Fitting the Bayesian AFT model with brms. This may take 3–12 minutes.")
  shinyjs::runjs(sprintf(
    "window.trackAnalysisStart('survival', 'bayes', '%s');",
    bayes_params$ref_group))
  promise <- .surv_bayes_promise(bayes_params)
  job_id  <- launch_job(sid, "bayes_surv", promise,
                        context = list(bayes_params = bayes_params))
  rv$surv_job_id <- job_id
  attach_bayes_poller(input, output, session, rv,
                      job_id_field          = "surv_job_id",
                      finalize_success      = .finalize_surv_success,
                      finalize_failure      = .finalize_surv_failure,
                      finalize_cancellation = .finalize_surv_cancellation)
  invisible(NULL)
}

#' @keywords internal
register_surv_resume <- function(input, output, session, rv) {
  register_bayes_resume(
    kind                  = "bayes_surv",
    tab_name              = "Survival",
    modal_title           = "Bayesian survival fit running…",
    modal_detail          = "Results will appear automatically when complete.",
    input = input, output = output, session = session, rv = rv,
    job_id_field          = "surv_job_id",
    finalize_success      = .finalize_surv_success,
    finalize_failure      = .finalize_surv_failure,
    finalize_cancellation = .finalize_surv_cancellation)
}

.surv_bayes_promise <- function(p) {
  promises::future_promise({
    mouseExperiment::bayesian_survival(
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
      return_model     = TRUE,
      plots            = TRUE)
  })
}

.finalize_surv_success <- function(result, context,
                                   input, output, session, rv) {
  bayes_params <- context$bayes_params
  result$model_type_used    <- "bayes"
  if (!is.null(bayes_params)) {
    result$event_src_used     <- bayes_params$event_src
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
  if (!is.null(rv$bayes_models)) rv$bayes_models$surv <- result$model
  else rv$bayes_models <- list(surv = result$model)
  result$model <- NULL
  gc()
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
      status_alert(paste("Bayesian analysis failed:", err_msg), "danger"))
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
