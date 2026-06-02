#' Tumor Growth — F.7 resume helpers (CODE_REVIEW.md F.7).
#'
#' Routes the Bayesian TG fit through the package-level job registry so
#' the fit survives tab disconnect / browser sleep / NPM proxy timeout.
#' On reconnect, the resume banner dispatches to the handler registered
#' below; the handler points the same poller at the existing job id and
#' the existing UI surfaces light up exactly as for an inline-completing
#' fit.
#'
#' Two code paths share one consumer:
#'   * Fresh "Run analysis" click → launch_job(), `rv$tg_job_id` set,
#'     poller picks up completion 2 s later.
#'   * Resume from banner → register_resume_handler("bayes_tg") fires,
#'     sets `rv$tg_job_id` to the existing job, attaches the same
#'     poller.
#'
#' Either path ends in `.finalize_tg_bayes_success()` or
#' `.finalize_tg_bayes_failure()` — both idempotent.

#' Launch the Bayesian TG fit through the job registry.
#'
#' Builds the promise, registers it with `launch_job()`, sets up the
#' poller, and shows the analysis modal. Returns invisibly so the
#' caller's `observeEvent` body can `return(NULL)` cleanly.
#'
#' @param sid          Session-identity cookie value for the user
#'                     starting the fit (`rv$session_id`).
#' @param bayes_params The bundled-up parameters list assembled by
#'                     `tumor_growth_server()` from the Bayesian sidebar.
#' @param input,output,session,rv  Module-scoped Shiny objects passed
#'                     through to the poller / finalizer.
#' @keywords internal
launch_tg_bayes_job <- function(sid, bayes_params,
                                input, output, session, rv) {

  if (is.null(sid) || identical(sid, SID_NONE)) {
    # No sid yet — fall back to the legacy inline path. This shouldn't
    # happen in practice because session-id.js sets the cookie before
    # the websocket opens.
    .legacy_inline_tg_bayes(bayes_params, input, output, session, rv)
    return(invisible(NULL))
  }

  show_analysis_modal(
    if (bayes_params$model_type == "gam")
      "Bayesian GAMM tumor growth fit running…"
    else
      "Bayesian tumor growth fit running…",
    if (bayes_params$model_type == "gam")
      "Fitting Bayesian GAMM with brms smoothers. This may take 3–12 minutes."
    else
      "Fitting Bayesian LMM with brms. This may take 3–12 minutes.")

  promise <- promises::future_promise({
    mouseExperiment::bayesian_tumor_growth(
      df                           = bayes_params$data,
      time_column                  = bayes_params$col_names$time_column,
      volume_column                = bayes_params$col_names$volume_column,
      id_column                    = bayes_params$col_names$id_column,
      treatment_column             = bayes_params$col_names$treatment_column,
      cage_column                  = bayes_params$col_names$cage_column,
      dose_column                  = bayes_params$dose_col,
      model_type                   = bayes_params$model_type,
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
      plots                        = FALSE,
      verbose                      = FALSE
    )
  })

  job_id <- launch_job(sid, "bayes_tg", promise,
                       context = list(bayes_params = bayes_params))
  rv$tg_job_id <- job_id
  attach_tg_job_poller(input, output, session, rv)
  invisible(NULL)
}

#' Attach the polling observer that consumes a TG bayes job.
#'
#' Polls `get_job(rv$session_id, rv$tg_job_id)` every 2 s; on terminal
#' status, dispatches to the success / failure / cancellation
#' finalizer, retires the job from the registry, and clears
#' `rv$tg_job_id` so the poller stops. Safe to attach more than once
#' for the same `tg_job_id`: each observer sees the registry transition
#' and races to call `retire_job()`, but `retire_job()` is itself a
#' no-op on a missing entry.
#'
#' @keywords internal
attach_tg_job_poller <- function(input, output, session, rv) {
  shiny::observe({
    job_id <- rv$tg_job_id
    sid    <- rv$session_id
    if (is.null(job_id) || is.null(sid) || identical(sid, SID_NONE)) return()
    shiny::invalidateLater(2000L, session)
    job <- get_job(sid, job_id)
    if (is.null(job)) {
      rv$tg_job_id <- NULL
      return()
    }
    if (job$status == "running") return()
    # Terminal state — consume and stop polling.
    if (job$status == "completed") {
      .finalize_tg_bayes_success(job$result, job$context$bayes_params,
                                 input, output, session, rv)
    } else if (job$status == "failed") {
      .finalize_tg_bayes_failure(job$error, output, rv)
    } else if (job$status == "cancelled") {
      .finalize_tg_bayes_cancellation(output, rv)
    }
    retire_job(sid, job_id)
    rv$tg_job_id <- NULL
  })
}

# ---- Resume handler registration ------------------------------------------

#' Register the bayes_tg resume handler.
#'
#' Called once from inside `tumor_growth_server()`. Captures the
#' module's input / output / session / rv via closure.
#' @keywords internal
register_tg_resume <- function(input, output, session, rv) {
  register_resume_handler("bayes_tg",
    function(job, app_session, app_rv) {
      # rv === app_rv (same reactiveValues, captured here from the
      # module's scope). Use the closured `rv` for module-scoped writes.
      rv$tg_job_id <- job$id
      attach_tg_job_poller(input, output, session, rv)
      shiny::updateNavbarPage(
        session  = app_session,
        inputId  = "main_nav",
        selected = "Tumor Growth")
      if (job$status == "running") {
        show_analysis_modal(
          "Bayesian tumor growth fit running…",
          paste0("Reconnected to a fit started in a previous session.",
                 " Results will appear automatically when complete."))
      } else {
        # Already finished — the poller's next tick will populate rv.
        # Don't show the modal because there's nothing to wait for.
      }
    })
}

# ---- Finalizers (package-private; idempotent so the poller can
#      safely fire more than once on the same registry entry) --------------

.finalize_tg_bayes_success <- function(result, bayes_params,
                                       input, output, session, rv) {
  result$model_type_used <- "bayes"
  if (!is.null(rv$bayes_models)) rv$bayes_models$tg <- result$model
  else rv$bayes_models <- list(tg = result$model)
  result$model <- NULL
  gc()
  # CODE_REVIEW.md Round 2 F.3 — preserve previous result for compare
  # tab when the model_type just changed (e.g. user ran LME4 then bayes).
  if (!is.null(rv$tg_results) &&
      !identical(rv$tg_results$model_type_used %||% "",
                 result$model_type_used)) {
    rv$tg_results_prev <- shiny::isolate(rv$tg_results)
  }
  rv$tg_results  <- result
  shiny::updateTabsetPanel(session, "tg_results_tabs", selected = "Summary")
  rv$tg_running  <- FALSE
  if (!is.null(bayes_params)) {
    rv$tg_plot_data <- bayes_params$data
    rv$tg_col_names <- bayes_params$col_names
  }
  shiny::removeModal()
  output$tg_status_msg <- shiny::renderUI(
    status_alert("Bayesian analysis completed successfully!", "success"))
  shinyjs::runjs("window.trackAnalysisComplete('tumor_growth', true);")
}

.finalize_tg_bayes_failure <- function(err_msg, output, rv) {
  rv$tg_running <- FALSE
  shiny::removeModal()
  output$tg_status_msg <- shiny::renderUI(
    status_alert(paste("Bayesian analysis failed:", err_msg), "danger"))
  shinyjs::runjs(sprintf(
    "window.trackAnalysisComplete('tumor_growth', false, %s);",
    js_string(err_msg)))
}

.finalize_tg_bayes_cancellation <- function(output, rv) {
  rv$tg_running <- FALSE
  shiny::removeModal()
  output$tg_status_msg <- shiny::renderUI(
    status_alert("Bayesian analysis was cancelled.", "warning"))
}

# ---- Legacy fallback ------------------------------------------------------
# Preserves the prior behavior when no sid is available (e.g. running
# the dashboard locally without HTTPS so the Secure cookie is refused).
# The fit still runs but isn't resumable.

.legacy_inline_tg_bayes <- function(bayes_params, input, output, session, rv) {
  show_analysis_modal(
    "Bayesian tumor growth fit running…",
    "Fitting Bayesian model. This may take 3–12 minutes.")
  promises::future_promise({
    mouseExperiment::bayesian_tumor_growth(
      df                           = bayes_params$data,
      time_column                  = bayes_params$col_names$time_column,
      volume_column                = bayes_params$col_names$volume_column,
      id_column                    = bayes_params$col_names$id_column,
      treatment_column             = bayes_params$col_names$treatment_column,
      cage_column                  = bayes_params$col_names$cage_column,
      dose_column                  = bayes_params$dose_col,
      model_type                   = bayes_params$model_type,
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
      plots                        = FALSE,
      verbose                      = FALSE)
  }) %...>% (function(result) {
    .finalize_tg_bayes_success(result, bayes_params,
                               input, output, session, rv)
  }) %...!% (function(err) {
    .finalize_tg_bayes_failure(conditionMessage(err), output, rv)
  })
  invisible(NULL)
}
