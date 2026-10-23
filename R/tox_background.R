# Toxicity: finishing a run, and the body-weight model on the worker
#
# CODE_REVIEW.md R20.D6. The body-weight model (about 15 s) ran in the Shiny
# process and froze every session. It now runs on the background worker;
# when it returns, the weight-loss threshold and the trajectory data (fast)
# are computed here, and `.tox_finish()` stores the run and reports it, as
# the therapeutic-window path does synchronously.

#' Store a toxicity run and report its status
#' @keywords internal
.tox_finish <- function(results, model_type, data, tx_col, n_err, n_total,
                        output, session, rv) {
  rv$tox_results <- results
  shiny::updateTabsetPanel(session, "tox_results_tabs", selected = "Summary")

  # Store treatment groups for colour/order pickers
  rv$tox_plot_groups <- if (!is.null(results$trajectory_data) && nrow(results$trajectory_data) > 0)
    unique(as.character(results$trajectory_data$Treatment))
  else if (!is.null(tx_col) && has_col(data, tx_col))
    unique(as.character(data[[tx_col]]))
  else NULL

  # Data order (so colours match the other tabs), limited to the arms the
  # therapeutic-window table reports.
  rv$tox_et_plot_groups <- if (model_type == "et" && !is.null(tx_col) &&
                               has_col(data, tx_col)) {
    grps    <- unique(as.character(data[[tx_col]]))
    twm_tbl <- results$twm$window_table
    if (!isTRUE(results$twm$error) && !is.null(twm_tbl))
      grps <- grps[grps %in% twm_tbl$Treatment]
    grps
  } else NULL

  # R20.D34: the backend's warnings are shown, not dropped; a failed
  # analysis is reported as a failure with its message.
  analyses <- if (model_type == "bw") results[c("mixed_model", "weight_loss")] else list(results$twm)
  failed   <- Filter(function(x) isTRUE(x$error), analyses)
  if (n_err == n_total) {
    msgs <- unique(vapply(failed, function(x) x$message %||% "", character(1L)))
    output$tox_status_msg <- shiny::renderUI({
      status_alert(paste("Analysis failed:", paste(msgs, collapse = "; ")), "danger")
    })
  } else if (n_err > 0) {
    msgs <- unique(vapply(failed, function(x) x$message %||% "", character(1L)))
    output$tox_status_msg <- shiny::renderUI(shiny::tagList(
      status_alert(sprintf("%d of %d analyses failed: %s", n_err, n_total,
                           paste(msgs, collapse = "; ")), "warning"),
      completion_status(Filter(function(x) !isTRUE(x$error), analyses),
                        "The other analysis completed successfully!")))
  } else {
    output$tox_status_msg <- shiny::renderUI(completion_status(analyses))
  }
}

#' The rest of a body-weight run, once the model is back
#' @keywords internal
.tox_bw_rest <- function(mixed_model, pending, output, session, rv) {
  results <- pending$results
  data    <- pending$data
  cl      <- pending$cols
  results$mixed_model <- mixed_model
  results$weight_loss <- safe_analysis(module = "tox_threshold",
    do.call(mouseExperiment::weight_loss_threshold, pending$wl_args))
  results$call_records <- pending$call_records   # R20-AA: for the R script

  # R20.D34 / R20.D31: plot the body-weight model's own net weight rather
  # than re-deriving it here with a hard-coded mm3 conversion. Lines are
  # drawn per animal (the backend's treatment + ID + cage key).
  mm <- results$mixed_model
  results$trajectory_data <- if (!isTRUE(mm$error) && !is.null(mm$weight_data)) {
    wd <- mm$weight_data
    td <- data.frame(ID = as.character(wd$.MouseKey),
                     Treatment = as.character(wd$Treatment),
                     Day = as.numeric(wd$Day), Weight = as.numeric(wd$Net_Weight),
                     stringsAsFactors = FALSE)
    td <- td[!is.na(td$Weight) & !is.na(td$Day), ]
    attr(td, "weight_label") <- if (isTRUE(mm$model_info$adjust_tumor))
      "Net body weight (g, tumour mass removed)" else "Body weight (g)"
    td
  } else tryCatch({
    # The model failed: show gross weight, labelled as such.
    td <- data.frame(
      ID        = animal_key(data[[cl$tx]], data[[cl$id]],
                             if (has_col(data, cl$cage)) data[[cl$cage]]),
      Treatment = as.character(data[[cl$tx]]),
      Day       = as.numeric(data[[cl$time]]),
      Weight    = as.numeric(data[[cl$wt]]),
      stringsAsFactors = FALSE
    )
    td <- td[!is.na(td$Weight) & !is.na(td$Day), ]
    attr(td, "weight_label") <- "Gross body weight (g)"
    td
  }, error = function(e) NULL)

  n_err <- sum(vapply(results[c("mixed_model", "weight_loss")],
                      function(x) isTRUE(x$error), logical(1)))
  .tox_finish(results, "bw", data, cl$tx, n_err, 2L, output, session, rv)
}

.finalize_bw_success <- function(result, context, input, output, session, rv) {
  shiny::removeModal()
  pending <- context$pending %||% shiny::isolate(rv$tox_pending)
  rv$tox_pending <- NULL
  if (is.null(pending)) return(invisible(NULL))
  log_info(sprintf("analysis_complete | module=tox_mixed_model | elapsed_s=%.2f",
                   round(result$elapsed_seconds %||% NA_real_, 2)))
  .tox_bw_rest(result, pending, output, session, rv)
}

# A failed model still gets the weight-loss analysis and gross-weight plots.
.finalize_bw_failure <- function(err_msg, input, output, session, rv) {
  shiny::removeModal()
  pending <- shiny::isolate(rv$tox_pending)
  rv$tox_pending <- NULL
  res <- background_failure(err_msg, "tox_mixed_model")
  if (is.null(pending)) {   # resumed after a reload: no inputs to go on with
    output$tox_status_msg <- shiny::renderUI(
      status_alert(paste("Analysis failed:", res$message), "danger"))
    return(invisible(NULL))
  }
  .tox_bw_rest(background_failure(err_msg, "tox_mixed_model"), pending, output, session, rv)
}

.finalize_bw_cancellation <- function(input, output, session, rv) {
  shiny::removeModal()
  rv$tox_pending <- NULL
  output$tox_status_msg <- shiny::renderUI(status_alert("The analysis was cancelled.", "warning"))
}

.finalize_bw_lost <- function(input, output, session, rv) {
  shiny::removeModal()
  rv$tox_pending <- NULL
  output$tox_status_msg <- shiny::renderUI(status_alert(paste(
    "This analysis was dismissed or expired before its result was shown.",
    "Run it again."), "warning"))
}
