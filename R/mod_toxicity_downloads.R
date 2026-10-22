#' Toxicity — all download handlers (CODE_REVIEW.md C.1/C.3).
#'
#' Bundles the toxicity module's `downloadHandler` bindings (CSV, active plot
#' and ZIP) in one helper so additions and behaviour changes happen in one
#' place.
#'
#' `download_tox_active_plot` reaches into `output[[plot_id]]()` to obtain
#' the currently-selected plot. That cross-binding still works after the
#' split because every handler is bound to the same `output` environment
#' the main server uses.
#'
#' @param input,output,session Standard moduleServer arguments.
#' @param rv Shared reactiveValues.
#' @keywords internal
add_tox_download_handlers <- function(input, output, session, rv, plots = list()) {

  # ── Active plot PNG (trajectory / KM / TWM) ───────────────────────────
  output$download_tox_active_plot <- downloadHandler(
    filename = function() {
      if ((input$tox_model_type %||% "bw") == "bw") {
        switch(input$tox_bw_plot_selector %||% "trajectory",
          trajectory = "body_weight_trajectory.png",
          km = "weight_loss_km.png", "toxicity_plot.png")
      } else {
        "therapeutic_window.png"
      }
    },
    content = function(file) {
      w   <- input$tox_plot_width  %||% 800
      h   <- input$tox_plot_height %||% 500
      dpi <- input$tox_plot_dpi    %||% 100
      model <- input$tox_model_type %||% "bw"
      if (model == "bw" && (input$tox_bw_plot_selector %||% "") == "km") {
        p <- tryCatch(output$tox_km_plot(), error = function(e) NULL)
        if (is.null(p)) { showNotification("Plot is not available.", type = "warning"); return(invisible(NULL)) }
        png(file, width = w, height = h, units = "px", res = dpi)
        print(p); dev.off(); return()
      }
      # R20.D50 (therapeutic window only): the plot is built directly; an
      # output cannot be read back inside a download handler.
      p <- if (model == "bw") {
        tryCatch(output[["tox_trajectory_plot"]](), error = function(e) NULL)
      } else if (!is.null(rv$tox_results$twm) && !isTRUE(rv$tox_results$twm$error)) {
        tryCatch(tox_window_plot(rv$tox_results$twm, input, rv$tox_et_plot_groups,
                                 text_size = input$tox_plot_text_size %||% 14),
                 error = function(e) NULL)
      }
      if (is.null(p)) {
        showNotification("Plot is not available.", type = "warning")
        return(invisible(NULL))
      }
      png(file, width = w, height = h, units = "px", res = dpi)
      print(p)
      dev.off()
    }
  )

  # ── CSV downloads ─────────────────────────────────────────────────────
  output$download_tox_model <- downloadHandler("body_weight_model.csv", function(file) {
    req(rv$tox_results, rv$tox_results$mixed_model)
    if (isTRUE(rv$tox_results$mixed_model$error)) return()
    utils::write.csv(rv$tox_results$mixed_model$fixed_effects, file, row.names = FALSE)
  })
  output$download_tox_emmeans_csv <- downloadHandler("body_weight_rates.csv", function(file) {
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$emmeans_table)) return()
    utils::write.csv(as.data.frame(mm$emmeans_table), file, row.names = FALSE)
  })
  output$download_tox_rate_comparisons_csv <- downloadHandler("body_weight_rate_comparisons.csv", function(file) {
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$pairwise_comparisons)) return()
    utils::write.csv(as.data.frame(mm$pairwise_comparisons), file, row.names = FALSE)
  })
  output$download_tox_events <- downloadHandler("weight_loss_events.csv", function(file) {
    req(rv$tox_results, rv$tox_results$weight_loss)
    if (isTRUE(rv$tox_results$weight_loss$error)) return()
    utils::write.csv(rv$tox_results$weight_loss$event_data, file, row.names = FALSE)
  })
  output$download_et_twm_csv <- downloadHandler("therapeutic_window.csv", function(file) {
    req(rv$tox_results, rv$tox_results$twm)
    if (isTRUE(rv$tox_results$twm$error)) return()
    utils::write.csv(rv$tox_results$twm$window_table, file, row.names = FALSE)
  })

  # ---- "Download Results" (CODE_REVIEW.md R20-D.Y) ----
  # Replaced "Download All", which zipped the CSV tables only.
  add_download_results(
    input, output, session, prefix = "tox", module = "toxicity",
    items = reactive(tox_dl_items(
      rv$tox_results,
      traj_plot   = function() plots$traj(),
      km_plot     = function() plots$km(),
      window_plot = function() tox_window_plot(rv$tox_results$twm, input,
                                               rv$tox_et_plot_groups,
                                               text_size = input$tox_plot_text_size %||% 14))),
    readme_header = function() dl_readme_header("Toxicity", rv$tox_results, tox_run_settings),
    plot_size = function() list(width  = input$tox_plot_width  %||% 800,
                                height = input$tox_plot_height %||% 500,
                                dpi    = input$tox_plot_dpi    %||% 100))
}

#' The run's settings, for the README
#' @keywords internal
tox_run_settings <- function(res) {
  list(analysis        = if (identical(res$model_type, "et")) "therapeutic window" else "body weight",
       reference_group = res$reference_group,
       estimation      = res$mixed_model$estimation_method %||% res$estimation_method,
       threshold       = res$weight_loss$threshold)
}

