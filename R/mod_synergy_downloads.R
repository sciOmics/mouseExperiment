#' Synergy — all download handlers (CODE_REVIEW.md C.1/C.3).
#'
#' The downloadHandler bindings split out of `mod_synergy.R`: one CSV
#' (`download_syn_metrics`) and the plot PNGs (`download_syn_plot_main`,
#' `_trend`, `_ci`, `_tgi`).
#'
#' The plot PNGs each call the reactive that builds that plot from the
#' current settings (R23.14; they held the last rendered plot until 0.56.0),
#' passed in explicitly so the cross-file dependency is visible at the call
#' site.
#'
#' @param input,output,session Standard moduleServer arguments.
#' @param rv Shared reactiveValues.
#' @param current_syn_main_plot,current_syn_trend_plot,current_syn_tgi_plot
#'   Reactives defined in `synergy_server()` that build each panel's plot.
#' @keywords internal
add_syn_download_handlers <- function(input, output, session, rv,
                                      current_syn_main_plot,
                                      current_syn_trend_plot,
                                      current_syn_tgi_plot) {

  output$download_syn_metrics <- downloadHandler(
    filename = function() paste0("synergy_metrics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      df <- if (isTRUE(rv$syn_over_time) && !is.null(rv$syn_results$synergy_summary)) {
        rv$syn_results$synergy_summary
      } else if (!is.null(rv$syn_results$synergy_metrics)) {
        as.data.frame(rv$syn_results$synergy_metrics)
      } else NULL
      if (!is.null(df)) write.csv(df, file, row.names = FALSE)
    }
  )

  output$download_syn_plot_main <- downloadHandler(
    filename = function() paste0("synergy_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      p <- current_syn_main_plot()
      if (is.null(p)) {
        showNotification("Plot is not available.", type = "warning")
        return(invisible(NULL))
      }
      png(file, width = input$syn_plot_width %||% 900,
                height = input$syn_plot_height %||% 600,
                units = "px", res = input$syn_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_syn_plot_trend <- downloadHandler(
    filename = function() paste0("synergy_trend_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      p <- current_syn_trend_plot()
      if (is.null(p)) {
        showNotification("Plot is not available.", type = "warning")
        return(invisible(NULL))
      }
      png(file, width = input$syn_plot_width %||% 900,
                height = input$syn_plot_height %||% 600,
                units = "px", res = input$syn_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )

  output$download_syn_plot_tgi <- downloadHandler(
    filename = function() paste0("synergy_tgi_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      p <- current_syn_tgi_plot()
      if (!is.null(p)) {
        png(file, width = input$syn_plot_width %||% 900,
                  height = input$syn_plot_height %||% 600,
                  units = "px", res = input$syn_plot_dpi %||% 100)
        print(p); dev.off()
      }
    }
  )
}
