#' Survival — all download handlers (CODE_REVIEW.md C.1/C.3).
#'
#' Six downloadHandler bindings: two Bayesian (`download_surv_pp_check`,
#' `download_surv_prior_posterior`) previously at `mod_survival.R:1039-1054`
#' and four frequentist (`download_surv_plot_km`, `download_surv_cox`,
#' `download_surv_risk_table`, `download_surv_forest_plot`) previously at
#' `mod_survival.R:1481-1551`.
#'
#' `download_surv_plot_km` references the `km_plot_r` reactive defined in
#' `survival_server()`; it's passed in explicitly so the cross-file
#' contract is visible.
#'
#' @param input,output,session Standard moduleServer arguments.
#' @param rv Shared reactiveValues.
#' @param km_plot_r Reactive returning the Kaplan-Meier plot, defined in
#'   `survival_server()`.
#' @keywords internal
add_surv_download_handlers <- function(input, output, session, rv, km_plot_r) {

  output$download_surv_pp_check <- downloadHandler(
    filename = function() paste0("surv_posterior_check_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$surv_results$pp_check_plot))
      ggplot2::ggsave(file, plot = rv$surv_results$pp_check_plot, width = 8,
                      height = 5, dpi = 150)
    }
  )

  output$download_surv_prior_posterior <- downloadHandler(
    filename = function() paste0("surv_prior_posterior_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$surv_results$prior_posterior_plot))
      ggplot2::ggsave(file, plot = rv$surv_results$prior_posterior_plot, width = 8, height = 5, dpi = 150)
    }
  )

  output$download_surv_plot_km <- downloadHandler(
    filename = function() paste0("kaplan_meier_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      shinyjs::runjs("window.trackDownload('survival', 'plot', 'png');")
      km_obj <- isolate(km_plot_r())
      req(!is.null(km_obj))
      png(file, width  = input$surv_plot_width  %||% 800,
                height = input$surv_plot_height %||% 600,
                units  = "px", res = input$surv_plot_dpi %||% 100)
      print(km_obj)
      dev.off()
    }
  )

  output$download_surv_cox <- downloadHandler(
    filename = function() paste0("survival_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      shinyjs::runjs("window.trackDownload('survival', 'table', 'csv');")
      if (!is.null(rv$surv_results$results)) {
        write.csv(rv$surv_results$results, file, row.names = FALSE)
      }
    }
  )

  output$download_surv_risk_table <- downloadHandler(
    filename = function() paste0("survival_risk_table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      shinyjs::runjs("window.trackDownload('survival', 'table', 'csv');")
      req(rv$surv_results$km_fit)
      km_sum <- summary(rv$surv_results$km_fit)
      risk_df <- data.frame(
        Group     = sub("^Treatment=", "", as.character(km_sum$strata)),
        Time      = km_sum$time,
        N_Risk    = km_sum$n.risk,
        N_Events  = km_sum$n.event,
        N_Censored = km_sum$n.censor,
        Survival  = round(km_sum$surv, 4)
      )
      write.csv(risk_df, file, row.names = FALSE)
    }
  )

  output$download_surv_forest_plot <- downloadHandler(
    filename = function() paste0("forest_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      shinyjs::runjs("window.trackDownload('survival', 'plot', 'png');")
      req(rv$surv_results, !isTRUE(rv$surv_results$error), rv$surv_results$results)
      p <- build_surv_forest_plot(
        res      = rv$surv_results$results,
        ref_grp  = rv$surv_results$reference_group,
        pt_color = input$surv_forest_pt_color %||% "steelblue",
        show_hr  = isTRUE(input$surv_forest_show_hr %||% TRUE),
        show_ci  = isTRUE(input$surv_forest_show_ci %||% TRUE),
        show_p   = isTRUE(input$surv_forest_show_p  %||% TRUE)
      )
      req(!is.null(p))
      p <- apply_plot_labels(p,
                             input$surv_plot_title        %||% "",
                             input$surv_plot_subtitle     %||% "",
                             input$surv_plot_title_align  %||% "center",
                             input$surv_opt_text_size     %||% 16,
                             xlab       = input$surv_plot_xaxis_label  %||% "",
                             ylab       = input$surv_plot_yaxis_label  %||% "",
                             legend_lab = input$surv_plot_legend_label %||% "")
      p <- apply_legend_pos(p, input$surv_opt_legend_pos %||% "right")
      png(file, width = input$surv_plot_width %||% 900,
                height = input$surv_plot_height %||% 600,
                units = "px", res = input$surv_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )
}
