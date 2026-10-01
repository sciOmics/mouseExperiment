#' Tumor Growth — render handlers for the necrosis table + every download
#' handler. CODE_REVIEW.md Round 2 C.1 / C.3 — splits ~230 LOC of download /
#' export wiring out of `mod_tumor_growth.R` (which was 2,659 LOC). The
#' separation is by concern, not by tab: the main server still owns the
#' lifecycle (observers, fit dispatch, status messages) and the per-tab
#' renderers; this file holds only the per-output download glue plus the
#' single render that shares a tab with its download.
#'
#' Reactive plot functions (`tg_growth_plot`, `tg_rates_plot`, etc.) are
#' passed in as arguments rather than captured by closure so the function
#' contract is explicit. The implementation is otherwise unchanged from
#' the pre-split inline version.
#'
#' @param input,output,session Standard moduleServer arguments.
#' @param rv Shared reactiveValues.
#' @param tg_growth_plot,tg_rates_plot,tg_auc_plot,tg_forest_plot,tg_wf_plot
#'   Reactive plot expressions defined in `tumor_growth_server()`.
#' @keywords internal
add_tg_download_handlers <- function(input, output, session, rv,
                                     tg_growth_plot, tg_rates_plot,
                                     tg_auc_plot, tg_forest_plot,
                                     tg_wf_plot) {

  output$download_tg_pp_check <- downloadHandler(
    filename = "posterior_predictive_check.png",
    content = function(file) {
      req(!is.null(rv$tg_results$pp_check_plot))
      ggplot2::ggsave(file, plot = rv$tg_results$pp_check_plot, width = 8,
                      height = 5, dpi = 150)
    }
  )

  # --- Necrosis summary renderer and download ---

  output$tg_necrosis_table <- renderDT({
    req(rv$tg_results)
    req(!isTRUE(rv$tg_results$error))
    ns_df <- rv$tg_results$necrosis_summary
    req(!is.null(ns_df))
    datatable(
      format_df_for_display(as.data.frame(ns_df)),
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE)
    )
  })

  output$download_tg_necrosis <- downloadHandler(
    filename = "necrosis_summary.csv",
    content = function(file) {
      req(!is.null(rv$tg_results$necrosis_summary))
      write.csv(rv$tg_results$necrosis_summary, file, row.names = FALSE)
    }
  )

  # Download handlers (placeholder implementations)
  output$download_tg_anova <- downloadHandler(
    filename = "anova_results.csv",
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'table', 'csv');")
      if (!is.null(rv$tg_results$anova)) {
        write.csv(as.data.frame(rv$tg_results$anova), file, row.names = TRUE)
      }
    }
  )

  output$download_tg_comparisons <- downloadHandler(
    filename = "pairwise_comparisons.csv",
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'table', 'csv');")
      if (!is.null(rv$tg_results$pairwise_comparisons)) {
        write.csv(as.data.frame(rv$tg_results$pairwise_comparisons), file, row.names = FALSE)
      }
    }
  )

  output$download_tg_effects <- downloadHandler(
    filename = "treatment_effects.csv",
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'table', 'csv');")
      if (!is.null(rv$tg_results$treatment_effects)) {
        write.csv(as.data.frame(rv$tg_results$treatment_effects), file, row.names = FALSE)
      }
    }
  )

  output$download_tg_rates <- downloadHandler(
    filename = "growth_rates.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$growth_rates)) {
        df <- as.data.frame(rv$tg_results$growth_rates)
        if ("growth_rate" %in% colnames(df)) {
          df$daily_fold_change  <- exp(df$growth_rate)
          df$pct_change_per_day <- (exp(df$growth_rate) - 1) * 100
        }
        write.csv(df, file, row.names = FALSE)
      }
    }
  )

  output$download_tg_auc_individual <- downloadHandler(
    filename = "auc_individual_data.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$auc_analysis$individual)) {
        write.csv(as.data.frame(rv$tg_results$auc_analysis$individual), file, row.names = FALSE)
      }
    }
  )

  output$download_tg_auc_summary <- downloadHandler(
    filename = "auc_summary.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$auc_analysis$summary)) {
        write.csv(as.data.frame(rv$tg_results$auc_analysis$summary), file, row.names = FALSE)
      }
    }
  )

  output$download_tg_data_summary <- downloadHandler(
    filename = "data_summary.csv",
    content = function(file) {
      if (!is.null(rv$tg_results$data_summary)) {
        write.csv(as.data.frame(rv$tg_results$data_summary), file, row.names = FALSE)
      }
    }
  )

  # Plot download handlers - save as PNG
  output$download_tg_plot_growth <- downloadHandler(
    filename = "tumor_growth_curves.png",
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'plot', 'png');")
      p <- tg_growth_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )

  output$download_tg_plot_rates <- downloadHandler(
    filename = "growth_rate_distribution.png",
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'plot', 'png');")
      p <- tg_rates_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )

  output$download_tg_plot_auc <- downloadHandler(
    filename = "treatment_effects.png",
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'plot', 'png');")
      p <- tg_auc_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 800,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )

  output$download_tg_plot_effect_forest <- downloadHandler(
    filename = "effect_size_forest_plot.png",
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'plot', 'png');")
      p <- tg_forest_plot()
      if (!is.null(p)) {
        png(file, width = input$tg_plot_width %||% 900,
                  height = input$tg_plot_height %||% 600,
                  units = "px", res = input$tg_plot_dpi %||% 100)
        print(p)
        dev.off()
      }
    }
  )

  output$download_tg_plot_waterfall <- downloadHandler(
    filename = "waterfall_plot.png",
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'plot', 'png');")
      p <- tg_wf_plot()
      if (!is.null(p)) {
        png(file, width  = input$tg_plot_width  %||% 800,
                  height = input$tg_plot_height %||% 550,
                  units  = "px", res = input$tg_plot_dpi %||% 100)
        tryCatch(print(p), finally = dev.off())
      }
    }
  )

  # ---- "Download All Results" — ZIP of all available CSV tables + plots ------
  output$download_tg_all <- downloadHandler(
    filename = function() {
      paste0("mouseExperiment_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      shinyjs::runjs("window.trackDownload('tumor_growth', 'all_results', 'zip');")
      # R20.D52: this was file.path(tempdir(), "mouseExperiment_results"),
      # one directory for every session in the process, cleaned only after a
      # successful export -- so a file left by another user's failed export
      # was zipped into the next user's download. One directory per call now,
      # removed however the call ends.
      outdir <- tempfile("tg_results_")
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
      on.exit(unlink(outdir, recursive = TRUE), add = TRUE)

      # Helper: write a CSV if the data exists
      write_if <- function(data, name) {
        if (!is.null(data)) {
          write.csv(as.data.frame(data), file.path(outdir, name), row.names = FALSE)
        }
      }

      res <- rv$tg_results
      if (!is.null(res)) {
        write_if(res$anova,                "anova_results.csv")
        write_if(res$pairwise_comparisons, "pairwise_comparisons.csv")
        write_if(res$treatment_effects,    "treatment_effects.csv")
        write_if(res$growth_rates,         "growth_rates.csv")
        write_if(res$data_summary,         "data_summary.csv")
        write_if(res$auc_analysis$individual, "auc_individual.csv")
        write_if(res$auc_analysis$summary,    "auc_summary.csv")
      }

      # Save available plots as PNG
      save_plot <- function(plot_fn, name, w = 800, h = 600) {
        p <- tryCatch(plot_fn(), error = function(e) NULL)
        if (!is.null(p)) {
          png(file.path(outdir, name), width = w, height = h, units = "px", res = 100)
          tryCatch(print(p), finally = dev.off())
        }
      }
      save_plot(tg_growth_plot, "growth_curves.png")
      save_plot(tg_rates_plot,  "growth_rates.png")
      save_plot(tg_auc_plot,    "treatment_effects.png")
      save_plot(tg_forest_plot, "effect_forest.png", w = 900)
      save_plot(tg_wf_plot,     "waterfall.png")

      # Create ZIP
      files <- list.files(outdir, full.names = TRUE)
      zip(file, files, flags = "-j")   # -j: junk directory paths

      log_info("download_all | module=tumor_growth")
    },
    contentType = "application/zip"
  )
}
