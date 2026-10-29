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
                                     tg_wf_plot, endpoints = NULL, growth_record = NULL) {

  output$download_tg_pp_check <- downloadHandler(
    filename = "posterior_predictive_check.png",
    content = function(file) {
      req(!is.null(rv$tg_results$pp_check_plot))
      ggplot2::ggsave(file, plot = rv$tg_results$pp_check_plot, width = 8,
                      height = 5, dpi = 150)
    }
  )

  # --- Necrosis summary renderer and download ---

  output$tg_necrosis_table <- DT::renderDT({
    req(rv$tg_results)
    req(!isTRUE(rv$tg_results$error))
    ns_df <- rv$tg_results$necrosis_summary
    req(!is.null(ns_df))
    dt_table(
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
      df <- tg_growth_rates_table(rv$tg_results)
      if (!is.null(df)) write.csv(df, file, row.names = FALSE)
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

  # ---- "Download Results" (CODE_REVIEW.md R20-D.X) -------------------------
  # Replaced "Download All", which zipped a fixed set of tables and PNG plots.
  # The checklist is built from what the result holds; plots are drawn only
  # when the ZIP is written.
  tg_dl_items <- reactive({
    res <- rv$tg_results
    if (is.null(res) || isTRUE(res$error)) return(NULL)
    has   <- function(x) !is.null(x) && (!is.data.frame(x) || nrow(x) > 0L)
    get   <- function(x) function() x
    drawn <- function(r) function() tryCatch(r(), error = function(e) NULL)
    mt    <- res$model_type_used %||% "lme4"
    bayes <- identical(mt, "bayes")
    auc   <- identical(mt, "auc")
    draws <- res$posterior_draws
    # F2: the Growth Delay tab's results, as shown, with their call.
    ep    <- function(r) if (is.function(r)) tryCatch(r(), error = function(e) NULL)
    ok    <- function(x) !is.null(x) && !inherits(x$result, "error")
    delay <- ep(endpoints$delay)
    # F9: the plots the script redraws, with the dashboard's styling.
    growth <- ep(growth_record)
    items <- list(
      # Tables
      if (has(res$anova))
        dl_item("anova", "tables", "ANOVA", "anova", "table", get(res$anova),
                "ANOVA table of the fitted model", row_names = TRUE),
      if (has(res$pairwise_comparisons))
        dl_item("comparisons", "tables", "Comparisons", "comparisons", "table",
                get(res$pairwise_comparisons),
                "Comparisons between arms, with the correction selected"),
      if (has(res$treatment_effects))
        dl_item("effects", "tables", if (auc) "AUC by arm" else "Treatment effects",
                "treatment_effects", "table", get(res$treatment_effects),
                if (auc) "Each arm's model-based AUC with its 95% interval"
                else "Estimated effects by arm"),
      if (has(res$growth_rates))
        dl_item("rates", "tables", "Growth rates", "growth_rates", "table",
                function() tg_growth_rates_table(res),
                "Each animal's growth rate, with daily fold change"),
      if (has(res$data_summary))
        dl_item("summary", "tables", "Data summary", "data_summary", "table",
                get(res$data_summary), "Volume by arm and day, on the volume scale"),
      if (has(res$auc_analysis$individual))
        dl_item("auc_individual", "tables", "Per-animal AUC", "auc_individual", "table",
                get(res$auc_analysis$individual),
                "Each animal's trapezoidal AUC, for description only"),
      if (has(res$auc_analysis$summary))
        dl_item("auc_summary", "tables", "Per-animal AUC summary", "auc_summary", "table",
                get(res$auc_analysis$summary),
                "Group summaries of the per-animal AUCs, for description only"),
      if (has(res$necrosis_summary))
        dl_item("necrosis", "tables", "Necrosis summary", "necrosis_summary", "table",
                get(res$necrosis_summary), "Necrotic observations by arm"),
      if (ok(delay))
        dl_item("growth_delay", "tables", "Growth delay", "growth_delay", "table",
                get(delay$result$delay_table),
                "Each arm's median time to the volume, T - C with its interval, the specific growth delay and the log-rank test"),
      if (ok(delay))
        dl_item("growth_delay_animals", "tables", "Growth delay by animal", "growth_delay_by_animal",
                "table", get(delay$result$per_animal),
                "Each animal's time to the volume (Reached = 0: censored at its last measurement)"),
      if (bayes && has(res$posterior_summary))
        dl_item("posterior_summary", "tables", "Posterior summary", "posterior_summary",
                "table", get(res$posterior_summary), "Posterior summaries of the parameters"),
      # Plots
      dl_item("plot_growth", "plots", "Growth curves", "growth_curves", "plot",
              drawn(tg_growth_plot), "Tumour volume over time by arm"),
      if (has(res$growth_rates))
        dl_item("plot_rates", "plots", "Growth rates", "growth_rates", "plot",
                drawn(tg_rates_plot), "Growth-rate distribution by arm"),
      if (has(res$treatment_effects))
        dl_item("plot_effects", "plots", "Treatment effects", "treatment_effects", "plot",
                drawn(tg_auc_plot), "Treatment effects by arm"),
      if (has(res$pairwise_comparisons))
        dl_item("plot_forest", "plots", "Effect-size forest plot", "effect_forest", "plot",
                drawn(tg_forest_plot), "Comparisons with their intervals"),
      if (has(res$growth_rates))
        dl_item("plot_waterfall", "plots", "Waterfall", "waterfall", "plot",
                drawn(tg_wf_plot), "Each animal's growth rate, sorted"),
      # Diagnostics
      if (!bayes && has(res$diag_qq_plot))
        dl_item("diag_qq", "diagnostics", "Residual Q-Q", "residual_qq", "plot",
                get(res$diag_qq_plot), "Q-Q plot of the standardised residuals"),
      if (!bayes && has(res$diag_resid_fitted_plot))
        dl_item("diag_resid", "diagnostics", "Residuals vs fitted", "residuals_vs_fitted",
                "plot", get(res$diag_resid_fitted_plot), "Residuals against fitted values"),
      if (!bayes && has(res$diag_scale_location_plot))
        dl_item("diag_scale", "diagnostics", "Scale-location", "scale_location", "plot",
                get(res$diag_scale_location_plot), "Spread of the residuals against fitted values"),
      if (!bayes && has(res$diag_re_qq_plot))
        dl_item("diag_re_qq", "diagnostics", "Random-effects Q-Q", "random_effects_qq",
                "plot", get(res$diag_re_qq_plot), "Q-Q plot of the per-animal random effects"),
      if (!bayes && has(res$diag_cooks_distance))
        dl_item("diag_cooks", "diagnostics", "Cook's distance (CSV)", "cooks_distance",
                "table", get(res$diag_cooks_distance), "Influence of each animal on the fit"),
      if (!bayes && has(res$diag_dfbetas))
        dl_item("diag_dfbetas", "diagnostics", "DFBETAS (CSV)", "dfbetas", "table",
                get(res$diag_dfbetas), "Influence of each animal on each coefficient"),
      if (bayes && has(res$mcmc_diagnostics))
        dl_item("diag_mcmc", "diagnostics", "MCMC diagnostics (CSV)", "mcmc_diagnostics",
                "table", get(res$mcmc_diagnostics), "Rhat and effective sample sizes"),
      if (bayes && has(draws))
        dl_item("diag_trace", "diagnostics", "MCMC trace plots", "mcmc_trace", "plot",
                function() bayesplot::mcmc_trace(draws), "Each chain's draws, after warmup"),
      if (bayes && has(draws))
        dl_item("diag_rank", "diagnostics", "MCMC rank plots", "mcmc_rank", "plot",
                function() build_mcmc_rank_plot(draws), "Rank-normalised draws by chain"),
      if (bayes)
        dl_item("diag_ppc", "diagnostics", "Posterior predictive coverage", "ppc_coverage",
                "plot", function() build_ppc_coverage_plot(res),
                "Share of observations inside the posterior predictive intervals"),
      if (bayes)
        dl_item("diag_pareto", "diagnostics", "Pareto-k (LOO)", "pareto_k", "plot",
                function() build_pareto_k_plot(res), "Influential observations by PSIS-LOO"),
      if (has(res$pp_check_plot))
        dl_item("diag_pp_check", "diagnostics", "Posterior predictive check", "pp_check",
                "plot", get(res$pp_check_plot), "Observed data against posterior predictions"),
      # Data: the frame this run analysed, recorded with its call (R20-Z).
      # rv$processed_data is rewritten when another tab prepares the data.
      if (has(res$call_record$data %||% rv$processed_data))
        dl_item("data", "data", "Analysed data", "analysed_data", "table",
                get(res$call_record$data %||% rv$processed_data),
                "The data as analysed: one row per measurement, with the derived Day and Volume"),
      # R script, with the data it reads (R20-Z)
      if (!is.null(res$call_record))
        script_dl_item(function(up) tg_script_bundle(res, up, extra = c(
          tg_endpoint_steps(delay), tg_plot_steps(growth))))
    )
    Filter(Negate(is.null), items)
  })

  add_download_results(
    input, output, session, prefix = "tg", module = "tumor_growth",
    items = tg_dl_items,
    readme_header = function() dl_readme_header("Tumor Growth", rv$tg_results,
                                                tg_report_settings),
    plot_size = function() list(width  = input$tg_plot_width  %||% 800,
                                height = input$tg_plot_height %||% 550,
                                dpi    = input$tg_plot_dpi    %||% 100))
}

#' Growth-rate table as downloaded: the backend's table plus the daily fold
#' change and the percentage change per day.
#' @keywords internal
tg_growth_rates_table <- function(res) {
  if (is.null(res$growth_rates)) return(NULL)
  df <- as.data.frame(res$growth_rates)
  if ("growth_rate" %in% colnames(df)) {
    df$daily_fold_change  <- exp(df$growth_rate)
    df$pct_change_per_day <- (exp(df$growth_rate) - 1) * 100
  }
  df
}

