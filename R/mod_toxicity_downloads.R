#' Toxicity — all download handlers (CODE_REVIEW.md C.1/C.3).
#'
#' Bundles the 17 `downloadHandler` bindings previously scattered across
#' `mod_toxicity.R` into a single helper. The Bayesian BW downloads were
#' adjacent to their renderers (lines 961-977); the Bayesian TWM downloads
#' were just below the TWM plot factory (1026-1043); the bulk of the
#' frequentist + active-plot + ZIP downloads were at the bottom of the
#' module (1999-2172). All three groups are co-located here so additions
#' and behavior changes happen in one place.
#'
#' `download_tox_active_plot` reaches into `output[[plot_id]]()` to obtain
#' the currently-selected plot. That cross-binding still works after the
#' split because every handler is bound to the same `output` environment
#' the main server uses.
#'
#' @param input,output,session Standard moduleServer arguments.
#' @param rv Shared reactiveValues.
#' @keywords internal
add_tox_download_handlers <- function(input, output, session, rv) {

  # ── Bayesian BW plot downloads (formerly mod_toxicity.R:961-977) ───────
  output$download_tox_bw_pp_check <- downloadHandler(
    filename = function() paste0("bw_posterior_check_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$bayes_models$bw))
      p <- brms::pp_check(rv$bayes_models$bw, type = "dens_overlay", ndraws = 50)
      ggplot2::ggsave(file, plot = p, width = 8, height = 5, dpi = 150)
    }
  )

  output$download_tox_bw_prior_posterior <- downloadHandler(
    filename = function() paste0("bw_prior_posterior_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$tox_results$mixed_model$prior_posterior_plot))
      ggplot2::ggsave(file, plot = rv$tox_results$mixed_model$prior_posterior_plot,
                      width = 8, height = 5, dpi = 150)
    }
  )

  # ── Bayesian TWM downloads (formerly mod_toxicity.R:1026-1043) ─────────
  output$download_bayes_twm_csv <- downloadHandler(
    filename = function() paste0("bayesian_twm_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      req(!is.null(rv$tox_results$twm))
      tbl <- rv$tox_results$twm$twm_table %||% rv$tox_results$twm$treatment_effects
      req(!is.null(tbl))
      write.csv(as.data.frame(tbl), file, row.names = FALSE)
    }
  )

  output$download_bayes_twm_plot <- downloadHandler(
    filename = function() paste0("bayesian_twm_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      req(!is.null(rv$tox_results$twm))
      p <- rv$tox_results$twm$twm_plot %||% rv$tox_results$twm$plot
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 8, height = 5, dpi = 150)
    }
  )

  # ── Active plot PNG (cycles among trajectory / AUC / KM / ET subplots) ─
  output$download_tox_active_plot <- downloadHandler(
    filename = function() {
      if ((input$tox_model_type %||% "bw") == "bw") {
        switch(input$tox_bw_plot_selector %||% "trajectory",
          trajectory = "body_weight_trajectory.png", auc = "body_weight_auc.png",
          km = "weight_loss_km.png", "toxicity_plot.png")
      } else {
        switch(input$tox_et_plot_selector %||% "twm",
          twm = "therapeutic_window.png", scatter = "safety_efficacy_scatter.png",
          benefit = "total_benefit.png", wctgi = "weight_corrected_tgi.png",
          "efficacy_toxicity_plot.png")
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
      plot_id <- if (model == "bw") {
        switch(input$tox_bw_plot_selector %||% "trajectory",
          trajectory = "tox_trajectory_plot", auc = "tox_auc_plot", "tox_trajectory_plot")
      } else {
        switch(input$tox_et_plot_selector %||% "twm",
          twm = "tox_et_twm_plot", scatter = "tox_et_scatter_plot",
          benefit = "tox_et_benefit_plot", wctgi = "tox_et_wctgi_plot",
          "tox_et_twm_plot")
      }
      p <- tryCatch(output[[plot_id]](), error = function(e) NULL)
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
  output$download_tox_emmeans_csv <- downloadHandler("body_weight_emmeans.csv", function(file) {
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$emmeans_table)) return()
    utils::write.csv(mm$emmeans_table, file, row.names = FALSE)
  })
  output$download_tox_auc_csv <- downloadHandler("body_weight_auc.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bw_auc)
    if (isTRUE(rv$tox_results$bw_auc$error)) return()
    utils::write.csv(rv$tox_results$bw_auc$auc_per_mouse, file, row.names = FALSE)
  })
  output$download_tox_nadir_csv <- downloadHandler("body_weight_nadir.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bw_auc)
    if (isTRUE(rv$tox_results$bw_auc$error)) return()
    utils::write.csv(rv$tox_results$bw_auc$nadir_data, file, row.names = FALSE)
  })
  output$download_tox_comparisons_csv <- downloadHandler("bw_auc_comparisons.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc      <- rv$tox_results$bw_auc
    cmp_type <- isolate(input$tox_comparison_type) %||% "reference"
    if (isTRUE(auc$error) || is.null(auc$comparisons)) return()
    df  <- auc$comparisons
    ref <- rv$tox_ref_group %||% ""
    if (cmp_type == "reference" && ref != "")
      df <- df[df$Group1 == ref | df$Group2 == ref, , drop = FALSE]
    if (nrow(df) > 0) df$P_Adjusted_BH <- stats::p.adjust(df$P_Value, method = "BH")
    utils::write.csv(df, file, row.names = FALSE)
  })
  output$download_tox_events <- downloadHandler("weight_loss_events.csv", function(file) {
    req(rv$tox_results, rv$tox_results$weight_loss)
    if (isTRUE(rv$tox_results$weight_loss$error)) return()
    utils::write.csv(rv$tox_results$weight_loss$event_data, file, row.names = FALSE)
  })
  output$download_et_twm_csv <- downloadHandler("therapeutic_window.csv", function(file) {
    req(rv$tox_results, rv$tox_results$twm)
    if (isTRUE(rv$tox_results$twm$error)) return()
    utils::write.csv(rv$tox_results$twm$twm_table, file, row.names = FALSE)
  })
  output$download_et_bivariate_csv <- downloadHandler("safety_efficacy_per_mouse.csv", function(file) {
    req(rv$tox_results, rv$tox_results$bivariate)
    if (isTRUE(rv$tox_results$bivariate$error)) return()
    utils::write.csv(rv$tox_results$bivariate$per_mouse, file, row.names = FALSE)
  })
  output$download_et_benefit_csv <- downloadHandler("total_benefit.csv", function(file) {
    req(rv$tox_results, rv$tox_results$benefit)
    if (isTRUE(rv$tox_results$benefit$error)) return()
    utils::write.csv(rv$tox_results$benefit$benefit_table, file, row.names = FALSE)
  })
  output$download_et_wctgi_csv <- downloadHandler("weight_corrected_tgi.csv", function(file) {
    req(rv$tox_results, rv$tox_results$wctgi)
    if (isTRUE(rv$tox_results$wctgi$error)) return()
    utils::write.csv(rv$tox_results$wctgi$comparison, file, row.names = FALSE)
  })
  output$download_et_report <- downloadHandler("efficacy_toxicity_report.csv", function(file) {
    req(rv$tox_results)
    res    <- rv$tox_results
    report <- NULL
    if (!isTRUE(res$twm$error))    report <- res$twm$twm_table
    if (!isTRUE(res$benefit$error) && !is.null(report))
      report <- merge(report, res$benefit$benefit_table[, c("Treatment","Benefit")],
                      by = "Treatment", all.x = TRUE)
    if (!isTRUE(res$wctgi$error)   && !is.null(report))
      report <- merge(report,
                      res$wctgi$comparison[, c("Treatment","Corrected_TGI","N_Excluded")],
                      by = "Treatment", all.x = TRUE)
    if (is.null(report)) return()
    utils::write.csv(report, file, row.names = FALSE)
  })

  # ── Download All (ZIP) ────────────────────────────────────────────────
  output$download_tox_all <- downloadHandler(
    filename = "toxicity_results.zip",
    content  = function(file) {
      req(rv$tox_results)
      res    <- rv$tox_results
      tmpdir <- tempfile(pattern = "tox_export_")
      dir.create(tmpdir, recursive = TRUE)
      on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

      write_if <- function(df, fname) {
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0)
          utils::write.csv(df, file.path(tmpdir, fname), row.names = FALSE)
      }

      if (res$model_type == "bw") {
        mm <- res$mixed_model
        if (!isTRUE(mm$error)) {
          write_if(mm$fixed_effects,  "bw_fixed_effects.csv")
          write_if(mm$emmeans_table,  "bw_emmeans.csv")
        }
        bwa <- res$bw_auc
        if (!isTRUE(bwa$error)) {
          write_if(bwa$auc_per_mouse, "bw_auc_per_mouse.csv")
          write_if(bwa$auc_summary,   "bw_auc_summary.csv")
          write_if(bwa$nadir_data,    "bw_nadir.csv")
          if (!is.null(bwa$comparisons)) {
            df  <- bwa$comparisons
            df$P_Adjusted_BH <- stats::p.adjust(df$P_Value, method = "BH")
            write_if(df, "bw_auc_comparisons.csv")
          }
        }
        wl <- res$weight_loss
        if (!isTRUE(wl$error)) write_if(wl$event_data, "bw_weight_loss_events.csv")
      } else {
        if (!isTRUE(res$twm$error))      write_if(res$twm$twm_table,          "et_twm.csv")
        if (!isTRUE(res$bivariate$error)) {
          write_if(res$bivariate$per_mouse, "et_per_mouse.csv")
          write_if(res$bivariate$per_group, "et_per_group.csv")
        }
        if (!isTRUE(res$benefit$error))  write_if(res$benefit$benefit_table,  "et_benefit.csv")
        if (!isTRUE(res$wctgi$error)) {
          write_if(res$wctgi$comparison,     "et_corrected_tgi.csv")
          write_if(res$wctgi$excluded_summary, "et_excluded_mice.csv")
        }
      }

      old_wd <- getwd()
      setwd(tmpdir)
      on.exit(setwd(old_wd), add = TRUE)
      csv_files <- list.files(".", pattern = "\\.csv$")
      if (length(csv_files) == 0) return()
      utils::zip(zipfile = file, files = csv_files)
    }
  )
}
