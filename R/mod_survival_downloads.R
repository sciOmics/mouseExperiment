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
      ef <- surv_effects(rv$surv_results)   # R22.2: Bayesian runs too
      if (!is.null(ef)) write.csv(ef$table, file, row.names = FALSE)
    }
  )

  output$download_surv_risk_table <- downloadHandler(
    filename = function() paste0("survival_risk_table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      shinyjs::runjs("window.trackDownload('survival', 'table', 'csv');")
      req(rv$surv_results$km_fit)
      write.csv(surv_risk_table(rv$surv_results$km_fit), file, row.names = FALSE)
    }
  )

  output$download_surv_forest_plot <- downloadHandler(
    filename = function() paste0("forest_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      shinyjs::runjs("window.trackDownload('survival', 'plot', 'png');")
      req(rv$surv_results, !isTRUE(rv$surv_results$error), surv_effects(rv$surv_results))
      p <- surv_forest_plot_styled(rv$surv_results, input)
      req(!is.null(p))
      png(file, width = input$surv_plot_width %||% 900,
                height = input$surv_plot_height %||% 600,
                units = "px", res = input$surv_plot_dpi %||% 100)
      print(p); dev.off()
    }
  )

  # ---- "Download Results" (CODE_REVIEW.md R20-D.Y) ----
  add_download_results(
    input, output, session, prefix = "surv", module = "survival",
    items = reactive(surv_dl_items(
      rv$surv_results,
      km_plot     = function() km_plot_r(),
      forest_plot = function() surv_forest_plot_styled(rv$surv_results, input))),
    readme_header = function() dl_readme_header("Survival", rv$surv_results, surv_run_settings),
    plot_size = function() list(width  = input$surv_plot_width  %||% 800,
                                height = input$surv_plot_height %||% 600,
                                dpi    = input$surv_plot_dpi    %||% 100))
}

#' Risk table from a Kaplan-Meier fit: at risk, events and censored by arm and day
#' @keywords internal
surv_risk_table <- function(km_fit) {
  km_sum <- summary(km_fit)
  data.frame(
    Group      = sub("^Treatment=", "", as.character(km_sum$strata)),
    Time       = km_sum$time,
    N_Risk     = km_sum$n.risk,
    N_Events   = km_sum$n.event,
    N_Censored = km_sum$n.censor,
    Survival   = round(km_sum$surv, 4)
  )
}

#' The forest plot as the Plots tab styles it (labels, legend position)
#' @keywords internal
surv_forest_plot_styled <- function(res, input) {
  ef <- surv_effects(res)
  if (is.null(ef)) return(NULL)
  p <- build_surv_forest_plot(
    res      = ef$forest,
    ref_grp  = res$reference_group,
    pt_color = input$surv_forest_pt_color %||% "steelblue",
    show_hr  = isTRUE(input$surv_forest_show_hr %||% TRUE),
    show_ci  = isTRUE(input$surv_forest_show_ci %||% TRUE),
    show_p   = isTRUE(input$surv_forest_show_p  %||% TRUE),
    measure  = ef$measure, abbrev = ef$abbrev
  )
  if (is.null(p)) return(NULL)
  p <- apply_plot_labels(p,
                         input$surv_plot_title        %||% "",
                         input$surv_plot_subtitle     %||% "",
                         input$surv_plot_title_align  %||% "center",
                         input$surv_opt_text_size     %||% 16,
                         xlab       = input$surv_plot_xaxis_label  %||% "",
                         ylab       = input$surv_plot_yaxis_label  %||% "",
                         legend_lab = input$surv_plot_legend_label %||% "")
  apply_legend_pos(p, input$surv_opt_legend_pos %||% "right")
}

#' A survival run's effects by arm, and what they measure
#'
#' CODE_REVIEW.md R22.2: a Cox run reports hazard ratios in `results`; a
#' Bayesian AFT run reports time ratios with credible intervals in
#' `treatment_effects` and has no `results`. The download, the CSV button,
#' the forest plot and the R script read `results`, so after a Bayesian run
#' they offered nothing, wrote no file, drew nothing and printed NULL.
#' @param res `rv$surv_results`.
#' @return NULL, or list(table, forest, measure, abbrev): the table as shown,
#'   the same in the forest plot's columns (Group, HR, CI_Lower, CI_Upper,
#'   P_Value), and the measure's name and abbreviation.
#' @keywords internal
surv_effects <- function(res) {
  if (is.null(res) || isTRUE(res$error)) return(NULL)
  # R23.11: an AFT run's time ratios, with confidence intervals and
  # corrected p-values.
  if (identical(res$model_type_used, "aft")) {
    tr <- res$time_ratios
    if (!is.data.frame(tr) || !nrow(tr)) return(NULL)
    return(list(table = tr,
                forest = data.frame(Group = tr$Treatment, HR = tr$Time_Ratio, CI_Lower = tr$Lower_CL,
                                    CI_Upper = tr$Upper_CL, P_Value = tr$P_Adjusted,
                                    stringsAsFactors = FALSE),
                measure = "Time Ratio", abbrev = "TR"))
  }
  if (identical(res$model_type_used, "bayes")) {
    te <- res$treatment_effects
    if (!is.data.frame(te) || !nrow(te)) return(NULL)
    return(list(table = te,
                forest = data.frame(Group = te$Group, HR = te$Time_Ratio, CI_Lower = te$Lower_CrI,
                                    CI_Upper = te$Upper_CrI, P_Value = NA_real_,
                                    stringsAsFactors = FALSE),
                measure = "Time Ratio", abbrev = "TR"))
  }
  r <- res$results
  if (!is.data.frame(r) || !nrow(r)) return(NULL)
  list(table = r, forest = r, measure = "Hazard Ratio", abbrev = "HR")
}

#' The run's settings, for the README
#' @keywords internal
surv_run_settings <- function(res) {
  list(model           = res$model_type_used %||% "frequentist",
       method          = res$method_used,
       reference_group = res$reference_group,
       events_from     = res$event_src_used,
       volume_limit    = res$vol_threshold_used,
       deaths_as_events = res$death_events_added,
       p_adjust        = if (is.data.frame(res$results)) unique(res$results$P_Adjust_Method)
                         else res$p_adjust_method,
       distribution    = res$distribution)
}

