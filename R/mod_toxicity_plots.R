#' Toxicity — plot renderers + diagnostics (CODE_REVIEW.md C.1/C.3).
#'
#' Holds the `renderPlot` bindings (residual diagnostics, BW trajectory,
#' weight-loss KM, TWM), the influence tables, the two `renderUI`
#' group-order/colour wrappers and the palette-preset observers. Split out of
#' `mod_toxicity.R` so the main server owns lifecycle and table renderers and
#' this file owns plot rendering.
#'
#' @param input,output,session Standard moduleServer arguments.
#' @param rv Shared reactiveValues.
#' @param ordered_tox_grps,ordered_et_grps Reactives returning user-ordered
#'   group vectors for the BW and ET sub-modules, defined in
#'   `toxicity_server()`.
#' @keywords internal
add_tox_plots <- function(input, output, session, rv,
                          ordered_tox_grps, ordered_et_grps) {

  # ---- Diagnostics ----
  # Prefer the ggplot2 diagnostic plots produced by analyze_body_weight()
  # (consistent with TG and AUC paths). Fall back to a base-R qqnorm /
  # scatter if the backend hasn't produced the ggplots yet (e.g., older
  # cached result).

  output$tox_diag_qq <- renderPlot({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$model)) return(NULL)
    if (!is.null(mm$diag_qq_plot)) return(mm$diag_qq_plot)
    # Fallback (pre-v0.4.8 results)
    r <- tryCatch(stats::residuals(mm$model), error = function(e) NULL)
    req(!is.null(r))
    stats::qqnorm(r, main = "Q-Q Plot of Mixed Model Residuals", pch = 16,
                  col = "#2c7bb6")
    stats::qqline(r, col = "#d7191c", lwd = 2)
  })

  output$tox_diag_resid_fitted <- renderPlot({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$model)) return(NULL)
    if (!is.null(mm$diag_resid_fitted_plot)) return(mm$diag_resid_fitted_plot)
    r <- tryCatch(stats::residuals(mm$model), error = function(e) NULL)
    f <- tryCatch(stats::fitted(mm$model),    error = function(e) NULL)
    req(!is.null(r), !is.null(f))
    plot(f, r, main = "Residuals vs Fitted",
         xlab = "Fitted Values", ylab = "Residuals",
         pch = 16, col = "#2c7bb6")
    graphics::abline(h = 0, col = "#d7191c", lty = 2, lwd = 2)
  })

  output$tox_diag_scale_location <- renderPlot({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$model)) return(NULL)
    if (!is.null(mm$diag_scale_location_plot)) return(mm$diag_scale_location_plot)
    plot(1, 1, type = "n",
         main = "Scale-Location plot not available (re-run analysis)")
  })

  output$tox_diag_re_qq <- renderPlot({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$model)) return(NULL)
    if (!is.null(mm$diag_re_qq_plot)) return(mm$diag_re_qq_plot)
    plot(1, 1, type = "n",
         main = "Random-effects Q-Q not available (re-run analysis)")
  })

  # CODE_REVIEW.md DIAGNOSTICS gap (13) — LMM influence for BW.
  output$tox_show_influence <- shiny::reactive({
    mm <- rv$tox_results$mixed_model
    !is.null(mm) && !isTRUE(mm$error) &&
      (!is.null(mm$diag_cooks_distance) || !is.null(mm$diag_dfbetas))
  })
  shiny::outputOptions(output, "tox_show_influence", suspendWhenHidden = FALSE)

  output$tox_cooks_distance <- renderDT({
    req(rv$tox_results, rv$tox_results$mixed_model)
    cd <- rv$tox_results$mixed_model$diag_cooks_distance
    if (is.null(cd) || nrow(cd) == 0L)
      return(datatable(data.frame(Message = "Cook's distance not available")))
    datatable(format_df_for_display(as.data.frame(cd)),
              options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
  })

  output$tox_dfbetas <- renderDT({
    req(rv$tox_results, rv$tox_results$mixed_model)
    db <- rv$tox_results$mixed_model$diag_dfbetas
    if (is.null(db) || nrow(db) == 0L)
      return(datatable(data.frame(Message = "DFBETAS not available")))
    datatable(format_df_for_display(as.data.frame(db)),
              options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
  })

  # ---- BW Plots ----

  output$tox_group_order_color_ui <- renderUI({
    req(rv$tox_plot_groups)
    group_order_and_color_ui("tox", rv$tox_plot_groups, ns = session$ns)
  })
  observe_palette_preset(input, session, "tox", reactive(rv$tox_plot_groups))

  output$tox_et_group_order_color_ui <- renderUI({
    req(rv$tox_et_plot_groups)
    group_order_and_color_ui("tox_et", rv$tox_et_plot_groups, ns = session$ns)
  })
  observe_palette_preset(input, session, "tox_et", reactive(rv$tox_et_plot_groups))

  output$tox_trajectory_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$trajectory_data)
    td <- rv$tox_results$trajectory_data
    req(nrow(td) > 0)
    traj_type    <- input$tox_traj_type %||% "spaghetti"
    text_sz      <- input$tox_plot_text_size %||% 14
    line_w       <- input$tox_traj_line_width %||% 1
    ordered_grps <- ordered_tox_grps()
    colors       <- build_color_vector(input, "tox_color", ordered_grps)
    y_lab        <- attr(td, "weight_label") %||% "Weight (g)"

    if (traj_type %in% c("spaghetti", "mean_spaghetti")) {
      td$.__LineGroup <- paste(td[["Treatment"]], td[["ID"]], sep = "·")
      if (!is.null(ordered_grps))
        td$Treatment <- factor(td$Treatment, levels = ordered_grps)
      p <- ggplot2::ggplot(td, ggplot2::aes(.data[["Day"]], .data[["Weight"]],
                                            group = .data[[".__LineGroup"]], color = .data[["Treatment"]])) +
        ggplot2::geom_line(alpha = if (traj_type == "mean_spaghetti") 0.3 else 0.6,
                           linewidth = line_w) +
        ggplot2::geom_point(size = 1.5, alpha = if (traj_type == "mean_spaghetti") 0.3 else 0.6) +
        ggplot2::labs(title = "Body Weight Trajectories",
                     x = "Day", y = y_lab, color = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
      if (traj_type == "mean_spaghetti") {
        grp_mean <- stats::aggregate(Weight ~ Treatment + Day, data = td, FUN = mean, na.rm = TRUE)
        if (!is.null(ordered_grps))
          grp_mean$Treatment <- factor(grp_mean$Treatment, levels = ordered_grps)
        p <- p + ggplot2::geom_line(
          data = grp_mean,
          ggplot2::aes(.data[["Day"]], .data[["Weight"]], group = .data[["Treatment"]],
                       color = .data[["Treatment"]]),
          linewidth = line_w * 2.2, inherit.aes = FALSE)
      }
      if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)
      p
    } else {
      agg <- stats::aggregate(Weight ~ Treatment + Day, data = td, FUN = function(x)
        c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = length(x)))
      af  <- do.call(data.frame, agg)
      names(af) <- c("Treatment", "Day", "Mean", "SD", "N")
      af$SEM   <- af$SD / sqrt(af$N)
      af$Error <- if (traj_type == "mean_sd") af$SD else af$SEM
      err_lbl  <- if (traj_type == "mean_sd") "SD" else "SEM"
      if (!is.null(ordered_grps))
        af$Treatment <- factor(af$Treatment, levels = ordered_grps)
      p <- ggplot2::ggplot(af, ggplot2::aes(.data[["Day"]], .data[["Mean"]],
                                            color = .data[["Treatment"]], fill = .data[["Treatment"]])) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["Mean"]] - .data[["Error"]],
                                          ymax = .data[["Mean"]] + .data[["Error"]]),
                             alpha = 0.2, color = NA) +
        ggplot2::geom_line(linewidth = line_w) +
        ggplot2::geom_point(size = 2) +
        ggplot2::labs(title = paste0("Body Weight Over Time (Mean ± ", err_lbl, ")"),
                     x = "Day", y = y_lab, color = "Treatment", fill = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
      if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors) +
                                     ggplot2::scale_fill_manual(values = colors)
      p
    }
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_km_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$weight_loss)
    wl <- rv$tox_results$weight_loss
    if (isTRUE(wl$error)) return(NULL)
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- ordered_tox_grps()
    colors       <- build_color_vector(input, "tox_color", ordered_grps)
    # Apply group order to event_data so ggsurvplot respects it
    event_df <- wl$event_data
    if (!is.null(ordered_grps) && has_col(event_df, "Treatment"))
      event_df$Treatment <- factor(event_df$Treatment, levels = ordered_grps)
    km_fit <- tryCatch(
      survival::survfit(survival::Surv(Time, Event) ~ Treatment, data = event_df),
      error = function(e) wl$km_fit
    )
    pal <- if (!is.null(colors)) unname(colors) else "hue"
    km_plot <- tryCatch(
      survminer::ggsurvplot(km_fit, data = event_df,
        palette    = pal,
        risk.table = TRUE, pval = km_pval_label(event_df), conf.int = TRUE,  # R20.D29
        title = paste0("Time to ", round(wl$threshold * 100), "% Body Weight Loss"),
        xlab = "Day", ylab = "Proportion Without Threshold Weight Loss",
        ggtheme = ggplot2::theme_classic(base_size = text_sz)),
      error = function(e) NULL)
    req(!is.null(km_plot))
    print(km_plot)
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  # ---- ET Plots ----
  output$tox_et_twm_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$twm)
    twm <- rv$tox_results$twm
    if (isTRUE(twm$error)) return(NULL)
    tox_window_plot(twm, input, ordered_et_grps(),
                    text_size = input$tox_plot_text_size %||% 14)
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)
}


#' The therapeutic window plot: TGI against mean worst weight loss
#'
#' R20-K (package v0.27.0). Each arm is a point with its two 95% intervals;
#' the dashed line is the tolerability threshold, and the shaded side of it is
#' "not tolerated". It replaces the TWM bar chart, which ranked arms by a ratio.
#' @param twm The backend's therapeutic-window result.
#' @param input Module input, for the colour pickers.
#' @param ordered_grps Arms in display order, or NULL.
#' @keywords internal
tox_window_plot <- function(twm, input, ordered_grps = NULL, text_size = 14) {
  w   <- twm$window_table
  thr <- twm$tolerability_threshold %||% 20
  lv  <- if (!is.null(ordered_grps)) intersect(ordered_grps, w$Treatment) else w$Treatment
  lv  <- c(lv, setdiff(w$Treatment, lv))
  w$Treatment <- factor(w$Treatment, levels = lv)
  colors <- build_color_vector(input, "tox_et_color", lv)
  p <- ggplot2::ggplot(w, ggplot2::aes(x = .data[["Worst_Loss"]], y = .data[["TGI"]],
                                       colour = .data[["Treatment"]])) +
    ggplot2::annotate("rect", xmin = thr, xmax = Inf, ymin = -Inf, ymax = Inf,
                      fill = "#c62828", alpha = 0.06) +
    ggplot2::geom_vline(xintercept = thr, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75") +
    ggplot2::geom_segment(ggplot2::aes(x = .data[["Worst_Loss_Lower"]],
                                       xend = .data[["Worst_Loss_Upper"]],
                                       yend = .data[["TGI"]]), na.rm = TRUE) +
    ggplot2::geom_segment(ggplot2::aes(y = .data[["TGI_Lower"]], yend = .data[["TGI_Upper"]],
                                       xend = .data[["Worst_Loss"]]), na.rm = TRUE) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = .data[["Treatment"]]), vjust = -1.1,
                       size = text_size / 4, show.legend = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
    ggplot2::labs(
      title = "Therapeutic window",
      subtitle = sprintf(paste0("TGI on day %s against each arm's mean worst weight loss. ",
                                "Bars: 95%% intervals. Dashed line: tolerability ",
                                "threshold (%s%%)."), twm$endpoint_day, format(thr)),
      x = "Mean worst weight loss (% of each animal's baseline)",
      y = "Tumour growth inhibition (%)", colour = "Treatment") +
    ggplot2::theme_classic(base_size = text_size)
  if (!is.null(colors)) p <- p + ggplot2::scale_colour_manual(values = colors)
  p
}
