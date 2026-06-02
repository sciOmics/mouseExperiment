#' Toxicity — frequentist plot renderers + diagnostics (CODE_REVIEW.md C.1/C.3).
#'
#' Holds all 10 frequentist `renderPlot` bindings (Q-Q, residuals-vs-fitted,
#' BW AUC forest, BW trajectory, BW AUC, KM, ET TWM, ET scatter, ET benefit,
#' ET weight-corrected TGI) plus the two `renderUI` group-order/colour
#' wrappers and the palette-preset observers. Split out of `mod_toxicity.R`
#' (was 2,175 LOC) so the main server owns lifecycle / Bayesian renderers /
#' table renderers and this file owns frequentist plot rendering.
#'
#' Bayesian plot renderers (`tox_bw_pp_check_plot`, `tox_bw_mcmc_trace_plot`,
#' `tox_bw_prior_posterior_plot`, `bayes_twm_plot`, and the six
#' `twm_*` aliased plots) intentionally remain in `mod_toxicity.R` because
#' they're interleaved with the Bayesian tab-insertion observer.
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

  output$tox_diag_qq <- renderPlot({
    req(rv$tox_results, rv$tox_results$mixed_model)
    mm <- rv$tox_results$mixed_model
    if (isTRUE(mm$error) || is.null(mm$model)) return(NULL)
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
    r <- tryCatch(stats::residuals(mm$model), error = function(e) NULL)
    f <- tryCatch(stats::fitted(mm$model),    error = function(e) NULL)
    req(!is.null(r), !is.null(f))
    plot(f, r, main = "Residuals vs Fitted",
         xlab = "Fitted Values", ylab = "Residuals",
         pch = 16, col = "#2c7bb6")
    graphics::abline(h = 0, col = "#d7191c", lty = 2, lwd = 2)
  })

  # ---- BW AUC Forest Plot ----

  output$tox_bw_forest_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc      <- rv$tox_results$bw_auc
    if (isTRUE(auc$error) || is.null(auc$comparisons)) return(NULL)
    cmp_type <- input$tox_comparison_type %||% "reference"
    text_sz  <- input$tox_plot_text_size %||% 14
    df  <- auc$comparisons
    ref <- rv$tox_ref_group %||% ""
    if (cmp_type == "reference" && ref != "")
      df <- df[df$Group1 == ref | df$Group2 == ref, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    df$P_Adjusted_BH <- stats::p.adjust(df$P_Value, method = "BH")
    df$Comparison    <- paste(df$Group1, "vs", df$Group2)
    df$Comparison    <- factor(df$Comparison, levels = rev(df$Comparison))
    df$Significant   <- df$P_Adjusted_BH < 0.05
    lbl <- if (cmp_type == "reference") "vs Reference Group" else "All Pairwise"
    ggplot2::ggplot(df,
      ggplot2::aes(x = .data[["Mean_Diff"]], y = .data[["Comparison"]],
                   color = .data[["Significant"]])) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data[["CI_Lower"]], xmax = .data[["CI_Upper"]]),
        height = 0.25, linewidth = 0.8
      ) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(
        values = c("TRUE" = "#d7191c", "FALSE" = "#999999"),
        labels = c("TRUE" = "BH p < 0.05", "FALSE" = "n.s."),
        name = ""
      ) +
      ggplot2::labs(
        title    = "BW AUC Comparison Forest Plot",
        subtitle = sprintf("%s comparisons; Benjamini-Hochberg correction", lbl),
        x = "Mean Difference in AUC (% weight change)", y = NULL
      ) +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "top")
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

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
                     x = "Day", y = "Weight (g)", color = "Treatment") +
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
                     x = "Day", y = "Weight (g)", color = "Treatment", fill = "Treatment") +
        ggplot2::theme_classic(base_size = text_sz)
      if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors) +
                                     ggplot2::scale_fill_manual(values = colors)
      p
    }
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_auc_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bw_auc)
    auc <- rv$tox_results$bw_auc
    if (isTRUE(auc$error)) return(NULL)
    text_sz      <- input$tox_plot_text_size %||% 14
    center       <- input$tox_auc_center_bar %||% "mean"
    err_bars     <- input$tox_auc_error_bar  %||% "SEM"
    pt_sz        <- input$tox_auc_pt_size    %||% 2
    jitter_w     <- input$tox_auc_jitter_width %||% 0.15
    ordered_grps <- ordered_tox_grps()
    colors       <- build_color_vector(input, "tox_color", ordered_grps)
    auc_df       <- auc$auc_per_mouse
    if (!is.null(ordered_grps))
      auc_df$Treatment <- factor(auc_df$Treatment, levels = ordered_grps)

    # Per-group summary for center bar and error bars
    gs <- stats::aggregate(AUC_Pct_Change ~ Treatment, data = auc_df,
      FUN = function(x) c(mean   = mean(x, na.rm = TRUE),
                          median = stats::median(x, na.rm = TRUE),
                          sd     = stats::sd(x, na.rm = TRUE),
                          n      = sum(!is.na(x))))
    gs <- do.call(data.frame, gs)
    names(gs) <- c("Treatment", "Mean", "Median", "SD", "N")
    gs$SEM  <- gs$SD / sqrt(gs$N)
    gs$CI95 <- gs$SEM * stats::qt(0.975, pmax(gs$N - 1, 1))
    if (!is.null(ordered_grps))
      gs$Treatment <- factor(gs$Treatment, levels = ordered_grps)
    gs$center_y <- if (center == "median") gs$Median else gs$Mean

    p <- ggplot2::ggplot(auc_df,
      ggplot2::aes(x = .data[["Treatment"]], y = .data[["AUC_Pct_Change"]])) +
      ggplot2::geom_jitter(ggplot2::aes(color = .data[["Treatment"]]),
                           width = jitter_w, alpha = 0.8, size = pt_sz)

    if (center != "none") {
      p <- p + ggplot2::geom_crossbar(
        data = gs,
        ggplot2::aes(x = .data[["Treatment"]], y = .data[["center_y"]],
                     ymin = .data[["center_y"]], ymax = .data[["center_y"]]),
        width = 0.45, linewidth = 0.8, color = "black", fill = NA,
        inherit.aes = FALSE)
      if (err_bars != "none") {
        err_half <- switch(err_bars,
          SEM = gs$SEM, SD = gs$SD, CI = gs$CI95, NULL)
        if (!is.null(err_half)) {
          gs$ymin_err <- gs$center_y - err_half
          gs$ymax_err <- gs$center_y + err_half
          p <- p + ggplot2::geom_errorbar(
            data = gs,
            ggplot2::aes(x = .data[["Treatment"]], ymin = .data[["ymin_err"]],
                         ymax = .data[["ymax_err"]]),
            width = 0.2, linewidth = 0.8, color = "black",
            inherit.aes = FALSE)
        }
      }
    }

    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Body Weight AUC (% Change from Baseline)",
                   x = "Treatment", y = "AUC of % Weight Change (·days)") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "none")
    if (!is.null(colors)) p <- p +
      ggplot2::scale_color_manual(values = colors)
    p
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
        risk.table = TRUE, pval = TRUE, conf.int = TRUE,
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
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- ordered_et_grps()
    colors       <- build_color_vector(input, "tox_et_color", ordered_grps)
    tbl <- twm$twm_table
    if (!is.null(ordered_grps))
      tbl$Treatment <- factor(tbl$Treatment, levels = ordered_grps)
    else
      tbl$Treatment <- factor(tbl$Treatment, levels = tbl$Treatment[order(tbl$TWM)])
    p <- ggplot2::ggplot(tbl,
           ggplot2::aes(.data[["Treatment"]], .data[["TWM"]], fill = .data[["Treatment"]])) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["TWM"]], 2)), vjust = -0.5, size = 4) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = "Therapeutic Window Metric by Treatment",
                   subtitle = "TWM = TGI(%) / Mean Weight Loss(%). Dashed line = TWM of 1.",
                   x = "Treatment", y = "TWM") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
    if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
    p
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_scatter_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$bivariate)
    biv <- rv$tox_results$bivariate
    if (isTRUE(biv$error)) return(NULL)
    pm           <- biv$per_mouse
    text_sz      <- input$tox_plot_text_size %||% 14
    pt_size      <- input$tox_et_scatter_pt_size %||% 3
    cent_size    <- input$tox_et_centroid_size   %||% 6
    ell_lwd      <- input$tox_et_ellipse_lwd     %||% 0.8
    ell_lty      <- input$tox_et_ellipse_lty     %||% "solid"
    show_ell     <- isTRUE(input$tox_et_show_ellipses)
    show_cent    <- isTRUE(input$tox_et_show_centroids)
    shape_by_grp <- isTRUE(input$tox_et_shape_by_group)
    ordered_grps <- ordered_et_grps()
    colors       <- build_color_vector(input, "tox_et_color", ordered_grps)
    if (!is.null(ordered_grps))
      pm$Treatment <- factor(pm$Treatment, levels = ordered_grps)
    eff_label <- switch(biv$efficacy_metric,
      tgi = "TGI (%)", tumor_auc = "Tumor AUC", "Unknown metric")
    aes_base <- if (shape_by_grp)
      ggplot2::aes(.data[["Efficacy"]], .data[["Max_Pct_Weight_Loss"]],
                   color = .data[["Treatment"]], shape = .data[["Treatment"]])
    else
      ggplot2::aes(.data[["Efficacy"]], .data[["Max_Pct_Weight_Loss"]],
                   color = .data[["Treatment"]])
    p <- ggplot2::ggplot(pm, aes_base) +
      ggplot2::geom_point(size = pt_size, alpha = 0.7)
    if (show_ell && nrow(pm) > 3)
      p <- p + ggplot2::stat_ellipse(level = 0.68, linewidth = ell_lwd,
                                     linetype = ell_lty, show.legend = FALSE)
    if (show_cent && !is.null(biv$per_group)) {
      pg <- biv$per_group
      if (!is.null(ordered_grps))
        pg$Treatment <- factor(pg$Treatment, levels = ordered_grps)
      p <- p + ggplot2::geom_point(data = pg,
                 ggplot2::aes(.data[["Efficacy_Mean"]], .data[["Toxicity_Mean"]],
                              color = .data[["Treatment"]]),
                 shape = 18, size = cent_size, inherit.aes = FALSE)
    }
    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(title = "Safety-Efficacy Scatter",
                   subtitle = paste0("Each point is one mouse. Efficacy: ", eff_label),
                   x = eff_label, y = "Max Weight Loss (%)", color = "Treatment") +
      ggplot2::theme_classic(base_size = text_sz)
    if (!is.null(colors)) p <- p + ggplot2::scale_color_manual(values = colors)
    p
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_benefit_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$benefit)
    ben <- rv$tox_results$benefit
    if (isTRUE(ben$error)) return(NULL)
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- ordered_et_grps()
    bt <- ben$benefit_table
    if (!is.null(ordered_grps))
      bt$Treatment <- factor(bt$Treatment, levels = ordered_grps)
    else
      bt$Treatment <- factor(bt$Treatment, levels = bt$Treatment[order(bt$Benefit_Score)])
    ggplot2::ggplot(bt, ggplot2::aes(.data[["Treatment"]], .data[["Benefit_Score"]],
                                      fill = .data[["Benefit_Score"]] > 0)) +
      ggplot2::geom_col(alpha = 0.8, width = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = round(.data[["Benefit_Score"]], 1)),
                         vjust = ifelse(bt$Benefit_Score >= 0, -0.5, 1.5), size = 4) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.8) +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "#d7191c"),
                                 labels = c("TRUE" = "Net Benefit", "FALSE" = "Net Harm"), name = "") +
      ggplot2::labs(title = "Total Benefit Score by Treatment",
                   subtitle = sprintf("B = AUC(efficacy) - lambda x AUC(toxicity), lambda = %.1f",
                                      ben$lambda),
                   x = "Treatment", y = "Benefit Score") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)

  output$tox_et_wctgi_plot <- renderPlot({
    req(rv$tox_results, rv$tox_results$wctgi)
    wc <- rv$tox_results$wctgi
    if (isTRUE(wc$error)) return(NULL)
    text_sz      <- input$tox_plot_text_size %||% 14
    ordered_grps <- ordered_et_grps()
    comp <- wc$comparison
    if (!is.null(ordered_grps))
      comp$Treatment <- factor(comp$Treatment, levels = ordered_grps)
    plot_df <- data.frame(
      Treatment = rep(comp$Treatment, 2),
      Type      = rep(c("Uncorrected TGI", "Weight-Corrected TGI"), each = nrow(comp)),
      TGI       = c(comp$TGI_Uncorrected, comp$TGI_Corrected),
      stringsAsFactors = FALSE
    )
    ggplot2::ggplot(plot_df,
      ggplot2::aes(.data[["Treatment"]], .data[["TGI"]], fill = .data[["Type"]])) +
      ggplot2::geom_col(position = ggplot2::position_dodge(0.7), width = 0.6, alpha = 0.8) +
      ggplot2::labs(title = "Weight-Corrected vs Standard TGI",
                   subtitle = sprintf("Safety threshold: %d%% weight loss",
                                      round(wc$safety_threshold * 100)),
                   x = "Treatment", y = "TGI (%)", fill = "") +
      ggplot2::theme_classic(base_size = text_sz) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  }, width  = function() input$tox_plot_width  %||% 800,
     height = function() input$tox_plot_height %||% 500)
}
