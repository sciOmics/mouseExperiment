#' Tumor Growth — plot reactives + plot renderers (CODE_REVIEW.md C.1/C.3).
#'
#' All five plot reactives (`tg_growth_plot`, `tg_rates_plot`, `tg_auc_plot`,
#' `tg_wf_plot`, `tg_forest_plot`), their associated `renderPlot`/`renderUI`
#' bindings, the group-order/colour UI, the palette-preset observer, and the
#' frequentist diagnostic plots (`tg_diag_qq`, `tg_diag_resid_fitted`) live in
#' this file. Split out of `mod_tumor_growth.R` (was 2,659 LOC) so the main
#' server owns lifecycle (observers, fit dispatch, status) and this file owns
#' rendering.
#'
#' The five reactives are returned in a named list so the caller can hand them
#' to `add_tg_download_handlers()` — capturing them explicitly keeps the cross-
#' file contract visible instead of relying on lexical scope.
#'
#' @param input,output,session Standard moduleServer arguments.
#' @param rv Shared reactiveValues.
#' @param ordered_tg_grps Reactive returning the user-ordered group vector,
#'   defined in `tumor_growth_server()`.
#' @return Named list of the five plot reactives.
#' @keywords internal
add_tg_plots <- function(input, output, session, rv, ordered_tg_grps) {

  # Group order and color combined control for Plots tab
  output$tg_group_order_color_ui <- renderUI({
    req(rv$tg_plot_groups)
    group_order_and_color_ui("tg", rv$tg_plot_groups, ns = session$ns)
  })

  # Apply colorblind-friendly palette presets when selector changes
  observe_palette_preset(input, session, "tg", reactive(rv$tg_plot_groups))

  # Growth curves plot.
  # R20.D53: in two layers. The data layer -- plot_tumor_growth(), colours and
  # shapes -- is rebuilt only when the data or those options change. The
  # cosmetic layer (titles, axis labels, legend, scales) is applied on top,
  # from debounced text inputs: typing a 33-character title used to rebuild
  # the whole plot 35 times (19.5 s on the shared process).
  tg_growth_base <- reactive({
    tryCatch({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    req(rv$tg_plot_data, rv$tg_col_names)
    data      <- rv$tg_plot_data
    col_names <- rv$tg_col_names
    ordered_grps <- ordered_tg_grps()
    colors    <- build_color_vector(input, "tg_color", ordered_grps)
    # For dose data, apply factor level order on the original treatment column;
    # plot_tumor_growth will build Group = 'Treatment - Dose: X' internally.
    # For non-dose data, apply order on the treatment column directly.
    if (isTRUE(length(ordered_grps) > 0)) {
      if (!is.null(rv$dose_col) && rv$dose_col %in% colnames(data) &&
          !is.null(rv$settings$col_treatment) && rv$settings$col_treatment %in% colnames(data)) {
        # Derive treatment-level order from the combined group order
        tx_order <- unique(sub(" - Dose: .*$", "", ordered_grps))
        data[[rv$settings$col_treatment]] <- factor(data[[rv$settings$col_treatment]], levels = tx_order)
      } else if (isTRUE(col_names$treatment_column %in% colnames(data))) {
        data[[col_names$treatment_column]] <- factor(data[[col_names$treatment_column]],
                                                     levels = ordered_grps)
      }
    }

    surv_col_mapped <- local({
      col <- rv$settings$col_survival_event
      is.character(col) && length(col) == 1L &&
        nzchar(col) && col != "None" && isTRUE(col %in% colnames(data))
    })
    extra_pts <- local({
      ep <- input$tg_extrapolation_points
      if (is.numeric(ep) && length(ep) == 1L && !is.na(ep)) as.integer(ep) else 0L
    })
    do_extrapolate <- isTRUE(surv_col_mapped) && extra_pts > 0L && isTRUE(input$tg_opt_show_extrap)

    p <- tryCatch({
      mouseExperiment::plot_tumor_growth(
        df                  = data,
        volume_column       = col_names$volume_column,
        day_column          = col_names$time_column,
        # When a dose column is present, pass the original treatment column so
        # plot_tumor_growth can build 'Group = Treatment - Dose: X' itself,
        # matching the names used in rv$tg_plot_groups and the colour pickers.
        treatment_column    = if (!is.null(rv$dose_col)) rv$settings$col_treatment else col_names$treatment_column,
        dose_column         = rv$dose_col,
        ID_column           = col_names$id_column,
        cage_column         = col_names$cage_column,
        survival_column     = if (isTRUE(surv_col_mapped)) rv$settings$col_survival_event else NULL,
        extrapolate_volumes = do_extrapolate,
        extrapolation_points = if (extra_pts > 0L) extra_pts else "all",
        group_summary_line  = isTRUE(input$tg_opt_group_line),
        point_size          = input$tg_opt_growth_pt_size %||% 2,
        # R20.D55: the displayed plot dropped the necrotic x markers, which
        # only the observer's (fallback) plot had.
        necrotic_column     = rv$necrotic_col,
        volume_units        = backend_volume_units(rv)
      ) + ggplot2::theme_classic()
    }, error = function(e) {
      if (!is.null(rv$tg_results$tumor_growth_plot)) return(rv$tg_results$tumor_growth_plot)
      NULL
    })
    req(!is.null(p))
    shapes <- build_shape_vector(input, "tg_shape", ordered_grps)
    if (!is.null(colors)) {
      p <- p + ggplot2::scale_color_manual(values = colors)
      # scale_fill_manual keeps dual-fill shapes (e.g. pch 25, down-triangle)
      # rendered with the correct group colour interior
      p <- p + ggplot2::scale_fill_manual(values = colors)
    }
    if (!is.null(shapes)) {
      p <- p + ggplot2::scale_shape_manual(values = shapes)
    }
    attr(p, "do_extrapolate") <- do_extrapolate
    p
    }, error = function(e) {
      msg <- conditionMessage(e)
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
          label = paste0("Plot error:\n", msg),
          hjust = 0.5, vjust = 0.5, size = 4, color = "red") +
        ggplot2::theme_void()
      attr(p, "plot_error") <- TRUE
      p
    })
  })

  tg_label_text <- debounce(reactive(list(
    title    = input$tg_plot_title        %||% "",
    subtitle = input$tg_plot_subtitle     %||% "",
    xlab     = input$tg_plot_xaxis_label  %||% "",
    ylab     = input$tg_plot_yaxis_label  %||% "",
    legend   = input$tg_plot_legend_label %||% ""
  )), 600)

  tg_growth_plot <- reactive({
    p <- tg_growth_base()
    req(p)
    if (isTRUE(attr(p, "plot_error"))) return(p)
    do_extrapolate <- isTRUE(attr(p, "do_extrapolate"))
    tryCatch({
    # Apply text size, title/subtitle, and legend position
    text_size  <- input$tg_opt_text_size %||% 16
    legend_pos <- input$tg_opt_legend_pos %||% "right"
    lbl <- tg_label_text()
    p <- apply_plot_labels(p, lbl$title, lbl$subtitle,
                           input$tg_plot_title_align %||% "center", text_size,
                           lbl$xlab, lbl$ylab, lbl$legend)
    p <- apply_legend_pos(p, legend_pos)
    # Apply y-axis scale
    yscale <- input$tg_opt_yscale %||% "linear"
    if (yscale == "log10") p <- p + ggplot2::scale_y_log10(expand = ggplot2::expansion(mult = c(0, 0.05)))
    # Apply x-axis start and ensure 5% right padding.
    # coord_cartesian is always applied so the explicit x_upper value
    # (data max + 5%) is honoured regardless of scale expand settings.
    # Pre-treatment (negative) days are kept unless a later start is chosen.
    x_start <- local({
      v <- input$tg_opt_xstart
      first <- suppressWarnings(min(0, rv$tg_plot_data[[rv$tg_col_names$time_column]], na.rm = TRUE))
      if (is.numeric(v) && length(v) == 1L && !is.na(v) && v > 0) v else first
    })
    x_upper <- tryCatch({
      time_col <- rv$tg_col_names$time_column
      x_max <- max(rv$tg_plot_data[[time_col]], na.rm = TRUE)
      x_max + 0.05 * (x_max - x_start)
    }, error = function(e) NA)
    if (!is.na(x_upper))
      p <- suppressMessages(p + ggplot2::coord_cartesian(xlim = c(x_start, x_upper)))
    # Annotate hollow circles = extrapolated when extrapolation is active
    # Added to the package's captions (necrotic markers, mean lines), not
    # in place of them.
    if (isTRUE(do_extrapolate))
      p <- p + ggplot2::labs(caption = paste(c(p$labels$caption,
        "Hollow circles / dashed lines = extrapolated data points after mouse death"),
        collapse = "\n"))
    p
    }, error = function(e) {
      msg <- conditionMessage(e)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
          label = paste0("Plot error:\n", msg),
          hjust = 0.5, vjust = 0.5, size = 4, color = "red") +
        ggplot2::theme_void()
    })
  })

  # R20.D55: extrapolation needs a mapped Event/Censor column; without one the
  # setting did nothing, silently.
  output$tg_extrap_note <- renderUI({
    ep <- input$tg_extrapolation_points
    col <- rv$settings$col_survival_event
    mapped <- is.character(col) && length(col) == 1L && nzchar(col) && col != "None"
    if (!is.numeric(ep) || is.na(ep) || ep <= 0 || mapped) return(NULL)
    helpText(icon("circle-info"), " No Event/Censor column is mapped on Data Upload, ",
             "so no curve is extrapolated: the plot cannot tell an animal that ",
             "died from one that was simply not measured.")
  })

  output$tg_plot_growth <- renderPlot({
    tg_growth_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)

  tg_rates_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    req(rv$tg_results$growth_rates)
    ordered_grps <- ordered_tg_grps()
    colors     <- remap_colors_to_treatment(build_color_vector(input, "tg_color", ordered_grps))
    shapes     <- remap_shapes_to_treatment(build_shape_vector(input, "tg_shape", ordered_grps))
    center_bar <- input$tg_opt_center_bar %||% "mean"
    err_type   <- input$tg_opt_error_bar  %||% "SEM"
    pt_size    <- input$tg_opt_point_size %||% 2.5
    jit_width  <- input$tg_opt_jitter_width %||% 0.2
    text_size  <- input$tg_opt_text_size %||% 12
    legend_pos <- input$tg_opt_legend_pos %||% "auto"
    rates      <- rv$tg_results$growth_rates
    grp_col_r  <- intersect(c("Treatment", "Group"), colnames(rates))[1]
    req(!is.na(grp_col_r))
    val_col_r  <- intersect(c("growth_rate", "GrowthRate"), colnames(rates))[1]
    req(!is.na(val_col_r))
    # Derive treatment-level order (convert " - Dose: " → " " to match TreatmentGroup column)
    gr_order <- if (length(ordered_grps) > 0) {
      raw <- unique(gsub(" - Dose: ", " ", ordered_grps))
      raw[raw %in% as.character(unique(rates[[grp_col_r]]))]
    } else NULL
    if (!is.null(gr_order) && length(gr_order) > 0)
      rates[[grp_col_r]] <- factor(rates[[grp_col_r]], levels = gr_order)

    # Helper: build summary stats
    make_gr_summary <- function(df, grp, val) {
      grps <- if (is.factor(df[[grp]])) levels(df[[grp]]) else unique(as.character(df[[grp]]))
      do.call(rbind, lapply(grps, function(g) {
        x <- df[[val]][as.character(df[[grp]]) == g]
        x <- x[!is.na(x)]
        n <- length(x)
        ctr <- switch(center_bar,
          mean = mean(x), median = median(x),
          # R20.D57: a growth rate is a slope on the log scale, so its mean is
          # already the geometric mean of the daily growth factor; the old
          # exp(mean(log(pmax(x, 1e-9)))) showed -0.044 as 0.000.
          geometric_mean = mean(x), mean(x))
        half <- switch(err_type,
          SEM = sd(x) / sqrt(n), SD = sd(x),
          CI = qt(0.975, df = max(n - 1L, 1L)) * sd(x) / sqrt(n), 0)
        data.frame(grp = g, ctr = ctr, lo = ctr - half, hi = ctr + half,
                   stringsAsFactors = FALSE)
      }))
    }
    summ_r <- make_gr_summary(rates, grp_col_r, val_col_r)
    summ_r$grp <- factor(summ_r$grp, levels = if (is.factor(rates[[grp_col_r]])) levels(rates[[grp_col_r]]) else unique(as.character(rates[[grp_col_r]])))

    # Growth rates are slopes of log(Volume) ~ Day, so always on the log scale.
    # Back-transform to daily fold change when the user has asked for original
    # units; CIs stay symmetric on the log scale (asymmetric in fold-change space).
    do_back <- isTRUE(input$tg_opt_backtransform)
    if (do_back) {
      rates[[val_col_r]] <- exp(rates[[val_col_r]])
      summ_r$ctr <- exp(summ_r$ctr)
      summ_r$lo  <- exp(summ_r$lo)
      summ_r$hi  <- exp(summ_r$hi)
    }
    y_label    <- if (do_back) "Daily Fold Change" else "Growth Rate (log-volume / day)"
    null_yint  <- if (do_back) 1 else 0

    # Build plot from scratch — ensures individual points are drawn LAST (on top)
    p <- ggplot2::ggplot(rates,
           ggplot2::aes(x = .data[[grp_col_r]], y = .data[[val_col_r]],
                        colour = .data[[grp_col_r]], fill = .data[[grp_col_r]],
                        shape = .data[[grp_col_r]])) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Treatment Group", y = y_label,
                    title = "Growth Rates by Treatment Group") +
      ggplot2::geom_hline(yintercept = null_yint, linetype = "dashed", colour = "red", linewidth = 0.5)

    # Error bars
    if (err_type != "none")
      p <- p + ggplot2::geom_errorbar(
        data = summ_r,
        ggplot2::aes(x = grp, ymin = lo, ymax = hi),
        inherit.aes = FALSE,
        colour = "black", width = 0.25, linewidth = 0.9)

    # Center bar
    if (center_bar != "none")
      p <- p + ggplot2::geom_crossbar(
        data = summ_r,
        ggplot2::aes(x = grp, y = ctr, ymin = ctr, ymax = ctr),
        inherit.aes = FALSE,
        colour = "black", width = 0.35, linewidth = 0.9)

    # Individual points — drawn LAST so they appear on top of bars
    p <- p + ggplot2::geom_point(
      position = ggplot2::position_jitter(width = jit_width, height = 0, seed = 42),
      size = pt_size, alpha = 0.85)

    if (!is.null(colors)) p <- p + ggplot2::scale_colour_manual(values = colors)
    if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
    if (!is.null(shapes)) p <- p + ggplot2::scale_shape_manual(values = shapes)
    p <- apply_plot_labels(p, input$tg_plot_title %||% "", input$tg_plot_subtitle %||% "",
                           input$tg_plot_title_align %||% "center", text_size,
                           input$tg_plot_xaxis_label %||% "", input$tg_plot_yaxis_label %||% "",
                           input$tg_plot_legend_label %||% "")
    p <- apply_legend_pos(p, legend_pos)
    p
  })

  output$tg_plot_rates <- renderPlot({
    tg_rates_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)

  # AUC/Treatment effects plot — reactive
  tg_auc_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    # Shared controls — remap colours/shapes to treatment-level names
    # so scale_*_manual values match the Treatment column, even when
    # rv$tg_plot_groups contains dose-composite names.
    ordered_grps <- ordered_tg_grps()
    colors     <- remap_colors_to_treatment(build_color_vector(input, "tg_color", ordered_grps))
    shapes     <- remap_shapes_to_treatment(build_shape_vector(input, "tg_shape", ordered_grps))
    center_bar <- input$tg_opt_center_bar %||% "mean"
    err_type   <- input$tg_opt_error_bar  %||% "SEM"
    pt_size    <- input$tg_opt_point_size %||% 2.5
    jit_width  <- input$tg_opt_jitter_width %||% 0.1
    text_size  <- input$tg_opt_text_size %||% 16
    legend_pos <- input$tg_opt_legend_pos %||% "right"

    # Helper: build a summary df (per group: center, lower, upper)
    make_summary <- function(df, grp_col, val_col) {
      groups <- if (is.factor(df[[grp_col]])) levels(df[[grp_col]]) else unique(as.character(df[[grp_col]]))
      do.call(rbind, lapply(groups, function(g) {
        x <- df[[val_col]][as.character(df[[grp_col]]) == g]
        x <- x[!is.na(x)]
        n <- length(x)
        ctr <- switch(center_bar,
          mean           = mean(x),
          median         = median(x),
          geometric_mean = exp(mean(log(pmax(x, 1e-9)))),
                           mean(x))
        half <- switch(err_type,
          SEM = sd(x) / sqrt(n),
          SD  = sd(x),
          CI  = qt(0.975, df = max(n - 1L, 1L)) * sd(x) / sqrt(n),
                0)
        data.frame(grp = g, ctr = ctr, lo = ctr - half, hi = ctr + half,
                   stringsAsFactors = FALSE)
      }))
    }

    # Helper: y-axis limits — floor at 0 (unless data dips below), 10% top buffer
    y_limits <- function(point_vals, lo_vals = NULL, hi_vals = NULL) {
      all_vals <- c(point_vals, lo_vals, hi_vals)
      all_vals <- all_vals[!is.na(all_vals)]
      ymin <- min(0, min(all_vals))
      ymax <- max(all_vals) * 1.10
      c(ymin, ymax)
    }

    # Helper: apply colours + labels + legend
    finish_plot <- function(p) {
      if (!is.null(colors)) p <- p + ggplot2::scale_colour_manual(values = colors)
      if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)
      if (!is.null(shapes)) p <- p + ggplot2::scale_shape_manual(values = shapes)
      p <- apply_plot_labels(p, input$tg_plot_title %||% "", input$tg_plot_subtitle %||% "",
                             input$tg_plot_title_align %||% "center", text_size = text_size,
                             xlab = input$tg_plot_xaxis_label %||% "",
                             ylab = input$tg_plot_yaxis_label %||% "",
                             legend_lab = input$tg_plot_legend_label %||% "")
      apply_legend_pos(p, legend_pos)
    }

    # ── Path A: AUC model ── per-animal AUC values available ──────────────
    if (!is.null(rv$tg_results$auc_analysis) &&
        !is.null(rv$tg_results$auc_analysis$individual)) {

      auc_ind   <- rv$tg_results$auc_analysis$individual
      grp_col_a <- intersect(c("Treatment", "Group"), colnames(auc_ind))[1]
      req(!is.na(grp_col_a))

      # Factor ordering
      if (length(ordered_grps) > 0) {
        auc_levels <- unique(gsub(" - Dose: ", " ", ordered_grps))
        auc_levels <- auc_levels[auc_levels %in% as.character(unique(auc_ind[[grp_col_a]]))]
        if (length(auc_levels) > 0)
          auc_ind[[grp_col_a]] <- factor(auc_ind[[grp_col_a]], levels = auc_levels)
      }

      summ <- make_summary(auc_ind, grp_col_a, "AUC")
      summ$grp <- factor(summ$grp, levels = if (is.factor(auc_ind[[grp_col_a]])) levels(auc_ind[[grp_col_a]]) else unique(as.character(auc_ind[[grp_col_a]])))

      p <- ggplot2::ggplot(auc_ind,
             ggplot2::aes(x = .data[[grp_col_a]], y = .data[["AUC"]],
                          colour = .data[[grp_col_a]], fill = .data[[grp_col_a]],
                          shape = .data[[grp_col_a]])) +
        ggplot2::theme_classic() +
        ggplot2::labs(x = "Treatment", y = "AUC", title = "Area Under the Curve") +
        ggplot2::theme(legend.title = ggplot2::element_blank())

      if (err_type != "none")
        p <- p + ggplot2::geom_errorbar(
          data = summ,
          ggplot2::aes(x = grp, ymin = lo, ymax = hi),
          inherit.aes = FALSE,
          colour = "black", width = 0.25, linewidth = 0.9)

      if (center_bar != "none")
        p <- p + ggplot2::geom_crossbar(
          data = summ,
          ggplot2::aes(x = grp, y = ctr, ymin = ctr, ymax = ctr),
          inherit.aes = FALSE,
          colour = "black", width = 0.25, linewidth = 0.9)

      p <- p + ggplot2::geom_point(
        position = ggplot2::position_jitter(width = jit_width, height = 0, seed = 42),
        size = pt_size, alpha = 0.85)

      ylims <- y_limits(auc_ind[["AUC"]], summ$lo, summ$hi)
      p <- p + ggplot2::coord_cartesian(ylim = ylims)
      return(finish_plot(p))
    }

    # ── Path B: LME4 model ── same visual style as Path A, using growth rates ──
    rates <- rv$tg_results$growth_rates
    if (!is.null(rates) && nrow(rates) > 0) {
      grp_col_b <- intersect(c("Treatment", "Group"), colnames(rates))[1]
      val_col_b <- intersect(c("growth_rate", "GrowthRate"), colnames(rates))[1]
      if (!is.na(grp_col_b) && !is.na(val_col_b)) {
        gr_order <- if (length(ordered_grps) > 0) {
          raw <- unique(gsub(" - Dose: ", " ", ordered_grps))
          raw[raw %in% as.character(unique(rates[[grp_col_b]]))]
        } else NULL
        if (!is.null(gr_order) && length(gr_order) > 0)
          rates[[grp_col_b]] <- factor(rates[[grp_col_b]], levels = gr_order)

        summ_b <- make_summary(rates, grp_col_b, val_col_b)
        summ_b$grp <- factor(summ_b$grp, levels = if (is.factor(rates[[grp_col_b]])) levels(rates[[grp_col_b]]) else unique(as.character(rates[[grp_col_b]])))

        # Growth rates are always log-volume slopes; back-transform to daily
        # fold change when the user has asked for original units.
        do_back_b <- isTRUE(input$tg_opt_backtransform)
        if (do_back_b) {
          rates[[val_col_b]] <- exp(rates[[val_col_b]])
          summ_b$ctr <- exp(summ_b$ctr)
          summ_b$lo  <- exp(summ_b$lo)
          summ_b$hi  <- exp(summ_b$hi)
        }
        y_label_b <- if (do_back_b) "Daily Fold Change" else "Growth Rate (log-volume / day)"

        p <- ggplot2::ggplot(rates,
               ggplot2::aes(x = .data[[grp_col_b]], y = .data[[val_col_b]],
                            colour = .data[[grp_col_b]], fill = .data[[grp_col_b]],
                            shape = .data[[grp_col_b]])) +
          ggplot2::theme_classic() +
          ggplot2::labs(x = "Treatment", y = y_label_b,
                        title = "Growth Rate Comparison") +
          ggplot2::theme(legend.title = ggplot2::element_blank())

        if (err_type != "none")
          p <- p + ggplot2::geom_errorbar(
            data = summ_b,
            ggplot2::aes(x = grp, ymin = lo, ymax = hi),
            inherit.aes = FALSE,
            colour = "black", width = 0.25, linewidth = 0.9)

        if (center_bar != "none")
          p <- p + ggplot2::geom_crossbar(
            data = summ_b,
            ggplot2::aes(x = grp, y = ctr, ymin = ctr, ymax = ctr),
            inherit.aes = FALSE,
            colour = "black", width = 0.25, linewidth = 0.9)

        p <- p + ggplot2::geom_point(
          position = ggplot2::position_jitter(width = jit_width, height = 0, seed = 42),
          size = pt_size, alpha = 0.85)

        ylims <- y_limits(rates[[val_col_b]], summ_b$lo, summ_b$hi)
        p <- p + ggplot2::coord_cartesian(ylim = ylims)
        return(finish_plot(p))
      }
    }

    # ── Path C: last resort — stored plot ─────────────────────────────────
    p <- rv$tg_results$auc_plot
    if (is.null(p)) return(NULL)
    finish_plot(p)
  })

  output$tg_plot_auc <- renderPlot({
    tg_auc_plot()
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)

  # ── Waterfall Plot ────────────────────────────────────────────────────────
  tg_wf_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    rates <- rv$tg_results$growth_rates
    req(!is.null(rates), nrow(rates) > 0)

    tryCatch({
      wf_data <- rates |>
        dplyr::filter(!is.na(.data$growth_rate)) |>
        dplyr::mutate(
          animal_label = if ("Cage" %in% colnames(rates))
                           paste0(.data$Cage, "/", .data$ID)
                         else
                           as.character(.data$ID),
          pct_per_day = (exp(.data$growth_rate) - 1) * 100
        ) |>
        dplyr::arrange(dplyr::desc(.data$pct_per_day)) |>
        dplyr::mutate(rank = dplyr::row_number())

      req(nrow(wf_data) > 0)

      trt_col      <- "Treatment"
      ordered_grps <- ordered_tg_grps()
      colors       <- remap_colors_to_treatment(build_color_vector(input, "tg_color", ordered_grps))
      text_size    <- input$tg_opt_text_size %||% 16

      p <- ggplot2::ggplot(
             wf_data,
             ggplot2::aes(x = .data$rank,
                          y = .data$pct_per_day,
                          fill = .data[[trt_col]])) +
        ggplot2::geom_col(width = 0.8, colour = "grey30", linewidth = 0.2) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.6) +
        { if (!is.null(colors)) ggplot2::scale_fill_manual(values = colors, name = trt_col)
          else ggplot2::scale_fill_discrete(name = trt_col) } +
        ggplot2::scale_x_continuous(
          breaks = wf_data$rank,
          labels = wf_data$animal_label,
          expand = c(0.01, 0)) +
        ggplot2::labs(
          x = "Animal",
          y = "Tumor growth rate (% change per day)") +
        ggplot2::theme_classic(base_size = text_size) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))

      p <- apply_legend_pos(p, input$tg_opt_legend_pos %||% "right")
      apply_plot_labels(p,
        input$tg_plot_title    %||% "",
        input$tg_plot_subtitle %||% "",
        input$tg_plot_title_align %||% "center",
        text_size,
        input$tg_plot_xaxis_label %||% "",
        input$tg_plot_yaxis_label %||% "",
        input$tg_plot_legend_label %||% "")

    }, error = function(e) {
      NULL
    })
  })

  output$tg_plot_waterfall <- renderPlot({
    p <- tg_wf_plot()
    req(!is.null(p))
    p
  }, width  = function() input$tg_plot_width  %||% 800,
     height = function() input$tg_plot_height %||% 550)

  # Effect size forest plot — rebuilt from raw comparison data each render
  # so group order changes take effect without re-running analysis
  tg_forest_plot <- reactive({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)

    # Pick best comparison source (Dunnett preferred)
    using_dunnett <- !is.null(rv$tg_results$reference_comparisons_dunnett)
    pc <- if (using_dunnett) rv$tg_results$reference_comparisons_dunnett
          else rv$tg_results$all_pairwise_comparisons
    if (is.null(pc)) return(NULL)

    text_size    <- input$tg_opt_text_size %||% 12
    ordered_grps <- ordered_tg_grps()
    grps_clean   <- unique(gsub(" - Dose: ", " ", ordered_grps))
    colors       <- remap_colors_to_treatment(build_color_vector(input, "tg_color", ordered_grps))
    shapes       <- remap_shapes_to_treatment(build_shape_vector(input, "tg_shape", ordered_grps))

    pc_df     <- as.data.frame(pc)
    label_col <- intersect(c("contrast", "Contrast", "comparison", "pair"), colnames(pc_df))[1]
    est_col   <- intersect(c("estimate", "Estimate"), colnames(pc_df))[1]
    se_col    <- intersect(c("SE"), colnames(pc_df))[1]
    # R20.D43: the backend's adjusted intervals first (0.31.0, R20.85), so the
    # interval and the adjusted p-value beside it agree.
    lower_col <- intersect(c("CI_Adjusted_Lower", "lower.CL", "Lower_CrI", "Lower", "lower_CL", "lower"), colnames(pc_df))[1]
    upper_col <- intersect(c("CI_Adjusted_Upper", "upper.CL", "Upper_CrI", "Upper", "upper_CL", "upper"), colnames(pc_df))[1]
    # R14.6: "p_value" held the ADJUSTED p-value on the GAM path and the RAW
    # one on the AUC path, so this sniff marked significance from unadjusted
    # p-values on AUC while using adjusted ones on GAM. Backend R14.4 now
    # supplies P_Value_Adjusted on every path; prefer it and keep the old
    # list only as a fallback for results from an older backend.
    p_col     <- intersect(c("P_Value_Adjusted", "p_adjusted", "p.value",
                             "p_value", "p.adj", "adj.p.value", "P_Value"),
                           colnames(pc_df))[1]
    if (is.na(label_col) || is.na(est_col)) return(NULL)

    if (is.na(lower_col) && !is.na(se_col)) {
      df_col <- intersect(c("df", "DF"), colnames(pc_df))[1]
      t_mult <- if (!is.na(df_col)) qt(0.975, pmax(as.numeric(pc_df[[df_col]]), 1)) else qnorm(0.975)
      pc_df$lower.CL <- pc_df[[est_col]] - t_mult * pc_df[[se_col]]
      pc_df$upper.CL <- pc_df[[est_col]] + t_mult * pc_df[[se_col]]
      lower_col <- "lower.CL"; upper_col <- "upper.CL"
    }

    p_label <- if (!is.na(p_col)) vapply(pc_df[[p_col]], format_pval, character(1)) else NULL
    ref_grp <- rv$tg_reference_used  # the reference of the displayed run

    contrast_labels <- as.character(pc_df[[label_col]])
    if (isTRUE(using_dunnett) && !is.null(ref_grp)) {
      suffix_pat <- paste0("\\s*[-\u2013]\\s*", gsub("([+()])", "\\\\\\1", ref_grp), "$")
      clean_tmp  <- sub(suffix_pat, "", contrast_labels)
      if (all(clean_tmp != contrast_labels)) contrast_labels <- clean_tmp
    }

    forest_df <- data.frame(
      Contrast = contrast_labels,
      Estimate = as.numeric(pc_df[[est_col]]),
      Lower    = if (!is.na(lower_col)) as.numeric(pc_df[[lower_col]]) else as.numeric(pc_df[[est_col]]),
      Upper    = if (!is.na(upper_col)) as.numeric(pc_df[[upper_col]]) else as.numeric(pc_df[[est_col]]),
      P_Label  = if (is.null(p_label)) "" else paste0("p=", p_label),
      is_ref   = FALSE,
      stringsAsFactors = FALSE
    )
    # R20.D57: all-pairs comparisons have no single reference row.
    if (!is.null(ref_grp) && isTRUE(using_dunnett)) {
      forest_df <- rbind(forest_df, data.frame(
        Contrast = paste0(ref_grp, "\n(Reference)"),
        Estimate = 0, Lower = NA_real_, Upper = NA_real_,
        P_Label = "", is_ref = TRUE,
        stringsAsFactors = FALSE
      ))
    }

    # Sort non-reference rows by user group order; reference always at bottom of plot
    non_ref <- forest_df[!forest_df$is_ref, , drop = FALSE]
    ref_row <- forest_df[ forest_df$is_ref, , drop = FALSE]
    if (length(grps_clean) > 0) {
      rank_by_grp <- match(non_ref$Contrast, grps_clean)
      rank_by_grp[is.na(rank_by_grp)] <- .Machine$integer.max
      non_ref <- non_ref[order(rank_by_grp), , drop = FALSE]
    } else {
      non_ref <- non_ref[order(-non_ref$Estimate), , drop = FALSE]
    }
    # Factor levels bottom-to-top: reference at bottom, first group-order item at top
    forest_df <- rbind(ref_row, non_ref)
    forest_df$Contrast <- factor(
      forest_df$Contrast,
      levels = c(as.character(ref_row$Contrast), rev(as.character(non_ref$Contrast)))
    )

    show_pval   <- isTRUE(input$tg_opt_forest_show_pval   %||% TRUE)
    show_effect <- isTRUE(input$tg_opt_forest_show_effect %||% FALSE)
    show_ci     <- isTRUE(input$tg_opt_forest_show_ci     %||% FALSE)
    any_ann     <- show_pval || show_effect || show_ci

    # Back-transform to fold change only when the model fit a log transform.
    # Sqrt back-transform of a difference of sqrt(volume) has no clean
    # biological interpretation, so we deliberately skip it.
    transform_used <- rv$tg_results$transform_used %||% "none"
    do_back        <- isTRUE(input$tg_opt_backtransform) && transform_used == "log"

    if (do_back) {
      forest_df$Estimate <- exp(forest_df$Estimate)
      forest_df$Lower    <- exp(forest_df$Lower)
      forest_df$Upper    <- exp(forest_df$Upper)
      # Reference row's "Estimate = 0" on log scale becomes 1 (no effect) on
      # fold-change scale. is_ref rows keep NA CIs so they aren't drawn.
      forest_df$Estimate[forest_df$is_ref] <- 1
      non_ref <- forest_df[!forest_df$is_ref, , drop = FALSE]
    }

    # Effect-size label format depends on scale
    eff_prefix <- if (do_back) "FC: " else "est: "

    # Build per-row annotation label from whichever elements are enabled
    forest_df$Ann_Label <- vapply(seq_len(nrow(forest_df)), function(j) {
      if (forest_df$is_ref[j]) return("")
      parts <- character(0)
      if (show_effect)
        parts <- c(parts, paste0(eff_prefix, signif(forest_df$Estimate[j], 3)))
      if (show_ci && !is.na(forest_df$Lower[j]))
        parts <- c(parts,
          paste0("[", signif(forest_df$Lower[j], 3), ", ", signif(forest_df$Upper[j], 3), "]"))
      if (show_pval && nzchar(forest_df$P_Label[j]))
        parts <- c(parts, forest_df$P_Label[j])
      paste(parts, collapse = "\n")
    }, character(1L))

    # Build per-contrast colour/shape maps (contrast labels correspond to group names)
    contrast_lvls <- setdiff(as.character(levels(forest_df$Contrast)),
                             as.character(ref_row$Contrast))
    forest_colours <- if (!is.null(colors)) {
      vals <- ifelse(contrast_lvls %in% names(colors), unname(colors[contrast_lvls]), "steelblue")
      setNames(vals, contrast_lvls)
    } else setNames(rep("steelblue", length(contrast_lvls)), contrast_lvls)
    forest_shapes <- if (!is.null(shapes)) {
      vals <- ifelse(contrast_lvls %in% names(shapes), unname(shapes[contrast_lvls]), 16L)
      setNames(as.integer(vals), contrast_lvls)
    } else setNames(rep(16L, length(contrast_lvls)), contrast_lvls)

    # x-axis limits & breaks differ between linear (difference) and log10
    # (fold-change) scales. On fold-change scale, range is centred on 1 in
    # log space so 0.5x and 2x are equidistant from 1.
    if (do_back) {
      all_vals  <- c(non_ref$Estimate, non_ref$Lower, non_ref$Upper)
      all_vals  <- all_vals[is.finite(all_vals) & all_vals > 0]
      max_log10 <- max(ceiling(max(abs(log10(all_vals)), na.rm = TRUE)), 1L)
      # Widen right margin when annotations are shown so text is not clipped
      right_pad <- if (any_ann) 0.5 else 0.1
      x_lim     <- c(10^(-max_log10), 10^(max_log10 + right_pad))
      x_breaks  <- 10^seq(-max_log10, max_log10, by = 1)
      x_lab     <- if (isTRUE(using_dunnett))
                     "Fold Change (treatment / reference)"
                   else "Fold Change"
    } else {
      all_vals  <- c(non_ref$Estimate, non_ref$Lower, non_ref$Upper)
      max_abs   <- ceiling(max(abs(all_vals), na.rm = TRUE))
      # Widen right margin when annotations are shown so text is not clipped
      x_lim     <- c(-max_abs, max_abs * if (any_ann) 1.55 else 1.20)
      x_breaks  <- pretty(c(-max_abs, max_abs), n = 6)
      x_breaks  <- x_breaks[x_breaks >= -max_abs & x_breaks <= max_abs]
      scale_tag <- scale_suffix(transform_used)
      x_lab     <- if (isTRUE(using_dunnett))
                     paste0("Mean Difference (treatment \u2212 reference)", scale_tag)
                   else paste0("Mean Difference", scale_tag)
    }

    null_xint <- if (do_back) 1 else 0

    p <- ggplot2::ggplot(forest_df, ggplot2::aes(x = Estimate, y = Contrast)) +
      ggplot2::geom_vline(xintercept = null_xint, linetype = "dashed", colour = "red", linewidth = 0.6) +
      ggplot2::geom_errorbarh(
        data = ~ subset(., !is_ref),
        ggplot2::aes(xmin = Lower, xmax = Upper, colour = Contrast),
        height = 0.25) +
      ggplot2::geom_point(data = ~ subset(., !is_ref),
        ggplot2::aes(colour = Contrast, shape = Contrast),
        size = input$tg_opt_forest_pt_size %||% 3) +
      ggplot2::geom_point(data = ~ subset(., is_ref),
        size = input$tg_opt_forest_pt_size %||% 3, colour = "grey50", shape = 18) +
      ggplot2::geom_text(
        data = ~ subset(., !is_ref & nzchar(Ann_Label) & !is.na(Upper)),
        ggplot2::aes(label = Ann_Label, x = Upper),
        hjust = -0.15, vjust = 0.5,
        size = (input$tg_opt_forest_ann_size %||% 12) / ggplot2::.pt,
        lineheight = 0.9) +
      ggplot2::scale_colour_manual(values = forest_colours, guide = "none") +
      ggplot2::scale_shape_manual(values = forest_shapes, guide = "none") +
      (if (do_back)
         ggplot2::scale_x_log10(limits = x_lim, breaks = x_breaks)
       else
         ggplot2::scale_x_continuous(limits = x_lim, breaks = x_breaks)) +
      ggplot2::labs(x = x_lab, y = NULL, title = "Effect Size Forest Plot",
                    caption = forest_ci_note(pc_df, lower_col)) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
    apply_plot_labels(p,
      input$tg_plot_title    %||% "",
      input$tg_plot_subtitle %||% "",
      input$tg_plot_title_align %||% "center",
      text_size,
      input$tg_plot_xaxis_label  %||% "",
      input$tg_plot_yaxis_label  %||% "",
      input$tg_plot_legend_label %||% "")
  })

  output$tg_plot_effect_forest <- renderPlot({
    p <- tg_forest_plot()
    if (!is.null(p)) {
      p
    } else {
      plot(1, 1, type="n", main="Effect size forest plot not available\n(requires pairwise comparison data)", xlab="", ylab="")
    }
  }, width = function() input$tg_plot_width %||% 800,
     height = function() input$tg_plot_height %||% 550)

  # Diagnostics tab: Q-Q plot
  output$tg_diag_qq <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    if (isTRUE(rv$tg_results$model_type_used == "bayes")) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Frequentist diagnostics not applicable in Bayesian mode.\nSee the Posterior Check tab for model validation.")
      return(invisible(NULL))
    }
    if (!is.null(rv$tg_results$diag_qq_plot)) {
      return(rv$tg_results$diag_qq_plot)
    }
    plot(1, 1, type = "n", main = "Q-Q plot not available", xlab = "", ylab = "")
  })

  # Diagnostics tab: Residuals vs Fitted
  output$tg_diag_resid_fitted <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    if (isTRUE(rv$tg_results$model_type_used == "bayes")) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Frequentist diagnostics not applicable in Bayesian mode.\nSee the Posterior Check tab for model validation.")
      return(invisible(NULL))
    }
    if (!is.null(rv$tg_results$diag_resid_fitted_plot)) {
      return(rv$tg_results$diag_resid_fitted_plot)
    }
    plot(1, 1, type = "n", main = "Residuals vs Fitted not available", xlab = "", ylab = "")
  })

  # Scale-Location (heteroscedasticity check). Falls back to a "not
  # available" plot for backwards compat with old result lists.
  output$tg_diag_scale_location <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    if (isTRUE(rv$tg_results$model_type_used == "bayes")) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Not applicable in Bayesian mode.")
      return(invisible(NULL))
    }
    if (!is.null(rv$tg_results$diag_scale_location_plot)) {
      return(rv$tg_results$diag_scale_location_plot)
    }
    plot(1, 1, type = "n", main = "Scale-Location plot not available",
         xlab = "", ylab = "")
  })

  # Random-effects Q-Q. Only meaningful for LME4 fits with random
  # intercept (or slope); AUC path has no random effects.
  output$tg_diag_re_qq <- renderPlot({
    req(rv$tg_results)
    if (isTRUE(rv$tg_results$error)) return(NULL)
    if (isTRUE(rv$tg_results$model_type_used == "bayes")) {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Not applicable in Bayesian mode.")
      return(invisible(NULL))
    }
    if (!is.null(rv$tg_results$diag_re_qq_plot)) {
      return(rv$tg_results$diag_re_qq_plot)
    }
    plot(1, 1, type = "n",
         main = "Random-effects Q-Q not available (AUC path has no random effects)",
         xlab = "", ylab = "")
  })

  # CODE_REVIEW.md DIAGNOSTICS gap (13) — LMM influence diagnostics.
  output$tg_show_influence <- shiny::reactive({
    !is.null(rv$tg_results) &&
      !isTRUE(rv$tg_results$error) &&
      (!is.null(rv$tg_results$diag_cooks_distance) ||
       !is.null(rv$tg_results$diag_dfbetas))
  })
  shiny::outputOptions(output, "tg_show_influence", suspendWhenHidden = FALSE)

  output$tg_cooks_distance <- DT::renderDT({
    req(rv$tg_results)
    cd <- rv$tg_results$diag_cooks_distance
    if (is.null(cd) || nrow(cd) == 0L)
      return(DT::datatable(data.frame(Message = "Cook's distance not available")))
    DT::datatable(format_df_for_display(as.data.frame(cd)),
              options = list(scrollX = TRUE, pageLength = 10),
              rownames = FALSE)
  })

  output$tg_dfbetas <- DT::renderDT({
    req(rv$tg_results)
    db <- rv$tg_results$diag_dfbetas
    if (is.null(db) || nrow(db) == 0L)
      return(DT::datatable(data.frame(Message = "DFBETAS not available")))
    DT::datatable(format_df_for_display(as.data.frame(db)),
              options = list(scrollX = TRUE, pageLength = 10),
              rownames = FALSE)
  })

  list(
    growth    = tg_growth_plot,
    rates     = tg_rates_plot,
    auc       = tg_auc_plot,
    waterfall = tg_wf_plot,
    forest    = tg_forest_plot
  )
}
