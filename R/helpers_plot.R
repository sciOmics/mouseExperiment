# Plot aesthetics: colours, shapes, group ordering and export
#
# Split out of the former `helpers.R` under CODE_REVIEW.md D.1.
#
# Lower stakes than the data and format helpers -- a wrong colour misleads
# nobody -- with one exception worth stating: `remap_colors_to_treatment()` and
# `remap_shapes_to_treatment()` re-key aesthetics onto treatment names, so a
# mismatch there mislabels which arm a curve belongs to. `get_ordered_groups()`
# feeds the reference-group choice, which is a statistical input, not a display
# one.



#' Apply custom title, subtitle, text size, and alignment to a ggplot
#'
#' When title/subtitle are empty strings the existing labs() are preserved
#' (i.e. the hardcoded plot title is not removed).
#'
#' @param p A ggplot2 object.
#' @param title Character; if non-empty overrides the plot title.
#' @param subtitle Character; if non-empty adds a subtitle.
#' @param align One of "left", "center", "right" — alignment for title/subtitle.
#' @param text_size Numeric; base font size in pt.
#' @param xlab Character; if non-empty overrides the X-axis label.
#' @param ylab Character; if non-empty overrides the Y-axis label.
#' @param legend_lab Character; if non-empty overrides the legend title.
#' @return A ggplot2 object with the theme layers applied.
#' @keywords internal
apply_plot_labels <- function(p, title = "", subtitle = "", align = "center",
                              text_size = 12, xlab = "", ylab = "", legend_lab = "") {
  # F9: the exported style_plot(), which the R script calls too.
  style_plot(p, title = title, subtitle = subtitle,
             title_align = if (align %in% c("left", "center", "right")) align else "center",
             text_size = text_size, x_label = xlab, y_label = ylab, legend_title = legend_lab)
}

#' Build a named colour vector from colourInput values
#'
#' Reads shiny `input[[paste0(prefix, "_", make.names(group))]]` for each group
#' and returns a named character vector suitable for ggplot2's
#' `scale_color_manual(values = ...)` / `scale_fill_manual(values = ...)`.
#' Entries whose input value is NULL or empty are dropped.
#'
#' @param input  The Shiny `input` list.
#' @param prefix Character prefix used when naming the colourInput widgets,
#'   e.g. `"tg_color"` → input ID `"tg_color_<make.names(group)>"`.
#' @param groups Character vector of group names.
#' @return Named character vector of hex colours, or NULL if all are missing.
#' @keywords internal
build_color_vector <- function(input, prefix, groups) {
  if (is.null(groups) || length(groups) == 0L) return(NULL)
  colors <- vapply(groups, function(g) {
    val <- input[[paste0(prefix, "_", make.names(g))]]
    if (is.null(val) || !nzchar(val)) NA_character_ else val
  }, character(1L))
  # Return NULL unless ALL groups have values — partial coverage causes
  # ggplot2 scale_colour_manual to map missing groups as grey/invisible.
  if (any(is.na(colors))) return(NULL)
  colors
}

#' Remap colour vector names from dose-composite to treatment-group names
#'
#' When groups are stored as \code{"Treatment - Dose: X"}, colour names
#' include the \code{" - Dose: "} separator.  Plot data (growth_rates, AUC) uses
#' the \code{TreatmentGroup} column which pastes treatment and dose with a
#' single space (e.g. \code{"Control 0"}).  This function converts
#' \code{" - Dose: "} to \code{" "} and deduplicates (keeping first).
#'
#' @param colors Named character vector (or NULL).
#' @return Named character vector with treatment-group names, or NULL.
#' @keywords internal
remap_colors_to_treatment <- function(colors) {
  if (is.null(colors)) return(NULL)
  nms <- gsub(" - Dose: ", " ", names(colors))
  dups <- duplicated(nms)
  out <- colors[!dups]
  names(out) <- nms[!dups]
  out
}

#' Remap shape vector names from dose-composite to treatment-group names
#'
#' Same logic as \code{remap_colors_to_treatment} but for the integer shape vector.
#'
#' @param shapes Named integer vector (or NULL).
#' @return Named integer vector with treatment-level names, or NULL.
#' @keywords internal
remap_shapes_to_treatment <- function(shapes) {
  if (is.null(shapes)) return(NULL)
  nms <- gsub(" - Dose: ", " ", names(shapes))
  dups <- duplicated(nms)
  out <- shapes[!dups]
  names(out) <- nms[!dups]
  out
}

#' Build a named integer vector of per-group point shapes from Shiny inputs
#'
#' Reads \code{input[[paste0(prefix, "_", make.names(group))]]} for each group.
#' Expected values are pch integers as character strings ("16", "15", etc.).
#'
#' @param input  Shiny \code{input} object.
#' @param prefix Input ID prefix, e.g. \code{"tg_shape"}.
#' @param groups Character vector of group names.
#' @return Named integer vector of pch codes, or NULL if none set.
#' @keywords internal
build_shape_vector <- function(input, prefix, groups) {
  if (is.null(groups) || length(groups) == 0L) return(NULL)
  shapes <- vapply(groups, function(g) {
    val <- input[[paste0(prefix, "_", make.names(g))]]
    if (is.null(val) || !nzchar(val)) NA_integer_ else as.integer(val)
  }, integer(1L))
  # Return NULL unless ALL groups have values — partial coverage causes
  # ggplot2 scale_shape_manual to map missing groups as NA shape (invisible).
  if (any(is.na(shapes))) return(NULL)
  shapes
}

#' Generate n evenly-spaced ggplot2-style default colours for n groups
#'
#' @param n Integer number of colours.
#' @return Character vector of hex colour strings.
#' @keywords internal
default_group_colors <- function(n) {
  if (n <= 0L) return(character(0L))
  scales::hue_pal()(n)
}

#' Named palette presets (colorblind-friendly options)
#'
#' Returns a named character vector of hex colours for the requested palette.
#' Falls back to the ggplot2 default hue palette when \code{name} is
#' \code{"default"} or unrecognised.
#'
#' @param name Character palette name: \code{"default"}, \code{"viridis"},
#'   \code{"okabe_ito"}, \code{"set2"}, or \code{"wes_anderson"}.
#' @param n    Number of colours required.
#' @return Character vector of hex colour strings of length \code{n}.
#' @keywords internal
palette_colors <- function(name, n) {
  if (n <= 0L) return(character(0L))
  switch(name,
    viridis  = substring(scales::viridis_pal(option = "D")(n), 1L, 7L),
    okabe_ito = {
      oi <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
              "#0072B2", "#D55E00", "#CC79A7", "#999999")
      rep_len(oi, n)
    },
    set2 = {
      s2 <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
              "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
      rep_len(s2, n)
    },
    wes_anderson = {
      wa <- c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")
      rep_len(wa, n)
    },
    default_group_colors(n)
  )
}

#' Set up an observer that applies palette presets to group colour pickers
#'
#' Call this once inside a module server.
#' When the user selects a non-custom preset, every per-group \code{colourInput}
#' is updated to the preset colour.  Selecting \code{"custom"} is a no-op so
#' the user can freely pick colours after choosing a starting palette.
#'
#' @param input   Shiny \code{input} object (module-scoped).
#' @param session Shiny \code{session} object (module-scoped).
#' @param prefix  Same \code{prefix} passed to \code{group_order_and_color_ui}.
#' @param groups_rv A reactive expression returning the current group names.
#' @keywords internal
observe_palette_preset <- function(input, session, prefix, groups_rv) {
  pal_id <- paste0(prefix, "_palette")
  observeEvent(input[[pal_id]], {
    pal_name <- input[[pal_id]]
    if (is.null(pal_name) || pal_name == "custom") return()
    groups <- groups_rv()
    if (is.null(groups) || length(groups) == 0L) return()
    cols <- palette_colors(pal_name, length(groups))
    for (i in seq_along(groups)) {
      color_id <- paste0(prefix, "_color_", make.names(groups[i]))
      colourpicker::updateColourInput(session, color_id, value = cols[i])
    }
  }, ignoreInit = TRUE)
}


#' Return groups in the user-specified drag-and-drop order
#'
#' Reads \code{input[[paste0(prefix, "_order")]]} and validates it.
#' Falls back to \code{groups} (original order) if the input is absent or
#' does not contain exactly the same set of labels.
#'
#' @param input  Shiny \code{input} object.
#' @param prefix Input-ID prefix used when calling \code{group_order_ui}.
#' @param groups Character vector of group names in their default order.
#' @return Character vector of group names in the current user-selected order.
#' @keywords internal
get_ordered_groups <- function(input, prefix, groups) {
  if (is.null(groups)) return(groups)
  ord <- input[[paste0(prefix, "_order")]]
  if (!is.null(ord)) ord <- as.character(unlist(ord))
  if (!is.null(ord) &&
      length(ord) == length(groups) &&
      isTRUE(setequal(ord, groups))) {
    ord
  } else {
    groups
  }
}

#' Apply a legend position value (including in-plot coordinates) to a ggplot
#'
#' Translates the dashboard's legend position selector values into proper
#' ggplot2 theme settings, including the coordinate-based in-plot options.
#'
#' @param p   A ggplot2 object.
#' @param pos Character; one of "none", "top", "bottom", "left",
#'   "right", "inplot_tl" (top-left), "inplot_tr" (top-right),
#'   "inplot_bl" (bottom-left), or "inplot_br" (bottom-right).
#' @return The modified ggplot2 object.
#' @keywords internal
apply_legend_pos <- function(p, pos) {
  if (is.null(pos)) return(p)
  style_plot(p, legend_position = pos)   # F9
}

#' Render a drag-and-drop rank list with embedded color pickers
#'
#' Combines group ordering (drag-and-drop) with color selection in a single
#' control. Groups can be reordered by dragging, and colors are adjusted via
#' inline colour pickers. The order is stored in \code{input[[paste0(prefix, "_order")]]}
#' and colors are read from \code{input[[paste0(prefix, "_color_", make.names(group))]]}
#' for each group.
#'
#' @param prefix    Input-ID prefix (e.g. \code{"tg"}).
#' @param groups    Character vector of group names in their default order.
#' @return A \code{tagList} containing HTML structure and colour inputs, or
#'   \code{NULL} when fewer than 2 groups.
#' @keywords internal
group_order_and_color_ui <- function(prefix, groups, ns = NS(NULL)) {
  if (is.null(groups) || length(groups) < 2L) return(NULL)

  defaults <- default_group_colors(length(groups))
  list_id  <- ns(paste0(prefix, "_drag_list"))
  input_id <- ns(paste0(prefix, "_order"))
  pal_id   <- ns(paste0(prefix, "_palette"))

  # Closed pch shapes available for per-group selection.
  # Open shapes (0, 1, 2, ...) are reserved for extrapolated points.
  # The symbols alone: with their names ("Triangle \u2191") the menu needed
  # 125 px and overflowed the sidebar.
  shape_choices <- c(
    "\u25cf" = "16",   # circle
    "\u25a0" = "15",   # square
    "\u25b2" = "17",   # triangle, point up
    "\u25bc" = "25",   # triangle, point down
    "\u25c6" = "18"    # diamond
  )
  default_shapes <- c(16L, 15L, 17L, 25L, 18L)

  # Each group is two lines: drag handle + full group name, then the colour
  # swatch + shape menu. One line overflowed the sidebar: at a 1400 px window
  # the panel is 209 px wide, the row needed 275 px, and the shape menus ran
  # out of the box while the names were cut to "Drug_..." (2026-10-22).
  # data-value holds the exact group name string — read directly by JS,
  # no innerText parsing, no custom Shiny input type.
  items <- lapply(seq_along(groups), function(i) {
    g        <- groups[i]
    color_id <- ns(paste0(prefix, "_color_", make.names(g)))
    shape_id <- ns(paste0(prefix, "_shape_", make.names(g)))
    def_shape <- as.character(default_shapes[((i - 1L) %% length(default_shapes)) + 1L])
    shiny::tags$li(
      `data-value` = g,
      class = "group-style-row",
      style = paste0(
        "padding:5px 8px; margin-bottom:4px;",
        "background:#f8f9fa; border:1px solid #dee2e6;",
        "border-radius:4px; cursor:default;"
      ),
      shiny::tags$div(
        style = "display:flex; align-items:baseline; gap:6px;",
        shiny::tags$span(
          class = "drag-handle",
          style = "cursor:grab; color:#adb5bd; font-size:14px; user-select:none; flex-shrink:0;",
          "\u2630"   # grip / hamburger icon
        ),
        shiny::tags$span(
          style = "flex:1; min-width:0; font-size:13px; overflow-wrap:anywhere;",
          g
        )
      ),
      shiny::tags$div(
        class = "group-style-controls",
        style = "display:flex; align-items:center; gap:6px; margin-top:4px; padding-left:20px;",
        shiny::tags$div(
          style = "flex:0 0 44px;",
          colourpicker::colourInput(
            inputId          = color_id,
            label            = NULL,
            value            = defaults[i],
            showColour       = "background",
            palette          = "square",
            allowTransparent = FALSE,
            width            = "44px"
          )
        ),
        shiny::tags$div(
          style = "flex:1; min-width:0;",
          shiny::selectInput(
            inputId  = shape_id,
            label    = NULL,
            choices  = shape_choices,
            selected = def_shape,
            width    = "100%"
          )
        )
      )
    )
  })

  # JS helper: read data-value attributes from the sortable <li> elements
  # and push the current order to Shiny.  Uses getElementById for reliability
  # (avoids depending on SortableJS `this.el` inside async callbacks).
  send_order_js <- htmlwidgets::JS(sprintf(
    "function(evt) {
       var el = document.getElementById('%s');
       if (el && typeof Shiny !== 'undefined' && Shiny.setInputValue) {
         var order = Array.prototype.slice.call(el.children).map(
           function(li) { return li.getAttribute('data-value'); }
         );
         Shiny.setInputValue('%s', order, {priority: 'event'});
       }
     }",
    list_id, input_id
  ))

  # Fire once after the sortable is rendered so the initial order reaches Shiny
  # (onLoad is NOT a standard SortableJS event and is silently ignored).
  init_js <- sprintf(
    "$(function(){
       var el = document.getElementById('%s');
       if (el && typeof Shiny !== 'undefined' && Shiny.setInputValue) {
         var order = Array.prototype.slice.call(el.children).map(
           function(li) { return li.getAttribute('data-value'); }
         );
         Shiny.setInputValue('%s', order, {priority: 'event'});
       }
     });",
    list_id, input_id
  )

  shiny::tagList(
    shiny::selectInput(
      pal_id, shiny::tags$span(shiny::icon("eye"), " Color Palette Presets"),
      choices  = c("Default" = "custom",
                   "Viridis" = "viridis",
                   "Okabe-Ito" = "okabe_ito",
                   "Set2" = "set2",
                   "Wes Anderson" = "wes_anderson"),
      selected = "custom",
      width    = "100%"
    ),
    shiny::tags$p(
      style = "font-size:11px; color:#6c757d; margin-bottom:6px;",
      shiny::icon("grip-vertical"), " Drag the \u2630 handle to reorder"
    ),
    shiny::tags$ul(
      id    = list_id,
      style = "padding:0; margin:0; list-style:none;",
      items
    ),
    sortable::sortable_js(
      css_id  = list_id,
      options = sortable::sortable_options(
        handle    = ".drag-handle",
        animation = 150L,
        onSort    = send_order_js
      )
    ),
    shiny::tags$script(shiny::HTML(init_js))
  )
}

#' Hazard-ratio forest plot for the Survival tab
#'
#' Shared by the on-screen plot and its PNG download. CODE_REVIEW.md R20.D19:
#' this was a local function inside `survival_server()`, so the download handler
#' in `mod_survival_downloads.R` -- a separate top-level function -- could not
#' find it, and "Download Plot" returned HTTP 500 while the plot showed fine.
#'
#' @param res Cox results table with `Group`, `HR`, `CI_Lower`, `CI_Upper`,
#'   `P_Value`.
#' @param ref_grp Reference group name.
#' @param pt_color Point and error-bar colour.
#' @param show_hr,show_ci,show_p Which annotations to print beside each row.
#' @return A ggplot, or NULL when there is nothing to plot.
#' @keywords internal
build_surv_forest_plot <- function(res, ref_grp, pt_color, show_hr, show_ci, show_p,
                                   measure = "Hazard Ratio", abbrev = "HR") {
  plot_df <- res[!is.na(res$HR), ]
  if (!is.null(ref_grp)) plot_df <- plot_df[plot_df$Group != ref_grp, ]
  if (nrow(plot_df) == 0) return(NULL)

  # Add reference group row at HR = 1 (no CI)
  ref_row <- plot_df[1L, , drop = FALSE]
  ref_row[1L, ] <- NA
  ref_row$Group    <- paste0(ref_grp, "\n(Reference)")
  ref_row$HR       <- 1
  ref_row$CI_Lower <- NA_real_
  ref_row$CI_Upper <- NA_real_
  plot_df <- rbind(plot_df, ref_row)

  # Ensure numeric columns stay numeric after rbind
  for (nc in c("HR", "CI_Lower", "CI_Upper", "P_Value")) {
    if (nc %in% colnames(plot_df)) plot_df[[nc]] <- as.numeric(plot_df[[nc]])
  }
  plot_df$is_ref <- plot_df$Group == paste0(ref_grp, "\n(Reference)")

  # Symmetric log10 axis (exclude reference row)
  all_hr_vals <- c(plot_df$HR[!plot_df$is_ref],
                   plot_df$CI_Lower[!plot_df$is_ref],
                   plot_df$CI_Upper[!plot_df$is_ref])
  all_hr_vals <- all_hr_vals[!is.na(all_hr_vals) & all_hr_vals > 0]
  max_log10 <- max(ceiling(max(abs(log10(all_hr_vals)), na.rm = TRUE)), 1L)
  x_lim    <- c(10^(-max_log10), 10^max_log10)
  x_breaks <- 10^seq(-max_log10, max_log10, by = 1)
  x_labels <- ifelse(x_breaks >= 1,
                     as.character(round(x_breaks, 0)),
                     as.character(x_breaks))

  # Build per-row annotation labels (vapply avoids the apply(sapply) dim error)
  hr_vals  <- as.numeric(plot_df$HR)
  ci_lo    <- as.numeric(plot_df$CI_Lower)
  ci_hi    <- as.numeric(plot_df$CI_Upper)
  p_vals   <- as.numeric(plot_df$P_Value)
  plot_df$label_text <- vapply(seq_len(nrow(plot_df)), function(i) {
    if (isTRUE(plot_df$is_ref[i])) return("")
    parts <- character(0)
    if (show_hr && !is.na(hr_vals[i]))  parts <- c(parts, paste0(abbrev, ": ", round(hr_vals[i], 2)))
    if (show_ci && !is.na(ci_lo[i]))    parts <- c(parts, paste0("[", round(ci_lo[i], 2), "\u2013", round(ci_hi[i], 2), "]"))
    if (show_p  && !is.na(p_vals[i]))   parts <- c(parts, paste0("p = ", format_pval(p_vals[i])))
    paste(parts, collapse = "  ")
  }, character(1))

  ggplot2::ggplot(plot_df,
                 ggplot2::aes(x = HR, y = stats::reorder(Group, HR),
                              xmin = CI_Lower, xmax = CI_Upper)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    ggplot2::geom_errorbarh(
      data = ~ subset(., !is_ref),
      height = 0.25, linewidth = 0.8, color = pt_color
    ) +
    ggplot2::geom_point(
      data = ~ subset(., !is_ref), size = 5, color = pt_color
    ) +
    ggplot2::geom_point(
      data = ~ subset(., is_ref), size = 5, color = "grey50", shape = 18
    ) +
    { if (any(nzchar(plot_df$label_text[!plot_df$is_ref])))
        ggplot2::geom_text(
          data = ~ subset(., !is_ref),
          ggplot2::aes(label = label_text),
          vjust = 3, hjust = -0.1, size = 4, color = "grey30")
      else list() } +
    ggplot2::scale_x_continuous(
      trans  = "log10",
      limits = x_lim,
      breaks = x_breaks,
      labels = x_labels
    ) +
    ggplot2::labs(
      title = paste0("Forest Plot: ", measure, "s vs ", ref_grp),
      # Plain "log10": the PNG device font has no subscript digits, so the
      # label rendered as "log" plus two boxes, on screen and in the download.
      x     = paste(measure, "(log10 scale)"),
      y     = NULL
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.text.y    = ggplot2::element_text(size = 11),
      legend.position = "none"
    )
}

#' The over-time synergy TGIs and their intervals, one row per arm and day
#'
#' CODE_REVIEW.md R20.D37: the Synergy TGI plot draws the backend's TGIs (the
#' endpoint model on each evaluable day) instead of recomputing TGI from the
#' animals present each day.
#'
#' @param ss `synergy_summary` of `analyze_drug_synergy_over_time()`.
#' @param grps List with `drug_a`, `drug_b` and `combo` labels.
#' @return Data frame: Day, Treatment (factor), mean_tgi, lo, hi.
#' @keywords internal
syn_tgi_long <- function(ss, grps) {
  long <- function(label, est, lo, hi) data.frame(
    Day = ss$Time_Point, Treatment = label, mean_tgi = ss[[est]],
    lo = if (lo %in% names(ss)) ss[[lo]] else NA_real_,
    hi = if (hi %in% names(ss)) ss[[hi]] else NA_real_,
    stringsAsFactors = FALSE)
  out <- rbind(long(grps$drug_a, "TGI_Drug_A", "TGI_Drug_A_Lower", "TGI_Drug_A_Upper"),
               long(grps$drug_b, "TGI_Drug_B", "TGI_Drug_B_Lower", "TGI_Drug_B_Upper"),
               long(grps$combo,  "TGI_Combo",  "TGI_Combo_Lower",  "TGI_Combo_Upper"))
  out$Treatment <- factor(out$Treatment, levels = c(grps$drug_a, grps$drug_b, grps$combo))
  out
}

#' The KM plot's p-value, labelled -- R20.D29
#'
#' survminer prints the omnibus log-rank p as a bare "p = ...", next to Cox
#' p-values that are adjusted for multiple comparisons. It now says what it is.
#'
#' @param data Data with `Time`, `Event` and `Treatment`.
#' @return A label for `ggsurvplot(pval = )`, or TRUE when it cannot be computed.
#' @keywords internal
km_pval_label <- function(data) {
  sd <- tryCatch(survival::survdiff(survival::Surv(Time, Event) ~ Treatment, data = data),
                 error = function(e) NULL)
  if (is.null(sd)) return(TRUE)
  p <- stats::pchisq(sd$chisq, length(sd$n) - 1L, lower.tail = FALSE)
  paste0("Log-rank, all arms (unadjusted): p ", format_pval_label(p))
}

format_pval_label <- function(p) {
  if (!is.finite(p)) return("= NA")
  if (p < 1e-4) "< 0.0001" else paste("=", formatC(p, digits = 2, format = "g"))
}

