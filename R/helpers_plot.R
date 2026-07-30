#' Plot aesthetics: colours, shapes, group ordering and export
#'
#' Split out of the former `helpers.R` under CODE_REVIEW.md D.1.
#'
#' Lower stakes than the data and format helpers -- a wrong colour misleads
#' nobody -- with one exception worth stating: `remap_colors_to_treatment()` and
#' `remap_shapes_to_treatment()` re-key aesthetics onto treatment names, so a
#' mismatch there mislabels which arm a curve belongs to. `get_ordered_groups()`
#' feeds the reference-group choice, which is a statistical input, not a display
#' one.

#' Safely generate a plot, returning an error placeholder on failure
#'
#' Wraps a plot expression in tryCatch. On success returns the plot object.
#' On failure logs via structured logger and returns a simple ggplot with the
#' error text, so users see a message instead of a blank panel.
#'
#' @param expr    Expression that produces a ggplot / base plot grob.
#' @param label   Short label for the plot (used in the error message).
#' @return A plot object.
#' @keywords internal
safe_plot <- function(expr, label = "Plot") {
  tryCatch(
    expr,
    error = function(e) {
      msg <- paste0(label, " generation failed: ", conditionMessage(e))
      log_caught_error("plot", label, e)
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                          size = 5, hjust = 0.5) +
        ggplot2::theme_void()
    }
  )
}

#' Plot export dimension controls
#'
#' Generates Width (px), Height (px), and Resolution (DPI) numeric inputs that
#' are placed at the bottom of each analysis sidebar. Download handlers read
#' these with \code{input[[paste0(prefix, "_plot_width")]]} etc.
#'
#' @param prefix Input ID prefix, e.g. \code{"tg"}, \code{"surv"}.
#' @keywords internal
plot_export_ui <- function(prefix) {
  tagList(
    hr(style = "margin-top: 14px; margin-bottom: 8px;"),
    p(style = "font-weight: 600; margin-bottom: 4px; font-size: 0.88em;",
      "\U0001F4BE  Plot Export Settings"),
    fluidRow(
      column(6,
        numericInput(paste0(prefix, "_plot_width"),  "Width (px)",
                     value = 800, min = 200, max = 5000, step = 50)
      ),
      column(6,
        numericInput(paste0(prefix, "_plot_height"), "Height (px)",
                     value = 600, min = 200, max = 5000, step = 50)
      )
    ),
    numericInput(paste0(prefix, "_plot_dpi"), "Resolution (DPI)",
                 value = 100, min = 72, max = 600, step = 10)
  )
}

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
  hjust <- switch(align, left = 0, center = 0.5, right = 1, 0.5)
  if (nzchar(trimws(title)))    p <- p + ggplot2::labs(title = title)
  if (nzchar(trimws(subtitle))) p <- p + ggplot2::labs(subtitle = subtitle)
  if (nzchar(trimws(xlab)))     p <- p + ggplot2::labs(x = xlab)
  if (nzchar(trimws(ylab)))     p <- p + ggplot2::labs(y = ylab)
  if (nzchar(trimws(legend_lab))) p <- p + ggplot2::labs(fill = legend_lab, color = legend_lab, shape = legend_lab)
  p + ggplot2::theme(
    text          = ggplot2::element_text(size = text_size),
    plot.title    = ggplot2::element_text(hjust = hjust),
    plot.subtitle = ggplot2::element_text(hjust = hjust)
  )
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

#' Render colourInput widgets for a set of treatment groups
#'
#' Returns a \code{tagList} of \code{colourpicker::colourInput} widgets,
#' one per group, with auto-assigned ggplot2-style default colours.
#'
#' @param prefix    Input-ID prefix (e.g. \code{"tg_color"}).
#' @param groups    Character vector of group names.
#' @param header    Optional header string. Default \code{"Group Colors"}.
#' @return A \code{tagList}, or \code{NULL} when \code{groups} is empty.
#' @keywords internal
group_color_pickers_ui <- function(prefix, groups, header = "Group Colors") {
  if (is.null(groups) || length(groups) == 0L) return(NULL)
  defaults <- default_group_colors(length(groups))
  shiny::tagList(
    shiny::tags$h6(style = "font-size: 14px; font-weight: 700; text-decoration: underline; margin-top: 10px; margin-bottom: 4px;",
                   shiny::icon("palette"), " ", header),
    lapply(seq_along(groups), function(i) {
      colourpicker::colourInput(
        inputId = paste0(prefix, "_", make.names(groups[i])),
        label   = groups[i],
        value   = defaults[i]
      )
    })
  )
}

#' Render a drag-and-drop rank list for reordering treatment groups
#'
#' Returns a \code{sortable::rank_list} widget.  The current order is
#' read back from \code{input[[paste0(prefix, "_order")]]}.
#'
#' @param prefix    Input-ID prefix shared with colour pickers (e.g. \code{"tg"}).
#' @param groups    Character vector of group names in their default order.
#' @return A \code{sortable::rank_list} tag, or \code{NULL} when fewer than 2 groups.
#' @keywords internal
group_order_ui <- function(prefix, groups) {
  if (is.null(groups) || length(groups) < 2L) return(NULL)
  sortable::rank_list(
    text     = "Group order",
    labels   = as.list(groups),
    input_id = paste0(prefix, "_order"),
    class    = "default-sortable"
  )
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
  if (pos == "inplot_tl") {
    p + ggplot2::theme(
      legend.position    = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.background  = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), colour = NA)
    )
  } else if (pos == "inplot_tr") {
    p + ggplot2::theme(
      legend.position    = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background  = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), colour = NA)
    )
  } else if (pos == "inplot_bl") {
    p + ggplot2::theme(
      legend.position    = c(0.02, 0.02),
      legend.justification = c("left", "bottom"),
      legend.background  = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), colour = NA)
    )
  } else if (pos == "inplot_br") {
    p + ggplot2::theme(
      legend.position    = c(0.98, 0.02),
      legend.justification = c("right", "bottom"),
      legend.background  = ggplot2::element_rect(
        fill = scales::alpha("white", 0.8), colour = NA)
    )
  } else {
    p + ggplot2::theme(legend.position = pos)
  }
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
  shape_choices <- c(
    "\u25cf Circle"        = "16",
    "\u25a0 Square"        = "15",
    "\u25b2 Triangle \u2191" = "17",
    "\u25bc Triangle \u2193" = "25",
    "\u25c6 Diamond"       = "18"
  )
  default_shapes <- c(16L, 15L, 17L, 25L, 18L)

  # Each row: drag handle + group label + colour picker + shape selector
  # data-value holds the exact group name string — read directly by JS,
  # no innerText parsing, no custom Shiny input type.
  items <- lapply(seq_along(groups), function(i) {
    g        <- groups[i]
    color_id <- ns(paste0(prefix, "_color_", make.names(g)))
    shape_id <- ns(paste0(prefix, "_shape_", make.names(g)))
    def_shape <- as.character(default_shapes[((i - 1L) %% length(default_shapes)) + 1L])
    shiny::tags$li(
      `data-value` = g,
      style = paste0(
        "display:flex; align-items:center; gap:6px;",
        "padding:5px 8px; margin-bottom:4px;",
        "background:#f8f9fa; border:1px solid #dee2e6;",
        "border-radius:4px; cursor:default;"
      ),
      shiny::tags$span(
        class = "drag-handle",
        style = "cursor:grab; color:#adb5bd; font-size:14px; user-select:none; flex-shrink:0;",
        "\u2630"   # grip / hamburger icon
      ),
      shiny::tags$span(
        style = "flex:1; font-size:13px; min-width:55px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap;",
        g
      ),
      colourpicker::colourInput(
        inputId          = color_id,
        label            = NULL,
        value            = defaults[i],
        showColour       = "background",
        palette          = "square",
        allowTransparent = FALSE
      ),
      shiny::tags$div(
        style = "flex-shrink:0;",
        shiny::selectInput(
          inputId  = shape_id,
          label    = NULL,
          choices  = shape_choices,
          selected = def_shape,
          width    = "125px"
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

# ---- Dosing-schedule annotation (CODE_REVIEW.md D.15) ------------------------
#
# Replaces the removed `plot_treatments()`, which needed a separately uploaded
# schedule file and a second panel kept aligned by hand. Dosing timing is not
# recoverable from measurement data — `Dose` is constant within an animal, and
# dosing days are typically not on the measurement grid — so the user supplies it
# as text and it is drawn into the plot's lower margin.

#' Read a module's dosing-schedule inputs
#'
#' @param input,prefix Shiny input object and module prefix (e.g. `"tg"`).
#' @param groups Treatment groups present in the data.
#' @return A list with `schedule` (or `NULL`) and `messages`, a character vector
#'   of parser warnings. The messages are returned rather than swallowed: a
#'   mistyped arm name silently dropping its row is the §K.2 failure this review
#'   has found five times.
#' @keywords internal
resolve_dose_schedule <- function(input, prefix, groups) {
  if (!isTRUE(input[[paste0(prefix, "_show_dosing")]])) {
    return(list(schedule = NULL, messages = character(0)))
  }
  txt <- input[[paste0(prefix, "_dosing_days")]]
  msgs <- character(0)
  sched <- withCallingHandlers(
    tryCatch(mouseExperiment::parse_dose_schedule(txt, groups = groups),
             error = function(e) { msgs <<- c(msgs, conditionMessage(e)); NULL }),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(schedule = sched, messages = msgs)
}

#' Add a dosing-schedule annotation to a time-axis plot
#'
#' @param p A ggplot with study day on x.
#' @param input,prefix,groups As for [resolve_dose_schedule()].
#' @param colors Named colour vector so marks match their curves.
#' @return The annotated plot, or `p` unchanged when no schedule is in force.
#'   Annotation failure never costs the user their plot.
#' @keywords internal
apply_dose_schedule <- function(p, input, prefix, groups, colors = NULL) {
  res <- resolve_dose_schedule(input, prefix, groups)
  if (is.null(res$schedule)) return(p)
  style <- input[[paste0(prefix, "_dosing_style")]] %||% "auto"
  tryCatch(
    mouseExperiment::annotate_dose_schedule(p, res$schedule, style = style,
                                            colors = colors),
    error = function(e) {
      log_caught_error("plot", paste0(prefix, " dose-schedule annotation"), e)
      p
    })
}

#' Sidebar controls for the dosing-schedule annotation
#'
#' @param ns Module namespace function.
#' @param prefix Module prefix.
#' @return A `tagList` of the toggle, the days field, the style selector and a
#'   message slot.
#' @keywords internal
dose_schedule_ui <- function(ns, prefix) {
  shiny::tagList(
    shiny::checkboxInput(ns(paste0(prefix, "_show_dosing")),
                         "Show dosing schedule on plots", value = FALSE),
    shiny::conditionalPanel(
      ns = ns,
      condition = paste0("input.", prefix, "_show_dosing == true"),
      shiny::textAreaInput(
        ns(paste0(prefix, "_dosing_days")), "Dosing days",
        placeholder = "1, 5, 9, 13", rows = 2, resize = "vertical"),
      shiny::helpText(
        "Days on which doses were given. Enter a list — ",
        shiny::tags$code("1, 5, 9, 13"), " — to mark every treated arm, or one ",
        "line per arm when schedules differ: ",
        shiny::tags$code("DrugA: 1,5,9,13"), ". Marks are drawn below the plot ",
        "so they never cover the data."),
      shiny::selectInput(
        ns(paste0(prefix, "_dosing_style")), "Mark style",
        choices  = c("Automatic"          = "auto",
                     "Individual doses"   = "rug",
                     "Dosing window"      = "band"),
        selected = "auto"),
      shiny::helpText(
        "Automatic draws one mark per dose, switching to a shaded window when ",
        "dosing is daily or there are more than eight doses — twenty marks is ",
        "noise, a window reads at a glance."),
      shiny::uiOutput(ns(paste0(prefix, "_dosing_msg")))
    )
  )
}

#' Render parser feedback for the dosing-schedule field
#'
#' @param input,prefix,groups As above.
#' @return A UI element, or `NULL` when the schedule parsed cleanly.
#' @keywords internal
dose_schedule_msg_ui <- function(input, prefix, groups) {
  res <- resolve_dose_schedule(input, prefix, groups)
  if (!length(res$messages)) return(NULL)
  status_alert(paste(res$messages, collapse = " "), "warning")
}
