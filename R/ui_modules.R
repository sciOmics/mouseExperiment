# Shared UI helpers for the dashboard: info boxes, help icons, plot notes and
# collapsible sections. The Data tab's UI is in ui_data.R, Tumor Growth's in
# ui_tumor_growth.R, Help and About in ui_help.R (CODE_REVIEW.md R22-K); the
# other tabs' are with their modules.

# Shared helper: creates a styled info box for tab descriptions
tab_info <- function(...) {
  div(
    class = "alert alert-info",
    style = "font-size:0.88em; margin:10px 0 8px 0; padding:8px 14px;",
    ...
  )
}

# Shared helper: inline help icon. It shows `text` in a tooltip on hover, on
# keyboard focus, and on a click or tap (www/js/help-tips.js). A bare title
# attribute showed nothing on a click or a tap. The icon sits inside input
# labels, and that script also keeps its clicks from reaching the label, which
# passed them on to the input: the Control group menu opened instead.
help_tip <- function(text) {
  tags$span(
    class = "help-tip text-muted",
    style = "cursor:help; margin-left:4px;",
    tabindex = "0", role = "button", `aria-label` = paste("Help:", text),
    title = text,
    icon("circle-question", class = "fa-sm")
  )
}

# Shared helper: a "How to read this plot" note under a plot's title, collapsed
# until clicked. It used to be a small grey line below each plot, often a
# screen further down, and was easy to miss (2026-10-22).
plot_help <- function(..., summary = "How to read this plot") {
  tags$details(
    class = "plot-help",
    tags$summary(icon("circle-info"), " ", summary),
    tags$div(class = "plot-help-body", ...)
  )
}

# Helper: collapsible section for the plot control panel.
# Starts collapsed; click the header to expand.
ctrl_section <- function(id, icon_nm, title, ...) {
  collapse_id <- paste0("ctrl-", id)
  tagList(
    tags$a(
      href = paste0("#", collapse_id),
      `data-toggle` = "collapse",
      class = "collapsed",
      style = paste0(
        "display:block; font-size:1.2rem; font-weight:700; ",
        "color:#333; text-decoration:none; ",
        "padding:5px 2px; margin-top:4px; user-select:none; ",
        "cursor:pointer;"
      ),
      icon(icon_nm, style = "margin-right:5px;"),
      title,
      tags$span(
        class = "pull-right text-muted",
        style = "font-size:0.75rem; padding-top:3px;",
        HTML("&#9660;")
      )
    ),
    tags$div(id = collapse_id, class = "collapse", ...),
    hr(style = "margin:5px 0;")
  )
}

