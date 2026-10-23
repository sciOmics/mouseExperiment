# =============================================================================
# "Download Results": the user chooses which tables, plots, diagnostics, data
# and report go into one ZIP, with a README.txt listing every file, the run's
# settings and the package versions (2026-10-22, CODE_REVIEW.md R20-D.X).
#
# It replaced "Download All", which zipped a fixed set of tables and PNG plots
# and left out the diagnostics, the analysed data and the report.
#
# A module describes what it can export as a list of dl_item()s built from its
# result. The checklist is built from that list alone; plots are drawn only when
# the ZIP is written, so changing a colour on the Plots tab does not rebuild the
# checklist and reset the user's ticks.
# =============================================================================

# The checklist's groups, in display order. Each is also a folder in the ZIP,
# except "code" and "report", whose files sit at the top level.
DL_GROUPS <- c(tables      = "Tables (CSV)",
               plots       = "Plots",
               diagnostics = "Diagnostics",
               data        = "Data (CSV)",
               code        = "R script",
               report      = "Report")
DL_TOP_LEVEL <- c("code", "report")

#' One file the "Download Results" tab can put in its ZIP
#'
#' @param id Unique within the module; the checkbox value.
#' @param group One of `names(DL_GROUPS)`.
#' @param label What the checkbox says.
#' @param file File name; tables and plots without the extension.
#' @param kind "table" (a data frame), "plot" (a ggplot or other printable
#'   plot), "file" (a function that writes the file to the path it is given)
#'   or "files" (a function that writes several files into the directory it is
#'   given and returns their paths relative to it, e.g. a script and its data).
#' @param value For "table" and "plot", a function returning the object (or
#'   NULL when it cannot be made); for "file", a function(path); for "files",
#'   a function(dir).
#' @param about One line describing the file, for README.txt.
#' @param row_names Write the data frame's row names (e.g. an ANOVA table).
#' @keywords internal
dl_item <- function(id, group, label, file, kind = c("table", "plot", "file", "files"),
                    value, about = label, row_names = FALSE) {
  kind <- match.arg(kind)
  stopifnot(group %in% names(DL_GROUPS), is.function(value))
  list(id = id, group = group, label = label, file = file, kind = kind,
       value = value, about = about, row_names = row_names)
}

#' The "Download Results" tab's contents
#' @param ns The module's namespace function.
#' @param prefix Input-ID prefix for this module (e.g. "tg").
#' @param upload Offer to include the uploaded file, so the R script starts
#'   from it (R20-AC); not on Power, which reads no data.
#' @keywords internal
download_results_ui <- function(ns, prefix, upload = TRUE) {
  id <- function(x) ns(paste0(prefix, "_dl_", x))
  tagList(
    tags$p("Choose what to download. It all goes into one ZIP file, with a README.txt ",
           "that lists every file, the settings of the run and the package versions."),
    uiOutput(id("choices")),
    if (upload) checkboxInput(id("include_upload"), tagList(
      "Include my original file, so the R script starts from it",
      help_tip(paste(
        "The ZIP then holds the file you uploaded, as it was, and the R script rebuilds",
        "the analysed data from it with the package's preparation functions: your",
        "exclusions, the study days, the volumes and the arm labels. It checks the result",
        "against the data the dashboard analysed before running the analysis. The file",
        "includes any animals you excluded and columns no analysis uses."))),
      value = FALSE),
    fluidRow(
      column(6,
        radioButtons(id("plot_format"),
          tagList("Plot format", help_tip(paste(
            "PNG for slides and documents. PDF or SVG scale without losing quality and",
            "can be edited in Illustrator or Inkscape. Plots use the size and resolution",
            "set under Plot Size / Export on the Plots tab."))),
          choices = c(PNG = "png", PDF = "pdf", SVG = "svg"), selected = "png",
          inline = TRUE)),
      column(6, style = "padding-top:28px;",
        actionLink(id("all"), "Select all"), " \u00b7 ",
        actionLink(id("none"), "Select none"))
    ),
    downloadButton(id("zip"), "Download selected (ZIP)",
                   class = "btn-lg btn-primary", icon = icon("download"))
  )
}

#' Server side of the "Download Results" tab
#'
#' @param input,output,session The module's.
#' @param prefix Input-ID prefix, as in [download_results_ui()].
#' @param module Short module name, for the ZIP's file name and analytics.
#' @param items A reactive returning the list of [dl_item()]s for the current
#'   result, or NULL before a successful run. It must not draw plots.
#' @param readme_header A function returning the README's opening lines.
#' @param plot_size A function returning list(width, height, dpi): pixels and
#'   pixels per inch.
#' @keywords internal
add_download_results <- function(input, output, session, prefix, module,
                                 items, readme_header, plot_size) {
  ns  <- session$ns
  id  <- function(x) paste0(prefix, "_dl_", x)
  gid <- function(g) id(paste0("group_", g))
  current <- function() items() %||% list()

  output[[id("choices")]] <- renderUI({
    it <- current()
    if (!length(it)) {
      return(div(class = "alert alert-info", icon("circle-info"),
                 " Run an analysis first; what it produced will be listed here."))
    }
    groups <- intersect(names(DL_GROUPS), vapply(it, `[[`, "", "group"))
    fluidRow(lapply(groups, function(g) {
      gi <- Filter(function(x) x$group == g, it)
      ids <- vapply(gi, `[[`, "", "id")
      column(4, checkboxGroupInput(ns(gid(g)), DL_GROUPS[[g]],
        choices = stats::setNames(ids, vapply(gi, `[[`, "", "label")),
        selected = ids))
    }))
  })

  select <- function(all) {
    it <- current()
    for (g in names(DL_GROUPS)) {
      ids <- vapply(Filter(function(x) x$group == g, it), `[[`, "", "id")
      if (length(ids))
        updateCheckboxGroupInput(session, gid(g), selected = if (all) ids else character(0))
    }
  }
  observeEvent(input[[id("all")]],  select(TRUE))
  observeEvent(input[[id("none")]], select(FALSE))

  output[[id("zip")]] <- downloadHandler(
    filename = function() {
      paste0("mouseExperiment_", module, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      chosen <- unlist(lapply(names(DL_GROUPS), function(g) input[[gid(g)]]))
      it <- Filter(function(x) x$id %in% chosen, current())
      write_results_zip(file, it, fmt = input[[id("plot_format")]] %||% "png",
                        size = plot_size(), header = readme_header(),
                        opts = list(include_upload = isTRUE(input[[id("include_upload")]])))
      shinyjs::runjs(sprintf("window.trackDownload(%s, 'selected_results', 'zip');",
                             js_string(module)))
      log_info(sprintf("download_results | module=%s | files=%d", module, length(it)))
    },
    contentType = "application/zip"
  )
}

#' Write the chosen items and a README into a ZIP
#'
#' An item whose object cannot be made (NULL or an error) is left out and
#' listed in the README as not available, rather than failing the download.
#' `opts` goes to "files" items: `include_upload`, the "Include my original
#' file" tick (R20-AC).
#' @return The ZIP's path, invisibly.
#' @keywords internal
write_results_zip <- function(zipfile, items, fmt = c("png", "pdf", "svg"),
                              size = list(width = 800, height = 600, dpi = 100),
                              header = character(0), opts = list()) {
  fmt <- match.arg(fmt)
  outdir <- tempfile("results_")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)

  written <- character(0)
  listing <- character(0)
  missing <- character(0)
  for (x in items) {
    folder <- if (x$group %in% DL_TOP_LEVEL) "" else x$group
    dir <- if (nzchar(folder)) file.path(outdir, folder) else outdir
    dir.create(dir, showWarnings = FALSE)
    f <- tryCatch(dl_write(x, dir, fmt, size, opts), error = function(e) NULL)
    if (!length(f)) {
      missing <- c(missing, x$label)
      next
    }
    # A "files" item writes several; the first is the item, the rest its data.
    rel <- if (nzchar(folder)) file.path(folder, f) else f
    new <- !rel %in% written
    written <- c(written, rel[new])
    listing <- c(listing, sprintf("  %-40s %s", rel[new],
                                  c(x$about, rep(paste("read by", rel[1L]), length(rel) - 1L))[new]))
  }

  readme <- c(header, "", "Files", listing)
  if (length(missing))
    readme <- c(readme, "", "Chosen but not available for this run",
                paste0("  ", missing))
  writeLines(readme, file.path(outdir, "README.txt"))

  old <- setwd(outdir)
  on.exit(setwd(old), add = TRUE, after = FALSE)
  utils::zip(zipfile, c("README.txt", written), flags = "-q9X")
  invisible(zipfile)
}

# Write one item into `dir`; returns its file name (several for "files"), or
# NULL when the object is not available.
dl_write <- function(x, dir, fmt, size, opts = list()) {
  if (identical(x$kind, "files")) {
    rel <- x$value(dir, opts)
    return(rel[file.exists(file.path(dir, rel))])
  }
  if (identical(x$kind, "file")) {
    x$value(file.path(dir, x$file))
    return(if (file.exists(file.path(dir, x$file))) x$file else NULL)
  }
  v <- x$value()
  if (is.null(v)) return(NULL)
  if (identical(x$kind, "table")) {
    v <- as.data.frame(v)
    if (!nrow(v)) return(NULL)
    f <- paste0(x$file, ".csv")
    utils::write.csv(v, file.path(dir, f), row.names = x$row_names)
    return(f)
  }
  f <- paste0(x$file, ".", fmt)
  path <- file.path(dir, f)
  w <- size$width %||% 800
  h <- size$height %||% 600
  dpi <- size$dpi %||% 100
  # grDevices, not ggsave: ggplot2 saves SVG only through svglite, which the
  # image does not have; R's own svg() device uses Cairo.
  switch(fmt,
    png = grDevices::png(path, width = w, height = h, units = "px", res = dpi),
    pdf = grDevices::pdf(path, width = w / dpi, height = h / dpi),
    svg = grDevices::svg(path, width = w / dpi, height = h / dpi))
  # A plot object is printed; a function draws (base graphics, e.g. cox.zph).
  tryCatch(if (is.function(v)) v() else print(v), finally = grDevices::dev.off())
  f
}

#' The README's opening lines: module, time, versions and the run's settings
#' @param module_title E.g. "Tumor Growth".
#' @param res The module's stored result.
#' @param settings_fn A function(res) returning the run's settings as a named
#'   list (as for the HTML report), or NULL.
#' @keywords internal
dl_readme_header <- function(module_title, res, settings_fn = NULL) {
  # As written in DESCRIPTION: packageVersion() prints 26.10.008 as 26.10.8.
  ver <- function(p) tryCatch(utils::packageDescription(p)$Version %||% "?",
                              error = function(e) "?")
  s <- if (is.function(settings_fn)) tryCatch(settings_fn(res), error = function(e) NULL)
  s <- Filter(function(v) length(v) && !all(is.na(v)) && any(nzchar(as.character(v))), s %||% list())
  c(paste("mouseExperiment dashboard:", module_title, "results"),
    paste("Downloaded:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
    paste0("Versions: mouseExperiment ", ver("mouseExperiment"), ", R ", getRversion()),
    if (length(s)) c("", "Settings of the run",
                     sprintf("  %s: %s", names(s),
                             vapply(s, function(v) paste(as.character(v), collapse = ", "), ""))))
}

