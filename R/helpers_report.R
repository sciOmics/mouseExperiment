# Self-contained HTML reports (CODE_REVIEW.md C.1/C.3, R22 proposal F4).

# ---- Reports (CODE_REVIEW.md Round 2 F.1; R22 proposal F4) -------------------
#
# Until 0.46.0 only Tumor Growth had a report, of five tables. Every tab's
# report is now built from what its Download Results offers: the run's
# settings, its tables, its plots (embedded) and its diagnostics. A study
# report puts every tab with a result in one file.

#' HTML-escape a value (vectors are joined with ", ")
#' @keywords internal
report_esc <- function(x) {
  if (is.null(x)) return("")
  xc <- as.character(x)
  xc <- gsub("&", "&amp;", xc, fixed = TRUE)
  xc <- gsub("<", "&lt;",  xc, fixed = TRUE)
  xc <- gsub(">", "&gt;",  xc, fixed = TRUE)
  paste(xc, collapse = ", ")
}

#' A data frame as an HTML table, formatted as on screen
#' @param df A data frame.
#' @param caption Optional caption.
#' @keywords internal
report_table <- function(df, caption = NULL) {
  df <- tryCatch(as.data.frame(df), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) return("<p><em>(empty)</em></p>")
  df <- format_df_for_display(df)
  # R20.D15: one header cell per column.
  header <- paste0("<tr>", paste0("<th>", vapply(colnames(df), report_esc, character(1L)), "</th>",
                                  collapse = ""), "</tr>")
  rows <- vapply(seq_len(nrow(df)), function(i) {
    cells <- vapply(seq_len(ncol(df)), function(j) report_esc(df[i, j]), character(1L))
    paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  }, character(1L))
  paste0(if (!is.null(caption)) paste0("<p class=\"caption\">", report_esc(caption), "</p>"),
         "<table><thead>", header, "</thead><tbody>", paste(rows, collapse = ""),
         "</tbody></table>")
}

#' One Download Results item as a piece of the report: a table, or a plot
#' drawn to PNG and embedded; NULL when it is not available
#' @param x A [dl_item()].
#' @param size list(width, height, dpi) for plots.
#' @keywords internal
report_item <- function(x, size) {
  if (x$kind == "table") {
    v <- tryCatch(x$value(), error = function(e) NULL)
    if (is.null(v)) return(NULL)
    return(paste0("<h3>", report_esc(sub(" \\(CSV\\)$", "", x$label)), "</h3>", report_table(v, x$about)))
  }
  if (x$kind != "plot" || !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  dir <- tempfile("report_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  f <- tryCatch(dl_write(x, dir, "png", size), error = function(e) NULL)
  if (!length(f) || !file.exists(file.path(dir, f))) return(NULL)
  png <- readBin(file.path(dir, f), "raw", n = file.size(file.path(dir, f)))
  paste0("<h3>", report_esc(x$label), "</h3><p class=\"caption\">", report_esc(x$about), "</p>",
         "<img alt=\"", report_esc(x$label), "\" src=\"data:image/png;base64,",
         jsonlite::base64_enc(png), "\">")
}

#' One tab's part of a report: its settings, tables, plots and diagnostics
#' @param title The tab, e.g. "Tumor Growth".
#' @param items The tab's Download Results items.
#' @param header The README's opening lines (see [dl_readme_header()]).
#' @param size Plot size: list(width, height, dpi).
#' @keywords internal
report_section <- function(title, items, header = character(0), size = list()) {
  part <- function(group, heading) {
    got <- Filter(Negate(is.null), lapply(Filter(function(x) x$group == group, items),
                                          report_item, size = size))
    if (length(got)) paste0("<h2>", heading, "</h2>", paste(got, collapse = ""))
  }
  paste0("<section><h1>", report_esc(title), "</h1>",
         if (length(header)) paste0("<pre class=\"meta\">", report_esc(paste(header, collapse = "\n")), "</pre>"),
         part("tables", "Tables"), part("plots", "Plots"), part("diagnostics", "Diagnostics"),
         "</section>")
}

#' Write a report: sections in one self-contained HTML file
#' @param file Output path.
#' @param title The document's title.
#' @param sections HTML from [report_section()].
#' @keywords internal
report_document <- function(file, title, sections) {
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>", report_esc(title), "</title>",
    "<style>",
    "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif; ",
    "max-width: 1100px; margin: 20px auto; padding: 0 16px; color: #222; }",
    "h1 { color: #2166AC; border-bottom: 2px solid #2166AC; padding-bottom: 6px; margin-top: 36px; }",
    "h2 { color: #444; margin-top: 28px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }",
    "h3 { margin: 20px 0 4px 0; font-size: 1.05em; }",
    "table { border-collapse: collapse; margin: 8px 0; font-size: 0.85em; width: 100%; }",
    "th, td { border: 1px solid #ccc; padding: 5px 8px; text-align: left; }",
    "thead { background: #f5f5f5; } tbody tr:nth-child(even) { background: #fafafa; }",
    "img { max-width: 100%; border: 1px solid #eee; }",
    ".meta { background: #f8f9fa; border-left: 4px solid #2166AC; padding: 8px 14px; ",
    "font-size: 0.85em; white-space: pre-wrap; }",
    ".caption { font-style: italic; color: #666; margin: 2px 0 4px 0; }",
    "footer { color: #999; margin-top: 40px; font-size: 0.8em; border-top: 1px solid #ddd; padding-top: 10px; }",
    "@media print { section { page-break-after: always; } }",
    "</style></head><body>",
    paste(sections, collapse = ""),
    "<footer>Generated by the mouseExperiment dashboard on ",
    format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    ". Print to PDF from your browser to keep it with the study.</footer></body></html>")
  writeLines(html, file)
  invisible(file)
}

#' The "HTML report" item a tab's Download Results adds for itself (F4)
#' @param title The tab.
#' @param items A function returning the tab's other items.
#' @param header,size Functions returning the README header and plot size.
#' @keywords internal
report_dl_item <- function(title, items, header, size) {
  dl_item("report", "report", "HTML report", "report.html", "file",
          function(path) report_document(path, paste(title, "report"),
                                         report_section(title, items(), header(), size())),
          "This run's settings, tables, plots and diagnostics in one page; open it in a browser")
}

#' The tabs registered for the study report, with a result now
#' @param registry `session$userData$dl_registry`: per tab, list(title,
#'   items, header, size), each a function.
#' @return HTML sections, in the registry's order.
#' @keywords internal
study_report_sections <- function(registry) {
  Filter(Negate(is.null), lapply(registry, function(r) {
    it <- tryCatch(r$items(), error = function(e) NULL)
    it <- Filter(function(x) x$group %in% c("tables", "plots", "diagnostics"), it %||% list())
    if (!length(it)) return(NULL)
    report_section(r$title(), it, tryCatch(r$header(), error = function(e) character(0)),
                   tryCatch(r$size(), error = function(e) list()))
  }))
}

#' The settings that produced a tumour-growth result -- R20.D46
#'
#' The report read the sidebar when it was downloaded, so a dropdown changed
#' after the run put "model_type: auc, transform: sqrt, reference_group:
#' Drug_B" above LME4 / log / Vehicle tables. The settings now come from what
#' the backend reports about the run, then from the snapshot taken when the
#' run started (`result$run_settings`).
#'
#' @param res A tumour-growth result as stored in `rv$tg_results`.
#' @return A named list for the README header and the report.
#' @keywords internal
tg_report_settings <- function(res) {
  rs <- res$run_settings %||% list()
  re <- res$random_effects
  list(
    model_type      = res$model_type_used %||% rs$model_type %||% "",
    transform       = res$transform_used %||% rs$transform %||% "",
    random_effects  = if (!is.null(re$used)) {
      if (identical(re$used, re$requested)) re$used
      else paste0(re$used, " (requested: ", re$requested, ")")
    } else rs$random_effects %||% "",
    reference_group = res$summary$data_description$reference_group %||%
                      rs$reference_group %||% "",
    cage_effect     = res$summary$model_specification$cage_placement %||%
                      rs$cage_effect %||% "",
    comparisons     = if (!is.null(res$comparison_family)) {
      paste0(res$comparison_family, ", ", pretty_adjust(res$p_adjust_method_used))
    } else "",
    zero_floor      = trimws(paste(zero_floor_line(res), collapse = "")),
    n_rows          = rs$n_rows %||% NA
  )
}

