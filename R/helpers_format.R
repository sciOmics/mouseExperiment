# Scale, formatting and result-provenance helpers
#
# Split out of the former `helpers.R` under CODE_REVIEW.md D.1.
#
# These decide how a number is presented, which is not a cosmetic question when
# the number lives on a transformed scale. `back_transform()` and
# `scale_suffix()` are the pair that turn a model estimate back into a volume; a
# disagreement between them mislabels every figure on screen. `comparisons_title()`
# and `resolve_interval_cols()` report multiplicity and interval provenance, and
# the latter exists because the previous column-name sniffing silently returned
# NA and blanked the panel (R4.D5).

#' Format numeric values for display
#'
#' @param x Numeric value
#' @param digits Number of decimal places
#'
#' @return Formatted character string
#' @keywords internal
format_number <- function(x, digits = 3) {
  if (is.na(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

#' Format a p-value for display
#'
#' Very small p-values are shown in scientific notation rather than being
#' rounded to "0". Values that underflow to exactly 0 are reported as
#' "< 2.2e-16" (the smallest representable positive double). Uses
#' \code{formatC} with \code{format="g"} which automatically switches to
#' scientific notation for small values.
#'
#' @param x Numeric p-value(s)
#' @param digits Significant digits (default 3)
#' @return Character string(s)
#' @keywords internal
format_pval <- function(x, digits = 3) {
  limit <- .Machine$double.eps  # ~2.2e-16
  ifelse(
    is.na(x), NA_character_,
    ifelse(
      x == 0 | (!is.na(x) & x < limit),
      paste0("< ", formatC(limit, digits = 2, format = "e")),
      formatC(x, digits = digits, format = "g")
    )
  )
}

#' Format a data frame for tabular display
#'
#' Rounds ordinary numeric columns to \code{digits} decimal places.
#' Columns whose names match common p-value patterns are formatted with
#' \code{\link{format_pval}} so very small values appear in scientific
#' notation rather than being displayed as "0".
#'
#' @param df Data frame
#' @param digits Decimal places for regular numeric columns (default 4)
#' @return Data frame with formatted columns
#' @keywords internal
format_df_for_display <- function(df, digits = 4) {
  # Flatten any array/table/matrix columns to plain atomic vectors first.
  # Columns returned by survival_statistics (Events, Total, Event_Rate) can be
  # table- or array-class objects; round() and datatable() both have trouble with
  # these non-standard types.
  for (col in names(df)) {
    v <- df[[col]]
    if (!is.null(dim(v)) || inherits(v, "table")) {
      df[[col]] <- as.vector(v)
    }
  }
  pval_pattern <- "p[._]?value|p[._]?adj|padj|^p$|Pr\\("
  for (col in names(df)) {
    if (!is.numeric(df[[col]])) next
    tryCatch({
      if (grepl(pval_pattern, col, ignore.case = TRUE)) {
        df[[col]] <- format_pval(df[[col]])
      } else {
        df[[col]] <- round(df[[col]], digits)
      }
    }, error = function(e) {
      # Log and leave column as-is if formatting fails
      log_caught_error("format", paste0("format_df_for_display col '", col, "'"), e)
    })
  }
  df
}

#' The zero-floor line of the tumour-growth summary
#'
#' Backend 0.30.0 (R20.32): on the log scale a zero volume takes a fixed
#' floor, 2 mm³ or 0.002 cm³, reported in `result$zero_floor`.
#'
#' @param res A `tumor_growth_statistics()` result.
#' @return One summary line, or `character(0)` when no floor applies.
#' @keywords internal
zero_floor_line <- function(res) {
  zf <- res$zero_floor
  if (is.null(zf) || !identical(res$transform_used, "log")) return(character(0))
  unit <- switch(zf$units %||% "mm3", cm3 = "cm\u00b3", "mm\u00b3")
  sprintf("  Zero volumes: set to %s %s before the log (%d value%s)",
          format(zf$value), unit, as.integer(zf$n_floored %||% 0L),
          if (identical(as.integer(zf$n_floored), 1L)) "" else "s")
}

#' The tumour-growth data summary, labelled for display -- R20.D45
#'
#' @param ds The backend's `data_summary` (volume scale since 0.31.0).
#' @param unit The declared unit's label, from [volume_unit_label()].
#' @return A data frame with display headers.
#' @keywords internal
data_summary_display <- function(ds, unit) {
  ds <- as.data.frame(ds)
  labels <- c(Treatment = "Treatment", Day = "Day", n = "N",
              mean_volume = paste0("Mean (", unit, ")"),
              sd_volume = "SD", sem_volume = "SEM",
              geometric_mean_volume = paste0("Geometric mean (", unit, ")"))
  keep <- intersect(names(labels), names(ds))
  ds <- ds[, keep, drop = FALSE]
  names(ds) <- unname(labels[keep])
  ds
}

#' What the forest plot's intervals are -- R20.D43
#'
#' @param pc_df The comparison table the plot is drawn from.
#' @param lower_col The column the lower limits came from.
#' @return One caption line.
#' @keywords internal
forest_ci_note <- function(pc_df, lower_col) {
  if (identical(lower_col, "Lower_CrI")) {           # Bayesian runs (R20.D59)
    return("Intervals: 95% posterior credible intervals")
  }
  adj <- if ("CI_Adjustment" %in% names(pc_df)) as.character(pc_df$CI_Adjustment[1]) else NA_character_
  if (is.na(lower_col) || !lower_col %in% c("CI_Adjusted_Lower", "lower.CL") || is.na(adj)) {
    return("Intervals: 95%, not adjusted for multiple comparisons")
  }
  if (startsWith(adj, "none:")) {          # BH: no simultaneous interval exists
    return("Intervals: 95% per comparison; the FDR correction has no matching interval")
  }
  if (startsWith(adj, "none")) return("Intervals: 95% per comparison (no correction applied)")
  paste0("Intervals adjusted for multiple comparisons (", adj, "), like the p-values")
}

#' Back-transform values from analysis scale to original scale
#'
#' Inverse of the transform applied in `tumor_growth_statistics()` and
#' `bayesian_tumor_growth()`, which both use `log(Volume)` (not
#' `log1p`). For `"sqrt"`, returns `x^2`. For `"none"`, returns `x` unchanged.
#'
#' @param x Numeric vector on the transformed scale
#' @param transform One of "log", "sqrt", or "none"
#' @return Numeric vector on the original scale
#' @keywords internal
back_transform <- function(x, transform) {
  x <- suppressWarnings(as.numeric(x))
  switch(transform,
    "log"  = exp(x),
    "sqrt" = x ^ 2,
    "none" = x,
    x
  )
}

#' Axis-label suffix for an analysis-scale plot
#'
#' Returns " (log scale)", " (sqrt scale)", or "" so plot y-axis labels can
#' indicate that values are on the transformed scale without overriding any
#' user-supplied label.
#'
#' @param transform One of "log", "sqrt", or "none"
#' @return Character scalar suitable for appending to an axis label
#' @keywords internal
scale_suffix <- function(transform) {
  switch(transform,
    "log"  = " (log scale)",
    "sqrt" = " (sqrt scale)",
    ""
  )
}


#' A correction's display name
#' @keywords internal
pretty_adjust <- function(adj) {
  switch(adj %||% "",
    bonferroni = "Bonferroni", holm = "Holm",
    fdr = "Benjamini-Hochberg FDR",
    dunnett = "Dunnett (exact, multivariate-t)",
    tukey = "Tukey HSD", none = "unadjusted", adj %||% "")
}

#' The p-value line of the Comparisons help, from the correction used -- R20.D44
#'
#' It said "Dunnett ... Tukey HSD" whatever the run had used, including no
#' correction at all.
#'
#' @param res A tumour-growth result, or NULL before a run.
#' @return One sentence.
#' @keywords internal
p_value_help <- function(res) {
  adj <- res$p_adjust_method_used
  if (is.null(adj)) return("p-value, adjusted by the correction chosen in the sidebar.")
  if (identical(adj, "none")) {
    return("p-value, not adjusted for multiple comparisons (the sidebar's correction is None).")
  }
  paste0("p-value adjusted for multiple comparisons (", pretty_adjust(adj),
         "). Values below 0.05 indicate a significant difference.")
}

#' Build a comparisons-table title that states what was actually applied
#'
#' CODE_REVIEW.md R3.D1 / R3.D2 — the headers used to assert "Dunnett-Adjusted"
#' / "Tukey HSD Adjusted" regardless of what had been computed, and for the GAM
#' path the table silently held per-day contrasts under a header that did not say
#' so. Read the adjustment from the result instead of asserting it.
#'
#' @param res A tumour-growth result list.
#' @param family_label Human-readable family name.
#' @return A character title.
#' @keywords internal
comparisons_title <- function(res, family_label) {
  adj   <- res$p_adjust_method_used %||% NULL
  scope <- if (!is.null(res$pairwise_comparisons) &&
               is.data.frame(res$pairwise_comparisons) &&
               "Adjust_Scope" %in% names(res$pairwise_comparisons)) {
    unique(res$pairwise_comparisons$Adjust_Scope)[1]
  } else NULL

  adj_txt <- if (is.null(adj)) "" else paste0(" (", pretty_adjust(adj), ")")
  scope_txt <- if (is.null(scope)) "" else paste0(" \u2014 ", scope)
  paste0("Comparisons: ", family_label, adj_txt, scope_txt)
}

#' Resolve a result's interval columns, preferring the backend's declaration
#'
#' CODE_REVIEW.md R4.D5. Column sniffing —
#' `intersect(c("lower.CL", "lower_CL", "Lower", "lower"), colnames(df))[1]` —
#' is the §K.2 pattern: a miss yields `NA`, the calling code takes a silent
#' fallback, and nothing fails. It had already gone wrong. The backend's
#' `treatment_effects` carries `Lower_CL`, which matches none of the four names
#' above, so the Treatment Effects plot had been drawing SE-derived error bars
#' instead of the model's confidence interval.
#'
#' Backend v0.12.0 attaches `meta$interval_columns` declaring the names it used.
#' Prefer that; fall back to sniffing only for results from an older backend, and
#' widen the sniff list to include the names that were being missed.
#'
#' @param df The data frame to read.
#' @param result Optional full result object carrying `$meta`.
#' @return A list with `lower` / `upper` column names, either `NA` when absent.
#' @keywords internal
resolve_interval_cols <- function(df, result = NULL) {
  meta <- result$meta
  if (!is.null(meta) && !identical(meta$interval_type, "none")) {
    lo <- meta$interval_columns[["lower"]]
    hi <- meta$interval_columns[["upper"]]
    if (!is.na(lo) && all(c(lo, hi) %in% names(df))) {
      return(list(lower = lo, upper = hi))
    }
    # Some tables in a result use different names from the headline one (e.g.
    # survival's CI_Lower/CI_Upper); honour an override when present.
    ov <- meta$interval_columns_override
    if (!is.null(ov) && all(ov %in% names(df))) {
      return(list(lower = ov[["lower"]], upper = ov[["upper"]]))
    }
  }
  # Legacy fallback. Note "Lower_CL"/"Upper_CL" and "Lower_CrI"/"Upper_CrI" —
  # both were absent from the original list.
  lo <- intersect(c("Lower_CL", "Lower_CrI", "lower.CL", "lower_CL",
                    "Lower", "lower"), colnames(df))[1]
  hi <- intersect(c("Upper_CL", "Upper_CrI", "upper.CL", "upper_CL",
                    "Upper", "upper"), colnames(df))[1]
  list(lower = lo, upper = hi)
}


#' A 95% interval as text, e.g. `"[12.3, 18.9]"`; a dash when there is none
#' @keywords internal
fmt_interval <- function(lo, hi, digits = 1) {
  ifelse(is.finite(lo) & is.finite(hi),
         sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), lo, hi),
         "\u2014")
}

#' The therapeutic window table as displayed
#'
#' R20-K (package v0.27.0): one row per arm with its TGI and mean worst
#' weight loss, each with a 95% interval, the animals over the threshold and
#' the tolerability flag. The CSV download keeps the backend's columns.
#' @param w The backend's `window_table`.
#' @param reference The reference arm, whose TGI is 0 by definition.
#' @keywords internal
window_display_table <- function(w, reference = NULL) {
  is_ref <- !is.null(reference) & w$Treatment %in% reference
  data.frame(
    Treatment                  = w$Treatment,
    Animals                    = w$N_Animals,
    `TGI (%)`                  = ifelse(is_ref, NA_real_, round(w$TGI, 1)),
    `TGI 95% interval`         = ifelse(is_ref, "reference",
                                        fmt_interval(w$TGI_Lower, w$TGI_Upper)),
    `Worst weight loss (%)`    = round(w$Worst_Loss, 1),
    `Weight loss 95% interval` = fmt_interval(w$Worst_Loss_Lower, w$Worst_Loss_Upper),
    `Over threshold`           = sprintf("%d of %d", as.integer(w$N_Over_Threshold),
                                         as.integer(w$N_Animals)),
    Tolerability               = ifelse(is.na(w$Tolerability), "No interval",
                                        w$Tolerability),
    check.names = FALSE, stringsAsFactors = FALSE)
}
