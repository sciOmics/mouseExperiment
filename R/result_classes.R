# Result classes (CODE_REVIEW.md R22-I, R22 proposal F8). The analysis
# functions return lists, which they still are; each now also carries a class,
# so print() shows what was found instead of every field, and summary()
# returns the headline tables.

# What each result's print() and summary() show: a title, and the tables that
# carry its answer, by field. `about` gives a one-line description of the run.
.me_result_spec <- list(
  me_tumor_growth = list(
    title = function(x) paste0("Tumour growth (", switch(x$model_type_used %||% "lme4",
                                 lme4 = "LME4", auc = "AUC", x$model_type_used), ")"),
    tables = c(anova = "ANOVA", pairwise_comparisons = "Comparisons",
               treatment_effects = "Treatment effects"),
    about = function(x) c(if (!is.null(x$transform_used)) paste("Transform:", x$transform_used),
                          if (!is.null(x$p_adjust_method_used))
                            paste("Comparisons:", x$comparison_family, "with", x$p_adjust_method_used))),
  me_bayes_tumor_growth = list(
    title = function(x) "Tumour growth (Bayesian LMM)",
    tables = c(treatment_effects = "Treatment effects", pairwise_comparisons = "Comparisons",
               mcmc_diagnostics = "MCMC diagnostics")),
  me_survival = list(
    title = function(x) "Survival",
    tables = c(results = "Hazard ratios against the reference"),
    about = function(x) c(paste("Reference:", x$reference_group), paste("Method:", x$method_used))),
  me_bayes_survival = list(
    title = function(x) paste0("Survival (Bayesian AFT, ", x$family_used %||% "", ")"),
    tables = c(treatment_effects = "Time ratios against the reference",
               mcmc_diagnostics = "MCMC diagnostics")),
  me_body_weight = list(
    title = function(x) "Body weight (LME4)",
    tables = c(emmeans_table = "Rate of weight change by arm",
               pairwise_comparisons = "Rates against the reference")),
  me_weight_loss = list(
    title = function(x) paste0("Time to ", round(100 * x$threshold), " % weight loss"),
    tables = c(censoring_summary = "How each animal left the study"),
    about = function(x) c(paste("Reference:", x$reference_group),
                          if (!is.null(x$cox_skipped)) x$cox_skipped
                          else if (!is.na(x$cox_method %||% NA)) paste("Cox model:", x$cox_method))),
  me_therapeutic_window = list(
    title = function(x) "Therapeutic window",
    tables = c(window_table = "TGI, %T/C and worst weight loss by arm"),
    about = function(x) c(paste("Endpoint day:", x$endpoint_day), x$tolerability_rule)),
  me_synergy = list(
    title = function(x) "Drug synergy (Bliss independence)",
    tables = c(synergy_metrics = "Verdict", synergy_ci = "Intervals"),
    about = function(x) paste0("Day ", x$evaluation_time_point, "; ", x$drug_a_name, " + ",
                               x$drug_b_name, " as ", x$combo_name, " against ", x$control_name)),
  me_synergy_over_time = list(
    title = function(x) "Drug synergy over time (Bliss independence)",
    tables = c(synergy_summary = "Bliss excess by day")),
  me_dose_response = list(
    title = function(x) "Dose response",
    tables = c(tgi_table = "TGI and %T/C by dose", summary_table = "Volume by dose"),
    about = function(x) paste("Endpoint day:", x$endpoint_day)),
  me_power = list(
    title = function(x) "A priori power (analytic, volume endpoint)",
    tables = c(scenario_table = "Scenarios"),
    about = function(x) x$method_note),
  me_survival_power = list(
    title = function(x) "A priori power (survival endpoint, log-rank)",
    tables = c(scenario_table = "Scenarios"),
    about = function(x) x$method_note),
  me_power_simulation = list(
    title = function(x) "A priori power (LMM simulation)",
    tables = c(power_table = "Power by N per group")),
  me_permutation = list(
    title = function(x) "Randomisation test of the treatment x time interaction",
    tables = character(0),
    about = function(x) c(sprintf("p = %s (%s, %d permutations of %d %s units)",
                                  format(signif(x$p_value, 3)), x$method %||% "",
                                  as.integer(x$n_perm %||% NA), as.integer(x$n_units %||% NA),
                                  x$unit %||% ""),
                          sprintf("Smallest attainable p: %s", format(signif(x$min_attainable_p, 3))))),
  me_responses = list(
    title = function(x) "mRECIST responses",
    tables = c(summary = "Animals in each category, and the response rate, by arm"),
    about = function(x) paste("Responses count from day", x$min_day, "after baseline")),
  me_growth_delay = list(
    title = function(x) paste("Tumour growth delay to", format(x$volume)),
    tables = c(delay_table = "Median time to the volume, and T - C"),
    about = function(x) paste("Reference:", x$reference_group)),
  me_data_check = list(
    title = function(x) "Data check",
    tables = c(summary = "What was found"),
    about = function(x) if (!nrow(x$summary)) "Nothing found." else
      paste0(sum(x$summary$Severity == "Problem"), " problem(s), ",
             sum(x$summary$Severity == "Check"), " thing(s) to check, ",
             sum(x$summary$Severity == "Note"), " note(s); the rows are in $flagged.")))

#' Give an analysis result its class
#' @param x The result list.
#' @param class One of the names of `.me_result_spec`.
#' @return `x`, with class `c(class, "me_result")`.
#' @noRd
me_result <- function(x, class) {
  stopifnot(is.list(x), class %in% names(.me_result_spec))
  structure(x, class = c(class, "me_result"))
}

# Row names worth printing: an ANOVA's terms, not 1, 2, 3.
me_named_rows <- function(df) {
  rn <- rownames(df)
  !is.null(rn) && !identical(rn, as.character(seq_len(nrow(df))))
}

# The spec for a result's class.
me_result_spec <- function(x) .me_result_spec[[intersect(class(x), names(.me_result_spec))[1L]]]

#' Print and summarise analysis results
#'
#' The analysis functions return lists with a class (since 0.47.0). `print()`
#' shows what the analysis found: a title, a line or two about the run, and the
#' first rows of its headline tables. Every field is still there: see
#' `names(result)`, or index the result as a list. `summary()` returns the
#' headline tables, whole.
#'
#' @param x,object A result from one of the analysis functions, such as
#'   [tumor_growth_statistics()] or [growth_delay()].
#' @param n The number of rows of each table to print.
#' @param ... Not used.
#' @return `print()` returns `x` invisibly; `summary()` a list of data frames
#'   with class `summary.me_result`.
#' @name result-printing
#' @examples
#' gd <- growth_delay(master_synthetic_data, volume = 1000, cage_column = "Cage",
#'                    reference_group = "Vehicle", n_boot = 100, seed = 1)
#' gd
#' summary(gd)$delay_table
NULL

#' @rdname result-printing
#' @export
print.me_result <- function(x, n = 10L, ...) {
  spec <- me_result_spec(x)
  cat("<mouseExperiment: ", spec$title(x), ">\n", sep = "")
  about <- if (is.function(spec$about)) tryCatch(spec$about(x), error = function(e) NULL)
  about <- about[nzchar(about %||% "")]
  if (length(about)) cat(paste0("  ", about, "\n"), sep = "")
  for (f in names(spec$tables)) {
    df <- x[[f]]
    if (!is.data.frame(df) || !nrow(df)) next
    cat("\n", spec$tables[[f]], "  ($", f, ")\n", sep = "")
    print(utils::head(as.data.frame(df), n), digits = 4, row.names = me_named_rows(df))
    if (nrow(df) > n) cat("  ... and", nrow(df) - n, "more rows\n")
  }
  cat("\nEvery field: names(x). The tables, whole: summary(x).\n")
  invisible(x)
}

#' @rdname result-printing
#' @export
summary.me_result <- function(object, ...) {
  spec <- me_result_spec(object)
  tabs <- Filter(is.data.frame, lapply(names(spec$tables), function(f) object[[f]]))
  names(tabs) <- names(spec$tables)[vapply(names(spec$tables), function(f) is.data.frame(object[[f]]), NA)]
  structure(tabs, class = "summary.me_result", title = spec$title(object))
}

#' @rdname result-printing
#' @export
print.summary.me_result <- function(x, ...) {
  cat("<mouseExperiment summary: ", attr(x, "title"), ">\n", sep = "")
  for (f in names(x)) {
    cat("\n$", f, "\n", sep = "")
    print(as.data.frame(x[[f]]), digits = 4, row.names = me_named_rows(x[[f]]))
  }
  invisible(x)
}
