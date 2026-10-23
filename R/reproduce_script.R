# =============================================================================
# "R script" in Download Results: a script that repeats the run with the
# package, from the data exactly as the dashboard analysed it
# (CODE_REVIEW.md R20-Z).
#
# The call is recorded when it is made -- the function and every argument the
# module sent, explicit NULLs included -- and the script is written from that
# record. A script rebuilt from the sidebar afterwards would drift from what
# actually ran, as the HTML report did before it recorded the run's settings.
# =============================================================================

#' Record a package call as the dashboard makes it
#'
#' @param fun The function's name in mouseExperiment.
#' @param args The argument list passed to it, exactly as passed.
#' @param data_arg The argument that holds the data. Its value is kept as
#'   `data` (the frame the script's CSV is written from), not deparsed.
#' @param exprs Named calls to write in place of argument values that do not
#'   deparse readably, e.g. `permutation_test = quote(perm_spec(...))`.
#' @return A list: fun, args (without the data), data_arg, data, exprs and the
#'   time of the run.
#' @keywords internal
record_call <- function(fun, args, data_arg = "df", exprs = list()) {
  stopifnot(is.character(fun), length(fun) == 1L, data_arg %in% names(args))
  list(fun = fun, args = args[setdiff(names(args), data_arg)], data_arg = data_arg,
       data = args[[data_arg]], exprs = exprs, when = Sys.time())
}

#' One argument value as R source
#'
#' Single line, with integers, NAs and names kept, so the script passes the
#' same values: `999L` stays an integer and `NULL` stays `NULL`.
#' @keywords internal
script_value <- function(x) {
  if (is.language(x)) return(paste(deparse(x, width.cutoff = 500L), collapse = " "))
  paste(deparse(x, width.cutoff = 500L, control = c("keepNA", "keepInteger", "niceNames")),
        collapse = " ")
}

#' A recorded call as R source: `name <- fun(data = df, arg = value, ...)`
#'
#' Every recorded argument is written, NULLs included: a NULL the dashboard
#' passed can differ from the default (`cage_column` defaults to "Cage").
#' @param rec A [record_call()] record.
#' @param name The variable the result is assigned to.
#' @param data_name The variable holding the data in the script.
#' @param override Named values that replace or add arguments.
#' @param suffix Text after the call, e.g. "$permutation_test".
#' @keywords internal
script_call <- function(rec, name, data_name = "df", override = list(), suffix = "") {
  args <- rec$args
  args[names(override)] <- override
  vals <- vapply(names(args), function(a) {
    if (a %in% names(rec$exprs)) script_value(rec$exprs[[a]]) else script_value(args[[a]])
  }, "")
  lines <- c(sprintf("  %s = %s", rec$data_arg, data_name), sprintf("  %s = %s", names(args), vals))
  c(sprintf("%s <- %s(", name, rec$fun),
    paste0(lines, c(rep(",", length(lines) - 1L), "")),
    paste0(")", suffix))
}

#' The opening lines every script shares
#'
#' @param title What was run, e.g. "Tumor Growth (LME4)".
#' @param rec The main call's record, for the time of the run.
#' @param data_file The data's path inside the ZIP, without the extension:
#'   the script reads the `.rds` (exact) and falls back to the `.csv`.
#' @keywords internal
script_header <- function(title, rec, data_file) {
  ver <- tryCatch(utils::packageDescription("mouseExperiment")$Version, error = function(e) NULL) %||% "?"
  sha <- tryCatch(utils::packageDescription("mouseExperiment")$RemoteSha, error = function(e) NULL)
  install <- if (length(sha) && grepl("^[0-9a-f]{40}$", sha)) {
    sprintf('#   remotes::install_github("sciOmics/mouseExperiment@%s")', sha)
  } else {
    '#   remotes::install_github("sciOmics/mouseExperiment")   # then check the version below'
  }
  c(sprintf("# Repeats the %s analysis run in the mouseExperiment dashboard", title),
    sprintf("# on %s, with mouseExperiment %s%s.", format(rec$when, "%Y-%m-%d %H:%M %Z"), ver,
            if (length(sha) && nzchar(sha)) paste0(" (commit ", substr(sha, 1, 12), ")") else ""),
    "#",
    sprintf("# It reads %s.rds, the data exactly as the dashboard analysed it:", data_file),
    "# after the subset filter, the day and volume calculations and the arm labels.",
    sprintf("# %s.csv is the same table as text, which keeps 15 significant digits;", data_file),
    "# from it, results can differ in about the 6th digit.",
    "# Run it from the folder the ZIP was unzipped into (setwd() there first).",
    "#",
    "# The same package version gives the same results. To install it:",
    install,
    "#",
    "# Plots are not repeated with the dashboard's styling (colours, labels, group",
    "# order); the plots in this ZIP are the ones the dashboard drew.",
    "",
    "library(mouseExperiment)",
    sprintf('if (as.character(packageVersion("mouseExperiment")) != "%s")', ver),
    sprintf('  warning("This script was written by mouseExperiment %s; this is ",', ver),
    '          packageVersion("mouseExperiment"), ". Results may differ.")',
    "",
    sprintf('df <- if (file.exists("%s.rds")) readRDS("%s.rds") else', data_file, data_file),
    sprintf('  utils::read.csv("%s.csv", check.names = FALSE, encoding = "UTF-8")', data_file),
    "")
}

#' The Tumor Growth run as an R script
#'
#' @param res The stored Tumor Growth result, with its `call_record` (and
#'   `perm_call_record` when the randomisation test ran).
#' @param data_file The analysed data's path inside the ZIP, without the
#'   extension.
#' @return Character lines, or NULL when the run left no record.
#' @keywords internal
tg_reproduce_script <- function(res, data_file = "data/analysed_data") {
  rec <- res$call_record
  if (is.null(rec)) return(NULL)
  bayes <- identical(rec$fun, "bayesian_tumor_growth")
  title <- paste0("Tumor Growth (", switch(res$model_type_used %||% "lme4",
                                           lme4 = "LME4", auc = "AUC", bayes = "Bayesian LMM",
                                           res$model_type_used), ")")
  shown <- c(anova = "The model's ANOVA table",
             pairwise_comparisons = "Comparisons between arms, with the correction selected",
             treatment_effects = "Estimated effects by arm",
             posterior_summary = "Posterior summaries of the parameters")
  shown <- shown[vapply(names(shown), function(f) is.data.frame(res[[f]]) && nrow(res[[f]]) > 0L, TRUE)]
  out <- c(script_header(title, rec, data_file),
           if (bayes) c("# MCMC: the seed is the package default (42), so the same version on the",
                        "# same platform draws the same chains; another compiler or Stan version can",
                        "# move the draws slightly.", ""),
           script_call(rec, "result"),
           "",
           unlist(lapply(names(shown), function(f) c(paste0("# ", shown[[f]]),
                                                     sprintf("print(result$%s)", f)))))
  prec <- res$perm_call_record
  if (!is.null(prec) && !is.null(res$permutation_test)) {
    out <- c(out, "",
             "# The randomisation test of the treatment x time interaction, as the dashboard",
             "# ran it: the same call with the test, without the model and diagnostics.",
             script_call(prec, "perm", suffix = "$permutation_test"),
             'print(perm[c("method", "p_value", "n_perm", "min_attainable_p")])')
  } else if (!is.null(prec)) {
    out <- c(out, "", "# The randomisation test was requested but had not finished when this was",
             "# downloaded; download again once its result is shown.")
  }
  c(out, "", "sessionInfo()")
}
