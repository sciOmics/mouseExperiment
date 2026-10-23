# =============================================================================
# "R script" in Download Results: a script that repeats the run with the
# package, from the data exactly as the dashboard analysed it
# (CODE_REVIEW.md R20-Z, R20-AA).
#
# The call is recorded when it is made -- the function and every argument the
# module sent, explicit NULLs included -- and the script is written from that
# record. A script rebuilt from the sidebar afterwards would drift from what
# actually ran, as the HTML report did before it recorded the run's settings.
#
# A script is a list of steps, one per recorded call. Each step reads its data
# from a file in data/: an .rds that holds the frame exactly, and a .csv of
# the same table for other software (CSV text keeps 15 significant digits, and
# R's parser is not correctly rounded on every platform).
# =============================================================================

#' Record a package call as the dashboard makes it
#'
#' @param fun The function's name in mouseExperiment.
#' @param args The argument list passed to it, exactly as passed.
#' @param data_arg The argument that holds the data, or NULL when the call
#'   takes none (power analysis). Its value is kept as `data` (the frame the
#'   script's data file is written from), not deparsed.
#' @param exprs Named calls to write in place of argument values that do not
#'   deparse readably, e.g. `permutation_test = quote(perm_spec(...))`.
#' @return A list: fun, args (without the data), data_arg, data, exprs and the
#'   time of the run.
#' @keywords internal
record_call <- function(fun, args, data_arg = "df", exprs = list()) {
  stopifnot(is.character(fun), length(fun) == 1L,
            is.null(data_arg) || data_arg %in% names(args))
  list(fun = fun,
       args = if (is.null(data_arg)) args else args[setdiff(names(args), data_arg)],
       data_arg = data_arg,
       data = if (!is.null(data_arg)) args[[data_arg]],
       exprs = exprs, when = Sys.time())
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
  lines <- c(if (!is.null(rec$data_arg)) sprintf("  %s = %s", rec$data_arg, data_name),
             sprintf("  %s = %s", names(args), vals))
  c(sprintf("%s <- %s(", name, rec$fun),
    paste0(lines, c(rep(",", length(lines) - 1L), "")),
    paste0(")", suffix))
}

#' One recorded call in a script
#'
#' @param rec A [record_call()] record.
#' @param name The variable the result is assigned to.
#' @param data_file The step's data in the ZIP, without the extension. Steps
#'   that share a file must have the same data.
#' @param data_var The variable that holds that data in the script.
#' @param comment Lines written above the call (without "# ").
#' @param prints Named: result field -> what it is; each is printed.
#' @param after Lines written after the prints.
#' @param suffix,override As for [script_call()].
#' @keywords internal
script_step <- function(rec, name, data_file = "data/analysed_data", data_var = "df",
                        comment = character(0), prints = character(0), after = character(0),
                        suffix = "", override = list()) {
  list(rec = rec, name = name, data_file = data_file, data_var = data_var,
       comment = comment, prints = prints, after = after, suffix = suffix,
       override = override)
}

#' The opening lines every script shares
#'
#' @param title What was run, e.g. "Tumor Growth (LME4)".
#' @param when The time of the run.
#' @param data_files The data files the script reads, without extensions.
#' @keywords internal
script_header <- function(title, when, data_files) {
  ver <- tryCatch(utils::packageDescription("mouseExperiment")$Version, error = function(e) NULL) %||% "?"
  sha <- tryCatch(utils::packageDescription("mouseExperiment")$RemoteSha, error = function(e) NULL)
  install <- if (length(sha) && grepl("^[0-9a-f]{40}$", sha)) {
    sprintf('#   remotes::install_github("sciOmics/mouseExperiment@%s")', sha)
  } else {
    '#   remotes::install_github("sciOmics/mouseExperiment")   # then check the version below'
  }
  data_lines <- if (length(data_files)) {
    c(sprintf("# It reads %s, the data exactly as the dashboard analysed %s:",
              paste0(data_files, ".rds", collapse = " and "),
              if (length(data_files) > 1L) "them" else "it"),
      "# after the subset filter, the day and volume calculations and the arm labels.",
      "# The .csv beside each is the same table as text, which keeps 15 significant",
      "# digits; from it, results can differ in about the 6th digit.",
      "# Run it from the folder the ZIP was unzipped into (setwd() there first).")
  } else {
    "# It needs no data: the calculation uses only the design values below."
  }
  c(sprintf("# Repeats the %s analysis run in the mouseExperiment dashboard", title),
    sprintf("# on %s, with mouseExperiment %s%s.", format(when, "%Y-%m-%d %H:%M %Z"), ver,
            if (length(sha) && nzchar(sha)) paste0(" (commit ", substr(sha, 1, 12), ")") else ""),
    "#",
    data_lines,
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
    "")
}

#' A script and the data it reads
#'
#' @param title What was run.
#' @param steps A list of [script_step()]s, in the order they are written;
#'   NULLs are dropped.
#' @param notes Lines written after the header (without "# ").
#' @return list(lines, data): the script, and the data frames to write,
#'   named by their file (without extension).
#' @keywords internal
script_bundle <- function(title, steps, notes = character(0)) {
  steps <- Filter(Negate(is.null), steps)
  stopifnot(length(steps) > 0L)
  with_data <- Filter(function(s) !is.null(s$rec$data_arg), steps)
  first_for <- function(f) Filter(function(s) identical(s$data_file, f), with_data)[[1L]]
  files <- unique(vapply(with_data, `[[`, "", "data_file"))
  data <- lapply(stats::setNames(files, files), function(f) first_for(f)$rec$data)
  reads <- if (length(files)) {
    c("# The .rds holds the data exactly; the .csv is the fallback.",
      "read_analysed <- function(path) {",
      '  if (file.exists(paste0(path, ".rds"))) return(readRDS(paste0(path, ".rds")))',
      '  utils::read.csv(paste0(path, ".csv"), check.names = FALSE, encoding = "UTF-8")',
      "}",
      sprintf('%s <- read_analysed("%s")', vapply(files, function(f) first_for(f)$data_var, ""), files),
      "")
  }
  body <- unlist(lapply(steps, function(s) {
    c(if (length(s$comment)) paste0("# ", s$comment),
      script_call(s$rec, s$name, s$data_var, s$override, s$suffix),
      unlist(lapply(names(s$prints), function(f)
        c(paste0("# ", s$prints[[f]]), sprintf("print(%s$%s)", s$name, f)))),
      s$after,
      "")
  }))
  lines <- c(script_header(title, steps[[1L]]$rec$when, files),
             if (length(notes)) c(paste0("# ", notes), ""),
             reads, body, "sessionInfo()")
  list(lines = lines, data = data)
}

#' Write a script bundle into a directory
#' @param bundle A [script_bundle()].
#' @param dir The ZIP's top level.
#' @return The files written, relative to `dir`.
#' @keywords internal
write_script_bundle <- function(bundle, dir) {
  writeLines(bundle$lines, file.path(dir, "analysis.R"))
  out <- "analysis.R"
  for (f in names(bundle$data)) {
    dir.create(dirname(file.path(dir, f)), recursive = TRUE, showWarnings = FALSE)
    saveRDS(bundle$data[[f]], file.path(dir, paste0(f, ".rds")))
    utils::write.csv(bundle$data[[f]], file.path(dir, paste0(f, ".csv")), row.names = FALSE)
    out <- c(out, paste0(f, ".rds"), paste0(f, ".csv"))
  }
  out
}

#' The Download Results item for a module's R script
#' @param bundle_fn A function returning the [script_bundle()], or NULL.
#' @keywords internal
script_dl_item <- function(bundle_fn) {
  dl_item("script", "code", "R script, with the data it reads", "analysis.R", "files",
          function(dir) {
            b <- bundle_fn()
            if (is.null(b)) character(0) else write_script_bundle(b, dir)
          },
          "Repeats this run with the package; reads the files in data/ written with it")
}

#' The Tumor Growth run as a script bundle
#'
#' @param res The stored Tumor Growth result, with its `call_record` (and
#'   `perm_call_record` when the randomisation test ran).
#' @return A [script_bundle()], or NULL when the run left no record.
#' @keywords internal
tg_script_bundle <- function(res) {
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
  prec <- res$perm_call_record
  steps <- list(
    script_step(rec, "result", prints = shown),
    if (!is.null(prec) && !is.null(res$permutation_test))
      script_step(prec, "perm", suffix = "$permutation_test",
                  comment = c("The randomisation test of the treatment x time interaction, as the dashboard",
                              "ran it: the same call with the test, without the model and diagnostics."),
                  after = 'print(perm[c("method", "p_value", "n_perm", "min_attainable_p")])'))
  notes <- c(
    if (bayes) c("MCMC: the seed is the package default (42), so the same version on the",
                 "same platform draws the same chains; another compiler or Stan version can",
                 "move the draws slightly."),
    if (!is.null(prec) && is.null(res$permutation_test))
      c("The randomisation test was requested but had not finished when this was",
        "downloaded; download again once its result is shown."))
  script_bundle(title, steps, notes)
}

#' The Survival run as a script bundle
#' @param res The stored Survival result, with its `call_record`.
#' @return A [script_bundle()], or NULL when the run left no record.
#' @keywords internal
surv_script_bundle <- function(res) {
  rec <- res$call_record
  if (is.null(rec)) return(NULL)
  bayes <- identical(rec$fun, "bayesian_survival")
  title <- if (bayes) "Survival (Bayesian AFT)" else "Survival"
  notes <- c(
    if (identical(res$event_src_used, "volume"))
      c(sprintf("The event and time columns (%s, %s) were derived by the dashboard:",
                rec$args$censor_column %||% rec$args$event_column, rec$args$time_column),
        sprintf("an event is the first day a tumour reached the volume limit of %s.",
                format(res$vol_threshold_used)),
        "They are in the data."),
    if (isTRUE(res$death_events_added > 0))
      c(sprintf("%d deaths from other causes were counted as events, from the removal", res$death_events_added),
        "reasons; the event column in the data already includes them."),
    if (bayes) c("MCMC: the seed is the package default, so the same version on the same",
                 "platform draws the same chains; another compiler or Stan version can",
                 "move the draws slightly."))
  shown <- c(results = if (bayes) "Posterior treatment effects by arm"
                       else "Each arm against the reference, with the correction selected")
  script_bundle(title, list(script_step(rec, "result", prints = shown)), notes)
}

#' The Toxicity run as a script bundle
#'
#' A body-weight run makes two calls on the same data (the model and the
#' weight-loss threshold); a therapeutic-window run makes one.
#' @param res The stored Toxicity result, with its `call_records`.
#' @return A [script_bundle()], or NULL when the run left no record.
#' @keywords internal
tox_script_bundle <- function(res) {
  recs <- res$call_records
  if (is.null(recs)) return(NULL)
  if (identical(res$model_type, "et")) {
    if (is.null(recs$twm)) return(NULL)
    return(script_bundle("Toxicity (therapeutic window)", list(
      script_step(recs$twm, "twm", prints = c(
        window_table = "TGI and worst weight loss by arm, with tolerability",
        weight_loss_data = "Each animal's worst weight loss"))),
      c("The intervals come from a bootstrap with the seed in the call (boot_seed).")))
  }
  ok <- function(x) !is.null(x) && !isTRUE(x$error)
  steps <- list(
    if (!is.null(recs$bw) && ok(res$mixed_model))
      script_step(recs$bw, "bw", comment = "The body-weight mixed model", prints = c(
        fixed_effects = "Fixed effects of the body-weight model",
        emmeans_table = "Each arm's rate of weight change",
        pairwise_comparisons = "Each arm's rate against the reference's")),
    if (!is.null(recs$wl) && ok(res$weight_loss))
      script_step(recs$wl, "wl", comment = "Time to the weight-loss threshold", prints = c(
        event_data = "Each animal's time to the weight-loss threshold")))
  steps <- Filter(Negate(is.null), steps)
  if (!length(steps)) return(NULL)
  script_bundle("Toxicity (body weight)", steps)
}

# The prints a step makes: only the fields this run's result holds.
script_prints <- function(res, prints) {
  prints[vapply(names(prints), function(f) {
    x <- res[[f]]
    !is.null(x) && (!is.data.frame(x) || nrow(x) > 0L)
  }, TRUE)]
}

#' The Drug Synergy run as a script bundle
#' @param res The stored Drug Synergy result, with its `call_record`.
#' @return A [script_bundle()], or NULL when the run left no record.
#' @keywords internal
syn_script_bundle <- function(res) {
  rec <- res$call_record
  if (is.null(rec)) return(NULL)
  over_time <- identical(rec$fun, "analyze_drug_synergy_over_time")
  prints <- script_prints(res, if (over_time) c(
    synergy_summary = "Synergy metrics for each evaluable day") else c(
    synergy_metrics = "Bliss expectation, observed TGI and the difference",
    summary = "Each arm's tumour growth inhibition",
    synergy_ci = "95% intervals of the synergy metrics"))
  script_bundle(if (over_time) "Drug Synergy (every evaluable day)" else "Drug Synergy",
                list(script_step(rec, "result", prints = prints)),
                "The intervals come from a bootstrap with the seed in the call (boot_seed).")
}

#' The Dose Response run as a script bundle
#' @param res The stored Dose Response result, with its `call_record`.
#' @return A [script_bundle()], or NULL when the run left no record.
#' @keywords internal
dr_script_bundle <- function(res) {
  rec <- res$call_record
  if (is.null(rec)) return(NULL)
  prints <- script_prints(res, c(
    summary_table = "Tumour volume by dose group",
    tgi_table = "Each dose's TGI with its interval"))
  script_bundle("Dose Response", list(script_step(
    rec, "result", prints = prints,
    after = c("# Trend tests, the EC50 with its interval, and the curve's parameters",
              "str(result$statistics, max.level = 1)"))))
}

#' The Power Analysis run as a script bundle
#'
#' The analytic mode calls apriori_power_analysis() once per effect size and
#' stacks the scenario tables, as the dashboard does; the LMM mode is one
#' seeded simulation. Neither reads data.
#' @param res The stored Power result, with its `call_records`.
#' @return A [script_bundle()], or NULL when the run left no record.
#' @keywords internal
pw_script_bundle <- function(res) {
  recs <- res$call_records
  if (!length(recs)) return(NULL)
  if (identical(res$mode, "apriori_lmm")) {
    return(script_bundle("Power Analysis (LMM simulation)", list(script_step(
      recs[[1L]], "sim", prints = c(power_table = "Power by N per group, from the simulations"))),
      "The simulation is seeded (seed in the call), so the same version gives the same power."))
  }
  names <- if (length(recs) == 1L) "power" else paste0("power_", seq_along(recs))
  steps <- Map(function(r, n) script_step(r, n), recs, names)
  combine <- function(field) sprintf("%s <- do.call(rbind, lapply(list(%s), `[[`, \"%s\"))",
                                     field, paste(names, collapse = ", "), field)
  after <- if (length(recs) == 1L) {
    c("# Required N or achieved power for each scenario", "print(power$scenario_table)",
      if (!is.null(res$sensitivity_table))
        c("# How the answer moves with the pooled SD", "print(power$sensitivity_table)"))
  } else {
    c("# One call per effect size, stacked as the dashboard shows them",
      combine("scenario_table"), "print(scenario_table)")
  }
  steps[[length(steps)]]$after <- after
  script_bundle("Power Analysis (analytic)", steps)
}
