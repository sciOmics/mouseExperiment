# The internal animal key ("aPD1|||L|||3") showed as the first column of the
# Growth Rates table (user report 2026-10-25, 0.49.1). Per-animal results
# built with split() and rbind() took the keys as row names, and DT shows row
# names. No result carries keys as row names now, and the dashboard's tables
# show row names only when they say something.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

keyed <- function(x, path = "result", depth = 0) {
  if (depth > 4 || is.null(x)) return(character(0))
  if (is.data.frame(x)) return(if (any(grepl("|||", rownames(x), fixed = TRUE))) path else character(0))
  if (is.list(x) && !inherits(x, c("ggplot", "merMod", "lm")))
    return(unlist(lapply(names(x), function(n) keyed(x[[n]], paste0(path, "$", n), depth + 1))))
  character(0)
}

test_that("tumour growth results name no row by its animal key", {
  q <- function(e) suppressWarnings(suppressMessages(e))
  for (mt in c("lme4", "auc")) {
    r <- q(tumor_growth_statistics(master_synthetic_data, model_type = mt, reference_group = "Vehicle",
                                   cage_column = "Cage", plots = FALSE))
    expect_identical(keyed(unclass(r)), character(0), info = mt)
    expect_identical(rownames(r$growth_rates), as.character(seq_len(nrow(r$growth_rates))), info = mt)
    expect_null(names(r$growth_rates$growth_rate))
  }
  # The random effects, whose rows were named only by the key, now carry the
  # animal as columns.
  re <- r <- q(tumor_growth_statistics(master_synthetic_data, reference_group = "Vehicle",
                                       cage_column = "Cage", plots = FALSE))$diagnostics$random_effects
  expect_identical(names(re$intercepts)[1:3], c("Treatment", "ID", "Cage"))
  expect_identical(nrow(re$intercepts), 48L)
  expect_null(names(re$slopes))
  # The therapeutic window's per-animal weight loss had the key as a column.
  tw <- q(therapeutic_window_metric(master_synthetic_data, reference_group = "Vehicle",
                                    cage_column = "Cage", volume_units = "mm3", n_boot = 0))
  expect_identical(names(tw$weight_loss_data)[1:3], c("Treatment", "ID", "Cage"))
})

test_that("row names show when they say something: terms yes, numbers and keys no", {
  named <- g("me_named_rows")
  expect_true(named(data.frame(x = 1:2, row.names = c("Treatment:Day", "Day"))))
  expect_false(named(data.frame(x = 1:3)))
  expect_false(named(data.frame(x = 1:3, row.names = c("6", "1", "2"))))
  expect_false(named(data.frame(x = 1:2, row.names = c("aPD1|||L|||3", "aPD1|||R|||3"))))
  expect_false(named(data.frame()))
})

test_that("the Growth Rates table shows no animal key and no row-number column", {
  rv <- shiny::reactiveValues(
    uploaded_data = master_synthetic_data, data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference", tg_p_adjust = "holm")
    session$setInputs(run_tg_analysis = 1)
    # The tables are server-side: the widget holds the header, the rows come
    # later. The first column is the treatment, not a row-name column.
    first_th <- function(widget) regmatches(widget, regexpr("<th>[^<]*<", widget))
    widget <- paste(as.character(output$tg_growth_rates_table), collapse = "")
    expect_identical(first_th(widget), "<th>Treatment<")
    expect_false(grepl("|||", widget, fixed = TRUE))
    expect_identical(rownames(rv$tg_results$growth_rates)[1:2], c("1", "2"))
    # The ANOVA keeps its row-name column: its terms.
    anova <- paste(as.character(output$tg_anova_table), collapse = "")
    expect_identical(first_th(anova), "<th> <")
  }))
})

# ---- a Stan compile killed for lack of memory (user report 2026-10-25) -------------

test_that("a compile killed for lack of memory reads as that, not as compiler output", {
  log <- paste("stan::mcmc::dense_e_point&) [with Model = model437ca635f5__namespace::model437ca635f5_ ...",
               "warning: ignoring attributes on template argument ... g++: fatal error: Killed signal",
               "terminated program cc1plus compilation terminated. make: *** [Makeconf:211: file43.o] Error 1")
  msg <- g("translate_error")(log)
  expect_match(msg, "ran out of memory", fixed = TRUE)
  expect_false(grepl("cc1plus", msg, fixed = TRUE))
  # Both Bayesian failure paths translate.
  for (f in c("tg_resume.R", "surv_resume.R"))
    expect_false(grepl('paste("Bayesian analysis failed:", err_msg)', dash_source("R", f), fixed = TRUE),
                 info = f)
})
