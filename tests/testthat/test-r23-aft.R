# Accelerated failure time models (CODE_REVIEW.md Round 23, R23.11, R23-F):
# aft_survival() and the Survival tab's AFT option.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

# Weibull AFT times: an arm's times are its time ratio times the reference's,
# censored at day 40.
sim_aft <- function(n = 40, ratios = c(Vehicle = 1, Slow = 2, Same = 1), shape = 3, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(names(ratios), function(a) {
    t <- ratios[[a]] * 12 * stats::rweibull(n, shape = shape, scale = 1)
    data.frame(ID = paste0(a, seq_len(n)), Treatment = a, Time = pmin(t, 40),
               Event = as.integer(t <= 40), stringsAsFactors = FALSE)
  }))
}

test_that("R23.11: the time ratios are survreg's, with Wald intervals and corrected p-values", {
  d <- sim_aft()
  a <- aft_survival(d, reference_group = "Vehicle", distribution = "weibull")
  expect_s3_class(a, c("me_aft", "me_result"))
  fit <- survival::survreg(survival::Surv(Time, Event) ~ stats::relevel(factor(Treatment), "Vehicle"),
                           data = d, dist = "weibull")
  co <- summary(fit)$table
  tr <- a$time_ratios
  expect_identical(tr$Treatment, c("Vehicle", "Same", "Slow"))
  expect_identical(tr$Time_Ratio[1], 1)
  expect_equal(tr$Time_Ratio[2:3], unname(exp(co[2:3, "Value"])))
  expect_equal(tr$Lower_CL[3], unname(exp(co[3, "Value"] - stats::qnorm(0.975) * co[3, "Std. Error"])))
  expect_equal(tr$P_Value[2:3], unname(co[2:3, "p"]))
  expect_equal(tr$P_Adjusted[2:3], stats::p.adjust(tr$P_Value[2:3], "holm"))
  # The truth is recovered: Slow takes about twice as long; Same is not different.
  expect_true(tr$Lower_CL[3] < 2 && 2 < tr$Upper_CL[3])
  expect_true(tr$Lower_CL[2] < 1 && 1 < tr$Upper_CL[2])
  expect_lt(tr$P_Adjusted[3], 0.001)
  # The predicted median of a Weibull with shape 3 and scale 12 x ratio.
  expect_equal(tr$Median_Predicted[1], 12 * log(2)^(1 / 3), tolerance = 0.1)
  # The fitted curves start near 1 and fall.
  fc <- a$fitted_curves
  expect_setequal(unique(fc$Treatment), c("Vehicle", "Same", "Slow"))
  expect_true(all(diff(fc$Survival[fc$Treatment == "Vehicle"]) <= 0))
  expect_s3_class(g("aft_fit_plot")(a), "ggplot")
  expect_match(paste(capture.output(print(a)), collapse = "\n"), "Accelerated failure time model (weibull)",
               fixed = TRUE)
})

test_that("R23.11: 'auto' uses the distribution with the lowest AIC", {
  set.seed(3)
  d <- do.call(rbind, lapply(c(Vehicle = 1, Drug = 1.6), function(r) data.frame(Time = r * stats::rlnorm(60, 2.5, 0.6))))
  d$Treatment <- rep(c("Vehicle", "Drug"), each = 60)
  d$Event <- 1L
  d$ID <- seq_len(nrow(d))
  a <- aft_survival(d, reference_group = "Vehicle")
  dc <- a$distribution_comparison
  expect_identical(dc$Distribution, c("weibull", "lognormal", "loglogistic"))
  expect_identical(a$distribution, dc$Distribution[which.min(dc$AIC)])
  expect_identical(dc$Used, dc$Distribution == a$distribution)
  expect_identical(nrow(aft_survival(d, reference_group = "Vehicle", distribution = "lognormal")$distribution_comparison), 1L)
})

test_that("R23.11: no events in an arm, in the reference, doses, cages, repeated rows", {
  d <- sim_aft(n = 10)
  d$Event[d$Treatment == "Slow"] <- 0L
  a <- aft_survival(d, reference_group = "Vehicle")
  r <- a$time_ratios[a$time_ratios$Treatment == "Slow", ]
  expect_true(is.na(r$Time_Ratio))
  expect_match(r$Note, "no events", fixed = TRUE)
  expect_false(is.na(a$time_ratios$Time_Ratio[a$time_ratios$Treatment == "Same"]))
  d2 <- d; d2$Event[d2$Treatment == "Vehicle"] <- 0L
  expect_error(aft_survival(d2, reference_group = "Vehicle"), "reference arm .* has no events")
  # One treatment at two doses: two arms.
  d3 <- sim_aft(n = 8, ratios = c(Vehicle = 1, Drug = 1.5, Drug2 = 2))
  d3$Dose <- ifelse(d3$Treatment == "Vehicle", 0, ifelse(d3$Treatment == "Drug", 10, 30))
  d3$Treatment[d3$Treatment == "Drug2"] <- "Drug"
  a3 <- suppressWarnings(aft_survival(d3, dose_column = "Dose", reference_group = "Vehicle"))
  expect_identical(a3$time_ratios$Treatment, c("Vehicle", "Drug 10", "Drug 30"))
  expect_warning(aft_survival(d, reference_group = "Vehicle", randomisation_unit = "cage"),
                 "treats animals as independent")
  expect_error(aft_survival(rbind(d, d[1, ]), reference_group = "Vehicle"), "one row per")
  d4 <- d; d4$Time[1] <- 0
  expect_error(aft_survival(d4, reference_group = "Vehicle"), "positive times")
})

test_that("R23.11: the Survival tab fits AFT, shows and downloads it, and its R script reproduces it", {
  rv <- shiny::reactiveValues(
    uploaded_data = utils::read.csv(system.file("sample_data", "all_features_synthetic_data.csv",
                                                package = "mouseExperiment"), stringsAsFactors = FALSE),
    data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_weight = "Weight", col_dose = "Dose", col_survival_event = "Survival_Censor",
                    control_group = "Vehicle", volume_units = "mm3", randomisation_unit = "mouse"))
  shiny::isolate(g("ensure_preprocessed")(rv))
  suppressWarnings(shiny::testServer(g("survival_server"), args = list(rv = rv), {
    session$setInputs(surv_model_type = "aft", surv_aft_dist = "auto", surv_aft_p_adjust = "bonferroni")
    session$setInputs(run_surv_analysis = 1)
    res <- rv$surv_results
    expect_false(isTRUE(res$error))
    expect_identical(res$model_type_used, "aft")
    expect_identical(res$p_adjust_method, "bonferroni")
    expect_false(is.null(res$km_fit))                  # the Kaplan-Meier plot and risk table still work
    head <- paste(as.character(output$surv_cox_table), collapse = "")
    expect_match(head, "Time_Ratio", fixed = TRUE)
    expect_match(head, "Lower 95% CI", fixed = TRUE)
    expect_match(output$surv_summary, "Time ratios against the reference", fixed = TRUE)
    expect_match(paste(as.character(output$surv_aft_dist_table), collapse = ""), "<th>AIC", fixed = TRUE)
    expect_identical(res$distribution_comparison$Distribution, c("weibull", "lognormal", "loglogistic"))
    ef <- g("surv_effects")(res)
    expect_identical(ef$measure, "Time Ratio")
    groups <- c("tables", "plots", "diagnostics", "data", "report")
    ins <- stats::setNames(rep(list(character(0)), length(groups)), paste0("surv_dl_group_", groups))
    ins$surv_dl_group_tables <- "results"
    ins$surv_dl_group_diagnostics <- c("aft_dists", "aft_fit")
    ins$surv_dl_group_code <- "script"
    do.call(session$setInputs, ins)
    dir <- tempfile("aft_")
    utils::unzip(output$surv_dl_zip, exdir = dir)
    files <- list.files(dir, recursive = TRUE)
    expect_true(all(c("tables/survival_results.csv", "diagnostics/aft_distributions.csv",
                      "diagnostics/aft_fit.png", "analysis.R") %in% files))
    script <- readLines(file.path(dir, "analysis.R"))
    expect_true(any(grepl("result <- aft_survival(", script, fixed = TRUE)))
    env <- new.env(parent = globalenv())
    old <- setwd(dir); on.exit(setwd(old))
    suppressMessages(suppressWarnings(utils::capture.output(
      sys.source("analysis.R", envir = env, keep.source = FALSE))))
    expect_equal(env$result$time_ratios$Time_Ratio, res$time_ratios$Time_Ratio)
    expect_equal(env$result$time_ratios$P_Adjusted, res$time_ratios$P_Adjusted)
  }))
})

test_that("R23.11: a failed proportional-hazards check points to the AFT model", {
  rv <- shiny::reactiveValues(settings = list())
  shiny::testServer(g("survival_server"), args = list(rv = rv), {
    tab <- matrix(c(9.1, 3, 0.028), nrow = 1, dimnames = list("GLOBAL", c("chisq", "df", "p")))
    rv$surv_results <- list(error = FALSE, ph_test = list(table = tab), method_used = "cox")
    expect_match(paste(as.character(output$surv_ph_test_banner$html), collapse = ""),
                 "Statistical Model: AFT", fixed = TRUE)
    tab[1, "p"] <- 0.4
    rv$surv_results <- list(error = FALSE, ph_test = list(table = tab), method_used = "cox")
    session$flushReact()
    expect_identical(paste(as.character(output$surv_ph_test_banner$html), collapse = ""), "")
  })
})
