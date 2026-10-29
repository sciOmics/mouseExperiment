# Round 22 proposals F2-F3 (CODE_REVIEW.md R22-E, 0.43.0): %T/C beside TGI and
# tumour growth delay, with its test and specific growth delay (R23.2,
# 0.52.0). The mRECIST response categories (F1) were removed in 0.52.0
# (R23.1).

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

# ---- F3: %T/C --------------------------------------------------------------------

test_that("F3: the therapeutic window and dose-response report %T/C = 100 - TGI, interval mirrored", {
  m <- master_synthetic_data
  tw <- suppressWarnings(suppressMessages(therapeutic_window_metric(m, reference_group = "Vehicle",
    cage_column = "Cage", volume_units = "mm3", n_boot = 100, seed = 1)))
  w <- tw$window_table
  expect_equal(w$TC, 100 - w$TGI)
  expect_equal(w$TC_Lower, 100 - w$TGI_Upper)
  expect_equal(w$TC_Upper, 100 - w$TGI_Lower)
  expect_identical(match(c("TGI_Upper", "TC", "TC_Lower", "TC_Upper"), names(w)),
                   match("TGI_Upper", names(w)) + 0:3)
  dr <- suppressWarnings(suppressMessages(dose_response_statistics(m, cage_column = "Cage",
    treatments = c("Drug_A Low", "Drug_A Mid", "Drug_A High"), reference_group = "Vehicle",
    verbose = FALSE)))
  tt <- dr$tgi_table
  expect_equal(tt$TC, 100 - tt$TGI)
  expect_equal(tt$TC_Lower, 100 - tt$TGI_Upper)
  # The dashboard's table shows it.
  shown <- g("window_display_table")(w, "Vehicle")
  expect_true(all(c("T/C (%)", "T/C 95% interval") %in% names(shown)))
})

# ---- F1: mRECIST -------------------------------------------------------------------

# One animal per row of `paths`: volumes on days 0, 4, ..., 24.
resp_data <- function(paths) {
  days <- seq(0, 24, 4)
  do.call(rbind, lapply(seq_along(paths), function(i)
    data.frame(ID = i, Treatment = names(paths)[i], Day = days[seq_along(paths[[i]])],
               Volume = paths[[i]])))
}

test_that("R23.1: the mRECIST response categories and their waterfall are gone", {
  ns <- getNamespaceExports("mouseExperiment")
  expect_false(any(c("response_categories", "plot_response_waterfall") %in% ns))
  expect_false(exists("response_categories", envir = .pkg, inherits = FALSE))
  ui <- paste(as.character(g("tumorGrowthUI")("tg")), collapse = "")
  expect_false(grepl("mRECIST", ui, fixed = TRUE))
  expect_false(grepl("tg_resp_min_day", ui, fixed = TRUE))
  expect_match(ui, 'data-value="Growth Delay"', fixed = TRUE)
})

test_that("F2: growth delay is T - C of Kaplan-Meier medians, times interpolated on the log scale", {
  days <- seq(0, 40, 5)
  mk <- function(arm, rates) do.call(rbind, lapply(seq_along(rates), function(i)
    data.frame(ID = i, Treatment = arm, Day = days, Volume = 100 * exp(rates[i] * days))))
  # Control doubles every ~5 days; the treated arm grows half as fast.
  d <- rbind(mk("Vehicle", c(0.13, 0.14, 0.15, 0.16)), mk("Drug", c(0.065, 0.07, 0.075, 0.08)))
  gd <- growth_delay(d, volume = 1000, reference_group = "Vehicle", n_boot = 200, seed = 1)
  # Each animal reaches 1000 at log(10) / rate exactly (log-linear growth).
  pa <- gd$per_animal
  expect_equal(pa$Time[pa$Treatment == "Vehicle"], log(10) / c(0.13, 0.14, 0.15, 0.16))
  t <- gd$delay_table
  expect_identical(t$Treatment, c("Vehicle", "Drug"))
  # Kaplan-Meier median of four events: survival sits at 0.5 between the
  # second and third, and survfit() takes the midpoint.
  km <- function(r) mean(sort(log(10) / r)[2:3])
  expect_equal(t$Median_Days, c(km(c(0.13, 0.14, 0.15, 0.16)), km(c(0.065, 0.07, 0.075, 0.08))))
  expect_equal(t$Delay_Days[2], t$Median_Days[2] - t$Median_Days[1])
  expect_equal(t$Delay_Pct[2], 100 * t$Delay_Days[2] / t$Median_Days[1])
  expect_true(t$Delay_Lower[2] <= t$Delay_Days[2] && t$Delay_Days[2] <= t$Delay_Upper[2])
  # Reproducible with a seed.
  expect_identical(growth_delay(d, volume = 1000, reference_group = "Vehicle", n_boot = 200,
                                seed = 1)$delay_table, t)
})

test_that("R23.2: the log-rank test against the reference, and the specific growth delay", {
  days <- seq(0, 40, 5)
  mk <- function(arm, rates) do.call(rbind, lapply(seq_along(rates), function(i)
    data.frame(ID = i, Treatment = arm, Day = days, Volume = 100 * exp(rates[i] * days))))
  ctl <- c(0.13, 0.14, 0.15, 0.16)
  d <- rbind(mk("Vehicle", ctl), mk("Slow", c(0.065, 0.07, 0.075, 0.08)),
             mk("Same", c(0.125, 0.135, 0.155, 0.165)))
  gd <- growth_delay(d, volume = 1000, reference_group = "Vehicle", n_boot = 100, seed = 1)
  t <- gd$delay_table
  expect_identical(t$Treatment, c("Vehicle", "Same", "Slow"))
  pa <- gd$per_animal
  lr <- function(arm) {
    s <- pa[pa$Treatment %in% c("Vehicle", arm), ]
    stats::pchisq(survival::survdiff(survival::Surv(Time, Reached) ~ Treatment, data = s)$chisq,
                  1, lower.tail = FALSE)
  }
  expect_true(is.na(t$P_Value[1]))
  expect_equal(t$P_Value[2:3], c(lr("Same"), lr("Slow")))
  expect_equal(t$P_Adjusted[2:3], stats::p.adjust(t$P_Value[2:3], "holm"))
  expect_lt(t$P_Value[3], 0.05)
  expect_gt(t$P_Value[2], 0.05)
  expect_equal(growth_delay(d, volume = 1000, reference_group = "Vehicle", n_boot = 0,
                            p_adjust_method = "none")$delay_table$P_Adjusted[2:3], t$P_Value[2:3])
  # The reference's doubling time: the median of its animals', log(2) / rate.
  expect_equal(gd$doubling_time_reference, stats::median(log(2) / ctl))
  expect_equal(t$SGD, t$Delay_Days / gd$doubling_time_reference)
  expect_equal(t$SGD_Lower[3], t$Delay_Lower[3] / gd$doubling_time_reference)
  expect_identical(gd$p_adjust_method, "holm")
  expect_error(growth_delay(d, volume = 1000, p_adjust_method = "dunnett"), "should be one of")
})

test_that("F2: an arm that never reaches the volume has a median of Inf (not reached)", {
  d <- data.frame(ID = rep(1:6, each = 4), Treatment = rep(c("Vehicle", "Cured"), each = 12),
                  Day = rep(c(0, 7, 14, 21), 6),
                  Volume = c(rep(c(100, 400, 900, 2000), 3), rep(c(100, 80, 60, 40), 3)))
  gd <- growth_delay(d, volume = 1000, reference_group = "Vehicle", n_boot = 50, seed = 1)
  t <- gd$delay_table
  expect_identical(t$Median_Days[2], Inf)
  expect_identical(t$N_Reached, c(3L, 0L))
  expect_identical(t$Delay_Days[2], Inf)
  expect_error(growth_delay(d, volume = -1), "one positive number")
})

# ---- the dashboard ------------------------------------------------------------------

test_that("F2/R23.2: the Growth Delay tab shows and downloads both tables; its R script reproduces them", {
  rv <- shiny::reactiveValues(
    uploaded_data = master_synthetic_data, data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference", tg_p_adjust = "holm",
                      tg_delay_volume = 1000, tg_delay_n_boot = 100, tg_delay_p_adjust = "bonferroni")
    session$setInputs(run_tg_analysis = 1)
    expect_false(isTRUE(rv$tg_results$error))
    head <- paste(as.character(output$tg_delay_table), collapse = "")
    expect_match(head, "p (log-rank)", fixed = TRUE)
    expect_match(head, "Specific growth delay", fixed = TRUE)
    expect_match(paste(as.character(output$tg_delay_animals), collapse = ""), "Outcome", fixed = TRUE)
    groups <- c("tables", "plots", "diagnostics", "data", "report")
    ins <- stats::setNames(rep(list(character(0)), length(groups)), paste0("tg_dl_group_", groups))
    ins$tg_dl_group_tables <- c("growth_delay", "growth_delay_animals")
    ins$tg_dl_group_code <- "script"
    do.call(session$setInputs, ins)
    dir <- tempfile("endpoints_")
    utils::unzip(output$tg_dl_zip, exdir = dir)
    files <- list.files(dir, recursive = TRUE)
    expect_true(all(c("tables/growth_delay.csv", "tables/growth_delay_by_animal.csv",
                      "analysis.R") %in% files))
    expect_false(any(grepl("response", files)))
    script <- readLines(file.path(dir, "analysis.R"))
    expect_false(any(grepl("response_categories", script, fixed = TRUE)))
    expect_true(any(grepl("volume = 1000", script, fixed = TRUE)))
    expect_true(any(grepl('p_adjust_method = "bonferroni"', script, fixed = TRUE)))
    env <- new.env(parent = globalenv())
    old <- setwd(dir); on.exit(setwd(old))
    suppressMessages(suppressWarnings(utils::capture.output(
      sys.source("analysis.R", envir = env, keep.source = FALSE))))
    delay <- utils::read.csv(file.path(dir, "tables", "growth_delay.csv"))
    expect_equal(env$delay$delay_table$Delay_Lower, delay$Delay_Lower)
    expect_equal(env$delay$delay_table$Median_Days, delay$Median_Days)
    expect_equal(env$delay$delay_table$P_Adjusted, delay$P_Adjusted)
    expect_equal(env$delay$delay_table$SGD, delay$SGD)
    animals <- utils::read.csv(file.path(dir, "tables", "growth_delay_by_animal.csv"))
    expect_equal(env$delay$per_animal$Time, animals$Time)
  }))
})

test_that("F2: the delay table says 'not reached', and gives no interval when neither end is reached", {
  t <- data.frame(Treatment = c("Vehicle", "A", "B"), N = 8L, N_Reached = c(8L, 0L, 5L),
                  Median_Days = c(11.7, Inf, 24.6), Delay_Days = c(0, Inf, 12.9),
                  Delay_Lower = c(NA, Inf, 10.5), Delay_Upper = c(NA, Inf, Inf),
                  Delay_Pct = c(0, Inf, 110))
  d <- g("delay_display_table")(t, "Vehicle")
  expect_identical(d$`Median days`, c("11.7", "not reached", "24.6"))
  expect_identical(d[["Delay 95% interval"]], c("\u2014", "\u2014", "[10.5, not reached]"))
  expect_identical(d[["Delay T \u2212 C (days)"]][1], "reference")  # R23.2: the test and the specific growth delay, when the result has them.
  t$SGD <- c(0, Inf, 4.12); t$SGD_Lower <- c(NA, Inf, 3.36); t$SGD_Upper <- c(NA, Inf, Inf)
  t$P_Value <- c(NA, 1e-5, 0.0004); t$P_Adjusted <- c(NA, 2e-5, 0.0004)
  d <- g("delay_display_table")(t, "Vehicle")
  expect_identical(d[["Specific growth delay"]], c("\u2014", "not reached", "4.12"))
  expect_identical(d[["SGD 95% interval"]], c("\u2014", "\u2014", "[3.36, not reached]"))
  expect_identical(d[["p (log-rank)"]][1], "\u2014")
  expect_identical(d[["p (adjusted)"]][3], g("format_pval")(0.0004))
  a <- g("delay_animal_table")(data.frame(Treatment = c("B", "A"), ID = c("1", "2"), Time = c(20.04, 9.96),
                                          Reached = c(0L, 1L)))
  expect_identical(a$Outcome, c("reached", "censored (never reached)"))
  expect_identical(a[["Day (reached or last measured)"]], c(10, 20))
})
