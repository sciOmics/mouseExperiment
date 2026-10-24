# Round 22 proposals F1-F3 (CODE_REVIEW.md R22-E, 0.43.0): %T/C beside TGI,
# mRECIST response categories with a waterfall, and tumour growth delay.

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

test_that("F1: the mRECIST categories follow Gao et al.'s thresholds", {
  d <- resp_data(list(
    CR = c(100, 60, 20, 4, 2, 1, 1),            # best -99 %, average below -40 %
    PR = c(100, 80, 60, 45, 40, 40, 40),        # best -60 %
    SD = c(100, 110, 120, 125, 120, 115, 110),  # best +10 %
    PD = c(100, 150, 200, 300, 400, 500, 600),
    late_drop = c(100, 400, 600, 800, 40, 30, 20),  # best -80 %, but average stays high
    short = c(100, 120)))                      # not measured 10 days after baseline
  rc <- response_categories(d)
  cat <- stats::setNames(rc$per_animal$Category, rc$per_animal$Treatment)
  expect_identical(unname(cat[c("CR", "PR", "SD", "PD", "late_drop", "short")]),
                   c("mCR", "mPR", "mSD", "mPD", "mPD", "Not evaluable"))
  pr <- rc$per_animal[rc$per_animal$Treatment == "PR", ]
  expect_equal(pr$Best_Response, -60)
  expect_identical(pr$Best_Response_Day, 16)
  # Running mean from baseline over days 0..t, its minimum over t >= 10.
  dv <- 100 * (c(100, 80, 60, 45, 40, 40, 40) - 100) / 100
  expect_equal(pr$Best_Avg_Response, min((cumsum(dv) / seq_along(dv))[-(1:3)]))
  s <- rc$summary
  expect_identical(s$Response_Rate[s$Treatment == "CR"], 100)
  expect_true(is.na(s$Response_Rate[s$Treatment == "short"]))
})

test_that("F1: animals with reused ear tags stay apart, and the waterfall draws", {
  d <- rbind(resp_data(list(A = c(100, 80, 60, 45, 40, 40, 40))),
             resp_data(list(B = c(100, 150, 200, 300, 400, 500, 600))))
  d$ID <- 1   # the same tag in both arms
  rc <- response_categories(d)
  expect_identical(nrow(rc$per_animal), 2L)
  p <- plot_response_waterfall(rc, cap = 150)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "1 animal above 150% drawn at 150% (largest 200%)", fixed = TRUE)
  expect_s3_class(plot_response_waterfall(rc, fill = "category"), "ggplot")
})

# ---- F2: growth delay ------------------------------------------------------------

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

test_that("F1/F2: Tumor Growth's Response & Delay tab downloads both, and its R script reproduces them", {
  rv <- shiny::reactiveValues(
    uploaded_data = master_synthetic_data, data_valid = TRUE, session_id = NULL,
    settings = list(col_id = "ID", col_day = "Day", time_method = "day", col_volume = "Volume",
                    volume_method = "precalculated", col_treatment = "Treatment", col_cage = "Cage",
                    col_dose = "", control_group = "Vehicle", volume_units = "mm3",
                    randomisation_unit = "mouse"))
  suppressWarnings(shiny::testServer(g("tumor_growth_server"), args = list(rv = rv), {
    session$setInputs(tg_model_type = "lme4", tg_transform = "log", tg_random_effects = "slope",
                      tg_cage_effect = "auto", tg_comparison_type = "reference", tg_p_adjust = "holm",
                      tg_resp_min_day = 10, tg_delay_volume = 1000, tg_delay_n_boot = 100,
                      tg_waterfall_fill = "treatment", tg_waterfall_cap = 200)
    session$setInputs(run_tg_analysis = 1)
    expect_false(isTRUE(rv$tg_results$error))
    expect_false(is.null(output$tg_response_summary))
    expect_false(is.null(output$tg_delay_table))
    groups <- c("tables", "plots", "diagnostics", "data", "report")
    ins <- stats::setNames(rep(list(character(0)), length(groups)), paste0("tg_dl_group_", groups))
    ins$tg_dl_group_tables <- c("responses", "responses_animals", "growth_delay")
    ins$tg_dl_group_plots <- "plot_response_waterfall"
    ins$tg_dl_group_code <- "script"
    do.call(session$setInputs, ins)
    dir <- tempfile("endpoints_")
    utils::unzip(output$tg_dl_zip, exdir = dir)
    files <- list.files(dir, recursive = TRUE)
    expect_true(all(c("tables/response_summary.csv", "tables/response_by_animal.csv",
                      "tables/growth_delay.csv", "plots/response_waterfall.png", "analysis.R") %in% files))
    script <- readLines(file.path(dir, "analysis.R"))
    expect_true(any(grepl("responses <- response_categories(", script, fixed = TRUE)))
    expect_true(any(grepl("volume = 1000", script, fixed = TRUE)))
    env <- new.env(parent = globalenv())
    old <- setwd(dir); on.exit(setwd(old))
    suppressMessages(suppressWarnings(utils::capture.output(
      sys.source("analysis.R", envir = env, keep.source = FALSE))))
    shown <- utils::read.csv(file.path(dir, "tables", "response_summary.csv"), check.names = FALSE)
    expect_equal(env$responses$summary, shown, ignore_attr = TRUE)
    delay <- utils::read.csv(file.path(dir, "tables", "growth_delay.csv"))
    expect_equal(env$delay$delay_table$Delay_Lower, delay$Delay_Lower)
    expect_equal(env$delay$delay_table$Median_Days, delay$Median_Days)
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
  expect_identical(d[["Delay T \u2212 C (days)"]][1], "reference")
})
