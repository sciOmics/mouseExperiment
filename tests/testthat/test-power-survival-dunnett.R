# Round 22 proposal F6 (CODE_REVIEW.md R22-F, 0.44.0): power for a survival
# endpoint, and Dunnett's correction in the analytic power.

skip_if_not_installed("shiny", minimum_version = "1.5.0")

.pkg <- tryCatch(getNamespace("mouseExperiment"), error = function(e) NULL)
skip_if(is.null(.pkg), "mouseExperiment namespace not loadable")
g <- function(nm) get(nm, envir = .pkg)

test_that("F6: Dunnett's critical values match the published tables", {
  crit <- function(m, df) {
    a <- g("me_dunnett_alpha")(0.05, m, df)
    if (is.finite(df)) stats::qt(1 - a / 2, df) else stats::qnorm(1 - a / 2)
  }
  # Two-sided, alpha 0.05 (Dunnett 1964; e.g. Hochberg & Tamhane, Table 4).
  expect_equal(crit(2, Inf), 2.21, tolerance = 0.005)
  expect_equal(crit(3, 20), 2.54, tolerance = 0.005)
  expect_equal(crit(5, Inf), 2.51, tolerance = 0.005)
  expect_equal(crit(4, 30), 2.58, tolerance = 0.005)   # between the grid's 24 and 40
  expect_identical(g("me_dunnett_alpha")(0.05, 1, 10), 0.05)
  a <- g("me_dunnett_alpha")(0.05, 3, 40)
  expect_true(a > 0.05 / 3 && a < 0.05)
})

test_that("F6: the analytic power's Dunnett N lies between the unadjusted and Bonferroni ones, and reaches the target", {
  run <- function(adj) apriori_power_analysis(effect_size = 1.2, n_groups = 4, alpha = 0.05,
                                              target_power = 0.8, p_adjust_method = adj)
  n <- vapply(c("none", "dunnett", "bonferroni"), function(a) run(a)$scenario_table$Required_N, 1)
  expect_true(n[["none"]] <= n[["dunnett"]] && n[["dunnett"]] <= n[["bonferroni"]])
  st <- run("dunnett")$scenario_table
  pw <- function(k) stats::power.t.test(n = k, delta = 1.2, sd = 1, type = "two.sample",
                                        sig.level = g("me_dunnett_alpha")(0.05, 3, 4 * (k - 1)))$power
  expect_gte(pw(st$Required_N), 0.8)
  expect_lt(pw(st$Required_N - 1), 0.8)
  expect_equal(st$Alpha_Per_Comparison, g("me_dunnett_alpha")(0.05, 3, 4 * (st$Required_N - 1)))
  fp <- apriori_power_analysis(effect_size = 1.2, n_groups = 4, alpha = 0.05, mode = "find_power",
                               n_per_group = st$Required_N, p_adjust_method = "dunnett")
  expect_equal(fp$scenario_table$Achieved_Power, pw(st$Required_N))
})

test_that("F6: survival power follows Schoenfeld: 66 events to halve the hazard at 80 % power", {
  s <- apriori_survival_power(hazard_ratio = 0.5, control_event_prob = 1, alpha = 0.05,
                              target_power = 0.8)
  st <- s$scenario_table
  expect_identical(st$Required_Events, 66)
  expect_identical(st$Required_N, 33L)   # every animal has the event: two per pair
  # With 90 % of controls, the treated arm's share follows proportional hazards.
  s2 <- apriori_survival_power(hazard_ratio = 0.5, control_event_prob = 0.9, alpha = 0.05,
                               target_power = 0.8, dropout_rate = 0.2)
  st2 <- s2$scenario_table
  expect_equal(st2$Treated_Event_Prob, 1 - 0.1^0.5)
  expect_identical(st2$Required_N, as.integer(ceiling(65.35 / (0.9 + 1 - 0.1^0.5))) )
  expect_identical(st2$Enroll_N, ceiling(st2$Required_N / 0.8))
  # The power at that N reaches the target; one fewer does not.
  fp <- apriori_survival_power(0.5, 0.9, alpha = 0.05, mode = "find_power",
                               n_per_group = c(st2$Required_N, st2$Required_N - 1))$scenario_table
  expect_gte(fp$Achieved_Power[1], 0.8)
  expect_lt(fp$Achieved_Power[2], 0.8)
  # A hazard ratio above 1 is as detectable as its reciprocal on the log scale.
  expect_identical(apriori_survival_power(2, 1, alpha = 0.05, target_power = 0.8)$scenario_table$Required_Events, 66)
  expect_error(apriori_survival_power(1, 0.9), "other than 1")
  expect_error(apriori_survival_power(0.5, 0), "in \\(0, 1\\]")
})

test_that("F6: the Power tab plans a survival endpoint, and its R script reproduces the table", {
  rv <- shiny::reactiveValues(uploaded_data = NULL, data_valid = FALSE, session_id = NULL)
  suppressWarnings(shiny::testServer(g("power_server"), args = list(rv = rv), {
    session$setInputs(pw_mode = "apriori_analytic", pw_ap_es_mode = "hazard_ratio",
                      pw_ap_hr = "0.3, 0.5", pw_ap_ctrl_events = 90, pw_ap_n_groups = "3",
                      pw_alpha = "0.05", pw_power = "0.80", pw_ap_p_adjust = "dunnett",
                      pw_ap_n_comparisons = 2, pw_ap_dropout = 10, pw_ap_submode = "find_n")
    session$setInputs(run_pw_analysis = 1)
    res <- rv$ap_results
    expect_identical(res$endpoint, "survival")
    st <- res$scenario_table
    expect_identical(st$Hazard_Ratio, c(0.3, 0.5))
    expect_true(all(st$Required_N > 0))
    expect_match(paste(capture.output(print(output$pw_ap_summary_ui)), collapse = " "),
                 "Dunnett over 2", fixed = TRUE)
    groups <- c("tables", "plots", "diagnostics", "data", "report")
    ins <- stats::setNames(rep(list(character(0)), length(groups)), paste0("pw_dl_group_", groups))
    ins$pw_dl_group_code <- "script"
    do.call(session$setInputs, ins)
    dir <- tempfile("pw_")
    utils::unzip(output$pw_dl_zip, exdir = dir)
    script <- readLines(file.path(dir, "analysis.R"))
    expect_true(any(grepl("apriori_survival_power(", script, fixed = TRUE)))
    env <- new.env(parent = globalenv())
    old <- setwd(dir); on.exit(setwd(old))
    suppressMessages(utils::capture.output(sys.source("analysis.R", envir = env, keep.source = FALSE)))
    expect_equal(env$scenario_table, st, ignore_attr = TRUE)
  }))
})
