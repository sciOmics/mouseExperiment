# =============================================================================
# Round 18 — namespace resolution, and known-answer checks for the surfaces the
# audit had never exercised.
#
# R CMD check found what testthat could not: five `pkg::fn` calls to objects that
# are visible inside a namespace but not exported, so `::` throws. Four were live.
# The most consequential was silent -- influence diagnostics wrapped in tryCatch,
# returning NULL on every run while the dashboard advertised them.
# =============================================================================

# ---- R18.1 namespace calls that could never resolve --------------------------

test_that("R18.1: the five previously-unresolvable calls are gone", {
  # `pkg::fn` requires fn to be EXPORTED. These were visible via imports only:
  #   stats::make.names  (make.names is base)
  #   brms::Gamma        (Gamma is stats)
  #   brms::update       (update is a stats generic)
  #   lme4::influence    (influence is a stats generic; lme4 registers the method)
  #   brms::ebfmi        (not exported by brms at all)
  code <- unlist(lapply(list.files(testthat::test_path("..", "..", "R"),
                                   pattern = "[.]R$", full.names = TRUE),
                        function(f) sub("#.*$", "", readLines(f, warn = FALSE))))
  for (bad in c("stats::make\\.names", "brms::Gamma", "brms::update",
                "lme4::influence", "brms::ebfmi")) {
    expect_false(any(grepl(bad, code)), info = bad)
  }
})

test_that("R18.1: influence diagnostics are produced, not silently NULL", {
  # Two faults kept these NULL on every lme4 run. Fixing the namespace was not
  # enough: influence.merMod refits by case deletion and re-evaluates the model
  # call, whose `data = analysis_df` names a local of the fitting function that
  # is gone by the time a caller holds the result.
  set.seed(5)
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) for (m in 1:8) {
    k <- k + 1L; d <- seq(0, 21, 3)
    rows[[k]] <- data.frame(
      ID = paste0(tx, m), Treatment = tx, Day = d,
      Volume = exp(log(150) + stats::rnorm(1, 0, 0.1) +
                     (if (tx == "Control") 0.12 else 0.07) * d +
                     stats::rnorm(length(d), 0, 0.05)),
      stringsAsFactors = FALSE)
  }
  r <- suppressWarnings(suppressMessages(
    tumor_growth_statistics(do.call(rbind, rows), plots = FALSE, verbose = FALSE)))

  expect_false(is.null(r$diag_cooks_distance))
  expect_false(is.null(r$diag_dfbetas))
  expect_s3_class(r$diag_cooks_distance, "data.frame")
  expect_gt(nrow(r$diag_cooks_distance), 0L)
  expect_true(all(c("Obs", "Cooks_D", "Threshold", "Is_Influential") %in%
                    names(r$diag_cooks_distance)))
  expect_true(all(is.finite(r$diag_cooks_distance$Cooks_D)))
})

# ---- R18.2 percent escaping under roxygen markdown mode ----------------------

test_that("R18.2: no Rd file carries the malformed double-escaped percent", {
  # DESCRIPTION sets Roxygen: list(markdown = TRUE). In markdown mode `\%` is a
  # markdown escape, so roxygen emits a literal backslash plus an escaped
  # percent -- `\\%` -- which Rd cannot parse. The error cascades: every
  # subsequent \item and section header in that file is lost, which is how
  # apriori_power_simulation lost the documentation of ten parameters.
  man <- testthat::test_path("..", "..", "man")
  skip_if_not(dir.exists(man), "man/ not locatable")
  bad <- character(0)
  for (f in list.files(man, pattern = "[.]Rd$", full.names = TRUE)) {
    if (any(grepl("\\\\\\\\%", readLines(f, warn = FALSE)))) bad <- c(bad, basename(f))
  }
  expect_length(bad, 0L)
})

test_that("R18.2: percent stays escaped inside literal Rd macros", {
  # The correction has a boundary: roxygen does NOT escape % inside \code{},
  # \eqn{} or \deqn{}, so there it must remain `\%`. A blanket unescape breaks
  # those -- the percent starts an Rd comment and eats the rest of the line.
  src <- unlist(lapply(list.files(testthat::test_path("..", "..", "R"),
                                  pattern = "[.]R$", full.names = TRUE),
                       function(f) readLines(f, warn = FALSE)))
  rox <- grep("^#'", src, value = TRUE)
  offenders <- grep("\\\\(code|eqn|deqn)\\{[^}]*[^\\\\]%", rox, value = TRUE)
  expect_length(offenders, 0L)
})

# ---- Known-answer verifications ---------------------------------------------

test_that("R18.3: all six volume formulas match hand calculation", {
  d <- data.frame(ID = "m1", Treatment = "A", Day = 0,
                  Length = 20, Width = 10, Height = 8)
  expected <- c(ellipsoid          = pi * 20 * 10^2 / 6,
                modified_ellipsoid = 20 * 10^2 / 2,
                ellipsoid_3axis    = pi * 20 * 10 * 8 / 6,
                cylinder           = pi * 10^2 * 20 / 4,
                sphere             = pi * 10^3 / 6,
                box                = 20 * 10 * 8)
  for (f in names(expected)) {
    got <- suppressWarnings(suppressMessages(calculate_volume(
      d, length_column = "Length", width_column = "Width",
      height_column = "Height", formula = f)))$Volume
    expect_equal(got, unname(expected[[f]]), tolerance = 1e-9, info = f)
  }
  # Height absent: the 3-axis form falls back to width, documented behaviour.
  d2 <- d[, setdiff(names(d), "Height")]
  got <- suppressWarnings(suppressMessages(
    calculate_volume(d2, formula = "ellipsoid_3axis")))$Volume
  expect_equal(got, pi * 20 * 10 * 10 / 6, tolerance = 1e-9)
})

test_that("R18.3: doubling time equals ln(2)/rate exactly", {
  for (r0 in c(0.05, 0.10, 0.20)) {
    rows <- list(); k <- 0L
    for (tx in c("Control", "DrugA")) for (m in 1:5) {
      k <- k + 1L; d <- seq(0, 28, 4)
      rr <- if (tx == "Control") r0 else r0 / 2
      rows[[k]] <- data.frame(ID = paste0(tx, m), Treatment = tx, Day = d,
                              Volume = 100 * exp(rr * d), stringsAsFactors = FALSE)
    }
    res <- suppressWarnings(suppressMessages(
      tumor_doubling_time(do.call(rbind, rows))))
    dt <- res$doubling_time[res$Treatment == "Control"][1]
    expect_equal(dt, log(2) / r0, tolerance = 1e-6)
    expect_equal(res$growth_rate[res$Treatment == "Control"][1], r0, tolerance = 1e-6)
  }
})

test_that("R18.3: AUC equals the exact trapezoid integral", {
  # V(t) = 100 + 10t on [0,20] integrates to 100*20 + 10*20^2/2 = 4000.
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) for (m in 1:5) {
    k <- k + 1L; d <- seq(0, 20, 5)
    slope <- if (tx == "Control") 10 else 5
    rows[[k]] <- data.frame(ID = paste0(tx, m), Treatment = tx, Day = d,
                            Volume = 100 + slope * d, stringsAsFactors = FALSE)
  }
  a <- suppressWarnings(suppressMessages(tumor_auc_analysis(do.call(rbind, rows))))
  s <- a$auc_summary
  expect_equal(s$Mean_AUC[s$Treatment == "Control"], 4000, tolerance = 1e-6)
  expect_equal(s$Mean_AUC[s$Treatment == "DrugA"],   3000, tolerance = 1e-6)
})

test_that("R18.3: volume-to-mass respects the declared units", {
  # 1 cm3 == 1000 mm3; at density 1.0 g/cm3 both are 1.0 g.
  expect_equal(as.numeric(volume_to_mass(1000, tumor_density = 1.0,
                                         volume_units = "mm3")), 1.0, tolerance = 1e-9)
  expect_equal(as.numeric(volume_to_mass(1, tumor_density = 1.0,
                                         volume_units = "cm3")), 1.0, tolerance = 1e-9)
  # Getting this wrong scales the tumour-weight correction by 1000 (R3.30).
  expect_identical(detect_volume_units(c(150, 800, 2000)), "mm3")
  expect_identical(detect_volume_units(c(0.15, 0.8, 2.0)), "cm3")
})

test_that("R18.3: cage ICC matches the variance components it is built from", {
  skip_if_not_installed("lme4")
  set.seed(11)
  n_cage <- 40; n_per <- 8
  ce <- stats::rnorm(n_cage, 0, 1)
  rows <- list(); k <- 0L
  for (ci in seq_len(n_cage)) for (m in seq_len(n_per)) {
    k <- k + 1L
    rows[[k]] <- data.frame(Cage = paste0("C", ci),
                            Treatment = if (ci <= n_cage / 2) "A" else "B",
                            y = ce[ci] + stats::rnorm(1, 0, 1),
                            stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  fit <- suppressMessages(lme4::lmer(y ~ Treatment + (1 | Cage), data = d))
  got <- cage_icc(fit, "Cage")
  vc <- as.data.frame(lme4::VarCorr(fit))
  manual <- vc$vcov[vc$grp == "Cage"] / sum(vc$vcov)
  expect_equal(as.numeric(got$ICC), manual, tolerance = 1e-8)
  expect_gt(got$ICC, 0); expect_lt(got$ICC, 1)
})

test_that("R18.3: GAMM recovers a plateauing trajectory", {
  skip_if_not_installed("gamm4")
  # Control grows throughout; treated plateaus after day 12. A single slope per
  # arm must average across the two phases.
  set.seed(67)
  days <- seq(0, 28, 2)
  tru <- function(tx, d) if (tx == "Control") log(150) + 0.10 * d else
    log(150) + ifelse(d <= 12, 0.10 * d, 0.10 * 12 + 0.005 * (d - 12))
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) for (m in 1:10) {
    k <- k + 1L
    rows[[k]] <- data.frame(ID = paste0(tx, m), Treatment = tx, Day = days,
      Volume = exp(tru(tx, days) + stats::rnorm(1, 0, 0.08) +
                     stats::rnorm(length(days), 0, 0.05)),
      stringsAsFactors = FALSE)
  }
  g <- suppressWarnings(suppressMessages(tumor_growth_statistics(
    do.call(rbind, rows), model_type = "gam", plots = FALSE, verbose = FALSE)))

  expect_false(is.null(g$gam_k_check))
  expect_false(is.null(g$gam_concurvity))
  pc <- g$pairwise_comparisons
  expect_true("Day" %in% names(pc))
  # The effect must be absent early and present late -- that ordering is the
  # whole point of fitting a smoother rather than one slope.
  early <- pc$P_Value_Adjusted[pc$Day == min(pc$Day)][1]
  late  <- pc$P_Value_Adjusted[pc$Day == max(pc$Day)][1]
  expect_gt(early, 0.05)
  expect_lt(late, 0.01)
})

test_that("R18.3: the bivariate toxicity metric recovers simulated weight loss", {
  set.seed(13)
  rows <- list(); k <- 0L
  for (tx in c("Control", "DrugA")) for (m in 1:8) {
    k <- k + 1L; d <- seq(0, 21, 3)
    drop <- if (tx == "DrugA") 0.12 else 0.02
    rows[[k]] <- data.frame(
      ID = paste0(tx, m), Treatment = tx, Day = d,
      Weight = 25 * (1 - drop * d / 21) + stats::rnorm(length(d), 0, 0.05),
      Volume = exp(log(150) + (if (tx == "Control") 0.12 else 0.07) * d),
      stringsAsFactors = FALSE)
  }
  e <- suppressWarnings(suppressMessages(efficacy_toxicity_bivariate(
    do.call(rbind, rows), volume_column = "Volume", weight_column = "Weight",
    time_column = "Day", treatment_column = "Treatment", id_column = "ID",
    reference_group = "Control", adjust_tumor_weight = FALSE)))
  pg <- e$per_group
  # Simulated losses were 2 % and 12 %.
  expect_equal(pg$Toxicity_Mean[pg$Treatment == "Control"], 2, tolerance = 0.5)
  expect_equal(pg$Toxicity_Mean[pg$Treatment == "DrugA"],  12, tolerance = 0.5)
  # TGI by hand: 1 - exp(0.07*21)/exp(0.12*21) = 64.9 %.
  expect_equal(pg$Efficacy_Mean[pg$Treatment == "DrugA"], 64.9, tolerance = 1)
})
