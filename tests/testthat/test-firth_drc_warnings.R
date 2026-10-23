# =============================================================================
# The two warnings in the maintainer's dashboard notes (2026-10-22; v0.33.1,
# CODE_REVIEW.md R20-X): Firth's profile-likelihood p-values that did not
# converge, and optim()'s L-BFGS-B tolerance warning on every dose-response fit.
# Each test fails on v0.33.0.
# =============================================================================

# Six arms of 8, shaped like the dashboard's All Features demo: every control
# animal has an event and two arms have none, which is when Firth is used. With
# coxphf's defaults (maxit = 50, maxstep = 0.5) the p-values of the two
# event-free arms do not converge.
firth_study <- function() {
  set.seed(4)
  arm <- rep(c("Vehicle", "A", "B", "C", "D", "E"), each = 8)
  ev  <- c(rep(1, 8), stats::rbinom(8, 1, 0.6), stats::rbinom(8, 1, 0.1), rep(0, 8),
           stats::rbinom(8, 1, 0.6), rep(0, 8))
  t   <- ifelse(ev == 1, round(stats::runif(48, 10, 30)), 35)
  t[arm == "Vehicle"] <- round(stats::runif(8, 12, 22))
  data.frame(ID = seq_along(arm), Treatment = arm, Day = t, Event = ev)
}
surv_fit <- function(d) {
  w <- character(0)
  res <- withCallingHandlers(
    survival_statistics(d, time_column = "Day", event_column = "Event",
                        treatment_column = "Treatment", id_column = "ID",
                        reference_group = "Vehicle", verbose = FALSE),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
  list(res = res, warnings = w)
}

test_that("R20-X: Firth's model converges for arms without events, without a warning", {
  skip_if_not_installed("coxphf")
  d <- firth_study()
  f <- surv_fit(d)
  expect_identical(f$res$method_used, "coxphf")
  # v0.33.0: "Convergence in estimating profile likelihood CI or p-values not
  # attained for all variables", and C and E fell back to log-rank p-values.
  expect_false(any(grepl("profile likelihood", f$warnings, fixed = TRUE)))
  r <- f$res$results
  non_ref <- r$Group != "Vehicle"
  expect_false(anyNA(r$P_Value_Unadjusted[non_ref]))
  expect_true(all(r$Note[non_ref] == ""))
  # The p-values are Firth's own profile-likelihood ones.
  a <- data.frame(time = d$Day, status = d$Event,
                  group = factor(d$Treatment, levels = c("Vehicle", "A", "B", "C", "D", "E")))
  m <- coxphf::coxphf(survival::Surv(time, status) ~ group, data = a, maxit = 1000, maxstep = 0.05)
  expect_equal(r$P_Value_Unadjusted[match(c("A", "B", "C", "D", "E"), r$Group)],
               unname(m$prob), tolerance = 1e-8)
})

test_that("R20-X: a fallback p-value keeps its precision, and the table says which test gave it", {
  skip_if_not_installed("coxphf")
  d <- firth_study()
  real <- coxphf::coxphf
  # Force one arm's p-value to fail, as it can for harder data.
  local_mocked_bindings(coxphf = function(...) { m <- real(...); m$prob[3] <- NA; m },
                        .package = "coxphf")
  r <- surv_fit(d)$res$results
  expect_match(r$Note[r$Group == "C"], "pairwise log-rank", fixed = TRUE)
  pair <- d[d$Treatment %in% c("Vehicle", "C"), ]
  lr <- survival::survdiff(survival::Surv(Day, Event) ~ Treatment, data = pair)
  # v0.33.0 rounded it to 4 decimals (here 0, and 1e-4 at most).
  expect_equal(r$P_Value_Unadjusted[r$Group == "C"],
               stats::pchisq(lr$chisq, df = 1, lower.tail = FALSE))
  expect_gt(r$P_Value_Unadjusted[r$Group == "C"], 0)
})

test_that("R20-X: a dose-response fit no longer warns about L-BFGS-B tolerances", {
  set.seed(3)
  doses <- c(0, 5, 10, 25, 50, 100, 200)
  df <- do.call(rbind, lapply(doses, function(dz) data.frame(
    ID = paste0("d", dz, "_", 1:8), Treatment = paste0("D", dz), Dose = dz, Day = 21,
    Volume = pmax(1, 100 + 900 / (1 + (max(dz, 1e-9) / 25)^1.5) + stats::rnorm(8, 0, 60)))))
  w <- character(0)
  withCallingHandlers(suppressMessages(dose_response_statistics(
    df, dose_column = "Dose", treatment_column = "Treatment", volume_column = "Volume",
    time_column = "Day", id_column = "ID", reference_group = "D0",
    treatments = setdiff(unique(df$Treatment), "D0"))),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("L-BFGS-B uses 'factr'", w, fixed = TRUE)))
  # Only that warning is muffled.
  expect_silent(me_muffle_lbfgsb_tol(warning(
    "method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'")))
  expect_warning(me_muffle_lbfgsb_tol(warning("something else")), "something else")
})
